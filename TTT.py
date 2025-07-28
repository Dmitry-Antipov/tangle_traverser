#!/usr/bin/env python3
import sys
import argparse
import logging
import math
import os
from src.path_optimizer import PathOptimizer
from src.path_supplementary import get_gaf_string
from src.alignment_scorer import AlignmentScorer
from src.logging_utils import setup_logging
from src.node_id_mapper import NodeIdMapper
from src.MIP_optimizer import MIPOptimizer
from src.input_parsing import (
    parse_gfa, parse_gaf, read_tangle_nodes, read_coverage_file, 
    coverage_from_graph, verify_coverage, calculate_median_coverage,
    identify_boundary_nodes, clean_tips, reverse_complement
)
from src.graph_transformation import (
    get_canonical_rc_vertex,
    create_dual_graph, create_multi_dual_graph,
    get_traversable_subgraph, write_multiplicities
)

# Create global instances
node_id_mapper = NodeIdMapper()
# Add this length (max) from the nodes neighboring tangle to fasta for better alignment
UNIQUE_BORDER_LENGTH = 200000

def optimize_paths(multi_graph, boundary_nodes, original_graph, num_initial_paths, max_iterations, early_stopping_limit, alignment_scorer: AlignmentScorer):
    """Optimize Eulerian paths."""
    best_path = None
    best_score = -1
    logging.info(f"Starting optimization with {num_initial_paths} initial paths, max {max_iterations} iterations per path.")
    rc_vertex_map = {}
    for vertex in multi_graph.nodes():
        rc_vertex_map[vertex] = get_canonical_rc_vertex(vertex, original_graph, node_id_mapper)
    for seed in range(num_initial_paths):
        logging.info(f"Generating initial path with seed {seed}.")
        #TODO: shouldn't it be outside cycle?
        subgraph_to_traverse, start_vertex = get_traversable_subgraph(multi_graph, boundary_nodes, original_graph, seed, node_id_mapper)

        if not subgraph_to_traverse:
            logging.warning(f"No Eulerian path found for seed {seed}.")
            continue
        pathOptimizer = PathOptimizer(subgraph_to_traverse, start_vertex, seed, node_id_mapper, rc_vertex_map)

        current_path = pathOptimizer.get_path()
        current_score = alignment_scorer.score_corasick(current_path)
        logging.info(f"Initial path score for seed {seed}: {current_score}.")

        iterations_since_improvement = 0
        for i in range(max_iterations):
            if iterations_since_improvement >= early_stopping_limit:
                logging.info(f"Early stopping for seed {seed} on iteration {i} after {iterations_since_improvement} iterations without improvement.")
                break
            
            new_path = pathOptimizer.get_random_change(seed * max_iterations + i)
            new_score = alignment_scorer.score_corasick(new_path)
            if new_score > current_score:
                logging.info(f"Improved score for seed {seed} at iteration {i}: {current_score} -> {new_score}.")
                current_path = new_path
                pathOptimizer.set_path(new_path)
                current_score = new_score
                iterations_since_improvement = 0
            else:
                iterations_since_improvement += 1
        
        logging.info(f"Final score for seed {seed}: {current_score}.")
        if current_score > best_score:
            logging.info(f"New best path found for seed {seed} with score {current_score}.")
            best_path = current_path
            best_score = current_score
    
    logging.info(f"Optimization completed. Best score: {best_score}.")
    pathOptimizer.get_synonymous_changes(best_path, alignment_scorer)
    return best_path, best_score

# Only UNIQUE_BORDER_LENGTH (=200K) suffix/prefix for border unique nodes used in fasta, but its HPC anyways...
def output_path(best_path, original_graph, output_fasta, output_gaf):
    aux = -1
    if node_id_mapper.has_name("AUX"):
        aux_id = node_id_mapper.get_id_for_name("AUX")
        for i in range(len(best_path)):
            if abs(best_path[i].original_node) == aux_id:
                aux = i
                logging.debug(f"Found AUX at position {aux}")
                break
    if aux > 0:
        paths = [best_path[:aux - 1], best_path[aux + 1:]]
    else:
        paths = [best_path]

    gaf_file = open(output_gaf, 'w')
    count = 0
    for path in paths:
        gaf_file.write(f"traversal_{count}\t{get_gaf_string(path, node_id_mapper)}\n")
        count += 1
    
    with open(output_fasta, 'w') as fasta_file:
        count = 0
        for path in paths:

            contig_sequence = ""
            last_node = None
            for i in range(len(path)):
                edge_id = path[i].original_node
                node_id = abs(edge_id)
                orientation = edge_id > 0

                # Retrieve the sequence from the graph
                node_sequence = original_graph.nodes[node_id]['sequence']
                if node_sequence == "*":
                    logging.warning("Provided noseq assembly graph, no HPC fasta output possible")
                    return
                if not orientation:
                    # Reverse complement the sequence if orientation is negative
                    node_sequence = reverse_complement(node_sequence)

                if i == 0:
                    overlap = max (0, len(node_sequence) - UNIQUE_BORDER_LENGTH)                
                else:                
                    overlap = original_graph.get_edge_data(last_node, edge_id)['overlap']

                if i == len(best_path) - 1:                            
                    contig_sequence += node_sequence[overlap:min(len(node_sequence),UNIQUE_BORDER_LENGTH)]
                else:
                    contig_sequence += node_sequence[overlap:]
                last_node = edge_id

            # Write the contig to the FASTA file
            fasta_file.write(f">traversal_{count}\n")
            fasta_file.write(f"{contig_sequence}\n")
            count += 1
    logging.info("Fasta output finished")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Solve for integer multiplicities in a GFA tangle graph based on coverage.")
    parser.add_argument("--gfa", required=False, dest="gfa_file", help="Path to the GFA file.")
    parser.add_argument("--alignment", required=False, help="Path to a file with graphaligner alignment")
    parser.add_argument("--outdir", required=True, type=str, help="Output directory for all result files (will be created if it doesn't exist)")
    parser.add_argument("--verkko-output", required=False, type=str, help="Path to dir with verkko results. Ovewrites --gfa, --alignment, --coverage-file with standart paths")
    
    parser.add_argument("--coverage-file", help="Path to a file with node coverages (node-id coverage). If not provided, coverage will be filled from the GFA file.")
    parser.add_argument("--median-unique", type=float, help="Median coverage for unique nodes.")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Set the logging level (default: INFO).")
    parser.add_argument("--tangle-file", help="Path to a file listing nodes in the tangle (one per line).")
    parser.add_argument("--tangle-node", type=str, help="Node ID to construct tangle as connected component of short edges around it")
    parser.add_argument("--boundary-nodes", type=str, help="Path to a file listing boundary node pairs, tab-separated (required for 2-2 tangles).")
    parser.add_argument("--tangle-length-cutoff", type=int, default=500000, help="Length cutoff for tangle detection, default 500K")
    parser.add_argument("--num-initial-paths", type=int, default=10, help="Number of initial paths to generate (default: 10).")
    parser.add_argument("--max-iterations", type=int, default=100000, help="Maximum iterations for path optimization (default: 100000).")
    parser.add_argument("--early-stopping-limit", type=int, default=15000, help="Early stopping limit for optimization (default: 15000).")
    #TODO: for quality 0 use random of alignments with highest score?
    parser.add_argument("--quality-threshold", type=int, default=20, help="Alignments with quality less than this will be filtered out, default 20")
    parser.add_argument("--basename", required=False, default="traversal", type=str, help="Basename for most of the output files, default `traversal`")
    args = parser.parse_args()
    if not args.verkko_output  and (not args.gfa_file  or not args.alignment ):
        #logging not initialized yet
        sys.stderr.write("Either --verkko-output OR both --gfa and --alignment are required options\n")
        return

    if args.verkko_output: 
        args.coverage_file = os.path.join(args.verkko_output, "2-processGraph", "unitig-unrolled-hifi-resolved.ont-coverage.csv")
        args.gfa_file = os.path.join(args.verkko_output, "2-processGraph", "unitig-unrolled-hifi-resolved.gfa")
        args.alignment = os.path.join(args.verkko_output, "3-align", "alns-ont.gaf")
    return args

def main():    
    args = parse_arguments()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args)
    logging.debug(f"args: {args}")
    logging.info("Reading files...")

    #TODO: separate function to clean Z connection, save them somewhere
    original_graph = parse_gfa(args.gfa_file, node_id_mapper)
    if args.coverage_file:
        cov = read_coverage_file(args.coverage_file, node_id_mapper)
    else:
        cov = coverage_from_graph(original_graph)
    # Verifying that coverage matches the graph
    # For rare cases coverage may be missing for some nodes, will update with median then 
    verify_coverage(cov, original_graph, node_id_mapper)

    tangle_nodes = read_tangle_nodes(args, original_graph, node_id_mapper)    
    clean_tips(tangle_nodes, original_graph, node_id_mapper)
    nor_nodes = {abs(node) for node in tangle_nodes}
    
    boundary_nodes = identify_boundary_nodes(args, original_graph, tangle_nodes, node_id_mapper)

    used_nodes = nor_nodes.copy()
    for b in boundary_nodes:
        used_nodes.add(abs(b))
        used_nodes.add(abs(boundary_nodes[b]))
   

    median_unique_range = calculate_median_coverage(args, nor_nodes, original_graph, cov, boundary_nodes, node_id_mapper)    
    median_unique = math.sqrt(median_unique_range[0] * median_unique_range[1])
    filtered_alignment_file = os.path.join(args.outdir, f"{args.basename}.q{args.quality_threshold}.used_alignments.gaf")
    alignments = parse_gaf(args.alignment, used_nodes, filtered_alignment_file, args.quality_threshold, node_id_mapper)
    alignment_scorer = AlignmentScorer(alignments)
    logging.info("Starting multiplicity counting...")
    # TODO: instead of a_values we just use coverage
    mip_optimizer = MIPOptimizer(node_id_mapper)
    equations, nonzeros, a_values, boundary_values = mip_optimizer.generate_MIP_equations(tangle_nodes, nor_nodes, cov, median_unique, original_graph, boundary_nodes, directed=True)
    best_solution = mip_optimizer.solve_MIP(equations, nonzeros, boundary_values, a_values, median_unique_range)
    
    # Define output filenames based on the output directory
    output_csv = os.path.join(args.outdir, args.basename + ".multiplicities.csv")
    output_fasta = os.path.join(args.outdir, args.basename + ".hpc.fasta")
    output_gaf = os.path.join(args.outdir, args.basename + ".gaf")

    write_multiplicities(output_csv, best_solution, cov, node_id_mapper)

    # Now all sequence is stored in edges, junctions are new vertices
    dual_graph = create_dual_graph(original_graph, node_id_mapper)
    multi_graph = create_multi_dual_graph(dual_graph, best_solution, tangle_nodes, boundary_nodes, original_graph, node_id_mapper)

    best_path, best_score = optimize_paths(multi_graph, boundary_nodes, original_graph, args.num_initial_paths, args.max_iterations, args.early_stopping_limit, alignment_scorer)
    logging.info("Path optimizing finished.")
    if best_path:
        best_path_str = get_gaf_string(best_path, node_id_mapper)
        logging.info(f"Found traversal\t{best_path_str}")   

        # Output FASTA file
        logging.info(f"Writing best path to {output_fasta} and gaf to {output_gaf}")
        output_path(best_path, original_graph, output_fasta, output_gaf)




if __name__ == "__main__":
    main()
