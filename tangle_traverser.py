#!/usr/bin/env python3
import sys
import pulp
import argparse
import logging
import statistics
import networkx as nx
import random
import math
import subprocess
import cProfile
import pstats
import time
import os
import dataclasses
from src.path_optimizer import PathOptimizer
from src.path_supplementary import EdgeDescription, get_gaf_path, get_gaf_string, rc_path
from src.alignment_scorer import AlignmentScorer
from src.logging_utils import setup_logging, log_assert
from src.node_mapper import NodeIdMapper
from src.MIP_optimizer import MIPOptimizer

# Create a global node mapper instance
node_mapper = NodeIdMapper()

# Create a global MIP optimizer instance
mip_optimizer = MIPOptimizer(node_mapper)

#allowed median coverage range in range [median_unique/MEDIAN_COVERAGE_VARIATION, median_unique * MEDIAN_COVERAGE_VARIATION]
DETECTED_MEDIAN_COVERAGE_VARIATION = 1.5
GIVEN_MEDIAN_COVERAGE_VARIATION = 1.2

#Add this length (max) from the nodes neighboring tangle to fasta for better alignment
UNIQUE_BORDER_LENGTH = 200000

#Multiplicities detection magic is here
#TODO: length-based weights?

       
def parse_gfa(file_path):
    """
    Parse a GFA file and construct a directed graph using networkx.DiGraph.
    Nodes are stored as integers (positive for '+' orientation, negative for '-' orientation).
    
    :param file_path: Path to the GFA file
    :return: networkx.DiGraph representing the graph
    """
    original_graph = nx.DiGraph()

    with open(file_path, 'r') as gfa_file:
        for line in gfa_file:
            # Skip header lines
            if line.startswith('#') or not line.strip():
                continue
            if line.startswith('S'):
                parts = line.strip().split('\t')
                node_id = node_mapper.parse_node_id(parts[1])
                node_mapper.add_mapping(node_id, parts[1])  # Map node ID to its string name
                #node_mapper.node_id_to_name[-node_id] = '-' + parts[1]  # Map reverse complement

                # Extract node length from LN:i:<length> or length of provided string
                length = len(parts[2])
                cov = -1
                for part in parts:
                    if part.startswith('LN:i:'):
                        length = int(part.split(':')[2])
                    elif part.startswith('ll:'):
                        cov = float(part.split(':')[2])
                original_graph.add_node(node_id, length=length, sequence=parts[2], coverage = cov)
                original_graph.add_node(-node_id, length=length, sequence=reverse_complement(parts[2]), coverage = cov)

    #links and segments can be interlaced
    with open(file_path, 'r') as gfa_file:
        for line in gfa_file:
            if line.startswith('L'):
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue  # Skip malformed lines

                from_node = node_mapper.parse_node_id(parts[1])
                if parts[2] == '-':
                    from_node = -from_node

                to_node = node_mapper.parse_node_id(parts[3])
                if parts[4] == '-':
                    to_node = -to_node

                # Extract overlap size if available
                overlap_size = 0
                if len(parts) > 5 and parts[5].endswith('M'):
                    try:
                        overlap_size = int(parts[5][:-1])
                    except ValueError:
                        logging.warning(f"Invalid overlap size in line: {line.strip()}")

                # Add to the graph with overlap size as an edge attribute
                original_graph.add_edge(from_node, to_node, overlap=overlap_size)
                original_graph.add_edge(rc_node(to_node), rc_node(from_node), overlap=overlap_size)

    # Symmetrization of junctions
    #TODO: should be reimagined.
    #with missing BD: A+D=B+C, A>= C, B >= D
    non_transitive_junctions = 0
    for start_node in list(original_graph.nodes):
        end_nodes = set(original_graph.successors(start_node))
        if len(end_nodes) == 0:
            continue
        start_nodes = set()
        for e in end_nodes:
            for n in original_graph.successors(-e):
                start_nodes.add(-n)
        #We may need multiple iterations of symmetrization, this should be enough
        for _ in range (10):
            for n1 in start_nodes:
                for n2 in original_graph.successors(n1):
                    if n2 not in original_graph.successors(start_node):
                        # Calculate the minimum overlap size for the symmetrized link
                        overlap_counted = False
                        for e in end_nodes:
                            if rc_node(n1) in original_graph.successors(rc_node(e)):
                                overlap_sizes = []
                                overlap_sizes.append(original_graph[start_node][e].get('overlap', 0))
                                overlap_sizes.append(original_graph[rc_node(e)][rc_node(n1)].get('overlap', 0))
                                overlap_sizes.append(original_graph[n1][n2].get('overlap', 0))
                                min_overlap = min(overlap_sizes) if overlap_sizes else 0
                                original_graph.add_edge(start_node, n2, overlap=min_overlap)
                                logging.debug(f"Non-transitive junction! Adding {node_mapper.node_id_to_name_safe(start_node)} -> {node_mapper.node_id_to_name_safe(n2)}, overlap {min_overlap}")
                                non_transitive_junctions += 1
                                overlap_counted = True
                                break
                        if not overlap_counted:
                            logging.warning(f"Non-transitive junction failed to get overlap. Adding {node_mapper.node_id_to_name_safe(start_node)} -> {node_mapper.node_id_to_name_safe(n2)}, overlap 0")
                            original_graph.add_edge(start_node, n2, overlap=0)
    logging.info(f"Tuned non-transitive junctions: {non_transitive_junctions}")
    return original_graph

#RC nodes stored as negative
def rc_node(node):
    return -node

#Possibly we need to remove dead-ends and very low covered nodes

#Detect unique nodes in tangle and count median
def is_forward_unique(oriented_node, original_graph):
    """Checks if an oriented node has exactly one outgoing edge."""
    return len(list(original_graph.successors(oriented_node))) == 1

def calculate_median_unique_coverage(nor_nodes, original_graph, cov, min_b):
    """
    Calculates the median coverage of nodes in nor_nodes that are structurally unique.
    A node is structurally unique if both its '+' and '-' orientations satisfy
    the conditions checked by is_forward_unique.
    """
    unique_coverages = []
    unique_node_ids = set()

    logging.debug(f"Identifying structurally unique nodes from {len(nor_nodes)} tangle nodes.")
    for node_id in nor_nodes:
        is_plus_unique = is_forward_unique(node_id, original_graph)
        is_minus_unique = is_forward_unique(-node_id, original_graph)

        if is_plus_unique and is_minus_unique:
            if node_id in cov and cov[node_id] < min_b * GIVEN_MEDIAN_COVERAGE_VARIATION:
                coverage = float(cov[node_id])
                if coverage == 0:
                    logging.info(f"Node {node_mapper.node_id_to_name_safe(node_id)} has zero coverage ,excluding from unique coverage calculation.")
                    continue
                unique_coverages.append([coverage, original_graph.nodes[node_id].get('length', 0), node_id])
                unique_node_ids.add(node_id)
                logging.debug(f"Node {node_mapper.node_id_to_name_safe(node_id)} is unique. Coverage: {coverage}")
            elif cov[node_id] >= min_b * GIVEN_MEDIAN_COVERAGE_VARIATION:
                logging.debug(f"Node {node_mapper.node_id_to_name_safe(node_id)} looks structurally unique but coverage {cov[node_id]} is higher than borders {min_b} * variation {GIVEN_MEDIAN_COVERAGE_VARIATION}")
            else:
                logging.warning(f"Structurally unique node {node_mapper.node_id_to_name_safe(node_id)} not found in coverage file.")
        #else:
            #logging.debug(f"Node {node_mapper.node_id_to_name_safe(node_id)} is not structurally unique (+ unique: {is_plus_unique}, - unique: {is_minus_unique}).")


    if not unique_coverages:
        logging.warning("No structurally unique nodes found in the tangle to calculate median coverage.")
        return None
    unique_coverages.sort()
    total_len = 0
    for coverage, length, node_id in unique_coverages:
        total_len += length
    cur_len = 0
    for ind in range(len(unique_coverages)):
        cur_len += unique_coverages[ind][1]
        if cur_len * 2 >= total_len:
            median_cov = unique_coverages[ind][0]
            break
    unique_debug = []
    #TODO: global mapping to restore original names
    for u in unique_node_ids:
        unique_debug.append(node_mapper.node_id_to_name[u])
    logging.info(f"Found {len(unique_node_ids)} structurally unique nodes: {sorted(unique_debug)}")
    logging.info(f"Calculated median coverage of unique nodes: {median_cov:.2f}")
    logging.debug(f"Unique coverages: {unique_coverages}")
    return median_cov
   
def get_canonical_nodepair(oriented_node1, oriented_node2, original_graph):
    """Returns a canonical lexicographically smallest link for a junction"""
    to_node = oriented_node2
    from_node = oriented_node1

    for node in original_graph.successors(oriented_node1):
        if node < to_node:
            to_node = node
    
    for node in original_graph.successors(rc_node(to_node)):
        if rc_node(node) < from_node:
            from_node = rc_node(node)
    if not (to_node in original_graph.successors(from_node)):
        logging.warning(f"Irregular junction {node_mapper.node_id_to_name_safe(from_node)} {node_mapper.node_id_to_name_safe(to_node)} {node_mapper.node_id_to_name_safe(oriented_node1)} {node_mapper.node_id_to_name_safe(oriented_node2)}")
        logging.warning(f"{node_mapper.node_id_to_name_safe(from_node)} -> {[node_mapper.node_id_to_name_safe(n) for n in original_graph.successors(from_node)]}")
        logging.warning(f"{node_mapper.node_id_to_name_safe(rc_node(to_node))} -> {[node_mapper.node_id_to_name_safe(n) for n in original_graph.successors(rc_node(to_node))]}")        
    return (from_node, to_node)

def get_canonical_rc_vertex(v, original_graph):
    """Returns a canonical reverse connection vertex for a given vertex"""
    return get_canonical_nodepair(-v[1], -v[0], original_graph)

#Transform to dual graph, vertices = junctions, edges = old nodes
#TODO: vertices should be not exactly junctions but junction+nodes since z-connections are legit

def create_dual_graph(original_graph:nx.MultiDiGraph):
    dual_graph = nx.MultiDiGraph()
    canonical_edges_set = set() # To store unique canonical edges (dual nodes)
    logging.info("Creating dual graph representation (nodes = canonical connections) using NetworkX...")

    # Iterate through each oriented node 'v' which acts as the junction
    for v in sorted(original_graph.nodes):
        # Determine the canonical edge C1 leading INTO node v
        predecessors = list(original_graph.successors(rc_node(v)))
        if not predecessors:
            # v is a source tip (no incoming edges for this orientation)
            C1 = ("TIP", v)
        else:
            # Pick any predecessor 'u' to define the canonical incoming edge
            # Note: get_canonical_edge expects (from, to), so use (pred, v)
            u = rc_node(predecessors[0])
            C1 = get_canonical_nodepair(u, v, original_graph)
        canonical_edges_set.add(C1)

        # Determine the canonical edge C2 leading OUT FROM node v
        successors = list(original_graph.successors(v))
        if not successors:
            # v is a sink tip (no outgoing edges for this orientation)
            C2 = (v, "TIP")
        else:
            # Pick any successor 'w' to define the canonical outgoing edge
            # Note: get_canonical_edge expects (from, to), so use (v, succ)
            w = successors[0]
            C2 = get_canonical_nodepair(v, w, original_graph)
        canonical_edges_set.add(C2)

        # Add the edge in the dual graph representing the transition through v
        # Ensure nodes C1 and C2 exist before adding the edge (add_nodes_from handles this later too)
        dual_graph.add_node(C1)
        dual_graph.add_node(C2)
        dual_graph.add_edge(C1, C2, original_node=v)
        logging.debug(f"Added edge from {C1} to {C2} in dual graph, original node: {node_mapper.node_id_to_name_safe(v)}")

    logging.info(f"Total oriented nodes in original GFA structure: {len(original_graph.nodes)}")
    logging.info(f"Total unique canonical edges (dual nodes): {len(canonical_edges_set)}")
    logging.info(f"Total nodes in final dual graph: {dual_graph.number_of_nodes()}") # Should match set size
    logging.info(f"Total edges (connections through original nodes) in dual graph: {dual_graph.number_of_edges()}")
    #logging.debug(f"Edges of dual graph: {list(dual_graph.edges())}")
    return dual_graph

def is_tangle_vertex (C, tangle_set):
    return C[0] in tangle_set or C[1] in tangle_set



#edge of multiplicity X -> X multiedges multiplicity 1
def create_multi_dual_graph(dual_graph: nx.DiGraph, multiplicities: dict, tangle_nodes: set, boundary_nodes: dict, G):
    multi_dual_graph = nx.MultiDiGraph()
    logging.info("Creating multi-dual graph from dual graph and multiplicities...")

    # Ensure all nodes from the original dual graph are present
    for node in dual_graph.nodes():
        if is_tangle_vertex(node, tangle_nodes):
            multi_dual_graph.add_node(node)

    #adding start and sink nodes
    incomings = set(boundary_nodes.keys())
    outgoings = set(boundary_nodes.values())
    for b in boundary_nodes:
        incomings.add(-boundary_nodes[b])
        outgoings.add(-b)
    edges_added = 0

    for b in boundary_nodes.keys():
        for next_node in G.successors(b):
            fw_edge = get_canonical_nodepair(b, next_node, G)                                
            logging.debug(f"start {node_mapper.node_id_to_name_safe(b)} {node_mapper.node_id_to_name_safe(next_node)}")
            start = (0, b)
            multi_dual_graph.add_node(start)
            multi_dual_graph.add_edge(start, fw_edge, original_node=b, key = f"{b}_{edges_added}")
            edges_added += 1
            break 
    for b in boundary_nodes.values():
        for next_node in G.predecessors(b):
            bw_edge = get_canonical_nodepair(next_node, b, G)                                
            logging.debug(f"start {node_mapper.node_id_to_name_safe(next_node)} {node_mapper.node_id_to_name_safe(b)}")
            start = (b, 0)
            multi_dual_graph.add_node(start)
            multi_dual_graph.add_edge(bw_edge, start, original_node=b, key = f"{b}_{edges_added}")
            edges_added += 1
            break
    
    # Iterate through edges of the dual graph and add them to the multi-graph
    # based on the multiplicity of the original node they represent.
    for u, v, data in dual_graph.edges(data=True):
        if (not u in multi_dual_graph.nodes()) or (not v in multi_dual_graph.nodes()):
            continue
        original_node_oriented = data.get('original_node')
        if not original_node_oriented:
            logging.warning(f"Edge ({u}, {v}) in dual graph is missing 'original_node' attribute. Skipping.")
            continue
        
        original_node_base = original_node_oriented #abs(original_node_oriented)

        multiplicity = multiplicities[original_node_base]

        if multiplicity < 0:
            logging.error(f"Negative multiplicity {multiplicity} for node {node_mapper.node_id_to_name_safe(original_node_base)}. Treating as 0 for edge ({u} -> {v}).")
            exit(0)
        # Add the edge 'multiplicity' times to the multi-dual graph
        logging.debug(f"Adding {multiplicity} multiedges for {node_mapper.node_id_to_name_safe(original_node_oriented)}")
        for _ in range(multiplicity):
            # Add edge with the original node attribute
            multi_dual_graph.add_edge(u, v, original_node=original_node_oriented, key = str(original_node_oriented) + "_"+str(edges_added))
            edges_added += 1

    logging.info(f"Created multi-dual graph with {multi_dual_graph.number_of_nodes()} nodes.")
    logging.info(f"Added {edges_added} edges to multi-dual graph based on multiplicities (original dual graph had {dual_graph.number_of_edges()} unique edges).")
    for v in multi_dual_graph.nodes():
        logging.debug(f"{v}: in {multi_dual_graph.in_degree[v] } out {multi_dual_graph.out_degree[v] }")    

    return multi_dual_graph

#Supplementary for Euler path search - moved to path_optimizer.py






#not_oriented border nodes, currently should be just one in one out
#TODO: non-random starts based on most support?
def get_traversable_subgraph(multi_dual_graph: nx.MultiDiGraph, border_nodes, original_graph, seed):
    random.seed(seed)
    start_vertices = []

    end_vertices = []
    for n in border_nodes:
        #border tips encoding
        start_vertices.append((0, n))
        end_vertices.append((n, 0))

    logging.info(f"Start and end vertices in the graph {start_vertices}, {end_vertices}")
    #adding fake links 
    

    border_nodes_count = len(border_nodes)
    # Only 1-1 or 2-2 tangles for now
    log_assert(border_nodes_count == 1 or border_nodes_count == 2, f"Only 1-1 or 2-2 tangles are supported")
    log_assert(len(start_vertices) == border_nodes_count and len(end_vertices) == border_nodes_count, f"Start and end vertices count mismatch: {len(start_vertices)} vs {border_nodes_count} or {len(end_vertices)} vs {border_nodes_count}")
    start_vertices.sort()
    start_vertex = start_vertices[0]
    start_vertex_node = start_vertex[1]
    matching_end_vertex_node = border_nodes[start_vertex_node]
    matching_end_vertex = (matching_end_vertex_node, 0)
    for e in multi_dual_graph.edges(keys=True):
        logging.debug(f"Edge {e}")

    # If there are 4 border nodes, we need to find the correct end vertex and add the auxiliary connection
    if border_nodes_count == 2:
        aux_node_str = "AUX"
        aux_int_id = node_mapper.parse_node_id(aux_node_str)
        multi_dual_graph.add_edge(matching_end_vertex, start_vertices[1], original_node=aux_int_id, key=f"{aux_int_id}_0")
        logging.info(f"Added auxiliary edge from {node_mapper.node_id_to_name_safe(matching_end_vertex_node)} to {node_mapper.node_id_to_name_safe(start_vertices[1][1])}")
        reachable_verts = nx.descendants(multi_dual_graph, start_vertex)
        # Add the start_node itself
        reachable_verts.add(start_vertex)

        # Recreate reachable subgraph after aux edge added
        reachable_subgraph = multi_dual_graph.subgraph(reachable_verts)

    unreachable_edges = set()
    for _ in range (100):
        reachable_verts = nx.descendants(multi_dual_graph, start_vertex)
        # Add the start_node itself
        reachable_verts.add(start_vertex)
        unreachable_verts = set(multi_dual_graph.nodes()) - reachable_verts
        unreachable_edges.clear()
        for e in multi_dual_graph.edges(keys=True):
            if e[0] in unreachable_verts:
                unreachable_edges.add(e)
                logging.debug(e)

        logging.info(f"Clearing unreachable edges: iteration {_} {len(unreachable_edges)} present")
        if len(unreachable_edges) == 0:
            break
        
        for e in unreachable_edges:
            logging.debug(f"Reversing edge: {e}")
            data = multi_dual_graph.get_edge_data(e[0], e[1], key=e[2])
            logging.debug(f"Data {data}")
            multi_dual_graph.remove_edge(e[0], e[1], key = e[2])
            if e[2][0] == '-':
                new_key = e[2][1:]
            else:
                new_key = '-' + e[2]
            multi_dual_graph.add_edge(get_canonical_rc_vertex(e[1], original_graph), get_canonical_rc_vertex(e[0], original_graph), original_node=-int(data['original_node']), key = new_key)

    logging.debug(f"After transformation")
    for e in multi_dual_graph.edges(keys=True):
        logging.debug(f"Edge {e}")
    log_assert(len(unreachable_edges) == 0, f"Unreachable edges are still present: {unreachable_edges}")
    reachable_subgraph = multi_dual_graph.subgraph(reachable_verts)
    return reachable_subgraph, start_vertex

def parse_gaf(gaf_file, interesting_nodes, filtered_file, quality_threshold):
    res = []
    if filtered_file:
        out_file = open(filtered_file, 'w')
    with open(gaf_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            path_id = parts[0]
            if path_id == "node":
                continue

            quality_score = int(parts[11])  # Assuming quality score is in the 11th column
            if quality_score < quality_threshold:
                continue

            nodes = []
            node = ""
            for char in parts[5].strip():
                if char in "<>":
                    if node:
                        nodes.append(node)
                    node = char
                else:
                    node += char
            if node:
                nodes.append(node)
            filtered_nodes = nodes.copy()
            nodes = []
            good = True
            for fnode in filtered_nodes:            
                int_node = node_mapper.parse_node_id(fnode[1:])
                if not (int_node in interesting_nodes):
                    good = False
                    break
                if fnode[0] == "<":
                    int_node = -int_node
                nodes.append(int_node)   
            #reverse_complement would be added in AlignmentScorer
            if good and len(nodes) > 1:
                res.append(nodes)
                if filtered_file:
                    out_file.write(line)
    return  res 

def node_to_tangle(directed_graph, length_cutoff, target_node):
    """
    tangle of short edges containing target_node
    """
    indirect_graph = nx.Graph()

    # Add edges between nodes that meet the length cutoff
    for u, v in directed_graph.edges():
        if directed_graph.nodes[u].get('length', float('inf')) < length_cutoff and \
           directed_graph.nodes[v].get('length', float('inf')) < length_cutoff:
            indirect_graph.add_edge(u, v)

    # Check if the target node exists in the indirect graph
    if (target_node not in indirect_graph) and (-target_node not in indirect_graph):
        raise ValueError(f"Target node {target_node} is not in the indirect graph.")

    # Get the connected component containing the target node        
    connected_component = nx.node_connected_component(indirect_graph, target_node)
    logging.info(f"Total {len(connected_component)} tangle nodes in connected component")
    for u in connected_component:
        logging.debug(f"Tangle node {node_mapper.node_id_to_name_safe(abs(u))}")
    rc_component = nx.node_connected_component(indirect_graph, -target_node)
    connected_component.update(rc_component)
    
    return connected_component

def clean_tips (tangle_nodes, directed_graph):
    #removing tips from tangle and graph
    changed= True
    while changed:
        changed = False
        to_erase = []
        for n in tangle_nodes:
            if directed_graph.out_degree(n) == 0 or directed_graph.in_degree(n) == 0:
                changed = True
                to_erase.append(n)
        for n in to_erase:
            logging.info(f"Cleaning {node_mapper.node_id_to_name_safe(n)} from tangle")
            tangle_nodes.remove(n)
            directed_graph.remove_node(n)
            


def read_tangle_nodes(args, original_graph):
    """Read or construct tangle nodes."""
    if args.tangle_file:
        tangle_nodes = set()
      #  nor_nodes = set()
        with open(args.tangle_file, 'r') as f:
            for line in f:
                node_str = line.strip()
                if node_str:
                    node_id = node_mapper.parse_node_id(node_str)
                    tangle_nodes.add(node_id)
                    tangle_nodes.add(-node_id)
                  #  nor_nodes.add(node_id)
    elif args.tangle_node:
        tangle_nodes = node_to_tangle(original_graph, args.tangle_length_cutoff, target_node=node_mapper.parse_node_id(args.tangle_node))
    else:
        logging.error("Either --tangle_file or --tangle_node must be provided.")
        sys.exit(1)

    with open(os.path.join(args.outdir, "nodes.txt"), 'w') as f:
        for node in tangle_nodes:
            #only forward nodes to file
            if node > 0:
                f.write(f"{node_mapper.node_id_to_name[node]}\n")

    return tangle_nodes

def read_coverage_file(coverage_file):
    """node_id\tnode_cov"""
    logging.info(f"Reading coverage from {coverage_file}")

    cov = {}
    with open(coverage_file, 'r') as f:
        for line in f:
            arr = line.strip().split()
            #Usual cov csv format starts with node*
            if len(arr) >= 2 and arr[0] != "node":
                node_id = node_mapper.parse_node_id(arr[0])
                cov[node_id] = float(arr[1])
    return cov

def coverage_from_graph(assembly_graph):
    logging.info(f"Using coverage provided with the assembly graph")
    cov = {}
    for id in assembly_graph.nodes():
        if id < 0:
            continue
        cov[id] = assembly_graph.nodes[id]['coverage']
        if cov[id] < 0:
            logging.error("Coverage file not provided, failed to get coverage from the assembly graph gfa")
            exit()
    return cov

def verify_coverage(cov, original_graph):
    graph_nodes = set(original_graph.nodes())
    for node_id in cov.keys():
        if abs(node_id) not in graph_nodes:
            logging.error(f"Node {node_mapper.node_id_to_name_safe(node_id)} from coverage file is not present in the assembly graph, please check that coverage corresponds to the same graph")
            exit(1)    
    missing_length = 0
    total_length = 0
    nodes = []
    for node_id in graph_nodes:
        if node_id < 0:
            continue

        if node_id not in cov:
            logging.info(f"Node {node_mapper.node_id_to_name_safe(node_id)} from assembly graph is not present in the coverage file")
            missing_length += original_graph.nodes[node_id].get('length', 0)
        else:
            nodes.append([original_graph.nodes[node_id].get('length', 0), cov[node_id]])
        total_length += original_graph.nodes[node_id].get('length', 0)
        
    if missing_length > 0:
        if missing_length * 10 >= total_length:
            logging.error(f"Too many missing nodes in coverage file: {missing_length} out of {total_length}, please check that coverage corresponds to the same graph")
            exit(1)
    else:
        nodes.sort(key=lambda x: x[0])
        cur_len = 0
        estimated_unique_coverage = 0
        for length, coverage in nodes:
            cur_len += length
            if cur_len * 2 >= total_length:
                estimated_unique_coverage = coverage
                break

        logging.info(f"Total length of missing nodes in whole graph: {missing_length} out of {total_length}, median coverage {estimated_unique_coverage} will be used")
        for node_id in graph_nodes:
            if node_id not in cov:
                cov[node_id] = estimated_unique_coverage
    

def calculate_median_coverage(args, nor_nodes, original_graph, cov, boundary_nodes):
    """Calculate or use provided median unique coverage."""
    min_b = 1000000
    all_boundary = list(boundary_nodes.keys()) + list(boundary_nodes.values())
    for b in all_boundary:
        logging.info (f"{node_mapper.node_id_to_name_safe(b)} {cov[abs(b)]}")
        min_b = min(min_b, cov[abs(b)])
    if args.median_unique is not None:
        return [args.median_unique/GIVEN_MEDIAN_COVERAGE_VARIATION, args.median_unique * GIVEN_MEDIAN_COVERAGE_VARIATION]
    calculated_median = calculate_median_unique_coverage(nor_nodes, original_graph, cov, min_b)
    if calculated_median is not None:
        return [calculated_median/DETECTED_MEDIAN_COVERAGE_VARIATION, calculated_median* DETECTED_MEDIAN_COVERAGE_VARIATION]
    else:
        res = [min_b/DETECTED_MEDIAN_COVERAGE_VARIATION, min_b* DETECTED_MEDIAN_COVERAGE_VARIATION]
        logging.warning(f"Failed to calculate median unique coverage for tangle. Using coverage based on neighbours {res} but better provide it manually with--median-unique.")
        return res

def write_multiplicities(output_file, solution, cov):
    with open(output_file, 'w') as out_file:
        out_file.write("node\tcoverage\tmult\n")
        if len (solution) > 0:
            for node_id in cov.keys():                    
                mult_value = solution.get(node_id, "X")
                rev_mult_value = solution.get(-node_id, "X")
                if mult_value != "X" and rev_mult_value != "X":
                    mult_value = int(mult_value) + int(rev_mult_value)
                cov_value = cov.get(node_id, "N/A")
                out_file.write(f"{node_mapper.node_id_to_name[node_id]}\t{cov_value}\t{mult_value}\n")
    logging.info(f"Wrote multiplicity solutions to {output_file}")


def identify_boundary_nodes(args, original_graph, tangle_nodes):
    if args.boundary_nodes:
        boundary_nodes = {}
        for line in open(args.boundary_nodes):
            # Parse the line and extract boundary node pairs
            # map incoming->matching outgoing
            parts = line.strip().split()
            if len(parts) == 2:
                node1 = node_mapper.parse_node_id(parts[0])
                is_incoming = False
                for next_node in original_graph.successors(node1):
                    if next_node in tangle_nodes:
                        is_incoming = True
                        break
                if not is_incoming:
                    node1 = -node1
                node2 = node_mapper.parse_node_id(parts[1])
                is_outgoing = False
                for prev in original_graph.predecessors(node2):
                    if prev in tangle_nodes:
                        is_outgoing = True
                        break
                if not is_outgoing:
                    node2 = -node2
                boundary_nodes[node1] = node2
        return boundary_nodes
    else:
        out_boundary_nodes = set()
        in_boundary_nodes = set()
        for first in original_graph.nodes:
            for second in original_graph.successors(first):
                if first in tangle_nodes and second not in tangle_nodes:
                    out_boundary_nodes.add(second)
                elif second in tangle_nodes and first not in tangle_nodes:
                    in_boundary_nodes.add(first)
        log_assert(len(out_boundary_nodes) == 2 and len(in_boundary_nodes) == 2, 
                    f"Autodetection works only for 1-1 tangles, detected out_boundary: {[node_mapper.node_id_to_name_safe(n) for n in out_boundary_nodes]}, in_boundary: {[node_mapper.node_id_to_name_safe(n) for n in in_boundary_nodes]}. Specify boundary node pairs manually")
        res = {}
        start = max(in_boundary_nodes)
        end = max(out_boundary_nodes)
        log_assert(abs(start) != abs(end), 
                    f"Start and end boundary nodes should be different, got {node_mapper.node_id_to_name_safe(start)} and {node_mapper.node_id_to_name_safe(end)}")
        res[start] = end
        return res


# Update the optimize_paths function to use command-line options for default parameters
def optimize_paths(multi_graph, boundary_nodes, original_graph, num_initial_paths, max_iterations, early_stopping_limit, alignment_scorer: AlignmentScorer):
    """Optimize Eulerian paths."""
    best_path = None
    best_score = -1
    logging.info(f"Starting optimization with {num_initial_paths} initial paths, max {max_iterations} iterations per path.")
    rc_vertex_map = {}
    for vertex in multi_graph.nodes():
        rc_vertex_map[vertex] = get_canonical_rc_vertex(vertex, original_graph)
    for seed in range(num_initial_paths):
        logging.info(f"Generating initial path with seed {seed}.")
        #TODO: shouldn't it be outside cycle?
        subgraph_to_traverse, start_vertex = get_traversable_subgraph(multi_graph, boundary_nodes, original_graph, seed)

        if not subgraph_to_traverse:
            logging.warning(f"No Eulerian path found for seed {seed}.")
            continue
        pathOptimizer = PathOptimizer(subgraph_to_traverse, start_vertex, seed, node_mapper, rc_vertex_map)

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

#Only UNIQUE_BORDER_LENGTH (=200K) suffix/prefix for border unique nodes used in fasta, but its HPC anyways...
def output_path(best_path, original_graph, output_fasta, output_gaf):

    aux = -1
    if node_mapper.has_name("AUX"):
        aux_id = node_mapper.get_id_for_name("AUX")
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
        gaf_file.write(f"traversal_{count}\t{get_gaf_string(path, node_mapper)}\n")
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

def reverse_complement(sequence):
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return sequence.translate(complement)[::-1]

    
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
    
    # Create output directory if it doesn't exist

    #TODO: save it somewhere, some connections added while parsing    
    original_graph = parse_gfa(args.gfa_file)
    if args.coverage_file:
        cov = read_coverage_file(args.coverage_file)
    else:
        cov = coverage_from_graph(original_graph)
    #verifying that coverage matches the graph
    #for rare cases coverage may be missing for some nodes, will update with median then 
    verify_coverage(cov, original_graph)

    tangle_nodes = read_tangle_nodes(args, original_graph)    
    #TODO: Do we really need it?
    clean_tips(tangle_nodes, original_graph)
    nor_nodes = {abs(node) for node in tangle_nodes}
    
    boundary_nodes = identify_boundary_nodes(args, original_graph, tangle_nodes)

    used_nodes = nor_nodes.copy()
    for b in boundary_nodes:
        used_nodes.add(abs(b))
        used_nodes.add(abs(boundary_nodes[b]))
   

    median_unique_range = calculate_median_coverage(args, nor_nodes, original_graph, cov, boundary_nodes)    
    median_unique = math.sqrt(median_unique_range[0] * median_unique_range[1])
    filtered_alignment_file = os.path.join(args.outdir, f"{args.basename}.q{args.quality_threshold}.used_alignments.gaf")
    alignments = parse_gaf(args.alignment, used_nodes, filtered_alignment_file, args.quality_threshold)
    alignment_scorer = AlignmentScorer(alignments)
    #Shit is hidden here
    logging.info("Starting multiplicity counting...")
    #TODO: instead of a_values we just use coverage
    equations, nonzeros, a_values, boundary_values = mip_optimizer.generate_MIP_equations(tangle_nodes, nor_nodes, cov, median_unique, original_graph, boundary_nodes, directed=True)
    #Failed to generate suboptimal solutions yet
    best_solution = mip_optimizer.solve_MIP(equations, nonzeros, boundary_values, a_values, median_unique_range)
    
    # Define output filenames based on the output directory
    output_csv = os.path.join(args.outdir, args.basename + ".multiplicities.csv")
    output_fasta = os.path.join(args.outdir, args.basename + ".hpc.fasta")
    output_gaf = os.path.join(args.outdir, args.basename + ".gaf")

    # Write multiplicities to CSV
    write_multiplicities(output_csv, best_solution, cov)

    
    #Now all sequence is stored in edges, junctions are new vertices
    dual_graph = create_dual_graph(original_graph)
    multi_graph = create_multi_dual_graph(dual_graph, best_solution, tangle_nodes, boundary_nodes, original_graph)

    best_path, best_score = optimize_paths(multi_graph, boundary_nodes, original_graph, args.num_initial_paths, args.max_iterations, args.early_stopping_limit, alignment_scorer)
    logging.info("Path optimizing finished.")
    if best_path:
        best_path_str = get_gaf_string(best_path, node_mapper)
        logging.info(f"Found traversal\t{best_path_str}")   

        
        
        # Output FASTA file
        logging.info(f"Writing best path to {output_fasta} and gaf to {output_gaf}")
        output_path(best_path, original_graph, output_fasta, output_gaf)


if __name__ == "__main__":
    main()
