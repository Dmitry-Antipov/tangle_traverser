#!/usr/bin/env python3

import sys
import os
import argparse
import logging
import networkx as nx
from .logging_utils import log_assert
from .node_mapper import NodeIdMapper

#allowed median coverage range in range [median_unique/MEDIAN_COVERAGE_VARIATION, median_unique * MEDIAN_COVERAGE_VARIATION]
DETECTED_MEDIAN_COVERAGE_VARIATION = 1.5
GIVEN_MEDIAN_COVERAGE_VARIATION = 1.2

def reverse_complement(sequence):
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return sequence.translate(complement)[::-1]

def rc_node(node):
    """RC nodes stored as negative"""
    return -node

def parse_gfa(file_path, node_mapper):
    """
    Parse a GFA file and construct a directed graph using networkx.DiGraph.
    Nodes are stored as integers (positive for '+' orientation, negative for '-' orientation).
    
    :param file_path: Path to the GFA file
    :param node_mapper: NodeIdMapper instance
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

def parse_gaf(gaf_file, interesting_nodes, filtered_file, quality_threshold, node_mapper):
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

def node_to_tangle(directed_graph, length_cutoff, target_node, node_mapper):
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

def clean_tips(tangle_nodes, directed_graph, node_mapper):
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

def read_tangle_nodes(args, original_graph, node_mapper):
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
        tangle_nodes = node_to_tangle(original_graph, args.tangle_length_cutoff, target_node=node_mapper.parse_node_id(args.tangle_node), node_mapper=node_mapper)
    else:
        logging.error("Either --tangle_file or --tangle_node must be provided.")
        sys.exit(1)

    with open(os.path.join(args.outdir, "nodes.txt"), 'w') as f:
        for node in tangle_nodes:
            #only forward nodes to file
            if node > 0:
                f.write(f"{node_mapper.node_id_to_name[node]}\n")

    return tangle_nodes

def read_coverage_file(coverage_file, node_mapper):
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

def verify_coverage(cov, original_graph, node_mapper):
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

def is_forward_unique(oriented_node, original_graph):
    """Checks if an oriented node has exactly one outgoing edge."""
    return len(list(original_graph.successors(oriented_node))) == 1

def calculate_median_unique_coverage(nor_nodes, original_graph, cov, min_b, node_mapper):
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

def calculate_median_coverage(args, nor_nodes, original_graph, cov, boundary_nodes, node_mapper):
    """Calculate or use provided median unique coverage."""
    min_b = 1000000
    all_boundary = list(boundary_nodes.keys()) + list(boundary_nodes.values())
    for b in all_boundary:
        logging.info (f"{node_mapper.node_id_to_name_safe(b)} {cov[abs(b)]}")
        min_b = min(min_b, cov[abs(b)])
    if args.median_unique is not None:
        return [args.median_unique/GIVEN_MEDIAN_COVERAGE_VARIATION, args.median_unique * GIVEN_MEDIAN_COVERAGE_VARIATION]
    calculated_median = calculate_median_unique_coverage(nor_nodes, original_graph, cov, min_b, node_mapper)
    if calculated_median is not None:
        return [calculated_median/DETECTED_MEDIAN_COVERAGE_VARIATION, calculated_median* DETECTED_MEDIAN_COVERAGE_VARIATION]
    else:
        res = [min_b/DETECTED_MEDIAN_COVERAGE_VARIATION, min_b* DETECTED_MEDIAN_COVERAGE_VARIATION]
        logging.warning(f"Failed to calculate median unique coverage for tangle. Using coverage based on neighbours {res} but better provide it manually with--median-unique.")
        return res

def identify_boundary_nodes(args, original_graph, tangle_nodes, node_mapper):
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


