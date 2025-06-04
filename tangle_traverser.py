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
import ahocorasick

# Add a global map for node ID to string node names
#only positive ids (unoriented)
node_id_to_name = {}
name_to_node_id = {}
last_enumerated_node = 0
#allowed median coverage range in range [median_unique/MEDIAN_COVERAGE_VARIATION, median_unique * MEDIAN_COVERAGE_VARIATION]
DETECTED_MEDIAN_COVERAGE_VARIATION = 1.5
GIVEN_MEDIAN_COVERAGE_VARIATION = 1.2

#Add this length (max) from the nodes neighboring tangle to fasta for better alignment
UNIQUE_BORDER_LENGTH = 200000

#Multiplicities detection magic is here
#TODO: length-based weights?
def solve_MIP(equations, nonzeros, boundary_values, coverages, unique_coverage_range):
    #last element in equation - amount of boundary nodes

    prob = pulp.LpProblem("Minimize_Deviation", pulp.LpMinimize)    
    # Define multiplicity variables    

    keys = set()
    abs_keys = set()
    for eq in equations:
        keys.update(eq[0])
        keys.update(eq[1])
        for ind in (0, 1):
            for j in eq[ind]:
                abs_keys.add(abs(j))
    for key in abs_keys:
        log_assert(key in keys and -key in keys, f"One of the keys {key} or {-key} is not in equations {equations}")
        log_assert(key in coverages , f"Key {key} is not in coverages")

    #TODO: we can check for hairpin presence just here.
    x_vars = {i: pulp.LpVariable(f"x_{i}", cat="Integer") for i in keys}

    # Constraints: x_i = x_j + x_k + x_l
    # for lhs, rhs in equations.items():
    # prob += x_vars[lhs] == sum(x_vars[j] for j in rhs)
    for eq in equations:
        prob += sum(x_vars[j] for j in eq[0])  == sum(x_vars[j] for j in eq[1]) 

    for x_var in x_vars:
        prob += x_vars[x_var] >= 0
        if x_var > 0 and x_var in nonzeros:
            prob += x_vars[x_var] + x_vars[-x_var] >= 1

    #boundary nodes set
    for x_var in boundary_values:
        prob += x_vars[x_var] == boundary_values[x_var]


    # Define continuous variables for absolute deviation |x_i - a_i|
    d_vars = {i: pulp.LpVariable(f"d_{i}", lowBound=0, cat="Continuous") for i in coverages}   
    inv_unique_coverage = pulp.LpVariable("uniqueness_multiplier", lowBound=1/unique_coverage_range[1], upBound=1/unique_coverage_range[0], cat="Continuous")
    for i in abs_keys:
        prob += d_vars[i] >= (x_vars[i] + x_vars[-i]) - coverages[i] * inv_unique_coverage
        prob += d_vars[i] >= coverages[i] * inv_unique_coverage - x_vars[i] - x_vars[-i]
    # Objective function: minimize sum of absolute deviations
    objective = pulp.lpSum(d_vars[i] for i in abs_keys)
    prob += objective

    logging.debug("Linear programming problem:")
    logging.debug(f"Objective: {prob.objective}")
    for constraint in prob.constraints.values():
        logging.debug(f"Constraint: {constraint}")
    for var in x_vars.values():
        logging.debug(f"Variable: {var.name}, LowBound: {var.lowBound}, Cat: {var.cat}")

    # MIP magic
    prob.solve()
    result = {}    
    
    if pulp.LpStatus[prob.status] != "Optimal":
        logging.warning("MIP did not find an optimal solution.")
    result = {i: int(pulp.value(x_vars[i])) for i in keys}
    
    logging.info(f"Found multiplicities from MIP {result}")
    score = pulp.value(objective)  
    detected_coverage = 1/pulp.value(inv_unique_coverage)
    logging.info (f"Unique coverage for the best MIP solution {detected_coverage}, solution score {score}")
    if abs(detected_coverage - unique_coverage_range[0]) < 0.01 or abs(detected_coverage - unique_coverage_range[1]) < 0.01:
        logging.info (f"Warning, detected best coverage is close to the allowed unique coverage borders{unique_coverage_range}")                                                                            
    return result
       
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
                node_id = parse_node_id(parts[1])
                node_id_to_name[node_id] = parts[1]  # Map node ID to its string name
                #node_id_to_name[-node_id] = '-' + parts[1]  # Map reverse complement

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

                from_node = parse_node_id(parts[1])
                if parts[2] == '-':
                    from_node = -from_node

                to_node = parse_node_id(parts[3])
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
                                logging.debug(f"Non-transitive junction! Adding {start_node} -> {n2}, overlap {min_overlap}")
                                non_transitive_junctions += 1
                                overlap_counted = True
                                break
                        if not overlap_counted:
                            logging.warning(f"Non-transitive junction failed to get overlap. Adding {start_node} -> {n2}, overlap 0")
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
                unique_coverages.append(coverage)
                unique_node_ids.add(node_id)
                logging.debug(f"Node {node_id} is unique. Coverage: {coverage}")
            elif cov[node_id] >= min_b * GIVEN_MEDIAN_COVERAGE_VARIATION:
                logging.debug(f"Node {node_id} looks structurally unique but coverage {cov[node_id]} is higher than borders {min_b} * variation {GIVEN_MEDIAN_COVERAGE_VARIATION}")
            else:
                logging.warning(f"Structurally unique node {node_id} not found in coverage file.")
        #else:
            #logging.debug(f"Node {node_id} is not structurally unique (+ unique: {is_plus_unique}, - unique: {is_minus_unique}).")


    if not unique_coverages:
        logging.warning("No structurally unique nodes found in the tangle to calculate median coverage.")
        return None
    unique_coverages.sort()
    median_cov = statistics.median(unique_coverages)
    unique_debug = []
    #TODO: global mapping to restore original names
    for u in unique_node_ids:
        unique_debug.append(node_id_to_name[u])
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
        logging.warning(f"Irregular junction {from_node} {to_node} {oriented_node1} {oriented_node2}")
        logging.warning(f"{from_node} -> {list(original_graph.successors(from_node))}")
        logging.warning(f"{rc_node(to_node)} -> {list(original_graph.successors(rc_node(to_node)))}")        
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
        logging.debug(f"Added edge from {C1} to {C2} in dual graph, original node: {v}")

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
            logging.debug(f"start {b} {next_node}")
            start = (0, b)
            multi_dual_graph.add_node(start)
            multi_dual_graph.add_edge(start, fw_edge, original_node=b, key = f"{b}_{edges_added}")
            edges_added += 1
            break 
    for b in boundary_nodes.values():
        for next_node in G.predecessors(b):
            bw_edge = get_canonical_nodepair(next_node, b, G)                                
            logging.debug(f"start {next_node} {b}")
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
            logging.error(f"Negative multiplicity {multiplicity} for node {original_node_base}. Treating as 0 for edge ({u} -> {v}).")
            exit(0)
        # Add the edge 'multiplicity' times to the multi-dual graph
        logging.debug(f"Adding {multiplicity} multiedges for {original_node_oriented}")
        for _ in range(multiplicity):
            # Add edge with the original node attribute
            multi_dual_graph.add_edge(u, v, original_node=original_node_oriented, key = str(original_node_oriented) + "_"+str(edges_added))
            edges_added += 1

    logging.info(f"Created multi-dual graph with {multi_dual_graph.number_of_nodes()} nodes.")
    logging.info(f"Added {edges_added} edges to multi-dual graph based on multiplicities (original dual graph had {dual_graph.number_of_edges()} unique edges).")
    for v in multi_dual_graph.nodes():
        logging.debug(f"{v}: in {multi_dual_graph.in_degree[v] } out {multi_dual_graph.out_degree[v] }")    

    return multi_dual_graph

#Supplementary for Euler path search
def next_element(temp_graph, current):
    edges = list(temp_graph.out_edges(current, data=True, keys=True))    
    # Randomly choose the next edge based on seed
    random.shuffle(edges)
    chosen_edge = edges[0]  # Take the first edge after shuffling    
    # Extract edge data
    u, v, key, data = chosen_edge
    temp_graph.remove_edge(u, v, key)
    return ((u, v, int(key.split("_")[0])), v)


#not_oriented border nodes, currently should be just one in one out
#TODO: non-random starts based on most support?
def get_traversing_eulerian_path(multi_dual_graph: nx.MultiDiGraph, border_nodes, original_graph, seed):
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

    unreachable_edges = set()
    for _ in range (2):
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
    reachable_end_vertices = []
    for v in reachable_subgraph.nodes():
        if v in end_vertices:
            reachable_end_vertices.append(v)

    # If there are 4 border nodes, we need to find the correct end vertex and add the auxiliary connection
    if border_nodes_count == 2:
        aux_node_str = "AUX"
        aux_int_id = parse_node_id(aux_node_str)  
        multi_dual_graph.add_edge(matching_end_vertex, start_vertices[1], original_node=aux_int_id, key = f"{aux_int_id}_0")


        logging.info(f"Added auxiliary edge {matching_end_vertex_node} -> {start_vertices[1]}")
        reachable_verts = nx.descendants(multi_dual_graph, start_vertex)
        # Add the start_node itself
        reachable_verts.add(start_vertex)

        # Recreate reachable subgraph after aux edge added
        reachable_subgraph = multi_dual_graph.subgraph(reachable_verts)

    # Eulerian path generation, manual to randomize
    if nx.has_eulerian_path(reachable_subgraph, start_vertex):
        # Initialize the path
        path = []
        current = start_vertex
        temp_graph = reachable_subgraph.copy()
        total_edges = reachable_subgraph.number_of_edges()
#        vertex_stack = [(start_vertex, None)]
#       last_key = None
        # While there are still outgoing edges from current node
        while total_edges > len(path):
            if temp_graph.out_degree(current) > 0:
                next_el, current = next_element(temp_graph, current)
                path.append(next_el)
                # Get all outgoing edges, can just work further
            else:
                for i in range (0, len(path)):
                    if temp_graph.out_degree(path[i][0]) > 0:
                        add_cycle = []
                        cycle_current = path[i][0]
                        while temp_graph.out_degree(cycle_current) > 0:
                            next_el, cycle_current = next_element(temp_graph, cycle_current)
                            add_cycle.append(next_el)
                        logging.debug (f"adding cycle of length {len(add_cycle)}")
                        logging.debug (f"{add_cycle}")
                        logging.debug(f"after {path[:i]}")
                        path = path[:i] + add_cycle + path[i:]                        

            # Move to next node
        logging.info(f"Randomized Eulerian path found with seed {seed} with {len(path)} nodes !")
        logging.debug (f"{get_gaf_string(path)}")
        return path
    else:
        logging.error(f"Path not found from {start_vertex}")
        for v in reachable_subgraph.nodes():
            logging.info(f"{v}: in {reachable_subgraph.in_degree[v]} out {reachable_subgraph.out_degree[v]}")
        for v in reachable_subgraph.nodes():
            if reachable_subgraph.in_degree[v] != reachable_subgraph.out_degree[v]:
                logging.warning(f"Not Eulerian {v}: in {reachable_subgraph.in_degree[v]} out {reachable_subgraph.out_degree[v]}")
        exit()
        return []

def get_gaf_path(path):
    res = []
    for i in range(len(path)):
        edge = path[i][2]
        #TODO: should be from hashes to original names        
        if edge > 0:            
            node = ">" + node_id_to_name[abs(edge)]
        else:
            node = "<" + node_id_to_name[abs(edge)]
        res.append(node)
    return res

def get_gaf_string(path):
    path_arr = get_gaf_path(path)
    res = "".join(path_arr)
    #logging.info(res)
    return res


class AlignmentScorer:
    #TODO: deprioritize based on lengths and not just node counts
    DEPRIORITIZE_RC_COEFFICIENT = 0.9
    def __init__(self, alignments):
        self.automaton = ahocorasick.Automaton()
        self.pattern_counts = {}
        
        #from lexicographical minimum of (pattern, rc_pattern) to max
        self.rc_patterns = {}
        for idx, alignment in enumerate(alignments):
            logging.debug(f"Adding alignment {idx}: {alignment}")
            pattern_str = self.aln_to_string(alignment)
            rc_nodes = [-n for n in alignment]
            rc_nodes.reverse()
            rc_pattern_str = self.aln_to_string(rc_nodes)
            
            #lexicographical minimum
            #if ",".join(nodes) > ",".join(rc_nodes):
             #   nodes = rc_nodes

            if pattern_str not in self.pattern_counts:
                self.pattern_counts[pattern_str] = 0
                self.pattern_counts[rc_pattern_str] = 0
                self.automaton.add_word(pattern_str, pattern_str)
                self.automaton.add_word(rc_pattern_str, rc_pattern_str)     
                self.rc_patterns[pattern_str] = rc_pattern_str
                self.rc_patterns[rc_pattern_str] = pattern_str   
            self.pattern_counts[pattern_str] += 1
            self.pattern_counts[rc_pattern_str] += 1

        self.automaton.make_automaton()
        logging.debug(f"automaton keys {list(self.automaton.keys())}")
        logging.info(f"Built automaton with {len(self.pattern_counts)} unique alignment patterns")

    def path_to_string(self, path):
        return "," + ",".join(str(edge[2]) for edge in path) + ","
    
    def aln_to_string(self, aln):
        return "," + ",".join(str(node) for node in aln) + ","
    
    def score_corasick(self, path):
        node_ids = [edge[2] for edge in path]
        path_length = len(node_ids)
        path_str = self.path_to_string(path)
        found = set()
        for item in self.automaton.iter(path_str):
            pattern = item[1]
            found.add(pattern)
        score = 0
        for item in found:
            #only looking for one of (pattern, rc_pattern)
            if self.rc_patterns[item] in found:
                score += self.pattern_counts[item] * self.DEPRIORITIZE_RC_COEFFICIENT
                #TODO: possibly add paths length 1 for better deprioritization?
                logging.debug(f"deprioritizing {item} because of rc")
            else:
                score += self.pattern_counts[item] * 2
        return score

def rc_path(path, rc_vertex_map):
    new_path = []
    for edge in path:
        u, v, edge_id = edge
        new_u = rc_vertex_map[u]
        new_v = rc_vertex_map[v]
        new_path.append((new_v, new_u, -edge_id))
    new_path.reverse()
    return new_path

def get_random_change(path, iter, rc_vertex_map):
    """
    Finds two non-overlapping intervals in the Eulerian path that start and end at the same vertices
    and swaps them to create a new path or inverts a random self-rc interval
       
    :param path: List of edges representing the Eulerian path.
    :param iter: Iteration number for random seed
    :return: Modified path with swapped intervals
    """
    
    random.seed(iter)
    path_length = len(path)
    
    # Pre-build index of matching start/end positions for faster lookup
    start_positions = {}
    end_positions = {}
    
    for idx, (u, v, _) in enumerate(path):
        if u not in start_positions:
            start_positions[u] = []
        start_positions[u].append(idx)
        
        if v not in end_positions:
            end_positions[v] = []
        end_positions[v].append(idx)
    
    max_tries = 10000  # Adjust based on path length
    
    for _ in range(max_tries):
        # Select i and j such that i < j
        i = random.randint(0, path_length - 3)
        j = random.randint(i + 1,  path_length - 2) 
        
        start_node = path[i][0]
        end_node = path[j][1]
        #We can do inversion
        if rc_vertex_map[start_node] == end_node and iter % 2 == 0:
            # not allowing to invert AUX node            
            if "AUX" in name_to_node_id:
                forbidden = False
                aux_id = name_to_node_id["AUX"]
                for ind in range(i, j + 1):
                    if path[ind][2] == aux_id:
                        forbidden = True
                        break
                if forbidden:
                    logging.debug(f"Skipping inversion due to AUX node presence between {i} and {j}")
                    continue
            # Invert the interval from i to j
            new_path = path[:i] + rc_path(path[i:j + 1], rc_vertex_map) + path[j + 1:]
            logging.debug(f"{_ + 1} attempts to generate random inversion")
            logging.debug(f"{get_gaf_string(new_path)}")
            return new_path
        
        # Find k where path[k][0] == start_node and k > j
        k_candidates = [k for k in start_positions.get(start_node, []) if k > j]
        if not k_candidates:
            continue
        # Choose k randomly without temperature weighting
        k = random.choice(k_candidates)

        # Find l such that path[l][1] == end_node and l > k
        l_candidates = [l for l in end_positions.get(end_node, []) if l > k]
        if not l_candidates:
            continue
        
        # Choose l randomly without temperature weighting
        l = random.choice(l_candidates)
        # Swap the intervals
        new_path = (
            path[:i]
            + path[k:l + 1]
            + path[j + 1:k]
            + path[i:j + 1]
            + path[l + 1:]
        )
        logging.debug(f"{_ + 1} attempts to generate random swap")
        logging.debug(f"{get_gaf_string(path)}")
        return new_path    
    logging.warning("Failed to find valid intervals to swap")
    logging.warning(f"{get_gaf_string(path)}")
    #exit(9)
    return path

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
                int_node = parse_node_id(fnode[1:])
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

#utig4-234 -> 234
def parse_node_id(node_str): 
    if node_str in name_to_node_id:
        return name_to_node_id[node_str]
    global last_enumerated_node
    parts = node_str.split('-')
    if len(parts) < 2 or not parts[1].isdigit():
        #ribotin graph case
        if node_str.isdigit():            
            node_id = int(node_str)
        else:
            last_enumerated_node += 1
            node_id = last_enumerated_node
            logging.debug(f"Assigned enumerated ID {node_id} to node {node_str}")
    else:
        node_id = int(parts[1])

    #utig4-0 same as its RC
    if node_id < last_enumerated_node:
        last_enumerated_node += 1
        node_id = last_enumerated_node
        logging.info(f"Assigned enumerated ID {node_id} to node {node_str}, utig4-0 special case")
    node_id_to_name[node_id] = node_str
    name_to_node_id[node_str] = node_id
    return node_id

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
        logging.debug(f"Tangle node {abs(u)}")
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
            logging.info(f"Cleaning {n} from tangle")
            tangle_nodes.remove(n)
            directed_graph.remove_node(n)
            


def setup_logging(args):
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    # Exclude date, include time, function name, and line number
    log_format = '%(asctime)s - %(levelname)s - [%(funcName)s:%(lineno)d] - %(message)s'
    datefmt = '%H:%M:%S'
        
    # Always log to both file and console
    log_file = os.path.join(args.outdir, f"{args.basename}.log")

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Clear any existing handlers
    root_logger.handlers = []
    
    # File handler
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setFormatter(logging.Formatter(log_format, datefmt=datefmt))
    file_handler.setLevel(log_level)
    root_logger.addHandler(file_handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(logging.Formatter(log_format, datefmt=datefmt))
    console_handler.setLevel(log_level)
    root_logger.addHandler(console_handler)
    
    # Log the GitHub commit hash if available
    try:
        commit_hash = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=sys.path[0]).strip().decode('utf-8')
        logging.info(f"GitHub Commit Hash: {commit_hash}")
    except Exception as e:
        logging.warning(f"Failed to retrieve GitHub commit hash: {e}")

    # Log the command-line arguments
    logging.info(f"Command-line arguments: {' '.join(sys.argv)}")
    logging.info(f"Logging to file: {log_file}")

def log_assert(condition, message, logger=None):
    """Assert a condition and log an error message if it fails."""
    if not condition:
        error_msg = f"Assertion failed: {message}"
        if logger:
            logger.error(error_msg)
        else:
            logging.error(error_msg)
        exit(1)
        #raise AssertionError(error_msg)

def read_tangle_nodes(args, original_graph):
    """Read or construct tangle nodes."""
    if args.tangle_file:
        tangle_nodes = set()
      #  nor_nodes = set()
        with open(args.tangle_file, 'r') as f:
            for line in f:
                node_str = line.strip()
                if node_str:
                    node_id = parse_node_id(node_str)
                    tangle_nodes.add(node_id)
                    tangle_nodes.add(-node_id)
                  #  nor_nodes.add(node_id)
    elif args.tangle_node:
        tangle_nodes = node_to_tangle(original_graph, args.tangle_length_cutoff, target_node=parse_node_id(args.tangle_node))
    else:
        logging.error("Either --tangle_file or --tangle_node must be provided.")
        sys.exit(1)

    with open(os.path.join(args.outdir, "nodes.txt"), 'w') as f:
        for node in tangle_nodes:
            #only forward nodes to file
            if node > 0:
                f.write(f"{node_id_to_name[node]}\n")

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
                node_id = parse_node_id(arr[0])
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

def calculate_median_coverage(args, nor_nodes, original_graph, cov, boundary_nodes):
    """Calculate or use provided median unique coverage."""
    min_b = 1000000
    all_boundary = list(boundary_nodes.keys()) + list(boundary_nodes.values())
    for b in all_boundary:
        logging.info (f"{b} {cov[abs(b)]}")
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
                out_file.write(f"{node_id_to_name[node_id]}\t{cov_value}\t{mult_value}\n")
    logging.info(f"Wrote multiplicity solutions to {output_file}")

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
        path = get_traversing_eulerian_path(multi_graph, boundary_nodes, original_graph, seed)
        if not path:
            logging.warning(f"No Eulerian path found for seed {seed}.")
            continue
        
        current_path = path
        current_score = alignment_scorer.score_corasick(current_path)
        logging.info(f"Initial path score for seed {seed}: {current_score}.")

        iterations_since_improvement = 0
        for i in range(max_iterations):
            if iterations_since_improvement >= early_stopping_limit:
                logging.info(f"Early stopping for seed {seed} on iteration {i} after {iterations_since_improvement} iterations without improvement.")
                break
            
            new_path = get_random_change(current_path, seed * max_iterations + i, rc_vertex_map)
            new_score = alignment_scorer.score_corasick(new_path)
            if new_score > current_score:
                logging.info(f"Improved score for seed {seed} at iteration {i}: {current_score} -> {new_score}.")
                current_path = new_path
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
    return best_path, best_score

#Messy stuff excluded from main
def generate_MIP_equations(tangle_nodes, nor_nodes, cov, median_unique, original_graph, boundary_nodes, directed = False):
    junction_equations = []    
    used = set()
    extended_tangle = tangle_nodes.copy()
    #format: used ins to corresponding used outs
    all_boundary_nodes = set()
    for b in boundary_nodes:
        all_boundary_nodes.add(abs(b))
        all_boundary_nodes.add(-abs(b))
        all_boundary_nodes.add(abs(boundary_nodes[b]))
        all_boundary_nodes.add(-abs(boundary_nodes[b]))
    extended_tangle.update(all_boundary_nodes)
    canonic_name = {}
    for node in extended_tangle:
        if directed:
            canonic_name[node] = node
        else:
            canonic_name[node] = abs(node)
    # Generate equations
    for from_node in extended_tangle:
        if from_node in used:
            continue
        arr = [[],[]]
        back_node = ""
        bad_extension = 0
        for to_node in original_graph.successors(from_node):
            back_node = to_node
            if not directed:
                used.add(-to_node)
            if to_node in extended_tangle:
                arr[1].append(canonic_name[to_node])
            else:
                bad_extension = to_node
                break
           
        if  bad_extension != 0:
            if not (from_node in all_boundary_nodes):
                logging.error(f"Somehow jumped over boundary nodes {from_node} {back_node} { boundary_nodes}")
                exit()
            else:
                continue
        if back_node != "":
            for alt_start_node in original_graph.predecessors(back_node):
                used.add(alt_start_node)
                if alt_start_node in extended_tangle:
                    arr[0].append(canonic_name[alt_start_node])
                else:
                    bad_extension = alt_start_node
                    break
        if bad_extension != 0:
            logging.error(f"Somehow jumped over boundary nodes (backwards) {alt_start_node} {bad_extension} { boundary_nodes}")
            exit()
        junction_equations.append(arr)

    must_use_nodes = []
    coverage = {}
    for node in nor_nodes:        
        coverage[node] = float(cov[node]) #/ median_unique
        logging.debug(f"Coverage of {node} : {coverage[node]}")
        if coverage[node] / median_unique >= 0.5:
            must_use_nodes.append(node)
    boundary_values = {}
    for b in boundary_nodes:
        boundary_values[b] = 1
        boundary_values[boundary_nodes[b]] = 1
        boundary_values[-b] = 0
        boundary_values[-boundary_nodes[b]] = 0
        coverage[abs(b)] = float(cov[abs(b)]) 
        coverage[abs(boundary_nodes[b])] = float(cov[abs(boundary_nodes[b])])
    return junction_equations, must_use_nodes, coverage, boundary_values


def identify_boundary_nodes(args, original_graph, tangle_nodes):
    if args.boundary_nodes:
        boundary_nodes = {}
        for line in open (args.boundary_nodes):
            # Parse the line and extract boundary node pairs
            #map incoming->matching outgoing
            parts = line.strip().split()
            if len(parts) == 2:
                node1 = parse_node_id(parts[0])
                is_incoming = False
                for next in original_graph.successors(node1):
                    if next in tangle_nodes:
                        is_incoming = True
                        break
                if not is_incoming:
                    node1 = -node1
                node2 = parse_node_id(parts[1])
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
        log_assert(len(out_boundary_nodes) == 2 and len(in_boundary_nodes) == 2, f"Autodetection works only for 1-1 tangles, detected out_boundary: {out_boundary_nodes}, in_boundary: {in_boundary_nodes}. Specify boundary node pairs manually")
        res = {}
        start = max(in_boundary_nodes)
        end = max(out_boundary_nodes)
        log_assert(abs(start) != abs(end), f"Start and end boundary nodes should be different, got {start} and {end}")
        res[start] = end
        return res

#Only UNIQUE_BORDER_LENGTH (=200K) suffix/prefix for border unique nodes used
def output_path_fasta(best_path, original_graph, output_file):
    with open(output_file, 'w') as fasta_file:
        contig_sequence = ""
        last_node = None
        for i in range (len(best_path)):
            edge_id = best_path[i][2]
            node_id = abs(edge_id)
            orientation = edge_id > 0

            # Retrieve the sequence from the graph
            node_sequence = original_graph.nodes[node_id]['sequence']
            if node_sequence == "*":
                logging.error("Provided noseq assembly graph, no fasta output possible")
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
        fasta_file.write(">contig\n")
        fasta_file.write(f"{contig_sequence}\n")
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
    if args.verkko_output is None and (args.gfa_file is None or args.alignment is None):
        logging.error("Either --verkko-output OR both --gfa and --alignment are required options")
        return
    if args.verkko_output: 
        args.coverage_file = os.path.join(args.verkko_output, "2-processGraph", "unitig-unrolled-hifi-resolved.ont-coverage.csv")
        args.gfa = os.path.join(args.verkko_output, "2-processGraph", "unitig-unrolled-hifi-resolved.gfa")
        args.alignment = os.path.join(args.verkko_output, "3-align", "alns-ont.gaf")
    return args

def main():    
    args = parse_arguments()
    
    os.makedirs(args.outdir, exist_ok=True)

    setup_logging(args)
    logging.info("Reading files...")

    # Create output directory if it doesn't exist

    #TODO: save it somewhere, some connections added while parsing    
    original_graph = parse_gfa(args.gfa_file)
    tangle_nodes = read_tangle_nodes(args, original_graph)    
    clean_tips(tangle_nodes, original_graph)
    nor_nodes = {abs(node) for node in tangle_nodes}
    
    boundary_nodes = identify_boundary_nodes(args, original_graph, tangle_nodes)

    used_nodes = nor_nodes.copy()
    for b in boundary_nodes:
        used_nodes.add(abs(b))
        used_nodes.add(abs(boundary_nodes[b]))
   
    if args.coverage_file:
        cov = read_coverage_file(args.coverage_file)
    else:
        cov = coverage_from_graph(original_graph)
    #TODO: some edges can be missing, fill them with ??? (longest node coverage?)
    median_unique_range = calculate_median_coverage(args, nor_nodes, original_graph, cov, boundary_nodes)    
    median_unique = math.sqrt(median_unique_range[0] * median_unique_range[1])
    filtered_alignment_file = os.path.join(args.outdir, f"{args.basename}.q{args.quality_threshold}.used_alignments.gaf")
    alignments = parse_gaf(args.alignment, used_nodes, filtered_alignment_file, args.quality_threshold)
    alignment_scorer = AlignmentScorer(alignments)
    #Shit is hidden here
    logging.info("Starting multiplicity counting...")
    #TODO: instead of a_values we just use coverage
    equations, nonzeros, a_values, boundary_values = generate_MIP_equations(tangle_nodes, nor_nodes, cov, median_unique, original_graph, boundary_nodes, directed=True)
    #Failed to generate suboptimal solutions yet
    best_solution = solve_MIP(equations, nonzeros, boundary_values, a_values, median_unique_range)
    
    # Define output filenames based on the output directory
    output_csv = os.path.join(args.outdir, args.basename + ".multiplicities.csv")
    output_fasta = os.path.join(args.outdir, args.basename + ".fasta")
    output_gaf = os.path.join(args.outdir, args.basename + ".gaf")

    # Write multiplicities to CSV
    write_multiplicities(output_csv, best_solution, cov)

    
    #Now all sequence is stored in edges, junctions are new vertices
    dual_graph = create_dual_graph(original_graph)
    multi_graph = create_multi_dual_graph(dual_graph, best_solution, tangle_nodes, boundary_nodes, original_graph)
    
    best_path, best_score = optimize_paths(multi_graph, boundary_nodes, original_graph, args.num_initial_paths, args.max_iterations, args.early_stopping_limit, alignment_scorer)
    logging.info("Path optimizing finished.")
    if best_path:
        best_path_str = get_gaf_string(best_path)
        logging.info(f"Found traversal\t{best_path_str}")
        
        # Output FASTA file
        logging.info(f"Writing best path to {output_fasta}")
        output_path_fasta(best_path, original_graph, output_fasta)
        
        # Output GAF file
        logging.info(f"Writing best path to {output_gaf}")
        with open(output_gaf, 'w') as outf:
            outf.write(f"tangle\t{best_path_str}\n")
    

if __name__ == "__main__":
    main()
