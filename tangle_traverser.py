#!/usr/bin/env python
import sys
import pulp
import argparse
import logging
import statistics
import networkx as nx
import random
import cov_support
import cProfile
import pstats
import time

# Add a global map for node ID to string node names
#only positive ids (unoriented)
node_id_to_name = {}

def solve_integer_system(equations, nonzeros, boundary_nodes, a_values, num_solutions=3):
    # Create a linear programming problem
    prob = pulp.LpProblem("Minimize_Deviation", pulp.LpMinimize)
    
    # Define integer variables
    x_vars = {i: pulp.LpVariable(f"x_{i}", cat="Integer") for i in a_values}
    
    # Define continuous variables for absolute deviation |x_i - a_i|
    d_vars = {i: pulp.LpVariable(f"d_{i}", lowBound=0, cat="Continuous") for i in a_values}
    
    # Objective function: minimize sum of absolute deviations
    objective = pulp.lpSum(d_vars[i] for i in a_values)
    prob += objective
    
    # Constraints: x_i = x_j + x_k + x_l
    #for lhs, rhs in equations.items():
     #   prob += x_vars[lhs] == sum(x_vars[j] for j in rhs)
    for eq in equations:
        boundary = 0
        real_nodes = []
        for j in eq[1:]:
            if j in boundary_nodes:
                boundary += 1
            else:
                real_nodes.append(j)
        prob += x_vars[eq[0]] == sum (x_vars[j] for j in real_nodes) + boundary
    for x_var in x_vars:
        if x_var in nonzeros:
            prob += x_vars[x_var] >= 1
        else:
            prob += x_vars[x_var] >= 0
    # Constraints for absolute deviation linearization: d_i >= |x_i - a_i|
    for i in a_values:
        prob += d_vars[i] >= x_vars[i] - a_values[i]
        prob += d_vars[i] >= a_values[i] - x_vars[i]
    
    solutions = []
    iter = 0
    all_inds = {}
    # python
    while len(solutions) < num_solutions:
        # Solve the problem
        prob.solve()
        
        # Extract results
        result = {i: pulp.value(x_vars[i]) for i in a_values}
        score = pulp.value(objective)  
        for ind in all_inds.keys():
            if pulp.value(ind) > 0:
                logging.info (f"nonzero {ind} {pulp.value(ind)} {all_inds[ind]}")
                
                break
        logging.info (f"iteration of MIP {iter}")
        # Check if the solver found a feasible solution
        if pulp.LpStatus[prob.status] != "Optimal" or result in solutions or None in result.values():
            break  # No more unique solutions found or no feasible solution
        solutions.append((result, score))
        #all suboptimal solution hacks fail.
        break
        diff_indicators = {i: pulp.LpVariable(f"diff_{iter}_{i}", cat="Binary") for i in a_values}
        for ind in diff_indicators:
            all_inds[diff_indicators[ind]] = ind
        M = 1000
        for i in a_values:
            # If x_vars[i] equals result[i], diff_indicators[i] can be 0 or 1
            # If x_vars[i] not equals result[i], diff_indicators[i] can be only 1
            prob += x_vars[i] - result[i] <= diff_indators[i] * M
            prob += result[i] - x_vars[i] <= diff_indicators[i] * M

            # If x_vars[i] equals result[i], diff_indicators[i] can be  only 0
            # If x_vars[i] not equals result[i], diff_indicators[i] can be  0  or 1
            prob += x_vars[i] - result[i] <  M - diff_indicators[i] * M
            prob += result[i] - x_vars[i] <  M - diff_indicators[i] * M
            print (f"{i} {x_vars[i]} {result[i]} {diff_indicators[i]}")
            #exit(0)
        # Ensure that not all variables are the same
        prob += pulp.lpSum(diff_indicators[i] for i in a_values) >= 1   # At least one variable must differ

        iter += 1   
    return solutions
       
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

            # Process only link lines (L)
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
                                logging.warning(f"Non-transitive junction! Adding {start_node} -> {n2}, overlap {min_overlap}")
                                overlap_counted = True
                                break
                        if not overlap_counted:
                            logging.error(f"Non-transitive junction failed to get overlap. Adding {start_node} -> {n2}, overlap 0")
                            original_graph.add_edge(start_node, n2, overlap=0)

    return original_graph

#RC nodes stored as negative
def rc_node(node):
    return -node

#Possibly we need to remove dead-ends and very low covered nodes

#Detect unique nodes in tangle and count median
def is_forward_unique(oriented_node, original_graph):
    """Checks if an oriented node has exactly one outgoing edge."""
    return len(list(original_graph.successors(oriented_node))) == 1

def calculate_median_unique_coverage(nor_nodes, original_graph, cov):
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
            if node_id in cov:
                coverage = float(cov[node_id])
                unique_coverages.append(coverage)
                unique_node_ids.add(node_id)
                logging.debug(f"Node {node_id} is unique. Coverage: {coverage}")
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
   
def get_canonical_edge(oriented_node1, oriented_node2, original_graph):
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
    
#Transform to dual graph, vertices = junctions, edges = old nodes
def create_dual_graph(original_graph:nx.MultiDiGraph):
    dual_graph = nx.MultiDiGraph()
    canonical_edges_set = set() # To store unique canonical edges (dual nodes)
    logging.info("Creating dual graph representation (nodes = canonical connections) using NetworkX...")

    # Iterate through each oriented node 'v' which acts as the junction
    for v in original_graph.nodes:
        # Determine the canonical edge C1 leading INTO node v
        predecessors = list(original_graph.successors(rc_node(v)))
        if not predecessors:
            # v is a source tip (no incoming edges for this orientation)
            C1 = ("TIP", v)
        else:
            # Pick any predecessor 'u' to define the canonical incoming edge
            # Note: get_canonical_edge expects (from, to), so use (pred, v)
            u = rc_node(predecessors[0])
            C1 = get_canonical_edge(u, v, original_graph)
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
            C2 = get_canonical_edge(v, w, original_graph)
        canonical_edges_set.add(C2)

        # Add the edge in the dual graph representing the transition through v
        # Ensure nodes C1 and C2 exist before adding the edge (add_nodes_from handles this later too)
        dual_graph.add_node(C1)
        dual_graph.add_node(C2)
        dual_graph.add_edge(C1, C2, original_node=v)

    logging.info(f"Total oriented nodes in original GFA structure: {len(original_graph.nodes)}")
    logging.info(f"Total unique canonical edges (dual nodes): {len(canonical_edges_set)}")
    logging.info(f"Total nodes in final dual graph: {dual_graph.number_of_nodes()}") # Should match set size
    logging.info(f"Total edges (connections through original nodes) in dual graph: {dual_graph.number_of_edges()}")

    return dual_graph

#as for now, multiplities are for unoriented nodes
#Transform dual graph to multi-dual graph based on multiplicities
#TODO: not supporting hairpins and unoriented multiplicities

def is_tangle_vertex (C, tangle_set):
    return C[0] in tangle_set or C[1] in tangle_set

def create_multi_dual_graph(dual_graph: nx.DiGraph, multiplicities: dict, tangle_nodes: set, boundary_nodes: set, G):
    multi_dual_graph = nx.MultiDiGraph()
    logging.info("Creating multi-dual graph from dual graph and multiplicities...")

    # Ensure all nodes from the original dual graph are present
    for node in dual_graph.nodes():
        if is_tangle_vertex(node, tangle_nodes):
            multi_dual_graph.add_node(node)

#Dirty hacks to add start and sink
    for b in boundary_nodes:
        for orientation in [1, -1]:
            id = b * orientation
            for s in G.successors(id):
                if s in tangle_nodes:
                    fw_edge = get_canonical_edge(id, s, G)                                
                    logging.debug(f"start {id} {s}")
                    start = (0, id)
                    multi_dual_graph.add_node(start)
                    multi_dual_graph.add_edge(start, fw_edge, original_node=id, key = str(id) + "_"+ "0")
                    break
            
            for s in G.predecessors(id):
                if s in tangle_nodes:         
                    bw_edge = get_canonical_edge(s, id, G)                
                    logging.debug(f"end {s} {id}")
                    end = (id, 0)
                    multi_dual_graph.add_node(end)
                    multi_dual_graph.add_edge(bw_edge, end, original_node=id, key = str(id) + "_"+ "0")
                    break

    
    
    # Iterate through edges of the dual graph and add them to the multi-graph
    # based on the multiplicity of the original node they represent.
    edges_added = 0
    for u, v, data in dual_graph.edges(data=True):
        if (not u in multi_dual_graph.nodes()) or (not v in multi_dual_graph.nodes()):
            continue
        original_node_oriented = data.get('original_node')
        if not original_node_oriented:
            logging.warning(f"Edge ({u}, {v}) in dual graph is missing 'original_node' attribute. Skipping.")
            continue
        
        original_node_base = abs(original_node_oriented)

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
def get_traversing_eulerian_path(multi_dual_graph: nx.MultiDiGraph, border_nodes, original_graph, seed):
    random.seed(seed)
    #TODO: start from ("Tip, +-border_node"), would simplify code
    start_vertices = []

    end_vertices = []
    for n in multi_dual_graph.nodes():
        if n[0] == 0:
            start_vertices.append(n)
        elif n[1] == 0:
            end_vertices.append(n)
    logging.info(f"Start and end vertices in the graph {start_vertices}, {end_vertices}")

    # Only 1-1 tangles for now
    if len(start_vertices) != 2 or len(end_vertices) != 2:
        logging.error(f"border vertices fail") 
        return []
    
    #TODO: possibly add option to always start from specified border node
    if abs(start_vertices[0][1]) < abs(start_vertices[1][1]):
        start_vertex = start_vertices[0]
    else:
        start_vertex = start_vertices[1]    
    logging.info(f"Using start vertex {start_vertex} for seed {seed}")
    reachable_verts = nx.descendants(multi_dual_graph, start_vertex)
    # Add the start_node itself
    reachable_verts.add(start_vertex)

    # Create the subgraph containing only the reachable nodes and edges between them
    reachable_subgraph = multi_dual_graph.subgraph(reachable_verts)
    reachable_end_vertices = []
    for v in reachable_subgraph.nodes():
        if v in end_vertices:
            reachable_end_vertices.append(v)
    if len(reachable_end_vertices) != 1:
        logging.error(f"Wrong amount of reachable end vertices {start_vertex} {reachable_end_vertices} all end vertices {end_vertices}")
        return []
    end_vertex = reachable_end_vertices[0]

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
        return path
    else:
        logging.info(f"Path not found from {start_vertex}")
        for v in reachable_subgraph.nodes():
            logging.info(f"{v}: in {reachable_subgraph.in_degree[v]} out {reachable_subgraph.out_degree[v]}")
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

def score_compressed_path(alignments, path):
    res = 0    
    for add in range (2, len(alignments)):
        possible_next = {}
        for i in range(len(path) - add):
            if not path[i][2] in possible_next:
                possible_next[path[i][2]] = set()
            possible_next[path[i][2]].add(path[i+add][2])
        for start in alignments[add]:
            for end in alignments[add][start]:
                if start in possible_next and end in possible_next[start]:
                    res += alignments[add][start][end]
        logging.debug(f"Score compressed after processing distance {add}: {res}")                
    return res

def get_random_change(path, iter, temperature=1.0):
    """
    Finds two non-overlapping intervals in the Eulerian path that start and end at the same vertices
    and swaps them to create a new path.
       
    :param path: List of edges representing the Eulerian path.
    :param iter: Iteration number for random seed
    :param temperature: Not used, need to understand local/global changes better
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

#TODO: should be hashes for general graph
def parse_gaf(gaf_file, interesting_nodes, filtered_file, quality_threshold=0):
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
            for char in parts[5]:
                if char in "<>":
                    if node:
                        nodes.append(node)
                    node = char
                else:
                    node += char
            if node:
                nodes.append(node)
            nodes = [node.split('_')[0] for node in nodes if node]
            filtered_nodes = []
            for i, node in enumerate(nodes):
                if i == 0 or node != nodes[i - 1]:
                    filtered_nodes.append(node)
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
            if good and len(nodes) > 2:
                res.append(nodes)
                if filtered_file:
                    out_file.write(line)
            #Not appending reverce-complement, this is processed in scoring!
    return  res 

#comppressed storing for graphaligner alignments.
def add_pair(add, start_a, end_a, arr):
    if not (start_a in arr[add].keys()):
        arr[add][start_a] = {}
    if not (end_a in arr[add][start_a]):
        arr[add][start_a][end_a] = 0
    arr[add][start_a][end_a] += 1

#align[i]: dict of node->node->multiplicity
def reform_alignments(alignments):   
    #TODO: do we need this constraint?
    MAX_NODES = 50
    perfect_score = 0
    different_pairs = 0
    longest = 0
    for al in alignments:
        longest = max(longest, len(al))
    
    res = []
    max_range = min(MAX_NODES, longest)
    for i in range (longest):
        res.append({})
    for al in alignments:
        lenal = len(al)        
        for i in range (lenal):
            for add in range (2, min(max_range, lenal - i)):  
                add_pair(add, al[i], al[i + add],res)                
                add_pair (add, -al[i+add], -al[i],res)
                perfect_score += 2
    for i in range(len(res)):
        for st in res[i]:
            different_pairs += len(res[i][st])
    #Not correct if hairpins
    logging.info (f"longest alignment {max_range} best possible score is {perfect_score} / 2, different scored pairs {different_pairs}")            
    return res

#utig4-234 -> 234
#TODO: hashes for non-verkko cases
def parse_node_id(node_str):    
    
    parts = node_str.split('-')
    if len(parts) < 2 or not parts[1].isdigit():
        #ribotin graph case
        if node_str.isdigit():            
            node_id = int(node_str)
        else:
            raise ValueError(f"Invalid node string: {node_str}")
    else:
        node_id = int(parts[1])

    #utig4-0 same as its RC
    if node_id == 0:
        node_id = 123123123
    node_id_to_name[node_id] = node_str
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
    rc_component = nx.node_connected_component(indirect_graph, -target_node)
    connected_component.update(rc_component)
    return connected_component

def setup_logging(args):
    """Setup logging configuration."""
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    if args.log_file:
        logging.basicConfig(filename=args.log_file, level=log_level, format=log_format, filemode='w')
    else:
        logging.basicConfig(level=log_level, format=log_format, stream=sys.stderr)

# Update the parse_arguments function to include new options
def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Solve for integer multiplicities in a GFA tangle graph based on coverage.")
    parser.add_argument("gfa_file", help="Path to the GFA file.")
    parser.add_argument("alignment", help="Path to a file with graphaligner alignment")
    parser.add_argument("--coverage-file", help="Path to a file with node coverages (node-id coverage). If not provided, coverage will be filled from the GFA file.")
    parser.add_argument("--median-unique", type=float, help="Median coverage for unique nodes.")
    parser.add_argument("--output-multiplicities", default="mip-mults.csv", help="Path to the output CSV file (default: mip-mults.csv).")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Set the logging level (default: INFO).")
    parser.add_argument("--log-file", help="Optional file to write logs to (default: stderr).")
    parser.add_argument("--tangle-file", help="Path to a file listing nodes in the tangle (one per line).")
    parser.add_argument("--tangle-node", type=str, help="Node ID to construct tangle as connected component of short edges around it")
    parser.add_argument("--tangle-length-cutoff", type=int, default=500000, help="Length cutoff for tangle detection, default 500K")
    parser.add_argument("--num-initial-paths", type=int, default=10, help="Number of initial paths to generate (default: 10).")
    parser.add_argument("--max-iterations", type=int, default=100000, help="Maximum iterations for path optimization (default: 100000).")
    parser.add_argument("--early-stopping-limit", type=int, default=15000, help="Early stopping limit for optimization (default: 15000).")
    parser.add_argument("--filtered-alignment-file", type=str, help="Filtered alignments to the tangle can be saved in this file for later reuse")
    parser.add_argument("--quality-threshold", type=int, default=0, help="Alignments with quality less than this will be filtered out, default 0 (no quality filtering)")
    parser.add_argument("--output-fasta", type=str, help="Path to output the best path as a contig in FASTA format.")
    parser.add_argument("--output-gaf", default="tangle.gaf", type=str, help="Path to output the best path in GAF format.")
    return parser.parse_args()

def read_tangle_nodes(args, original_graph):
    """Read or construct tangle nodes."""
    if args.tangle_file:
        tangle_nodes = set()
        nor_nodes = set()
        with open(args.tangle_file, 'r') as f:
            for line in f:
                node_str = line.strip()
                if node_str:
                    node_id = parse_node_id(node_str)
                    tangle_nodes.add(node_id)
                    tangle_nodes.add(-node_id)
                    nor_nodes.add(node_id)
    elif args.tangle_node:
        tangle_nodes = node_to_tangle(original_graph, args.tangle_length_cutoff, target_node=parse_node_id(args.tangle_node))
        nor_nodes = {abs(node) for node in tangle_nodes}
    else:
        logging.error("Either --tangle_file or --tangle_node must be provided.")
        sys.exit(1)
    return tangle_nodes, nor_nodes

def read_coverage_file(coverage_file):
    """Read coverage data from file."""
    logging.info(f"Reading coverage from {coverage_file}")
    cov = {}
    with open(coverage_file, 'r') as f:
        for line in f:
            arr = line.strip().split()
            #Usual cov csv format starts with node*
            if len(arr) == 2 and arr[0] != "node":
                node_id = parse_node_id(arr[0])
                cov[node_id] = arr[1]
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
def calculate_median_coverage(args, nor_nodes, original_graph, cov):
    """Calculate or use provided median unique coverage."""
    if args.median_unique is not None:
        return args.median_unique
    calculated_median = calculate_median_unique_coverage(nor_nodes, original_graph, cov)
    if calculated_median is not None:
        return calculated_median
    logging.error("Failed to calculate median unique coverage. Provide it manually with --median-unique.")
    sys.exit(1)

def write_solutions_to_file(output_file, solutions, cov):
    """Write solver solutions to the output file."""
    try:
        with open(output_file, 'w') as out_file:
            out_file.write("node\tcoverage\tmult\n")
            if solutions:
                for node in cov.keys():
                    mult_value = solutions[0].get(node, "X")
                    cov_value = cov.get(node, "N/A")
                    out_file.write(f"{node}\t{cov_value}\t{mult_value}\n")
    except IOError as e:
        logging.error(f"Failed to write output file {output_file}: {e}")
        sys.exit(1)

# Update the optimize_paths function to use command-line options for default parameters
def optimize_paths(multi_graph, boundary_nodes, original_graph, alignments, num_initial_paths, max_iterations, early_stopping_limit):
    """Optimize Eulerian paths."""
    best_path = None
    best_score = -1
    reformed_alignments = reform_alignments(alignments)
    logging.info(f"Starting optimization with {num_initial_paths} initial paths, max {max_iterations} iterations per path.")
    
    for seed in range(num_initial_paths):
        logging.info(f"Generating initial path with seed {seed}.")
        path = get_traversing_eulerian_path(multi_graph, boundary_nodes, original_graph, seed)
        if not path:
            logging.warning(f"No Eulerian path found for seed {seed}.")
            continue
        
        current_path = path
        current_score = score_compressed_path(reformed_alignments, path)
        logging.info(f"Initial path score for seed {seed}: {current_score}.")
        
        iterations_since_improvement = 0
        for i in range(max_iterations):
            if iterations_since_improvement >= early_stopping_limit:
                logging.info(f"Early stopping for seed {seed} on iteration {i} after {iterations_since_improvement} iterations without improvement.")
                break
            
            new_path = get_random_change(current_path, seed * max_iterations + i)
            new_score = score_compressed_path(reformed_alignments, new_path)
            
            if new_score > current_score:
                logging.debug(f"Improved score for seed {seed} at iteration {i}: {current_score} -> {new_score}.")
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
def generate_equations_and_values(tangle_nodes, nor_nodes, cov, median_unique, original_graph, boundary_nodes):
    """
    Generate equations, nonzeros, and a_values (observed multiplicities) for solving the integer system.
    """
    equations = []    

    # Generate equations
    for from_node in tangle_nodes:
        if len(list(original_graph.successors(from_node))) != 1:
            correct = True
            for to_node in original_graph.successors(from_node):
                rc_to_node = rc_node(to_node)
                if len(list(original_graph.successors(rc_to_node))) != 1:
                    correct = False
                    break
            if correct:
                arr = [abs(from_node)]
                for node in original_graph.successors(from_node):
                    if abs(node) in nor_nodes or abs(node) in boundary_nodes:
                        arr.append(abs(node))
                equations.append(arr)

    # Prepare nonzeros = nodes that we always want to include and a_values = "observed multiplicities"
    nonzeros = []
    a_values = {}
    for node in nor_nodes:
        a_values[node] = float(cov[node]) / median_unique
        if a_values[node] >= 0.5:
            nonzeros.append(node)

    return equations, nonzeros, a_values

def identify_boundary_nodes(original_graph, tangle_nodes):
    boundary_nodes = set()
    for first in original_graph.nodes:
        for second in original_graph.successors(first):
            if first in tangle_nodes and second not in tangle_nodes:
                boundary_nodes.add(abs(second))
            elif second in tangle_nodes and first not in tangle_nodes:
                boundary_nodes.add(abs(first))
    return boundary_nodes

#not including border nodes!
def output_path_fasta(best_path, original_graph, output_file):
    """
    Write the best path as a contig in FASTA format.

    :param best_path: List of tuples representing the best path [(u, v, edge_id), ...].
    :param original_graph: The original graph containing node sequences.
    :param output_file: Path to the output FASTA file.
    """
    with open(output_file, 'w') as fasta_file:
        contig_sequence = ""
        last_node = None
        for i in range (len(best_path)):
            if i == 0 or i == len(best_path) - 1:
                continue
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

            if i == 1:
                overlap = 0                
            else:                
                overlap = original_graph.get_edge_data(last_node, edge_id)['overlap']            
            contig_sequence += node_sequence[overlap:]
            last_node = edge_id

        # Write the contig to the FASTA file
        fasta_file.write(">contig\n")
        fasta_file.write(f"{contig_sequence}\n")
    logging.info("Fasta output finished")

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return sequence.translate(complement)[::-1]

def main():
    args = parse_arguments()
    setup_logging(args)

    logging.info("Reading files...")

    #TODO: save it somewhere, some connections added while parsing    
    original_graph = parse_gfa(args.gfa_file)
    tangle_nodes, nor_nodes = read_tangle_nodes(args, original_graph)    
    boundary_nodes = identify_boundary_nodes(original_graph, tangle_nodes)

    used_nodes = nor_nodes.copy()
    used_nodes.update(boundary_nodes)
    if args.coverage_file:
        cov = read_coverage_file(args.coverage_file)
    else:
        cov = coverage_from_graph(original_graph)
    median_unique = calculate_median_coverage(args, nor_nodes, original_graph, cov)    
    alignments = parse_gaf(args.alignment, used_nodes, args.filtered_alignment_file, args.quality_threshold)
    
    #Shit is hidden here
    logging.info("Starting multiplicity counting...")
    equations, nonzeros, a_values = generate_equations_and_values(tangle_nodes, nor_nodes, cov, median_unique, original_graph, boundary_nodes)
    #Failed to generate suboptimal solutions yet
    solutions = solve_integer_system(equations, nonzeros, boundary_nodes, a_values, num_solutions=10)
    
    # Extract the first solution (dictionary) from the tuple
    best_solution = solutions[0][0] if solutions else {}
    #Some implementations of pulp would give you 2.0 instead of 2
    multiplicities = {node: int(round(value)) for node, value in best_solution.items()}
    logging.info(f"Found multiplicities from MIP {multiplicities}")
    
    #Now all sequence is stored in edges, junctions are new vertices
    dual_graph = create_dual_graph(original_graph)
    multi_graph = create_multi_dual_graph(dual_graph, multiplicities, tangle_nodes, boundary_nodes, original_graph)
    
    best_path, best_score = optimize_paths(multi_graph, boundary_nodes, original_graph, alignments, args.num_initial_paths, args.max_iterations, args.early_stopping_limit)
    if best_path:
        best_path_str = get_gaf_string(best_path)
        logging.info(f"Found traversal\t{best_path_str}")
        if args.output_fasta:
            output_path_fasta(best_path, original_graph, args.output_fasta)
        if args.output_gaf:
            outf = open (args.output_gaf, 'w')
            outf.write(f"tangle\t{best_path_str}\n")
    
    #TODO: output to fasta
    logging.info("Path optimizing finished.")
    #TODO: id 0 for start/end replace?
if __name__ == "__main__":
    main()
