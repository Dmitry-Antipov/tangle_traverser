#!/usr/bin/env python3

import logging
import random
import networkx as nx
from .logging_utils import log_assert
from .node_id_mapper import NodeIdMapper

def rc_node(node):
    """RC nodes stored as negative"""
    return -node

def get_canonical_nodepair(oriented_node1, oriented_node2, original_graph, node_mapper):
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

def get_canonical_rc_vertex(v, original_graph, node_mapper):
    """Returns a canonical reverse connection vertex for a given vertex"""
    return get_canonical_nodepair(-v[1], -v[0], original_graph, node_mapper)

def create_dual_graph(original_graph: nx.MultiDiGraph, node_mapper):
    """
    Transform to dual graph, vertices = junctions, edges = old nodes
    TODO: vertices should be not exactly junctions but junction+nodes since z-connections are legit
    """
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
            C1 = get_canonical_nodepair(u, v, original_graph, node_mapper)
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
            C2 = get_canonical_nodepair(v, w, original_graph, node_mapper)
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

def is_tangle_vertex(C, tangle_set):
    return C[0] in tangle_set or C[1] in tangle_set

def create_multi_dual_graph(dual_graph: nx.DiGraph, multiplicities: dict, tangle_nodes: set, boundary_nodes: dict, G, node_mapper):
    """
    edge of multiplicity X -> X multiedges multiplicity 1
    """
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
            fw_edge = get_canonical_nodepair(b, next_node, G, node_mapper)                                
            logging.debug(f"start {node_mapper.node_id_to_name_safe(b)} {node_mapper.node_id_to_name_safe(next_node)}")
            start = (0, b)
            multi_dual_graph.add_node(start)
            multi_dual_graph.add_edge(start, fw_edge, original_node=b, key = f"{b}_{edges_added}")
            edges_added += 1
            break 
    for b in boundary_nodes.values():
        for next_node in G.predecessors(b):
            bw_edge = get_canonical_nodepair(next_node, b, G, node_mapper)                                
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

def get_traversable_subgraph(multi_dual_graph: nx.MultiDiGraph, border_nodes, original_graph, node_mapper):
    """
    Supplementary for Euler path search - search itself moved to path_optimizer.py
    """
    start_vertices = []
    end_vertices = []
    for n in border_nodes:
        #border tips encoding
        start_vertices.append((0, n))
        end_vertices.append((n, 0))

    logging.info(f"Start and end vertices in the graph {start_vertices}, {end_vertices}")      

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
            multi_dual_graph.add_edge(get_canonical_rc_vertex(e[1], original_graph, node_mapper), get_canonical_rc_vertex(e[0], original_graph, node_mapper), original_node=-int(data['original_node']), key = new_key)

    logging.debug(f"After transformation")
    for e in multi_dual_graph.edges(keys=True):
        logging.debug(f"Edge {e}")
    log_assert(len(unreachable_edges) == 0, f"Unreachable edges are still present: {unreachable_edges}")
    reachable_subgraph = multi_dual_graph.subgraph(reachable_verts)
    return reachable_subgraph, start_vertex
