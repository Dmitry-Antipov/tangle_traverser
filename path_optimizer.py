#!/usr/bin/env python3

import logging
import random
import dataclasses
import networkx as nx


@dataclasses.dataclass
class EdgeDescription:
    source: int
    target: int
    #possibly string is faster here?
    original_node: int


class PathOptimizer:
    def __init__(self, graph, start_vertex, seed, node_id_to_name_safe_func=None, get_gaf_string_func=None, name_to_node_id_dict=None):
        self.graph = graph
        self.start_vertex = start_vertex
        self.seed = seed
        self.node_id_to_name_safe = node_id_to_name_safe_func
        self.get_gaf_string = get_gaf_string_func
        self.name_to_node_id = name_to_node_id_dict or {}
        self.traversing_path = self.generate_random_eulerian_path()
        #Storing for debug only
        

    def generate_random_eulerian_path(self):
        #Supplementary for eulerian path construction: selects random outgoing edge, returns it, removes it from graph
        def next_element(temp_graph, current_vertex):
            edges = list(temp_graph.out_edges(current_vertex, data=True, keys=True))
            # Randomly choose the next edge based on seed
            random.shuffle(edges)
            chosen_edge = edges[0]  # Take the first edge after shuffling    
            # Extract edge data
            u, v, key, data = chosen_edge
            temp_graph.remove_edge(u, v, key)
            return (EdgeDescription(u, v, int(key.split("_")[0])), v)        
        
        random.seed(self.seed)
        if nx.has_eulerian_path(self.graph, self.start_vertex):
            # Initialize the path
            path = []
            current = self.start_vertex
            temp_graph = self.graph.copy()
            total_edges = self.graph.number_of_edges()
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
                        if temp_graph.out_degree(path[i].source) > 0:
                            add_cycle = []
                            cycle_current = path[i].source
                            while temp_graph.out_degree(cycle_current) > 0:
                                next_el, cycle_current = next_element(temp_graph, cycle_current)
                                add_cycle.append(next_el)
                            logging.debug (f"adding cycle of length {len(add_cycle)}")
                            logging.debug (f"{add_cycle}")
                            logging.debug(f"after {path[:i]}")
                            path = path[:i] + add_cycle + path[i:]                        

                # Move to next node
            logging.info(f"Randomized Eulerian path found with seed {self.seed} with {len(path)} nodes !")
            if self.get_gaf_string:
                logging.debug (f"{self.get_gaf_string(path)}")
            return path
        else:
            logging.error(f"Path not found from {self.start_vertex}")
            for v in self.graph.nodes():
                logging.info(f"{v}: in {self.graph.in_degree[v]} out {self.graph.out_degree[v]}")
            for v in self.graph.nodes():
                if self.graph.in_degree[v] != self.graph.out_degree[v]:
                    logging.warning(f"Not Eulerian {v}: in {self.graph.in_degree[v]} out {self.graph.out_degree[v]}")
            exit(1)
            return []

    def rc_path(self, path, rc_vertex_map):
        new_path = []
        for edge in path:
            new_u = rc_vertex_map[edge.source]
            new_v = rc_vertex_map[edge.target]
            new_path.append(EdgeDescription(new_v, new_u, -edge.original_node))
        new_path.reverse()
        return new_path

    def get_path(self):
        return self.traversing_path

    def set_path(self, new_path):
        self.traversing_path = new_path

    def get_random_change(self, iter, rc_vertex_map):
        """
        Finds two non-overlapping intervals in the Eulerian path that start and end at the same vertices
        and swaps them to create a new path or inverts a random self-rc interval
        
        :param path: List of edges representing the Eulerian path.
        :param iter: Iteration number for random seed
        :return: Modified path with swapped intervals
        """
        
        random.seed(iter)
        path_length = len(self.traversing_path)
        
        # Pre-build index of matching start/end positions for faster lookup
        start_positions = {}
        end_positions = {}

        for idx, edge_descr in enumerate(self.traversing_path):
            if edge_descr.source not in start_positions:
                start_positions[edge_descr.source] = []
            start_positions[edge_descr.source].append(idx)

            if edge_descr.target not in end_positions:
                end_positions[edge_descr.target] = []
            end_positions[edge_descr.target].append(idx)

        max_tries = 10000  # Adjust based on path length
        
        for _ in range(max_tries):
            # Select i and j such that i < j
            i = random.randint(0, path_length - 3)
            j = random.randint(i + 1,  path_length - 2) 

            start_node = self.traversing_path[i].source
            end_node = self.traversing_path[j].target
            #We can do inversion
            if rc_vertex_map[start_node] == end_node and iter % 2 == 0:
                # not allowing to invert AUX node            
                if "AUX" in self.name_to_node_id:
                    forbidden = False
                    aux_id = self.name_to_node_id["AUX"]
                    for ind in range(i, j + 1):
                        if self.traversing_path[ind].original_node == aux_id:
                            forbidden = True
                            break
                    if forbidden:
                        logging.debug(f"Skipping inversion due to AUX node presence between {i} and {j}")
                        continue
                # Invert the interval from i to j
                new_path = self.traversing_path[:i] + self.rc_path(self.traversing_path[i:j + 1], rc_vertex_map) + self.traversing_path[j + 1:]

                logging.debug(f"{_ + 1} attempts to generate random inversion")
                if self.get_gaf_string:
                    logging.debug(f"{self.get_gaf_string(new_path)}")
                return new_path

            # Find k where path[k].source == start_node and k > j
            k_candidates = [k for k in start_positions.get(start_node, []) if k > j]
            if not k_candidates:
                continue
            # Choose k randomly without temperature weighting
            k = random.choice(k_candidates)

            # Find l such that path[l].target == end_node and l > k
            l_candidates = [l for l in end_positions.get(end_node, []) if l > k]
            if not l_candidates:
                continue
            
            # Choose l randomly without temperature weighting
            l = random.choice(l_candidates)
            # Swap the intervals
            new_path = (
                self.traversing_path[:i]
                + self.traversing_path[k:l + 1]
                + self.traversing_path[j + 1:k]
                + self.traversing_path[i:j + 1]
                + self.traversing_path[l + 1:]
            )
            logging.debug(f"{_ + 1} attempts to generate random swap")
            if self.get_gaf_string:
                logging.debug(f"{self.get_gaf_string(new_path)}")
            return new_path    
        logging.warning("Failed to find valid intervals to swap")
        if self.get_gaf_string:
            logging.warning(f"{self.get_gaf_string(self.traversing_path)}")
        #exit(9)
        return self.traversing_path