#!/usr/bin/env python3

import logging
import random
import networkx as nx

from logging_utils import log_assert
from node_mapper import NodeIdMapper
from path_supplementary import EdgeDescription, get_gaf_path, get_gaf_string, rc_path


class PathOptimizer:
    def __init__(self, graph, start_vertex, seed, node_mapper=None):
        self.graph = graph
        self.start_vertex = start_vertex
        self.seed = seed
        self.node_mapper = node_mapper
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
            logging.debug (f"{get_gaf_string(path, self.node_mapper)}")
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
                if self.node_mapper and self.node_mapper.has_name("AUX"):
                    forbidden = False
                    aux_id = self.node_mapper.get_id_for_name("AUX")
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
                logging.debug(f"{get_gaf_string(new_path, self.node_mapper)}")
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
            logging.debug(f"{get_gaf_string(new_path, self.node_mapper)}")
            return new_path
        logging.warning("Failed to find valid intervals to swap")
        logging.warning(f"{get_gaf_string(self.traversing_path, self.node_mapper)}")
        #exit(9)
        return self.traversing_path

    def get_synonymous_changes(self, path, rc_vertex_map, alignment_scorer):
        """
        Find synonymous changes in the path - intervals that can be swapped or inverted 
        without changing the alignment score.
        """
        # Import necessary functions that might be needed
        
        #some copy-paste
        path_length = len(path)    
        # Pre-build index of matching start/end positions for faster lookup
        start_positions = {}
        end_positions = {}
        final_score = alignment_scorer.score_corasick(path)
        for idx, edge_descr in enumerate(path):
            if edge_descr.source not in start_positions:
                start_positions[edge_descr.source] = []
            start_positions[edge_descr.source].append(idx)

            if edge_descr.target not in end_positions:
                end_positions[edge_descr.target] = []
            end_positions[edge_descr.target].append(idx)

        swappable_intervals = set()
        for start_v in start_positions:
            for end_v in end_positions:
                for first_start_ind in range(len(start_positions[start_v])-1):
                    for second_start_ind in range(first_start_ind + 1, len(start_positions[start_v])):
                        for first_end_ind in range(len(end_positions[end_v]) - 1):
                            for second_end_ind in range(first_end_ind + 1, len(end_positions[end_v])):
                                #valid swap
                                start_path_first_ind = start_positions[start_v][first_start_ind]
                                end_path_first_ind = end_positions[end_v][first_end_ind]
                                start_path_second_ind = start_positions[start_v][second_start_ind]
                                end_path_second_ind = end_positions[end_v][second_end_ind]
                                #for an interval start and end can be same - one node paths. But different intervals should not overlap
                                if end_path_first_ind < start_path_first_ind or end_path_second_ind < start_path_second_ind or start_path_second_ind <= end_path_first_ind:                                
                                    continue
                                first_gaf_fragment = get_gaf_path(path[start_path_first_ind:end_path_first_ind+1], self.node_mapper)
                                second_gaf_fragment = get_gaf_path(path[start_path_second_ind:end_path_second_ind+1], self.node_mapper)
                                logging.debug(first_gaf_fragment)
                                logging.debug(second_gaf_fragment)
                                if first_gaf_fragment == second_gaf_fragment:
                                    logging.debug(f"Found synonymous change: {first_gaf_fragment} positions {start_path_first_ind}-{end_path_first_ind} and {start_path_second_ind}-{end_path_second_ind}, not checking")
                                    continue
                                new_path = (
                                path[:start_path_first_ind]
                                    + path[start_path_second_ind:end_path_second_ind + 1]
                                    + path[end_path_first_ind + 1:start_path_second_ind]
                                    + path[start_path_first_ind:end_path_first_ind + 1]
                                    + path[end_path_second_ind + 1:]
                                )
                                log_assert((len(new_path) == len(path)), f"New path length does not match original path len = {len(path)} indices {start_path_first_ind}-{end_path_first_ind} and {start_path_second_ind}-{end_path_second_ind}")
                                new_score = alignment_scorer.score_corasick(new_path)
                                log_assert(new_score <= final_score, "New path score is greater than original path score")
                                logging.debug(f"New path score: {new_score}, original path score: {final_score} positions {start_path_first_ind}-{end_path_first_ind} and {start_path_second_ind}-{end_path_second_ind}")
                                if new_score < final_score:
                                    logging.debug(f"Swap {start_path_first_ind}-{end_path_first_ind} with {start_path_second_ind}-{end_path_second_ind} decreases score, continuing")
                                elif new_score == final_score:

                                    left_shift = 0
                                    right_shift = 0
                                    while first_gaf_fragment[left_shift] == second_gaf_fragment[left_shift]:
                                        left_shift += 1
                                    while first_gaf_fragment[-right_shift-1] == second_gaf_fragment[-right_shift-1]:
                                        right_shift += 1
                                    logging.debug(f"Swap {start_path_first_ind}-{end_path_first_ind} with {start_path_second_ind}-{end_path_second_ind} does not change score")
                                    logging.debug(f"Tuned intervals {start_path_first_ind + left_shift}-{end_path_first_ind - right_shift} with {start_path_second_ind + left_shift}-{end_path_second_ind - right_shift}")
                                    swappable_intervals.add((start_path_first_ind + left_shift, end_path_first_ind - right_shift, start_path_second_ind + left_shift, end_path_second_ind - right_shift))
                                    logging.debug(f"Edge paths are {first_gaf_fragment} and {second_gaf_fragment}")
        invertable_intervals = set()
        for start_v in start_positions:
            rc_start_v = rc_vertex_map[start_v]
            if not rc_start_v in end_positions:
                continue
            for start_path_ind in start_positions[start_v]:
                for end_path_ind in end_positions[rc_start_v]:
                    #check if we can invert a self-rc interval
                    if start_path_ind < end_path_ind:
                        #invert the interval
                        # Use the static rc_path from tangle_traverser for this calculation
                        inverted_path = rc_path(path[start_path_ind:end_path_ind + 1], rc_vertex_map)
                        if get_gaf_string(inverted_path, self.node_mapper) == get_gaf_string(path[start_path_ind:end_path_ind + 1], self.node_mapper):
                            logging.debug(f"Found equivalent inverted path: {get_gaf_string(inverted_path, self.node_mapper)}")
                            continue
                        new_path = path[:start_path_ind] + rc_path(path[start_path_ind:end_path_ind + 1], rc_vertex_map) + path[end_path_ind + 1:]
                        new_score = alignment_scorer.score_corasick(new_path)
                        log_assert(new_score <= final_score, "New path score is greater than original path score")
                        logging.debug(f"New path score: {new_score}, original path score: {final_score} positions {start_path_ind}-{end_path_ind}")
                        if new_score < final_score:
                            logging.debug(f"Inversion {start_path_ind}-{end_path_ind} decreases score, continuing")
                        elif new_score == final_score:
                            invertable_intervals.add((start_path_ind, end_path_ind))
                            logging.debug(f"Inversion {start_path_ind}-{end_path_ind} does not change score")
                            logging.debug(f"Edge paths are {get_gaf_string(path[start_path_ind:end_path_ind + 1], self.node_mapper)}")
        #TODO: check for possible inversions
        if len(invertable_intervals) + len(swappable_intervals) > 0:
            if len(swappable_intervals) > 0:
                logging.warning(f"Total {len(swappable_intervals)} path swaps with the same score found!")
                for first_start, first_end, second_start, second_end in swappable_intervals:
                    logging.info(f"Swappable interval: {first_start}-{first_end} with {second_start}-{second_end}")
                    logging.info(f"Subpaths {get_gaf_string(path[first_start:first_end + 1], self.node_mapper)} and {get_gaf_string(path[second_start:second_end + 1], self.node_mapper)}")
            if len(invertable_intervals) > 0:
                logging.warning(f"Total {len(invertable_intervals)} path inversions with the same score found!")
                for start, end in invertable_intervals:
                    logging.info(f"Invertable interval: {start}-{end}")
                    logging.info(f"Subpath {get_gaf_string(path[start:end + 1], self.node_mapper)}")
        else:
            logging.info("No equivalent paths found")