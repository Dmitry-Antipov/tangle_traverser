#!/usr/bin/env python3

from importlib.resources import path
import logging
import random
import networkx as nx

from .logging_utils import log_assert
from .node_id_mapper import NodeIdMapper
from .path_supplementary import EdgeDescription, get_gaf_path, get_gaf_string, rc_path
from .input_parsing import reverse_complement

# Add this length (max) from the nodes neighboring tangle to fasta for better alignment
UNIQUE_BORDER_LENGTH = 200000


class PathOptimizer:
    def __init__(self, graph, start_vertex,  node_mapper, rc_vertex_map):
        self.graph = graph
        self.start_vertex = start_vertex
        self.node_mapper = node_mapper
        #map reverse complement vertices
        self.rc_vertex_map = rc_vertex_map
        self.generate_random_eulerian_path(0)

        #indices from vertices to path positions
        self.start_positions = {}
        self.end_positions = {}
    
    def set_seed(self, seed):
        self.seed = seed
        random.seed(seed)

    def generate_random_eulerian_path(self, seed):
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
        self.seed = seed        
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
            self.traversing_path = path
            self.__update_start_end_positions()
        else:
            logging.error(f"Path not found from {self.start_vertex}")
            for v in self.graph.nodes():
                logging.info(f"{v}: in {self.graph.in_degree[v]} out {self.graph.out_degree[v]}")
            for v in self.graph.nodes():
                if self.graph.in_degree[v] != self.graph.out_degree[v]:
                    logging.warning(f"Not Eulerian {v}: in {self.graph.in_degree[v]} out {self.graph.out_degree[v]}")
            exit(1)
            return []

    def get_path(self):
        return self.traversing_path

    def set_path(self, new_path):
        self.traversing_path = new_path
        self.__update_start_end_positions()    
    
    def __update_start_end_positions(self):
        self.start_positions = {}
        self.end_positions = {}

        for idx, edge_descr in enumerate(self.traversing_path):
            if edge_descr.source not in self.start_positions:
                self.start_positions[edge_descr.source] = []
            self.start_positions[edge_descr.source].append(idx)

            if edge_descr.target not in self.end_positions:
                self.end_positions[edge_descr.target] = []
            self.end_positions[edge_descr.target].append(idx)
        return

    def get_random_change(self, iter):
        """
        Finds two non-overlapping intervals in the Eulerian path that start and end at the same vertices
        and swaps them to create a new path or inverts a random self-rc interval
        
        :param iter: Iteration number for random seed
        :return: Modified path with swapped intervals
        """
        
        random.seed(iter)
        path_length = len(self.traversing_path)            
        self.__update_start_end_positions()
        max_tries = 10000  # Adjust based on path length
        
        for _ in range(max_tries):
            # Select i and j such that i < j
            i = random.randint(0, path_length - 3)
            start_vertex = self.traversing_path[i].source            
            
            rc_start_v = self.rc_vertex_map[start_vertex]
            if iter % 2 == 0 and rc_start_v in self.end_positions:
                # We can do inversion
                j_candidates = [j for j in self.end_positions[rc_start_v] if j > i]
                if not j_candidates:
                    continue
                j = random.choice(j_candidates)
                # not allowing to invert AUX node
                forbidden = False
                if self.node_mapper and self.node_mapper.has_name("AUX"):                
                    aux_id = self.node_mapper.get_id_for_name("AUX")
                    for ind in range(i, j + 1):
                        if self.traversing_path[ind].original_node == aux_id:
                            forbidden = True
                            break
                if forbidden:
                    logging.debug(f"Skipping inversion due to AUX node presence between {i} and {j}")
                    continue
                
                # Invert the interval from i to j
                new_path = self.traversing_path[:i] + rc_path(self.traversing_path[i:j + 1], self.rc_vertex_map) + self.traversing_path[j + 1:]
                logging.debug(f"{_ + 1} attempts to generate random inversion")                
                logging.debug(f"{get_gaf_string(new_path, self.node_mapper)}")
                return new_path
            else: 
                #Doing regular interval swap
                j = random.randint(i + 1,  path_length - 2) 
                end_vertex = self.traversing_path[j].target

                k_candidates = [k for k in self.start_positions.get(start_vertex, []) if k > j]
                if not k_candidates:
                    continue
                k = random.choice(k_candidates)

                l_candidates = [l for l in self.end_positions.get(end_vertex, []) if l > k]
                if not l_candidates:
                    continue
                
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
        return self.traversing_path

    #Only for logging
    def get_synonymous_changes(self, alignment_scorer):
        self.__update_start_end_positions()
        swappable_intervals = set()
        final_score = alignment_scorer.score_corasick(self.traversing_path)
        for start_v in self.start_positions:
            for end_v in self.end_positions:
                for first_start_ind in range(len(self.start_positions[start_v])-1):
                    for second_start_ind in range(first_start_ind + 1, len(self.start_positions[start_v])):
                        for first_end_ind in range(len(self.end_positions[end_v]) - 1):
                            for second_end_ind in range(first_end_ind + 1, len(self.end_positions[end_v])):
                                #valid swap
                                start_path_first_ind = self.start_positions[start_v][first_start_ind]
                                end_path_first_ind = self.end_positions[end_v][first_end_ind]
                                start_path_second_ind = self.start_positions[start_v][second_start_ind]
                                end_path_second_ind = self.end_positions[end_v][second_end_ind]
                                #for an interval start and end can be same - one node paths. But different intervals should not overlap
                                if end_path_first_ind < start_path_first_ind or end_path_second_ind < start_path_second_ind or start_path_second_ind <= end_path_first_ind:                                
                                    continue
                                first_gaf_fragment = get_gaf_path(self.traversing_path[start_path_first_ind:end_path_first_ind+1], self.node_mapper)
                                second_gaf_fragment = get_gaf_path(self.traversing_path[start_path_second_ind:end_path_second_ind+1], self.node_mapper)
                                logging.debug(first_gaf_fragment)
                                logging.debug(second_gaf_fragment)
                                if first_gaf_fragment == second_gaf_fragment:
                                    logging.debug(f"Found synonymous change: {first_gaf_fragment} positions {start_path_first_ind}-{end_path_first_ind} and {start_path_second_ind}-{end_path_second_ind}, not checking")
                                    continue
                                new_path = (
                                    self.traversing_path[:start_path_first_ind]
                                    + self.traversing_path[start_path_second_ind:end_path_second_ind + 1]
                                    + self.traversing_path[end_path_first_ind + 1:start_path_second_ind]
                                    + self.traversing_path[start_path_first_ind:end_path_first_ind + 1]
                                    + self.traversing_path[end_path_second_ind + 1:]
                                )
                                if get_gaf_string(new_path, self.node_mapper) == get_gaf_string(self.traversing_path, self.node_mapper):
                                    logging.debug(f"Synonymous change, different fragments but same path: {first_gaf_fragment} positions {start_path_first_ind}-{end_path_first_ind} and {start_path_second_ind}-{end_path_second_ind}, not checking")
                                    continue
                                log_assert((len(new_path) == len(self.traversing_path)), f"New path length does not match original path len = {len(self.traversing_path)} indices {start_path_first_ind}-{end_path_first_ind} and {start_path_second_ind}-{end_path_second_ind}")
                                new_score = alignment_scorer.score_corasick(new_path)
                                log_assert(new_score <= final_score, "New path score is greater than original path score")
                                logging.debug(f"New path score: {new_score}, original path score: {final_score} positions {start_path_first_ind}-{end_path_first_ind} and {start_path_second_ind}-{end_path_second_ind}")
                                if new_score < final_score:
                                    logging.debug(f"Swap {start_path_first_ind}-{end_path_first_ind} with {start_path_second_ind}-{end_path_second_ind} decreases score, continuing")
                                elif new_score == final_score:

                                    left_shift = 0
                                    right_shift = 0

                                    min_l = min(len(first_gaf_fragment), len(second_gaf_fragment))
                                    while left_shift < min_l and first_gaf_fragment[left_shift] == second_gaf_fragment[left_shift]:
                                        left_shift += 1
                                    if left_shift == min_l:
                                        left_shift -= 1

                                    while right_shift < min_l and first_gaf_fragment[-right_shift-1] == second_gaf_fragment[-right_shift-1]:
                                        right_shift += 1
                                    if right_shift == min_l:
                                        right_shift -= 1

                                    logging.debug(f"Swap {start_path_first_ind}-{end_path_first_ind} with {start_path_second_ind}-{end_path_second_ind} does not change score")
                                    logging.debug(f"Tuned intervals {start_path_first_ind + left_shift}-{end_path_first_ind - right_shift} with {start_path_second_ind + left_shift}-{end_path_second_ind - right_shift}")
                                    swappable_intervals.add((start_path_first_ind + left_shift, end_path_first_ind - right_shift, start_path_second_ind + left_shift, end_path_second_ind - right_shift))
                                    logging.debug(f"Edge paths are {first_gaf_fragment} and {second_gaf_fragment}")
        invertable_intervals = set()
        for start_v in self.start_positions:
            rc_start_v = self.rc_vertex_map[start_v]
            if not rc_start_v in self.end_positions:
                continue
            for start_path_ind in self.start_positions[start_v]:
                for end_path_ind in self.end_positions[rc_start_v]:
                    #check if we can invert a self-rc interval
                    if start_path_ind < end_path_ind:
                        #invert the interval
                        # Use the static rc_path from tangle_traverser for this calculation
                        inverted_path = rc_path(self.traversing_path[start_path_ind:end_path_ind + 1], self.rc_vertex_map)
                        if get_gaf_string(inverted_path, self.node_mapper) == get_gaf_string(self.traversing_path[start_path_ind:end_path_ind + 1], self.node_mapper):
                            logging.debug(f"Found equivalent inverted path: {get_gaf_string(inverted_path, self.node_mapper)}")
                            continue
                        new_path = self.traversing_path[:start_path_ind] + rc_path(self.traversing_path[start_path_ind:end_path_ind + 1], self.rc_vertex_map) + self.traversing_path[end_path_ind + 1:]
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
                    logging.info(f"Subpaths {get_gaf_string(self.traversing_path[first_start:first_end + 1], self.node_mapper)} and {get_gaf_string(self.traversing_path[second_start:second_end + 1], self.node_mapper)}")
            if len(invertable_intervals) > 0:
                logging.warning(f"Total {len(invertable_intervals)} path inversions with the same score found!")
                for start, end in invertable_intervals:
                    logging.info(f"Invertable interval: {start}-{end}")
                    logging.info(f"Subpath {get_gaf_string(self.traversing_path[start:end + 1], self.node_mapper)}")
        else:
            logging.info("No equivalent paths found")

    def output_path(self, original_graph, output_fasta, output_gaf):
        """Output the best path to FASTA and GAF files."""
        aux = -1
        if self.node_mapper.has_name("AUX"):
            aux_id = self.node_mapper.get_id_for_name("AUX")
            for i in range(len(self.traversing_path)):
                if abs(self.traversing_path[i].original_node) == aux_id:
                    aux = i
                    logging.debug(f"Found AUX at position {aux}")
                    break
        if aux > 0:
            paths = [self.traversing_path[:aux - 1], self.traversing_path[aux + 1:]]
        else:
            paths = [self.traversing_path]

        gaf_file = open(output_gaf, 'w')
        count = 0
        for path in paths:
            gaf_file.write(f"traversal_{count}\t{get_gaf_string(path, self.node_mapper)}\n")
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
                        logging.info("Provided noseq assembly graph, no HPC fasta output possible")
                        return
                    if not orientation:
                        # Reverse complement the sequence if orientation is negative
                        node_sequence = reverse_complement(node_sequence)

                    if i == 0:
                        overlap = max (0, len(node_sequence) - UNIQUE_BORDER_LENGTH)                
                    else:                
                        overlap = original_graph.get_edge_data(last_node, edge_id)['overlap']

                    if i == len(self.traversing_path) - 1:                            
                        contig_sequence += node_sequence[overlap:min(len(node_sequence),UNIQUE_BORDER_LENGTH)]
                    else:
                        contig_sequence += node_sequence[overlap:]
                    last_node = edge_id

                # Write the contig to the FASTA file
                fasta_file.write(f">traversal_{count}\n")
                fasta_file.write(f"{contig_sequence}\n")
                count += 1
        logging.info("Homopolymer-compressed fasta output finished")