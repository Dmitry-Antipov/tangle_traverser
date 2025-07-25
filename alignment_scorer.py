#!/usr/bin/env python3
import logging
import ahocorasick


class AlignmentScorer:
    #TODO: deprioritize based on lengths and not just node counts
    #With othervise equivaltent solution we belive more in one with smaller inversion
    DEPRIORITIZE_RC_COEFFICIENT = 0.9

    #TODO:
    #add coefficient to balance path lengths for two-haplotype tangles
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
        return "," + ",".join(str(edge.original_node) for edge in path) + ","
    
    def aln_to_string(self, aln):
        return "," + ",".join(str(node) for node in aln) + ","
    
    def score_corasick(self, path):
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
