#!/usr/bin/env python
import sys
import os
import itertools
import networkx as nx

class Node:
    def __init__(self, node_id):
        self.node_id = node_id
        self.edges = set()
        self.path = []
        self.adjusted_coverage = []
        self.utig1_length = []
        #length of maximum overlap in paths between different nodes
        self.max_inter = [0,0]

    def add_edge(self, to_node):
        self.edges.add(to_node)

    def add_path(self, path):
        self.path = path

    def get_coverage(self):        
        if self.max_inter[0]  < len (self.path) - self.max_inter[1]:
            start = self.max_inter[0] 
            end = len (self.path) - self.max_inter[1] 
        else: 
            # print(f"for node {self.node_id} of {len(self.path)} utig1 incident vertices of sizes {self.max_inter} overlaps")
            start = 0
            end = len(self.path)
        #to avoid 0 length issues
        total_l = 1
        total_cov = 0
        #we should use overlap
        for i in range (start, end):
            total_l += self.utig1_length[i]
            total_cov += self.utig1_length[i] * self.adjusted_coverage[i]        
        return total_cov/total_l
    
class AssemblyGraph:
    def __init__(self, gfa_file=None):
        self.nodes = {}
        self.edges = []
        self.intersections = {}
        self.utig1_multiplitcities = {}
        self.utig1_coverage = {}
        self.utig1_length = {}
        self.utig1_subgraphs = {}
        if gfa_file:
            self.parse_gfa(gfa_file)
        
    def add_node(self, node_id):
        if node_id not in self.nodes:
            self.nodes[node_id] = Node(node_id)

    def add_edge(self, from_node, to_node):        
        self.nodes[from_node].add_edge(to_node)        
        self.edges.append((from_node, to_node))

    def add_path(self, path_id, path):
        self.nodes[path_id].add_path(path)

    def add_path_intersections(self):
        for from_node, node in self.nodes.items():
            for to_node in node.edges:
                from_path = self.nodes[from_node].path
                to_path = self.nodes[to_node].path
                intersection = []
                for i in range (1, min(len(from_path), len(to_path)) + 1):
                    if i > 20:
                        break
                    if from_path[-i:] == to_path[:i]:
                        intersection = from_path[-i:]
                        break
                if intersection:
                    if from_node not in self.intersections:
                        self.intersections[from_node] = {}
                    self.intersections[from_node][to_node] = intersection
                    lens = len(intersection)
                    self.nodes[from_node].max_inter[1] = max(lens, self.nodes[from_node].max_inter[1])
                    self.nodes[to_node].max_inter[0] = max(lens, self.nodes[to_node].max_inter[0])

    def get_all_utig1(self):
        self.utig1_multiplicities = {}
        for node in self.nodes.values():
            for utig1 in node.path:
                if utig1 in self.utig1_multiplicities:
                    self.utig1_multiplicities[utig1] += 1
                else:
                    self.utig1_multiplicities[utig1] = 1
        

    def parse_gfa(self, gfa_file):
        rc = {'-':'+','+':'-'}
        to_prefs = {'-':'<','+':'>'}
        with open(gfa_file, 'r') as file:
            for line in file:
                if line.startswith('S'):
                    parts = line.strip().split()
                    node_id = parts[1]
                    self.add_node('>' + node_id)
                    self.add_node('<' + node_id)
                elif line.startswith('L'):
                    parts = line.strip().split()
                    from_node = to_prefs[parts[2]] + parts[1]
                    to_node = to_prefs[parts[4]] + parts[3]
                    rc_to_node = to_prefs[rc[parts[4]]] + parts[3]
                    rc_from_node = to_prefs[rc[parts[2]]] + parts[1]
                    self.add_edge(rc_to_node, rc_from_node)
                    self.add_edge(from_node, to_node)

    def parse_gaf(self, gaf_file):
        with open(gaf_file, 'r') as file:
            for line in file:
                parts = line.strip().split()
                path_id = parts[0]
                if path_id == "node":
                    continue
                nodes = []
                node = ""
                for char in parts[1]:
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
                nodes = filtered_nodes
                self.add_path(">" + path_id , nodes)
                reverse_nodes = ["<" + n[1:] if n.startswith(">") else ">" + n[1:] for n in nodes[::-1]]
                self.add_path("<" + path_id , reverse_nodes)

    def construct_subgraphs(self):
        instances = {}
        sys.stderr.write(f"{len(self.utig1_multiplicities)} different utig1s\n")
        added = 0
        for node in self.nodes:
            path = self.nodes[node].path
            
            for i in range (len(path)):
                utig1 = path[i]
                if self.utig1_multiplicities[utig1] > 1:
                    if not utig1 in instances:
                        instances[utig1] = []
                        self.utig1_subgraphs[utig1] = nx.Graph()
                    instances[utig1].append((node, i))
                    self.utig1_subgraphs[utig1].add_node(f"{node}_{i}")
        for utig1 in instances.keys():
            for i in range(len(instances[utig1])):
                for j in range(0, len(instances[utig1])):
                    node1, i1 = instances[utig1][i]
                    node2, i2 = instances[utig1][j]
                    if node1 in graph.intersections and node2 in graph.intersections[node1]:
                        intersection = graph.intersections[node1][node2]
                        lenpath1 = len(self.nodes[node1].path)
                        if (i2 < len(intersection) and i1 >= lenpath1 - len(intersection) and
                          i2 == i1 - (lenpath1 - len(intersection))):
                            self.utig1_subgraphs[utig1].add_edge(f"{node1}_{i1}", f"{node2}_{i2}")
                            added += 1
        sys.stderr.write(f"Added {added} edges to subgraphs\n")

    def read_utig1_coverage(self, csv_file):   
        with open(csv_file, 'r') as file:
            for line in file:
                parts = line.strip().split()
                if len(parts) == 3:
                    node, coverage, length = parts
                    if coverage.replace('.', '', 1).isdigit():
                        self.utig1_coverage[">" + node] = float(coverage)
                        self.utig1_coverage["<" + node] = float(coverage)
                        self.utig1_length["<" + node] = int(length)
                        self.utig1_length[">" + node] = int(length)
        for node in self.nodes.values():
            for utig1 in node.path:
                if utig1 in self.utig1_coverage:
                    node.adjusted_coverage.append(self.utig1_coverage[utig1])
                    node.utig1_length.append(self.utig1_length[utig1])
                else:
                    node.adjusted_coverage.append(10000.0)
                    node.utig1_length.append(0)
                

    def adjust_coverage(self):
        for utig1 in self.utig1_subgraphs:                    
            subgraph = self.utig1_subgraphs[utig1]
            sys.stderr.write(f"Subgraph edges for utig1 {utig1}: {subgraph.edges()}\n")
            connected_components = list(nx.connected_components(subgraph))
            component_sizes = [len(component) for component in connected_components]
            sys.stderr.write(f"utig1: {utig1}, component sizes: {component_sizes}\n")
            for component in connected_components:
                sys.stderr.write(f"Component: {component}\n")
            total_components = len(connected_components)
            for positioned_utig1 in subgraph.nodes():
                arr = positioned_utig1.split('_')
                node = arr[0]
                ind = int(arr[1])
                self.nodes[node].adjusted_coverage[ind] /= total_components
        
def get_utig_seq(gaf_file):
    utig_seq = {}
    with open(gaf_file, 'r') as file:
        for line in file:
            entities = line.split()
            utig4 = entities[0]
            for entity in entities:
                utig1_nodes = [node for node in entity.replace('<', '>').replace('_','>').split('>') if 'utig1-' in node]
                
                utig_seq[utig4] = utig1_nodes
    return utig_seq


def count_multiplicity(utig_seq):
    multiplicity = {}

    for node, seq in utig_seq.items():
        for element in seq:
            if element in multiplicity:
                multiplicity[element] += 1
            else:
                multiplicity[element] = 1
    return multiplicity

def get_contained_count(utig_seq):
    borders = {}
    for node, seq in utig_seq.items():
        if len(seq) == 1:
            border_nodes = [seq[0]]
        elif len(seq) > 1:
            border_nodes = [seq[0], seq[-1]]
        else:
            continue
        for node in border_nodes:
            if node in borders:
                borders[node] += 1
            else:
                borders[node] = 1
    contained_count = {}
    for node, seq in utig_seq.items():
        if len(seq) == 1 and seq[0] in borders:
            contained_count[seq[0]] = borders[seq[0]]
    return contained_count

def count_updated_coverage(multiplicity, coverage_data, contained_count, gaf_file):
    upd_cov = {}
    with open(gaf_file, 'r') as file:
        for line in file:
            entities = line.split()
            for entity in entities:
                utig1_nodes = [node for node in entity.replace('<', '>').replace('_','>').split('>') if 'utig1-' in node]
                mincov = get_lowest_coverage(utig1_nodes, coverage_data, contained_count, multiplicity)
                upd_cov[entities[0]] = mincov

    return upd_cov
# Example usage


def read_coverage(csv_file):
    coverage_data = {}
    with open(csv_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                node, coverage, length = parts
                if coverage.replace('.', '', 1).isdigit():
                    coverage_data[node] = float(coverage)
    return coverage_data

def get_lowest_coverage(utig1_nodes, coverage_data, contained_count, multiplicity):
    coverages = [100000]
    for node in utig1_nodes:
        if node in coverage_data:
            if node in contained_count:
                multiplier = contained_count[node]
               # print (f"Node {node} contained {multiplier}")
            else:
                multiplier = 1
            if not (node in multiplicity):
                sys.stderr.write (f"Node {node} not in multiplicity\n")
                exit(1)                
            coverages.append(coverage_data[node] * multiplier/multiplicity[node])
    return min(coverages) if coverages else None

def get_all_unique_strings(graph):
    unique_strings = set()
    for path in graph.paths.values():
        unique_strings.update(path)
    return list(unique_strings)

# Example usage
#TODO: split not evenly between components but with respect to the component's ins/outs.
if __name__ == "__main__":
    gaf_file = sys.argv[1]
    graph_file = sys.argv[2]
    csv_file = sys.argv[3]
    graph = AssemblyGraph(graph_file)
    graph.parse_gaf(gaf_file)
    graph.get_all_utig1()
    graph.add_path_intersections()
    graph.read_utig1_coverage(csv_file)
    graph.construct_subgraphs()
    graph.adjust_coverage()
    #print ("Printing adjusted coverage")
    print ("node\tcoverage\t")
    for node_id, node in graph.nodes.items():
        if node_id[0] == ">":
            print (f"{node_id[1:]}\t{node.get_coverage():.1f}")
            #print (f"{node_id} adjusted coverage {node.get_coverage()} {node.path} {node.max_inter}")
    exit(0)
    # Construct subgraphs for all utig1 with multiplicities larger than 3 and output connected component sizes
    utig1_multiplicities = graph.get_all_utig1()
    for utig1, count in utig1_multiplicities.items():
        if count > 3:
            subgraph = construct_subgraph(graph, utig1)
            print(f"Subgraph edges for utig1 {utig1}: {subgraph.edges()}")
            connected_components = list(nx.connected_components(subgraph))
            component_sizes = [len(component) for component in connected_components]
            print(f"utig1: {utig1}, component sizes: {component_sizes}")
            for component in connected_components:
                print(f"Component: {component}")

    coverage_data = read_coverage(csv_file)
    utig_seq = get_utig_seq(gaf_file)
    multiplicity = count_multiplicity(utig_seq)
    contained_count = get_contained_count(utig_seq)
    upd_cov = get_lowest_coverage(multiplicity, coverage_data, contained_count, gaf_file)
    for entity, cov in upd_cov.items():
        if cov is not None:
            print(f"{entity}\t{cov:.1f}")
        else:
            print(f"{entity}\t{0}")
