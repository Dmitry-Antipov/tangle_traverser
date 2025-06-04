#!/usr/bin/env python3


import sys
import re
import os
count = 0


def revnode(n):
	assert len(n) >= 2
	assert n[0] == "<" or n[0] == ">"
	return (">" if n[0] == "<" else "<") + n[1:]

def canon(left, right):
	if revnode(right) + revnode(left) < left + right:
		return (revnode(right), revnode(left))
	return (left, right)
#partial copypaste from verkko's  get_layout_from_mbg

#for the input utig4 outputs utig1 paths (without left/right shifts)

#6-layoutContigs/combined-nodemap.txt
#path_to_run
node_mapping = {}
with open(os.path.join(sys.argv[1], "6-layoutContigs","combined-nodemap.txt")) as f:
	for l in f:
		parts = l.strip().split('\t')
		assert parts[0] not in node_mapping
		path = parts[1].split(':')[0].replace('<', "\t<").replace('>', "\t>").strip().split('\t')
		node_mapping[parts[0]] = path

def get_utig1_paths(node_mappings):
    result = {}
    for node in node_mapping.keys():
        if node[:5] == "utig4":
            result[node] = [">" +node]
    while True:    
        changed = False
        for node in result.keys():
            new_path = []
            for pnode in result[node]:

                if (len(pnode) >5 and pnode[1:6] == "utig1") or not (pnode[1:] in node_mapping):
                    new_path.append(pnode)
                else:
                    changed = True
                    part = [n for n in node_mapping[pnode[1:]]]
                    if pnode[0] == "<":
                        part = [revnode(n) for n in part[::-1]]
                    new_path.extend(part)
            result[node] = new_path
        if not(changed):
            break
    return result
result = get_utig1_paths(node_mapping)
print("node\tutig1-label")
for node in sorted(result.keys()):
    print (f"{node}\t{''.join(result[node])}")                
