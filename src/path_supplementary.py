#!/usr/bin/env python3

import dataclasses

from .node_id_mapper import NodeIdMapper

@dataclasses.dataclass
class EdgeDescription:
    source: int
    target: int
    #possibly string is faster here?
    original_node: int


def get_gaf_path(path, node_mapper):
    """Convert a path of EdgeDescriptions to a list of GAF node representations."""
    # Import here to avoid circular imports   
    res = []
    for i in range(len(path)):
        edge = path[i].original_node
        node = node_mapper.node_id_to_name_safe(edge)
        res.append(node)
    return res


def get_gaf_string(path, node_mapper):
    """Convert a path of EdgeDescriptions to a GAF string representation."""
    path_arr = get_gaf_path(path, node_mapper)
    res = "".join(path_arr)
    #logging.info(res)
    return res


def rc_path(path, rc_vertex_map):
    """Create reverse complement path from given path."""
    new_path = []
    for edge in path:
        new_u = rc_vertex_map[edge.source]
        new_v = rc_vertex_map[edge.target]
        new_path.append(EdgeDescription(new_v, new_u, -edge.original_node))
    new_path.reverse()
    return new_path
