#!/usr/bin/env python3
import logging


class NodeIdMapper:
    """Class to handle node ID to name mapping and vice versa."""
    
    def __init__(self):
        # Map node ID to string node names (only positive ids, unoriented)
        self.node_id_to_name = {}
        self.name_to_node_id = {}
        self.last_enumerated_node = 0
    
    def node_id_to_name_safe(self, node_id):
        """Safely convert node ID to name, handling cases where name might not exist"""
        if node_id > 0:
            prefix = ">"
        else:
            prefix = "<"
        node_id = abs(node_id)
        if node_id in self.node_id_to_name:
            return prefix + self.node_id_to_name[node_id]
        else:
            logging.error(f"Node ID {node_id} not found in name mapping.")
            logging.error(f"map size: {len(self.node_id_to_name)}")
            exit(1)
    
    def parse_node_id(self, node_str): 
        """Parse node string to get unique integer ID."""
        if node_str in self.name_to_node_id:
            return self.name_to_node_id[node_str]
        
        parts = node_str.split('-')
        if len(parts) < 2 or not parts[1].isdigit():
            # ribotin graph case
            if node_str.isdigit():            
                node_id = int(node_str)
            else:
                self.last_enumerated_node += 1
                node_id = self.last_enumerated_node
                logging.debug(f"Assigned enumerated ID {node_id} to node {node_str}")
        else:
            node_id = int(parts[1])

        # utig4-0 same as its RC
        if node_id < self.last_enumerated_node:
            self.last_enumerated_node += 1
            node_id = self.last_enumerated_node
            logging.info(f"Assigned enumerated ID {node_id} to node {node_str}, utig4-0 special case")
        
        self.node_id_to_name[node_id] = node_str
        self.name_to_node_id[node_str] = node_id
        return node_id
    
    def add_mapping(self, node_id, node_str):
        """Manually add a mapping between node ID and name."""
        self.node_id_to_name[node_id] = node_str
        self.name_to_node_id[node_str] = node_id
    
    def get_name_to_id_dict(self):
        """Get the name to ID dictionary."""
        return self.name_to_node_id.copy()
    
    def has_name(self, name):
        """Check if a name exists in the mapping."""
        return name in self.name_to_node_id
    
    def get_id_for_name(self, name):
        """Get ID for a given name."""
        return self.name_to_node_id.get(name)
    
    def size(self):
        """Get the size of the mapping."""
        return len(self.node_id_to_name)
