"""
This script compares the pubsubs SubscribersAdjacencyList from ROOT-Sim
with the topology created by generate_1k_topology.py script
"""
import sys
import json


def read_sim_topology_as_adjacency_list(path):
    sim_adj = []
    with open(path, "r") as f:
        sim_adj = json.load(f)
    for s in sim_adj:
        s.sort()

    return sim_adj

def read_topology_file_as_adjacency_list(path):
    top_adj = []
    with open(path, "r") as f:
        lines = f.readlines()
        lines = lines[1:]
    
    for l in lines:
        vals = l.strip().strip(",").split(",")
        vals = vals[5:]
        vals = [int(v) for v in vals]
        vals.sort()
        top_adj.append(vals)
    
    return top_adj


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Two arguments required: ROOT-Sim/Nest Adjacency list, and the topology")
        exit(1)

    sim_adj_path = sys.argv[1]
    topology_path = sys.argv[2]

    sim_adj = read_sim_topology_as_adjacency_list(sim_adj_path)
    top_adj = read_topology_file_as_adjacency_list(topology_path)
    
    with open("sorted_topology.txt", "w") as f:
        json.dump(top_adj, f)
    
    with open("sorted_sim_topology.txt", "w") as f:
        json.dump(sim_adj, f)

    print(top_adj==sim_adj)
