import random
import math
import sys
"""
This is an ugly script with a lot of magic numbers that generates a topology
file with 1000 neurons. The topology graph presents no loops thanks to the 
connection probability table.
The populations are 6: Input, Exc1, Inh1, Exc2, Inh2, Output
"""

# NeuronID to population index
def n2pop(neuron_id, popsizes):
    for i, sz in enumerate(popsizes):
        if neuron_id < sz:
            return i
        neuron_id -= sz

# One input that fires to every other population with different connection probs
# 2 middle layers, with 2 populations each
# 1 output layer
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Must provide path to output file")
        exit(1)
    outpath = sys.argv[1]
    
    INPUT_CURRENT = 1800 # Input current to the first 100 neurons
    v_r = -65.0 # Reset potential
    v_th = -50.0 # Fire threshold
    syn_exc_weight = 200 # Excitatory synapse weight
    in_w_g = 3.0
    syn_inh_weight = -(syn_exc_weight * in_w_g) # Inhibitory synapse weight
    syn_exc_delay = 1.5
    syn_inh_delay = 0.8
    
    random.seed(42)
    
    #            In,  1e,  1i,  2e,  2i, out
    popsizes = [100, 200, 200, 200, 200, 100]
    poptypes = ["e", "e", "i", "e", "i", "e"]
    popweights = []
    popdelays = []
    for p in poptypes:
        if p=="e":
            popweights.append(syn_exc_weight)
            popdelays.append(syn_exc_delay)
        else:
            popweights.append(syn_inh_weight)
            popdelays.append(syn_inh_delay)

    pop_offset = [sum(popsizes[:i]) for i in range(len(popsizes))]
    pop_offset.append(sum(popsizes))

    table =[
        [0.0, 0.292, 0.192, 0.049, 0.237, 0.169],
        [0.0, 0.0  , 0.0  , 0.106, 0.254, 0.438],
        [0.0, 0.0  , 0.0  , 0.409, 0.250, 0.309],
        [0.0, 0.0  , 0.0  , 0.0  , 0.0  , 0.491],
        [0.0, 0.0  , 0.0  , 0.0  , 0.0  , 0.225],
        [0.0, 0.0  , 0.0  , 0.0  , 0.0  , 0.0]]

    """ table[i][j] contains connection probability from i to j In -> 
    1e, 1i, 2e, 2i, out 1e -> 2e, 2i, out 1i -> 2e, 2i, out 2e -> out 
    2i -> out """

    connections = []

    for i in range(len(popsizes)):
        print("Source population:", i)
        src_first = pop_offset[i]
        src_last  = src_first + popsizes[i] - 1
        for j in range(len(popsizes)):
            print(j)
            dst_first = pop_offset[j]
            dst_last  = dst_first + popsizes[j] - 1
            print("Dest first:", dst_first, "last:", dst_last)
            nsyns = math.floor(table[i][j]*popsizes[i]*popsizes[j])
            print("Nsyns", nsyns)
            for _ in range(nsyns):
                src = random.randint(src_first, src_last)
                dst = random.randint(dst_first, dst_last)
                connections.append([src, dst])

    topology = [[] for _ in range(sum(popsizes))]

    for c in connections:
        topology[c[0]].append(c[1])

    with open(outpath, "w") as f:
        f.write("Neuron, v0, Iext, syn_weight, syn_delay, [Targets]\n")
        for src, dests in enumerate(topology):
            pop = n2pop(src, popsizes)
            f.write(f"{src}, {round(random.uniform(v_r, v_th), 17)}, {INPUT_CURRENT if src < 100 else 0}, {popweights[pop]}, {popdelays[pop]}, ")
            # ~ f.write(str(src))
            # ~ f.write(", ")
            dst_string = ""
            for dst in dests:
                dst_string += str(dst) + ", "
            f.write(dst_string[:-2]+"\n")

