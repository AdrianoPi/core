import sys
import time
import random
import json
import datetime

ONLY_PRINT_CONNECTIONS = False

# This function reads a file from the topology path.
def read_topology_from_file(path):
    topology = []
    init_potentials = []
    iexts = []
    syn_weights = []
    syn_delays = []
    exc_syn_weight = 0
    exc_syn_delay = 0
    inh_syn_weight = 0
    inh_syn_delay = 0
    
    with open(path, "r") as f:
        iter_lines = iter(f.readlines())
        next(iter_lines)
        for l in iter_lines:
            toks = l.split(", ")
            nrn_id = int(toks[0])
            nrn_v0 = float(toks[1])
            nrn_ie = float(toks[2])
            syn_w = float(toks[3])
            syn_delay = float(toks[4])
            nrn_outs = [int(i.strip()) for i in toks[5:] if len(i.strip()) != 0]
            while nrn_id >= len(topology):
                topology.append(None)
                init_potentials.append(None)
                iexts.append(None)
            topology[nrn_id] = nrn_outs
            init_potentials[nrn_id] = nrn_v0
            # Only the first entry is read to set current to input population
            iexts[nrn_id] = nrn_ie
            
            if syn_w > 0:
                exc_syn_weight = syn_w
                exc_syn_delay = syn_delay
            else:
                inh_syn_weight = syn_w
                inh_syn_delay = syn_delay
    
    print("exc_syn_weight, exc_syn_delay, inh_syn_weight, inh_syn_delay\n", exc_syn_weight, exc_syn_delay, inh_syn_weight, inh_syn_delay)
    return topology, init_potentials, iexts, exc_syn_weight, exc_syn_delay, inh_syn_weight, inh_syn_delay

#TODO: use the correct neuron model and physical parameters!!
if __name__ == "__main__":
        
    if len(sys.argv)!=5:
        print("Four arguments needed: path/to/topology/file, number of threads, time to simulate, timestep (dt)",file=sys.stderr)
        exit(0)

    import nest

    t_file = sys.argv[1]
    topology, init_potentials, iexts, exc_weight, exc_delay, inh_weight, inh_delay = read_topology_from_file(t_file)
    print(f"Will use the following parameters:")
    print(f"Iext for input population = {iexts[0]}")
    print(f"Excitatory Weight {exc_weight}, delay {exc_delay}")
    print(f"Inhibitory Weight {inh_weight}, delay {inh_delay}")

    # Number of neurons
    nrns = len(topology)
    # Number of local threads
    threads = int(sys.argv[2])
    # Stop time
    tstop = int(sys.argv[3])
    # Simulation timestep (ms)
    # ~ dt = 0.1
    dt = float(sys.argv[4])
    print("Will simulate", nrns, "neurons for", tstop, "msec on", threads, "workers, with timestep ", dt)

    # Population sizes
    # In,  1e,  1i,  2e,  2i, out
    popsizes = [100, 200, 200, 200, 200, 100]
    pops = []

    min_delay = min(exc_delay, inh_delay)
    max_delay = max(exc_delay, inh_delay)

    kernel_params = {
        'resolution'		: dt,
        'min_delay'		    : min_delay,
        'max_delay'		    : max_delay,
        'local_num_threads'	: threads
    }

    nest.SetKernelStatus(kernel_params)

    # Neuron Physical parameters
    tau_m    = 10.  # (ms)
    tau_exc  = 0.5  # (ms)
    tau_inh  = 0.5  # (ms)
    E_leak   = -65.	# (mV)
    v_reset  = -65. # (mV)
    v_thresh = -50. # (mV)
    C_m      = 250. # (µF/cm²)
    t_refrac = 2.   # (ms) (clamped at v_reset)
    v_m      = -58. # (mV) starting membrane potential

    # Neuron type and parameters
    celltype = "iaf_psc_exp"    
    cell_params = {
        'tau_m' : tau_m,    'tau_syn_ex' : tau_exc, 'tau_syn_in' : tau_inh,
        'E_L' : E_leak,     'V_reset' : v_reset,    'V_th': v_thresh,
        'C_m' : C_m,        't_ref' : t_refrac,
        'I_e' : 0.0,        'V_m' : v_m
    }

    nest.SetDefaults(celltype, cell_params)

    # Now we can start building the network

    # Start the timer
    t0 = time.time()

    # Init a single big population, with random starting potential -> needs to be not random.
    for popsize in popsizes:
        pops.append(nest.Create(celltype, popsize))

    # Set input population's input current.
    pops[0].set(I_e=iexts[0])

    # Create connections with all_to_all as we provide 1 source neuron, and many dst
    # Set static synapse and the receptor type.
    conn_dict = {'rule': 'one_to_one', "allow_multapses":True}
    nest.CopyModel("static_synapse","exc_syn",{"weight": exc_weight, "delay": exc_delay})
    nest.CopyModel("static_synapse","inh_syn",{"weight": inh_weight, "delay": inh_delay})


    # This should connect everything the way we want: 
    for src, dsts in enumerate(topology):
        is_excitatory = src<300 or 500 <= src < 700 or 900 <= src
        if not len(dsts)>0:
            continue
        
        if is_excitatory:
            nest.Connect([src+1]*len(dsts),  dsts, conn_dict, syn_spec="exc_syn") 
        else:
            nest.Connect([src+1]*len(dsts),  dsts, conn_dict, syn_spec="inh_syn") 

    connections = nest.GetConnections()
    
    if ONLY_PRINT_CONNECTIONS:
        with open(f"nest_connections_{str(datetime.datetime.now()).replace(' ', '')}", "w") as f:
            connections_gotten = connections.get(["source", "target"])
            # ~ conn_list = list(zip(connections_gotten["source"], connections_gotten["target"]))
            adj_list = [[] for _ in range(nrns)]
            for c in zip(connections_gotten["source"], connections_gotten["target"]):
                adj_list[c[0]-1].append(c[1])
            json.dump(adj_list, f)
            exit()
    
    # Spike recorder.
    print("\nSpikes will be recorded to disk")
    nest.overwrite_files = True
    sr = nest.Create('spike_recorder', params={'record_to': 'ascii'})
    nest.Connect(range(901, 1001), sr, syn_spec='exc_syn')

    buildCPUTime = time.time() - t0

    # Run the simulation!
    nest.Simulate(tstop)

    simCPUTime = time.time() - (t0+buildCPUTime)

    # printa le stats
    print("\n--- 1000 Neurons Network Simulation ---"			)
    # ~ print("Number of Neurons      : %d" % nrns				    )
    # ~ print("Number of Synapses     : %s" % n_syns				)
    print("Build time             : %g s" % buildCPUTime			)
    print("Simulation time        : %g s" % simCPUTime  			)
