#include "cuba.h"
#include <stdio.h>

//~ #define M_DEBUG                // Undefine to decrease verbosity

#ifdef M_DEBUG

# define printdbg(...) \
    printf(__VA_ARGS__)
#else
# define printdbg(...)
#endif

#undef M_DEBUG

#define MODEL_MAX_SIMTIME global_config.termination_time

#define V_TOLERANCE 0.01	//mV
#define I_TOLERANCE 0.01	//pA
#define T_TOLERANCE 0.001	//mS

#define e_portion 0.8f

#define g_El -49.0f			//mV
#define g_we (60*0.27f/10)	//pA
#define g_wi (-20*4.5f/10)	//pA
#define g_syn_delay 1.0f	//mS

unsigned long long int syn_ct = 0;

// Parameters used in CUBA benchmark
const struct neuron_params_t n_params = {
	.inv_tau_m = 1/20.0,		// [1/ms]
	.inv_tau_e = 1/5.0, 		// [1/ms]
	.inv_tau_i = 1/10.0,		// [1/ms]
	.De = 1.0/(1.0-(20.0/5.0)),	// Inverse of (1-tau_m/tau_e)
	.Di = 1.0/(1.0-(20.0/10.0)),// Inverse of (1-tau_m/tau_i)
	.reset_potential = -60.0,	// [mV]
	.threshold = -50.0,			// [mV]
	.refractory_period = 5.0	// [ms]
};

/* Called at the start of the spike handling. Updates the state of the neuron from the last update to now */
void bring_to_present(neuron_state_t* state, simtime_t now);
/* Compute the next time at which the neuron will fire, if any */
simtime_t getNextFireTime(neuron_state_t* state);
/* Are a and b within epsilon from each other? */
inline bool doubleEquals(double a, double b, double epsilon);
/* Compute the time after which the spike will take place given T0=0, V0 and I0, on a self-spiking neuron. */
double findSpikeDeltaBinaryBlind(struct neuron_helper_t* helper, double V0, double Ge0, double Gi0);

/* Compute Ge(delta_t) given Ge0, delta_t, and Tau_e */
inline double get_Ge_f(double delta_t, double Ge0);
extern double get_Ge_f(double delta_t, double Ge0);
/* Compute Gi(delta_t) given Gi0, delta_t, and Tau_i */
inline double get_Gi_f(double delta_t, double Gi0);
extern double get_Gi_f(double delta_t, double Gi0);

/* Compute V(delta_t) given V0, I0, delta_t, and the parameters*/
inline double get_V_t(struct neuron_helper_t *helper, double V0, double Ge0, double Ge1, double Gi0, double Gi1, double delta_t);
extern double get_V_t(struct neuron_helper_t *helper, double V0, double Ge0, double Ge1, double Gi0, double Gi1, double delta_t);

/* Initializes the topology of the benchmark */
void CUBATopology();

void printNeuronState(neuron_state_t* state);
/* Generates an array a s.t. a[i] contains the index of the first neuron of population i */
void gen_indexes(unsigned int * popsizes, unsigned int *out, int size);
/* Get population index from neuron ID */
int n2pop(unsigned long int neuron_ID);
/* Compute the time it takes for a self-spiking neuron to spike with I=0 and V0=Vreset */
double getSelfSpikeTime(struct neuron_helper_t *params);

/* Initialize a LIF neuron with exponential synapses */
neuron_state_t* InitExpLIFNeuron(unsigned long int me);

static struct neuron_helper_t n_helper_p[2];

/* NETWORK PARAMETERS v*/
//~ Probability connection table
double table[2][2]={{0.02,  0.02},\
					{0.02,  0.02}};
/* END NETWORK PARAMETERS ^ */

FILE* outFile = NULL;


/* Initialize the neuron */
void* NeuronInit(unsigned long int me){
	
	printdbg("[N%lu] INITIALIZING\n", me);

	return InitExpLIFNeuron(me);
}


neuron_state_t* InitExpLIFNeuron(unsigned long int me){
	neuron_state_t* state = malloc(sizeof(neuron_state_t));
	
	int pop = n2pop(me);

	struct neuron_helper_t *n_helper = &n_helper_p[pop];

	if(!n_helper->self_spike_time){
		//INIZIALIZZA TUTTO DENTRO ST_PARAMS
		n_helper->El = g_El;
		n_helper->Q = n_helper->El * n_params.inv_tau_m * n_params.inv_tau_m; //El/(tau_m^2)
		printdbg("El: %lf, Q: %lf\n", n_helper_p->El, n_helper_p->Q);
		
		n_helper->self_spike_time = getSelfSpikeTime(n_helper);
		printf("Neuron population %d self spike time: %lf\n", pop, n_helper->self_spike_time);
	}
	
	state->helper = n_helper;
	state->times_fired = 0;
	state->last_updated = 0;
	state->Ge = 0.0;
	state->Gi = 0.0;
	state->last_fired = -(n_params.refractory_period + 1);
	
	
	state->membrane_potential = n_params.reset_potential + Random()*(n_params.threshold - n_params.reset_potential);
	
	return state;
}

/* neuron receives spike */
void NeuronHandleSpike(unsigned long int me, simtime_t now, double value, void* neuron_state){
	printdbg("[N%lu] Spike received at %lf\n", me, now);
	
	neuron_state_t* state = (neuron_state_t*) neuron_state;
	
	simtime_t fireTime = -1.0;
	
	printdbg("[N%lu] Spike received with value: %lf at %lf\n", me, value, now);
	printNeuronState(state);
	
	// This brings the neuron to present. Computes and updates Gi, Ge and V.
	bring_to_present(state, now);
	// Apply the spike
	if (value > 0.0) {
		state->Ge += value;
	} else {
		state->Gi+=value;
	}
	
	state->last_updated = now; // Unnecessary, bring_to_present already updates it
	
	printdbg("[N%lu] updated values: ", me);
	printNeuronState(state);
	
	/* Now we have to compute the point in time - if any - at which this neuron will spike with the current values */
	fireTime = getNextFireTime(state);
	printdbg("[N%lu] Next spike time is %lf\n", me, fireTime);
	
	if (fireTime > -0.000001){ // Spike might happen in future
		MaybeSpikeAndWake(me, fireTime);
	}
	
	return;
}


/*  Neuron has been woken by the MaybeSpikeAndWake */
void NeuronWake(unsigned long int me, simtime_t now, void* n_state){
	printdbg("[N%lu] Woken up from spiking at %lf**********************************************************************!\n", me, now);
	

	neuron_state_t* state = n_state;
	
	bring_to_present(state, now);
	
	state->membrane_potential = n_params.reset_potential;
	state->times_fired++;
	state->last_fired = now;
	state->last_updated = now;
	printNeuronState(state);
	
	/* Now we have to compute the point in time - if any - at which this neuron will spike with the current values */
	double fireTime = -1.0;
	fireTime = getNextFireTime(state);
	printdbg("[N%lu] Next spike time is %lf\n", me, fireTime);
	
	if (fireTime > -0.000001){ // Spike might happen in future
		MaybeSpikeAndWake(me, fireTime);
	}
	return;
}


/* synapse handles the spike by updating its state and *returning* the spike intensity */
double SynapseHandleSpike(simtime_t now, unsigned long int src_neuron, unsigned long int dest_neuron, synapse_t* state){
	//~ printdbg("[Synapse] Spike directed from N %lu to N %lu with value: %lf\n", src_neuron, dest_neuron, state->weight);
	
	return state->weight;
	//~ float w = *((float*) state);
	//~ return w;
}


/* is the neuron done? */
bool NeuronCanEnd(unsigned long int me, neuron_state_t* state){	
	//~ if(state->times_fired >= SPIKE_MAX){
		//~ printdbg("[N%lu] Spikes: %lu, CanEnd!!!\n", me, state->times_fired);
	//~ }
	
	//~ return (state->times_fired >= SPIKE_MAX);
	return false;
}


void ProbeRead(simtime_t now, unsigned long int monitored_neuron, const neuron_state_t* neuron_state){
	//~ printf("[Probe N%lu] Neuron has fired %lu times at time %lf\n", monitored_neuron, neuron_state->times_fired, now);
	return;
}


void GatherStatistics(simtime_t now, unsigned long int neuron_id, const neuron_state_t* state){	
	//TODO: Fill this in to gather statistics and print them to file
	double spike_frequency = 1000 * state->times_fired/now;
	//~ double spike_freq_till_last = state->times_fired/state->last_updated;
	
	if(outFile == NULL){//First neuron to write. Open file and write
		outFile = fopen("CUBARun", "w");
	}
	
	fprintf(outFile, "N%lu, f%lf\n", neuron_id, spike_frequency);
	
}


void SNNInitTopology(unsigned long int neuron_count){
	CUBATopology(neuron_count);
}


// Creates the network for the CUBA benchmark.
void CUBATopology(unsigned long int neuron_count){
	// Number of excitatory neurons
	unsigned long int e_count = e_portion * neuron_count;
	// Inhibitory neurons
	unsigned long int i_count = neuron_count - e_count;
	
	unsigned long int pop_sizes[2] = {e_count, i_count};
	unsigned long int nn_cum[2] = {0, e_count};
	
	unsigned long int src_neuron = 0;
	unsigned long int dst_neuron = 0;
	
	synapse_t* synapse;
	
	double we = g_we;
	double wi = g_wi;
	
	double delay = g_syn_delay;
	
	for(int r = 0; r<2; r++){
		
		double wt = (r == 0 ? we : wi);
		
		for(int c=0; c<2; c++){
			double p = table[r][c];
			
			unsigned long int nsyns = pop_sizes[r]*pop_sizes[c]*p;
			printdbg("Synapses from %d to %d: %lu\n", r, c, nsyns);
			syn_ct += nsyns;
			
			while(nsyns--){
				src_neuron = (unsigned long int) (Random()*(pop_sizes[r])) + nn_cum[r];
				dst_neuron = (unsigned long int) (Random()*(pop_sizes[c])) + nn_cum[c];
				
				synapse = NewSynapse(src_neuron, dst_neuron, sizeof(synapse_t), true, delay);
				if(synapse==NULL){continue;} // Synapse is on remote node
				synapse->weight = wt;
			}
		}
	}
	NewProbe(0);
	
	printdbg("Total number of synapses: %llu\n", syn_ct);
	
	
}


/* Generates an array a s.t. a[i] contains the index of the first neuron of population i */
void gen_indexes(unsigned int * popsizes, unsigned int *out, int size){
	if(size != 2) return;
	
	out[0] = 0;
	out[1] = (unsigned int) e_portion * n_lps;
}


/* Get population index from neuron ID */
int n2pop(unsigned long int neuron_ID){
	return (int) neuron_ID > (unsigned long int) e_portion * n_lps;
}


/* Called at the start of the spike handling. Updates the state of the neuron from the last update to now.
 * We are certain that the neuron did not spike from the last time it was updated to now.
 * So we don't need to do calculations to account for that */
void bring_to_present(neuron_state_t* state, simtime_t now){
	struct neuron_helper_t *n_helper = state->helper;
	
	/*
	 * t0 = state->last_updated
	 * tf = now
	 * Ge0 = state->Ge
	 * Gi0 = state->Gi
	 * V0 = state->membrane_potential
	 * */
	
	double delta_t;
	double Gef;
	double Gif;
	
	simtime_t ref_end = state->last_fired + n_params.refractory_period; // When will refractory period end
	
	if(ref_end > state->last_updated){ // The neuron is in refractory period at t0
		printdbg("BTP: neuron is still in refractory period: last_updated: %lf, ref_end: %lf\n", state->last_updated, ref_end);
		
		state->membrane_potential = n_params.reset_potential; // Potential is clamped to reset potential
		
		if (ref_end >= now){ // Refractory lasts at least till now time. Only update the current and don't touch membrane potential
			printdbg("BTP: refractory period lasts further than now: %lf\n", now);
			printNeuronState(state);
			
			delta_t = now - state->last_updated;
			
			Gef = get_Ge_f(delta_t, state->Ge);
			Gif = get_Gi_f(delta_t, state->Gi);
			state->Ge = Gef;
			state->Gi = Gif;
			
			state->last_updated = now;
			printdbg("BTP: computed new I\n");
			printNeuronState(state);
			return;
		
		} else { // Skip the remaining refractory period by updating I to the value it will have at the end of it
			printdbg("BTP: handling refractory period\n");
			printNeuronState(state);
			
			delta_t = ref_end - state->last_updated;
			
			Gef = get_Ge_f(delta_t, state->Ge);
			Gif = get_Gi_f(delta_t, state->Gi);
			state->Ge = Gef;
			state->Gi = Gif;
			
			// Remember to update the time!
			state->last_updated = ref_end;
			printdbg("BTP: computed new I\n");
			printNeuronState(state);
			
		}
	}
	
	printdbg("BTP: bringing neuron values to present\n");
	// Now refractory is out of the way. Compute the new values for I and Vm
	delta_t = now - state->last_updated;
	
	printdbg("BTP: delta_t: %lf\n", delta_t);
	
	Gef = get_Ge_f(delta_t, state->Ge);
	Gif = get_Gi_f(delta_t, state->Gi);
	
	printdbg("BTP: Gef: %lf, Gif: %lf\n", Gef, Gif);
	
	// 7 arguments. Can we drop one?
	double V = get_V_t(n_helper, state->membrane_potential, state->Ge, Gef, state->Gi, Gif, delta_t);
	
	printdbg("BTP: Vf: %lf\n", V);
	
	state->membrane_potential = V;
	state->Ge = Gef;
	state->Gi = Gif;
	state->last_updated = now;
}


// Only takes into account neurons that are self spiking.
simtime_t getNextFireTime(neuron_state_t* state){
	struct neuron_helper_t* n_helper = state->helper;
	
	// Down here only non-self-spiking neurons are considered.
	//~ double I0 = state->I;
	double Ge0 = state->Ge;
	double Gi0 = state->Gi;
	
	simtime_t t0 = state->last_updated;
	double V0 = state->membrane_potential;
	
	printNeuronState(state);
	
	double delta_t;
	double ref_end = state->last_fired + n_params.refractory_period;
	
	printdbg("getFT: Ref_end: %lf\n", ref_end);
	
	/* Now check if we are still in refractory period */
	if(ref_end > state->last_updated){
		
		t0 = ref_end; // When will refractory period end, new t0
		
		delta_t = t0 - state->last_updated;
		
		Ge0 = get_Ge_f(delta_t, state->Ge);
		Gi0 = get_Gi_f(delta_t, state->Gi);
		
		V0 = n_params.reset_potential; // Potential is clamped to reset potential
		printdbg("getFT: Ref period out of the way\n");
		printdbg("New V0: %lf, Ge0: %lf, Gi0: %lf, t0: %lf\n", V0, Ge0, Gi0, t0);
	}
	/* Refractory done */
	
	/* If El>Vth the neuron is self spiking and we need to
	 * check when it will spike regardless of the value of I*/
	if(n_helper->El > n_params.threshold){ // Neuron is self spiking
		delta_t = findSpikeDeltaBinaryBlind(n_helper, V0, Ge0, Gi0);
		if (delta_t<=0){
			printf("Something went wrong\n");
			abort();
		}
		return t0 + delta_t;
	}
	
	// THIS IS UNREACHABLE IN THIS IMPLEMENTATION
	printf("WARNING: UNREACHABLE CODE REACHED\n");
	abort();
	return -1;
	
	// However, to handle all cases this might need to be handled in the future.
	// Definitely can be done in a fashion similar to that using Icond, but taking into
	// 	account Ge+Gi rather than I.
}


inline double get_Ge_f(double delta_t, double Ge_i){
	return exp(-delta_t * n_params.inv_tau_e) * Ge_i;
}

inline double get_Gi_f(double delta_t, double Gi_i){
	return exp(-delta_t * n_params.inv_tau_i) * Gi_i;
}


// Might want to make another function that computes everything in one, given the state.
inline double get_V_t(struct neuron_helper_t *helper, double V0, double Ge0, double Ge1, double Gi0, double Gi1, double delta_t){
	/* Now compute V(delta_t) and return it */
	double aux = (V0 - Ge0*n_params.De - Gi0*n_params.Di - helper->Q) * exp(-delta_t * n_params.inv_tau_m);
	
	double expe = exp(-delta_t * n_params.inv_tau_e);
	double expi = exp(-delta_t * n_params.inv_tau_i);
	
	return (Ge0*expe*n_params.De + Gi0*expi*n_params.Di + helper->Q + aux);
}


/* Given a neuron, compute how long it takes to spike without external time-dependent inputs */
double getSelfSpikeTime(struct neuron_helper_t *n_helper){	
	// Check that dV/dt > 0 for V->Vth and ge=gi=0
	if(n_params.threshold > n_helper->El) return -1;
	
	double delta_t = 0;
	
	//~ double It;
	double Vt;
	double Ge_t;
	double Gi_t;
	double tmin = 0;
	double tmax = 50;
	double step = 50;
	
	// While Vt is below Vth
	while(tmin < MODEL_MAX_SIMTIME*50){// !isinf(delta_t)) {
		/* Now compute V(delta_t) */
		Vt = get_V_t(n_helper, n_params.reset_potential, 0, 0, 0, 0, tmax);
		
		if (Vt >= n_params.threshold){
			//~ printdbg("delta_t: %lf, Vt %lf > Vth %lf\n", tmax, Vt, n_params.threshold);
			break;
		} else {
			//~ printdbg("delta_t: %lf, Vt %lf < Vth %lf\n", tmax, Vt, n_params.threshold);
			tmin = tmax;
			//~ tmax = tmax + step;
			tmax *=2;
		}
	} // We now have a cap for the spike time
	
	// Might want to check this approach, could lead to errors with short simtimes.
	if(tmin > MODEL_MAX_SIMTIME*50){//isinf(delta_t)){
		return -1;
	}
	
	// While V is farther than epsilon from Vth. Binary search
	do {
		//~ printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		delta_t = (tmax + tmin)/2;
		/* Now compute V(delta_t) */
		Vt = get_V_t(n_helper, n_params.reset_potential, 0, 0, 0, 0, delta_t);
		
		
		if (Vt > n_params.threshold){
			//~ printdbg("Vt %lf > Vth %lf\n", Vt, n_params.threshold);
			tmax = delta_t;
		} else {
			//~ printdbg("Vt %lf < Vth %lf\n", Vt, n_params.threshold);
			tmin = delta_t;
		}
	} while(!doubleEquals(Vt, n_params.threshold, V_TOLERANCE) || !doubleEquals(tmin, tmax, T_TOLERANCE));
	
	printdbg("SelfSpikeTime spikedelta: %lf\n", delta_t);
	return delta_t;
}


/* Compute the time after which the spike will take place given T0=0 and V0
 * on a self spiking neuron */
double findSpikeDeltaBinaryBlind(struct neuron_helper_t* helper, double V0, double Ge0, double Gi0){
	double delta_t = 0;
	
	double Vt;
	double Ge_t;
	double Gi_t;
	double tmin = 0;
	double tmax = helper->self_spike_time;
	double step = helper->self_spike_time;
	
	if (step <= 0) {
		printf("FSDBB: SOMETHING WRONG HAPPENED\n");
		fflush(stdout);
		abort();
	}

	// While Vt is below Vth
	printdbg("FSDBB: Neuron is self spiking\n");
	while(true) {
		delta_t = tmax;
		/* Compute Ge(delta_t) and Gi(delta_t) */
		Ge_t = get_Ge_f(delta_t, Ge0);
		Gi_t = get_Gi_f(delta_t, Gi0);
		/* Now compute V(delta_t) */
		Vt = get_V_t(helper, V0, Ge0, Ge_t, Gi0, Gi_t, delta_t);
		
		if (Vt >= n_params.threshold){
			//~ printdbg("FSDBB: delta_t: %lf, Vt %lf > Vth %lf\n", delta_t, Vt, n_params.threshold);
			break;
		} else {
			//~ printdbg("FSDBB: delta_t: %lf, Vt %lf < Vth %lf\n", delta_t, Vt, n_params.threshold);
			tmin = tmax;
			tmax = tmax + step;
			//~ tmax *=2;
		}
	}
	
	// While V is farther than epsilon from Vth
	do {
		//~ printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		delta_t = (tmax + tmin)/2;
		
		/* Compute Ge(delta_t) and Gi(delta_t) */
		Ge_t = get_Ge_f(delta_t, Ge0);
		Gi_t = get_Gi_f(delta_t, Gi0);
		/* Now compute V(delta_t) */
		Vt = get_V_t(helper, V0, Ge0, Ge_t, Gi0, Gi_t, delta_t);
		
		if (Vt > n_params.threshold){
			//~ printdbg("Vt %lf > Vth %lf\n", Vt, n_params.threshold);
			tmax = delta_t;
		} else {
			//~ printdbg("Vt %lf < Vth %lf\n", Vt, n_params.threshold);
			tmin = delta_t;
		}
	} while(!doubleEquals(Vt, n_params.threshold, V_TOLERANCE) && !doubleEquals(tmin, tmax, T_TOLERANCE));
	
	printdbg("FSDBB spikedelta: %lf\n", delta_t);
	return delta_t;
}


inline bool doubleEquals(double a, double b, double epsilon){
	return fabs(a-b) < epsilon;
}
extern bool doubleEquals(double a, double b, double epsilon);


void printNeuronState(neuron_state_t* state){
	printdbg("V: %lf, Ge: %lf, Gi: %lf, last_up: %lf, last_spike: %lf\n", state->membrane_potential, state->Ge, state->Gi, state->last_updated, state->last_fired);
}