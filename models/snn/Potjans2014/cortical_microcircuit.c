#include "cortical_microcircuit.h"

//~ #define M_DEBUG                // Undefine to decrease verbosity

#ifdef M_DEBUG

# define printdbg(...) \
    printf(__VA_ARGS__)
#else
# define printdbg(...)
#endif

#undef M_DEBUG

#define MODEL_MAX_SIMTIME global_config.termination_time

#define V_TOLERANCE 0.01 //mV
#define T_TOLERANCE 0.01 //mS

#define doubleEquals(a, b, epsilon) (fabs((a)-(b)) < (epsilon))

__thread unsigned long long int syn_ct = 0;

// Parameters used in Potjans2014
static const struct neuron_params_t n_params = {
	.threshold = -50.0,
	.reset_potential = -65.0,
	.inv_tau_m = 1/10.0,
	.inv_tau_syn = 1/0.5,
	.inv_C_m = 1/250.0,
	.refractory_period = 2.0
};
// Synapses parameters
#define g_d_ex 			1.5 //* 0.001  			//# [ms] Excitatory delay
#define g_std_d_ex 		0.75 //* 0.001			//# [ms] Std. Excitatory delay
#define g_d_in			0.80 //* 0.001 			//# [ms] Inhibitory delay
#define g_std_d_in		0.4 //* 0.001			//# [ms] Std. Inhibitory delay
#define g_tau_syn		0.5 //* 0.001			//# [ms] Post-synaptic current time constant

#define g_w_ex			87.8 //* 0.000000001	//# [pA] Standard excitatory weight
#define g_std_w_ex		0.1*g_w_ex				//# standard deviation weight
#define in_g			4.0						//# Multiplicative factor for the inhibitory weights

// External current
#define g_I_ext			0.002 //* 0.001				//# [pA] External current

/* Called at the start of the spike handling. Updates the state of the neuron from the last update to now */
void bring_to_present(neuron_state_t* state, simtime_t delta_upd, simtime_t delta_spike);
/* Compute the next time at which the neuron will fire, if any */
simtime_t getNextFireTime(neuron_state_t* state, simtime_t delta_spike);
/* The neuron will spike. Find out how long it will take */
double findSpikeDeltaNewton(struct neuron_helper_t *params, double V0, double I0, double time_interval);
double findSpikeDeltaBinary(struct neuron_helper_t *params, double V0, double I0, double time_interval);
/* Compute the time after which the spike will take place given T0=0, V0 and I0, on a self-spiking neuron. */
double findSpikeDeltaBinaryBlind(struct neuron_helper_t *params, double V0, double I0);
/* Compute I(delta_t) given I0, delta_t, and Tau_s */
inline double get_I_f(double delta_t, double I0);
extern double get_I_f(double delta_t, double I0);

/* Compute V(delta_t) given V0, I0, delta_t, and the parameters*/
inline double get_V_t(struct neuron_helper_t *params, double V0, double I0, double I1, double delta_t);
extern double get_V_t(struct neuron_helper_t *params, double V0, double I0, double I1, double delta_t);

/* Initializes the topology of Potjans2014's Cortical microcircuit */
void MicrocircuitTopology();

void printNeuronState(neuron_state_t* state);
/* Generates an array a s.t. a[i] contains the index of the first neuron of population i */
void gen_indexes(unsigned int * popsizes, unsigned int *out, int size);
/* Computes and sets value of global var bg_layer */
void get_bg_layer(int bg_type);
/* Get population index from neuron ID */
int n2pop(unsigned long int neuron_ID);
/* Compute the time it takes for a self-spiking neuron to spike with I=0 and V0=Vreset */
double getSelfSpikeTime(struct neuron_helper_t *params);

/* Initialize a LIF neuron */
neuron_state_t* InitLIFNeuron(unsigned long int me);
/* Wake a Poisson Neuron */
void PoissonWake(unsigned long int me, simtime_t now);

static struct neuron_helper_t n_helper_p[9];

/* NETWORK PARAMETERS v*/
// Population size per layer
//          				2/3e   2/3i   4e    4i    5e    5i    6e     6i    Th
unsigned int n_layer[9] = {20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948, 902};
//Tot neurons: 78071 with thal=1, 77169 otherwise

//~ Probability connection table
double table[8][9]={{0.101,  0.169, 0.044, 0.082, 0.032, 0.,     0.008, 0.,     0.    },\
		{0.135,  0.137, 0.032, 0.052, 0.075, 0.,     0.004, 0.,     0.    },\
		{0.008,  0.006, 0.050, 0.135, 0.007, 0.0003, 0.045, 0.,     0.0983},\
		{0.069,  0.003, 0.079, 0.160, 0.003, 0.,     0.106, 0.,     0.0619},\
		{0.100,  0.062, 0.051, 0.006, 0.083, 0.373,  0.020, 0.,     0.    },\
		{0.055,  0.027, 0.026, 0.002, 0.060, 0.316,  0.009, 0.,     0.    },\
		{0.016,  0.007, 0.021, 0.017, 0.057, 0.020,  0.040, 0.225,  0.0512},\
		{0.036,  0.001, 0.003, 0.001, 0.028, 0.008,  0.066, 0.144,  0.0196}};
					
// Layer-specific background input
unsigned int bg_layer_specific[] = {1600, 1500 ,2100, 1900, 2000, 1900, 2900, 2100};

// Layer-independent background input
unsigned int bg_layer_independent[] = {2000, 1850 ,2000, 1850, 2000, 1850, 2000, 1850};
/* END NETWORK PARAMETERS ^ */

unsigned int protocol = 0;
unsigned int bg_type = 0;
unsigned int stim = 0;
unsigned int thal = 0;
// Init this in argp
unsigned int* bg_layer;

FILE* outFile = NULL;

static const double exp_adj[] = {
	1.036199565629158, 1.027277981818330, 1.019341151369058,
	1.012316468833528, 1.006139353234098, 1.000751255315155,
	0.996100136699148, 0.992140185571610, 0.988830196962552,
	0.986132544047804, 0.984013562521051, 0.982443170190283,
	0.981394102987059, 0.980841516073624, 0.980763042138659,
	0.981138430277712, 0.981957925989939, 0.983178546312421,
	0.984811614558841, 0.986834419454900, 0.989234236856811,
	0.991999812972351, 0.995121060808105, 0.998588377827231,
	1.002392956865277, 1.006527340270889, 1.010984919213803,
	1.015759089560808, 1.020843804511367, 1.026234345941742,
	1.031901110282296, 1.037916082604387
};

unsigned int v0type = 0;
double original_v0_mean[] = {-58.0, -58.0, -58.0, -58.0, -58.0, -58.0, -58.0, -58.0};
double original_v0_std[] = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
double optimized_v0_mean[] = {-68.28, -63.16, -63.33, -63.45, -63.11, -61.66, -66.72, -61.43};
double optimized_v0_std[] = {5.36, 4.57, 4.74, 4.94, 4.94, 4.55, 5.46, 4.48};
double* v0_mean = optimized_v0_mean;
double* v0_std = optimized_v0_std;


enum{
	OPT_PROTOCOL = 129,
	OPT_V0_TYPE = 130
};

struct ap_option model_options[] = {
	{"protocol", OPT_PROTOCOL, "VALUE",
		"Protocol to run the simulation with.\n\
\tprotocol = 0: spontaneous activity (figure 2)\n\
\tprotocol = 1: DC input and layer-independent experiments\n\
\tprotocol = 5: response to transient thalamic input\n"
	},
	{"v0", OPT_V0_TYPE, "VALUE",
		"Initial conditions for membrane potential.\n\
\tv0 = 0: Optimized - population-specific mean and standard deviation, allowing\
a reduction of the initial activity burst in the network (default)\n\
\tv0 = 1: Original - uniform mean and standard deviation for all populations as\
used in earlier implementations of the model\n"
	},
	{0}
};

void model_parse (int key, const char *arg){
	switch (key) {
		case AP_KEY_FINI:
			get_bg_layer(bg_type);
			break;
		case AP_KEY_INIT:
			break;
		case OPT_PROTOCOL:
		{
			if(sscanf(arg, "%u", &protocol) != 1) {
				printf("Could not parse protocol option\n");
				abort();
			}
			switch(protocol){
				case 0:
				{
					bg_type = 0;
					stim = 0;
					//TODO: fix poissonInput to make stim=0 viable
					printf("Protocol 0 not supported\n");
					abort();
					break;
				}
				case 1:
				{
					bg_type = 0;
					stim = 1;
					break;
				}
				case 5:
				{
					bg_type = 0;
					stim = 1;
					thal = 1;
					break;
				}
				default:
				{
					printf("Protocol not recognized\n");
					abort();
				}
			}
			
			break;
			
		}
		case OPT_V0_TYPE:
		{
			if(sscanf(arg, "%u", &v0type) != 1) {
				printf("Could not parse v0 option\n");
				abort();
			}
			switch(v0type){
				case 0:
				{
					// Optimized
					v0_mean = optimized_v0_mean;
					v0_std = optimized_v0_std;
					break;
				}
				case 1:
				{
					// Original
					v0_mean = original_v0_mean;
					v0_std = original_v0_std;
					break;
				}
				default:
				{
					printf("Initial membrane potential (V0) not recognized\n");
					abort();
				}
			}

			break;

		}
		default:
		{
			printf("Argument not recognized\n");
			abort();
		}
	}
}

static inline double csexp(double x)
{
	int64_t tmp = 1512775L * x + 1072632447L;
	int i = (tmp >> 15) & ((1 << 5) - 1);
	tmp <<= 32;
	double r;
	memcpy(&r, &tmp, sizeof(r));
	return r * exp_adj[i];
}

/* Initialize the neuron */
void* NeuronInit(unsigned long int me){
	
	printdbg("[N%lu] INITIALIZING\n", me);

	int pop = n2pop(me);
	if(pop < 8){
		return InitLIFNeuron(me);
	} else { // Pop = 8. i.e. Thalamic
		return NULL;
	}
}

neuron_state_t* InitLIFNeuron(unsigned long int me){
	neuron_state_t* state = malloc(sizeof(neuron_state_t));
	
	int pop = n2pop(me);

	struct neuron_helper_t *n_helper = &n_helper_p[pop];

	if(!n_helper->self_spike_time){
		//INIZIALIZZA TUTTO DENTRO ST_PARAMS
		if(stim==1){
			n_helper->Iext = (0.3512*bg_layer[pop]);
			printdbg("Iext: %lf\n", n_helper->Iext);
		} else {
			n_helper->Iext = 0;
		}
		
		n_helper->A0 = n_params.inv_C_m/(n_params.inv_tau_m - n_params.inv_tau_syn);
		n_helper->A2 = n_params.reset_potential + ((n_helper->Iext * n_params.inv_C_m)/n_params.inv_tau_m);// * 0.001;
		printdbg("A0: %lf, A2: %lf\n", n_helper_p->A0, n_helper_p->A2);
		
		n_helper->Icond = ((n_params.threshold - n_params.reset_potential) * (n_params.inv_tau_m)) / n_params.inv_C_m - (n_helper->Iext);
		printdbg("Neuron population %d Icond: %lf\n", pop, n_helper->Icond);
		n_helper->self_spike_time = getSelfSpikeTime(n_helper);
		printdbg("Neuron population %d self spike time: %lf\n", pop, n_helper->self_spike_time);
	}
	
	state->helper = n_helper;
	state->last_updated = 0;
	state->I = 0.0;
	state->last_fired = -(n_params.refractory_period + 1);
	
	state->membrane_potential = v0_mean[pop] + v0_std[pop] * Normal();
	if(state->membrane_potential > n_params.threshold){
		state->membrane_potential = n_params.threshold + 0.01;
	}
	return state;
}

/* neuron receives spike */
void NeuronHandleSpike(unsigned long int me, simtime_t now, double value, void* neuron_state){
	printdbg("[N%lu] Spike received at %lf\n", me, now);
	
	if(me >= 77169){
		printdbg("[N%lu] It's Poisson: waking\n", me);
		PoissonWake(me, now);
		return;
	}
	
	neuron_state_t* state = (neuron_state_t*) neuron_state;
	
	simtime_t fireTime;
	
	printdbg("[N%lu] Spike received with value: %lf at %lf\n", me, value, now);

	simtime_t delta_spike = now - state->last_fired;
	// This brings the neuron to present. Computes and updates I and V.
	bring_to_present(state, now - state->last_updated, delta_spike);
	// Apply the spike
	state->I += value;
	state->last_updated = now;
	
	printdbg("[N%lu] updated values: ", me);
	printNeuronState(state);
	
	/* Now we have to compute the point in time - if any - at which this neuron will spike with the current values */
	fireTime = getNextFireTime(state, delta_spike);
	printdbg("[N%lu] Next spike time is %lf\n", me, fireTime);
	
	if (fireTime > 0.0){ // Spike might happen in future
		MaybeSpikeAndWake(me, fireTime);
	}
	
	return;
}

#define THAL_RATE 120
#define THAL_START_DELAY 700
#define THAL_CYCLE 1000
#define THAL_ACTIVE_PERIOD 10

void PoissonWake(unsigned long int me, simtime_t now){
	printdbg("[Poisson %lu] Woken up from spiking at %lf**********************************************************************!\n", me, now);
	
	simtime_t delay = Expent(1000/THAL_RATE);
	simtime_t nextSpike = now;
	if(now == 0) {nextSpike += THAL_START_DELAY;}
	nextSpike += delay;

	simtime_t tmod = fmod(nextSpike - THAL_START_DELAY, THAL_CYCLE);
	unsigned jumps = delay / THAL_ACTIVE_PERIOD + (tmod > THAL_ACTIVE_PERIOD);
	
	nextSpike += jumps * (THAL_CYCLE - THAL_ACTIVE_PERIOD);
	MaybeSpikeAndWake(me, nextSpike);
}

/*  Neuron has been woken by the MaybeSpikeAndWake */
void NeuronWake(unsigned long int me, simtime_t now, void* n_state){
	printdbg("[N%lu] Woken up from spiking at %lf**********************************************************************!\n", me, now);
	
	if(me >= 77169){
		//~ printf("[N%lu] Spiker is a Poisson Neuron*****!\n", me);
		PoissonWake(me, now);
		return;
	}

	neuron_state_t* state = n_state;
	simtime_t t_since_last_eval = now - state->last_updated;
	state->I = get_I_f(t_since_last_eval, state->I);
	state->membrane_potential = n_params.reset_potential;
//	state->times_fired++;
	state->last_fired = now;
	state->last_updated = now;
	printNeuronState(state);
	
	/* Now we have to compute the point in time - if any - at which this neuron will spike with the current values */
	double fireTime = getNextFireTime(state, 0);
	printdbg("[N%lu] Next spike time is %lf\n", me, fireTime);
	
	if (fireTime > 0.0){ // Spike might happen in future
		MaybeSpikeAndWake(me, fireTime);
	}
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
	//~ return (state->times_fired >= SPIKE_MAX);
	return false;
}

void ProbeRead(simtime_t now, unsigned long int monitored_neuron, const neuron_state_t* neuron_state){
	//~ printf("[Probe N%lu] Neuron has fired %lu times at time %lf\n", monitored_neuron, neuron_state->times_fired, now);
	return;
}

void GatherStatistics(simtime_t now, unsigned long int neuron_id, const neuron_state_t* state){
	if(neuron_id >= 77169) {return;}//Poisson neuron

	if(outFile == NULL){//First neuron to write. Open file and write
		char* fname = malloc(strlen("Potjans2014Run_") + 10);
		sprintf(fname, "Potjans2014Run_%lu", neuron_id);
//		sprintf(fname + strlen("Potjans2014Run_"), "%lu", neuron_id);
		outFile = fopen(fname, "w");
		free(fname);
	}
}

void SNNInitTopology(unsigned long int neuron_count){
	MicrocircuitTopology();
}

/* Initializes the topology */
void MicrocircuitTopology(){//(stim, bg_type, w_ex, g, bg_freq, nsyn_type, thal)
	int nsyn_type = 0;

	// Index of the first neuron for every population
	unsigned int nn_cum[9];
	gen_indexes(n_layer, nn_cum, 9);
	
	//###########################################################################
	//# Connecting neurons
	//###########################################################################
	
	unsigned long int src_neuron;
	unsigned long int dst_neuron;
	synapse_t* synapse;
	double delay;
	double wt; // weight
	int nsyn = 0;
	
	// Connection between neuron populations
	for(int c=0; c<8; c++){
		for(int r=0; r<8; r++){

			if (nsyn_type==0){
				// number of synapses calculated with equation 3 from the article
				nsyn = log(1.0-table[r][c])/log(1.0 - (1.0/(n_layer[c]*n_layer[r])));
				printdbg("NSYN from %d to %d: %d\n", c, r, nsyn);
			} else if (nsyn_type==1) {
				// number of synapses calculated with equation 5 from the article
				nsyn = table[r][c]*n_layer[c]*n_layer[r];
				printdbg("NSYN from %d to %d: %d\n", c, r, nsyn);
			}

			syn_ct += nsyn;

			while(nsyn--){
				src_neuron = (unsigned long int) (Random()*(n_layer[c])) + nn_cum[c];
				dst_neuron = (unsigned long int) (Random()*(n_layer[r])) + nn_cum[r];
				
				if ((c % 2) == 0){ // Excitatory connections
					
					delay = g_d_ex + g_std_d_ex*Normal();
					if(delay < 0.1){delay=0.1;}
					// Make sure to always make the same number of calls to RNG regardless of conditions
					wt = g_w_ex + g_std_w_ex*Normal();
					synapse = NewSynapse(src_neuron, dst_neuron, sizeof(synapse_t), true, delay);
					if(synapse==NULL){continue;}
					synapse->weight = wt;
					if(synapse->weight < 0.0){synapse->weight = 0.0;}
					// Synaptic weight from L4e to L2/3e is doubled
					if (c == 2 && r == 0){
						synapse->weight *= 2.0;
					}
					
				} else { // Inhibitory connections
					
					delay = g_d_in + g_std_d_in*Normal();
					if(delay < 0.1){delay=0.1;}
					wt = -(g_w_ex + g_std_w_ex*Normal());
					synapse = NewSynapse(src_neuron, dst_neuron, sizeof(synapse_t), true, delay);
					if(synapse==NULL){continue;}
					synapse->weight = wt;
					if(synapse->weight > 0.0){synapse->weight = 0.0;}
					synapse->weight *= in_g;
					
				}
			}

			//printdbg("TotalSynapses %llu\n", syn_ct);
		}
	}
	
	if(thal==1){
		for (int r = 0; r<8; r++) {
			nsyn = (table[r][8]*n_layer[8]*n_layer[r]);
			printdbg("NSYN from 8 to %d: %d\n", r, nsyn);
			
			for(int i=0; i<nsyn; i++){
				src_neuron = (unsigned long int) (Random()*(n_layer[8])) + nn_cum[8];
				dst_neuron = (unsigned long int) (Random()*(n_layer[r])) + nn_cum[r];
				delay = 0.1;
				if(delay < 0.1){delay=0.1;}
				wt = g_w_ex + (g_w_ex*0.1*Normal());
				synapse = NewSynapse(src_neuron, dst_neuron, sizeof(synapse_t), true, delay);
				if(synapse==NULL){continue;}
				synapse->weight = wt;
				if(synapse->weight < 0.0){synapse->weight = 0.0;}
			}
		}
		syn_ct += nsyn;
		//printdbg("TotalSynapses %llu\n", syn_ct);
	}

	if(bg_type==2){free(bg_layer);}
}

/* Generates an array a s.t. a[i] contains the index of the first neuron of population i */
void gen_indexes(unsigned int * popsizes, unsigned int *out, int size){
	if(size <= 0) return;
	
	out[0] = 0;
	for (int i = 1; i < size; i++){
		out[i] = popsizes[i-1] + out[i-1];
	}
}

/* Get population index from neuron ID */
int n2pop(unsigned long int neuron_ID){
	for(int i=0; i<9; i++){
		if(neuron_ID < n_layer[i]){
			return i;
		}
		neuron_ID -= n_layer[i];
	}
	return -1;
}

void get_bg_layer(int bg_type){
	//# Background number per layer
	if (bg_type == 0) {
		//# layer-specific:
		bg_layer = bg_layer_specific;
	} else if (bg_type == 1) {
		//# layer-independent:
		bg_layer = bg_layer_independent;
	} else if (bg_type == 2) {
		bg_layer = malloc(sizeof(unsigned int) * 8);
		//#layer-independent-random:
		for (int i=0; i<8; i+=2) {
			//# range of the number of inputs given to an excitatory population:
			int exc_bound_A = bg_layer_specific[i];
			int exc_bound_B = bg_layer_independent[i];
			unsigned int diff_exc = abs(exc_bound_A-exc_bound_B);
			double exc_input;
			
			//# range of the number of inputs given to an inhibitory population:
			double inh_bound_A;
			double inh_bound_B;
			double diff_inh;
			double inh_input;
			
			//# randomly choosing a number for the external input to an excitatory population:
			if (exc_bound_A<=exc_bound_B) {
				exc_input = exc_bound_A + Random()*diff_exc;
			} else {//if (exc_bound_A>exc_bound_B){
				exc_input = exc_bound_B + Random()*diff_exc;
			}

			//# range of the number of inputs given to an inhibitory population:
			if (i!=6) {
				//# eq. 4 from the article
				inh_bound_A = ((1-0.1)/(1+0.1))*exc_input;
			} else {
				// # eq. 4 from the article
				inh_bound_A = ((1-0.2)/(1+0.2))*exc_input;
			}
			inh_bound_B = exc_input;
			diff_inh = fabs(inh_bound_A-inh_bound_B);

			//# randomly choosing a number for the external input to an inhibitory population:
			if (inh_bound_A<=inh_bound_B) {
				inh_input = inh_bound_A + Random()*diff_inh;
			} else {
				inh_input = inh_bound_B + Random()*diff_inh;
			}

			//# array created to save the values:
			bg_layer[i] = (unsigned int) exc_input;
			bg_layer[i+1] = (unsigned int) inh_input;
		}
	}
}


//OK
/* Called at the start of the spike handling. Updates the state of the neuron from the last update to now.
 * We are certain that the neuron did not spike from the last time it was updated to now.
 * So we don't need to do calculations to account for that */
void bring_to_present(neuron_state_t* state, simtime_t delta_update, simtime_t delta_spike){
	/*
	 * delta_update = now - state->last_updated
	 * delta_spike = now - state->last_spiked
	 * */

	struct neuron_helper_t *n_helper = state->helper;

	double If = get_I_f(delta_update, state->I);
	double Ient = state->I;
	state->I = If;

	if (delta_spike < n_params.refractory_period)
		return;

	// Amount of time of refractory period that we still have to take into account
	double remaining = delta_update + n_params.refractory_period - delta_spike;

	if (remaining > 0.0) {// computes current at end of refractory period
		Ient = get_I_f(remaining, state->I);
		delta_update -= remaining;
	}

	// Update membrane potential from end of refractory to now
	state->membrane_potential = get_V_t(state->helper, state->membrane_potential, Ient, If, delta_update);
}

simtime_t getNextFireTime(neuron_state_t* state, simtime_t delta_spike){
	struct neuron_helper_t* n_helper = state->helper;
	
	/* Necessary condition for the membrane to reach threshold potential:
	 * dV/dt>0 with V->Vth
	 * Indeed, as potential approaches the threshold value, the leakage will be stronger.
	 * We need an I that can overcome it */
	double Icond = state->helper->Icond;
	
	printdbg("getFT: Icond: %lf\n", Icond);	

	double I0 = state->I;
	simtime_t t0 = state->last_updated;
	double V0 = state->membrane_potential;
	
//	printNeuronState(state);
	
	/* Check that the condition for spiking is met.
	 * Since I decays towards 0 with the passing of time (regardless of whether I is positive or negative),
	 * if the spiking condition is not met at the start, it will never be met.
	 * Unless a spike is received. But that is a future event.*/
//	TODO: check that the changes do not break any assumption made for non-self-spiking neurons,
//	 then uncomment this and remove the abort below. The code should be ok though
//	if (n_helper->A2 < n_params.threshold && I0 < Icond){// Neuron is not self-spiking and I<Icond
//		printdbg("getFT: I %lf < Icond %lf. Will never spike\n", I0, Icond);
//		return -1; // Can never spike this way. Return
//	}

	// Time until end of refractory period
	simtime_t delta_t;
	/* Now check if we are still in refractory period */
	if(n_params.refractory_period > delta_spike){

		delta_t = n_params.refractory_period - delta_spike;
		// Current at end of refractory period
		I0 = get_I_f(delta_t, state->I);
	} else {
		delta_t = 0;
	}
	/* Refractory done */
	
	/* If A2>Vth the neuron is self spiking and we need to
	 * check when it will spike regardless of the value of I*/
	if(n_helper->A2 > n_params.threshold){ // Neuron is self spiking
		return t0 + findSpikeDeltaBinaryBlind(n_helper, V0, I0) + delta_t;
	}

	abort(); // TODO: check that the changes do not break any assumption made below, then remove the abort

	// Neuron is not self spiking
	if (I0 < Icond){
		printdbg("getFT: I %lf < Icond %lf. Will never spike\n", I0, Icond);
		return -1; // Can never spike this way. Return
	}
	
	/* Now let's find the t_upper(/delta_t), the time for(/after) which I = Icond.
	 * t_upper is an upper bound to the spike time. We have two cases:
	 * a) V(t_upper) >= Vth. The neuron spikes in a time between t0 and t_upper.
	 * b) V(t_upper) < Vth. The neuron did not spike and will not in the future (since from t_upper on, I<Icond)
	 * These are the only two cases: since Icond is needed to surpass Vth,
	 * any I > Icond is sufficient to keep V > Vth.
	 * As such, if there is a time t' in [t0, t_upper] s.t. V(t') > Vth,
	 * this will stay true for at least as long as I>Icond. i.e. for all t in [t', t_upper].
	 * Thus, if V(t_upper) >= Vth, we crossed Vth once in [t0,t_upper]. Otherwise, we never did. QED
	 * */
	delta_t = -(log(Icond/I0) / n_params.inv_tau_syn); //Divide by the inverse
	printdbg("getFT: delta_t for Icond: %lf\n", delta_t);
	printdbg("getFT: I at delta_t: %lf\n", get_I_f(delta_t, I0));
	
	/* Now compute V(t_upper) = V(t0+delta_t) */
	double V = get_V_t(n_helper, V0, I0, Icond, delta_t);
	
	/* Check which case we are in */
	if(V < n_params.threshold){ // Case b: neuron does not spike
		printdbg("V %lf < Vth %lf. Will not spike\n", V, n_params.threshold);
		return -1;
	}
	// Case a: V>=Vth. Neuron spiked at ts in [t0, t_upper]
	/* Here if we have V=Vth AND dV/dt = 0, we just touched Vth as a max. Check if this is the case */
	double dvdt = (-V + n_params.reset_potential)*(n_params.inv_tau_m) + (Icond + n_helper->Iext)*(n_params.inv_C_m);
	if(doubleEquals(V, n_params.threshold, V_TOLERANCE) && doubleEquals(dvdt, 0.0, V_TOLERANCE)){
		return t0 + delta_t;
	}
	
	/* Otherwise */
	return t0 + findSpikeDeltaBinary(n_helper, V0, I0, delta_t);
}

/* Compute the time after which the spike will take place given T0=0, V0 and I0, on a NON-self-spiking neuron. */
double findSpikeDeltaBinary(struct neuron_helper_t* n_helper, double V0, double I0, double time_interval){
	double It;
	double Vt;
	double tmin = 0;
	double tmax = time_interval;
	
	// While V is farther than epsilon from Vth
	do {
		printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		double delta_t = (tmax + tmin)/2;
		/* Compute I(delta_t) */
		It = get_I_f(delta_t, I0);
		/* Now compute V(delta_t) */
		Vt = get_V_t(n_helper, V0, I0, It, delta_t);
		
		if (Vt > n_params.threshold){
			printdbg("Vt %lf > Vth %lf\n", Vt, n_params.threshold);
			tmax = delta_t;
		} else {
			printdbg("Vt %lf < Vth %lf\n", Vt, n_params.threshold);
			tmin = delta_t;
		}
		if (tmax <= tmin + T_TOLERANCE){
			printdbg("BinarySearch spikedelta: %lf\n", delta_t);
			return  delta_t;
		}
	} while(1);
}

inline double get_I_f(double delta_t, double I_i){
	return csexp(-delta_t * n_params.inv_tau_syn) * I_i;
}

inline double get_V_t_acc(struct neuron_helper_t *helper, double V0, double I0, double I1, double delta_t){
	/* Now compute V(delta_t) and return it */
	double aux = (V0 - helper->A2 - I0 * helper->A0) * exp(-delta_t * n_params.inv_tau_m);
	return (I1 * helper->A0 + helper->A2 + aux);
}
extern double get_V_t_acc(struct neuron_helper_t *helper, double V0, double I0, double I1, double delta_t);

inline double get_V_t(struct neuron_helper_t *helper, double V0, double I0, double I1, double delta_t){
	/* Compute V(t0+delta_t) */
	double aux = (V0 - helper->A2 - I0 * helper->A0) * csexp(-delta_t * n_params.inv_tau_m);
	return (I1 * helper->A0 + helper->A2 + aux);
}

double getSelfSpikeTime(struct neuron_helper_t *n_helper){
	// Would hang
	if(n_helper->A2 < n_params.threshold) return -1;
	
	double delta_t;
	
	double Vt;
	double tmin = 0;
	double tmax = 50;
	
	// While Vt is below Vth
	while(1) {
		/* Now compute V(delta_t) */
		Vt = get_V_t_acc(n_helper, n_params.reset_potential, 0, 0, tmax);
		
		if (Vt >= n_params.threshold){
			printdbg("delta_t: %lf, Vt %lf > Vth %lf\n", tmax, Vt, n_params.threshold);
			break;
		}

		printdbg("delta_t: %lf, Vt %lf < Vth %lf\n", tmax, Vt, n_params.threshold);
		tmin = tmax;
		tmax *=2;
	} // We have a cap for the spike time
	
	// While V is farther than epsilon from Vth
	do {
		printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		delta_t = (tmax + tmin)/2;
		/* Now compute V(delta_t) */
		Vt = get_V_t_acc(n_helper, n_params.reset_potential, 0, 0, delta_t);
		
		if (Vt > n_params.threshold){
			printdbg("Vt %lf > Vth %lf\n", Vt, n_params.threshold);
			tmax = delta_t;
		} else {
			printdbg("Vt %lf < Vth %lf\n", Vt, n_params.threshold);
			tmin = delta_t;
		}
	} while(tmax > tmin + T_TOLERANCE);
	
	printdbg("SelfSpikeTime spikedelta: %lf\n", delta_t);
	return delta_t;
}

/* Compute the time after which the spike will take place given T0=0, V0 and I0,
 * on a neuron where Icond is < 0 and as such, no upper bound is available for the spike time. */
double findSpikeDeltaBinaryBlind(struct neuron_helper_t* helper, double V0, double I0){
	double Vt;
	double It;
	double tmin = 0;
	double tmax = helper->self_spike_time*2;

	// While Vt is below Vth
	printdbg("FSDBB: Neuron is self spiking\n");
	while(true) {
		/* Compute I(tmax) */
		It = get_I_f(tmax, I0);
		/* Now compute V(tmax) */
		Vt = get_V_t(helper, V0, I0, It, tmax);
		
		if (Vt >= n_params.threshold){
			printdbg("FSDBB: delta_t: %lf, Vt %lf > Vth %lf\n", tmax, Vt, n_params.threshold);
			break;
		}
		printdbg("FSDBB: delta_t: %lf, Vt %lf < Vth %lf\n", tmax, Vt, n_params.threshold);
		tmin = tmax;
		tmax *= 2;
	}
	
	// While V is farther than epsilon from Vth
	do {
		printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		double delta_t = (tmax + tmin)/2;
		/* Compute I(delta_t) */
		It = get_I_f(delta_t, I0);
		/* Now compute V(delta_t) */
		Vt = get_V_t(helper, V0, I0, It, delta_t);
		
		if (Vt > n_params.threshold){
			printdbg("Vt %lf > Vth %lf\n", Vt, n_params.threshold);
			tmax = delta_t;
		} else {
			printdbg("Vt %lf < Vth %lf\n", Vt, n_params.threshold);
			tmin = delta_t;
		}
		if(tmax <= tmin + T_TOLERANCE){
			printdbg("FSDBB spikedelta: %lf\n", delta_t);
			return delta_t;
		}
	} while(1);
}

void printNeuronState(neuron_state_t* state){
	printdbg("V: %lf, I: %lf, last_up: %lf, last_spike: %lf\n", state->membrane_potential, state->I, state->last_updated, state->last_fired);
}
