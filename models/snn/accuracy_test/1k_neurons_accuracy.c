#include "1k_neurons_accuracy.h"
#include <stdio.h>
#include <errno.h>
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
#define I_TOLERANCE 0.01 //pA
#define T_TOLERANCE 0.01 //mS

#define CONNECTIONS 1

const struct neuron_params_t n_params = {
	.threshold = -50.0,
	.reset_potential = -65.0,
	.inv_tau_m = 1/10.0,
	.inv_tau_syn = 1/0.5,
	.inv_C_m = 1/250.0,
	.refractory_period = 2.0
};
// Synapses parameters
double g_d_ex=1.5; //# [ms] Excitatory delay
double g_d_in=0.8; //# [ms] Inhibitory delay
#define g_tau_syn		0.5 	//# [ms] Post-synaptic current time constant

double g_w_ex=0; // Excitatory weight
double g_w_in=0; // Inhibitory weight
#define in_g			4.0						//# Multiplicative factor for the inhibitory weights

// External current
double g_I_ext=0; // External current for first population

/* Called at the start of the spike handling. Updates the state of the neuron from the last update to now */
void bring_to_present(neuron_state_t* state, simtime_t now);
/* Compute the next time at which the neuron will fire, if any */
simtime_t getNextFireTime(neuron_state_t* state);
/* Are a and b within epsilon from each other? */
inline bool doubleEquals(double a, double b, double epsilon);
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

/* Initializes the topology of Shimoura18 */
void ShimouraTopology();

/* Generates an array a s.t. a[i] contains the index of the first neuron of population i */
void gen_indexes(unsigned int * popsizes, unsigned int *out, int size);
/* Get population index from neuron ID */
int n2pop(unsigned long int neuron_ID);
/* Compute the time it takes for a self-spiking neuron to spike with I=0 and V0=Vreset */
double getSelfSpikeTime(struct neuron_helper_t *params);

/* Initialize a LIF neuron */
neuron_state_t* InitLIFNeuron(unsigned long int me);

static struct neuron_helper_t n_helper_p[7];

/* NETWORK PARAMETERS v*/
// Population size per layer
//							IN		L1e		L1i		L2e		L2i		Out
unsigned int n_layer[7] = {	100,	200,	200,	200,	200,	100};
//Tot neurons: 1000

//~ Probability connection table. src v\dest ->
double table[7][7] = {\
		{0.0, 0.292, 0.192, 0.049, 0.237, 0.169, 0.115},\
		{0.0, 0.224, 0.293, 0.106, 0.254, 0.438, 0.099},\
		{0.0, 0.135, 0.025, 0.409, 0.250, 0.309, 0.271},\
		{0.0, 0.165, 0.177, 0.122, 0.032, 0.491, 0.300},\
		{0.0, 0.448, 0.319, 0.080, 0.207, 0.225, 0.201},\
		{0.0, 0.395, 0.123, 0.265, 0.215, 0.476, 0.174},\
		{0.0, 0.223, 0.276, 0.358, 0.028, 0.065, 0.188}};
/* END NETWORK PARAMETERS ^ */

unsigned int protocol = 0;
unsigned int bg_type = 0;
unsigned int stim = 0;
unsigned int thal = 0;

FILE* topology_file = NULL;
// Make it stupid: just use some constants. FIXME: this needs to be parametrized
lp_id_t read_topology[1000][1000];
double init_potentials[1000];
bool topology_read_from_file = false;
void read_topology_from_file(FILE * file);

enum{
    OPT_TOPOLOGY_FILE = 129
};

struct ap_option model_options[] = {
	{"topology_file", OPT_TOPOLOGY_FILE, "path/to/file",
		"The file holding the topology, in format:\n"
		"Neuron, [Targets]\n"
        "1, 10, 14, 11, 50, 20\n"
        "2, 54, 214, 11\n"
        "...\n"
	},
    {0}
};

void model_parse (int key, const char *arg){
	switch (key) {
		case AP_KEY_FINI:
			break;
		case AP_KEY_INIT:
			break;
		case OPT_TOPOLOGY_FILE:
		{
			errno = 0;
			topology_file = fopen(arg, "r");
				if(topology_file != NULL) {
					printf("Reading the topology from file:\n%s\n", arg);
			read_topology_from_file(topology_file);
			topology_read_from_file=true;
			fclose(topology_file);
				} else {
			if(errno==2){
				fprintf(stderr, "\e[0;31m[ERROR]\e[0m Could not open file %s while loading topology, file does not exist.\n", arg);
			} else {
				fprintf(stderr, "\e[0;31m[ERROR]\e[0m Could not open file %s while loading topology.\n", arg);
				fprintf(stderr, "Error code: %d\n", errno);
			}
			abort();
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

void read_topology_from_file(FILE * file){
    // Read the file line by line and load the topology into the global array read_topology
    // read_topology[i] contains the list of postsynaptic neurons

    // Make it stupid: for now we suppose no line is bigger than 1MB
    // FIXME: this will need to be changed.
    const int BUF_SZ = 1<<20;
    char read_buf[BUF_SZ];
    read_buf[BUF_SZ - 1] = 'a';
    // First line is header
    fgets(read_buf, BUF_SZ, file);
    while(fgets(read_buf, BUF_SZ, file) != NULL) {
        // Overwrite last \n
        int lastchar = strlen(read_buf)-1;
        if(read_buf[lastchar]=='\n'){
            read_buf[lastchar] = '\0';
        }
        // For each token
        // Read src neuronID
        char* tk = strtok(read_buf, " ,");
        long src = strtol(tk, NULL, 10);
        // Starting potential
        tk = strtok(NULL, ", ");
        init_potentials[src] = strtod(tk, NULL);
        // Input current
        tk = strtok(NULL, ", ");
        if (g_I_ext == 0) {
            double aux = strtod(tk, NULL);
            if (aux > 0) {
                g_I_ext = aux;
                printf("Input current Iext: %lf\n", g_I_ext);
            }
        }
        // Synapse weight
        tk = strtok(NULL, ", ");
        bool skip_delay = false;
        if (g_w_ex == 0 || g_w_in == 0) {
            double aux = strtod(tk, NULL);
            tk = strtok(NULL, ", ");
            skip_delay = true;
            if (g_w_ex == 0 && aux > 0) {
                g_w_ex = aux;
                printf("Excitatory weight: %lf\n", aux);
                g_d_ex = strtod(tk, NULL);
                printf("Excitatory delay: %lf\n", g_d_ex);
            } else if (g_w_in == 0) {
                g_w_in = aux;
                printf("Inhibitory weight: %lf\n", aux);
                g_d_in = strtod(tk, NULL);
                printf("Inhibitory delay: %lf\n", g_d_in);
            }
        }
        // Synapse delay
        if (!skip_delay) strtok(NULL, ", ");
        // Postsynaptic neurons
        tk = strtok(NULL, ", ");
        int cur = 0;
        for(; tk!=NULL; cur++) {
            //~ int strl = strlen(tk);
            //~ printf("%s, len=%d\n", tk, strl);
            long dst = strtol(tk, NULL, 10);
            read_topology[src][cur] = dst;

            tk = strtok(NULL, ", ");
        }
        // End the table entries with a recognizable value
        read_topology[src][cur] = -1;
    }
}




/* Initialize the neuron */
void* NeuronInit(unsigned long int me){

	printdbg("[N%lu] INITIALIZING\n", me);
	return InitLIFNeuron(me);
}

neuron_state_t* InitLIFNeuron(unsigned long int me){
	neuron_state_t* state = malloc(sizeof(neuron_state_t));

	//Check which population I'm in, and check whether my helper state has been initialized already
	int pop = n2pop(me);

	struct neuron_helper_t *n_helper = &n_helper_p[pop];

	if (!n_helper->self_spike_time){
		//struct neuron_helper_t *n_helper = st_malloc(sizeof(struct neuron_helper_t));

		n_helper->Iext = pop ? 0.0 : g_I_ext; // Neurons in pop 0 (IDs 0-99) take Iext, others don't
		printdbg("[Pop %d] Iext: %lf\n", pop, n_helper->Iext);
		n_helper->A0 = n_params.inv_C_m/(n_params.inv_tau_m - n_params.inv_tau_syn);
		n_helper->A2 = n_params.reset_potential + ((n_helper->Iext * n_params.inv_C_m)/n_params.inv_tau_m);
		printdbg("A0: %lf, A2: %lf\n", n_helper_p->A0, n_helper_p->A2);

		n_helper->Icond = ((n_params.threshold - n_params.reset_potential) * (n_params.inv_tau_m)) / n_params.inv_C_m - (n_helper->Iext);
		printf("Neuron population %d Icond: %lf\n", pop, n_helper->Icond);
		n_helper->self_spike_time = getSelfSpikeTime(n_helper);
		printf("Neuron population %d self spike time: %lf\n", pop, n_helper->self_spike_time);

		n_helper_p[pop] = *n_helper;
	}
	state->helper = &(n_helper_p[pop]);
	state->last_updated = 0;
	state->I = 0.0;
	state->last_fired = -(n_params.refractory_period + 1);

	//~ state->membrane_potential = -58.0 + 10 * Normal();
	state->membrane_potential = -58.0;
	if(state->membrane_potential > n_params.threshold){
		state->membrane_potential = n_params.threshold + 0.01;
	}
	return state;
}

/* neuron receives spike */
void NeuronHandleSpike(unsigned long int me, simtime_t now, double value, void* neuron_state){
	printdbg("[N%lu] Spike received at %lf\n", me, now);

	neuron_state_t* state = (neuron_state_t*) neuron_state;

	simtime_t fireTime = -1.0;

	printdbg("[N%lu] Spike received with value: %lf at %lf\n", me, value, now);

	// This brings the neuron to present. Computes and updates I and V.
	bring_to_present(state, now);

	// Apply the spike
	state->I += value;
	state->last_updated = now; // Unnecessary, bring_to_present already updates it

	printdbg("[N%lu] updated values: ", me);

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
	state->last_fired = now;
	state->last_updated = now;

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
}

/* is the neuron done? */
bool NeuronCanEnd(unsigned long int me, neuron_state_t* state){
	return false;
}

void ProbeRead(simtime_t now, unsigned long int monitored_neuron, const neuron_state_t* neuron_state){
	return;
}

void GatherStatistics(simtime_t now, unsigned long int neuron_id, const neuron_state_t* state){
	return;
}

void ThousandNeuronsTopology();
void ReadTopologyInit();

/* Init the topology */
void SNNInitTopology(unsigned long int neuron_count){
    ReadTopologyInit();
}

// FIXME ugly magic numbers!
void ReadTopologyInit(){
    synapse_t* synapse;
	double delay;
    for(lp_id_t src_neuron=0; src_neuron<1000; src_neuron++){
        for(int cur = 0; read_topology[src_neuron][cur]!=-1; cur++){
            lp_id_t dst_neuron = read_topology[src_neuron][cur];

            int pop = n2pop(src_neuron);

            if ((pop % 2) || !pop){ // Excitatory connections
                delay = g_d_ex;
                if(delay < 0.1){delay=0.1;}
                synapse = NewSynapse(src_neuron, dst_neuron, sizeof(synapse_t), true, delay);
                if(synapse==NULL){continue;}
                synapse->weight = g_w_ex;
                if(synapse->weight < 0.0){synapse->weight = 0.0;}

            } else { // Inhibitory connections
                delay = g_d_in;
                if(delay < 0.1){delay=0.1;}
                synapse = NewSynapse(src_neuron, dst_neuron, sizeof(synapse_t), true, delay);
                if(synapse==NULL){continue;}
                synapse->weight = -(g_w_ex);
                if(synapse->weight > 0.0){synapse->weight = 0.0;}
                synapse->weight *= in_g;
            }
        }
    }
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
	for(int i=0; i<7; i++){
		if(neuron_ID < n_layer[i]){
			return i;
		}
		neuron_ID -= n_layer[i];
	}
	return -1;
}

//OK
/* Called at the start of the spike handling. Updates the state of the neuron from the last update to now.
 * We are certain that the neuron did not spike from the last time it was updated to now.
 * So we don't need to do calculations to account for that */
void bring_to_present(neuron_state_t* state, simtime_t now){
	struct neuron_helper_t *n_helper = state->helper;

	/*
	 * t0 = state->last_updated
	 * tf = now
	 * I0 = state->I
	 * V0 = state->membrane_potential
	 * */

	double delta_t;
	double If;

	simtime_t ref_end = state->last_fired + n_params.refractory_period; // When will refractory period end

	if(ref_end > state->last_updated){ // The neuron is in refractory period at t0
		printdbg("BTP: neuron is still in refractory period: last_updated: %lf, ref_end: %lf\n", state->last_updated, ref_end);

		state->membrane_potential = n_params.reset_potential; // Potential is clamped to reset potential

		if (ref_end >= now){ // Refractory lasts at least till now time. Only update the current and don't touch membrane potential
			printdbg("BTP: refractory period lasts further than now: %lf\n", now);
			delta_t = now - state->last_updated;
			If = get_I_f(delta_t, state->I);

			state->I = If;
			state->last_updated = now;
			printdbg("BTP: computed new I\n");
			return;

		} else { // Skip the remaining refractory period by updating I to the value it will have at the end of it
			printdbg("BTP: handling refractory period\n");
			delta_t = ref_end - state->last_updated;
			If = get_I_f(delta_t, state->I);

			//Update I
			state->I=If;
			state->last_updated = ref_end;
			printdbg("BTP: computed new I\n");
		}
	}
	printdbg("BTP: bringing neuron values to present\n");
	// Now refractory is out of the way. Compute the new values for I and Vm
	delta_t = now - state->last_updated;

	printdbg("BTP: delta_t: %lf\n", delta_t);

	If = get_I_f(delta_t, state->I);

	printdbg("BTP: If: %lf\n", If);

	double V = get_V_t(n_helper, state->membrane_potential, state->I, If, delta_t);

	printdbg("BTP: Vf: %lf\n", V);

	state->membrane_potential = V;
	state->I = If;
	state->last_updated = now;
}

simtime_t getNextFireTime(neuron_state_t* state){
	struct neuron_helper_t* n_helper = state->helper;

	/* Necessary condition for the membrane to reach threshold potential:
	 * dV/dt>0 with V->Vth
	 * Indeed, as potential approaches the threshold value, the leakage will be stronger.
	 * We need an I that can overcome it */
	double Icond = state->helper->Icond;

	printdbg("getFT: Icond: %lf\n", Icond);

	// Down here only non-self-spiking neurons are considered.
	double I0 = state->I;
	simtime_t t0 = state->last_updated;
	double V0 = state->membrane_potential;

	/* Check that the condition for spiking is met.
	 * Since I decays towards 0 with the passing of time (regardless of whether I is positive or negative),
	 * if the spiking condition is not met at the start, it will never be met.
	 * Unless a spike is received. But that is a future event.*/
	if (n_helper->A2 < n_params.threshold && I0 < Icond){// Neuron is not self-spiking and I<Icond
		printdbg("getFT: I %lf < Icond %lf. Will never spike\n", I0, Icond);
		return -1; // Can never spike this way. Return
	}

	double delta_t;
	double ref_end = state->last_fired + n_params.refractory_period;

	printdbg("getFT: Ref_end: %lf\n", ref_end);

	/* Now check if we are still in refractory period */
	if(ref_end > state->last_updated){

		t0 = ref_end; // When will refractory period end, new t0

		delta_t = t0 - state->last_updated;
		I0 = get_I_f(delta_t, I0);
		V0 = n_params.reset_potential; // Potential is clamped to reset potential
		printdbg("getFT: Ref period out of the way\n");
		printdbg("New V0: %lf, I0: %lf, t0: %lf\n", V0, I0, t0);
	}
	/* Refractory done */

	/* If A2>Vth the neuron is self spiking and we need to
	 * check when it will spike regardless of the value of I*/
	if(n_helper->A2 > n_params.threshold){ // Neuron is self spiking
		delta_t = findSpikeDeltaBinaryBlind(n_helper, V0, I0);
		if (delta_t<=0){
			printf("Something went wrong\n");
			abort();
		}
		return t0 + delta_t;
	}

	// Neuron is not self spiking and we proceed with the correct approach
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
	double delta_t = 0;

	double It;
	double Vt;
	double tmin = 0;
	double tmax = time_interval;
	unsigned int its = 0;

	// While V is farther than epsilon from Vth
	do {
		printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		delta_t = (tmax + tmin)/2;
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
		its++;
	} while(!doubleEquals(Vt, n_params.threshold, V_TOLERANCE) || !doubleEquals(tmin, tmax, T_TOLERANCE));


	printdbg("BinarySearch spikedelta: %lf\n", delta_t);
	return delta_t;
}

inline double get_I_f(double delta_t, double I_i){
	return exp(-delta_t * n_params.inv_tau_syn) * I_i;
}

inline double get_V_t(struct neuron_helper_t *helper, double V0, double I0, double I1, double delta_t){
	/* Now compute V(delta_t) and return it */
	double aux = (V0 - helper->A2 - I0 * helper->A0) * exp(-delta_t * n_params.inv_tau_m);
	return (I1 * helper->A0 + helper->A2 + aux);
}

double getSelfSpikeTime(struct neuron_helper_t *n_helper){
	// Would hang
	if(n_helper->A2 < n_params.threshold) return -1;

	double delta_t = 0;

	double Vt;
	double tmin = 0;
	double tmax = 50;
	double step = 50;

	// While Vt is below Vth
	while(tmin < MODEL_MAX_SIMTIME){// !isinf(delta_t)) {
		/* Now compute V(delta_t) */
		Vt = get_V_t(n_helper, n_params.reset_potential, 0, 0, tmax);

		if (Vt >= n_params.threshold){
			printdbg("delta_t: %lf, Vt %lf > Vth %lf\n", tmax, Vt, n_params.threshold);
			break;
		} else {
			printdbg("delta_t: %lf, Vt %lf < Vth %lf\n", tmax, Vt, n_params.threshold);
			tmin = tmax;
			//~ tmax = tmax + step;
			tmax *=2;
		}
	}

	if(tmin > MODEL_MAX_SIMTIME){//isinf(delta_t)){
		return -1;
	}

	// While V is farther than epsilon from Vth
	do {
		printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		delta_t = (tmax + tmin)/2;
		/* Now compute V(delta_t) */
		Vt = get_V_t(n_helper, n_params.reset_potential, 0, 0, delta_t);

		if (Vt > n_params.threshold){
			printdbg("Vt %lf > Vth %lf\n", Vt, n_params.threshold);
			tmax = delta_t;
		} else {
			printdbg("Vt %lf < Vth %lf\n", Vt, n_params.threshold);
			tmin = delta_t;
		}
	} while(!doubleEquals(Vt, n_params.threshold, V_TOLERANCE) || !doubleEquals(tmin, tmax, T_TOLERANCE));

	printdbg("SelfSpikeTime spikedelta: %lf\n", delta_t);
	return delta_t;
}

/* Compute the time after which the spike will take place given T0=0, V0 and I0,
 * on a neuron where Icond is < 0 and as such, no upper bound is available for the spike time. */
double findSpikeDeltaBinaryBlind(struct neuron_helper_t* helper, double V0, double I0){
	double delta_t = 0;

	double Vt;
	double It;
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
		/* Compute I(delta_t) */
		It = get_I_f(delta_t, I0);
		/* Now compute V(delta_t) */
		Vt = get_V_t(helper, V0, I0, It, delta_t);

		if (Vt >= n_params.threshold){
			printdbg("FSDBB: delta_t: %lf, Vt %lf > Vth %lf\n", delta_t, Vt, n_params.threshold);
			break;
		} else {
			printdbg("FSDBB: delta_t: %lf, Vt %lf < Vth %lf\n", delta_t, Vt, n_params.threshold);
			tmin = tmax;
			tmax = tmax + step;
			//~ tmax *=2;
		}
	}

	// While V is farther than epsilon from Vth
	do {
		printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		delta_t = (tmax + tmin)/2;
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
	} while(!doubleEquals(Vt, n_params.threshold, V_TOLERANCE) && !doubleEquals(tmin, tmax, T_TOLERANCE));

	printdbg("FSDBB spikedelta: %lf\n", delta_t);
	return delta_t;
}

inline bool doubleEquals(double a, double b, double epsilon){
	return fabs(a-b) < epsilon;
}
extern bool doubleEquals(double a, double b, double epsilon);
