#include "1k_neurons_accuracy.h"
#include <stdio.h>
#include <errno.h>
//~ #define M_DEBUG                // Undefine to decrease verbosity

#ifdef M_DEBUG

#define printdbg(...) printf(__VA_ARGS__)
#else
#define printdbg(...)
#endif

#undef M_DEBUG

#define T_TOLERANCE 0.0001 // mS

// Parameters used in Shimoura18
const struct neuron_params_t n_params = {.threshold = -50.0,
    .reset_potential = -65.0,
    .inv_tau_m = 1 / 10.0,
    .inv_tau_syn = 1 / 0.5,
    .inv_C_m = 1 / 250.0,
    .refractory_period = 2.0};
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
void bring_to_present(neuron_state_t *state, simtime_t delta_update, simtime_t delta_spike);

/* Compute the next time at which the neuron will fire, if any */
simtime_t getNextFireTime(neuron_state_t *state, simtime_t delta_spike);

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

/* Get population index from neuron ID */
int n2pop(unsigned long int neuron_ID);

/* Compute the time it takes for a self-spiking neuron to spike with I=0 and V0=Vreset */
double getSelfSpikeTime(struct neuron_helper_t *params);

/* Initialize a LIF neuron */
neuron_state_t *InitLIFNeuron(unsigned long int me);

static struct neuron_helper_t n_helper_p[7];

/* NETWORK PARAMETERS v*/
// Population size per layer
//							IN		L1e		L1i		L2e
// L2i Out
unsigned int n_layer[7] = {100, 200, 200, 200, 200, 100};
// Tot neurons: 1000

//~ Probability connection table. src v\dest ->
double table[7][7] = {{0.0, 0.292, 0.192, 0.049, 0.237, 0.169, 0.115}, {0.0, 0.224, 0.293, 0.106, 0.254, 0.438, 0.099},
    {0.0, 0.135, 0.025, 0.409, 0.250, 0.309, 0.271}, {0.0, 0.165, 0.177, 0.122, 0.032, 0.491, 0.300},
    {0.0, 0.448, 0.319, 0.080, 0.207, 0.225, 0.201}, {0.0, 0.395, 0.123, 0.265, 0.215, 0.476, 0.174},
    {0.0, 0.223, 0.276, 0.358, 0.028, 0.065, 0.188}};
/* END NETWORK PARAMETERS ^ */

unsigned int protocol = 0;
unsigned int bg_type = 0;
unsigned int stim = 0;
unsigned int thal = 0;

FILE *topology_file = NULL;
// Make it stupid: just use some constants. FIXME: this needs to be parametrized
lp_id_t read_topology[1000][1000];
double init_potentials[1000];
bool topology_read_from_file = false;

void read_topology_from_file(FILE *file);

static const double exp_adj[] = {1.036199565629158, 1.027277981818330, 1.019341151369058, 1.012316468833528,
    1.006139353234098, 1.000751255315155, 0.996100136699148, 0.992140185571610, 0.988830196962552, 0.986132544047804,
    0.984013562521051, 0.982443170190283, 0.981394102987059, 0.980841516073624, 0.980763042138659, 0.981138430277712,
    0.981957925989939, 0.983178546312421, 0.984811614558841, 0.986834419454900, 0.989234236856811, 0.991999812972351,
    0.995121060808105, 0.998588377827231, 1.002392956865277, 1.006527340270889, 1.010984919213803, 1.015759089560808,
    1.020843804511367, 1.026234345941742, 1.031901110282296, 1.037916082604387};

enum { OPT_TOPOLOGY_FILE = 129 };

struct ap_option model_options[] = {{"topology_file", OPT_TOPOLOGY_FILE, "path/to/file",
					"The file holding the topology, in format:\n"
					"Neuron, [Targets]\n"
					"1, 10, 14, 11, 50, 20\n"
					"2, 54, 214, 11\n"
					"...\n"},
    {0}};

void model_parse(int key, const char *arg)
{
	switch(key) {
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
					topology_read_from_file = true;
					fclose(topology_file);
				} else {
					if(errno == 2) {
						fprintf(stderr,
						    "\e[0;31m[ERROR]\e[0m Could not open file %s while loading topology, file does not exist.\n",
						    arg);
					} else {
						fprintf(stderr,
						    "\e[0;31m[ERROR]\e[0m Could not open file %s while loading topology.\n",
						    arg);
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

// FIXME OK FROM HERE ONWARDS!!

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
            } else if (g_w_in == 0 && aux < 0) {
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
void *NeuronInit(unsigned long int me)
{
	printdbg("[N%lu] INITIALIZING\n", me);
	return InitLIFNeuron(me);
}

neuron_state_t *InitLIFNeuron(unsigned long int me)
{
	neuron_state_t *state = malloc(sizeof(neuron_state_t));

	// Check which population I'm in, and check whether my helper state has been initialized already
	int pop = n2pop(me);

	struct neuron_helper_t *n_helper = &n_helper_p[pop];

	if(!n_helper->self_spike_time) {
		// Init physical params
		n_helper->Iext = pop ? 0.0 : g_I_ext; // Neurons in pop 0 (IDs 0-99) take Iext, others don't
		n_helper->A0 = n_params.inv_C_m / (n_params.inv_tau_m - n_params.inv_tau_syn);
		n_helper->A2 = n_params.reset_potential + ((n_helper->Iext * n_params.inv_C_m) / n_params.inv_tau_m);

		n_helper->Icond =
		    ((n_params.threshold - n_params.reset_potential) * (n_params.inv_tau_m)) / n_params.inv_C_m -
		    (n_helper->Iext);
		n_helper->self_spike_time = getSelfSpikeTime(n_helper);

		printf("Neuron population %d Icond: %lf\n", pop, n_helper->Icond);
		printf("Neuron population %d self spike time: %lf\n", pop, n_helper->self_spike_time);
	}
	state->helper = n_helper;
	state->last_updated = 0;
	state->I = 0.0;
	state->last_fired = -(n_params.refractory_period + 1);

	//~ state->membrane_potential = -58.0 + 10 * Normal();
	state->membrane_potential = init_potentials[me];
	if(state->membrane_potential > n_params.threshold) {
		state->membrane_potential = n_params.threshold + 0.01;
	}
	return state;
}

/* neuron receives spike */
void NeuronHandleSpike(unsigned long int me, simtime_t now, double value, void *neuron_state)
{
	neuron_state_t *state = (neuron_state_t *)neuron_state;

	simtime_t fireTime;

	simtime_t delta_spike = now - state->last_fired;
	// This brings the neuron to present. Computes and updates I and V.
	bring_to_present(state, now - state->last_updated, delta_spike);
	// Apply the spike
	state->I += value;
	state->last_updated = now;

	/* Given current state, when will the neuron spike next, if at all? */
	fireTime = getNextFireTime(state, delta_spike);

	if(fireTime > 0.0) { // Spike might happen in future
		MaybeSpikeAndWake(me, fireTime);
	}
}

/*  Neuron has been woken by the MaybeSpikeAndWake */
void NeuronWake(unsigned long int me, simtime_t now, void *n_state)
{
	neuron_state_t *state = n_state;
	simtime_t t_since_last_eval = now - state->last_updated;
	state->I = get_I_f(t_since_last_eval, state->I);
	state->membrane_potential = n_params.reset_potential;
	state->last_fired = now;
	state->last_updated = now;

	/* Given current state, when will the neuron spike next, if at all? */
	double fireTime = getNextFireTime(state, 0);

	if(fireTime > 0.0) { // Spike might happen in future
		MaybeSpikeAndWake(me, fireTime);
	}
}

/* synapse handles the spike by updating its state and *returning* the spike intensity */
double SynapseHandleSpike(simtime_t now, unsigned long int src_neuron, unsigned long int dest_neuron, synapse_t *state)
{
	return state->weight;
}

/* is the neuron done? */
bool NeuronCanEnd(unsigned long int me, neuron_state_t *state)
{
	return false;
}

void ProbeRead(simtime_t now, unsigned long int monitored_neuron, const neuron_state_t *neuron_state)
{
	return;
}

void GatherStatistics(simtime_t now, unsigned long int neuron_id, const neuron_state_t *state)
{
	return;
}

void ReadTopologyInit();

/* Init the topology */
void SNNInitTopology(unsigned long int neuron_count)
{
	ReadTopologyInit();
}

// FIXME ugly magic numbers!
void ReadTopologyInit()
{
	synapse_t *synapse;
	double delay;
	for(lp_id_t src_neuron = 0; src_neuron < 1000; src_neuron++) {
		for(int cur = 0; read_topology[src_neuron][cur] != -1; cur++) {
			lp_id_t dst_neuron = read_topology[src_neuron][cur];

			int pop = n2pop(src_neuron);

			if((pop % 2) || !pop) { // Excitatory connections
				delay = g_d_ex;
				if(delay < 0.1) {
					delay = 0.1;
				}
				synapse = NewSynapse(src_neuron, dst_neuron, sizeof(synapse_t), true, delay);
				if(synapse == NULL) {
					continue;
				}
				synapse->weight = g_w_ex;
				if(synapse->weight < 0.0) {
					synapse->weight = 0.0;
				}

			} else { // Inhibitory connections
				delay = g_d_in;
				if(delay < 0.1) {
					delay = 0.1;
				}
				synapse = NewSynapse(src_neuron, dst_neuron, sizeof(synapse_t), true, delay);
				if(synapse == NULL) {
					continue;
				}
                synapse->weight = g_w_in;
				if(synapse->weight > 0.0) {
					synapse->weight = 0.0;
				}
			}
		}
	}
}

/* Get population index from neuron ID */
int n2pop(unsigned long int neuron_ID)
{
	for(int i = 0; i < 7; i++) {
		if(neuron_ID < n_layer[i]) {
			return i;
		}
		neuron_ID -= n_layer[i];
	}
	return -1;
}

/* Called at the start of the spike handling. Updates the state of the neuron from the last update to now.
 * WARNING: does side effect on the state but does not update state->last_updated!
 * We are certain that the neuron did not spike from the last time it was updated to now.
 * So we don't need to do calculations to account for that */
void bring_to_present(neuron_state_t *state, simtime_t delta_update, simtime_t delta_spike)
{
	/*
	 * delta_update = now - state->last_updated
	 * delta_spike = now - state->last_spiked
	 * */

	double If = get_I_f(delta_update, state->I);
	double Ient = state->I;
	state->I = If;

	if(delta_spike < n_params.refractory_period)
		return;

	// Amount of time of refractory period that we still have to take into account
	double remaining = delta_update + n_params.refractory_period - delta_spike;

	if(remaining > 0.0) { // computes current at end of refractory period
		Ient = get_I_f(remaining, state->I);
		delta_update -= remaining;
	}

	// Update membrane potential from end of refractory to now
	state->membrane_potential = get_V_t(state->helper, state->membrane_potential, Ient, If, delta_update);
}

// delta_spike is how long has passed since the last time the neuron has spiked
simtime_t getNextFireTime(neuron_state_t *state, simtime_t delta_spike)
{
	struct neuron_helper_t *n_helper = state->helper;

	/* Necessary condition for the membrane to reach threshold potential:
	 * dV/dt>0 with V->Vth
	 * As potential approaches the threshold value, the leakage will be stronger.
	 * We need an I that can overcome it */
	double Icond = state->helper->Icond;
	double I0 = state->I;
	simtime_t t0 = state->last_updated;
	double V0 = state->membrane_potential;

	/* Check that the condition for possibly spiking is met.
	 * Either neuron is self-spiking, or I0 >= Icond */
	if(n_helper->A2 < n_params.threshold && I0 < Icond) { // Neuron is not self-spiking and I<Icond
		return -1;                                    // Can never spike this way. Return
	}

	// Time until end of refractory period
	simtime_t dt_ref_end;
	/* If we are still in refractory period, get I at the end of it */
	if(n_params.refractory_period > delta_spike) {
		// Current at end of refractory period. V is clamped at reset.
		dt_ref_end = n_params.refractory_period - delta_spike;
		I0 = get_I_f(dt_ref_end, state->I);
	} else {
		dt_ref_end = 0;
	}
	/* Refractory done */

	/* A2>Vth => neuron is self spiking */
	if(n_helper->A2 > n_params.threshold) { // Neuron is self spiking
		return t0 + dt_ref_end + findSpikeDeltaBinaryBlind(n_helper, V0, I0);
	}

	// Neuron is not self spiking
	if(I0 < Icond) {
		return -1; // Can never spike this way. Return
	}

	/* Now let's find the t_upper(/t_to_icond), the time for(/after) which I = Icond.
	 * t_upper is an upper bound to the spike time. We have two cases:
	 * a) V(t_upper) >= Vth. The neuron spikes in a time between t0 and t_upper.
	 * b) V(t_upper) < Vth. The neuron did not spike and will not in the future (since from t_upper on, I<Icond)
	 * These are the only two cases: since Icond is needed to surpass Vth,
	 * any I > Icond is sufficient to keep V > Vth.
	 * As such, if there is a time t' in [t0, t_upper] s.t. V(t') > Vth,
	 * this will stay true for at least as long as I>Icond. i.e. for all t in [t', t_upper].
	 * Thus, if V(t_upper) >= Vth, we crossed Vth once in [t0,t_upper]. Otherwise, we never did. QED
	 * */
	simtime_t t_to_icond = -(log(Icond / I0) / n_params.inv_tau_syn); // Divide by the inverse

	/* V(t_upper) = V(t0 + dt_ref_end + t_to_icond) */
	double V = get_V_t(n_helper, V0, I0, Icond, t_to_icond);

	if(V < n_params.threshold) { // Case b: neuron does not spike
		return -1;
	}
	// Case a: V(t_upper)>=Vth. Neuron spiked at ts in [t0, t_upper]
	return t0 + dt_ref_end + findSpikeDeltaBinary(n_helper, V0, I0, t_to_icond);
}

/* Compute the time after which the spike will take place given T0=0, V0 and I0, on a NON-self-spiking neuron. */
double findSpikeDeltaBinary(struct neuron_helper_t *n_helper, double V0, double I0, double time_interval)
{
	double It;
	double Vt;
	double tmin = 0;
	double tmax = time_interval;

	// While V is farther than epsilon from Vth
	do {
		double delta_t = (tmax + tmin) / 2;
		/* Compute I(t0 + delta_t) & V(t0 + delta_t) */
		It = get_I_f(delta_t, I0);
		Vt = get_V_t(n_helper, V0, I0, It, delta_t);

		if(Vt > n_params.threshold) {
			tmax = delta_t;
		} else {
			tmin = delta_t;
		}
		if(tmax <= tmin + T_TOLERANCE) {
			return delta_t;
		}
	} while(1);
}

inline double get_I_f(double delta_t, double I_i)
{
	return exp(-delta_t * n_params.inv_tau_syn) * I_i;
}

inline double get_V_t_acc(struct neuron_helper_t *helper, double V0, double I0, double I1, double delta_t)
{
	/* Now compute V(delta_t) and return it */
	double aux = (V0 - helper->A2 - I0 * helper->A0) * exp(-delta_t * n_params.inv_tau_m);
	return (I1 * helper->A0 + helper->A2 + aux);
}

extern double get_V_t_acc(struct neuron_helper_t *helper, double V0, double I0, double I1, double delta_t);

inline double get_V_t(struct neuron_helper_t *helper, double V0, double I0, double I1, double delta_t)
{
	/* Compute V(t0+delta_t) */
	double aux = (V0 - helper->A2 - I0 * helper->A0) * csexp(-delta_t * n_params.inv_tau_m);
	return (I1 * helper->A0 + helper->A2 + aux);
}

double getSelfSpikeTime(struct neuron_helper_t *n_helper)
{
	// Would hang
	if(n_helper->A2 < n_params.threshold)
		return -1;

	double delta_t;

	double Vt;
	double tmin = 0;
	double tmax = 50;

	// While Vt is below Vth
	while(1) {
		/* Now compute V(delta_t) */
		Vt = get_V_t_acc(n_helper, n_params.reset_potential, 0, 0, tmax);

		if(Vt >= n_params.threshold) {
			printdbg("delta_t: %lf, Vt %lf > Vth %lf\n", tmax, Vt, n_params.threshold);
			break;
		}

		printdbg("delta_t: %lf, Vt %lf < Vth %lf\n", tmax, Vt, n_params.threshold);
		tmin = tmax;
		tmax *= 2;
	} // We have a cap for the spike time

	// While V is farther than epsilon from Vth
	do {
		printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		delta_t = (tmax + tmin) / 2;
		/* Now compute V(delta_t) */
		Vt = get_V_t_acc(n_helper, n_params.reset_potential, 0, 0, delta_t);

		if(Vt > n_params.threshold) {
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
double findSpikeDeltaBinaryBlind(struct neuron_helper_t *helper, double V0, double I0)
{
	double Vt;
	double It;
	double tmin = 0;
	double tmax = helper->self_spike_time * 2;

	// While Vt is below Vth
	while(true) {
		/* Compute I(tmax) & V(tmax) */
		It = get_I_f(tmax, I0);
		Vt = get_V_t(helper, V0, I0, It, tmax);

		if(Vt >= n_params.threshold) {
			break;
		}
		tmin = tmax;
		tmax *= 2;
	}

	// While time window is bigger than tolerance, bisection
	do {
		double delta_t = (tmax + tmin) / 2;
		It = get_I_f(delta_t, I0);
		Vt = get_V_t(helper, V0, I0, It, delta_t);

		if(Vt > n_params.threshold) {
			tmax = delta_t;
		} else {
			tmin = delta_t;
		}
		if(tmax <= tmin + T_TOLERANCE) {
			return delta_t;
		}
	} while(1);
}

inline bool doubleEquals(double a, double b, double epsilon)
{
	return fabs(a - b) < epsilon;
}

extern bool doubleEquals(double a, double b, double epsilon);
