#include "cuba.h"
#include <stdio.h>
#include <assert.h>

//~ #define M_DEBUG                // Undefine to decrease verbosity

#ifdef M_DEBUG

#define printdbg(...) printf(__VA_ARGS__)
#else
#define printdbg(...)
#endif

#undef M_DEBUG

#define MODEL_MAX_SIMTIME global_config.termination_time

#define T_TOLERANCE 0.01 // mS

#define e_portion 0.8f

#define g_El (-49.0f)          // mV
#define g_we (60 * 0.27f / 10) // pA
#define g_wi (-20 * 4.5f / 10) // pA
// Possible syn delay #1
// #define g_syn_delay 1.0f	//mS
/* Syn delay used in many other models. Standard delay used in
 * Brette et al. 2007
 * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2638500/
 * (see Section 4.3.3 - Performance)
 * */
#define g_syn_delay 0.1f // mS

unsigned long long int syn_ct = 0;

// Parameters used in CUBA benchmark
static struct neuron_params_t n_params = {
    .inv_tau_m = 1 / 20.0,             // [1/ms]
    .inv_tau_e = 1 / 5.0,              // [1/ms]
    .inv_tau_i = 1 / 10.0,             // [1/ms]
    .De = 1.0 / (1.0 - (20.0 / 5.0)),  // Inverse of (1-tau_m/tau_e)
    .Di = 1.0 / (1.0 - (20.0 / 10.0)), // Inverse of (1-tau_m/tau_i)
    .reset_potential = -60.0,          // [mV]
    .threshold = -50.0,                // [mV]
    .refractory_period = 5.0           // [ms]
};

/* Called at the start of the spike handling. Updates the state of the neuron
 * from the last update to now */
void bring_to_present(neuron_state_t *state, simtime_t delta_upd,
    simtime_t delta_spike);

/* Compute the next time at which the neuron will fire, if any */
simtime_t getNextFireTime(neuron_state_t *state, simtime_t delta_spike);

/* Compute the time after which the spike will take place given T0=0, V0 and I0,
 * on a self-spiking neuron. */
double findSpikeDeltaBinaryBlind(struct neuron_helper_t *helper, double V0,
    double Ge0, double Gi0);

/* Compute Ge(delta_t) given Ge0, delta_t, and Tau_e */
inline double get_Ge_f(double delta_t, double Ge0);

extern double get_Ge_f(double delta_t, double Ge0);

/* Compute Gi(delta_t) given Gi0, delta_t, and Tau_i */
inline double get_Gi_f(double delta_t, double Gi0);

extern double get_Gi_f(double delta_t, double Gi0);

/* Compute V(delta_t) given V0, I0, delta_t, and the parameters*/
inline double get_V_t(struct neuron_helper_t *helper, double V0, double Ge0,
    double Ge1, double Gi0, double Gi1, double delta_t);

extern double get_V_t(struct neuron_helper_t *helper, double V0, double Ge0,
    double Ge1, double Gi0, double Gi1, double delta_t);

/* Initializes the topology of the benchmark */
void CUBATopology(unsigned long int neuron_count);

void printNeuronState(neuron_state_t *state);

/* Generates an array a s.t. a[i] contains the index of the first neuron of
 * population i */
void gen_indexes(unsigned int *popsizes, unsigned int *out, int size);

/* Get population index from neuron ID */
int n2pop(unsigned long int neuron_ID);

/* Extract a value from a binomial disribution */
unsigned random_binomial(unsigned trials, double p);
/* Is a population excitatory or inhibitory? */
bool is_excitatory(unsigned int population_id);

/* Compute the time it takes for a self-spiking neuron to spike with I=0 and
 * V0=Vreset */
double getSelfSpikeTime(struct neuron_helper_t *params);

/* Initialize a LIF neuron with exponential synapses */
neuron_state_t *InitExpLIFNeuron(unsigned long int me);

/* Connect presynaptic neurons to neuron my_id */
void ConnectPresynaptics(const unsigned long int my_id, const int my_population,
    const unsigned long int *pop_sizes, unsigned long int pop_count);

static struct neuron_helper_t n_helper_p[2];

/* NETWORK PARAMETERS v*/
//~ Probability connection table
double table[2][2] = {{0.02, 0.02}, {0.02, 0.02}};
/* END NETWORK PARAMETERS ^ */

FILE *outFile = NULL;

enum { OPT_SCALING = 129, OPT_TAU_M };

struct ap_option model_options[] = {
    {"scaling", OPT_SCALING, "VALUE",
	"Scaling factor applied to population connection probability."
	" conn_prob = conn_prob * scaling.\n"},
    {"tau_m", OPT_TAU_M, "VALUE",
	"Custom value for membrane time constant. Default: 20.0.\n"},
    {0}};

void model_parse(int key, const char *arg)
{
	switch(key) {
	case AP_KEY_FINI:
		break;
	case AP_KEY_INIT:
		break;
	case OPT_SCALING: {
		double conn_scaling;
		if(sscanf(arg, "%lf", &conn_scaling) != 1) {
			printf("Could not parse scaling option\n");
			abort();
		}
		if(conn_scaling <= 0.0) {
			printf("Option scaling has to be positive.\n");
			abort();
		}
		for(int i = 0; i < 2; ++i) {
			for(int j = 0; j < 2; ++j) {
				table[i][j] *= conn_scaling;
			}
		}
		break;
	}
	case OPT_TAU_M: {
		double tau_m;
		if(sscanf(arg, "%lf", &tau_m) != 1) {
			printf("Could not parse tau_m option\n");
			abort();
		}
		if(tau_m <= 0.0) {
			printf("Option tau_m has to be positive.\n");
			abort();
		}
		n_params.inv_tau_m = 1.0 / tau_m;
		break;
	}
	default: {
		printf("Argument not recognized\n");
		abort();
	}
	}
}

static const double exp_adj[] = {1.036199565629158, 1.027277981818330,
    1.019341151369058, 1.012316468833528, 1.006139353234098, 1.000751255315155,
    0.996100136699148, 0.992140185571610, 0.988830196962552, 0.986132544047804,
    0.984013562521051, 0.982443170190283, 0.981394102987059, 0.980841516073624,
    0.980763042138659, 0.981138430277712, 0.981957925989939, 0.983178546312421,
    0.984811614558841, 0.986834419454900, 0.989234236856811, 0.991999812972351,
    0.995121060808105, 0.998588377827231, 1.002392956865277, 1.006527340270889,
    1.010984919213803, 1.015759089560808, 1.020843804511367, 1.026234345941742,
    1.031901110282296, 1.037916082604387};

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

	neuron_state_t *state = InitExpLIFNeuron(me);

	// FIXME: This is obviously NOT the way to fill these fields. Just for
	// testing purposes.
	unsigned long int e_count = e_portion * n_lps;
	unsigned long int i_count = n_lps - e_count;
	unsigned long int pop_sizes[2] = {e_count, i_count};

	int my_pop = n2pop(me);
	if(populations[my_pop].record_spikes) NewProbe(me);

	// Init topology now
	ConnectPresynaptics(me, my_pop, pop_sizes, POPULATIONS_COUNT);

	return state;
}


neuron_state_t *InitExpLIFNeuron(unsigned long int me)
{
	neuron_state_t *state = malloc(sizeof(neuron_state_t));

	int pop = n2pop(me);

	struct neuron_helper_t *n_helper = &n_helper_p[pop];

	if(!n_helper->self_spike_time) {
		// INIZIALIZZA TUTTO DENTRO ST_PARAMS
		n_helper->El = g_El;
		n_helper->Q = n_helper->El * n_params.inv_tau_m *
		              n_params.inv_tau_m; // El/(tau_m^2)
		printdbg("El: %lf, Q: %lf\n", n_helper_p->El, n_helper_p->Q);

		n_helper->self_spike_time = getSelfSpikeTime(n_helper);
		printf("Neuron population %d self spike time: %lf\n", pop,
		    n_helper->self_spike_time);
	}

	state->helper = n_helper;
	state->last_updated = 0;
	state->Ge = 0.0;
	state->Gi = 0.0;
	state->last_fired = -(n_params.refractory_period + 1);


	state->membrane_potential =
	    n_params.reset_potential +
	    Random() * (n_params.threshold - n_params.reset_potential);

	return state;
}

/* neuron receives spike */
void NeuronHandleSpike(unsigned long int me, simtime_t now, double value,
    void *neuron_state)
{
	printdbg("[N%lu] Spike received at %lf\n", me, now);

	neuron_state_t *state = (neuron_state_t *)neuron_state;

	simtime_t fireTime;

	printdbg("[N%lu] Spike received with value: %lf at %lf\n", me, value,
	    now);
	printNeuronState(state);

	simtime_t delta_spike = now - state->last_fired;

	// This brings the neuron to present. Computes and updates Gi, Ge and V.
	bring_to_present(state, now - state->last_updated, delta_spike);
	// Apply the spike
	if(value > 0.0) {
		state->Ge += value;
	} else {
		state->Gi += value;
	}

	state->last_updated = now;

	printdbg("[N%lu] updated values: ", me);
	printNeuronState(state);

	/* Now we have to compute the point in time - if any - at which this
	 * neuron will spike with the current values */
	fireTime = getNextFireTime(state, delta_spike);
	printdbg("[N%lu] Next spike time is %lf\n", me, fireTime);

	if(fireTime > 0.0) { // Spike might happen in future
		MaybeSpikeAndWake(me, fireTime);
	}
}


/*  Neuron has been woken by the MaybeSpikeAndWake */
void NeuronWake(unsigned long int me, simtime_t now, void *n_state)
{
	printdbg(
	    "[N%lu] Woken up from spiking at %lf**********************************************************************!\n",
	    me, now);


	neuron_state_t *state = n_state;

	simtime_t t_since_last_eval = now - state->last_updated;
	state->Gi = get_Gi_f(t_since_last_eval, state->Gi);
	state->Ge = get_Ge_f(t_since_last_eval, state->Gi);
	state->membrane_potential = n_params.reset_potential;
	state->last_fired = now;
	state->last_updated = now;

	/* Now we have to compute the point in time - if any - at which this
	 * neuron will spike with the current values */
	double fireTime = getNextFireTime(state, now - state->last_fired);
	printdbg("[N%lu] Next spike time is %lf\n", me, fireTime);

	if(fireTime > 0.0) { // Spike might happen in future
		MaybeSpikeAndWake(me, fireTime);
	}
}


/* synapse handles the spike by updating its state and *returning* the spike
 * intensity */
double SynapseHandleSpike(simtime_t now, unsigned long int src_neuron,
    unsigned long int dest_neuron, synapse_t *state)
{
	//~ printdbg("[Synapse] Spike directed from N %lu to N %lu with value:
	//%lf\n", src_neuron, dest_neuron, state->weight);
	(void)now;
	(void)src_neuron;
	(void)dest_neuron;
	return state->weight;
}


/* is the neuron done? */
bool NeuronCanEnd(unsigned long int me, neuron_state_t *state)
{
	(void)me;
	(void)state;
	return false;
}


void ProbeRead(simtime_t now, unsigned long int monitored_neuron,
    const neuron_state_t *neuron_state)
{
	(void)now;
	(void)monitored_neuron;
	(void)neuron_state;
}


void GatherStatistics(simtime_t now, unsigned long int neuron_id,
    const neuron_state_t *state)
{
	(void)now;
	(void)state;

	// TODO: Fill this in to gather statistics and print them to file
	if(outFile == NULL) { // First neuron to write. Open file and write
		char *fname = malloc(strlen("CUBARun_") + 50);
		sprintf(fname, "CUBARun_%lu", neuron_id);
		outFile = fopen(fname, "w");
		free(fname);
	}
}


void SNNInitTopology(unsigned long int neuron_count)
{
	(void)neuron_count;
}

void ConnectPresynaptics(const unsigned long int my_id, const int my_population,
    const unsigned long int *pop_sizes, unsigned long int pop_count)
{
	unsigned long int cumulative_size = 0;

	unsigned long int src_neuron;
	synapse_t *synapse;

	double we = g_we;
	double wi = g_wi;
	double delay = g_syn_delay;

	unsigned long int total_synapses = 0;

	// For each presynaptic population
	for(unsigned int pre_pop = 0; pre_pop < pop_count; pre_pop++) {
		unsigned long int pre_size = pop_sizes[pre_pop];
		// connection_probability_table[pre_synaptic][post_synaptic] =
		// connection probability from pre to post
		double p = connection_probability_table[pre_pop][my_population];
		double weight = (is_excitatory(pre_pop) ? we : wi);

		// Then decide how many presynaptic neurons from that population
		// this neuron has -> binomial distribution
		unsigned long int n_pre = random_binomial(pre_size, p);
		total_synapses += n_pre;

		// Connect all the presynaptic neurons to this one. (remember to
		// sum the offset!)
		for(unsigned long int i = 0; i < n_pre; i++) {
			// FIXME: Does this never extract the last neuron of a
			// population?
			src_neuron =
			    ((unsigned long int)(Random() * pre_size)) +
			    cumulative_size;
			synapse = NewSynapse(src_neuron, my_id,
			    sizeof(synapse_t), true, delay);
			synapse->weight = weight;
		}

		cumulative_size += pre_size;
	}

	printdbg("Neuron %llu has %llu incoming connections\n", my_id,
	    total_synapses);
}

bool is_excitatory(unsigned int pre_pop)
{
	return populations[pre_pop].is_excitatory;
}

// From Luc Devroye's book "Non-Uniform Random Variate Generation." p. 522
unsigned random_binomial(unsigned trials, double p)
{
	if(p >= 1.0 || !trials) {
		return trials;
	}
	unsigned x = 0;
	double sum = 0,
	       log_q = log(1.0 - p); // todo cache those logarithm value
	while(1) {
		double r = Random();
		sum += log(r) / (trials - x);
		if(sum < log_q || trials == x) {
			return x;
		}
		x++;
	}
}

/* Generates an array a s.t. a[i] contains the index of the first neuron of
 * population i */
void gen_indexes(unsigned int *popsizes, unsigned int *out, int size)
{
	(void)popsizes;

	if(size != 2)
		return;

	out[0] = 0;
	out[1] = (unsigned int)e_portion * n_lps;
}


/* Get population index from neuron ID */
int n2pop(unsigned long int neuron_ID){
	for(int i=0; i<POPULATIONS_COUNT; i++){
		if(neuron_ID < population_sizes[i]){
			return i;
		}
		neuron_ID -= population_sizes[i];
	}
	fprintf(stderr, "ID %lu is too big. The neuron belongs to no declared population.\n", neuron_ID);
	abort();
}


/* Called at the start of the spike handling. Updates the state of the neuron
 * from the last update to now. We are certain that the neuron did not spike
 * from the last time it was updated to now. So we don't need to do calculations
 * to account for that */
void bring_to_present(neuron_state_t *state, simtime_t delta_update,
    simtime_t delta_spike)
{
	double Gef = get_Ge_f(delta_update, state->Ge);
	double Gif = get_Gi_f(delta_update, state->Gi);

	double Gent = state->Ge;
	double Gint = state->Gi;

	state->Ge = Gef;
	state->Gi = Gif;

	if(delta_spike < n_params.refractory_period)
		return;

	// Amount of time of refractory period that we still have to take into
	// account
	double remaining =
	    delta_update + n_params.refractory_period - delta_spike;

	if(remaining > 0.0) { // computes currents at end of refractory period
		Gent = get_Ge_f(remaining, state->Ge);
		Gint = get_Gi_f(remaining, state->Gi);
		// Time between refractory period and now
		delta_update -= remaining;
	}

	state->membrane_potential =
	    get_V_t(state->helper, state->membrane_potential, Gent, state->Ge,
		Gint, state->Gi, delta_update);
}


// Only takes into account neurons that are self spiking.
simtime_t getNextFireTime(neuron_state_t *state, simtime_t delta_spike)
{
	struct neuron_helper_t *n_helper = state->helper;

	double Ge0 = state->Ge;
	double Gi0 = state->Gi;

	simtime_t t0 = state->last_updated;
	double V0 = state->membrane_potential;

	simtime_t delta_t;
	if(n_params.refractory_period > delta_spike) {
		// Time till end of refractory period
		delta_t = n_params.refractory_period - delta_spike;
		// Currents at end of refractory period
		Ge0 = get_Ge_f(delta_t, state->Ge);
		Gi0 = get_Gi_f(delta_t, state->Gi);
	} else {
		delta_t = 0;
	}

	assert(n_helper->El > n_params.threshold); // Neuron is self spiking
	return t0 + findSpikeDeltaBinaryBlind(n_helper, V0, Ge0, Gi0) + delta_t;
}


inline double get_Ge_f(double delta_t, double Ge_i)
{
	return csexp(-delta_t * n_params.inv_tau_e) * Ge_i;
}

inline double get_Gi_f(double delta_t, double Gi_i)
{
	return csexp(-delta_t * n_params.inv_tau_i) * Gi_i;
}


// Might want to make another function that computes everything in one, given
// the state.
inline double get_V_t(struct neuron_helper_t *helper, double V0, double Ge0,
    double Ge1, double Gi0, double Gi1, double delta_t)
{
	(void)Ge1;
	(void)Gi1;
	/* Now compute V(delta_t) and return it */
	double aux = (V0 - Ge0 * n_params.De - Gi0 * n_params.Di - helper->Q) *
	             csexp(-delta_t * n_params.inv_tau_m);

	double expe = csexp(-delta_t * n_params.inv_tau_e);
	double expi = csexp(-delta_t * n_params.inv_tau_i);

	return (Ge0 * expe * n_params.De + Gi0 * expi * n_params.Di +
		helper->Q + aux);
}


// Might want to make another function that computes everything in one, given
// the state.
inline double get_V_t_acc(struct neuron_helper_t *helper, double V0, double Ge0,
    double Ge1, double Gi0, double Gi1, double delta_t)
{
	(void)Ge1;
	(void)Gi1;

	/* Now compute V(delta_t) and return it */
	double aux = (V0 - Ge0 * n_params.De - Gi0 * n_params.Di - helper->Q) *
	             exp(-delta_t * n_params.inv_tau_m);

	double expe = exp(-delta_t * n_params.inv_tau_e);
	double expi = exp(-delta_t * n_params.inv_tau_i);

	return (Ge0 * expe * n_params.De + Gi0 * expi * n_params.Di +
		helper->Q + aux);
}

extern double get_V_t_acc(struct neuron_helper_t *helper, double V0, double Ge0,
    double Ge1, double Gi0, double Gi1, double delta_t);


/* Given a neuron, compute how long it takes to spike without external
 * time-dependent inputs */
double getSelfSpikeTime(struct neuron_helper_t *n_helper)
{
	// Check that dV/dt > 0 for V->Vth and ge=gi=0
	if(n_params.threshold > n_helper->El)
		return -1;

	double delta_t;

	double Vt;
	double tmin = 0;
	double tmax = 50;

	// While Vt is below Vth
	while(1) {
		/* Now compute V(delta_t) */
		Vt = get_V_t_acc(n_helper, n_params.reset_potential, 0, 0, 0, 0,
		    tmax);

		if(Vt >= n_params.threshold) {
			//~ printdbg("delta_t: %lf, Vt %lf > Vth %lf\n", tmax,
			// Vt, n_params.threshold);
			break;
		}
		//~ printdbg("delta_t: %lf, Vt %lf < Vth %lf\n", tmax, Vt,
		// n_params.threshold);
		tmin = tmax;
		tmax *= 2;
	} // We now have a cap for the spike time

	// While V is farther than epsilon from Vth. Binary search
	do {
		//~ printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		delta_t = (tmax + tmin) / 2;
		/* Now compute V(delta_t) */
		Vt = get_V_t_acc(n_helper, n_params.reset_potential, 0, 0, 0, 0,
		    delta_t);


		if(Vt > n_params.threshold) {
			tmax = delta_t;
		} else {
			tmin = delta_t;
		}
	} while(tmax > tmin + T_TOLERANCE);

	printdbg("SelfSpikeTime spikedelta: %lf\n", delta_t);
	return delta_t;
}


/* Compute the time after which the spike will take place given T0=0 and V0
 * on a self spiking neuron */
double findSpikeDeltaBinaryBlind(struct neuron_helper_t *helper, double V0,
    double Ge0, double Gi0)
{
	double Vt;
	double Ge_t;
	double Gi_t;
	double tmin = 0;
	double tmax = helper->self_spike_time * 2;

	// While Vt is below Vth
	printdbg("FSDBB: Neuron is self spiking\n");
	while(true) {
		/* Compute Ge(delta_t) and Gi(delta_t) */
		Ge_t = get_Ge_f(tmax, Ge0);
		Gi_t = get_Gi_f(tmax, Gi0);
		/* Now compute V(delta_t) */
		Vt = get_V_t(helper, V0, Ge0, Ge_t, Gi0, Gi_t, tmax);

		if(Vt >= n_params.threshold) {
			//~ printdbg("FSDBB: delta_t: %lf, Vt %lf > Vth %lf\n",
			// delta_t, Vt, n_params.threshold);
			break;
		}
		//~ printdbg("FSDBB: delta_t: %lf, Vt %lf < Vth %lf\n", delta_t,
		// Vt, n_params.threshold);
		tmin = tmax;
		tmax *= 2;
	}

	// While V is farther than epsilon from Vth
	do {
		//~ printdbg("tmin: %lf, tmax: %lf\n", tmin, tmax);
		double delta_t = (tmax + tmin) / 2;

		/* Compute Ge(delta_t) and Gi(delta_t) */
		Ge_t = get_Ge_f(delta_t, Ge0);
		Gi_t = get_Gi_f(delta_t, Gi0);
		/* Now compute V(delta_t) */
		Vt = get_V_t(helper, V0, Ge0, Ge_t, Gi0, Gi_t, delta_t);

		if(Vt > n_params.threshold) {
			//~ printdbg("Vt %lf > Vth %lf\n", Vt,
			// n_params.threshold);
			tmax = delta_t;
		} else {
			//~ printdbg("Vt %lf < Vth %lf\n", Vt,
			// n_params.threshold);
			tmin = delta_t;
		}
		if(tmax <= tmin + T_TOLERANCE) {
			printdbg("FSDBB spikedelta: %lf\n", delta_t);
			return delta_t;
		}
	} while(1);
}

void printNeuronState(neuron_state_t *state)
{
	(void)state;
	printdbg("V: %lf, Ge: %lf, Gi: %lf, last_up: %lf, last_spike: %lf\n",
	    state->membrane_potential, state->Ge, state->Gi,
	    state->last_updated, state->last_fired);
}
