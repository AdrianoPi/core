#ifndef CORE_MODEL_TOPOLOGY_H
#define CORE_MODEL_TOPOLOGY_H

#include <stdbool.h>
#include <ROOT-Sim.h>

#define POPULATIONS_COUNT 2

typedef struct neuron_params_t{
	double inv_tau_m; // [1/ms]
	double inv_tau_e; // [1/ms]
	double inv_tau_i; // [1/ms]
	double De; // Inverse of (1-tau_m/tau_e)
	double Di; // Inverse of (1-tau_m/tau_i)
	double reset_potential; // [mV]
	double threshold; // [mV]
	double refractory_period; // [ms]
} neuron_params_t;

struct neural_population {
	struct neuron_params_t *parameters;
	unsigned population_id; // TODO: Remove?
	unsigned size; // TODO: Remove if we move this info to population_sizes?
	bool is_excitatory;
	bool record_spikes;
};

// connection_probability_table[pre][post] Contains the connection probability from pre to post
extern double connection_probability_table[POPULATIONS_COUNT][POPULATIONS_COUNT];
extern struct neural_population populations[POPULATIONS_COUNT];
extern unsigned population_sizes[POPULATIONS_COUNT];
extern struct simulation_configuration global_config_from_file;

#endif // CORE_MODEL_TOPOLOGY_H
