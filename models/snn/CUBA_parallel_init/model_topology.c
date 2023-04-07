#include "model_topology.h"
#include <stdbool.h>

#define POPULATIONS_COUNT 2

// connection_probability_table[pre][post] Contains the connection probability from pre to post
double connection_probability_table[POPULATIONS_COUNT][POPULATIONS_COUNT] = {
	{0.02, 0.02},
	{0.02, 0.02}
};

struct neural_population populations[POPULATIONS_COUNT] = {
    {.parameters =
	    {
		.inv_tau_m = 1 / 20.0,             // [1/ms]
		.inv_tau_e = 1 / 5.0,              // [1/ms]
		.inv_tau_i = 1 / 10.0,             // [1/ms]
		.De = 1.0 / (1.0 - (20.0 / 5.0)),  // Inverse of (1-tau_m/tau_e)
		.Di = 1.0 / (1.0 - (20.0 / 10.0)), // Inverse of (1-tau_m/tau_i)
		.reset_potential = -60.0,          // [mV]
		.threshold = -50.0,                // [mV]
		.refractory_period = 5.0           // [ms]
	    },
	.is_excitatory = true,
	.record_spikes = true},
    {.parameters =
	    {
		.inv_tau_m = 1 / 20.0,             // [1/ms]
		.inv_tau_e = 1 / 5.0,              // [1/ms]
		.inv_tau_i = 1 / 10.0,             // [1/ms]
		.De = 1.0 / (1.0 - (20.0 / 5.0)),  // Inverse of (1-tau_m/tau_e)
		.Di = 1.0 / (1.0 - (20.0 / 10.0)), // Inverse of (1-tau_m/tau_i)
		.reset_potential = -60.0,          // [mV]
		.threshold = -50.0,                // [mV]
		.refractory_period = 5.0           // [ms]
	    },
	.is_excitatory = false,
	.record_spikes = false}};

unsigned population_sizes[POPULATIONS_COUNT] = {800, 200};

struct simulation_configuration global_config_from_file = {
	.termination_time = 10,
	.gvt_period = 250,
	.verbosity = 1,
	.prng_seed = 1,
	.core_binding = true,
	.is_serial = false
};
