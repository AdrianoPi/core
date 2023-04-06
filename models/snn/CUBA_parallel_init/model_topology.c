#include "model_topology.h"
#include <stdbool.h>

#define POPULATIONS_COUNT 2

// connection_probability_table[pre][post] Contains the connection probability from pre to post
double connection_probability_table[POPULATIONS_COUNT][POPULATIONS_COUNT] = {
	{0.02, 0.02},
	{0.02, 0.02}
};

struct neural_population populations[POPULATIONS_COUNT] = {
	{
		.parameters = NULL,
		.is_excitatory = true,
		.record_spikes = true
	},
	{
		.parameters = NULL,
		.is_excitatory = false,
		.record_spikes = false
	}
};

unsigned population_sizes[POPULATIONS_COUNT] = {4200, 800};

struct simulation_configuration global_config_from_file = {
	.termination_time = 10,
	.gvt_period = 250,
	.verbosity = 1,
	.prng_seed = 1,
	.core_binding = true,
	.is_serial = false
};
