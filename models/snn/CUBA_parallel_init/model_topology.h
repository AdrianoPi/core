#ifndef CORE_MODEL_TOPOLOGY_H
#define CORE_MODEL_TOPOLOGY_H

#include <stdbool.h>
#include <ROOT-Sim.h>
#include "cuba.h"

#define POPULATIONS_COUNT 2

typedef struct neuron_params_t neuron_params_t;

struct neural_population {
	struct neuron_params_t *parameters;
	unsigned population_id;
	unsigned size;
	bool is_excitatory;
};

// connection_probability_table[pre][post] Contains the connection probability from pre to post
extern double connection_probability_table[POPULATIONS_COUNT][POPULATIONS_COUNT];
extern struct neural_population populations[POPULATIONS_COUNT];
extern struct simulation_configuration global_config_from_file;

#endif // CORE_MODEL_TOPOLOGY_H
