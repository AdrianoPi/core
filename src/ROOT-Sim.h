/**
 * @file ROOT-Sim.h
 *
 * @brief ROOT-Sim header for model development
 *
 * This header defines all the symbols which are needed to develop a model
 * to be simulated on top of ROOT-Sim.
 *
 * This header is the only file which should be included when developing
 * a simulation model. All function prototypes exposed to the application
 * developer are exposed and defined here.
 *
 * SPDX-FileCopyrightText: 2008-2021 HPDCS Group <rootsim@googlegroups.com>
 * SPDX-License-Identifier: GPL-3.0-only
 */
#pragma once

#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

enum rootsim_event {
	MODEL_FINI = 65532,
	LP_FINI,
	LP_INIT,
	MODEL_INIT
};

typedef double simtime_t;
typedef uint64_t lp_id_t;

struct ap_option {
	const char *name;
	int key;
	const char *arg;
	const char *doc;
};

enum ap_event_key {
	AP_KEY_INIT = 1 << 14,
	AP_KEY_FINI
};

__attribute((weak)) extern struct ap_option model_options[];
__attribute((weak)) extern void model_parse(int key, const char *arg);

extern lp_id_t n_lps;

extern void ScheduleNewEvent(lp_id_t receiver, simtime_t timestamp,
	unsigned event_type, const void *event_content, unsigned event_size);

#ifdef PUBSUB
// For publish-subscribe messages
extern void ProcessPublishedEvent(lp_id_t me, simtime_t msg_ts, unsigned int event, const void* msg_content, unsigned int size, const void* user_data);
/// To be called from within ProcessPublishedEvent to complete the unpacking of the message.
extern void PubsubDeliver(simtime_t timestamp, unsigned event_type, const void *event_content, unsigned event_size);
extern void PublishNewEvent(simtime_t timestamp, unsigned event_type, const void *event_content, unsigned event_size);
extern void Subscribe(lp_id_t subscriber_id, lp_id_t publisher_id);
extern void SubscribeAndHandle(lp_id_t subscriber_id, lp_id_t publisher_id, void* data);
#endif
// TODO: readd guards
// For retractable events
extern void ScheduleRetractableEvent(simtime_t timestamp, unsigned event_type);
extern void DescheduleRetractableEvent();
// </readd guards>

extern void SetState(void *new_state);

extern double Random(void);
extern uint64_t RandomU64(void);
extern double Expent(double mean);
extern double Normal(void);
extern int RandomRange(int min, int max);

enum _topology_geometry_t {
	TOPOLOGY_HEXAGON = 1,	//!< hexagonal grid topology
	TOPOLOGY_SQUARE,	//!< square grid topology
	TOPOLOGY_RING,		//!< a ring shaped topology walkable in a single direction
	TOPOLOGY_BIDRING,	//!< a ring shaped topology
	TOPOLOGY_TORUS,		//!< a torus shaped grid topology (a wrapping around square topology)
	TOPOLOGY_STAR,		//!< a star shaped topology
	TOPOLOGY_MESH,		//!< an arbitrary shaped topology
};

enum _direction_t {
	DIRECTION_N,	//!< North direction
	DIRECTION_S,	//!< South direction
	DIRECTION_E,	//!< East direction
	DIRECTION_W,	//!< West direction
	DIRECTION_NE,	//!< North-east direction
	DIRECTION_SW,	//!< South-west direction
	DIRECTION_NW,	//!< North-west direction
	DIRECTION_SE,	//!< South-east direction

	// FIXME this is bad if the n_lps is more than INT_MAX - 1
	DIRECTION_INVALID = INT_MAX	//!< A generic invalid direction
};

struct topology_t;

extern struct topology_t *TopologyInit(enum _topology_geometry_t geometry, unsigned int out_of_topology);
extern unsigned long long RegionsCount(const struct topology_t *topology);
extern unsigned long long DirectionsCount(const struct topology_t *topology);
extern lp_id_t GetReceiver(const struct topology_t *topology, lp_id_t from, enum _direction_t direction);
extern lp_id_t FindReceiver(const struct topology_t *topology);

/// A set of configurable values used by other modules
struct simulation_configuration {
	/// The target termination logical time
	simtime_t termination_time;
	/// The gvt period expressed in microseconds
	unsigned gvt_period;
	/// The log verbosity level
	int verbosity;
	/// The seed used to initialize the pseudo random numbers
	uint64_t prng_seed;
	/// If set, worker threads are bound to physical cores
	bool core_binding;
	/// If set, the simulation will run on the serial runtime
	bool is_serial;
};

extern struct simulation_configuration global_config;
