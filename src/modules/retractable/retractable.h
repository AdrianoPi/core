#pragma once

#include <datatypes/retractable_heap.h>
#include <lp/msg.h>
#include <lib/lib_internal.h>
#include <lib/lib.h>

#define rq_elem_update(rq, i) ((rq).m->pos = (i))

struct rq_elem {
	simtime_t t;
	struct lp_msg *m;
};

typedef rheap_declare(struct rq_elem) r_queue_t;

extern __thread r_queue_t r_queue;

void retractable_module_init(void);
void retractable_module_fini(void);

void retractable_module_lp_fini(void);
void retractable_module_lp_init(void);

extern void retractable_msg_schedule(simtime_t timestamp, unsigned event_type);

inline bool is_retractable(const struct lp_msg* msg)
{
	return atomic_load_explicit(&msg->flags, memory_order_relaxed) & MSG_FLAG_RETRACTABLE;
}

inline bool is_valid_retractable(const struct lp_msg* msg)
{
	return msg->dest_t == current_lp->lib_ctx_p->r_ts;
}

extern void retractable_rollback_handle(void);
extern void ScheduleRetractableEvent_pr(simtime_t timestamp, unsigned event_type);

inline void ScheduleRetractableEvent(simtime_t timestamp, unsigned event_type)
{
	ScheduleRetractableEvent_pr(timestamp, event_type);
}

inline void DescheduleRetractableEvent(void)
{
	// Just set the timestamp in LP memory to be invalid
	current_lp->lib_ctx_p->r_ts = -1.0;
}
