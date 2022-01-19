#include <modules/retractable/retractable.h>

__thread r_queue_t r_queue;

void retractable_module_init(void)
{
	rheap_init(r_queue);
}

void retractable_module_fini(void)
{
	rheap_fini(r_queue);
}

void retractable_module_lp_init()
{
	current_lp->r_msg = NULL;

	current_lp->lib_ctx_p->r_ts = -1.0;
	current_lp->lib_ctx_p->r_e_type = -1;
}

void retractable_module_lp_fini()
{
	struct lp_msg *msg = current_lp->r_msg;
	// If there is a message and it is not scheduled
	// (otherwise it gets freed by msg_queue_fini)
	if(msg && msg->dest_t < 0)
		msg_allocator_free(current_lp->r_msg);
}

/* This function populates the retractable message with the correct info
 * taken from the LP's lib state managed and then schedules it/moves it
 * to the right position in the input queue */
void msg_queue_insert_retractable()
{
	if(current_lp->lib_ctx_p->r_ts < 0) // The message is not to be scheduled
		return;

	// The message is to be scheduled
	struct lp_msg* msg = current_lp->r_msg;

	if(msg == NULL) { // Need to create a new retractable msg
		msg = msg_allocator_pack(current_lp - lps, current_lp->lib_ctx_p->r_ts,
			current_lp->lib_ctx_p->r_e_type, NULL, 0);
#if LOG_LEVEL <= LOG_DEBUG
		msg->send = current_lp - lps;
		msg->send_t = current_lp->p.last_t;
#endif
		msg->raw_flags = MSG_FLAG_RETRACTABLE;
		current_lp->r_msg = msg;

		struct rq_elem q_el = (struct rq_elem){.t = msg->dest_t, .m = msg};
		rheap_insert(r_queue, q_elem_is_before, rq_elem_update, q_el);
		return;
	}

	// Msg already exists

	// Is the message already in the incoming queue?
	bool already_in_Q = (msg->dest_t >= 0);

	// Set the correct values for the message
	bool lowered = current_lp->lib_ctx_p->r_ts > msg->dest_t;

	msg->dest_t = current_lp->lib_ctx_p->r_ts;
	msg->m_type = current_lp->lib_ctx_p->r_e_type;
#if LOG_LEVEL <= LOG_DEBUG
	msg->send_t = current_lp->p.last_t;
#endif

	// Schedule the message
	if(already_in_Q) {
		array_get_at(r_queue, msg->pos).t = msg->dest_t;
		struct rq_elem q_el = array_get_at(r_queue, msg->pos);

		if (lowered)
			rheap_priority_lowered(r_queue, q_elem_is_before, rq_elem_update, q_el);
		else
			rheap_priority_increased(r_queue, q_elem_is_before, rq_elem_update, q_el);
	} else {
		struct rq_elem q_el = (struct rq_elem){.t = msg->dest_t, .m = msg};
		rheap_insert(r_queue, q_elem_is_before, rq_elem_update, q_el);
	}
}

void retractable_rollback_handle(void)
{
	msg_queue_insert_retractable();
}

void ScheduleRetractableEvent_pr(simtime_t timestamp, unsigned event_type)
{
	assert(timestamp >= current_lp->p.last_t);

	current_lp->lib_ctx_p->r_ts = timestamp;
	current_lp->lib_ctx_p->r_e_type = event_type;

	if(silent_processing)
		return;

	msg_queue_insert_retractable();
}

extern bool is_retractable(const struct lp_msg* msg);
extern bool is_valid_retractable(const struct lp_msg* msg);
extern void DescheduleRetractableEvent(void);
extern void ScheduleRetractableEvent(simtime_t timestamp, unsigned event_type);
