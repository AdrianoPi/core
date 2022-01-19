#include <modules/publish_subscribe/pubsub.h>

#ifdef PUBSUB

#include <limits.h>		// for CHAR_BIT
#include <lp/process.h>

#define n_pub_info_ptr(msg) ((struct n_pubsub_info *)((msg)->pl + (msg)->pl_size - sizeof(struct n_pubsub_info)))
#define t_pub_info_ptr(msg) ((struct t_pubsub_info *)((msg)->pl + (msg)->pl_size - sizeof(struct t_pubsub_info)))

#define LP_ID_MSB (((lp_id_t) 1) << (sizeof(lp_id_t) * CHAR_BIT - 1))

#define PUBSUB_PRINT_TOPOLOGY true

#if LOG_LEVEL <= LOG_DEBUG
// File for logging pubsub messages
__thread FILE *pubsub_msgs_logfile;
#endif

#define mark_as_thread_lv(msg_p) ((struct lp_msg *)(((uintptr_t)(msg_p)) | 1U))
#define is_thread_lv(msg_p) (((uintptr_t)(msg_p)) & 1U)
#define unmark_msg(msg_p) \
	((struct lp_msg *)(((uintptr_t)(msg_p)) & (UINTPTR_MAX - 3)))

// Holds thread-level pubsub messages
__thread dyn_array(struct lp_msg *) past_pubsubs;
__thread dyn_array(struct lp_msg *) early_pubsub_antis;

typedef struct table_lp_entry_t{
	lp_id_t lid;
	void* data;
} table_lp_entry_t;

typedef dyn_array(table_lp_entry_t) lp_entry_arr;

typedef struct table_thread_entry_t{
	rid_t tid;
	lp_entry_arr lp_arr;
} table_thread_entry_t;

typedef dyn_array(table_thread_entry_t) t_entry_arr;

// Size of additional data needed by pubsub messages published locally
struct n_pubsub_info {
	size_t m_cnt;
	struct lp_msg **m_arr;
};

struct t_pubsub_info {
	size_t m_cnt;
	struct lp_msg **m_arr;
	lp_entry_arr *lp_arr_p;
};

t_entry_arr *subscribersTable;
spinlock_t *tableLocks;

/// Adds to past_pubsubs maintaining ordering
static void pubsub_insert_in_past(struct lp_msg *msg)
{
	simtime_t time = msg->dest_t;

	array_count_t i = array_count(past_pubsubs);
	if (!i || time >= array_get_at(past_pubsubs, i - 1)->dest_t) {
		array_push(past_pubsubs, msg);
		return;
	}

	do {
		--i;
	} while(i && time >= array_get_at(past_pubsubs, i - 1)->dest_t);

	array_add_at(past_pubsubs, i, msg);
}

static void thread_actually_antimessage(struct lp_msg *msg);
extern void pubsub_thread_msg_free(struct lp_msg* msg);

void pprint_subscribers_table();
void pprint_subscribers_adjacency_list();

// OK
void pubsub_module_lp_init()
{
	current_lp->subnodes = mm_alloc(bitmap_required_size(n_nodes));
	bitmap_initialize(current_lp->subnodes, n_nodes);
	current_lp->n_remote_sub_nodes = 0;
}

void pubsub_module_lp_fini(){
	mm_free(current_lp->subnodes);
}

// OK
/// Initializes structures needed for pubsub module
void pubsub_module_global_init(){
	// Right now, the hashtable is simply an array[n_lps]
	subscribersTable = mm_alloc(sizeof(t_entry_arr)*n_lps);
	tableLocks = mm_alloc(sizeof(spinlock_t)*n_lps);

	for(lp_id_t i=0; i<n_lps; i++){
		//~ subscribersTable[i].items = NULL;
		array_count(subscribersTable[i])=0; // Unorthodox. Works
		spin_init(&(tableLocks[i]));
	}
}

void print_pubsub_topology_to_file(){
	if(!PUBSUB_PRINT_TOPOLOGY) return;
	// Print two json files:
	// One is the subscribers table
	// The other is an array that at index i contains the list of subscribers of LP i.
	// Each MPI node prints its file containing its portion of the information
	FILE *f;
	char* fname = malloc(strlen("SubscribersTable_.json") + 100);
	sprintf(fname, "SubscribersTable_%d.json", nid);
	f = fopen(fname, "w");
	free(fname);
	pprint_subscribers_table(f);
	fclose(f);

	fname = malloc(strlen("SubscribersAdjacencyList_.json") + 100);
	sprintf(fname, "SubscribersAdjacencyList_%d.json", nid);
	f = fopen(fname, "w");
	free(fname);
	pprint_subscribers_adjacency_list(f);
	fclose(f);
}

/// Last call to the pubsub module
void pubsub_module_global_fini(){
#if LOG_LEVEL <= LOG_DEBUG
	log_log(LOG_DEBUG, "Starting pubsub_module_global_fini\n");

	// Free subscribersTable, its locks and its contents
	mm_free(tableLocks);
	for(lp_id_t pub_id=0; pub_id<n_lps; pub_id++){
		t_entry_arr t_arr = subscribersTable[pub_id];
		if(!array_count(t_arr)){
			continue;
		}
		// Free lp_arrs in t_arr.items
		for (array_count_t t_index = 0; t_index < array_count(t_arr); ++t_index) {
			lp_entry_arr lp_arr = array_get_at(t_arr, t_index).lp_arr;
			// free each allocated element of lp_arr.items
			for (array_count_t l_index = 0; l_index < array_count(lp_arr); ++l_index) {
				table_lp_entry_t lp_entry = array_get_at(lp_arr, l_index);

				if (lp_entry.data){
					mm_free(lp_entry.data);
				}
			}
			array_fini(lp_arr);
		}
		array_fini(t_arr);
	}
	mm_free(subscribersTable);

	log_log(LOG_DEBUG, "Ending pubsub_module_global_fini\n");
#endif
}

void pubsub_module_init(){
	array_init(past_pubsubs);
	array_init(early_pubsub_antis);
#if LOG_LEVEL <= LOG_DEBUG
	// Open logging file
	char* fname = malloc(strlen("pubsub_msgs_n_t.log") + 100);
	sprintf(fname, "pubsub_msgs_n%d_t%d.log", nid, rid);
	pubsub_msgs_logfile = fopen(fname, "w");
	setvbuf(pubsub_msgs_logfile, NULL, _IOFBF, 4096);
	fprintf(pubsub_msgs_logfile, "LP, send_time\n");
	free(fname);
#endif
}

#if LOG_LEVEL <= LOG_DEBUG
// FIXME: could this log a message twice if it is the last one of the array?
void log_pubsub_msgs_to_file(struct lp_msg** msg_array, array_count_t size, simtime_t max_time){
	// Since dyn_arrays are unnamed structs, we work with the items element
	for(array_count_t i = 0; i < size; i++){
		struct lp_msg* msg = unmark_msg(msg_array[i]);
		if (msg->dest_t > global_config.termination_time){
			return;
		}

		assert(!(msg->raw_flags & MSG_FLAG_ANTI));

		// Only log pubsubs that have not been undone
		if(is_pubsub_msg(msg)){
			fprintf(pubsub_msgs_logfile, "%lu, %lf\n", msg->send, msg->send_t);
		}
	}
}
#endif

void pubsub_module_fini(){
	for (array_count_t i = 0; i < array_count(past_pubsubs); ++i) {
		struct lp_msg* msg = array_get_at(past_pubsubs, i);
		if ((msg->raw_flags & MSG_FLAG_ANTI) && (msg->raw_flags & MSG_FLAG_PROCESSED)){
			// Is in msg_queue and will be freed there.
			continue;
		}
		pubsub_thread_msg_free(array_get_at(past_pubsubs, i));
	}
	array_fini(past_pubsubs);

	for (array_count_t i = 0; i < array_count(early_pubsub_antis); ++i) {
		msg_allocator_free(array_get_at(early_pubsub_antis, i));
	}
	array_fini(early_pubsub_antis);

#if LOG_LEVEL <= LOG_DEBUG
	fclose(pubsub_msgs_logfile);
#endif
}

/// This is called when a local LP publishes a message.
/// Sends message via MPI, and unpacks local copies.
void pub_node_handle_published_message(struct lp_msg* msg)
{
	/*
	 * The parent message keeps track of its children by keeping
	 * an array of pointers in the payload.
	 * Contents of msg->pl: [*(msg->pl) | size_t | lp_msg** ]
	 */

	t_entry_arr threads = subscribersTable[msg->dest];
	size_t stripped_payload_size = msg->pl_size - sizeof(struct n_pubsub_info);

	// One message per thread
	int n_ch_count = array_count(threads);

#ifdef ROOTSIM_MPI
	n_ch_count += current_lp->n_remote_sub_nodes;
#endif

	struct n_pubsub_info *pi = n_pub_info_ptr(msg);

	pi->m_cnt = n_ch_count;
	pi->m_arr = mm_alloc(sizeof(*pi->m_arr) * n_ch_count);

	// Index of next index to fill in n_children_ptr(msg)
	int it = 0;

#ifdef ROOTSIM_MPI
	// Create and send messages to other nodes
	int subbed_nodes = current_lp->n_remote_sub_nodes;
	if (subbed_nodes) {
		size_t stripped_msg_size = offsetof(struct lp_msg, pl) + stripped_payload_size;
		struct block_bitmap *subs = current_lp->subnodes;
		nid_t dest_nid = 0;
		do {
			while (!bitmap_check(subs, dest_nid))
				++dest_nid;

			struct lp_msg *mpi_msg = msg_allocator_alloc(
					stripped_payload_size);
			memcpy(mpi_msg, msg, stripped_msg_size);
			msg->pl_size = stripped_payload_size;

			n_pub_info_ptr(msg)->m_arr[it] = mpi_msg;
			++it;

			mpi_remote_msg_send(mpi_msg, dest_nid);
		} while (it < subbed_nodes);
	}

#endif

	// Size of the payload of thread-level messages
	size_t child_pl_size = stripped_payload_size + sizeof(struct t_pubsub_info);
	// Visualization of *child_payload once populated:
	// Offsets	:v-0       		v-original_pl_size
	// Contents	:[ 	*(msg->pl) 	| 	&lp_arr	|	childCount	| lp_msg**	]

	// Create and enqueue messages for each subscribed Thread
	for (int c = 0; it < n_ch_count; ++it, ++c) {

		table_thread_entry_t *t_entry = &array_get_at(threads, c);

		struct lp_msg *child_msg = msg_allocator_alloc(child_pl_size);
		memcpy(child_msg, msg, offsetof(struct lp_msg, pl_size));
		memcpy(child_msg->pl, msg->pl, stripped_payload_size);
		// Target is the target thread's tid
		child_msg->dest = t_entry->tid;
		child_msg->pl_size = child_pl_size;

		struct t_pubsub_info *pi = t_pub_info_ptr(child_msg);
		pi->lp_arr_p = &t_entry->lp_arr;
		pi->m_cnt = 0;
		pi->m_arr = NULL;

		// *child_payload right now:
		// Byte offsets	:	v-0       		v-user_pl_size
		// Contents	:		[ 	*(msg->pl) 	| 	lp_arr*	|	0	|	NULL	]
		// The info for pubsub will be filled out by the thread handling the messages

		child_msg->raw_flags = MSG_FLAG_PUBSUB;

		// Keep track of the child message for rollbacks!
		n_pub_info_ptr(msg)->m_arr[it] = child_msg;
		// Push child message into target thread's incoming queue
		pubsub_msg_queue_insert(child_msg);
	}
}

// Should be ok?
/// Called when MPI extracts a pubsub message
void sub_node_handle_published_message(struct lp_msg* msg)
{
	// On receiver nodes, the original message is only used to create thread-level ones

	t_entry_arr threads = subscribersTable[msg->dest];

	// One message per thread
	int n_ch_count = array_count(threads);

	assert(n_ch_count);

	size_t child_pl_size = msg->pl_size + sizeof(struct t_pubsub_info);
	// Visualization of child's payload once populated:
	// Offsets	:	v-0       		v-msg->pl_size
	// Contents	:	[ 	*(msg->pl) 	| 	lp_arr* |	childCount	| lp_msg**	]

	// Create and enqueue messages for each subscribed Thread
	for (int c = 0; c < n_ch_count; c++) {

		table_thread_entry_t *t_entry = &array_get_at(threads, c);

		// Create child message
		// Target holds the target thread's tid
		struct lp_msg *child_msg = msg_allocator_alloc(child_pl_size);
		memcpy(child_msg, msg, offsetof(struct lp_msg, pl_size));
		memcpy(child_msg->pl, msg->pl, msg->pl_size);
		// Target is the target thread's tid
		child_msg->dest = t_entry->tid;
		child_msg->pl_size = child_pl_size;

		struct t_pubsub_info *pi = t_pub_info_ptr(child_msg);
		pi->lp_arr_p = &t_entry->lp_arr;
		pi->m_cnt = 0;
		pi->m_arr = NULL;

		// *child_payload right now:
		// Byte offsets	:	v-0       		v-user_pl_size
		// Contents	:		[ 	*(msg->pl) 	| 	lp_arr*	|	0	|	NULL	]
		// The info for pubsub will be filled out by the thread handling the messages

		child_msg->raw_flags = MSG_FLAG_PUBSUB;

		// Push child message into target thread's incoming queue
		pubsub_msg_queue_insert(child_msg);

	}

	/* This msg is only a blueprint: Free it right now */
	msg_allocator_free(msg);
}

static inline bool check_early_anti_pubsub_messages(struct lp_msg *msg)
{
	uint32_t m_id;
	if (likely(!array_count(early_pubsub_antis) || !(m_id = msg->raw_flags >> MSG_FLAGS_BITS)))
		return false;
	uint32_t m_seq = msg->m_seq;
//	if (!m_id) // Already checked above
//		return false;
	for (array_count_t i = 0; i < array_count(early_pubsub_antis); ++i) {
		struct lp_msg *a_msg = array_get_at(early_pubsub_antis, i);
		if (a_msg->raw_flags == m_id && a_msg->m_seq == m_seq) {
			msg_allocator_free(msg);
			msg_allocator_free(a_msg);
			array_get_at(early_pubsub_antis, i) = array_peek(early_pubsub_antis);
			--array_count(early_pubsub_antis);
			return true;
		}
	}
	return false;
}

/* In ProcessPublishedEvent the user creates the message resulting from
 * the published one. The message is created using ScheduleNewEvent.
 * This will schedule the generated event and will push a pointer to it
 * in sent_msgs.
 * We pop that pointer from sent_msgs (it does not belong there) and add
 * it to the children of the msg we are unpacking here.
 * */
// This function is called when a thread extracts a pubsub message from its queue
void thread_handle_published_message(struct lp_msg* msg){

	if (unlikely(check_early_anti_pubsub_messages(msg)))
		return;

	uint32_t flags = atomic_fetch_add_explicit(&msg->flags, MSG_FLAG_PROCESSED, memory_order_relaxed);
	if (unlikely(flags & MSG_FLAG_ANTI)) {
		pubsub_thread_msg_free(msg);
		return;
	}

	/*
	 * The parent message keeps track of its children by keeping
	 * an array of pointers in the payload.
	 * Contents of msg->pl once populated:
	 * [ pl	| &lp_arr	| childCount	| lp_msg**	]
	 *
	 * Contents right now:
	 * [ pl	| &lp_arr	| 	0	| NULL		]
	 */

	struct t_pubsub_info *pi = t_pub_info_ptr(msg);
	lp_entry_arr lp_arr = *pi->lp_arr_p;

	assert(array_count(lp_arr));

	// One message per subscription
	pi->m_cnt = array_count(lp_arr);
	pi->m_arr = mm_alloc(sizeof(*pi->m_arr) * array_count(lp_arr));
	// Contents msg->pl now:
	// Byte offsets	:v-0    v-og_pl_size	v-(pl_size+sizeof(void*))
	// Contents	:[ pl	| &lp_arr	| childCount	| lp_msg**	]

	// Here calculate the size of the payload for the messages LPs will receive
	size_t original_pl_size = msg->pl_size - sizeof(struct t_pubsub_info);

	// For each subscribed LP
	for (array_count_t i = 0; i < array_count(lp_arr); i++) {
		struct lp_msg* child_msg;
		table_lp_entry_t c_lp_entry = array_get_at(lp_arr, i);

		lp_id_t target_lid = c_lp_entry.lid; // Still dirty
		bool justSub = target_lid & LP_ID_MSB;
		target_lid = target_lid & ~LP_ID_MSB; // Clear leftmost bit

		// Check: Subscribe or SubscribeAndHandle?
		if(justSub){ // Regular subscribe
			// Just enqueue a clone of the message!
			child_msg = msg_allocator_pack(
					target_lid,
					msg->dest_t,
					msg->m_type,
					msg->pl,
					original_pl_size
			);
#if LOG_LEVEL <= LOG_DEBUG
			child_msg->send = msg->send;
			child_msg->send_t = msg->send_t;
#endif

			atomic_store_explicit(&child_msg->flags, 0U, memory_order_relaxed);

			msg_queue_insert(child_msg);

		} else { // SubscribeAndHandle
			// Set the current LP to the target lp's id
			struct lp_ctx* this_lp = &lps[target_lid];
			current_lp = this_lp;

			struct process_data *proc_p = &current_lp->p;

			// Get user-provided data from entry
			void *usr_data = c_lp_entry.data;

			// > Invoke the handler with correct data and lp_id
			//ProcessPublishedEvent_pr(
			ProcessPublishedEvent(
				current_lp - lps,
				msg->dest_t,
				msg->m_type,
				msg->pl,
				original_pl_size,
				usr_data
			);

			// Flags and all initialized in ScheduleNewEvent

			// Pop the created message from sent_msgs
			child_msg = unmark_msg(array_pop(proc_p->p_msgs));
#if LOG_LEVEL <= LOG_DEBUG
			child_msg->send = msg->send;
			child_msg->send_t = msg->send_t;
#endif
		}

		// Keep track of the child message
		pi->m_arr[i] = child_msg;

	}
	pubsub_insert_in_past(msg);
}

// ok
inline void PublishNewEvent_pr(simtime_t timestamp, unsigned event_type, const void *payload, unsigned payload_size){
	PublishNewEvent(timestamp, event_type, payload, payload_size);
}

// ok?
// Here the logic that takes care to send the message to every node is implemented.
void PublishNewEvent(simtime_t timestamp, unsigned event_type, const void *payload, unsigned payload_size)
{
	if(silent_processing)
		return;

	struct process_data *proc_p = &current_lp->p;

	// A node-level pubsub message is created.
	struct lp_msg *msg = msg_allocator_alloc(payload_size + sizeof(struct n_pubsub_info));
	msg->dest = current_lp - lps;
	msg->dest_t = timestamp;
	msg->m_type = event_type;
#if LOG_LEVEL <= LOG_DEBUG
	msg->send = current_lp - lps;
	msg->send_t = proc_p->last_t;
#endif

	if (payload_size) {
		memcpy(msg->pl, payload, payload_size);
	}

	msg->raw_flags = 0;

	pub_node_handle_published_message(msg);

	// Do this after making copies, or it will dirty MPI messages!
	msg->raw_flags = MSG_FLAG_PUBSUB;
	array_push(proc_p->p_msgs, mark_msg_remote(msg));
}

// Ok?
/// This carries out the antimessaging when the node received the antimessage via MPI
void sub_node_handle_published_antimessage(struct lp_msg *msg)
{
	// This flag is only set here to propagate it to children messages
	msg->raw_flags += MSG_FLAG_PUBSUB;

	// On sub nodes, just create thread-level copies
	t_entry_arr threads = subscribersTable[msg->dest];

	// One message per thread
	int n_ch_count = array_count(threads);
	assert(n_ch_count);

	// For each subscribed Thread
	for(int c = 0; c < n_ch_count; c++){
		table_thread_entry_t *t_entry = &array_get_at(threads, c);

		// Child message is a clone of the message
		// Target holds the target thread's tid
		struct lp_msg *child_msg = msg_allocator_alloc(0);
		memcpy(child_msg, msg, msg_anti_size());
		child_msg->dest = t_entry->tid;
		// Flags are already set upon reception

		// Push child message into target thread's incoming queue
		pubsub_msg_queue_insert(child_msg);
	}

	// Remove the flag to make it be freed as a normal msg
	msg->raw_flags -= MSG_FLAG_PUBSUB;

	/* This msg is only a blueprint: Free it right now */
	msg_allocator_free(msg);
}

// Ok?
/// This carries out the antimessaging when the node is the one responsible for the publisher LP
void pub_node_handle_published_antimessage(struct lp_msg *msg)
{
	struct n_pubsub_info *pi = n_pub_info_ptr(msg);

	// Carry out the antimessaging
	// the message originated from the local node. More precisely, it came from
	// current_LP as this can only be called when rolling back
	struct lp_msg **children = pi->m_arr;
	int child_count = pi->m_cnt;
	int it = 0;

#ifdef ROOTSIM_MPI
	// Antimessage via MPI
	int subbed_nodes = current_lp->n_remote_sub_nodes;
	if (subbed_nodes) {
		struct block_bitmap* subs = current_lp->subnodes;
		nid_t dest_nid = 0;
		do {
			while (!bitmap_check(subs, dest_nid))
				++dest_nid;

			mpi_remote_anti_msg_send(children[it], dest_nid);
			++it;
		} while (it < subbed_nodes);
	}
#endif

	// Now antimessage every local child by setting flag + requeueing
	for (; it < child_count; ++it) {
		struct lp_msg *cmsg = children[it];

		// Set anti flag
		int cflags = atomic_fetch_add_explicit(
			&cmsg->flags, MSG_FLAG_ANTI,
			memory_order_relaxed);

		// Requeue if already processed
		if (cflags & MSG_FLAG_PROCESSED) {
			// cmsg->dest contains the tid of target thread
			pubsub_msg_queue_insert(cmsg);
		}
	}

	msg_allocator_free(msg);
}

/// Gets the positive message matching a_msg by removing it from past_pubsubs.
static inline struct lp_msg* get_positive_pubsub_msg(const struct lp_msg *a_msg)
{
	uint32_t m_id = a_msg->raw_flags >> MSG_FLAGS_BITS, m_seq = a_msg->m_seq;
	array_count_t past_i = 0;
	array_count_t ct = array_count(past_pubsubs);
	while (past_i<ct) {
		struct lp_msg *msg = unmark_msg(array_get_at(past_pubsubs, past_i));
		if ((msg->raw_flags >> MSG_FLAGS_BITS) == m_id && msg->m_seq == m_seq) {
			atomic_fetch_add_explicit(&msg->flags, MSG_FLAG_ANTI, memory_order_relaxed);

			// TODO: is it better to free while antimessaging
			//  or letting the fossil collector free it later?
			// Prevents double free
			return msg;
		}
		past_i++;
	}

	return NULL;
}

// This is called when a pubsub message with ANTI flag set is extracted
void thread_handle_published_antimessage(struct lp_msg *anti_msg){

	// Is the publisher local to node
	if (!(anti_msg->flags >> MSG_FLAGS_BITS)) {
		if (anti_msg->flags & MSG_FLAG_PROCESSED){
			thread_actually_antimessage(anti_msg);
		} else {
			// Was still not processed
			pubsub_thread_msg_free(anti_msg);
		}
		return;
	}

	// Did we process the positive copy of this message already?
	struct lp_msg *msg = get_positive_pubsub_msg(anti_msg);
	if (unlikely(!msg)) {
		// No, store for later
		anti_msg->raw_flags >>= MSG_FLAGS_BITS;
		array_push(early_pubsub_antis, anti_msg);
		return;
	}

	// We did
	msg_allocator_free(anti_msg);

	msg->raw_flags += MSG_FLAG_ANTI;
	thread_actually_antimessage(msg);
}

/// Carries out antimessaging of thread-level pubsub message
static void thread_actually_antimessage(struct lp_msg *msg)
{
	struct t_pubsub_info *pi = t_pub_info_ptr(msg);

	for (size_t i = 0; i < pi->m_cnt; i++) {
		struct lp_msg *cmsg = pi->m_arr[i];

		uint32_t cflags = atomic_fetch_add_explicit(
				&cmsg->flags, MSG_FLAG_ANTI, memory_order_relaxed);

		if (cflags & MSG_FLAG_PROCESSED)
			msg_queue_insert(cmsg);
	}
	// FIXME: aggiorna il resto per essere coerente con questo:
	// Tolgo flag processed. Se da past_pubsubs un messaggio non ha il flag processed lo posso freeare
	// Se ha processed e anti, allora sta ancora in msg_queue
	// La cosa si presenta solo in cleanup.
	msg->raw_flags -= MSG_FLAG_PROCESSED;
}

// OK
table_thread_entry_t new_thread_entry(void)
{
	// Create a new array to keep LPs of curr thread
	lp_entry_arr lparr;
	array_init(lparr);
	// Create thread subscription entry
	table_thread_entry_t entry = {
		.tid = rid,
		.lp_arr = lparr
	};
	return entry;
}

// OK
// return the index of the thread or that of its successor (i.e. where it should be)
int seekInSortedArray(t_entry_arr arr, rid_t tid){
	int size = array_count(arr);
	if(!size) return -1;
	//~ int min = array_get_at(arr, 0);
	//~ int max = array_get_at(arr, size-1);
	int min = 0;
	int max = size-1;

	int next;
	rid_t curr;
	do {
		next = (max+min)/2;
		curr = array_get_at(arr, next).tid;
		if (tid < curr){
			max = next-1;
		} else if (tid > curr) {
			min = next+1;
		} else { // Found it
			break;
		}
	} while( max >= min );
	if (max<min) return min; // Did not find it!
	return next;
}

extern void add_subbed_node(lp_id_t sub_id, lp_id_t pub_id);
// OK
// Crashes if LP pub_id is not node-local. Data race if LP pub_id is not owned by local thread
inline void add_subbed_node(lp_id_t sub_id, lp_id_t pub_id){
	nid_t sub_node_id = lid_to_nid(sub_id);
	assert(sub_id != pub_id);
	if (!bitmap_check(lps[pub_id].subnodes, sub_node_id)) {
		bitmap_set(lps[pub_id].subnodes, sub_node_id);

		// Newly subbed node is remote
		if (sub_node_id != nid) {
			lps[pub_id].n_remote_sub_nodes++;
			// Can also do this with bitmap_count_set after initialization
		}
	}
}

// OK
// Subscribe LP subscriber_id to LP publisher_id
void SubscribeAndHandle(lp_id_t dirty_subscriber_id, lp_id_t publisher_id, void* data){
	lp_id_t subscriber_id = dirty_subscriber_id & ~LP_ID_MSB; // Clear leftmost bit

	bool sub_is_local = lid_to_nid(subscriber_id)==nid && lid_to_rid(subscriber_id)==rid;
	bool pub_is_local = lid_to_nid(publisher_id)==nid && lid_to_rid(publisher_id)==rid;

	if(pub_is_local && !sub_is_local){ // If the publisher is managed by local thread
		add_subbed_node(subscriber_id, publisher_id);
		return;
	}

	// Acquire the lock
	spin_lock(&(tableLocks[publisher_id]));

	// dyn_array holding entries for threads subbed to pub_LP
	t_entry_arr *subbedThreads_p = &(subscribersTable[publisher_id]);

	int size = array_count(*subbedThreads_p);
	int pos = 0;

	if (!size){ // Uninitialized subscribersTable entry
		array_init(*subbedThreads_p);

		// Push a new thread entry in the array in table[pub_id]
		array_push(*subbedThreads_p, new_thread_entry());

		// pos is 0. Correct.

	} else { // Array has at least one thread entry
		// Is there an entry for the current thread?
		pos = seekInSortedArray(*subbedThreads_p, rid);

		if (pos == size){ // No entry - add at end of array
			// Add a new thread entry at the end of table[pub_id]
			array_push(*subbedThreads_p, new_thread_entry());
			// pos is size. Correct

		} else if (array_get_at(*subbedThreads_p, pos).tid!=rid){ // No entry - add in middle of array
			// Add a new thread entry at index pos of table[pub_id]
			array_add_at(*subbedThreads_p, pos, new_thread_entry());
			// pos is the correct index
		}
	}

	// Thread's entry is at index pos

	// Insert lp_entry in the thread's array
	table_lp_entry_t lp_entry = {
		.lid = dirty_subscriber_id,
		.data = data
	};

	array_push(array_get_at(*subbedThreads_p, pos).lp_arr, lp_entry);

	// Release the lock
	spin_unlock(&(tableLocks[publisher_id]));
}

//OK
void Subscribe(lp_id_t subscriber_id, lp_id_t publisher_id){
	// Set most significant bit of subscriber_id to indicate it is a normal Subscribe
	SubscribeAndHandle(
		subscriber_id | LP_ID_MSB,
		publisher_id,
		NULL);
}

// OK
// Free a (node level) pubsub msg
// If the message is node-level, then it will have local children and possibly
//  MPI children and will properly behave
// If the message is thread-level, then this was called from outside pubsub.c,
// meaning that it was called by msg_allocator_free
// message_allocator_free may only be called on a pubsub msg in one of two cases:
// 1. When fossil collecting => the message is node-level pubsub
// 2. When cleaning up (msg_queue_fini) => the message is thread level.
// In case 2, the thread-level message will not have any children
// the check on c_ptr will thus fail and behave exactly as pubsub_thread_msg_free
void pubsub_msg_free(struct lp_msg* msg)
{
	struct n_pubsub_info *pi = n_pub_info_ptr(msg);

	// If it has children
	size_t children_ct = pi->m_cnt;
	// We can use count of local subscribers to know
	// how many children were destined to MPI
	if (likely(children_ct)) {
		array_count_t local_children = array_count(subscribersTable[msg->dest]);
		children_ct -= local_children;
		// Free children for MPI
		for (size_t i = 0; i < children_ct; i++){
			msg_allocator_free(pi->m_arr[i]);
		}
	}
	// Free the array pointing to children
	// The local children will be freed independently
	mm_free(pi->m_arr);
}

void pubsub_thread_msg_free(struct lp_msg* msg)
{
	// Works for both node and thread-level
	// No race conditions: either GVT>dest_t, or already antimsgd
	msg->raw_flags &= ~MSG_FLAG_PUBSUB;

	if(msg->pl_size)
		mm_free(t_pub_info_ptr(msg)->m_arr);

	msg_allocator_free(msg);
}

// OK
// Insert a thread-level pubsub message in queue
void pubsub_msg_queue_insert(struct lp_msg* msg){

	// msg->dest contains the thread id
	unsigned dest_rid = msg->dest;
	struct msg_queue *mq = &queues[dest_rid];
	struct q_elem qe = {.t = msg->dest_t, .m = msg};

	spin_lock(&mq->q_lock);
	array_push(mq->b, qe);
	spin_unlock(&mq->q_lock);

}

void pubsub_on_gvt(simtime_t current_gvt)
{
	// Garbage collect past_pubsubs
	array_count_t ct = array_count(past_pubsubs);
	array_count_t i;

	for(i = 0; i < ct; ++i){
		struct lp_msg *msg = array_get_at(past_pubsubs, i);
		if (msg->dest_t >= current_gvt)
			break;

		pubsub_thread_msg_free(msg);
	}
	array_truncate_first(past_pubsubs, i);

}

#if LOG_LEVEL <= LOG_DEBUG
void pprint_table_lp_entry_t(table_lp_entry_t e, FILE *f){
	fprintf(f, "%lu", e.lid);
}

//typedef dyn_array(table_lp_entry_t) lp_entry_arr;
void pprint_lp_entry_arr(lp_entry_arr a, FILE *f){
	// è un dyn_array di table_lp_entry_t
	size_t ct = array_count(a);
	fprintf(f, "[");
	for(size_t i=0; i<ct; i++){
		pprint_table_lp_entry_t(array_get_at(a, i), f);
		if(i+1<ct) fprintf(f, ", ");
	}
	fprintf(f, "]\n");
}

void pprint_table_thread_entry_t(table_thread_entry_t e, FILE *f){
	fprintf(f, "\t\t");
	fprintf(f, "{\n");
	fprintf(f, "\t\t\t\"Tid\": %u,\n\t\t\t\"LP_arr\": ", e.tid);
	pprint_lp_entry_arr(e.lp_arr, f);
	fprintf(f, "\t\t");
	fprintf(f, "}");
}

void pprint_t_entry_arr(t_entry_arr a, FILE *f){
	fprintf(f, "\t");
	fprintf(f, "[\n");
	// è un dyn_array di table_thread_entry_t
	size_t ct = array_count(a);
	for(size_t i=0; i<ct; i++){
		pprint_table_thread_entry_t(array_get_at(a, i), f);
		if(i+1<ct) fprintf(f, ",");
		fprintf(f, "\n");
	}
	fprintf(f, "\t");
	fprintf(f, "]");
}

// Print the table to a file
void pprint_subscribers_table(FILE *f){
	t_entry_arr *t = subscribersTable;
	fprintf(f, "[\n");
	// Size of t: n_lps;
	for(lp_id_t i=0; i<n_lps; i++){
		//Entry dei sub dell'LP con id "i": t[i];
		pprint_t_entry_arr(t[i], f);
		if(i+1<n_lps) fprintf(f, ",");
		fprintf(f, "\n");
	}
	fprintf(f, "]\n");
}

void pprint_subscribers_adjacency_list(FILE *f){
	// Subscribers table
	t_entry_arr *t = subscribersTable;
	fprintf(f, "[\n");
	// Size of t: n_lps;
	for(lp_id_t pub_id=0; pub_id<n_lps; pub_id++){
		// Working with subscribers table entries, which are dyn_arrays of t_entry
		t_entry_arr t_entry_arr = t[pub_id];

		fprintf(f, "\t[");

		size_t t_ct = array_count(t_entry_arr);
		for(size_t t_i=0; t_i<t_ct; t_i++){
			//Now we get the lp_arr in t_entry, which are dyn_arrays of lp_entry
			lp_entry_arr lp_arr = array_get_at(t_entry_arr, t_i).lp_arr;

			size_t lp_ct = array_count(lp_arr);
			for(size_t l_i=0; l_i<lp_ct; l_i++){
				fprintf(f, "%lu", array_get_at(lp_arr, l_i).lid);
				if(t_i+1<t_ct || l_i+1<lp_ct) fprintf(f, ", ");
			}

		}

		fprintf(f, "]");
		if(pub_id+1<n_lps) fprintf(f, ",");
		fprintf(f, "\n");
	}
	fprintf(f, "]\n");
}
#endif

#endif
