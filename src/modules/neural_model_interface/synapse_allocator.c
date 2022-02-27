#include <modules/neural_model_interface/synapse_allocator.h>
#include <modules/neural_model_interface/neural_model_interface.h>
#include <core/core.h>
#include <datatypes/array.h>

/**
 * An utterly plain and possibly stupid synapse allocator.
 * Simply allocates bigger memory blocks at once.
 * Cannot deallocate (why would you).
 */

#define SYN_ALLOC_BLOCK_EXP 18
#define BASE_MAX_SYN_SIZE offsetof(__syn_t, data) + sizeof(double)

// dyn_array delle zone allocate
static __thread dyn_array(struct __syn_alloc_block *) allocd_list = {0};
__thread struct __syn_alloc_block *curr_block = NULL;
__thread size_t max_sin_size = BASE_MAX_SYN_SIZE;

void SetMaxSynSize(size_t size){
	max_sin_size = size + offsetof(__syn_t, data);
}

void move_to_new_block(){
	// Allocate the new block
	curr_block = mm_alloc(offsetof(struct __syn_alloc_block, data) + (max_sin_size << SYN_ALLOC_BLOCK_EXP));
	curr_block->size = max_sin_size << SYN_ALLOC_BLOCK_EXP;
	curr_block->curr_offset = 0;
	// Add it to array of allocd
	array_push(allocd_list, curr_block);
}

void synapse_allocator_init(){
	array_init(allocd_list);
	// Allocate the first block
	move_to_new_block();
}

void synapse_allocator_fini(){
	for (size_t i=0; i<array_count(allocd_list); ++i){
		mm_free(array_get_at(allocd_list, i));
	}
	array_fini(allocd_list);
}

__syn_t* synapse_alloc(size_t size){
	// Not enough space!
	if (curr_block->curr_offset + size >= curr_block->size) {
		move_to_new_block();
	}
	size_t offset = curr_block->curr_offset;
	curr_block->curr_offset += size;
	return (__syn_t*) &curr_block->data[offset];
}
