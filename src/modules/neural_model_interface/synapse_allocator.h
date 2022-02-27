#pragma once
#include <core/core.h>

/* <INTERFACES> */
void SetMaxSynSize(size_t size);
__syn_t* synapse_alloc(size_t size);

void synapse_allocator_init();
void synapse_allocator_fini();
/* </INTERFACES> */

struct __syn_alloc_block{
	size_t curr_offset; // Current offset inside data block
	size_t size;
	unsigned char data[];
};
