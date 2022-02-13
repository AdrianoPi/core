/**
 * @file test/core/sync.c
 *
 * @brief Test: synchronization primitives test
 *
 * SPDX-FileCopyrightText: 2008-2021 HPDCS Group <rootsim@googlegroups.com>
 * SPDX-License-Identifier: GPL-3.0-only
 */
#include <stdatomic.h>
#include <test.h>

#include <core/core.h>
#include <core/sync.h>

#define N_THREADS 16
#define THREAD_REP 100000

static atomic_uint counter;
static spinlock_t spin;
static mrswlock_t rw_lock;
static unsigned k = 0;
static unsigned d = 0;

_Static_assert(THREAD_REP % 5 == 0, "THREAD_REP must be a multiple of 5 for the spinlock test to work.");

static test_ret_t sync_barrier_test(__unused void *_)
{
	int ret = 0;
	unsigned j = THREAD_REP;
	while(j--) {
		atomic_fetch_add_explicit(&counter, 1U, memory_order_relaxed);
		sync_thread_barrier();
		unsigned val = atomic_load_explicit(&counter, memory_order_relaxed);
		sync_thread_barrier();
		ret += (int) (val % N_THREADS);
	}

	return ret;
}

static test_ret_t spinlock_test(__unused void *_)
{
	unsigned long long ret = 0;
	unsigned j = THREAD_REP;

	while(j--) {
		spin_lock(&spin);
		switch(k % 5) {
			case 0:
					__attribute__((fallthrough));
			case 2:
				k += 1;
					__attribute__((fallthrough));
			case 4:
				k += 1;
				break;
			case 3:
					__attribute__((fallthrough));
			case 1:
				k += 5;
				break;
			default:
				__builtin_unreachable();
		}
		spin_unlock(&spin);
	}

	sync_thread_barrier();
	ret += k != ((THREAD_REP * 5 * N_THREADS + 2) / 3);

	return ret;
}

static test_ret_t mrswlock_test(__unused void *_)
{
	unsigned long long ret = 0;
	unsigned j = THREAD_REP;

	while(j--) {
		mrswlock_wlock(&rw_lock, (int)N_THREADS);
		switch(d % 5) {
			case 0: __attribute__((fallthrough));
			case 2:
				d += 1;
				__attribute__((fallthrough));
			case 4:
				d += 1;
				break;
			case 3: __attribute__((fallthrough));
			case 1:
				d += 5;
				break;
			default:
				__builtin_unreachable();
		}
		mrswlock_wunlock(&rw_lock, N_THREADS);
	}

	sync_thread_barrier();
	ret += d != ((THREAD_REP * 5 * N_THREADS + 2) / 3);
	sync_thread_barrier();
	d = 0;
	k = 0;
	sync_thread_barrier();

	j = THREAD_REP;
	if(rid) {
		while(j--) {
			mrswlock_rlock(&rw_lock);
			ret += d != k;
			mrswlock_runlock(&rw_lock);
		}
	} else {
		while(j--) {
			mrswlock_wlock(&rw_lock, N_THREADS);
			d += 3;
			k += 2;
			k += d % 5;
			d += 2;
			mrswlock_wunlock(&rw_lock, N_THREADS);
		}
	}

	sync_thread_barrier();
	return ret;
}

void foo() {}

struct simulation_configuration conf = {
    .lps = N_THREADS,
    .n_threads = N_THREADS,
    .dispatcher = (ProcessEvent_t)foo,
    .committed = (CanEnd_t)foo,
};


int main(void)
{
	init(N_THREADS);

	RootsimInit(&conf);
	parallel_test("Testing synchronization barrier", sync_barrier_test, NULL);
	parallel_test("Testing spinlock", spinlock_test, NULL);
	mrswlock_init(&rw_lock, N_THREADS);
	parallel_test("Testing mrswlock", mrswlock_test, NULL);

	finish();
}
