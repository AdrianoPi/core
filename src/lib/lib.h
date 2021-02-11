/**
* @file lib/lib.h
*
* @brief Model library main header
*
* This is the main header to initialize core model libraries.
*
* @copyright
* Copyright (C) 2008-2021 HPDCS Group
* https://hpdcs.github.io
*
* This file is part of ROOT-Sim (ROme OpTimistic Simulator).
*
* ROOT-Sim is free software; you can redistribute it and/or modify it under the
* terms of the GNU General Public License as published by the Free Software
* Foundation; only version 3 of the License applies.
*
* ROOT-Sim is distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
* A PARTICULAR PURPOSE. See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with
* ROOT-Sim; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/
#pragma once

#include <lib/random/random.h>
#include <lib/state/state.h>
#include <lib/topology/topology.h>

struct lib_ctx {
	// random library
	uint64_t rng_s[4];
	double unif;
	// todo remove
	void *state_s;
};

extern void lib_global_init(void);
extern void lib_global_fini(void);

extern void lib_lp_init(void);
extern void lib_lp_fini(void);

extern void lib_lp_init_pr(void);
extern void lib_lp_fini_pr(void);

