/**
* @file gvt/gvt.h
*
* @brief Global Virtual Time
*
* This module implements the GVT reduction. The current implementation
* is non blocking for observable simulation plaftorms.
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

#include <core/core.h>
#include <lp/msg.h>

extern void gvt_global_init(void);
extern simtime_t gvt_phase_run(void);
extern void gvt_on_msg_process(simtime_t msg_t);

#ifdef ROOTSIM_MPI

extern __thread bool gvt_phase_green;
extern __thread unsigned remote_msg_sent[MAX_NODES];
extern atomic_int remote_msg_received[2];

extern void gvt_on_start_ctrl_msg(void);
extern void gvt_on_done_ctrl_msg(void);

#define gvt_on_remote_msg_send(dest_nid)				\
__extension__({ remote_msg_sent[dest_nid]++; })

#define gvt_on_remote_msg_receive(msg_phase)				\
__extension__({ atomic_fetch_add_explicit(remote_msg_received + 	\
	msg_phase, 1U, memory_order_relaxed); })

#define gvt_phase_get() __extension__({ gvt_phase_green;})
#endif
