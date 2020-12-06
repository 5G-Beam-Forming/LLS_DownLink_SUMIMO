/*
 *	LDPC_CC.h
 *
 *		C++ interface to LDPC C++ module.
 *		Programmed by PJW
 *		Last Updated: 2018-12-11
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */
 
#ifndef __LDPC_H__
#define __LDPC_H__

#include "sumimo.h"
#include "LDPC/method_ChannelCoding.h"


void LDPC_C_initialization(NR_ChannelCoding &, int, int, int, double, ChCoding &);
void LDPC_C_tx(NR_ChannelCoding &, ChCoding &);
void LDPC_C_rx(NR_ChannelCoding &, ChCoding &);



#endif	/* __LDPC_H__ : end of LDPC_CC.h */