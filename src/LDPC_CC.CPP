/*
 * LDPC_CC.cpp - Wrapper functions for LDPC C/C++ modules.
 *
 *		Wrapper functions interface to LDPC C/C++ module which is originally programmed by C/C++.
 *		Programmed by PJW
 *		Last Updated: 2018-12-11
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */

#include <cstdio>
#include <cstdlib>
#include <armadillo>
#include "sumimo.h"
#include "LDPC/stdafx.h"
#include "LDPC/SISO_PDSCH.h"
#include "LDPC/method_ChannelCoding.h"


using namespace arma;



/*
 * LDPC_C_initilization() - A wrapper function for Initialization LDPC method (Initialization method is renamed from NR_ChannelCoding() LDPC method.)
 */
void LDPC_C_initialization(NR_ChannelCoding &LDPC, int A, int G, int Qm, double Code_rate, ChCoding &ch_coding)
{
	/* local variable declarations */
	int			Mod_order;
	int			log2_Mod_order;
	NR_CODEC    nr_codec;
	int			B, C, L;
	int			CRC_Type, CRC_seg_Type;
	int			BG_index;
	int			K_cb, K_p, L_seg, B_p;
	int			CBGTI, CBGTI_present, C_p, N_L;

	Mod_order			 = Qm;
	nr_codec.data_length = A;
	nr_codec.code_rate   = Code_rate;
	nr_codec.G			 = G;
	nr_codec.ACK		 = 0;
	nr_codec.NACK        = 0;
	nr_codec.ldpc.LDPCMaxIter		 = ch_coding.Max_LDPC_itr;
	nr_codec.ldpc.LDPCDecOption		 = ch_coding.Decoding_type;
	nr_codec.ldpc.LDPC_MinSum_factor = ch_coding.MS_factor;


	log2_Mod_order = (int)(log((double)Mod_order) / log(2.0));

	/*-------- A and G  -------*/
	//int n_PRB = atoi(argv[2]); // total number of allocated PRBs for the UE
	int n_PRB     = 250;
	int N_symb_sh = 12;			// 12 or 14, 12: with pilots or 14: without pilots
								// limited buffer�� �̿��ϴ� ��� 13 OFDM symbols�� scheduling ��

	int *TBS_N_RE;
	TBS_N_RE = new int[2];
	// for LBRM
	int *TBS_LRBM;			
	TBS_LRBM = new int[2];	

	TBS_N_RE = TBS_determination(n_PRB, N_symb_sh, log2_Mod_order, &nr_codec);//TBS ;

	nr_codec.ldpc.LBRM_indicator = 1;							// 1:limited buffer rate matching, 0: othersize
	if (nr_codec.ldpc.LBRM_indicator == 1) {												
		TBS_LRBM = TBS_determination(n_PRB, N_symb_sh, log2_Mod_order, &nr_codec);	
		nr_codec.ldpc.TBS_LBRM = TBS_LRBM[0];											
	}																					

	nr_codec.harq.HARQ_indicator	 = 1;							// 1: HARQ on, 0: HAQR off
	nr_codec.harq.HARQ_type			 = 2;							// 1:CC-HARQ, 2:IR-HARQ

	if (nr_codec.harq.HARQ_indicator == 1)
		nr_codec.harq.HARQ_max_iter = 2;							// 1~16 ;
	else if (nr_codec.harq.HARQ_indicator == 0)
		nr_codec.harq.HARQ_max_iter = 1;

	nr_codec.harq.HARQ_iter = 0;
	nr_codec.ACK_cnt = 0;
	nr_codec.NACK_cnt = 0;

	nr_codec.harq.HARQ_index = new int[nr_codec.harq.HARQ_max_iter];
	int HARQ_index_temp[16] = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };	// can be manually modified

	for (int k = 0; k < nr_codec.harq.HARQ_max_iter; k++) {
		if (nr_codec.harq.HARQ_type == 1)
			nr_codec.harq.HARQ_index[k] = 0;
		else if (nr_codec.harq.HARQ_type == 2)
			nr_codec.harq.HARQ_index[k] = HARQ_index_temp[k];
	}

	/* Initialize NR_ChannelCoding class */
	LDPC.Initialization(Mod_order, &nr_codec);

	/* output parameters settings of LDPC_initialization() wrapper function */
	ch_coding.nr_codec.data_length = nr_codec.data_length;
	ch_coding.nr_codec.code_rate   = nr_codec.code_rate;
	ch_coding.nr_codec.G		   = nr_codec.G;
	ch_coding.nr_codec.ACK		   = nr_codec.ACK;
	ch_coding.nr_codec.NACK		   = nr_codec.NACK;
	ch_coding.nr_codec.ldpc.LDPCMaxIter		   = nr_codec.ldpc.LDPCMaxIter;
	ch_coding.nr_codec.ldpc.LDPCDecOption	   = nr_codec.ldpc.LDPCDecOption;
	ch_coding.nr_codec.ldpc.LDPC_MinSum_factor = nr_codec.ldpc.LDPC_MinSum_factor;

	/* CRC */
	if (A > 3824) {
		L = 24;			// CRC sequence length
		CRC_Type = 0;	// CRC24A
	}
	else {
		L = 16;
		CRC_Type = 2;	// CRC16
	}

	/* Base graph selection */
	if (A <= 292 || (A <= 3824 && Code_rate <= 0.67) || Code_rate <= 0.25)
		BG_index = 2;
	else
		BG_index = 1;

	/* Code block segmentation */
	B = A + L;

	if (BG_index == 1)
		K_cb = 8448;
	else
		K_cb = 3840;

	if (B <= K_cb) {
		L_seg = 0;
		C = 1;			// Total number of coded block
		B_p = B;
	}
	else {
		L_seg = 24;
		C = (int) ceil(B / (K_cb - L_seg));
		B_p = B + C*L_seg;

		CRC_seg_Type = 1;		// CRC24B
	}

	K_p = B_p / C;

	/* Rate matching */
	CBGTI = 1;			// 0 or 1, TBD
	CBGTI_present = 0;	// 0 or 1, TBD
	if (CBGTI_present == 0)
		C_p = C;
	N_L = 1;

	irowvec E = zeros<irowvec>(C);
	for (int r = 0, jj = 0; r < C; r++) {
		if (CBGTI == 0)
			E(r) = 0;
		else {
			if (jj <= C_p - mod(G / N_L * Qm, C_p) - 1)
				E(r) = N_L * Qm * (int)floor(G / (N_L * Qm * C_p));
			else
				E(r) = N_L * Qm * (int)ceil(G / (N_L * Qm * C_p));
			jj++;
		}
	}

	/* CRC */
	ch_coding.A = A;
	ch_coding.L = L;
	ch_coding.CRC_Type = CRC_Type;

	/* Base graph selection */
	ch_coding.BG_index = BG_index;

	/* Code block segmentation */
	ch_coding.B = B;
	ch_coding.B_p = B_p;
	ch_coding.L_seg = L_seg;
	ch_coding.C = C;
	ch_coding.K_p = K_p;
	
	if (C > 1)
		ch_coding.CRC_seg_Type = CRC_seg_Type;
	

	/* rate matching */
	ch_coding.E = E;
	ch_coding.Q = Qm;

	/* memory deallocations */
	E.reset();

}   /* end of LDPC_initialization() */





/*
 * LDPC_tx() - A LDPC_rx wrapper function for Tx_ChannelEncoding() fucntion.
 */
void LDPC_C_tx(NR_ChannelCoding &LDPC, ChCoding &ch_coding) {

	int			*data;
	NR_CODEC	 nr_codec;
	int			*tx_codeword;

	/* parameter initializations from ch_coding */
	nr_codec.data_length = ch_coding.nr_codec.data_length;
	nr_codec.code_rate   = ch_coding.nr_codec.code_rate;
	nr_codec.G			 = ch_coding.nr_codec.G;
	nr_codec.ACK		 = ch_coding.nr_codec.ACK;
	nr_codec.NACK		 = ch_coding.nr_codec.NACK;
	nr_codec.ldpc.LDPCMaxIter		 = ch_coding.nr_codec.ldpc.LDPCMaxIter;
	nr_codec.ldpc.LDPCDecOption		 = ch_coding.nr_codec.ldpc.LDPCDecOption;
	nr_codec.ldpc.LDPC_MinSum_factor = ch_coding.nr_codec.ldpc.LDPC_MinSum_factor;
	nr_codec.ldpc.LBRM_indicator	 = ch_coding.nr_codec.ldpc.LBRM_indicator;

	/* parameter initialization */
	int  Mod_order = (int)pow(2.0, 2.0);
	int  log2_Mod_order = (int)(log((double)Mod_order) / log(2.0));
	int *TBS_LRBM = new int[2];
	int  n_PRB = 1;			// total number of allocated PRBs for the UE
	TBS_LRBM = TBS_determination(n_PRB, 12, log2_Mod_order, &nr_codec);
	nr_codec.ldpc.TBS_LBRM = TBS_LRBM[0];

	/* the following setting values drived from MAIN_NR_PDSCH.cpp */
	nr_codec.harq.HARQ_indicator = 1;
	nr_codec.harq.HARQ_type		 = 2;
	nr_codec.harq.HARQ_iter		 = 0;
	if (nr_codec.harq.HARQ_indicator == 1)
		nr_codec.harq.HARQ_max_iter = 2;
	else /* = 0 */
		nr_codec.harq.HARQ_max_iter = 1;

	nr_codec.harq.HARQ_index = new int[nr_codec.harq.HARQ_max_iter];
	int HARQ_index_temp[16] = { 0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3 };
	for (int k = 0; k < nr_codec.harq.HARQ_max_iter; k++)
		nr_codec.harq.HARQ_index[k] = HARQ_index_temp[k];

	data = new int[nr_codec.data_length];
	for (uword i = 0; i < (uword)nr_codec.data_length; i++)
		data[i] = (int) ch_coding.tx_a(i);

	tx_codeword = Tx_ChannelEncoding(data, LDPC, &nr_codec);

	ch_coding.tx_g = zeros<urowvec>(nr_codec.G);
	for (uword i = 0; i < (uword)nr_codec.G; i++)
		ch_coding.tx_g(i) = tx_codeword[i];

	/* output parameters settings of LDPC_initialization() wrapper function */
	ch_coding.nr_codec.data_length = nr_codec.data_length;
	ch_coding.nr_codec.code_rate   = nr_codec.code_rate;
	ch_coding.nr_codec.G		   = nr_codec.G;
	ch_coding.nr_codec.ACK		   = nr_codec.ACK;
	ch_coding.nr_codec.NACK		   = nr_codec.NACK;
	ch_coding.nr_codec.ldpc.LDPCMaxIter		   = nr_codec.ldpc.LDPCMaxIter;
	ch_coding.nr_codec.ldpc.LDPCDecOption	   = nr_codec.ldpc.LDPCDecOption;
	ch_coding.nr_codec.ldpc.LDPC_MinSum_factor = nr_codec.ldpc.LDPC_MinSum_factor;

	ch_coding.nr_codec.ldpc.LBRM_indicator = nr_codec.ldpc.LBRM_indicator;
	ch_coding.nr_codec.ldpc.TBS_LBRM	   = nr_codec.ldpc.TBS_LBRM;

	ch_coding.nr_codec.harq.HARQ_indicator = nr_codec.harq.HARQ_indicator;
	ch_coding.nr_codec.harq.HARQ_type	   = nr_codec.harq.HARQ_type;
	ch_coding.nr_codec.harq.HARQ_iter	   = nr_codec.harq.HARQ_iter;
	ch_coding.nr_codec.harq.HARQ_max_iter  = nr_codec.harq.HARQ_max_iter;
	ch_coding.nr_codec.harq.HARQ_index	   = nr_codec.harq.HARQ_index;

	delete []data;
	delete []tx_codeword;

}	/* end of LDPC_tx() function */





/*
 * LDPC_rx() - LDPC_rx wrapper function
 */
void LDPC_C_rx(NR_ChannelCoding &LDPC, ChCoding &ch_coding) {

	double		*rx_LLR;
	NR_CODEC	 nr_codec;
	int			*rx_decoded_data;

	/* parameter initilizations */
	nr_codec.data_length = ch_coding.nr_codec.data_length;
	nr_codec.code_rate   = ch_coding.nr_codec.code_rate;
	nr_codec.G			 = ch_coding.nr_codec.G;
	nr_codec.ACK		 = ch_coding.nr_codec.ACK;
	nr_codec.NACK		 = ch_coding.nr_codec.NACK;
	nr_codec.ldpc.LDPCMaxIter		 = ch_coding.nr_codec.ldpc.LDPCMaxIter;
	nr_codec.ldpc.LDPCDecOption		 = ch_coding.nr_codec.ldpc.LDPCDecOption;
	nr_codec.ldpc.LDPC_MinSum_factor = ch_coding.nr_codec.ldpc.LDPC_MinSum_factor;
	nr_codec.harq.HARQ_iter			 = ch_coding.nr_codec.harq.HARQ_iter;
	nr_codec.harq.HARQ_index		 = ch_coding.nr_codec.harq.HARQ_index;

	rx_LLR = new double[ch_coding.symbol_LLR.n_elem];
	for (uword i = 0; i < ch_coding.symbol_LLR.n_elem; i++)
		rx_LLR[i] = ch_coding.symbol_LLR(i);

	rx_decoded_data = Rx_ChannelDecoding(rx_LLR, LDPC, &nr_codec);

	ch_coding.rx_a = zeros<urowvec>(nr_codec.data_length);
	for (uword i = 0; i < (uword)nr_codec.data_length; i++)
		ch_coding.rx_a(i) = (unsigned int)rx_decoded_data[i];

	/* output parameters settings of LDPC_initialization() wrapper function */
	ch_coding.nr_codec.data_length = nr_codec.data_length;
	ch_coding.nr_codec.code_rate   = nr_codec.code_rate;
	ch_coding.nr_codec.G		   = nr_codec.G;
	ch_coding.nr_codec.ACK		   = nr_codec.ACK;
	ch_coding.nr_codec.NACK		   = nr_codec.NACK;
	ch_coding.nr_codec.ldpc.LDPCMaxIter		   = nr_codec.ldpc.LDPCMaxIter;
	ch_coding.nr_codec.ldpc.LDPCDecOption	   = nr_codec.ldpc.LDPCDecOption;
	ch_coding.nr_codec.ldpc.LDPC_MinSum_factor = nr_codec.ldpc.LDPC_MinSum_factor;

	ch_coding.ACK_for_check = nr_codec.ACK;
	
	uword smCntA = 0, diffCntA = 0;
	for (uword i = 0; i < ch_coding.rx_a.n_elem; i++)
		if (ch_coding.rx_a(i) == ch_coding.tx_a(i))
			smCntA++;
		else
			diffCntA++;

	delete []rx_LLR;
	delete []rx_decoded_data;

}	/* end of LDPC_rx() function */



/* end of LDPC_initialization.cpp */