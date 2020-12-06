/*
 * module_BS_SUSISO.h
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *		Last Updated: 2018-10-01
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */
#ifndef __MODULE_BS_SUMIMO_H__
#define __MODULE_BS_SUMIMO_H__

#include <armadillo>
#include "sumimo.h"
#include "module_Parameter_MIMO.h"
#include "LDPC_CC.h"


using namespace arma;



struct HARQBuffer{
	uword					num_process;
	uword					max_retrans;
	field<rowvec>			coded_bits;
	uvec					process;
	uvec					timing;
	field<uvec>				retrans;
};   /* end of HARQ_buffer structure definition */



class ModuleBS_SUMIMO {
public:
	/* Related to data stream */
	field<urowvec>		bit_stream;
	field<urowvec>		coded_bit_stream;
	field<urowvec>		scrambled_bit_stream;
	field<cx_rowvec>	modulated_signal_stream;
	cx_mat				OFDM_signal_stream;
	cx_cube				resource_grid;

	/* RS */
	cx_cube				resource_grid_DMRS;
	cx_mat				resource_grid_DMRS_subframe;
	cx_cube				resource_grid_CSIRS;
	ucube				DMRS_position;
	umat				DMRS_position_total;
	cx_mat				DMRS_signal;
	cx_mat				CSIRS_signal;
	
	ucube				CSIRS_position;
	umat				CSIRS_position_total;
	umat				data_position_total;

	int					mu;

	/* Related to signaling */
	Genie				genie;

	/* Related to resource allocation */
	umat				PCFICH_position;
	umat				Data_position;
	imat				PSS_signal;
	imat				SSS_signal;
	urowvec				SS_PBCH_block_start;

	umat				PSS_position;
	umat				SSS_position;
	umat				PSS_5G_position;
	umat				SSS_5G_position;
	umat				PBCH_DMRS_5G_position;
	umat				PBCH_5G_position;
	umat				PDCCH_5G_position;
	umat				Control_position;
	umat				RS_position;
	umat				num_Data_RB;

	/* Related to MCS */
	irowvec				M_order;
	rowvec				Coding_rate;
	irowvec				CQI;
	int					CQI_table_num;
	field<ChCoding>		ch_coding;
	int					N_RE;

	/* Related to Power allocation */
	rowvec				power_data;
	rowvec				power_DMRS;
	rowvec				power_CSIRS;
	
	/* Related to HARQ */
	HARQBuffer			HARQ_buffer;
	uvec				HARQ_process_index;
	uvec				HARQ_newdata_indicator;

	/* otherwise */
	int					cell_ID;
	int					n_RNTI = 1;
	int					PMI;
	int					RI;
	int					num_layers;
	uword				ind_WB;

	int					num_codewords;

	field<cx_mat>		Narrowbeam;
	field<cx_mat>		Widebeam;
	field<cx_mat>		Precoder;
	cx_mat				Precoder_PDCCH;
	cx_mat				Precoder_CSIRS;
	mat					layered_PDCCH_stream;
	mat					layered_PBCH_stream;
	cx_mat				layered_signal;
	
	cx_mat				precoded_signal_stream;
	cx_mat				precoded_PDCCH_stream;
	cx_mat				precoded_PBCH_stream;


	/* class method definitions */
	ModuleBS_SUMIMO();
	~ModuleBS_SUMIMO();
	void module_BS_SUMIMO(ModuleParameterMIMO &);
	void Process(ModuleParameterMIMO &, NR_ChannelCoding &, int &, feedback &, ChannelOutput &);
	void method_feedback_reception(feedback &);
	void method_CSI_RS_allocation(ModuleParameterMIMO &, int &);
	void method_resource_allocation(ModuleParameterMIMO &, int &);
	void method_beam_management(ModuleParameterMIMO &, ChannelOutput &, int &);
	void method_DMRS_allocation(ModuleParameterMIMO &, int &);
	void method_data_allocation(ModuleParameterMIMO &, int &);
	void method_sync_generation(ModuleParameterMIMO &, int &);
	void method_HARQ_process(ModuleParameterMIMO &, int &, feedback &);
	void method_bit_generation(ModuleParameterMIMO &, int &);
	void method_channel_coding(ModuleParameterMIMO &, NR_ChannelCoding &, int &, feedback &);
	void method_scrambling(ModuleParameterMIMO &, int &);
	void method_modulation_mapping(void);
	void method_layer_mapping(ModuleParameterMIMO &);
	void method_precoding(ModuleParameterMIMO &);
	void method_OFDM_generation(ModuleParameterMIMO &);
	void method_genie(void);

};	/* end of class ModuleBS_SUMIMO */


/* local function definitions */
void func_pseudo_sequence_generation(int, double, mat &);



#endif		 /* end of module_BS_SUMIMO.h */