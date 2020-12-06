/*
 * module_UE_MIMO.h
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *		Last Updated: 2019-01-10
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */
#ifndef __MODULE_UE_MIMO_H__
#define __MODULE_UE_MIMO_H__

#include <armadillo>
#include "sumimo.h"
#include "module_Parameter_MIMO.h"
#include "LDPC_CC.h"


using namespace arma;


struct Result {
	irowvec			ACK;
	field<rowvec>	received_SNR;
	double			received_SNR_for_cali;
	irowvec			suc_data;
	irowvec			suc_bit;
	irowvec			fail_bit;
	irowvec			first_fail_bit;
	irowvec			suc_bit_data;
	irowvec			Error_count;
	irowvec			Totalbit;
};   /* end of Result struct definition */



/* HARQ Buffer specific to UE SISO module */
struct HARQBufferUE {
	int				num_process;
	irowvec			retransmission;
	unsigned int	max_retransmission;
	field<rowvec>	LLR;
	unsigned int	num_TB;
	field<ChCoding>	BS_signaling;
};   /* end of HARQ_buffer struct definition */



struct ChOutput {
	int				time;
	colvec			fft;
};   /* end of ch_Output struct definition */



struct EstFreqOffset {
	mat			frac;
	mat			res;
	mat			intt;
};



class ModuleUE_SUMIMO {
public:
	/* Related to data stream */
	cx_mat				received_signal_stream;
	cx_mat				received_signal_stream_sync;
	cx_cube				received_resource_grid;
	umat				Data_stream_indexf;
	cx_mat				Data_stream;

	/* Related to channel estimation */
	field<cx_cube>		H_estimated;
	EstFreqOffset		est_freq_offset;

	/* Related to received SNR */
	double				sigma_n_time;
	double				sigma_n_freq;

	/* Related to channel */
	cx_cube				HV_perfect;
	cx_cube				HV_estimated;
	cx_cube				HV_estimated_data;
	cx_cube				HV_estimated_power;
	field<cx_cube>		channel_buffer;

	/* Related to detection */
	mat					LLR;
	field<mat>			descrambled_LLR;
	field<mat>			descrambled_LLR_HARQ;
	field<mat>			delayered_LLR;


	/* Related to result calculation */
	field<umat>			coded_bit;
	field<umat>			decoded_bit;
	Result				Result;

	/* Related to HARQ & feedback */
	HARQBufferUE		HARQ_buffer;
	feedback			feedback;

	/* otherwise */
	uword				cell_ID;


	/* class method definitions */
	ModuleUE_SUMIMO();
	~ModuleUE_SUMIMO();
	void module_UE_SUMIMO(ModuleParameterMIMO &);
	void module_UE_SUSISO(ModuleParameterMIMO &);
	void process(ModuleParameterMIMO &, NR_ChannelCoding &, Genie &, ChannelOutput &, int);
	void method_channel_filtering(ModuleParameterMIMO &, Genie &, ChannelOutput &);
	void method_synchronization(ModuleParameterMIMO &, Genie &, int);
	void method_OFDM_demodulation(ModuleParameterMIMO &);
	void method_channel_estimation(ModuleParameterMIMO &, Genie &, ChannelOutput &, int);
	void method_Data_extraction(ModuleParameterMIMO &, Genie &);
	void method_detection(ModuleParameterMIMO &, Genie &);
	void method_delayer_mapping(Genie &);
	void method_descrambling(ModuleParameterMIMO &, Genie &, int);
	void method_HARQ_process(ModuleParameterMIMO &, Genie &);
	void method_decoding(ModuleParameterMIMO &, NR_ChannelCoding &, Genie &, int);
	void method_Result(ModuleParameterMIMO &, Genie &, int);

};   /* end of ModuleUE_SISO class definition */



#endif	/* end of module_UE_MIMO.h */