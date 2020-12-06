/*
 * module_UE_MIMO.cpp
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *		Last Updated: 2019-01-10
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */

#include <armadillo>
#include "sumimo.h"
#include "module_UE_MIMO.h"
#include "LDPC_CC.h"
#include "getCQITable.h"

#undef		max;
#undef		min;


using namespace arma;




/*
 * module_UE_MIMO() - custom constructor
 */
void ModuleUE_SUMIMO::module_UE_SUMIMO(ModuleParameterMIMO &para)
{
	
	this->module_UE_SUSISO(para);	/* obj@module_UE_SISO(para) */

	int num_slot_total = para.num_slot_in_subframe * para.num_subframe;

	this->result.ACK		  = zeros<irowvec>(num_slot_total);
	this->result.received_SNR = field<rowvec>(num_slot_total);
	this->result.suc_data	  = zeros<irowvec>(num_slot_total);
	this->result.Error_count  = zeros<irowvec>(num_slot_total);
	this->result.Totalbit	  = zeros<irowvec>(num_slot_total);
	this->result.suc_bit	  = zeros<irowvec>(num_slot_total);

	this->feedBack.CQI.reset();
	this->feedBack.PMI = 1;
	this->feedBack.RI  = 1;

	this->HARQ_buffer.num_process = para.num_HARQ_process;
	this->HARQ_buffer.retransmission = zeros<irowvec>(this->HARQ_buffer.num_process);
	this->HARQ_buffer.max_retransmission = para.num_HARQ_retransmission;
	this->HARQ_buffer.LLR = field<rowvec>(this->HARQ_buffer.num_process);
	this->HARQ_buffer.num_TB = 0;
	this->HARQ_buffer.BS_signaling = field<ChCoding>(this->HARQ_buffer.num_process);

}	/* end of module_UE_MIMO() method definition */





/*
 * Module_UE_MIMO::module_UE_SISO()
 */
void ModuleUE_SUMIMO::module_UE_SUSISO(ModuleParameterMIMO &para) {

	int num_slot_total = para.num_slot_in_subframe * para.num_subframe;

	this->result.ACK = zeros<irowvec>(num_slot_total);
	this->result.received_SNR = field<rowvec>(num_slot_total);
	this->result.suc_data = zeros<irowvec>(num_slot_total);
	this->result.Error_count = zeros<irowvec>(num_slot_total);
	this->result.Totalbit = zeros<irowvec>(num_slot_total);
	this->result.suc_bit = zeros<irowvec>(num_slot_total);
	this->result.suc_bit_data = zeros<irowvec>(num_slot_total);

	this->feedBack.CQI.clear();
	this->feedBack.ACK = zeros<irowvec>(1);
	this->feedBack.ACK(0) = 1;

	this->channel_buffer = field<cx_cube>(para.num_RB * para.num_sc, para.num_symb);
	for (uword i = 0; i < para.num_RB * para.num_sc; i++)
		for (uword j = 0; j < para.num_symb; j++)
			this->channel_buffer(i, j) = zeros<cx_cube>(para.num_Rx_antenna, para.num_port_CSIRS, num_slot_total);

	this->HARQ_buffer.num_process = para.num_HARQ_process;
	this->HARQ_buffer.retransmission = zeros<irowvec>(this->HARQ_buffer.num_process);
	this->HARQ_buffer.max_retransmission = para.num_HARQ_retransmission;
	this->HARQ_buffer.LLR = field<rowvec>(this->HARQ_buffer.num_process);
	this->HARQ_buffer.num_TB = 0;
	this->HARQ_buffer.BS_signaling = field<ChCoding>(this->HARQ_buffer.num_process);

}	/* end of module_UE_SUSISO() method definition */





/*
 * process() - Gateway function of UE module.
 */
void ModuleUE_SUMIMO::process(ModuleParameterMIMO &para, NR_ChannelCoding &LDPC, Genie &genie, ChannelOutput &ch_Output, int ind_slot) {

	method_channel_filtering(para, genie, ch_Output);

	method_synchronization(para, genie, ind_slot);

	method_OFDM_demodulation(para);

	method_channel_estimation(para, genie, ch_Output, ind_slot);

	if (genie.selected) {
		
		method_Data_extraction(para, genie);

		method_detection(para, genie);

		method_delayer_mapping(genie);

		method_descrambling(para, genie, ind_slot);

		method_decoding(para, LDPC, genie, ind_slot);

		method_Result(para, genie, ind_slot);

	}

}	/* end of process() method definition */





/*
 * method_channel_filtering and noise addition() - Channel filtering and noise addition
 */
void ModuleUE_SUMIMO::method_channel_filtering(ModuleParameterMIMO &para, Genie &genie, ChannelOutput &ch_Output)
{
	/* local variable declaration */
	cx_vec		temp_conv;
	cx_vec		sq_ch_time;
	rowvec		temp_vec;
	cx_mat		noise_received;

	switch (para.channel.time_fading) {
	case Ch_mode::Block_Fading:
		this->received_signal_stream = zeros<cx_mat>(genie.OFDM_signal_stream.n_rows, para.num_Rx_antenna);
		for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++) {
			temp_conv = zeros<cx_vec>(genie.OFDM_signal_stream.n_rows + ch_Output.time.n_slices - 1);
			for (uword ind_tx = 0; ind_tx < para.num_Tx_antenna; ind_tx++) {
				sq_ch_time = ch_Output.time(span(ind_rx), span(ind_tx), span::all);
				temp_conv += conv(genie.OFDM_signal_stream.col(ind_tx), sq_ch_time);
			}
			this->received_signal_stream.col(ind_rx) = temp_conv(span(0, genie.OFDM_signal_stream.n_rows - 1));
		}
	}

	/* add noise */
	this->received_signal_stream = this->received_signal_stream(span(0, genie.OFDM_signal_stream.n_rows - 1), span::all);
	temp_vec = zeros<rowvec>(this->received_signal_stream.n_rows);
	for (uword i = 0; i < this->received_signal_stream.n_rows; i++)
		temp_vec(i) = (double)i;
	
	this->received_signal_stream = this->received_signal_stream % repmat(exp(cx_double(0, 1) * 2.0 * datum::pi * ((double)para.sync_freq_offset * temp_vec).st() / (double)para.size_fft), 1, para.num_Rx_antenna);
	this->sigma_n_freq = 1.0 / pow(10.0, para.SNR / 10.0);
	this->sigma_n_time = (double)para.size_fft / (double)(para.num_RB * para.num_sc) * 1.0 / pow(10.0, para.SNR / 10.0);
	noise_received = sqrt(this->sigma_n_time) * (randn(arma::size(this->received_signal_stream)) + cx_double(0, 1) * randn(arma::size(this->received_signal_stream))) / sqrt(2.0);
	this->received_signal_stream += noise_received;

	/* local memory deallocation */
	temp_conv.reset();
	sq_ch_time.reset();
	temp_vec.reset();
	noise_received.reset();

}	/* end of method_channel filtering() method definition */





/*
 * method_synchronization() - Synchronization
 */
void ModuleUE_SUMIMO::method_synchronization(ModuleParameterMIMO &para, Genie &genie, int ind_slot)
{
	/* local varaible definition */
	rowvec		temp_vec;
	uvec		ind_ofdm_symbol(4);
	field<umat>	index_received_temp(4);
	urowvec		index_temp_2, index_temp_3, index_temp_4;
	urowvec		temp_vec_1, temp_vec_2, temp_vec_3, temp_vec_4;
	urowvec		Index_CP;
	cx_mat		received_CP1, received_CP2;
	cx_rowvec	received_CRS_sub, CRS_sub, RFO_epsilon;
	uword		index_argmax;

	switch (para.sync_freq_mode) {
	case Sync_freq_mode::Perfect:
		temp_vec = zeros<rowvec>(this->received_signal_stream.n_rows);
		for (uword i = 0; i < this->received_signal_stream.n_rows; i++)
			temp_vec(i) = (double)i;
		this->received_signal_stream_sync = this->received_signal_stream % repmat(exp(cx_double(0, -1) * 2.0 * datum::pi * ((double)para.sync_freq_offset * temp_vec).t() / (double)para.size_fft), 1, para.num_Rx_antenna);
		
		temp_vec.reset();	// local memeory deallocation
		break;
	case Sync_freq_mode::estimated:
		/* Estimate the fractional frequencey offset (FFO) */
		ind_ofdm_symbol(0) = para.size_fft + para.length_CP(0);
		ind_ofdm_symbol(1) = para.size_fft + para.length_CP(0) + (para.size_fft + para.length_CP(1)) * (para.num_symb / 2 - 1);
		ind_ofdm_symbol(2) = 2 * (para.size_fft + para.length_CP(0)) + (para.size_fft + para.length_CP(1)) * (para.num_symb / 2 - 1);
		ind_ofdm_symbol(3) = 2 * (para.size_fft + para.length_CP(0)) + 2 * (para.size_fft + para.length_CP(1)) + (para.num_symb / 2 - 1);

		index_received_temp(0) = find(ones<uvec>(ind_ofdm_symbol(0)));
		
		index_temp_2 = zeros<urowvec>(ind_ofdm_symbol(1) - (ind_ofdm_symbol(0)) + 1);
		for (uword i = 0, j = ind_ofdm_symbol(0) + 1; i < ind_ofdm_symbol(1) - (ind_ofdm_symbol(0) + 1); i++, j++)
			index_temp_2(i) = j;
		index_received_temp(1) = zeros<urowvec>(index_temp_2.n_rows * para.size_fft + para.length_CP(1), index_temp_2.n_cols * (para.num_symb / 2 - 1));
		index_received_temp(1) = reshape(index_temp_2, para.size_fft + para.length_CP(1), para.num_symb / 2 - 1);
		
		index_received_temp(2) = zeros<urowvec>(ind_ofdm_symbol(2) - (ind_ofdm_symbol(1)) + 1);
		for (uword i = 0, j = ind_ofdm_symbol(1) + 1; i < ind_ofdm_symbol(1) - (ind_ofdm_symbol(0) + 1); i++, j++)
			index_received_temp(2)(i) = j;

		index_temp_4 = zeros<urowvec>(ind_ofdm_symbol(3) - (ind_ofdm_symbol(2) + 1));
		for (uword i = 0, j = ind_ofdm_symbol(2) + 1; i < ind_ofdm_symbol(3) - (ind_ofdm_symbol(2) + 1); i++, j++)
			index_temp_4(i) = j;
		index_received_temp(3) = zeros<urowvec>(index_temp_4.n_rows * para.size_fft + para.length_CP(1), index_temp_4.n_cols * (para.num_symb / 2 - 1));
		index_received_temp(3) = reshape(index_temp_4, para.size_fft + para.length_CP(1), para.num_symb / 2 - 1);

		temp_vec_1 = index_received_temp(0)(0, para.length_CP(0) - 1);
		temp_vec_2 = reshape(index_received_temp(1).rows(0, para.length_CP(1)), 1, index_received_temp(1).n_elem);
		temp_vec_3 = index_received_temp(2)(0, para.length_CP(0));
		temp_vec_4 = reshape(index_received_temp(3).rows(1, para.length_CP(1)), 1, index_received_temp(3).n_elem);
		Index_CP = join_horiz(join_horiz(temp_vec_1, temp_vec_2), join_horiz(temp_vec_3, temp_vec_4));

		received_CP1 = this->received_signal_stream(Index_CP);
		received_CP2 = this->received_signal_stream(Index_CP + para.size_fft);

		this->est_freq_offset.frac = -1.0 / (2.0 * datum::pi) * arg(sum(received_CP1 % conj(received_CP2)));

		temp_vec = zeros<rowvec>(this->received_signal_stream.n_rows);
		for (uword i = 0; i < this->received_signal_stream.n_rows; i++)
			temp_vec(i) = (double)i;
		this->received_signal_stream_sync = this->received_signal_stream % exp(cx_double(0, -1) * 2.0 * datum::pi * this->est_freq_offset.frac * temp_vec).t() / para.size_fft;

		// Estimate the integer frequency offset (IFO)
		// We should obtain genie.received_resource_grid after FFO
		this->method_OFDM_demodulation(para);

		cx_mat IFO_epsilon;
		if (mod(ind_slot, 10) == 0 || mod(ind_slot, 10) == 5) {
			for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++)
				for (uword ind_tx = 0; ind_tx < para.num_Tx_antenna; ind_tx++) {
					IFO_epsilon = zeros<cx_mat>(63, 1);
					for (uword ind_IFO = -31; ind_IFO < 31; ind_IFO++) {
						IFO_epsilon(32 + ind_IFO) = real(exp(cx_double(0, 1) * 2.0 * datum::pi * (double)ind_IFO * (double)para.length_CP(1) / (double)para.size_fft)
							* sum(conj(this->received_resource_grid(shift(genie.SSS_position, -1 * ind_IFO))) % this->received_resource_grid(shift(genie.PSS_position, -1 * ind_IFO))
							% conj(genie.PSS_position_signal) % genie.DMRS_signal));
					}
					index_argmax = IFO_epsilon.index_max();
				}
			this->est_freq_offset.intt = 32 - index_argmax;
			temp_vec = zeros<rowvec>(this->received_signal_stream.n_rows);
			for (uword i = 0; i < this->received_signal_stream.n_rows; i++)
				temp_vec(i) = (double)i;
			this->received_signal_stream_sync = this->received_signal_stream_sync % exp(cx_double(0, -1) * 2.0 * datum::pi * this->est_freq_offset.intt * temp_vec.t() / (double)para.size_fft); 
		}
		else {
			temp_vec = zeros<rowvec>(this->received_signal_stream.n_rows);
			for (uword i = 0; i < this->received_signal_stream.n_rows; i++)
				temp_vec(i) = (double)i;
			this->received_signal_stream_sync = this->received_signal_stream_sync % exp(cx_double(0, -1) * 2.0 * datum::pi * this->est_freq_offset.intt * temp_vec.t() / (double)para.size_fft);
		}

		/* Estimate the residual frequency offset (RFO) */
		this->method_OFDM_demodulation(para);

		for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++)
			for (uword ind_tx = 0; ind_tx < para.num_Tx_antenna; ind_tx++) {
				received_CRS_sub = this->received_resource_grid(find(genie.CRS_position.slice(ind_tx)));
				CRS_sub = genie.resource_grid(find(genie.CRS_position.slice(ind_tx)));
				RFO_epsilon = sum(received_CRS_sub(span(0, 1)) % conj(received_CRS_sub(span(2, 3)) % CRS_sub(span(2, 3)) % conj(CRS_sub(span(0, 1)))));
			}

		this->est_freq_offset.res = -arg(RFO_epsilon) / (2.0 * datum::pi) * (double)para.size_fft / (this->received_signal_stream_sync.n_elem / 2);
		this->received_signal_stream_sync = this->received_signal_stream_sync % exp(cx_double(0, -1) * 2.0 * datum::pi * this->est_freq_offset.res * temp_vec).t() / para.size_fft;

		/* local memory deallocation */
		ind_ofdm_symbol.reset();
		for (uword i = 0; i < index_received_temp.n_elem; i++)
			index_received_temp(i).reset();
		index_received_temp.reset();
		index_temp_2.reset();
		index_temp_3.reset();
		index_temp_4.reset();
		temp_vec_1.reset();
		temp_vec_2.reset();
		temp_vec_3.reset();
		temp_vec_4.reset();
		Index_CP.reset();
		received_CP1.reset();
		received_CP2.reset();
		received_CRS_sub.reset();
		CRS_sub.reset();
		RFO_epsilon.reset();
		break;
	}

}	/* end of method_synchronization() method deifinition */





/*
 * method_OFDM_demodulation() - OFDM demodulation
 */
void ModuleUE_SUMIMO::method_OFDM_demodulation(ModuleParameterMIMO &para)
{
	/* local variable definition */
	rowvec			ind_ofdm_symbol(3);
	field<cx_mat>	received_signal_stream_temp(4, 1);
	cx_mat			resource_grid_temp(para.size_fft, para.num_symb);
	cx_mat			resource_grid_fft(para.size_fft, para.num_symb);
	uword			cidx, ridx, i;

	this->received_resource_grid = zeros<cx_cube>(para.num_RB * para.num_sc, para.num_symb, para.num_Rx_antenna);

	ind_ofdm_symbol(0) = para.size_fft + para.length_CP[0];
	ind_ofdm_symbol(1) = para.size_fft + para.length_CP[0] + (para.size_fft + para.length_CP[1]) * ((double)para.num_symb / 2.0 - 1.0);
	ind_ofdm_symbol(2) = ind_ofdm_symbol(1) + para.size_fft + para.length_CP[0];

	for (uword ind_Rx = 0; ind_Rx < para.num_Rx_antenna; ind_Rx++) {
		received_signal_stream_temp(0) = this->received_signal_stream_sync(span(0, ind_ofdm_symbol(0) - 1), span(ind_Rx));
		received_signal_stream_temp(1) = reshape(this->received_signal_stream_sync(span(ind_ofdm_symbol(0), ind_ofdm_symbol(1)), span(ind_Rx)), para.size_fft + para.length_CP[1], para.num_symb / 2 - 1);
		received_signal_stream_temp(2) = this->received_signal_stream_sync(span(ind_ofdm_symbol(1), ind_ofdm_symbol(2) - 1), span(ind_Rx));
		received_signal_stream_temp(3) = reshape(this->received_signal_stream_sync(span(ind_ofdm_symbol(2), this->received_signal_stream_sync.n_rows - 1), span(ind_Rx)), para.size_fft + para.length_CP[1], para.num_symb / 2 - 1);

		/* Remove CP */
		for (i = 0, cidx = 0; i < received_signal_stream_temp(0).n_cols; i++, cidx++)
			resource_grid_temp.col(cidx) = received_signal_stream_temp(0)(span(para.length_CP[0], received_signal_stream_temp(0).n_rows - 1), span(i));
		for (i = 0; i < received_signal_stream_temp(1).n_cols; i++, cidx++)
			resource_grid_temp.col(cidx) = received_signal_stream_temp(1)(span(para.length_CP[1], received_signal_stream_temp(1).n_rows - 1), span(i));
		for (i = 0; i < received_signal_stream_temp(2).n_cols; i++, cidx++)
			resource_grid_temp.col(cidx) = received_signal_stream_temp(2)(span(para.length_CP[0], received_signal_stream_temp(2).n_rows - 1), span(i));
		for (i = 0; i < received_signal_stream_temp(3).n_cols; i++, cidx++)
			resource_grid_temp.col(cidx) = received_signal_stream_temp(3)(span(para.length_CP[1], received_signal_stream_temp(3).n_rows - 1), span(i));

		/* FFT : can be a bottle neck */
		clock_t UE_fft1start = clock();
		resource_grid_fft = sqrt((double)para.num_RB * para.num_sc) / (double)para.size_fft * fft(resource_grid_temp);
		clock_t UE_fft1end = clock();

		/* Remove DC carrier, zero padding */
		for (ridx = 0, i = resource_grid_fft.n_rows - para.num_RB * para.num_sc / 2; i < resource_grid_fft.n_rows; i++, ridx++)
			this->received_resource_grid(span(ridx), span::all, span(ind_Rx)) = resource_grid_fft(span(i), span::all);
		for (i = 0; i < para.num_RB * para.num_sc / 2; i++, ridx++)
			this->received_resource_grid(span(ridx), span::all, span(ind_Rx)) = resource_grid_fft(span(i), span::all);
	}

	/* local memory deallocation */
	ind_ofdm_symbol.reset();
	for (uword i = 0; i < received_signal_stream_temp.n_elem; i++)
		received_signal_stream_temp(i).reset();
	received_signal_stream_temp.reset();
	resource_grid_temp.reset();
	resource_grid_fft.reset();

}	/* end of method_OFDM_demodulation() method definition */





/*
 * method_channel_estimation() - Channel estimation
 */
void ModuleUE_SUMIMO::method_channel_estimation(ModuleParameterMIMO &para, Genie &genie, ChannelOutput &ch_Output, int ind_slot)
{
	/* local variable decalaration */
	field<cx_cube>	H_hat, H_hat_temp;
	bool			interpolation_DMRS;
	cx_mat			sq_H_hat_temp, temp_H_hat_temp;
	cx_cube			HV_hat, HV_temp;
	uvec			this_RB_freq;
	cx_cube			received_grid_temp;
	cx_cube			resource_grid_temp;
	cx_cube			DMRS_signal;
	uword			est_freq_part;
	uvec			posi_freq_time;
	umat			time_cov;
	rowvec			time_interval;
	uvec			this_RB_freq_part;
	uvec			posi_vec, posi_freq, posi_time, posi;


	H_hat = ch_Output.fft;
	this->H_estimated = ch_Output.fft;

	/* Channel Estimation */
	if (!genie.selected)
		return;			/*  If this user is not served at this time, just a return */

	interpolation_DMRS = true;
	H_hat_temp = ch_Output.fft;
	this->HV_perfect = zeros<cx_cube>(para.num_RB * para.num_sc, para.num_symb, para.num_Rx_antenna);

	/* Effective channel (Channel X Precoder) calculation for CSIR */
	sq_H_hat_temp = zeros<cx_mat>(para.num_sc, H_hat_temp(0).n_slices);
	for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++) {
		for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {

			for (uword i = 0, j = para.num_sc * ind_RB; j < para.num_sc * (ind_RB + 1); i++, j++)		// sqeeze H_hat_temp
				for (uword n = 0; n < H_hat_temp(0).n_slices; n++)
					sq_H_hat_temp(i, n) = H_hat_temp(j)(0, ind_rx, n);
			
			if (para.num_Tx_antenna == 1 && para.num_Rx_antenna == 1) {
				for (uword i = 0, j = para.num_sc * ind_RB; j < para.num_sc * (ind_RB + 1); i++, j++)
					for (uword n = 0; n < H_hat_temp.n_slices; n++)
						this->HV_perfect(j, n, ind_rx) = sq_H_hat_temp(j, n);
			}
			else {
				temp_H_hat_temp = sq_H_hat_temp * genie.Precoder(ind_RB);

				for (uword i = 0, j = para.num_sc * ind_RB; j < para.num_sc * (ind_RB + 1); i++, j++)
					for (uword n = 0; n < temp_H_hat_temp.n_cols; n++)
						this->HV_perfect(j, n, ind_rx) = temp_H_hat_temp(i, n);
			}

			for (uword i = para.num_sc * ind_RB; i < para.num_sc * (ind_RB + 1); i++)	// repmat
				for (uword j = 0; j < para.num_symb; j++)
					this->HV_perfect(i, j, ind_rx) = this->HV_perfect(i, 0, ind_rx);

		}
	}

	HV_hat = zeros<cx_cube>(para.num_RB * para.num_sc, para.num_symb, para.num_Rx_antenna);
	received_grid_temp = zeros<cx_cube>(arma::size(this->received_resource_grid));
	resource_grid_temp = zeros<cx_cube>(arma::size(genie.resource_grid_DMRS));
	DMRS_signal = zeros<cx_cube>(arma::size(resource_grid_temp));

	/* Channel estimation for detection (LS, MMSE) */
	switch (para.ch_est_mode_DMRS) {
	case Ch_est_mode_DMRS::Perfect:
		HV_hat = this->HV_perfect;
		interpolation_DMRS = false;
		break;
	case Ch_est_mode_DMRS::LS:
		/* For the phase2 calibration */
		HV_temp = zeros<cx_cube>(para.num_RB * para.num_sc, para.num_symb, para.num_Rx_antenna);

		if (para.cal_scenario) {	// MMSE estimation
			for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++)
				for (uword ind_port = 0; ind_port < genie.RI; ind_port++)
					for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
						this_RB_freq = ind_RB * para.num_sc + find(ones<uvec>(para.num_sc));
						for (uword i = 0; i < para.num_sc; i++) {
							HV_temp(span(this_RB_freq(i)), span::all, span(ind_rx)) = HV_hat(span(this_RB_freq(i)), span::all, span(ind_rx));
							received_grid_temp(span(this_RB_freq(i)), span::all, span(ind_rx)) = this->received_resource_grid(span(this_RB_freq(i)), span::all, span(ind_rx));
						}

						for (uword i = 0; i < this_RB_freq.n_elem; i++)
							for (uword j = 0; j < genie.DMRS_position.n_cols; j++)
								if (genie.DMRS_position(i, j, ind_port) != 0)
									HV_temp(i, j, ind_port) = received_grid_temp(i, j, ind_port);


						for (uword i = 0; i < this_RB_freq.n_elem; i++)
							for (uword j = 0; j < genie.resource_grid_DMRS.n_cols; j++)
								resource_grid_temp(i, j, ind_port) = genie.resource_grid_DMRS(this_RB_freq(i), j, ind_port);

						for (uword i = 0; i < this_RB_freq.n_elem; i++)
							for (uword j = 0; j < genie.DMRS_position.n_cols; j++)
								if (genie.DMRS_position(i, j, ind_port) != 0) {
									DMRS_signal(i, j, ind_port) = resource_grid_temp(i, j, ind_port);
									HV_temp(i, j, ind_port) = HV_temp(i, j, ind_port) / DMRS_signal(i, j, ind_port);
								}

						for (uword i = 0; i < this_RB_freq.n_elem; i++)
							HV_hat(span(this_RB_freq(i)), span::all, span(ind_port)) = HV_temp;
					}
		}
		else {		// DMRS periodicity (CDM per # of est_freq_part subcarriers)
			switch (para.DMRS_type) {
			case 1:
				est_freq_part = 4;
				break;
			case 2:
				est_freq_part = 6;
				break;
			}

			for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++)
				for (uword ind_port = 0; ind_port < genie.RI; ind_port++)
					for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++)
						for (uword ind_part = 0; ind_part < para.num_sc / est_freq_part; ind_part++) {
							this_RB_freq = ind_RB * para.num_sc + find(ones<uvec>(para.num_sc)) * est_freq_part + find(ones<uvec>(est_freq_part));
							for (uword i = 0; i < para.num_sc; i++) {
								HV_temp(span(this_RB_freq(i)), span::all, span(ind_rx)) = HV_hat(span(this_RB_freq(i)), span::all, span(ind_rx));
								received_grid_temp(span(this_RB_freq(i)), span::all, span(ind_rx)) = this->received_resource_grid(span(this_RB_freq(i)), span::all, span(ind_rx));
							}

							for (uword i = 0; i < this_RB_freq.n_elem; i++)
								for (uword j = 0; j < genie.DMRS_position.n_cols; j++)
									if (genie.DMRS_position(i, j, ind_port) != 0)
										HV_temp(i, j, ind_port) = received_grid_temp(i, j, ind_port);


							for (uword i = 0; i < this_RB_freq.n_elem; i++)
								for (uword j = 0; j < genie.resource_grid_DMRS.n_cols; j++)
									resource_grid_temp(i, j, ind_port) = genie.resource_grid_DMRS(this_RB_freq(i), j, ind_port);

							for (uword i = 0; i < this_RB_freq.n_elem; i++)
								for (uword j = 0; j < genie.DMRS_position.n_cols; j++)
									if (genie.DMRS_position(i, j, ind_port) != 0)
										DMRS_signal(i, j, ind_port) = resource_grid_temp(i, j, ind_port);
									
							posi_vec = find(HV_temp);
							posi_time = posi_vec / HV_temp.n_rows;
							posi_freq = posi_vec - posi_time * HV_temp.n_rows;
							//HV_temp(posi_freq, posi_time) = repmat(DMRS_signal.t() * nonezero(HV_temp), posi_freq.n_elem, posi_time.n_elem/DMRS_signal.size());

							for (uword i = 0; i < this_RB_freq.n_elem; i++)
								HV_hat(span(this_RB_freq(i)), span::all, span(ind_port)) = HV_temp;
						}
		}

		interpolation_DMRS = true;
		break;

	case Ch_est_mode_DMRS::MMSE:	/* 2D(freq-time) MMSE channel estimation */
		cx_mat	HV_temp;
		cx_mat	received_grid_temp;
		cx_mat	resource_grid_temp;
		cx_mat	DMRS_signal;
		umat	DMRS_position_temp;
		uvec	posi_vec, posi_freq, posi_time, posi;
		cx_mat	HV_temp_vec;
		cx_mat	freq_cov, total_cov;

		switch (para.channel.time_fading) {
		case Ch_mode::Block_Fading:
			time_cov = toeplitz(ones<urowvec>(para.num_symb));
			break;
		case Ch_mode::Fast_Fading:
			time_interval = zeros<rowvec>(para.num_symb);
			double temp = 0.0;
			for (uword ind = 0; ind < para.num_symb; ind++) {
				time_interval(ind) = temp;
				if (ind != 7)
					temp = temp + (double)para.length_CP(1) * (double)para.size_fft;
				else
					temp = temp + (double)para.length_CP(0) * (double)para.size_fft;
			}
			// time_cov = toeplitz(besselj(0, 2.0 * PI * para.user_speed * para.f / para.light_speed * time_interval * 1.0 / para.Fs));

			break;
		}	/* end of switch (para.channel.time_fading) */


		if (para.cal_scenario) {	/* MMSE estimation */
			for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++)
				for (uword ind_port = 0; ind_port < genie.RI; ind_port++)
					for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
						this_RB_freq = ind_RB * para.num_sc + find(ones<urowvec>(para.num_sc));
						for (uword i = 0; i < this_RB_freq_part.n_elem; i++)		// squeeze
							for (uword j = 0; j < para.num_symb; j++) {
								HV_temp(i, j) = HV_hat(this_RB_freq_part(i), j, ind_rx);
								received_grid_temp(i, j) = this->received_resource_grid(this_RB_freq_part(i), j, ind_rx);
								DMRS_position_temp(i, j) = genie.DMRS_position(this_RB_freq_part(i), j, ind_port);
							}

						for (uword i = 0; i < this_RB_freq_part.n_elem; i++)
							for (uword j = 0; j < genie.DMRS_position.n_cols; j++)
								if (genie.DMRS_position(this_RB_freq_part(i), j, ind_port) != 0)
									HV_temp(i, j) = received_grid_temp(i, j);

						for (uword i = 0; i < this_RB_freq_part.n_elem; i++)
							for (uword j = 0; j < genie.resource_grid_DMRS.n_cols; j++)
								resource_grid_temp(i, j) = genie.resource_grid_DMRS(this_RB_freq_part(i), j, ind_port);

						umat rowsz = max(sum(DMRS_position_temp, 1), 0);
						umat colsz = max(sum(DMRS_position_temp, 0), 1);
						DMRS_signal = zeros<cx_mat>(rowsz(0, 0), colsz(0, 0));
						DMRS_signal = resource_grid_temp(find(DMRS_position_temp));

						posi_vec = find(HV_temp);
						posi_time = posi_vec / HV_temp.n_rows;
						posi_freq = posi_vec - posi_time * HV_temp.n_rows;
						posi = posi_freq + para.num_sc * posi_time;
						HV_temp_vec = nonzeros(HV_temp);

						freq_cov = para.Cov_Matrix_freq(this_RB_freq, this_RB_freq);

						mat db_time_cov(arma::size(time_cov));
						for (uword i = 0; i < time_cov.n_elem; i++)
							db_time_cov(i) = (double)time_cov(i);	// type conversion

						total_cov = kron(db_time_cov, freq_cov);

						HV_temp(posi_freq, posi_time) = reshape(DMRS_signal.t() * nonzeros(HV_temp), posi_freq.n_elem, posi_time.n_elem) / DMRS_signal.n_elem;
						for (uword i = 0; i < this_RB_freq_part.n_elem; i++)
							for (uword j = 0; j < HV_hat.n_cols; j++)
								HV_hat(this_RB_freq_part(i), j, ind_rx) = HV_temp(i, j);
					}
		}
		else {		// if para.cal_scenario = false
			switch (para.DMRS_type) {
			case 1:
				est_freq_part = 4;
				break;
			case 2:
				est_freq_part = 6;
				break;
			}

			this_RB_freq_part = zeros<uvec>(est_freq_part);
			HV_temp = zeros<cx_mat>(est_freq_part, para.num_symb);
			received_grid_temp = zeros<cx_mat>(est_freq_part, para.num_symb);
			resource_grid_temp = zeros<cx_mat>(est_freq_part, para.num_symb);
			DMRS_position_temp = zeros<umat>(est_freq_part, para.num_symb);

			for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++)
				for (uword ind_port = 0; ind_port < genie.RI; ind_port++)
					for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
						HV_temp = zeros<cx_mat>(est_freq_part, para.num_symb);
						for (uword ind_part = 0; ind_part < para.num_sc / est_freq_part; ind_part++) {
							this_RB_freq_part = ind_RB * para.num_sc + ind_part * est_freq_part + find(ones<uvec>(est_freq_part));

							for (uword i = 0; i < this_RB_freq_part.n_elem; i++)		// squeeze
								for (uword j = 0; j < para.num_symb; j++) {
									HV_temp(i, j) = HV_hat(this_RB_freq_part(i), j, ind_rx);
									received_grid_temp(i, j) = this->received_resource_grid(this_RB_freq_part(i), j, ind_rx);
									DMRS_position_temp(i, j) = genie.DMRS_position(this_RB_freq_part(i), j, 0);
								}

							for (uword i = 0; i < this_RB_freq_part.n_elem; i++)
								for (uword j = 0; j < genie.DMRS_position.n_cols; j++)
									if (genie.DMRS_position(this_RB_freq_part(i), j, 0) != 0)
										HV_temp(i, j) = received_grid_temp(i, j);

							for (uword i = 0; i < this_RB_freq_part.n_elem; i++)
								for (uword j = 0; j < genie.resource_grid_DMRS.n_cols; j++)
									resource_grid_temp(i, j) = genie.resource_grid_DMRS(this_RB_freq_part(i), j, 0);

							umat rowsz = max(sum(DMRS_position_temp, 1), 0);
							umat colsz = max(sum(DMRS_position_temp, 0), 1);
							DMRS_signal = zeros<cx_mat>(rowsz(0, 0), colsz(0, 0));
							DMRS_signal = resource_grid_temp(find(DMRS_position_temp));

							posi_vec = find(HV_temp);
							posi_time = posi_vec / HV_temp.n_rows;
							posi_freq = posi_vec - posi_time * HV_temp.n_rows;

							HV_temp(posi_freq, posi_time) = repmat(DMRS_signal.t() * nonzeros(HV_temp), posi_freq.n_elem, posi_time.n_elem) / DMRS_signal.n_elem;

							for (uword i = 0; i < this_RB_freq_part.n_elem; i++)
								for (uword j = 0; j < HV_hat.n_cols; j++)
									HV_hat(this_RB_freq_part(i), j, ind_rx) = HV_temp(i, j);

							HV_temp.fill(0);
							received_grid_temp.fill(0);
							resource_grid_temp.fill(0);
							DMRS_position_temp.fill(0);
						}

						HV_temp = zeros<cx_mat>(para.num_sc, para.num_symb);

						this_RB_freq = ind_RB * para.num_sc + find(ones<uvec>(para.num_sc));

						for (uword i = 0; i < para.num_sc; i++)
							for (uword j = 0; j < para.num_symb; j++)
								HV_temp(i, j) = HV_hat(this_RB_freq(i), j, ind_rx);

						posi_vec = find(HV_temp);
						posi_time = posi_vec / HV_temp.n_rows;
						posi_freq = posi_vec - posi_time * HV_temp.n_rows;

						posi = posi_freq + para.num_sc * posi_time;
						HV_temp_vec = nonzeros(HV_temp);

						freq_cov = para.Cov_Matrix_freq(this_RB_freq, this_RB_freq);
						
						mat dbl_time_cov(arma::size(time_cov));
						for (uword i = 0; i < time_cov.n_elem; i++)
							dbl_time_cov(i) = (double)time_cov(i);	// type conversion

						total_cov = kron(dbl_time_cov, freq_cov);

						HV_temp = reshape(total_cov.cols(posi) * solve(total_cov(posi, posi) + para.sigma_n_freq / pow(abs(DMRS_signal(0)), 2) * eye<mat>(posi.n_elem, posi.n_elem), HV_temp_vec), para.num_sc, para.num_symb);
						
						for (uword i = 0; i < this_RB_freq.n_elem; i++)
							for (uword j = 0; j < HV_hat.n_cols; j++)
								HV_hat(this_RB_freq(i), j, ind_rx) = HV_temp(i, j);

					}

		}	/* end of if (para.cal_scenario) */
		
		interpolation_DMRS = false;

		break;
	}	/* end of switch (para.ch_est_mode_DMRS) */


	/* Channel interpolation after LS channel estimation */
	if (interpolation_DMRS) {
		cx_vec	HV_temp, H_temp, temp;
		uvec	inter_point, tempf, posi, posi_freq, posi_time;

		switch (para.channel.time_fading) {
		case Ch_mode::Block_Fading:
			switch (para.ch_interp_mode_DMRS) {
			case Ch_interp_mode_DMRS::linear:
				cx_vec	HV_temp, temp;
				uvec	inter_point;

				for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++)
					for (uword ind_port = 0; ind_port < genie.RI; ind_port)
						for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
							this_RB_freq = ind_RB * para.num_sc + find(ones<uvec>(para.num_sc));
							for (uword i = 0; i < this_RB_freq.n_elem; i++)
								HV_temp(i) = HV_hat(this_RB_freq(i), 3, ind_rx);
							temp = nonzeros(HV_temp);
							HV_temp(0) = temp(0);
							HV_temp(HV_temp.n_elem - 1) = temp(temp.n_elem - 1);
							inter_point = find(HV_temp);
							//interp1(inter_point, nonzeros(HV_temp), find(ones<uvec>(HV_temp.n_elem)), HV_temp, para.ch_interp_mode_DMRS);
							HV_hat(span(this_RB_freq(0), this_RB_freq(this_RB_freq.n_elem - 1)), span::all, span(ind_rx)) = repmat(HV_temp, 1, para.num_symb);
						}
			}
			break;
		case Ch_mode::Fast_Fading:
			switch (para.ch_interp_mode) {
			case Ch_interp_mode::linear:
				for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++)
					for (uword ind_port = 0; ind_port < para.num_port; ind_port)
						for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
							for (uword ind_symb = 0; ind_symb < para.num_symb; ind_symb++) {
								this_RB_freq = ind_RB * para.num_sc + find(ones<uvec>(para.num_sc));
								for (uword i = 0; i < this_RB_freq.n_elem; i++)
									HV_temp(i) = HV_hat(this_RB_freq(i), ind_symb, ind_rx);
								temp = nonzeros(HV_temp);
								HV_temp(0) = temp(0);
								HV_temp(HV_temp.n_elem - 1) = temp(temp.n_elem - 1);
								inter_point = find(HV_temp);
								//interp1(inter_point, nonzeros(HV_temp), find(ones<uvec>(HV_temp.n_elem)), HV_temp, para.ch_interp_mode_DMRS);
								HV_hat(span(this_RB_freq(0), this_RB_freq(this_RB_freq.n_elem - 1)), span::all, span(ind_rx)) = repmat(HV_temp, 1, para.num_symb);
							}
							tempf = find(HV_hat(span(0), span::all, span(ind_rx)));
							HV_hat(span::all, span(0), span(ind_rx)) = HV_hat(span::all, span(tempf(0)), span(ind_rx));
							HV_hat(span::all, span(HV_hat.n_cols - 1), span(ind_rx)) = HV_hat(span::all, span(tempf(tempf.n_elem - 1)), span(ind_rx));
							for (uword ind_sc = ind_RB * para.num_sc; ind_sc < ind_RB * para.num_sc; ind_sc++) {
								inter_point = find(HV_hat(span(ind_sc), span::all, span(ind_rx)));
								//HV_hat(span(ind_sc), span::all, span(ind_rx)) = interp1(inter_point, nonzeros(HV_hat(span(ind_sc), span::all, span(ind_rx))), max(HV_hat.n_cols, ind_sc) , para.ch_interp_mode);
							}
						}
				break;
			case Ch_interp_mode::TDI:
				cx_vec	H_temp = zeros<cx_vec>(H_hat.n_rows);
				cx_vec	g_temp;
				uword	threshold, last;

				for (uword ind_rx = 0; ind_rx < para.num_Rx_antenna; ind_rx++)
					for (uword ind_port = 0; ind_port < para.num_port; ind_port) {
						for (uword ind_symb = 0; ind_symb < para.num_symb; ind_symb++) {
							for (uword i = 0; i < H_hat.n_rows; i++)
								H_temp = H_hat(i, ind_symb, ind_rx);
							
							posi = find(H_temp);
							posi_time = posi / HV_temp.n_rows;
							posi_freq = posi - posi_time * HV_temp.n_rows;

							g_temp = ifft(nonzeros(H_temp));
							threshold = 4;
							H_temp(span(0, threshold)) = g_temp(span(0, threshold));
							last = H_temp.n_elem - 1;
							H_temp(span(last - threshold, last)) = g_temp(span(last - threshold, last));
							H_temp == shift(fft(H_temp), posi_freq(0));
							for (uword i = 0; i < H_hat.n_rows; i++)
								H_hat(i, ind_symb, ind_rx) = H_temp(i);
						}
						tempf = find(HV_hat(span(0), span::all, span(ind_rx)));
						HV_hat(span::all, span(0), span(ind_rx)) = HV_hat(span::all, span(tempf(0)), span(ind_rx));
						HV_hat(span::all, span(HV_hat.n_cols - 1), span(ind_rx)) = HV_hat(span::all, span(tempf(tempf.n_elem - 1)), span(ind_rx));
						for (uword ind_sc = ind_port * para.num_sc; ind_sc < ind_port * para.num_sc; ind_sc++) {
							inter_point = find(HV_hat(span(ind_sc), span::all, span(ind_rx)));
							//HV_hat(span(ind_sc), span::all, span(ind_rx)) = interp1(inter_point, nonzeros(HV_hat(span(ind_sc), span::all, span(ind_rx))), max(HV_hat.n_cols, ind_sc) , para.ch_interp_mode);
						}
					}
				break;
			}
			break;
		}	/* end of switch (para.channel.time_fading) */
	}

	this->HV_estimated = HV_hat;
	if (para.power_DMRS.is_empty())
		this->HV_estimated_power = HV_hat;
	else
		this->HV_estimated_power = HV_hat / pow(para.power_DMRS(0), 2.0);


	/* local memory deallocation */
	for (uword i = 0; i < H_hat.n_elem; i++)
		H_hat(i).reset();
	H_hat.reset();
	for (uword i = 0; i < H_hat_temp.n_elem; i++)
		H_hat_temp(i).reset();
	H_hat_temp.reset();
	sq_H_hat_temp.reset();
	temp_H_hat_temp.reset();
	HV_hat.reset();
	HV_temp.reset();
	this_RB_freq.reset();
	received_grid_temp.reset();
	resource_grid_temp.reset();
	DMRS_signal.reset();
	posi_freq_time.reset();
	time_cov.reset();
	time_interval.reset();
	this_RB_freq_part.reset();
	posi_vec.reset();
	posi_freq.reset();
	posi_time.reset();
	posi.reset();

}	/* end of method_channel_estimation method definition */





/*
 * method_Data_extraction() - Data extraction
 */
void ModuleUE_SUMIMO::method_Data_extraction(ModuleParameterMIMO &para, Genie &genie)
{
	/* local variable declaration */
	uword		data_length, Data_index, length_temp;
	ucube		Data_stream_indexf_temp;
	umat		Data_position_temp, Data_indexf_temp;
	cx_mat		Data_temp;
	uword		temp_position, temp_length;
	cx_cube		temp, HV_temp;
	umat		temp_position_matrix;


	data_length = sum(sum(genie.Data_position));
	this->Data_stream = zeros<cx_mat>(data_length, para.num_Rx_antenna);
	this->Data_stream_indexf = zeros<umat>(data_length, para.num_Rx_antenna);
	Data_stream_indexf_temp = zeros<ucube>(para.num_RB * para.num_sc, para.num_symb, para.num_Rx_antenna);
	for (uword i = 0; i < para.num_RB * para.num_sc; i++)
		for (uword j = 0; j < para.num_symb; j++)
			for (uword k = 0; k < para.num_Rx_antenna; k++)
				Data_stream_indexf_temp(i, j, k) = i;
	
	/* Data_stream */
	for (uword ind_port = 0; ind_port < para.num_Rx_antenna; ind_port++) {
		Data_index = 0;
		for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
			Data_position_temp = genie.Data_position(span(ind_RB * para.num_sc, (ind_RB + 1) * para.num_sc - 1), span::all);
			length_temp = genie.num_Data_RB(ind_RB);
			Data_temp = this->received_resource_grid(span(ind_RB * para.num_sc, (ind_RB + 1) * para.num_sc - 1), span::all, span(ind_port));
			Data_indexf_temp = Data_stream_indexf_temp(span(ind_RB * para.num_sc, (ind_RB + 1) * para.num_sc - 1), span::all, span(ind_port));

			this->Data_stream(span(Data_index, Data_index + length_temp - 1), span(ind_port)) = Data_temp(find(Data_position_temp));
			this->Data_stream_indexf(span(Data_index, Data_index + length_temp - 1), span(ind_port)) = Data_indexf_temp(find(Data_position_temp));
			Data_index = Data_index + length_temp;
		}
	}

	/* HV_estimated_data */
	this->HV_estimated_data = zeros<cx_cube>(para.num_Rx_antenna, genie.num_layers, data_length);
	temp_position = 0;
	for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
		temp_length = sum(sum(genie.Data_position(span(para.num_sc * ind_RB, para.num_sc * (ind_RB + 1) - 1), span(0, para.num_symb - 1))));
		HV_temp = this->HV_estimated_power(span(para.num_sc * ind_RB, para.num_sc *(ind_RB + 1) - 1), span(0, para.num_symb - 1), span::all);
		temp = zeros <cx_cube>(HV_temp.n_slices, HV_temp.n_rows, HV_temp.n_cols);
		for (uword i = 0; i < HV_temp.n_rows; i++)
			for (uword j = 0; j < HV_temp.n_cols; j++)
				for (uword k = 0; k < HV_temp.n_slices; k++)
					temp(k, i, j) = HV_temp(i, j, k);

		temp_position_matrix = genie.Data_position(span(para.num_sc * ind_RB, para.num_sc * (ind_RB + 1) - 1), span(0, para.num_symb - 1));

		for (uword i = 0; i < temp.n_rows; i++)
			for (uword k = 0, n = 0; k < temp.n_slices; k++)
				for (uword j = 0; j < temp.n_cols; j++)
					if (temp_position_matrix(j, k) == 1)
						this->HV_estimated_data(i, 0, temp_position + n++) = temp(i, j, k);

		temp_position = temp_position + temp_length;
	}

	/* local memeory deallocation */
	Data_stream_indexf_temp.reset();
	Data_position_temp.reset();
	Data_indexf_temp.reset();
	Data_temp.reset();
	temp.reset();
	HV_temp.reset();
	temp_position_matrix.reset();

}	/* end of method_Data_extraction() method definition */





/*
 * method_detection() - Detection
 */
void ModuleUE_SUMIMO::method_detection(ModuleParameterMIMO &para, Genie &genie)
{
	/* local variable declaration */
	irowvec			M_order = genie.m_order;
	irowvec			codeword1, codeword2;
	field<imat>		modulation_bittable;
	field<cx_mat>	modulation_alphabet;
	umat			indexf_temp;
	uword			mapping_index;
	cx_cube			H_estimated_temp_effective;
	cx_mat			H_estimated_temp_effective_temp;
	cx_mat			Detection_matrix;
	cx_mat			Data_stream_detected;
	cx_mat			H_g;
	rowvec			noise_enhancement;
	uword			temp_length, M_order_temp;
	mat				LLR_temp, LLR_temp_sub;
	mat				LLR_0, LLR_1;
	uword			temp_idx = 0;

	M_order = genie.m_order;
	codeword1 = zeros<irowvec>(genie.num_layers / 2);
	codeword2 = zeros<irowvec>(genie.num_layers - (genie.num_layers / 2));	// 1 - 0 + 1 = 2
	for (uword i = 0; i < codeword1.n_elem; i++)
		codeword1(i) = i;
	for (uword i = 0, n = floor(genie.num_layers / 2) + 1; i < codeword2.n_elem; i++, n++)
		codeword2(i) = n;

	switch (genie.num_codewords) {
	case 1:
		LLR_temp = zeros<mat>(M_order(0) * genie.num_layers, max(this->Data_stream.n_rows, this->Data_stream.n_cols));
		break;
	case 2:
		LLR_temp = zeros<mat>(M_order(0) * codeword1.n_elem + M_order(1) * codeword2.n_elem * codeword2.n_elem, max(Data_stream.n_rows, Data_stream.n_cols));
		break;
	}

	/* get CQI Table section :: Already done in module_Parameter */
	modulation_bittable = field<imat>(genie.num_codewords, 1);
	modulation_alphabet = field<cx_mat>(genie.num_codewords, 1);
	for (uword ind_codeword = 0; ind_codeword < genie.num_codewords; ind_codeword++) {
		modulation_bittable(ind_codeword) = modulation_mapping_bittable(M_order(ind_codeword) - 1);
		modulation_alphabet(ind_codeword) = modulation_mapping(M_order(ind_codeword) - 1);
	}

	for (uword ind_f = 0; ind_f < para.num_RB * para.num_sc; ind_f++) {
		indexf_temp = this->Data_stream_indexf.col(0) - ind_f;
		indexf_temp = matlnot(indexf_temp);
		H_estimated_temp_effective = this->HV_estimated_data.slices(find(indexf_temp));
		
		if (!H_estimated_temp_effective.is_empty()) {
			H_estimated_temp_effective_temp = H_estimated_temp_effective.slice(0);

			switch (para.ue_detection_mode) {
			case Ue_detection_mode::ZF:					/* Zero-forcing */
				Detection_matrix = pinv(H_estimated_temp_effective_temp);
				Data_stream_detected = Detection_matrix * this->Data_stream.rows(find(indexf_temp)).t();
				H_g = Detection_matrix * H_estimated_temp_effective_temp;
				noise_enhancement = sum(pow(abs(Detection_matrix), 2), 1);

				for (uword ind_layer = 0; ind_layer < genie.num_layers; ind_layer++) {
					LLR_temp_sub.clear();
					temp_length = Data_stream_detected.row(ind_layer).n_elem;
						
					if (genie.num_codewords == 1) {
						M_order_temp = M_order(0);
						mapping_index = 0;
					}
					else if (ind_layer <= floor(genie.num_layers / 2)) {
						M_order_temp = M_order(0);
						mapping_index = 0;
					}
					else {
						M_order_temp = M_order(1);
						mapping_index = 1;
					}

					for (uword ind_m = 0; ind_m < M_order_temp; ind_m++) {
						LLR_0 = exp(-1 / this->sigma_n_freq / noise_enhancement(ind_layer) * abs(repmat(Data_stream_detected.row(ind_layer), pow(2, (M_order_temp - 1)), 1) - repmat(modulation_alphabet(mapping_index)(find(modulation_bittable(mapping_index).row(ind_m) == 0)), 1, temp_length)));
						LLR_0 = sum(LLR_0, 0);
						LLR_1 = exp(-1 / this->sigma_n_freq / noise_enhancement(ind_layer) * abs(repmat(Data_stream_detected.row(ind_layer), pow(2, (M_order_temp - 1)), 1) - repmat(modulation_alphabet(mapping_index)(find(modulation_bittable(mapping_index).row(ind_m) == 1)), 1, temp_length)));
						LLR_1 = sum(LLR_1, 0);

						for (uword i = 0; i < LLR_0.n_elem; i++)
							if (LLR_0(i) == 0)
								LLR_0(i) = 1.0e-300;
						for (uword i = 0; i < LLR_1.n_elem; i++)
							if (LLR_1(i) == 0)
								LLR_1(i) = 1.0e-300;

						LLR_temp_sub = join_vert(LLR_temp_sub, log(LLR_1 / LLR_0));
					}

					uvec indexf_temp_temp = find(indexf_temp);
					for (uword i = M_order_temp * ind_layer, m = 0; i < M_order_temp * (ind_layer + 1) - 1; i++, m++)
						for (uword j = 0, n = 0; j < indexf_temp_temp.n_elem; j++, n++)
							LLR_temp(i, indexf_temp_temp(j)) = LLR_temp_sub(m, n);

				}
				break;
			case Ue_detection_mode::MMSE:				/* MMSE */
				Detection_matrix = H_estimated_temp_effective_temp.t() * pinv(H_estimated_temp_effective_temp * H_estimated_temp_effective_temp.t() + this->sigma_n_freq * eye<cx_mat>(para.num_Rx_antenna, para.num_Rx_antenna));
				
				if (norm(Detection_matrix) < 0.01)
					Detection_matrix = Detection_matrix / norm(Detection_matrix);

				Data_stream_detected = Detection_matrix * this->Data_stream.rows(find(indexf_temp)).st();
				H_g = Detection_matrix * H_estimated_temp_effective_temp;
				noise_enhancement = sum(pow(abs(Detection_matrix), 2), 1);

				for (uword ind_codeword = 0; ind_codeword < genie.num_codewords; ind_codeword++) {
					for (uword ind_layer = 0; ind_layer < genie.num_layers; ind_layer++) {
						LLR_temp_sub.clear();
						temp_length = Data_stream_detected.row(ind_layer).n_elem;
						if (genie.num_codewords == 1)
							M_order_temp = M_order(0);
						else if (ind_layer <= floor(genie.num_layers / 2))
							M_order_temp = M_order(0);
						else
							M_order_temp = M_order(1);

						for (uword ind_m = 0; ind_m < M_order_temp; ind_m++) {
							LLR_0 = exp(-1.0 / this->sigma_n_freq / noise_enhancement(ind_layer) * pow(abs(repmat(Data_stream_detected.row(ind_layer), pow(2, (M_order_temp - 1)), 1) - repmat(modulation_alphabet(ind_codeword)(find(modulation_bittable(ind_codeword).row(ind_m) == 0)), 1, temp_length)), 2.0));
							LLR_0 = sum(LLR_0, 0);
							LLR_1 = exp(-1.0 / this->sigma_n_freq / noise_enhancement(ind_layer) * pow(abs(repmat(Data_stream_detected.row(ind_layer), pow(2, (M_order_temp - 1)), 1) - repmat(modulation_alphabet(ind_codeword)(find(modulation_bittable(ind_codeword).row(ind_m) == 1)), 1, temp_length)), 2.0));
							LLR_1 = sum(LLR_1, 0);

							for (uword i = 0; i < LLR_0.n_elem; i++)
								if (LLR_0(i) == 0)
									LLR_0(i) = 1.0e-300;
							for (uword i = 0; i < LLR_1.n_elem; i++)
								if (LLR_1(i) == 0)
									LLR_1(i) = 1.0e-300;

							LLR_temp_sub = join_vert(LLR_temp_sub, log(LLR_1 / LLR_0));
						}
						
						uvec indexf_temp_temp = find(indexf_temp);
						for (uword i = M_order_temp * ind_layer, m = 0; i < M_order_temp * (ind_layer + 1); i++, m++)
							for (uword j = 0, n = 0; j < indexf_temp_temp.n_elem; j++, n++)
								LLR_temp(i, indexf_temp_temp(j)) = LLR_temp_sub(m, n);
							
					}
				}
				break;
			case Ue_detection_mode::Sphere_decoding:	/* Sphere decoding */
				cx_mat		Q, R;
				cx_mat		y_tilde_f, y_tilde_f_ZF;
				umat		argmin_ZF;
				cx_mat		ZF_x;
				cx_mat		LLR_temp2_sub2;
				cx_mat		y_tilde;
				mat			lambda_ML_bar;
				mat			ZF_distance, lambda_ML;
				imat		x_lambda_ML;
				uword		ind_layer_sd, index;
				umat		ind_order_sd;
				mat			DI;
				imat		x_sub;
				cx_mat		s_sub;
				cx_double	A_sub_max;


				qr(Q, R, H_estimated_temp_effective_temp);

				if (R.n_cols < R.n_rows) {
					R = R.rows(0, R.n_cols - 1);
					Q = Q.cols(0, R.n_cols - 1);
				}

				y_tilde_f = this->Data_stream.rows(find(indexf_temp));
				y_tilde_f_ZF = pinv(H_estimated_temp_effective_temp * y_tilde_f);

				ZF_distance = zeros<mat>(1, y_tilde_f.n_cols);
				argmin_ZF = zeros<umat>(genie.num_layers, y_tilde_f.n_cols);
				for (uword ind_data_f = 0; ind_data_f < y_tilde_f.n_cols; ind_data_f++)
					for (uword ind_layer = 0; ind_layer < genie.num_layers; ind_layer++)
						argmin_ZF(ind_layer, ind_data_f) = (pow(abs(y_tilde_f_ZF(ind_layer, ind_data_f) * ones<cx_mat>(1, pow(2, M_order(0))) - modulation_alphabet(0,0)), 2.0)).index_min();

				ZF_x = modulation_alphabet(0, 0)(find(argmin_ZF));
				ZF_distance = sum(pow(abs(y_tilde_f_ZF - ZF_x), 2.0), 0);

				LLR_temp2_sub2.reset();
				for (uword ind_iter = 0; ind_iter < max(indexf_temp.n_elem, this->Data_stream.n_cols); ind_iter++) {
					y_tilde = Q.t() * y_tilde_f.col(ind_iter);
					lambda_ML = ZF_distance(ind_iter);

					x_lambda_ML = reshape(modulation_bittable(0,0).cols(argmin_ZF.col(ind_iter)), modulation_bittable.n_elem, 1);
					lambda_ML_bar = 100.0 * ones<mat>(M_order(0) * genie.num_layers, 1);

					ind_layer_sd = genie.num_layers;
					ind_order_sd = ones<umat>(1, genie.num_layers);
					DI = zeros<mat>(1, genie.num_layers);

					x_sub = zeros<imat>(M_order(0) * genie.num_layers, 1);
					s_sub = zeros<cx_mat>(genie.num_layers, 1);

					while (true) {
						if (ind_order_sd(ind_layer_sd) <= pow(2, M_order(0))) {
							x_sub(span(M_order(0) * ind_layer_sd, M_order(0) * ind_layer_sd - 1), 0) = modulation_bittable(0, 0).col(ind_order_sd(ind_layer_sd));
							s_sub(ind_layer_sd) = modulation_alphabet(0, 0)(ind_order_sd(ind_layer_sd));

							DI(ind_layer_sd) = DI(ind_layer_sd) + pow(abs(y_tilde(ind_layer_sd) - sum(R(span(ind_layer_sd), span(ind_layer_sd, R.n_elem - 1)) * s_sub(span(ind_layer_sd, s_sub.n_elem - 1), 0))), 2);

							index = 0;
							A_sub_max = cx_double(0, 0);

							for (uword ind_layer = 0; ind_layer < genie.num_layers; ind_layer++) {
								for (uword ind_order = 0; ind_order < M_order(0); ind_order++) {
									if (ind_layer_sd > ind_order) {
										if (lambda_ML_bar(index) > A_sub_max.real())
											A_sub_max = lambda_ML_bar(index);
									}
									else {
										if (lambda_ML_bar(index) > A_sub_max.real() && x_lambda_ML(index) != x_sub(index))
											A_sub_max = lambda_ML_bar(index);
									}
									index += 1;
								}
							}

							if (A_sub_max.real() < DI(ind_layer_sd))
								ind_order_sd(ind_layer_sd) += 1;
							else {
								if (ind_layer_sd == 1) {
									if (DI(ind_layer_sd) < lambda_ML(0,0)) {
										lambda_ML_bar(x_lambda_ML != x_sub) = lambda_ML;
										lambda_ML = DI(ind_layer_sd);
										x_lambda_ML = x_sub;
									}
									else {
										for (uword i = 0; i < lambda_ML_bar.n_elem; i++)
											if (x_lambda_ML(i) != x_sub(i) && DI(ind_layer_sd) < lambda_ML_bar(i))
												lambda_ML_bar(i) = DI(ind_layer_sd);
									}
									ind_order_sd(ind_layer_sd) += 1;
								}
								else {
									ind_layer_sd -= 1;
								}
							}
						}
						else {
							if (ind_layer_sd == genie.num_layers)
								break;

							ind_order_sd(ind_layer_sd) = 1;
							ind_layer_sd += 1;
							ind_order_sd(ind_layer_sd) += 1;
						}

						if (ind_order_sd(0) == pow(2, M_order(0)))
							break;

					}	/* end of while (true) loop */

					LLR_temp_sub = zeros<mat>(M_order(0) * genie.num_layers, 1);
					LLR_temp_sub(x_lambda_ML == 0) = lambda_ML -lambda_ML_bar(x_lambda_ML == 0);
					LLR_temp_sub(x_lambda_ML == 1) = lambda_ML_bar(x_lambda_ML == 1) - lambda_ML;
					LLR_temp_sub = join_horiz(LLR_temp_sub, LLR_temp_sub);

				}	/* end of for() loop */

				LLR_temp.cols(indexf_temp) == LLR_temp_sub;

				/* local memory deallocation */
				Q.reset();
				R.reset();
				y_tilde_f.reset();
				y_tilde_f_ZF.reset();
				argmin_ZF.reset();
				ZF_x.reset();
				LLR_temp2_sub2.reset();
				y_tilde.reset();
				lambda_ML_bar.reset();
				lambda_ML.reset();
				ZF_distance.reset();
				x_lambda_ML.reset();
				ind_order_sd.reset();
				DI.reset();
				x_sub.reset();
				s_sub.reset();

				break;
			}	/* end of switch (para.ue_detection_mode) */
		}
	}
	this->LLR = LLR_temp;

	/* local memory deallocation */
	M_order.reset();
	codeword1.reset();
	codeword2.reset();
	for (uword i = 0; i < modulation_bittable.n_elem; i++)
		modulation_bittable(i).reset();
	modulation_bittable.reset();
	for (uword i = 0; i < modulation_alphabet.n_elem; i++)
		modulation_alphabet(i).reset();
	modulation_alphabet.reset();
	indexf_temp.reset();
	H_estimated_temp_effective.reset();
	H_estimated_temp_effective_temp.reset();
	Detection_matrix.reset();
	Data_stream_detected.reset();
	H_g.reset();
	noise_enhancement.reset();
	LLR_temp.reset();
	LLR_temp_sub.reset();
	LLR_0.reset();
	LLR_1.reset();

}	/* end of method_detection() method definition */





/*
 * method_delayer_mapping() - Delayer mapping
 *                            % TS. 36.211, V13.3.0, Table 6.3.3.2-1
 */
void ModuleUE_SUMIMO::method_delayer_mapping(Genie &genie)
{
	/* local variable declaration */
	urowvec		codeword1, codeword2;

	this->delayered_LLR = field<mat>(genie.num_codewords);
	switch (genie.num_codewords) {
	case 1:
		this->delayered_LLR(0) = reshape(this->LLR, genie.m_order(0), this->LLR.n_elem / genie.m_order(0));
		break;
	case 2:
		codeword1 = find(ones<urowvec>(floor(genie.num_layers / 2)));
		codeword2 = find(ones<urowvec>(genie.num_layers - genie.num_layers / 2));
		this->delayered_LLR(0) = reshape(this->LLR.rows(0, codeword1.n_elem * genie.m_order(0) - 1), 1, codeword1.n_elem * genie.m_order(0));
		this->delayered_LLR(1) = reshape(this->LLR.rows(codeword1.n_elem * genie.m_order(0), codeword1.n_elem * genie.m_order(0) + codeword2.n_elem * genie.m_order(1)), 1, (codeword1.n_elem * genie.m_order(0) + codeword2.n_elem * genie.m_order(1)) - codeword1.n_elem * genie.m_order(0));
		break;
	}

	/* local memory deallocation */
	codeword1.reset();
	codeword2.reset();

}	/* end of method_delayer_mapping method definition */





/*
 * method_descrambling() - Descrambling
 *                         TS 38.211, 7.3.1.1.,  (TS 36.211, 6.3.1)
 */
void ModuleUE_SUMIMO::method_descrambling(ModuleParameterMIMO &para, Genie &genie, int ind_slot)
{
	/* local variable declaration */
	uword		n_RNTI = 1;
	uword		c_init;
	mat			scrambling_seq;

	this->cell_ID = para.cell_ID;
	this->descrambled_LLR = field<mat>(genie.num_codewords);

	for (uword ind_codeword = 0; ind_codeword < genie.num_codewords; ind_codeword++) {
		this->descrambled_LLR(ind_codeword) = reshape(this->delayered_LLR(ind_codeword), 1, this->delayered_LLR(ind_codeword).n_elem);
		c_init = n_RNTI * pow(2, 14) + (ind_codeword + 1) * pow(2, 13) + ind_slot * pow(2, 9) + this->cell_ID;
		func_pseudo_sequence_generation(c_init, this->descrambled_LLR(ind_codeword).n_elem, scrambling_seq);
		
		for (uword i = 0; i < scrambling_seq.n_elem; i++)
			if (scrambling_seq(i) == 0)
				scrambling_seq(i) = -1.0;

		this->descrambled_LLR(ind_codeword) = -1.0 * this->descrambled_LLR(ind_codeword) % scrambling_seq;
	}

	/* local memory deallocation */
	scrambling_seq.reset();

}	/* end of method_descrambling() method definition */






/*
 * method_HARQ_process() - HARQ process
 */
void ModuleUE_SUMIMO::method_HARQ_process(ModuleParameterMIMO &para, Genie &genie)
{
	if (para.MUMIMO)
		this->descrambled_LLR_HARQ = this->descrambled_LLR;
	else {
		for (uword ind_codeword = 0; ind_codeword < genie.num_codewords; ind_codeword++) {
			if (genie.newdata_indicator(ind_codeword) == 0) {
				if (this->descrambled_LLR(ind_codeword).n_elem > this->HARQ_buffer.LLR(ind_codeword, genie.HARQ_process_index(0)).n_elem)
					this->descrambled_LLR_HARQ(ind_codeword) = this->HARQ_buffer.LLR(ind_codeword, genie.HARQ_process_index(0)) + this->descrambled_LLR(ind_codeword)(span(0, this->HARQ_buffer.LLR(genie.HARQ_process_index(0)).n_elem), span::all);
				else
					this->descrambled_LLR_HARQ(ind_codeword) = this->HARQ_buffer.LLR(ind_codeword, genie.HARQ_process_index(0)) + join_horiz(this->descrambled_LLR(ind_codeword), zeros(1, this->HARQ_buffer.LLR(ind_codeword, genie.HARQ_process_index(0)).n_elem - this->descrambled_LLR(ind_codeword).n_elem));

				this->HARQ_buffer.retransmission(ind_codeword, genie.HARQ_process_index(0)) = this->HARQ_buffer.retransmission(ind_codeword, genie.HARQ_process_index(0));
				this->HARQ_buffer.LLR(ind_codeword, genie.HARQ_process_index(0)) = this->descrambled_LLR_HARQ(ind_codeword);
			}
			else {
				this->HARQ_buffer.LLR(ind_codeword, genie.HARQ_process_index(0)) = this->descrambled_LLR(ind_codeword);
				this->descrambled_LLR_HARQ(ind_codeword) = this->descrambled_LLR(ind_codeword);
			}
		}
	}

}	/* end of method_HARQ_process() method definition */





/*
 * method_decoding() - Decoding
 *                     For uncoded BER
 */
void ModuleUE_SUMIMO::method_decoding(ModuleParameterMIMO &para, NR_ChannelCoding &LDPC, Genie &genie, int ind_slot)
{
	/* local variable declaration */
	int		ACK;

	this->coded_bit = field<umat>(para.num_codewords);
	this->decoded_bit = field<umat>(para.num_codewords);

	for (uword ind_codeword = 0; ind_codeword < para.num_codewords; ind_codeword++) {
		this->coded_bit(ind_codeword) = this->descrambled_LLR(ind_codeword) > 0;
		genie.ch_coding(ind_codeword).symbol_LLR = zeros<rowvec>(this->descrambled_LLR(ind_codeword).n_elem);
		genie.ch_coding(ind_codeword).symbol_LLR = this->descrambled_LLR(ind_codeword);
		
		LDPC_C_rx(LDPC, genie.ch_coding(ind_codeword));

		this->decoded_bit(ind_codeword) = genie.ch_coding(ind_codeword).rx_a;
		ACK = genie.ch_coding(ind_codeword).ACK_for_check;
		this->feedBack.ACK = ACK;
		this->feedBack.NACK_a = genie.ch_coding(ind_codeword).NACK_a;
		this->feedBack.NACK_b = genie.ch_coding(ind_codeword).NACK_b;
		this->result.ACK(ind_slot) = this->result.ACK(ind_slot) + ACK;

		if (ACK || this->HARQ_buffer.retransmission(genie.HARQ_process_index(0)) == this->HARQ_buffer.max_retransmission) {
			this->HARQ_buffer.num_TB = this->HARQ_buffer.num_TB + 1;
			this->HARQ_buffer.retransmission(genie.HARQ_process_index(0)) = 0;
			this->HARQ_buffer.LLR(genie.HARQ_process_index(0)).reset();
		}
		else {
			this->HARQ_buffer.BS_signaling(genie.HARQ_process_index(0)) = genie.ch_coding(ind_codeword);
			this->HARQ_buffer.LLR(genie.HARQ_process_index(0)) = this->descrambled_LLR_HARQ(0);
		}
	}

}	/* end of method_decoding() method definition */





/*
 * method_Result() - Result
 */
void ModuleUE_SUMIMO::method_Result(ModuleParameterMIMO &para, Genie &genie, int ind_slot) {

	/* local variable definition */
	mat			sig_power;
	mat			interf_power;
	colvec		noise_enhancement;
	cx_mat		H_estimated_temp, sq_H_estimated_temp;
	cx_mat		Detection_matrix;
	cx_mat		H_effective_temp;

	if (para.channel.time_fading == Ch_mode::Block_Fading) {
		sig_power = zeros<mat>(para.num_sc * para.num_RB, 1);
		interf_power = zeros<mat>(para.num_sc * para.num_RB, 1);
		noise_enhancement = zeros<colvec>(para.num_sc * para.num_RB, 1);

		for (uword ind_freq = 0; ind_freq < para.num_sc * para.num_RB; ind_freq++) {
			sq_H_estimated_temp = this->H_estimated(ind_freq)(span(0), span::all, span::all);
			if (para.num_Rx_antenna == 1)
				H_estimated_temp = sq_H_estimated_temp.st() * genie.Precoder(floor(ind_freq / para.num_sc));
			else
				H_estimated_temp = sq_H_estimated_temp * genie.Precoder(floor(ind_freq / para.num_sc));

			Detection_matrix = pinv(H_estimated_temp);
			H_effective_temp = Detection_matrix * H_estimated_temp;
			sig_power(ind_freq) = pow(abs(H_effective_temp(0, 0)), 2);
			interf_power(ind_freq) = sum(pow(abs(H_effective_temp.row(0)), 2)) - sig_power(ind_freq);
			noise_enhancement(ind_freq) = sum(pow(abs(Detection_matrix.row(0)), 2));
		}

		/* Result */
		this->result.received_SNR_for_cali = mean(sig_power / (interf_power + para.sigma_n_freq * noise_enhancement));
		this->result.received_SNR(ind_slot) = (10.0 * log10((sig_power / (interf_power + para.sigma_n_freq * noise_enhancement)))).t();
	}

	int same_bit_cnt = 0, diff_bit_cnt = 0;
	for (uword ind_codeword = 0; ind_codeword < genie.num_codewords; ind_codeword++) {

		if (this->result.ACK(ind_slot))
			this->result.suc_data(ind_slot) = para.num_data_bits(0);
		else
			this->result.suc_data(ind_slot) = 0;

		same_bit_cnt = 0;
		for (uword i = 0; i < genie.coded_bit_stream(ind_codeword).n_elem; i++)
			if (genie.coded_bit_stream(ind_codeword)(i) == this->coded_bit(ind_codeword)(i))
				same_bit_cnt++;
		this->result.suc_bit(ind_slot) = same_bit_cnt;

		same_bit_cnt = 0;
		for (uword i = 0; i < genie.bit_stream(ind_codeword).n_elem; i++)
			if (genie.bit_stream(ind_codeword)(i) == this->decoded_bit(ind_codeword)(i))
				same_bit_cnt++;
			else
				diff_bit_cnt++;

		this->result.suc_bit_data(ind_slot) = same_bit_cnt;

		this->result.Totalbit(ind_slot) = genie.coded_bit_stream(ind_codeword).n_elem;
	}

	/* local memory deallocation */
	sig_power.reset();
	interf_power.reset();
	noise_enhancement.reset();
	H_estimated_temp.reset();
	sq_H_estimated_temp.reset();
	Detection_matrix.reset();
	H_effective_temp.reset();

}	/* end of method_Result() method definition */




/* default constructor & destructor */
ModuleUE_SUMIMO::ModuleUE_SUMIMO()
{
}

ModuleUE_SUMIMO::~ModuleUE_SUMIMO()
{
}



/* end of module_UE_MIMO.cpp */