/*
 * module_Parameter_MIMO.cpp
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *		Last Updated: 2018-09-19
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */

#include <armadillo>
#include "module_Parameter_MIMO.h"
#include "getCQITable.h"


using namespace std;
using namespace arma;



/* default constructor and destructor */
ModuleParameterMIMO::ModuleParameterMIMO(void) {}
ModuleParameterMIMO::~ModuleParameterMIMO(void) {}



/*
 * module_Parameter_MIMO() - Custom constructor of module_Parameter class.
s */
void ModuleParameterMIMO::module_Parameter_MIMO(const ParameterInit &init)
{
	/* local variable declaration */
	uword	i;				// general index
	int		L;				// variable for this->num_SS_block

	rowvec	PDP_Delay;
	rowvec	New_PDP_Delay_samples;
	uword	diag_index;

	/* default value settings */
	this->CSIRS_bitmat_value << 1 << 2 << 3 << 4 << 5 << 6 << endr;
	this->CSIRS_symbol_start_idx << 0 << 5 << endr;

	this->pol_slant_angle << -45 << 45 << endr;
	this->subband_set << 6 << 4 << 4 << 6 << 8 << 8 << endr;

	/* channel parameters */
	this->theta_v = 45.0;
	this->phi_v = 45.0;

	/* Parameters determined by scenario */
	this->Sim_Case = init.Sim_Case;
	this->BW = init.BW;

	if (init.SNR_range_isDirty)
		this->SNR_range = init.SNR_range;

	this->num_subframe = init.num_subframe;
	this->f = init.f;
	this->mu = init.mu;
	this->subcarrier_spacing = pow(2, this->mu) * 15e3;
	this->num_symb_in_subframe = (int)(pow(2, this->mu) * 14);
	this->num_slot_in_subframe = (int)pow(2, this->mu);

	switch (init.ch_mode) {
	case ch_mode::Block_Fading:
		this->t = zeros<rowvec>(1);
		this->t(0) = 0.0;
		break;
	case ch_mode::Fast_Fading:
		this->t = zeros<rowvec>(this->num_symb);
		for (i = 0; i < this->num_symb; i++)
			this->t(i) = i;
		this->t = this->t * (1e-3 / pow(2.0, this->mu));	// Time for 14 OFDM slots of one subframe(1ms); subject to change acoording to numerology.
		break;
	default:
		;
	}   /* end of switch(init.ch_mode) */

	/* TS 38.104 table 5.3.2 */
	switch (this->mu) {
	case 0:
		this->BW_set << 5 << 10 << 15 << 20 << 25 << 30 << 40 << 50 << endr;
		this->BW_set = this->BW_set * 1e6;
		this->RB_set << 25 << 52 << 79 << 106 << 133 << 160 << 216 << 270 << endr;
		break;
	case 1:
		this->BW_set << 5 << 10 << 15 << 20 << 25 << 30 << 40 << 50 << 60 << 70 << 80 << 90 << 100 << endr;
		this->BW_set *= 1e6;
		this->RB_set << 11 << 24 << 38 << 51 << 65 << 78 << 106 << 133 << 162 << 189 << 217 << 245 << 273 << endr;
		break;
	case 2:
		if (this->f < 6e9) {	// FR 1
			this->BW_set << 10 << 15 << 20 << 25 << 30 << 40 << 50 << 60 << 70 << 80 << 90 << 100 * 1e6 << endr;
			this->BW_set *= 1e6;
			this->RB_set << 11 << 18 << 24 << 31 << 38 << 51 << 65 << 79 << 93 << 107 << 121 << 135 << endr;
		}
		else {					// FR 2
			this->BW_set << 50 << 100 << 200 << endr;
			this->BW_set *= 1e6;
			this->RB_set << 66 << 132 << 264 << endr;
		}
		break;
	case 3:
		this->BW_set << 50 << 100 << 200 << 400 << endr;
		this->BW_set *= 1e6;
		this->RB_set << 32 << 66 << 132 << 264 << endr;
		break;
	default:
		;
	}   /* end of switch(this->mu) */

	for (i = 0; i < this->BW_set.n_elem; i++)			// Allowed values: Depending on numerology
		if (this->BW_set(i) == this->BW)
			break;
	this->num_RB_maximum = this->RB_set(i);
	this->num_RB_maximum_isDirty = true;				// set dirty bit of num_RB_maximum

	for (i = 0; this->BW_set.n_elem; i++)				// Depending on numerology
		if (this->BW_set(i) == this->BW)
			break;
	this->num_RB = this->RB_set(i);
	this->num_RB_isDirty = true;

	if (init.num_RB_isDirty) {
		this->num_RB = init.num_RB;
		this->num_RB_isDirty = true;
	}
	else if (this->use_fullband) {
		this->num_RB = this->num_RB_maximum;
		this->num_RB_isDirty = true;
	}
	else
		;	// empty statement

	if (!this->num_RB_isDirty)
		printf("ERROR: Input number of resource block for data transmission\n");


	switch (this->mu) {
	case 0:
		if (this->f <= 3e9) {
			this->num_SS_block = 2 * 2;
			this->slotset_SS_block << 1 << 2 << endr;
		}
		else if (this->f > 3e9 && this->f <= 6e9) {
			this->num_SS_block = 2 * 4;
			this->slotset_SS_block << 1 << 2 << 3 << 4 << endr;
		}
		break;
	case 1:
		if (this->f <= 3e9) {
			this->num_SS_block = 4 * 1;
			this->slotset_SS_block << 1 << 2 << endr;
		}
		else if (this->f > 3e9 && this->f <= 6e9) {
			this->num_SS_block = 4 * 2;
			this->slotset_SS_block << 1 << 2 << 3 << 4 << endr;
		}
		break;
	case 2:
		if (this->f <= 3e9) {
			this->num_SS_block = 2 * 2;
			this->slotset_SS_block << 1 << 2 << endr;
		}
		else if (this->f > 3e9 && this->f <= 6e9) {
			this->num_SS_block = 2 * 4;
			this->slotset_SS_block << 1 << 2 << 3 << 4 << endr;
		}
		break;
	case 3:
		if (this->f >= 6e9) {
			this->num_SS_block = 4 * 16;
			this->slotset_SS_block << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 11 << 12 << 13 << 14 << 15 << 16 << 17 << 18
				<< 21 << 22 << 23 << 24 << 25 << 26 << 27 << 28 << 31 << 32 << 33 << 34 << 35 << 36 << 37 << 38 << endr;
		}
		break;
	case 4:
		if (this->f >= 6e9) {
			this->num_SS_block = 8 * 8;
			this->slotset_SS_block << 1 << 2 << 3 << 6 << 7 << 8 << endr;
		}
	default:
		;
	}   /* end of switch (this->mu) */

	if (this->BW == 15e6 && this->subcarrier_spacing == 15e3)
		this->size_fft = 1536;
	else if (this->BW == 15e6 && this->subcarrier_spacing == 7.5e3)
		this->size_fft = 1536 * 2;
	else if (this->num_RB == 1 && init.ch_type == ch_type::CDL_A)
		this->size_fft = 31;
	else
		this->size_fft = pow(2, ceil(log2(this->num_RB * this->num_sc)));


	/* Sampling frequency */
	this->Fs = this->size_fft * this->subcarrier_spacing;

	this->Tx_pol = init.Tx_pol;		// Tx polarization
	this->Rx_pol = init.Rx_pol;		// Rx polarization

	if (init.num_Tx_antenna_horizontal_isExist && init.num_Tx_antenna_vertical_isExist) {
		this->N1 = init.num_Tx_antenna_horizontal;
		this->N2 = init.num_Tx_antenna_vertical;
		this->num_Tx_antenna = init.num_Tx_antenna_horizontal * init.num_Tx_antenna_vertical * this->Tx_pol;
		this->TxArrayType = TxArrayType::URA;
	}
	else if (init.num_Tx_antenna_isExist) {
		this->num_Rx_antenna = init.num_Rx_antenna;
		this->TxArrayType = TxArrayType::ULA;
	}
	else
		perror("Input the number of transmit antennas.\n");


	if (init.num_Tx_antenna_horizontal_isExist && init.num_Tx_antenna_vertical_isExist) {
		this->M1 = init.num_Rx_antenna_horizontal;
		this->M2 = init.num_Rx_antenna_vertical;
		this->num_Rx_antenna = init.num_Rx_antenna_horizontal * init.num_Rx_antenna_vertical * this->Rx_pol;
		this->RxArrayType = RxArrayType::URA;
	}
	else if (init.num_Rx_antenna_isDirty) {
		this->num_Rx_antenna = init.num_Rx_antenna * this->Rx_pol;
		this->RxArrayType = RxArrayType::ULA;
	}
	else
		perror("Input the number of Receive antennas\n");


	if (init.cal_scenario_isExist && init.cal_scenario) {
		this->cal_scenario = true;
		this->N1 = 4;
		this->N2 = 4;
		this->num_Tx_antenna = this->Tx_pol * this->N1 * this->N2;
		this->num_Rx_antenna = this->Rx_pol;
		this->scenario_beam = 1;
		this->DMRS_type = 1;
		this->Tx_downtilt = 0;
	}
	else {
		L = this->num_SS_block;

		if (this->N2 > 1) {
			this->L1 = L / 2;
			this->L2 = 2;
		}
		else {
			this->L1 = L;
			this->L2 = 1;
		}

		if (this->N1 > 4)
			this->O1 = 4;
		else
			this->O1 = 2;

		if (this->N2 > 2)
			this->O2 = 2;
		else
			this->O2 = 1;

		this->S1 = this->N1 * this->O1 / this->L1;
		this->S2 = this->N2 * this->O2 / this->L2;

		this->P1 = max(this->S1, 2 * this->O1);
		this->P2 = max(this->S2, 2 * this->O2);

		if ((this->P1 * this->P2) > 16)
			perror("This antenna configuration is not supported.\n");

		if ((this->P1 * this->P2) <= 2)
			this->P1 = 2 * this->P1;

	}   /* end of if (init.cal_scenario_isExist && init.cal_scenario) */

	switch (init.ch_type) {
	case ch_type::TDL_A:
	case ch_type::TDL_B:
	case ch_type::TDL_C:
	case ch_type::TDL_D:
	case ch_type::TDL_E:
	case ch_type::CDL_A:
	case ch_type::CDL_B:
	case ch_type::CDL_C:
	case ch_type::CDL_D:
	case ch_type::CDL_E:
		this->DS = 100.0 * 1e-9;
	}

	/* RS initial setting */
	this->PDSCH_mapping_type = PDSCH_mapping_type::A;
	this->DL_DMRS_config_type = DL_DMRS_config_type::type_2;
	this->DL_DMRS_typeA_pos = 3;
	this->DL_DMRS_max_len = 2;
	this->pos_last_PDSCH_symbol = 12;
	this->DL_DMRS_add_pos = 1;
	this->CSIRS_periodicity = 5;

	switch (this->Sim_Case) {
	case Sim_Case::SUSISO:
		this->ch_est_mode_DMRS = init.ch_est_mode_DMRS;
		this->row_CSIRS = init.CSIRS_row;
		this->CSIRS_density = init.CSIRS_density;
		this->CSIRS_periodicity = init.CSIRS_periodicity;
		this->DMRS_type = init.DMRS_config;
		this->PDSCH_mapping_type = init.PDSCH_type;

		this->port_set << 1000 << endr;
		this->port_set_DMRS << 1000 << endr;
		this->port_set_CSIRS << 3000 << endr;
										
		this->CSIRS_CDMType = CSIRS_CDMType::NO_CDM;
		this->row_CSIRS = 1;		// 1 ~ 19, TS 38.211 'Table 7.4.1.5.2-1: CSI-RS locations within a slot.'

		/* power */
		this->power_DMRS << 1.3 << endr;
		this->power_CSIRS << 1.4 << endr;
		this->power_data << 1.0 << endr;

		break;
	case Sim_Case::SUMIMO:
		this->port_set_DMRS << 1000 << 1001 << 1002 << 1003 << 1004 << 1005 << 1006 << 1007 << 1008 << 1009 << 1010 << 1011 << endr;
		this->port_set_CSIRS = zeros<irowvec>(this->Tx_pol * this->P1 * this->P2);
		for (i = 0; i < this->Tx_pol * this->P1 * this->P2; i++)
			this->port_set_CSIRS(i) = i + 3000;

		switch (this->port_set_CSIRS.n_elem) {
		case 4:
			this->row_CSIRS = 4;
			break;
		case 8:
			this->row_CSIRS = 7;
			break;
		case 12:
			this->row_CSIRS = 10;
			break;
		case 16:
			this->row_CSIRS = 12;
			break;
		case 24:
			this->row_CSIRS = 14;
			break;
		case 32:
			this->row_CSIRS = 17;
			break;
		}

		break;
	}   /* end of switch(this->Sim_Case) */

	this->num_port_DMRS  = this->port_set_DMRS.n_elem;
	this->num_port_CSIRS = this->port_set_CSIRS.n_elem;

	/* Calculation */
	method_calculation();			// invoke class member method

	/* Setting for channel model */
	this->Channel.type = init.ch_type;
	this->Channel.time_fading = init.ch_mode;
	this->Channel.time_correlation = time_correlation::independent;		// corelated or independent

	this->Cov_Matrix_time = zeros<mat>(this->size_fft, this->size_fft);
	switch (this->Channel.type) {
	case ch_type::TDL_A:	// Table 7.7.2-1
		TDL << 0.0000 << -13.4 << endr
			<< 0.3819 <<   0.0 << endr
			<< 0.4025 <<  -2.2 << endr
			<< 0.5868 <<  -4.0 << endr
			<< 0.4610 <<  -6.0 << endr
			<< 0.5375 <<  -8.2 << endr
			<< 0.6708 <<  -9.9 << endr
			<< 0.5750 << -10.5 << endr
			<< 0.7618 <<  -7.5 << endr
			<< 1.5375 << -15.9 << endr
			<< 1.8978 <<  -6.6 << endr
			<< 2.2242 << -16.7 << endr
			<< 2.1718 << -12.4 << endr
			<< 2.4942 << -15.2 << endr
			<< 2.5119 << -10.8 << endr
			<< 3.0582 << -11.3 << endr
			<< 4.0810 << -12.7 << endr
			<< 4.4579 << -16.2 << endr
			<< 4.5695 << -18.3 << endr
			<< 4.7966 << -18.9 << endr
			<< 5.0066 << -16.6 << endr
			<< 5.3043 << -19.9 << endr
			<< 9.6586 << -29.7 << endr;
		this->Channel.PDP_dB = zeros(arma::size(TDL.t()));
		this->Channel.PDP_dB = join_horiz(TDL.col(1), (TDL.col(0) * this->DS)).t();

		PDP_Delay = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		PDP_Delay = this->Channel.PDP_dB.row(1);
		New_PDP_Delay_samples = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		New_PDP_Delay_samples = round(PDP_Delay * this->Fs);
		for (uword ind = 0; ind < New_PDP_Delay_samples.n_elem; ind++) {
			diag_index = New_PDP_Delay_samples(ind);
			this->Cov_Matrix_time(diag_index, diag_index) = this->Cov_Matrix_time(diag_index, diag_index) + pow(10.0, this->Channel.PDP_dB(0, ind) / 10.0);
		}
		break;
	case ch_type::TDL_B:	// Table 7.7.2-2
		TDL << 0.0000 <<   0.0 << endr
			<< 0.1072 <<  -2.2 << endr
			<< 0.2155 <<  -4.0 << endr
			<< 0.2095 <<  -3.2 << endr
			<< 0.2870 <<  -9.8 << endr
			<< 0.2986 <<  -1.2 << endr
			<< 0.3752 <<  -3.4 << endr
			<< 0.5055 <<  -5.2 << endr
			<< 0.3681 <<  -7.6 << endr
			<< 0.3697 <<  -3.0 << endr
			<< 0.5700 <<  -8.9 << endr
			<< 0.5283 <<  -9.0 << endr
			<< 1.1021 <<  -4.8 << endr
			<< 1.2756 <<  -5.7 << endr
			<< 1.5474 <<  -7.5 << endr
			<< 1.7842 <<  -1.9 << endr
			<< 2.0169 <<  -7.6 << endr
			<< 2.8294 << -12.2 << endr
			<< 3.0219 <<  -9.8 << endr
			<< 3.6187 << -11.4 << endr
			<< 4.1067 << -14.9 << endr
			<< 4.2790 <<  -9.2 << endr
			<< 4.7834 << -11.3 << endr;
		this->Channel.PDP_dB = zeros(arma::size(TDL.t()));
		this->Channel.PDP_dB = join_horiz(TDL.col(1), (TDL.col(0) * this->DS)).t();

		PDP_Delay = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		PDP_Delay = this->Channel.PDP_dB.row(1);
		New_PDP_Delay_samples = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		New_PDP_Delay_samples = round(PDP_Delay * this->Fs);
		for (uword ind = 0; ind < New_PDP_Delay_samples.n_elem; ind++) {
			diag_index = New_PDP_Delay_samples(ind);
			this->Cov_Matrix_time(diag_index, diag_index) = this->Cov_Matrix_time(diag_index, diag_index) + pow(10.0, this->Channel.PDP_dB(0, ind) / 10.0);
		}
		break;
	case ch_type::TDL_C:		// Table 7.7.2-3
		TDL << 0.0000 <<  -4.4 << endr
			<< 0.2099 <<  -1.2 << endr
			<< 0.2219 <<  -3.5 << endr
			<< 0.2329 <<  -5.2 << endr
			<< 0.2176 <<  -2.5 << endr
			<< 0.6366 <<   0.0 << endr
			<< 0.6448 <<  -2.2 << endr
			<< 0.6560 <<  -3.9 << endr
			<< 0.6584 <<  -7.4 << endr
			<< 0.7935 <<  -7.1 << endr
			<< 0.8213 << -10.7 << endr
			<< 0.9336 << -11.1 << endr
			<< 1.2285 <<  -5.1 << endr
			<< 1.3083 <<  -6.8 << endr
			<< 2.1704 <<  -8.7 << endr
			<< 2.7105 << -13.2 << endr
			<< 4.2589 << -13.9 << endr
			<< 4.6003 << -13.9 << endr
			<< 5.4902 << -15.8 << endr
			<< 5.6077 << -17.1 << endr
			<< 6.3065 << -16.0 << endr
			<< 6.6374 << -15.7 << endr
			<< 7.0427 << -21.6 << endr
			<< 8.6523 << -22.8 << endr;
		this->Channel.PDP_dB = zeros(arma::size(TDL.t()));
		this->Channel.PDP_dB = join_horiz(TDL.col(1), (TDL.col(0) * this->DS)).t();

		PDP_Delay = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		PDP_Delay = this->Channel.PDP_dB.row(1);
		New_PDP_Delay_samples = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		New_PDP_Delay_samples = round(PDP_Delay * this->Fs);
		for (size_t ind = 0; ind < New_PDP_Delay_samples.n_elem; ind++) {
			diag_index = New_PDP_Delay_samples(ind);
			this->Cov_Matrix_time(diag_index, diag_index) = this->Cov_Matrix_time(diag_index, diag_index) + pow(10.0, this->Channel.PDP_dB(0, ind) / 10.0);
		}
		break;
	case ch_type::TDL_D:		// Table 7.7.2-4
		TDL << 0.0000 <<  -0.2 << endr
			<< 0.0350 << -18.8 << endr
			<< 0.6120 << -21.0 << endr
			<< 1.3630 << -22.8 << endr
			<< 1.4050 << -17.9 << endr
			<< 1.8040 << -20.1 << endr
			<< 2.5960 << -21.9 << endr
			<< 1.7750 << -22.9 << endr
			<< 4.0420 << -27.8 << endr
			<< 7.9370 << -23.6 << endr
			<< 9.4240 << -24.8 << endr
			<< 9.7080 << -30.0 << endr
			<< 12.5250 << -27.7 << endr;
		this->Channel.PDP_dB = zeros(arma::size(TDL.t()));
		this->Channel.PDP_dB = join_horiz(TDL.col(1), (TDL.col(0) * this->DS)).t();

		PDP_Delay = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		PDP_Delay = this->Channel.PDP_dB.row(1);
		New_PDP_Delay_samples = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		New_PDP_Delay_samples = round(PDP_Delay * this->Fs);
		for (size_t ind = 0; ind < New_PDP_Delay_samples.n_elem; ind++) {
			diag_index = New_PDP_Delay_samples(ind);
			this->Cov_Matrix_time(diag_index, diag_index) = this->Cov_Matrix_time(diag_index, diag_index) + pow(10.0, this->Channel.PDP_dB(0, ind) / 10.0);
		}
		break;
	case ch_type::TDL_E:		// Table 7.7.2-5
		TDL << 0.0000  << -0.03 << endr
			<< 0.5133  << -15.8 << endr
			<< 0.5440  << -18.1 << endr
			<< 0.5630  << -19.8 << endr
			<< 0.5440  << -22.9 << endr
			<< 0.7112  << -22.4 << endr
			<< 1.9092  << -18.6 << endr
			<< 1.9293  << -20.8 << endr
			<< 1.9589  << -22.6 << endr
			<< 2.6426  << -22.3 << endr
			<< 3.7136  << -25.6 << endr
			<< 5.4524  << -20.2 << endr
			<< 12.0034 << -29.8 << endr
			<< 20.6519 << -29.2 << endr;
		this->Channel.PDP_dB = zeros(arma::size(TDL.t()));
		this->Channel.PDP_dB = join_horiz(TDL.col(1), (TDL.col(0) * this->DS)).t();

		PDP_Delay = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		PDP_Delay = this->Channel.PDP_dB.row(1);
		New_PDP_Delay_samples = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		New_PDP_Delay_samples = round(PDP_Delay * this->Fs);
		for (size_t ind = 0; ind < New_PDP_Delay_samples.n_elem; ind++) {
			diag_index = New_PDP_Delay_samples(ind);
			this->Cov_Matrix_time(diag_index, diag_index) = this->Cov_Matrix_time(diag_index, diag_index) + pow(10.0, this->Channel.PDP_dB(0, ind) / 10.0);
		}
		break;
	case ch_type::PedA:
		this->Channel.PDP_dB << 0.0 << -9.7 << -19.2 << -22.8 << endr
			<< 0.0 << 110.0 * pow(10.0, -9.0) << 190.0 * pow(10.0, -9.0) << 410.0 * pow(10.0, -9.0) << endr;
		this->Cov_Matrix_time(0, 0) = 1.0 + pow(pow(10.0, -2.28), 2.0) + pow(pow(10.0, -1.92), 2.0);
		this->Cov_Matrix_time(1, 1) = pow(pow(10.0, -2.28), 2.0);
		break;
	case ch_type::flat_Rayleigh:
	case ch_type::Rayleigh2:
		this->Channel.PDP_dB << 0.0 << endr
			<< 0.0 << endr;
		this->Cov_Matrix_time(0, 0) = 1.0;
		break;
	case ch_type::PedB:
	case ch_type::PedBcorr:
		this->Channel.PDP_dB << 0.0 << -0.9 << -4.9 << -8.0 << -7.8 << -23.9 << endr
			<< 0.0 << 200.0 * pow(10.0, -9.0) << 800.0 * pow(10.0, -9.0) << 1200.0 * pow(10.0, -9.0) << 2300.0 * pow(10.0, -9.0) << 3700.0 * pow(10.0, -9.0) << endr;
		this->Cov_Matrix_time(0, 0) = 1.0 + pow(10.0, -0.09);
		this->Cov_Matrix_time(2, 2) = pow(10.0, -0.49) + pow(10.0, -0.8);
		this->Cov_Matrix_time(4, 4) = pow(10.0, -0.78);
		this->Cov_Matrix_freq(7, 7) = pow(10.0, -2.39);
		break;
	case ch_type::VehA:
		this->Channel.PDP_dB << 0.0 << -1.0 << -9.0 << -10.0 << -15.0 << -20.0 << endr
			<< 0.0 << 310.0 * pow(10.0, -9.0) << 710.0 * pow(10.0, -9.0) << 1090.0 * pow(10.0, -9.0) << 1730.0 * pow(10.0, -9.0) << 2510.0 * pow(10.0, -9.0) << endr;
		this->Cov_Matrix_time(0, 0) = 1.0;
		this->Cov_Matrix_time(1, 1) = pow(10.0, -0.1) + pow(10.0, -0.9);
		this->Cov_Matrix_time(2, 2) = pow(10.0, -1.0);
		this->Cov_Matrix_time(3, 3) = pow(10.0, -1.5);
		this->Cov_Matrix_time(5, 5) = pow(10.0, -2.0);
		break;
	case ch_type::VehB:
		this->Channel.PDP_dB << -2.5 << 0.0 << -12.8 << -10.0 << -25.2 << -16.0 << endr
			<< 0.0 << 300.0 * pow(10.0, -9.0) << 8900.0 * pow(10.0, -9.0) << 12900.0 * pow(10.0, -9.0) << 17100.0 * pow(10.0, -9.0) << 20000.0 * pow(10.0, -9.0) << endr;
		this->Cov_Matrix_time(0, 0) = pow(10.0, -0.25);
		this->Cov_Matrix_time(1, 1) = 1.0;
		this->Cov_Matrix_time(17, 17) = pow(10.0, -1.28);
		this->Cov_Matrix_time(25, 25) = pow(10.0, -1.0);
		this->Cov_Matrix_time(33, 33) = pow(10.0, -2.51);
		this->Cov_Matrix_time(38, 38) = pow(10.0, -1.6);
		break;
	case ch_type::EPA5Hz:		// extended multipath channel model
		this->Channel.PDP_dB << 0.0 << -1.0 << -2.0 << -3.0 << -8.0 << -17.2 << -20.8 << endr
			<< 0.0 << 30.0 * pow(10.0, -9.0) << 70.0 * pow(10.0, -9.0) << 90.0 * pow(10.0, -9.0) << 110.0 * pow(10.0, -9.0) << 190.0 * pow(10.0, -9.0) << 410.0 * pow(10.0, -9.0) << endr;
		this->Cov_Matrix_time(0, 0) = 1.0 + pow(10.0, -1.0) + pow(10.0, -0.2) + pow(10.0, -0.3) + pow(10.0, -0.8) + pow(10.0, -1.72);
		this->Cov_Matrix_time(1, 1) = pow(10.0, -2.08);
		break;
	case ch_type::EVA5Hz:		// extended multipath channel model
		this->Channel.PDP_dB << 0.0 << -1.5 << -1.4 << -3.6 << -0.6 << -9.1 << -7.0 << -12.0 << -16.9 << endr
			<< 0.0 << 30.0 * pow(10.0, -9.0) << 150.0 * pow(10.0, -9.0) << 310.0 * pow(10.0, -9.0) << 370.0 * pow(10.0, -9.0) << 710.0 * pow(10.0, -9.0) << 1090.0 * pow(10.0, -9.0) << 1730.0 * pow(10.0, -9.0) << 2510.0 * pow(10.0, -9.0) << endr;
		this->Cov_Matrix_time(0, 0) = 1.0 + pow(10.0, -0.15) + pow(10.0, -0.14);
		this->Cov_Matrix_time(1, 1) = pow(10.0, -0.36) + pow(10.0, -0.06) + pow(10.0, -0.91);
		this->Cov_Matrix_time(2, 2) = pow(10.0, -0.7);
		this->Cov_Matrix_time(3, 3) = pow(10.0, -1.2);
		this->Cov_Matrix_time(5, 5) = pow(10.0, -1.69);
		break;
	case ch_type::EVA70Hz:		// extended multipath channel model
		this->Channel.PDP_dB << 0.0 << -1.5 << -1.4 << -3.6 << -0.6 << -9.1 << -7.0 << -12.0 << -16.9 << endr
			<< 0.0 << 30.0 * pow(10.0, -9.0) << 150.0 * pow(10.0, -9.0) << 310.0 * pow(10.0, -9.0) << 370.0 * pow(10.0, -9.0) << 710.0 * pow(10.0, -9.0) << 1090.0 * pow(10.0, -9.0) << 1730.0 * pow(10.0, -9.0) << 20510.0 * pow(10.0, -9.0) << endr;
		this->Cov_Matrix_time(0, 0) = 2.0 + 3.0 * pow(10.0, -0.1);
		this->Cov_Matrix_time(1, 1) = 1.0;
		this->Cov_Matrix_time(3, 3) = pow(10.0, -0.3);
		this->Cov_Matrix_time(4, 4) = pow(10.0, -0.5);
		this->Cov_Matrix_time(10, 10) = pow(10.0, -0.7);
		break;
	case ch_type::ETU70Hz:		// extended multipath channel model
		this->Channel.PDP_dB << -1.0 << -1.0 << -1.0 << 0.0 << 0.0 << 0.0 << -3.0 << -5.0 << -7.0 << endr
			<< 0.0 << 50.0 * pow(10.0, -9.0) << 120.0 * pow(10.0, -9.0) << 200.0 * pow(10.0, -9.0) << 230.0 * pow(10.0, -9.0) << 500.0 * pow(10.0, -9.0) << 1600.0 * pow(10.0, -9.0) << 2300.0 * pow(10.0, -9.0) << 5000.0 * pow(10.0, -9.0) << endr;
		break;
	case ch_type::ETU300Hz:
		this->Channel.PDP_dB << -1.0 << -1.0 << -1.0 << 0.0 << 0.0 << 0.0 << -3.0 << -5.0 << -7.0 << endr
			<< 0.0 << 50.0 * pow(10.0, -9.0) << 120.0 * pow(10.0, -9.0) << 200.0 * pow(10.0, -9.0) << 230.0 * pow(10.0, -9.0) << 500.0 * pow(10.0, -9.0) << 1600.0 * pow(10.0, -9.0) << 2300.0 * pow(10.0, -9.0) << 5000.0 * pow(10.0, -9.0) << endr;
		break;
	case ch_type::AWGN:
	case ch_type::winner_II:
		this->Cov_Matrix_time(0, 0) = 1.0;
		break;
	case ch_type::CDL_A:
		this->CDL
			<< 0.0000 << -13.4 << -178.1 <<   51.3 <<  50.2 << 125.4 << endr
			<< 0.3819 <<   0.0 <<   -4.2 << -152.7 <<  93.2 <<  91.3 << endr
			<< 0.4025 <<  -2.2 <<   -4.2 << -152.7 <<  93.2 <<  91.3 << endr
			<< 0.5868 <<  -4.0 <<   -4.2 << -152.7 <<  93.2 <<  91.3 << endr
			<< 0.4610 <<  -6.0 <<   90.2 <<   76.6 << 122.0 <<  94.0 << endr
			<< 0.5375 <<  -8.2 <<   90.2 <<   76.6 << 122.0 <<  94.0 << endr
			<< 0.6708 <<  -9.9 <<   90.2 <<   76.6 << 122.0 <<  94.0 << endr
			<< 0.5750 << -10.5 <<  121.5 <<   -1.8 << 150.2 <<  47.1 << endr
			<< 0.7618 <<  -7.5 <<  -81.7 <<  -41.9 <<  55.2 <<  56.0 << endr
			<< 1.5375 << -15.9 <<  158.4 <<   94.2 <<  26.4 <<  30.1 << endr
			<< 1.8978 <<  -6.6 <<  -83.0 <<   51.9 << 126.4 <<  58.8 << endr
			<< 2.2242 << -16.7 <<  134.8 << -115.9 << 171.6 <<  26.0 << endr
			<< 2.1718 << -12.4 << -153.0 <<   26.6 << 151.4 <<  49.2 << endr
			<< 2.4942 << -15.2 << -172.0 <<   76.6 << 157.2 << 143.1 << endr
			<< 2.5119 << -10.8 << -129.9 <<   -7.0 <<  47.2 << 117.4 << endr
			<< 3.0582 << -11.3 << -136.0 <<  -23.0 <<  40.4 << 122.7 << endr
			<< 4.0810 << -12.7 <<  165.4 <<  -47.2 <<  43.3 << 123.2 << endr
			<< 4.4579 << -16.2 <<  148.4 <<  110.4 << 161.8 <<  32.6 << endr
			<< 4.5695 << -18.3 <<  132.7 <<  144.5 <<  10.8 <<  27.2 << endr
			<< 4.7966 << -18.9 << -118.6 <<  155.3 <<  16.7 <<  15.2 << endr
			<< 5.0066 << -16.6 << -154.1 <<  102.0 << 171.7 << 146.0 << endr
			<< 5.3043 << -19.9 <<  126.5 << -151.8 <<  22.7 << 150.7 << endr
			<< 9.6586 << -29.7 <<  -56.2 <<   55.2 << 144.9 << 156.1 << endr;

		this->c_ASD =   5.0;
		this->c_ASA =  11.0;
		this->c_ZSD =   3.0;
		this->c_ZSA =   3.0;
		this->XPR_dB = 10.0;

		this->Channel.PDP_dB = zeros(arma::size(CDL.t()));
		this->Channel.PDP_dB = join_horiz(this->CDL.col(1), this->CDL.col(0) * this->DS).t();

		this->Cov_Matrix_time(0,  0)  = pow(10.0, this->Channel.PDP_dB(0,  0) / 10.0);
		this->Cov_Matrix_time(1,  1)  = pow(10.0, this->Channel.PDP_dB(0,  1) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 2) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 4) / 10.0);
		this->Cov_Matrix_time(2,  2)  = pow(10.0, this->Channel.PDP_dB(0,  3) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 5) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 6) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 7) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 8) / 10.0);
		this->Cov_Matrix_time(5,  5)  = pow(10.0, this->Channel.PDP_dB(0,  9) / 10.0);
		this->Cov_Matrix_time(6,  6)  = pow(10.0, this->Channel.PDP_dB(0, 10) / 10.0);
		this->Cov_Matrix_time(7,  7)  = pow(10.0, this->Channel.PDP_dB(0, 11) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 12) / 10.0);
		this->Cov_Matrix_time(8,  8)  = pow(10.0, this->Channel.PDP_dB(0, 13) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 14) / 10.0);
		this->Cov_Matrix_time(9,  9)  = pow(10.0, this->Channel.PDP_dB(0, 15) / 10.0);
		this->Cov_Matrix_time(13, 13) = pow(10.0, this->Channel.PDP_dB(0, 16) / 10.0);
		this->Cov_Matrix_time(14, 14) = pow(10.0, this->Channel.PDP_dB(0, 17) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 18) / 10.0);
		this->Cov_Matrix_time(15, 15) = pow(10.0, this->Channel.PDP_dB(0, 19) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 20) / 10.0);
		this->Cov_Matrix_time(16, 16) = pow(10.0, this->Channel.PDP_dB(0, 21) / 10.0);
		this->Cov_Matrix_time(30, 30) = pow(10.0, this->Channel.PDP_dB(0, 22) / 10.0);
		break;
	case ch_type::CDL_B:
		this->CDL
			<< 0.0000 <<   0.0 <<    9.3 << -173.3 << 105.8 << 78.9 << endr
			<< 0.1072 <<  -2.2 <<    9.3 << -173.3 << 105.8 << 78.9 << endr
			<< 0.2155 <<  -4.0 <<    9.3 << -173.3 << 105.8 << 78.9 << endr
			<< 0.2095 <<  -3.2 <<  -34.1 <<  125.5 << 115.3 << 63.3 << endr
			<< 0.2870 <<  -9.8 <<  -65.4 <<  -88.0 << 119.3 << 59.9 << endr
			<< 0.2986 <<  -1.2 <<  -11.4 <<  155.1 << 103.2 << 67.5 << endr
			<< 0.3752 <<  -3.4 <<  -11.4 <<  155.1 << 103.2 << 67.5 << endr
			<< 0.5055 <<  -5.2 <<  -11.4 <<  155.1 << 103.2 << 67.5 << endr
			<< 0.3681 <<  -7.6 <<  -67.2 <<  -89.8 << 118.2 << 82.6 << endr
			<< 0.3697 <<  -3.0 <<   52.5 <<  132.1 << 102.0 << 66.3 << endr
			<< 0.5700 <<  -8.9 <<  -72.0 <<  -83.6 << 100.4 << 61.6 << endr
			<< 0.5283 <<  -9.0 <<   74.3 <<   95.3 <<  98.3 << 58.0 << endr
			<< 1.1021 <<  -4.8 <<  -52.2 <<  103.7 << 103.4 << 78.2 << endr
			<< 1.2756 <<  -5.7 <<  -50.5 <<  -87.8 << 102.5 << 82.0 << endr
			<< 1.5474 <<  -7.5 <<   61.4 <<  -92.5 << 101.4 << 62.4 << endr
			<< 1.7842 <<  -1.9 <<   30.6 << -139.1 << 103.0 << 78.0 << endr
			<< 2.0169 <<  -7.6 <<  -72.5 <<  -90.6 << 100.0 << 60.9 << endr
			<< 2.8294 << -12.2 <<  -90.6 <<   58.6 << 115.2 << 82.9 << endr
			<< 3.0219 <<  -9.8 <<  -77.6 <<  -79.0 << 100.5 << 60.8 << endr
			<< 3.6187 << -11.4 <<  -82.6 <<   65.8 << 119.6 << 57.3 << endr
			<< 4.1067 << -14.9 << -103.6 <<   52.7 << 118.7 << 59.9 << endr
			<< 4.2790 <<  -9.2 <<   75.6 <<   88.7 << 117.8 << 60.1 << endr
			<< 4.7834 << -11.3 <<  -77.6 <<  -60.4 << 115.7 << 62.3 << endr;

		this->c_ASD = 10.0;
		this->c_ASA = 22.0;
		this->c_ZSD =  3.0;
		this->c_ZSA =  8.0;
		this->XPR_dB = 8.0;

		this->Channel.PDP_dB = zeros(arma::size(CDL.t()));
		this->Channel.PDP_dB = join_horiz(this->CDL.col(1), this->CDL.col(0) * this->DS).t();

		this->Cov_Matrix_time(0, 0) = pow(10.0, this->Channel.PDP_dB(0, 0) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 1) / 10.0);
		this->Cov_Matrix_time(1, 1) = pow(10.0, this->Channel.PDP_dB(0, 2) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 3) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 4) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 5) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 6) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 8) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 9) / 10.0);
		this->Cov_Matrix_time(2, 2) = pow(10.0, this->Channel.PDP_dB(0, 7) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 10) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 11) / 10.0);
		this->Cov_Matrix_time(3, 3) = pow(10.0, this->Channel.PDP_dB(0, 12) / 10.0);
		this->Cov_Matrix_time(4, 4) = pow(10.0, this->Channel.PDP_dB(0, 13) / 10.0);
		this->Cov_Matrix_time(5, 5) = pow(10.0, this->Channel.PDP_dB(0, 14) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 15) / 10.0);
		this->Cov_Matrix_time(6, 6) = pow(10.0, this->Channel.PDP_dB(0, 16) / 10.0);
		this->Cov_Matrix_time(9, 9) = pow(10.0, this->Channel.PDP_dB(0, 17) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 18) / 10.0);
		this->Cov_Matrix_time(11, 11) = pow(10.0, this->Channel.PDP_dB(0, 19) / 10.0);
		this->Cov_Matrix_time(13, 13) = pow(10.0, this->Channel.PDP_dB(0, 20) / 10.0) + pow(10.0, this->Channel.PDP_dB(0, 21) / 10.0);
		this->Cov_Matrix_time(15, 15) = pow(10.0, this->Channel.PDP_dB(0, 12) / 10.0);

		break;
	case ch_type::CDL_C:
		this->CDL
			<< 0.0000 <<  -4.4 <<  -46.6 << -101.0 <<  97.2 <<  87.6 << endr
			<< 0.2099 <<  -1.2 <<  -22.8 <<  120.0 <<  98.6 <<  72.1 << endr
			<< 0.2219 <<  -3.5 <<  -22.8 <<  120.0 <<  98.6 <<  72.1 << endr
			<< 0.2329 <<  -5.2 <<  -22.8 <<  120.0 <<  98.6 <<  72.1 << endr
			<< 0.2176 <<  -2.5 <<  -40.7 << -127.5 << 100.6 <<  70.1 << endr
			<< 0.6366 <<   0.0 <<    0.3 <<  170.4 <<  99.2 <<  75.3 << endr
			<< 0.6448 <<  -2.2 <<    0.3 <<  170.4 <<  99.2 <<  75.3 << endr
			<< 0.6560 <<  -3.9 <<    0.3 <<  170.4 <<  99.2 <<  75.3 << endr
			<< 0.6584 <<  -7.4 <<   73.1 <<   55.4 << 105.2 <<  67.4 << endr
			<< 0.7935 <<  -7.1 <<  -64.5 <<   66.5 <<  95.3 <<  63.8 << endr
			<< 0.8213 << -10.7 <<   80.2 <<  -48.1 << 106.1 <<  71.4 << endr
			<< 0.9336 << -11.1 <<  -97.1 <<   46.9 <<  93.5 <<  60.5 << endr
			<< 1.2285 <<  -5.1 <<  -55.3 <<   68.1 << 103.7 <<  90.6 << endr
			<< 1.3083 <<  -6.8 <<  -64.3 <<  -68.7 << 104.2 <<  60.1 << endr
			<< 2.1704 <<  -8.7 <<  -78.5 <<   81.5 <<  93.0 <<  61.0 << endr
			<< 2.7105 << -13.2 <<  102.7 <<   30.7 << 104.2 << 100.7 << endr
			<< 4.2589 << -13.9 <<   99.2 <<  -16.4 <<  94.9 <<  62.3 << endr
			<< 4.6003 << -13.9 <<   88.8 <<    3.8 <<  93.1 <<  66.7 << endr
			<< 5.4902 << -15.8 << -101.9 <<  -13.7 <<  92.2 <<  52.9 << endr
			<< 5.6077 << -17.1 <<   92.2 <<    9.7 << 106.7 <<  61.8 << endr
			<< 6.3065 << -16.0 <<   93.3 <<    5.6 <<  93.0 <<  51.9 << endr
			<< 6.6374 << -15.7 <<  106.6 <<    0.7 <<  92.9 <<  61.7 << endr
			<< 7.0427 << -21.6 <<  119.5 <<  -21.9 << 105.2 <<  58.0 << endr
			<< 8.6523 << -22.8 << -123.8 <<   33.6 << 107.8 <<  57.0 << endr;

		this->Channel.PDP_dB = zeros(arma::size(CDL.t()));
		this->Channel.PDP_dB = join_horiz(CDL.col(1), (CDL.col(0) * this->DS)).t();

		this->c_ASD =  2.0;
		this->c_ASA = 15.0;
		this->c_ZSD =  3.0;
		this->c_ZSA =  7.0;
		this->XPR_dB = 7.0;

		/* calculate cov_matrix_time */
		PDP_Delay = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		PDP_Delay = this->Channel.PDP_dB.row(1);
		New_PDP_Delay_samples = zeros<rowvec>(this->Channel.PDP_dB.n_cols);
		New_PDP_Delay_samples = round(PDP_Delay * this->Fs);
		for (size_t ind = 0; ind < New_PDP_Delay_samples.size(); ind++) {
			diag_index = New_PDP_Delay_samples(ind);
			this->Cov_Matrix_time(diag_index, diag_index) = this->Cov_Matrix_time(diag_index, diag_index) + pow(10.0, this->Channel.PDP_dB(0, ind) / 10.0);
		}

		break;
	case ch_type::CDL_D:
		this->CDL
			<< 0.0    <<  -0.2 <<    0.0 << -180.0 <<  98.5 << 81.5 << endr
			<< 0.0    << -13.5 <<    0.0 << -180.0 <<  98.5 << 81.5 << endr
			<< 0.035  << -18.8 <<   89.2 <<   89.2 <<  85.5 << 86.9 << endr
			<< 0.612  << -21.0 <<   89.2 <<   89.2 <<  85.5 << 86.9 << endr
			<< 1.363  << -22.8 <<   89.2 <<   89.2 <<  85.5 << 86.9 << endr
			<< 1.405  << -17.9 <<   13.0 <<  163.0 <<  97.5 << 79.4 << endr
			<< 1.804  << -20.1 <<   13.0 <<  163.0 <<  97.5 << 79.4 << endr
			<< 2.596  << -21.9 <<   13.0 <<  163.0 <<  97.5 << 79.4 << endr
			<< 1.775  << -22.9 <<   34.6 << -137.0 <<  98.5 << 78.2 << endr
			<< 4.042  << -27.8 <<  -64.5 <<   74.5 <<  88.4 << 73.6 << endr
			<< 7.937  << -23.6 <<  -32.9 <<  127.7 <<  91.3 << 78.3 << endr
			<< 9.424  << -24.8 <<   52.6 << -119.6 << 103.8 << 87.0 << endr
			<< 9.708  << -30.0 << -132.1 <<   -9.1 <<  80.3 << 70.6 << endr
			<< 12.525 << -27.7 <<   77.2 <<  -83.8 <<  86.5 << 72.9 << endr;

		this->Channel.PDP_dB = zeros(arma::size(CDL.t()));
		this->Channel.PDP_dB = join_horiz(CDL.col(1), (CDL.col(0) * this->DS)).t();

		this->c_ASD  =  5.0;
		this->c_ASA  =  8.0;
		this->c_ZSD  =  3.0;
		this->c_ZSA  =  3.0;
		this->XPR_dB = 11.0;

		break;
	case ch_type::CDL_E:
		this->CDL
			<<  0.000  <<  -0.03 <<   0.0 << -180.0 <<  99.6 << 80.4 << endr
			<<  0.000  << -22.03 <<   0.0 << -180.0 <<  99.6 << 80.4 << endr
			<<  0.5133 << -15.8  <<  57.5 <<   18.2 << 104.2 << 80.4 << endr
			<<  0.5440 << -18.1  <<  57.5 <<   18.2 << 104.2 << 80.4 << endr
			<<  0.5630 << -19.8  <<  57.5 <<   18.2 << 104.2 << 80.4 << endr
			<<  0.5440 << -22.9  << -20.1 <<  101.8 <<  99.4 << 80.8 << endr
			<<  0.7112 << -22.4  <<  16.2 <<  112.9 << 100.8 << 86.3 << endr
			<<  1.9092 << -18.6  <<   9.3 << -155.5 <<  98.8 << 82.7 << endr
			<<  1.9293 << -20.8  <<   9.3 << -155.5 <<  98.8 << 82.7 << endr
			<<  1.9589 << -22.6  <<   9.3 << -155.5 <<  98.8 << 82.7 << endr
			<<  2.6426 << -22.3  <<  19.0 << -143.3 << 100.8 << 82.9 << endr
			<<  3.7136 << -25.6  <<  32.7 <<  -94.7 <<  96.4 << 88.0 << endr
			<<  5.4524 << -20.2  <<   0.5 <<  147.0 <<  98.9 << 81.0 << endr
			<< 12.0034 << -29.8  <<  55.9 <<  -36.2 <<  95.6 << 88.6 << endr
			<< 20.6419 << -29.2  <<  57.6 <<  -26.0 << 104.6 << 78.3 << endr;

		this->Channel.PDP_dB = zeros(arma::size(CDL.t()));
		this->Channel.PDP_dB = join_horiz(CDL.col(1), (CDL.col(0) * this->DS)).t();

		this->c_ASD  =  5.0;
		this->c_ASA  = 11.0;
		this->c_ZSD  =  3.0;
		this->c_ZSA  =  7.0;
		this->XPR_dB =  8.0;

		break;
	default:
		cout << "ERROR: Channel not supported" << endl;
	}   /* end of switch (this->Channel.type) */

	this->Cov_Matrix_time = this->Cov_Matrix_time / sum(this->Cov_Matrix_time.diag());
	cx_mat fft_matrix = fft(eye(this->size_fft, this->size_fft));

	int front_start_val = this->size_fft - this->num_RB * this->num_sc / 2 + 1;
	int front_end_val = this->size_fft;
	vector<int> freq_permute;
	for (int i = 0; i <= front_end_val - front_start_val; i++)
		freq_permute.push_back(front_start_val + i);
	int back_start_val = 2;
	int back_end_val = this->num_RB * this->num_sc / 2 + 1;
	for (int i = 0; i <= back_end_val - back_start_val; i++)
		freq_permute.push_back(back_start_val + i);

	this->Cov_Matrix_freq = zeros<cx_mat>(this->size_fft, this->size_fft);

	cx_mat Cov_Matrix_time_Comp = zeros<cx_mat>(size_fft, size_fft);
	cx_mat Cov_Matrix_time_imag = zeros<cx_mat>(size_fft, size_fft);
	for (int i = 0; i < size_fft; i++) {
		for (int j = 0; j < size_fft; j++) {
			Cov_Matrix_time_imag(i, j) = cx_double(0.0, 1.0);
			Cov_Matrix_time_Comp(i, j) = cx_double(Cov_Matrix_time(i, j), 1.0);
		}
	}
	this->Cov_Matrix_freq = ((fft_matrix * Cov_Matrix_time_Comp) - (fft_matrix * Cov_Matrix_time_imag)) * fft_matrix.t();

	cx_mat temp_Cov_Matrix_freq(freq_permute.size(), freq_permute.size(), fill::zeros);
	for (int i = 0; i < freq_permute.size(); i++)
		for (int j = 0; j < freq_permute.size(); j++)
			temp_Cov_Matrix_freq(i, j) = this->Cov_Matrix_freq(freq_permute[i] - 1, freq_permute[j] - 1);
	this->Cov_Matrix_freq = temp_Cov_Matrix_freq;

	this->user_speed = init.user_speed;
	this->HARQ_switch = init.HARQ_switch;

	/* get CQI Table */
	getCQITable();

	/* initialize para.temp_ch_coding() field */
	this->temp_ch_coding = field<ChCoding>(10 * pow(2, this->mu));

	/* memory deallocations */
	freq_permute.clear();			// vector container deallocation
	Cov_Matrix_time_Comp.reset();
	Cov_Matrix_time_Comp.reset();
	temp_Cov_Matrix_freq.reset();

}	/* end of module_Parameter_MIMO() method definition */





/*
 * ModuleParameterMIMO::method_calculation() - class method.
 */
void ModuleParameterMIMO::method_calculation(void)
{
	this->time_CP = zeros<rowvec>(2);
	this->length_CP = zeros<irowvec>(2);

	if (!this->num_RB_maximum_isDirty)
		cout << "ERROR: Input bandwidth is not allowed. The allowed values are " << this->BW_set / 1e6 << "MHz" << endl;

	this->Fs = this->size_fft * this->subcarrier_spacing;
	this->kappa = 64;
	this->Tc = this->Ts / this->kappa;
	this->time_CP(0) = (144.0 * this->kappa * pow(2.0, -this->mu) + 16.0 * this->kappa) * this->Tc;
	this->time_CP(1) = (144.0 * this->kappa * pow(2.0, -this->mu)) * this->Tc;
	this->length_CP(0) = this->time_CP[0] * this->Fs;
	this->length_CP(1) = this->time_CP[1] * this->Fs;

	this->num_port = this->port_set.n_elem;
	this->w_d = 2.0 * datum::pi * this->user_speed * this->f / this->light_speed;

}   /* end of method_calculation(void) method implementation */




/*
 * ModuleParameter::Reset - class method
 *	- input: double SNR, int ind_SNR
 *	- output: direct access to the class member variables
 */
void ModuleParameterMIMO::Reset(double SNR, int ind_SNR)
{
	this->SNR = SNR;
	this->ind_SNR = ind_SNR;
	this->sigma_n_freq = 1.0 / pow(10.0, this->SNR / 10.0);
	this->sigma_n_time = this->size_fft / (this->num_RB * this->num_sc) * 1.0 / pow(10.0, this->SNR / 10.0);

}   /* end of Reset(double SNR, int ind_SNR) method implementation */



/* end of module_Parameter_MIMO.cpp */