/*
 * module_Paramter_MIMO.h - 
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */
#ifndef __MODULE_PARAMETER_MIMO_H__
#define __MODULE_PARAMETER_MIMO_H__

#include <vector>
#include <armadillo>
#include "sumimo.h"


using namespace std;
using namespace arma;



enum class CSIRS_CDMType {
	NO_CDM,
	FD_CDM2,
	CDM4,
	CDM8,
	UNSUPPORTED
};


enum class DL_DMRS_config_type {
	type_1,
	type_2,
	UNSUPPORTED
};


enum class TxArrayType {
	ULA,
	URA,
	UNSUPPORTED
};


enum class RxArrayType {
	ULA,
	URA,
	UNSUPPORTED
};


enum class Tx_pattern_type {
	Pattern,
	Omni_directional,
	UNSUPPORTED
};


enum class Rx_pattern_type {
	Pattern,
	Omni_directional,
	UNSUPPORTED
};



struct Channel {
	ch_type				type;
	ch_mode				time_fading;
	time_correlation	time_correlation;
	mat					PDP_dB;
};



/*
 * ModuleParameterMIMO class definition
 */
class ModuleParameterMIMO {
public:
	/* LDPC init */
	field<ChCoding>		temp_ch_coding;

	/* numerology, setting */
	int					mu = 0;							// numerology, default : 0, value: 0 1 2 3 4
	int					subcarrier_spacing;				// TS 38.211, Table 4.2-1
	int					num_sc = 12;
	int					num_symb = 14;
	int					num_slot_in_subframe;
	double				Ts = 1.0 / (15e3 * 2048.0);		// Reference time units

	int					size_fft = 4096;				// 5G, maximum default
	double				light_speed = 299792458;
	double				Fs;								// Sample frequency
	rowvec				time_CP;						// CP time
	irowvec				length_CP;						// CP length

	/* RS----------------*/
	irowvec				num_port_DMRS;
	int					num_port_CSIRS;
	irowvec				port_set_DMRS;					// = ["1000"];  %port_set -> DMRS_port_set
	irowvec				port_set_CSIRS;					// = ["3000", "3001", "3002", "3003", "3004", "3005", "3006", "3007", "3008", "3009", "3010", "3011", "3012", "3013", "3014", "3015", "3016", "3017", "3018", "3019", "3020", "3021", "3022", "3023", "3024", "3025", "3026", "3027", "3028", "3029", "3030", "3031"];
	
	int					CSIRS_density = 1;
	CSIRS_CDMType		CSIRS_CDMType;					// = 'CDM8' 'no CDM', 'FD-CDM2', 'CDM4', 'CDM8'
	irowvec				CSIRS_bitmat_value;				// Fix this for now related with starting point of CSI-RS in freq. domain
	irowvec				CSIRS_symbol_start_idx;			// Fix this for now related with starting point of CSI-RS in time domain
	int					row_CSIRS;						// 1~ 19, TS 38.211 'Table 7.4.1.5.2-1: CSI-RS locations within a slot.'

	int					num_symb_in_slot = 14;
	int					num_symb_in_subframe;

	PDSCH_mapping_type	PDSCH_mapping_type;
	DL_DMRS_config_type	DL_DMRS_config_type;
	int					DMRS_type = 1;
	int					DL_DMRS_add_pos;
	int					DL_DMRS_typeA_pos;
	int					DL_DMRS_max_len;
	int					pos_last_PDSCH_symbol;
	int					dur_PDSCH_transmission;
	int					CSIRS_periodicity;

	int					kappa;
	double				Tc;

	int					num_SS_block;
	irowvec				slotset_SS_block;

	/* Related to Resource allocation */
	bool				use_PDCCH = false;
	bool				use_PBCH = false;
	bool				use_Power_allocation = false;		// default
	bool				use_Sync = false;

	int					CFI = 1;

	/* power */
	rowvec				power_DMRS;							// = [1.3], for MIMO, [1.1, 1.2, 1.1] layer ....
	rowvec				power_CSIRS;						// = [1.4, 1.4, 1.4, ..],
	rowvec				power_data;							// = [1]. for MIMO [1, 1, 1.1] layer...

	/* Related to Scenario */
	bool				cal_scenario = false;
	int					scenario_beam = 1;
	Sim_Case			Sim_Case;
	double				BW;
	double				f;									// Carrier frequency 2110e6 / 2.5e9

	irowvec				BW_set;
	irowvec				RB_set;
	int					num_RB_maximum;						// # of maximum number of RB for bandwidth - jyseo -
	bool				num_RB_maximum_isDirty = false;		//dirty bit of num_RB_maximum variable
	int					num_subframe;
	int					num_port;
	irowvec				port_set;							// Set of antenna port numbers used for transmission
	int					num_Tx_antenna;
	int					num_Rx_antenna;

	int					N1;									// # of horizontal Tx antennas
	int					N2;									// # of vertical Tx antennas

	int					M1;									// # of horizontal Rx antennas
	int					M2;									// # of vertical Rx antennas

	int					O1;									// # of DFT oversampling of horizontal axis
	int					O2;									// # of DFT oversampling of vertical axis

	int					P1;									// # of horizontal narrow beams per a wide beam
	int					P2;									// # of vertical narrow beams per a wide beam
															// L = L1 * L2, # of wide beams
	int					L1;									// # of horizontal wide beams
	int					L2;									// # of vertical wide beams
															// Actual wide beam is kronecker product of horizontal wide beam and vertical wide beam
	int					S1;									// Narrow beam interval that composed of different wide beam (horizontal)
	int					S2;									// Narrow beam interval that composed of different wide beam (vertical)

	irowvec				pol_slant_angle;					// polarization slant angle

	TxArrayType			TxArrayType = TxArrayType::URA;		// Transmit array type (ULA or URA)
	RxArrayType			RxArrayType = RxArrayType::URA;		// Receive array type (ULA or URA)

	double				Tx_d_lambda = 0.5;					// Tx antenna spacing
	double				Rx_d_lambda = 0.5;					// Rx antenna spacing
	int					Tx_pol;
	int					Rx_pol;
	int					Tx_downtilt = 0;

	Tx_pattern_type		Tx_pattern_type = Tx_pattern_type::Pattern;				// 'Omni-directional' or 'Pattern'
	Rx_pattern_type		Rx_pattern_type = Rx_pattern_type::Omni_directional;	// 'Omni-directional' or 'Pattern'

	irowvec				TM;
	irowvec				fixed_RI;
	int					RI_max = 2;

	mat					SNR_range;
	bool				SNR_range_isExist = true;
	bool				SNR_range_isDirty = false;
	double				SNR;
	int					ind_SNR;
	double				sigma_n_time;
	double				sigma_n_freq;
	int					cell_ID = 0;
	irowvec				num_user;
	bool				MUMIMO = false;

	bool				use_fullband = true;
	int					num_RB;								// # of RB for PDSCH
	bool				num_RB_isExist = true;
	bool				num_RB_isDirty = false;

	/* Related to Channel model */
	double				DS = 100 * 1e-9;					// delay spread
	mat					CDL;
	mat					TDL;
	double				user_speed = 0;						// (double) 5.0 / 6.0;	// [m/s], SUSISO

	Channel				Channel;							// Channel.type, Channel.PDP_dB
	mat					Cov_Matrix_time;
	cx_mat				Cov_Matrix_freq;
	double				Channel_coefRX;
	double				Channel_coefTx;

	method_interpolation	method_interpolation = method_interpolation::nearest_neighbor; // "sinc_interpolation" or "nearest_neighbor"

	int					sin_num = 10;						// Number of sin realizations
	double				w_d;

	double				theta_v;
	double				phi_v;
	rowvec				t;

	int					AS_ratio = 1;						// Channel realization from MIMO 5G channel model, default: 1 (AS_desired/AS_model)

	double				c_ASD;								// degrees
	double				c_ASA;
	double				c_ZSD;
	double				c_ZSA;
	double				XPR_dB;

	/* Relasted to channel estimation */
	ch_est_mode			ch_est_mode = ch_est_mode::MMSE;
	ch_est_mode_DMRS	ch_est_mode_DMRS = ch_est_mode_DMRS::MMSE;
	ch_interp_mode		ch_interp_mode = ch_interp_mode::linear;
	ch_interp_mode_DMRS	ch_interp_mode_DMRS = ch_interp_mode_DMRS::linear;

	/* Related to channel codeing */
	int					MCS_table_index = 1;				// whether use table 1 or 2 TS 38.214, Table 5.1.3.1-1, 5.1.3.1-2
	double				ch_coding;
	HARQ_switch			HARQ_switch;
	int					num_HARQ_retransmission = 0;		// 0, 1, 2, 3
	int					num_HARQ_process = 8;
	bool				hard_decision = false;

	/* Related to transmission */
	int					num_codewords = 1;					// Default : 1
	irowvec				num_data_bits;
	irowvec				num_coded_bits;

	int					M_order = 2;
	double				Coding_rate = (double) 1.0 / 3.0;

	/* Related to frequency syncronization */
	int					sync_freq_offset = 0;						// Normalized frequency offset
	sync_freq_mode		sync_freq_mode = sync_freq_mode::Perfect;	// "perfect", "estimated"

	/* Reated to Precoding */
	int					subband;
	irowvec				subband_set;

	bool					zero_TTI_feedback = true;
	zero_TTI_feedback_case	zero_TTI_feedback_case = zero_TTI_feedback_case::Perfect;

	precoder_case			precoder_case;
	string					codebook;						// codebook file name include file extension

	/* Related to receiver detection */
	ue_detection_mode		ue_detection_mode = ue_detection_mode::MMSE;
	int						num_feed_bit = 10;
	double					CB;
	codebook_mode			codebook_mode = codebook_mode::Random;


	/* class member method definition */
	ModuleParameterMIMO(void);
	~ModuleParameterMIMO(void);
	void module_Parameter_MIMO(const ParameterInit &);
	void method_calculation(void);
	void Reset(double, int);

};	/* end of ModuleParamterMIMO class definition */



#endif		/* end of module_Paramter_MIMO.h */