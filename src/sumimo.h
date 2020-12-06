/*
 * sumimo.h - init parameters and global variable definition.
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */
#ifndef __SUMIMO_H__
#define __SUMIMO_H__

#include <cmath>
#include <armadillo>


#if defined(_MSC_BUILD)
#define	sprintf		sprintf_s
#define strcpy		strcpy_s
#define stricmp		_stricmp
#endif

#define CQI_LENGTH			 1
#define SNR_LENGTH		    14
#define SNR_STEPSIZE_NUM	 4

#define LDPC_CPP			 0
#define LDPC_MATLAB			 1

using namespace std;
using namespace arma;

const double light_speed = 299792458;



enum class Sim_Case {
	SUSISO,
	SUMIMO,
	SUMIMO_CLSM,
	SUMIMO_Noncodebook,
	MUMIMO,
	MUMIMO_CLSM,
	MUMIMO_Noncodebook,
	UNSUPPORTED
};


enum class Ch_type {
	PedA,
	flat_Rayleigh,
	Rayleigh2,
	PedB,
	PedBcorr,
	VehA,
	VehB,
	EPA5Hz,
	EVA5Hz,
	EVA70Hz,
	ETU70Hz,
	ETU300Hz,
	AWGN,
	TDL_A,
	TDL_B,
	TDL_C,
	TDL_D,
	TDL_E,
	CDL_A,
	CDL_B,
	CDL_C,
	CDL_D,
	CDL_E,
	winner_II,
	UNSUPPORTED
};



enum class Ch_mode {
	Block_Fading,
	Fast_Fading,
	flat_Rayleigh,
	UNSUPPORTED
};


enum class Ch_est_mode {
	Perfect,
	LS,
	MMSE,
	UNSUPPORTED
};

enum class Ch_est_mode_DMRS {
	Perfect,
	LS,
	MMSE,
	UNSUPPORTED
};

enum class Ch_interp_mode {
	linear,
	nearest,
	spline,
	TDI,
	UNSUPPORTED
};

enum class Ch_interp_mode_DMRS {
	linear,
	UNSUPPORTED
};

enum class Precoder_case {
	ZF,
	UNSUPPORTED
};

enum class Time_correlation {
	independent,
	correlated,
	UNSUPPORTED
};


enum class Method_interpolation {
	nearest_neighbor,
	sinc_interpolation,
	UNSUPPORTED
};

enum class Sync_freq_mode {
	Perfect,
	estimated,
	UNSUPPORTED
};

enum class Zero_TTI_feedback_case {
	Perfect,
	UNSUPPORTED
};

enum class Ue_detection_mode {
	ZF,
	MMSE,
	Sphere_decoding,
	UNSUPPORTED
};

enum class Codebook_mode {
	Random,
	DFT,
	UNSUPPORTED
};

enum class Pol_tx {
	cross_pol,
	UNSUPPORTED
};

enum class Pol_rx {
	cross_pol,
	UNSUPPORTED
};

enum class HARQ_switch {
	on,
	off,
	UNSUPPORTED
};

enum class PDSCH_mapping_type {
	A,
	B,
	UNSUPPORTED
};



struct feedback {
	irowvec			CQI;
	irowvec			ACK;
	int				NACK_a;
	urowvec			NACK_b;
	int				PMI;
	int				RI;
};


struct ChCoding {
	int					A;
	int					B;
	int					B_p;
	int					C;
	imat				DD;
	irowvec				E;
	int					K;
	int					K_p;
	int					L;
	int					L_seg;
	int					M_H;
	int					N;
	int					N_cb;
	int					N_H;
	int					Q;
	int					Z_c;
	int					CRC_Type;
	int					CRC_seg_Type;
	int					BG_index;
	int					null_cnt;
	urowvec				null_position;
	irowvec				col_weight;
	irowvec				row_weight;
	int					max_col;
	int					max_row;
	irowvec				INDEX_col;
	irowvec				INDEX_row;
	imat				invBB_AA;
	int					rv_id;
	int					k0;
	int					Decoding_type;
	int					Max_LDPC_itr;
	int					MS_factor;
	int					LLR_Threshold;
	int					NACK_a;
	urowvec				NACK_b;
	urowvec				tx_a;
	urowvec				tx_g;
	urowvec				rx_a;
	urowvec				rx_g;
	rowvec				symbol_LLR;
	int					ACK_for_check;
	
	struct NR_CODEC {
		int		data_length;
		double	code_rate;
		int		G;
		int		ACK;
		int		ACK_cnt;
		int		NACK;
		int		NACK_cnt;
		struct LDPC {
			int		LDPCMaxIter;
			int		LDPCDecOption;
			double	LDPC_MinSum_factor;
			int		LBRM_indicator;
			int		TBS_LBRM;
		} ldpc;
		struct HARQ {
			int		HARQ_indicator;
			int		HARQ_type;
			int		HARQ_iter;
			int		HARQ_max_iter;
			int	   *HARQ_index;
		} harq;

	} nr_codec;

};   /* end of ChCoding struct definition */



struct angleFirst {
	double		ZOD;
	double		AOD;
};


struct ChannelOutput {
	cx_cube				time;
	field<cx_mat>		time4D;
	cx_colvec			fft1D;
	cx_mat				fft2D;
	field<cx_cube>		fft;
	angleFirst			angle_first;
};   /* end of struct output */



struct Genie {
	int					num_codewords;
	int					PMI;
	int					RI;
	int					CQI_table_num;
	irowvec				CQI;
	irowvec				m_order;
	
	/* Related to position */
	umat				PDCCH_position;
	umat				PSS_position;
	umat				SSS_position;
	umat				PCFICH_position;
	umat				PBCH_position;
	ucube				DMRS_position;
	ucube				CSIRS_position;
	umat				CSIRS_position_total;
	umat				Data_position;
	ucube				CRS_position;

	/* Related to signal */
	cx_mat				DMRS_signal;
	cx_mat				CSIRS_signal;
	cx_mat				PSS_position_signal;

	cx_cube				resource_grid;
	cx_cube				resource_grid_DMRS;
	cx_mat				resource_grid_CSIRS;
	
	/* HARQ and channel coding */
	uvec				HARQ_process_index;
	uvec				newdata_indicator;
	field<ChCoding>		ch_coding;
	
	/* power allocation */
	rowvec				power_data;
	rowvec				power_DMRS;
	rowvec				power_CSIRS;

	/* bit stream */
	field<urowvec>		bit_stream;
	field<urowvec>		coded_bit_stream;
	field<urowvec>		scrambled_bit_stream;
	cx_mat				OFDM_signal_stream;

	/* MIMO */
	umat				num_Data_RB;
	int					num_layers;
	field<cx_mat>		Precoder;

	bool				selected;

};   /* end of Genie structure definition */



enum class LDPC_MODULE {
	MATLAB,
	CPLUS,
	UNSUPPORTED
};


enum class RANDOM_GEN {
	RANDOM_FIX,
	RANDOM_ENABLE,
	UNSUPPORTED
};





class ParameterInit {
public:
	Sim_Case		sim_Case;			// SUSISO, SUMIMO
	unsigned short	mu;					// numerology: [0 1 2 3] 
	double			BW;					// Allowed values: [ 5 10 15 20 25 30 40 50 60 70 80 90 100 200 400 ] MHz
	unsigned short	num_subframe;		// The number of subframe
	Ch_type			ch_type;			// Channel model type: 'AWGN', 'CDL_A', 'CDL_B', 'CDL_C', 'CDL_D',  'CDL_E', ..., 'PedA', 'PedB', 'VehA', 'flat Rayleigh' ...
	Ch_mode			ch_mode;			// 'Block_Fading' or 'Fast_Fading' % 'flat Rayleigh' - Block_Fading only / 'Rayleigh2' - Fast_fading only
	Ch_est_mode		ch_est_mode;
	Ch_est_mode_DMRS	ch_est_mode_DMRS;
	double			f;					// Carrier frequencys
	Pol_tx			pol_tx;
	Pol_rx			pol_rx;

	unsigned short	num_Tx_antenna;
	bool			num_Tx_antenna_isExist = true;
	bool			num_Tx_antenna_isDirty = false;
	unsigned short	num_Rx_antenna;
	bool			num_Rx_antenna_isExist = true;
	bool			num_Rx_antenna_isDirty = false;
	int				num_Tx_antenna_horizontal;
	bool			num_Tx_antenna_horizontal_isExist = true;
	bool			num_Tx_antenna_horizontal_isDirty = false;
	int				num_Tx_antenna_vertical;
	bool			num_Tx_antenna_vertical_isExist = true;
	bool			num_Tx_antenna_vertical_isDirty = false;
	int				num_Rx_antenna_horizontal;
	int				num_Rx_antenna_vertical;

	double			user_speed;
	HARQ_switch		harq_switch;
	mat				SNR_range;
	bool			SNR_range_isExist = true;
	bool			SNR_range_isDirty = false;
	int				num_RB;
	bool			num_RB_isExist = true;
	bool			num_RB_isDirty = false;

	int				Tx_pol;
	int				Rx_pol;

	bool			cal_scenario = false;
	bool			cal_scenario_isExist = true;
	bool			cal_scenario_isDirty = false;

	irowvec			CQI;

	int				CSIRS_row;
	int				CSIRS_density;
	int				CSIRS_periodicity;
	int				DMRS_config;
	PDSCH_mapping_type	PDSCH_type;


};	/* end of ParameterInit class definition */



struct execTime {	// Release verision ������ ���� ����
	clock_t				start;
	clock_t				end;
	clock_t				initTime;
	clock_t				cqiTime;
	clock_t				snrTime;
	clock_t				slotTime;
	clock_t				chTime;
	clock_t				BSTime;
	clock_t				UETime;
	clock_t				exeTime;
	unsigned short		ms;
	unsigned short		sec;
	unsigned short		min;
	unsigned short		hr;
};



/*
 * global function definition
 */

void func_pseudo_sequence_generation(int, double, mat &);

int		mod(int, int);
int		imod(int, int);
double	dmod(double, double);
umat	matlnot(umat &);
double	fix(double);



#endif		/* end of sumimo.h */