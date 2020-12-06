/*
 * module_BS_SUMIMO.cpp
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *		Last Updated: 2018-12-05 (started on Octover, 10)
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */

#include <armadillo>
#include "sumimo.h"
#include "module_BS_SUMIMO.h"
#include "getCQITable.h"
#include "LDPC_MATLAB\/LDPC_initialization.h"
#include "LDPC_MATLAB\/LDPC_tx.h"
#include "LDPC_CCv1.h"


#undef max
#undef min


using namespace std;
using namespace arma;


extern bool		fixRandom;
extern char		LDPC_MODULE;
extern int		SLOT;


/* local function declarations */
void func_wide_beam_gen(ModuleParameterMIMO &, field<cx_mat> &, field<cx_mat> &);
void func_beam_selection(ModuleParameterMIMO &, cx_cube &, cx_mat &, int &, field<cx_mat> &);
void func_CDM(CSIRS_CDMType, imat &, imat &, int &, int &);




/*
 * module_BS_SUMIMO() - User constructor of module_BS_SUMIMO class
 */
void ModuleBS_SUMIMO::module_BS_SUMIMO(ModuleParameterMIMO &para)
{

	this->HARQ_buffer.num_process = para.num_HARQ_process;
	this->HARQ_buffer.max_retrans = para.num_HARQ_retransmission;
	this->HARQ_buffer.coded_bits  = field<rowvec>(HARQ_buffer.num_process);
	this->HARQ_buffer.process	  = zeros<uvec>(HARQ_buffer.num_process);
	this->HARQ_buffer.timing	  = zeros<uvec>(HARQ_buffer.num_process);
	this->HARQ_buffer.retrans	  = field<uvec>(HARQ_buffer.num_process);
	this->HARQ_process_index = zeros<uvec>(1);					//  method_data_allocation()에서 상수 1로 설정이 됨.
	this->HARQ_process_index(0)	  = 0;

	this->M_order     = para.M_order;
	this->Coding_rate = para.Coding_rate;

}	/* end of module_BS_SUMIMO() method */




/*
 * Process() - main gateway method of BS_SUMIMO module.
 */
void ModuleBS_SUMIMO::Process(ModuleParameterMIMO &para, NR_ChannelCoding &LDPC,  int &ind_slot, feedback &feedback, ChannelOutput &ch_Output)
{
	//cout << "DEBUG: BS.process enterance." << endl;

	method_feedback_reception(feedback);

	//cout << "DEBUG: BS.method_feedback_reception." << endl;

	if (!para.cal_scenario) {
		method_CSI_RS_allocation(para, ind_slot);
	
		//cout << "DEBUG: BS.method_CSI_RS_allocation." << endl;
		
		method_resource_allocation(para, ind_slot);

		//cout << "DEBUG: BS.method_resource_allocation." << endl;
	}

	method_beam_management(para, ch_Output, ind_slot);

	//cout << "DEBUG: BS.method_beam_management." << endl;

	method_DMRS_allocation(para, ind_slot);

	//cout << "DEBUG: BS.method_DMRS_allocation." << endl;

	method_data_allocation(para, ind_slot);

	//cout << "DEBUG: BS.method_data_allocation." << endl;

	method_sync_generation(para, ind_slot);

	//cout << "DEBUG: BS.method_sync_generation." << endl;

	method_HARQ_process(para, ind_slot, feedback);

	//cout << "DEBUG: BS.method_HARQ_process." << endl;

	method_bit_generation(para, ind_slot);

	//cout << "DEBUG: BS.method_bit_generation." << endl;

	method_channel_coding(para, LDPC, ind_slot, feedback);

	//cout << "DEBUG: BS.method_channel_coding." << endl;

	method_scrambling(para, ind_slot);

	//cout << "DEBUG: BS.method_scrambling." << endl;

	method_modulation_mapping();

	//cout << "DEBUG: BS.method_modulation_mapping." << endl;

	method_layer_mapping(para);

	//cout << "DEBUG: BS.method_layer_mapping." << endl;

	method_precoding(para);

	//cout << "DEBUG: BS.method_precoding." << endl;

	method_OFDM_generation(para);

	//cout << "DEBUG: BS.method_OFDM_generation." << endl;

	method_genie();

	//cout << "DEBUG: BS.method_genie. (BS module END)" << endl;
	
}	/* end of Process() method definition */




/*
 * method_feedback_receotion()
 */
void ModuleBS_SUMIMO::method_feedback_reception(feedback &feedback)
{

	this->CQI = feedback.CQI;

}	/* end of method_feedback_reception() method definition */




/*
 * method_CSI_RS_allocation() method - CSI-RS allocation.
 */
void ModuleBS_SUMIMO::method_CSI_RS_allocation(ModuleParameterMIMO &para, int &ind_slot)
{
	/* local variables definition */
	irowvec			k, l;
	irowvec			freq_pos, time_pos;
	imat			CSIRS_w_f, CSIRS_w_t;
	int				CDM_length, CDM_length_f, CDM_length_t;
	icube			RB_position;
	int				CSIRS_density;
	imat			RB_position_tot;
	field<cx_mat>	pseudo_sequence_CSIRS;
	int				c_init_CSIRS;
	mat				c;
	uvec			index_for_time;
	int				ind_symb;
	cx_colvec		temp;
	uword			ridx;
	umat			CSIRS_position_temp;


	this->cell_ID = para.cell_ID;
	this->CSIRS_position = zeros<ucube>(para.num_RB * para.num_sc, para.num_symb, para.num_port_CSIRS);

	if (ind_slot % para.CSIRS_periodicity == 0) {
		this->resource_grid_CSIRS = zeros<cx_cube>(para.num_RB * para.num_sc, para.num_symb, para.num_port_CSIRS);

		k << 4 << 8 << 0 << 2 << 6 << 10 << endr;		// bitmap value
		l << 5 << 9 << endr;							// bitmap value

		switch (para.row_CSIRS) {		// determine CSIRS setting with row numbering
		case 1:
			func_CDM(CSIRS_CDMType::NO_CDM, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << 2 << 6 << 10 << endr;
			time_pos << l(0) << endr;
			break;
		case 2:
			func_CDM(CSIRS_CDMType::NO_CDM, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << 6 << endr;
			time_pos << l(0) << endr;
			break;
		case 3:
			func_CDM(CSIRS_CDMType::FD_CDM2, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << endr;
			time_pos << l(0) << endr;
			break;
		case 4:
			func_CDM(CSIRS_CDMType::FD_CDM2, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << 4 << 6 << endr;
			time_pos << l(0) << endr;
			break;
		case 5:
			func_CDM(CSIRS_CDMType::FD_CDM2, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << endr;
			time_pos << l(0) << l(0) + 1 << endr;
			break;
		case 7:
			func_CDM(CSIRS_CDMType::FD_CDM2, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << k(3) << endr;
			time_pos << l(0) << endr;
			break;
		case 8:
			func_CDM(CSIRS_CDMType::FD_CDM2, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << endr;
			time_pos << l(0) << l(0) + 1 << endr;
			break;
		case 9:
			func_CDM(CSIRS_CDMType::CDM4, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << endr;
			time_pos << l(0) << endr;
			break;
		case 10:
			func_CDM(CSIRS_CDMType::FD_CDM2, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << k(3) << k(4) << k(5) << endr;
			time_pos << l(0) << endr;
			break;
		case 11:
			func_CDM(CSIRS_CDMType::CDM4, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << endr;
			time_pos << l(0) << endr;
			break;
		case 12:
			func_CDM(CSIRS_CDMType::FD_CDM2, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << k(3) << endr;
			time_pos << l(0) << l(0) + 1 << endr;
			break;
		case 13:
			func_CDM(CSIRS_CDMType::CDM4, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << k(3) << endr;
			time_pos << l(0) << endr;
			break;
		case 14:
			func_CDM(CSIRS_CDMType::FD_CDM2, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << endr;
			time_pos << l(0) << l(0) + 1 << l(1) << l(1) + 1 << endr;
			break;
		case 15:
			func_CDM(CSIRS_CDMType::CDM4, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << endr;
			time_pos << l(0) << l(1) << endr;
			break;
		case 16:
			func_CDM(CSIRS_CDMType::CDM8, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << endr;
			time_pos << l(0) << endr;
			break;
		case 17:
			func_CDM(CSIRS_CDMType::FD_CDM2, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << k(3) << endr;
			time_pos << l(0) << l(0) + 1 << l(1) << l(1) + 1 << endr;
			break;
		case 18:
			func_CDM(CSIRS_CDMType::CDM4, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << k(3) << endr;
			time_pos << l(0) << l(1) << endr;
			break;
		case 19:
			func_CDM(CSIRS_CDMType::CDM8, CSIRS_w_f, CSIRS_w_t, CDM_length_f, CDM_length_t);
			freq_pos << k(0) << k(1) << k(2) << k(3) << endr;
			time_pos << l(0) << endr;
			break;
		}   /* end of switch (para->row_CSIRS) */

		CDM_length = CDM_length_f * CDM_length_t;		// calculating length of CDM

		RB_position = zeros<icube>(para.num_sc, para.num_symb, freq_pos.n_elem * time_pos.n_elem * CDM_length);
		for (uword ind_freq = 0; ind_freq < freq_pos.n_elem; ind_freq++)
			for (uword ind_symb = 0; ind_symb < time_pos.n_elem; ind_symb++)
				for (uword ind_CDM = 0; ind_CDM < CDM_length; ind_CDM++) {
					RB_position(span(freq_pos(ind_freq), freq_pos(ind_freq) + CDM_length_f - 1),
						span(time_pos(ind_symb), time_pos(ind_symb) + CDM_length_t - 1),
						span((ind_freq * time_pos.n_elem + ind_symb) * CDM_length + ind_CDM))
						= CSIRS_w_f.row(ind_CDM).t() * CSIRS_w_t.row(ind_CDM);
				}

		CSIRS_density = para.CSIRS_density;				// Density of CSIRS
		if (para.row_CSIRS == 1) {
			RB_position_tot = zeros<imat>(para.num_sc, para.num_symb);
			for (uword i = 0; i < freq_pos.n_elem * time_pos.n_elem * CDM_length; i++)
				RB_position_tot = RB_position_tot + RB_position.slice(i);

			CSIRS_density = 1;
		}
		
		/* pseudo sequence for CSIRS */
		pseudo_sequence_CSIRS = field<cx_mat>(1, para.num_symb);
		for (uword ind_symb = 0; ind_symb < time_pos.n_elem; ind_symb++)
			for (uword ind_time = 0; ind_time < CDM_length_t; ind_time++) {
				c_init_CSIRS = mod((pow(2, 10) * (14 * (ind_slot) + time_pos(ind_symb) + ind_time + 1) * (2 * para.cell_ID + 1) + para.cell_ID), pow(2, 31));
				func_pseudo_sequence_generation(c_init_CSIRS, 2 * CDM_length_f * para.num_RB * para.CSIRS_density, c);
				pseudo_sequence_CSIRS(time_pos(ind_symb) + ind_time) = zeros<cx_mat>(1, c.n_elem / 2);
				for (uword i = 0, ridx = 0, iidx = 1; ridx < c.n_elem && iidx < c.n_elem; i++, ridx += 2, iidx += 2)
					pseudo_sequence_CSIRS(time_pos(ind_symb) + ind_time)(0, i) = ((1 - 2 * c(0, ridx)) + cx_double(0, 1) * (1 - 2 * c(0, iidx))) / sqrt(2);
			}

		/* Inserting pseudo sequence for CSIRS int the right position */
		this->CSIRS_position_total = zeros<umat>(para.num_RB * para.num_sc, para.num_symb);
		for (uword ind_port = 0; ind_port < para.num_port_CSIRS; ind_port++) {
			for (uword ind_RB = 0; ind_RB < para.num_RB * CSIRS_density; ind_RB++) {	// For the case of CSIRS_desity = 0.5
				if (para.row_CSIRS == 1) {
					for (uword i = 0; i < para.num_sc; i++) {
						ridx = i + ind_RB * para.num_sc;
						for (uword j = 0; j < para.num_symb; j++)
							this->CSIRS_position(ridx, j, ind_port) = abs(RB_position_tot(i, j));
					}
				}
				else {
					for (uword i = 0; i < para.num_sc; i++) {
						ridx = i + (ind_RB / CSIRS_density) * para.num_sc;
						for (uword j = 0; j < para.num_symb; j++)
							this->CSIRS_position(ridx, j, ind_port) = abs(RB_position(i, j, ind_port));
					}
				}
			}

			index_for_time = find(RB_position.slice(ind_port));
			ind_symb = ceil(index_for_time(0) / para.num_sc);

			temp = zeros<cx_colvec>(this->CSIRS_position.n_rows);
			for (uword ind_time = 0; ind_time < CDM_length_t; ind_time++) {
				for (uword i = 0; i < this->resource_grid_CSIRS.n_rows; i++)
					temp(i) = this->resource_grid_CSIRS(i, ind_symb + ind_time, ind_port);
				for (uword i = 0, j = 0; i < this->CSIRS_position.n_rows; i++)
					if (this->CSIRS_position(i, ind_symb + ind_time, ind_port) == 1)
						temp(i) = pseudo_sequence_CSIRS(ind_symb + ind_time)(j++);
				for (uword i = 0; i < temp.n_rows; i++)
					this->resource_grid_CSIRS(i, ind_symb + ind_time, ind_port) = temp(i);
			}

			if (para.row_CSIRS == 1)
				this->resource_grid_CSIRS.slice(ind_port) = this->resource_grid_CSIRS.slice(ind_port) % repmat(RB_position_tot, para.num_RB, 1);
			else
				this->resource_grid_CSIRS.slice(ind_port) = this->resource_grid_CSIRS.slice(ind_port) % repmat(RB_position.slice(ind_port), para.num_RB, 1);

			CSIRS_position_temp = zeros<umat>(this->CSIRS_position.n_rows, this->CSIRS_position.n_cols);
			for (uword i = 0; i < CSIRS_position_temp.n_rows; i++)		// squeeze this->CSIRS_postion
				for (uword j = 0; j < CSIRS_position_temp.n_cols; j++)
					CSIRS_position_temp(i, j) = this->CSIRS_position(i, j, ind_port);

			this->CSIRS_position_total = this->CSIRS_position_total || CSIRS_position_temp;

		}	/* end of for (uword ind_port = 0; ind_port < para.num_port_CSIRS; ind_port++) loop */

	}
	else {
		this->resource_grid_CSIRS = zeros<cx_cube>(para.num_RB * para.num_sc, para.num_symb, para.num_port_CSIRS);
		this->CSIRS_position = zeros<ucube>(para.num_RB * para.num_sc, para.num_symb, para.num_port_CSIRS);
	}

	/* local matrix memory deallocation */
	k.reset();
	l.reset();
	freq_pos.reset();
	time_pos.reset();
	CSIRS_w_f.reset();
	CSIRS_w_t.reset();
	RB_position.reset();
	RB_position_tot.reset();
	for (uword i = 0; i < pseudo_sequence_CSIRS.n_elem; i++)
		pseudo_sequence_CSIRS(i).reset();
	pseudo_sequence_CSIRS.reset();
	c.reset();
	index_for_time.reset();
	temp.reset();
	CSIRS_position_temp.reset();

}	/* end of method_CSI_RS_allocation() method definition */




/*
 * method_resource_allocation() -resource allocation
 */
void ModuleBS_SUMIMO::method_resource_allocation(ModuleParameterMIMO &para, int &ind_slot)
{
	/* local variable declarations */
	int		num_layers;

	/* load CQI_table from MATLAB MAT-File */
	this->CQI_table_num = para.MCS_table_index;
	switch (this->CQI_table_num) {
	case 1:
		if (!this->CQI.empty()) {
			this->M_order = CQI_Table_1.order(this->CQI[0] - 1);
			this->Coding_rate = CQI_Table_1.code_rate(this->CQI[0] - 1) / 1024.0;
		}
		break;
	case 2:
		if (!this->CQI.empty()) {
			this->M_order = CQI_Table_2.order(this->CQI[0] - 1);
			this->Coding_rate = CQI_Table_2.code_rate(this->CQI[0] - 1) / 1024.0;
		}
		break;
	default:
		perror("MCS table index not supported");
	}

	/* PSS and SSS position, PBCH, and PDCCH position */
	this->PSS_5G_position		= zeros<umat>(para.num_RB * para.num_sc, para.num_symb);
	this->SSS_5G_position		= zeros<umat>(para.num_RB * para.num_sc, para.num_symb);
	this->PBCH_5G_position		= zeros<umat>(para.num_RB * para.num_sc, para.num_symb);
	this->PDCCH_5G_position		= zeros<umat>(para.num_RB * para.num_sc, para.num_symb);

	/* PDCCH for 5G	  http://niviuk.free.fr/nr_frame.php */
	if (para.use_PDCCH)		// code comment: para.use_PDCCH = 0 => if 문장을 실행하지 않고 다음으로 넘어감.
		for (uword ll = 0; ll < para.num_slot_in_subframe; ll++)
			this->PDCCH_5G_position(span::all, span(ll * 14)) = true;	// MATLAB code: obj.PDCCH_5G_position(:,1+(ll-1)*14,ind_port) => 변수 ind_port 지정이 없음. (문의해 볼 것)

	/* PSS and SSS for 5G: 38.211 */
	this->SS_PBCH_block_start.reset();
	switch (para.mu) {
	case 0:		/* 38.213 4.1 Case A */
		if (para.f <= 3e9)
			switch (mod(ind_slot + 1, 10 * pow(2, para.mu))) {
			case 1:
			case 2:
				this->SS_PBCH_block_start << 2 << 8 << endr;
			}
		else if (para.f > 3e9 && para.f <= 6e9)
			switch (mod(ind_slot + 1, 10 * pow(2, para.mu))) {
			case 1:
			case 2:
			case 3:
			case 4:
				this->SS_PBCH_block_start << 2 << 8 << endr;
			}
		break;
	case 1:		/* 38.213 Case 4.1 Case B */
		if (para.f <= 3e9)
			switch (mod(ind_slot + 1, 10 * pow(2, para.mu))) {
			case 1:
				this->SS_PBCH_block_start << 4 << 8 << endr;
				break;
			case 2:
				this->SS_PBCH_block_start << 16 - 14 << 20 - 14 << endr;
				break;
			}
		else if (para.f > 3e9 && para.f < 6e9)
			switch (mod(ind_slot + 1, 10 * pow(2, para.mu))) {
			case 1:
			case 3:
				this->SS_PBCH_block_start << 4 << 8 << endr;
				break;
			case 2:
			case 4:
				this->SS_PBCH_block_start << 16 - 14 << 20 - 14 << endr;
				break;
			}
		break;
	case 2:		/* 38.213 4.1 Case C */
		if (para.f <= 3e9)
			switch (mod(ind_slot + 1, 10 * pow(2, para.mu))) {
			case 1:
			case 2:
				this->SS_PBCH_block_start << 2 << 8 << endr;
			}
		else if (para.f > 3e9 && para.f <= 6e9)
			switch (mod(ind_slot + 1, 10 * pow(2, para.mu))) {
			case 1:
			case 2:
			case 3:
			case 4:
				this->SS_PBCH_block_start << 2 << 8 << endr;
			}
		break;
	case 3:		/* 38.213 4.1 Case D */
		if (para.f > 6e9)
			switch (mod(ind_slot + 1, pow(10 * 2, para.mu))) {
			case 1:
			case 3:
			case 5:
			case 7:
			case 11:
			case 13:
			case 15:
			case 17:
			case 21:
			case 23:
			case 25:
			case 27:
			case 31:
			case 33:
			case 35:
			case 37:
				this->SS_PBCH_block_start << 4 << 8 << endr;
				break;
			case 2:
			case 4:
			case 6:
			case 8:
			case 12:
			case 14:
			case 16:
			case 18:
			case 22:
			case 24:
			case 26:
			case 28:
			case 32:
			case 34:
			case 36:
			case 38:
				this->SS_PBCH_block_start << 16 - 14 << 20 - 14 << endr;
			}
	}	/* end of switch (para->mu) */

	/* PSS, SSS, PBCH, DM-RS for PBCH */
	if (this->SS_PBCH_block_start.n_elem != 0) {
		uvec sn(this->SS_PBCH_block_start.n_elem);
		for (uword i = 0; i < sn.n_elem; i++)
			sn(i) = (unsigned int) this->SS_PBCH_block_start(i);		// type conversion

		for (uword i = 56; i < 183; i++)
			for (uword j = 0; j < sn.n_elem; j++) {
				this->PSS_5G_position(i, sn(j) - 1 + 1) = true;			// l = 0
				this->SSS_5G_position(i, sn(j) - 1 + 3) = true;			// l = 2
			}

		for (uword i = 0; i < 240; i++) 								// l = 1
			for (uword j = 0; j < sn.n_elem; j++) {
				this->PBCH_5G_position(i, sn(j) - 1 + 2) = true;
				this->PBCH_5G_position(i, sn(j) - 1 + 4) = true;
				if (i < 48 || i >= 192)
					this->PBCH_5G_position(i, sn(j) - 1 + 3) = true;	// l = 2
			}

		int v = mod(this->cell_ID, 4);

		this->PBCH_DMRS_5G_position = zeros<umat>(237 + v, sn(sn.n_elem - 1) + 4);
		for (uword i = 0 + v; i < 237 + v; i += 4)
			for (uword j = 0; j < sn.n_elem; j++) {
				this->PBCH_DMRS_5G_position(i, sn(j) - 1 + 2) = true;	// l = 1
				this->PBCH_DMRS_5G_position(i, sn(j) - 1 + 4) = true;	// l = 3
			}

		for (uword i = 0 + v; i < 45 + v; i += 4)
			for (uword k = 0; k < sn.n_elem; k++)
				this->PBCH_DMRS_5G_position(i, sn(k) - 1 + 3) = true;	// l = 2

		for (uword i = 192 + v; i < 237 + v; i += 4)
			for (uword k = 0; k < sn.n_elem; k++)
				this->PBCH_DMRS_5G_position(i, sn(k) - 1 + 3) = true;	// l = 2

		/* local memory deallocation */
		sn.reset();
	}

}	/* end of method_resource_allocation() method definition */




/*
 * method_beam_management() - Beam management
 */
void ModuleBS_SUMIMO::method_beam_management(ModuleParameterMIMO &para, ChannelOutput &ch_Output, int &ind_slot)
{
	/* local variable decalaration */
	int				O1;
	mat				Tx_theta_ver;
	cx_mat			TxB_ver;
	rowvec			TxOpA, TxOpB;
	cx_cube			rsp;
	cx_mat			squeezed_fft;
	double			max_ind_beam;
	cx_cube			sumRes;
	cx_rowvec		squeezedRes;
	field<cx_cube>	H_effective;
	cx_mat			H_effective_RB;
	cx_mat			U;
	vec				s;
	cx_mat			V;
	uword			ws;
	umat			WB_pos;
	uword			pos_size;
	mat				rsv_norm;

	cx_cube			ts_CSIRS_signal;
	cx_mat			temp_ts_signal;
	mat				UE_Beam;
	cx_cube			rv_CSIRS_signal;
	cx_colvec		squeezed_ts_sig;
	cx_cube			ch_est;
	uvec			pos_freq;
	cx_mat			temp_grid;
	cx_cube			temp_rv_CSIRS_signal;
	cx_mat			temp_CSIRS_signal;
	cx_mat			ttmp_CSIRS_signal;


	/* for Calibration */
	if (para.cal_scenario) {
		switch (para.scenario_beam) {
		case 1:
			O1 = 1;
			TxOpA = zeros<rowvec>(para.N2);
			for (uword i = 0; i < para.N2; i++)
				TxOpA(i) = i;
			TxOpB = zeros<rowvec>(para.N2 * O1);
			for (uword i = 0; i < para.N2 * O1; i++)
				TxOpB(i) = i + 1;
			Tx_theta_ver = TxOpA.t() * (datum::pi * TxOpB - 1.0 / 2.0) / (double)para.N2 / (double)O1;
			TxB_ver = (exp(cx_double(0, -1) * datum::pi * Tx_theta_ver) / sqrt((double)para.N2));
			break;
		case 2:
			TxOpA = zeros<rowvec>(para.N2);
			for (uword i = 0; i < para.N2; i++)
				TxOpA(i) = i;
			TxOpB = zeros<rowvec>(para.N2);
			for (uword i = 0; i < para.N2; i++)
				TxOpB(i) = i + 1;
			Tx_theta_ver = TxOpA.t() * (datum::pi * TxOpB - 1.0 / 2.0) / (double)para.N2;
			TxB_ver = (exp(cx_double(0, -1) * datum::pi * Tx_theta_ver) / sqrt((double)para.N2));
			break;
		}

		/* Initialize narrow beam */
		this->Narrowbeam = field<cx_mat>(TxB_ver.n_cols, 1);
		for (uword ind_ver = 0; ind_ver < TxB_ver.n_cols; ind_ver++) 
			this->Narrowbeam(ind_ver) = kron(eye(para.Tx_pol * para.N1, para.Tx_pol * para.N1), TxB_ver(span::all, ind_ver));
		
		/* Calculate rsp(received signal power) per beam */	
		rsp = zeros<cx_cube>(para.num_RB * para.num_sc, para.num_symb, TxB_ver.n_cols);
		for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {	// squeeze ch_Output.fft 4D-matrix to 2D-matrix
			squeezed_fft = zeros<cx_mat>(ch_Output.fft(0).n_rows, ch_Output.fft(0).n_slices);
			for (uword i = 0; i < ch_Output.fft(0).n_cols; i++)
				for (uword j = 0; j < ch_Output.fft(0).n_slices; j++)
					squeezed_fft(i, j) = ch_Output.fft(ind_RB * para.num_sc + 5)(6, i, j);
			
			for (uword ind_beam = 0; ind_beam < this->Narrowbeam.n_elem; ind_beam++)
				if (para.num_Rx_antenna == 1)
					rsp((ind_RB) * para.num_sc + 6, 7, ind_beam) = pow(norm(squeezed_fft.t() * this->Narrowbeam(ind_beam), 2.0), 2.0);
				else
					rsp((ind_RB) * para.num_sc + 6, 7, ind_beam) = pow(norm(squeezed_fft * this->Narrowbeam(ind_beam), 2.0), 2.0);
		}

		sumRes = sum(sum(rsp, 0), 1);
		squeezedRes = zeros<cx_rowvec>(sumRes.n_slices);
		squeezedRes = sumRes(span(0), span(0), span::all);
		max_ind_beam = index_max(squeezedRes);

		/* H_effective */
		H_effective = field<cx_cube>(para.num_RB * para.num_sc);
		for (uword i = 0; i < para.num_RB * para.num_sc; i++)
			H_effective(i) = zeros<cx_cube>(para.num_symb, para.num_Rx_antenna, para.Tx_pol * para.N1);

		for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
			squeezed_fft = zeros<cx_mat>(ch_Output.fft(0).n_cols, ch_Output.fft(0).n_slices);
			for (uword i = 0; i < ch_Output.fft(0).n_cols; i++)
				for (uword j = 0; j < ch_Output.fft(0).n_slices; j++)
					squeezed_fft(i, j) = ch_Output.fft(ind_RB * para.num_sc + 5)(6, i, j);

			cx_mat temp;
			if (para.num_Rx_antenna == 1) {
				//H_effective(ind_RB * para.num_sc + 5)(span(6), span::all, span::all) = squeezed_fft.t() * this->Narrowbeam(max_ind_beam);
				temp = squeezed_fft.t() * this->Narrowbeam(max_ind_beam);
			}
			else {
				//H_effective(ind_RB * para.num_sc + 5)(span(6), span::all, span::all) = squeezed_fft * this->Narrowbeam(max_ind_beam);// (span::all, span::all);
				temp = squeezed_fft * this->Narrowbeam(max_ind_beam);
			}
			for (uword i = 0; i < 2; i++)
				for (uword j = 0; j < 4; j++)
					H_effective(ind_RB * para.num_sc + 5)(6, i, j) = temp(i, j);
		}

		/* Precoder */
		this->Precoder = field<cx_mat>(para.num_RB);		// 2018-10-23 17:25
		//H_effective_RB = zeros<cx_mat>(H_effective(0).n_cols, H_effective(0).n_slices);
		for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
			if (para.num_Rx_antenna == 1) {
				H_effective_RB = H_effective(para.num_sc * ind_RB + 5)(span(6), span::all, span::all);
				H_effective_RB = H_effective_RB.t();
			}
			else
				H_effective_RB = H_effective(para.num_sc * ind_RB + 5)(span(6), span::all, span::all);
			
			svd(U, s, V, H_effective_RB);
			this->Precoder(ind_RB) = this->Narrowbeam(max_ind_beam) * V.col(0);
			this->Precoder(ind_RB) = this->Precoder(ind_RB) / norm(this->Precoder(ind_RB));
		}
		this->RI = 1;
	}
	else {	/* Wide Beam Selection */
		/* Wide beam generation */
		func_wide_beam_gen(para, this->Widebeam, this->Narrowbeam);

		ws = this->Widebeam.n_rows;
		WB_pos = this->PSS_5G_position || this->SSS_5G_position || this->PBCH_5G_position;
		pos_size = sum(sum(WB_pos));
		rsp = zeros<cx_cube>(para.num_Rx_antenna, pos_size, ws);
		rsv_norm = zeros<mat>(ws, 1);

		/* Calculate recevied signal power per WB */
		field<cx_cube> ch_Output_temp(ch_Output.fft(0).n_cols);
		for (uword i = 0; i < ch_Output_temp.n_elem; i++)
			ch_Output_temp(i) = zeros<cx_cube>(ch_Output.fft(0).n_slices, ch_Output.fft.n_elem, ch_Output.fft(0).n_rows);

		for (uword i = 0; i < ch_Output.fft.n_elem; i++)
			for (uword j = 0; j < ch_Output.fft(0).n_rows; j++)
				for (uword k = 0; k < ch_Output.fft(0).n_cols; k++)
					for (uword m = 0; m < ch_Output.fft(0).n_slices; m++)
						ch_Output_temp(k)(m, i, j) = ch_Output.fft(i)(j, k, m);

		uword n;
		cx_cube ch_Output_pos;
		ch_Output_pos = zeros<cx_cube>(ch_Output_temp.n_elem, ch_Output_temp(0).n_rows, sum(sum(WB_pos)));
		for (uword i = 0; i < ch_Output_pos.n_rows; i++)
			for (uword j = 0; j < ch_Output_pos.n_cols; j++) {
				n = 0;
				for (uword m = 0; m < WB_pos.n_cols; m++)
					for (uword k = 0; k < WB_pos.n_rows; k++)
						if (WB_pos(k, m) != 0)
							ch_Output_pos(i, j, n++) = ch_Output_temp(i)(j, k, m);
			}
		
		cx_colvec res(rsp.n_rows);
		cx_colvec rspvec(rsp.n_rows * rsp.n_cols);
		for (uword ind_wi = 0; ind_wi < ws; ind_wi++) {
			for (uword ind_pos = 0; ind_pos < pos_size; ind_pos++) {
				res = ch_Output_pos.slice(ind_pos) * this->Widebeam(ind_wi);
				for (uword i = 0; i < rsp.n_rows; i++)
					rsp(i, ind_pos, ind_wi) = res(i);
			}
			n = 0;
			for (uword i = 0; i < rsp.n_rows; i++)
				for (uword j = 0; j < rsp.n_cols; j++)
					rspvec(n++) = rsp(i, j, ind_wi);
			rsv_norm(ind_wi) = norm(rspvec, 2);
		}

		colvec maxval;
		maxval = max(rsv_norm, 0);
		for (uword i = 0; i < rsv_norm.n_elem; i++)
			if (maxval(0) == rsv_norm(i)) {
				this->ind_WB = i;
				break;
			}

		/* Narrow Beam Selection */
		uword s1, s2;
		s1 = this->Narrowbeam(this->ind_WB).n_rows;
		s2 = this->Narrowbeam(this->ind_WB).n_cols;
		this->Precoder_CSIRS = zeros<cx_mat>(s1, s2);
		this->Precoder_CSIRS = this->Narrowbeam(this->ind_WB);

		/* Virtual transmission of CSI-RS signal (ts_CSIRS_signal) */
		ts_CSIRS_signal = zeros<cx_cube>(s1, para.num_RB * para.num_sc, para.num_symb);
		for (uword ind_Nt = 0; ind_Nt < para.num_Tx_antenna; ind_Nt++) {
			temp_ts_signal = zeros<cx_mat>(para.num_RB * para.num_sc, para.num_symb);
			for (uword ind_port = 0; ind_port < s2; ind_port++)
				temp_ts_signal = temp_ts_signal + this->Precoder_CSIRS(ind_Nt, ind_port) * this->resource_grid_CSIRS.slice(ind_port);
			ts_CSIRS_signal(span(ind_Nt), span::all, span::all) = temp_ts_signal(span::all, span::all);
		}
		temp_ts_signal.reset();

		/* Virtual reception of CSI-RS signal (rv_CSIRS_signal) */
		UE_Beam = eye(para.num_Rx_antenna, para.num_Rx_antenna);
		rv_CSIRS_signal = zeros<cx_cube>(UE_Beam.n_cols, para.num_RB * para.num_sc, para.num_symb);
		squeezed_fft = zeros<cx_mat>(ch_Output.fft(0).n_cols, ch_Output.fft(0).n_slices);
		squeezed_ts_sig = zeros<cx_colvec>(ts_CSIRS_signal.n_rows);
		for (uword ind_sc = 0; ind_sc < para.num_RB * para.num_sc; ind_sc++)
			for (uword ind_symb = 0; ind_symb < para.num_symb; ind_symb++) {
				for (uword i = 0; i < squeezed_fft.n_rows; i++)
					for (uword j = 0; j < squeezed_fft.n_cols; j++)
						squeezed_fft(i, j) = ch_Output.fft(ind_sc)(ind_symb, i, j);

				for (uword i = 0; i < ts_CSIRS_signal.n_rows; i++)
					squeezed_ts_sig(i) = ts_CSIRS_signal(i, ind_sc, ind_symb);

				if (para.num_Rx_antenna == 1)
					rv_CSIRS_signal(span::all, span(ind_sc), span(ind_symb)) = UE_Beam * squeezed_fft.t() * squeezed_ts_sig;
				else
					rv_CSIRS_signal(span::all, span(ind_sc), span(ind_symb)) = UE_Beam * squeezed_fft * squeezed_ts_sig;
			}

		/* Beam channel estimation (ch_est) : 2018-10-31 01:06 MATLAB Code line 470 */
		ch_est = zeros<cx_cube>(para.num_Rx_antenna, para.num_RB, s2);
		temp_grid = zeros<cx_mat>(para.num_sc, this->resource_grid_CSIRS.n_cols);
		pos_freq = zeros<uvec>(para.num_sc);
		temp_rv_CSIRS_signal = zeros<cx_cube>(rv_CSIRS_signal.n_rows, pos_freq.n_elem, rv_CSIRS_signal.n_slices);

		umat temp_CSIRS_position = zeros<umat>(pos_freq.n_elem, this->CSIRS_position.n_cols);
		
		for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
			for (uword i = 0; i < para.num_sc; i++)
				pos_freq(i) = para.num_sc * ind_RB + i;

			for (uword ind_port = 0; ind_port < s2; ind_port++) {
				for (uword i = 0; i < pos_freq.n_elem; i++)
					for (uword j = 0; j < this->resource_grid_CSIRS.n_cols; j++)
						temp_grid(i, j) = this->resource_grid_CSIRS(pos_freq(i), j, ind_port);

				for (uword i = 0; i < rv_CSIRS_signal.n_rows; i++)
					for (uword j = 0; j < pos_freq.n_elem; j++)
						for (uword k = 0; k < rv_CSIRS_signal.n_slices; k++)
							temp_rv_CSIRS_signal(i, j, k) = rv_CSIRS_signal(i, pos_freq(j), k);

				for (uword i = 0; i < pos_freq.n_elem; i++)
					for (uword j = 0; j < this->CSIRS_position.n_cols; j++)
						temp_CSIRS_position(i, j) = this->CSIRS_position(pos_freq(i), j, ind_port);

				temp_CSIRS_signal = zeros<cx_mat>(arma::size(find(temp_CSIRS_position)));

				umat	sq_CSIRS_position = zeros<umat>(pos_freq.n_elem, this->CSIRS_position.n_cols);
				for (uword i = 0; i < sq_CSIRS_position.n_rows; i++)
					for (uword j = 0; j < sq_CSIRS_position.n_cols; j++)
						sq_CSIRS_position(i, j) = this->CSIRS_position(pos_freq(i), j, ind_port);

				cx_mat sq_temp_rv_CSIRS_signal(size(find(sq_CSIRS_position)));
				cx_mat CTemp(1, 1);
				for (uword ind_Nr = 0; ind_Nr < para.num_Rx_antenna; ind_Nr++) {
					temp_CSIRS_signal = temp_grid(find(temp_CSIRS_position));
					for (uword j = 0, n = 0; j < sq_CSIRS_position.n_cols; j++)
						for (uword i = 0; i < sq_CSIRS_position.n_rows; i++)
							if (sq_CSIRS_position(i, j) == 1)
								sq_temp_rv_CSIRS_signal(n++) = temp_rv_CSIRS_signal(ind_Nr, i, j);
					
					CTemp = temp_CSIRS_signal.t() * sq_temp_rv_CSIRS_signal;
					ch_est(ind_Nr, ind_RB, ind_port) = CTemp(0, 0);
				}
			}
		}

		/* Beam selection & feedback */
		func_beam_selection(para, ch_est, this->Precoder_CSIRS, this->RI, this->Precoder);

		if (imod(ind_slot, para.CSIRS_periodicity) != 0)
			this->CSIRS_position_total = zeros<umat>(para.num_RB * para.num_sc, para.num_symb);

	}

	/* local memory deallocation */
	Tx_theta_ver.reset();
	TxB_ver.reset();
	TxOpA.reset();
	TxOpB.reset();
	rsp.reset();
	squeezed_fft.reset();
	sumRes.reset();
	squeezedRes.reset();
	for (uword i = 0; i < H_effective.n_elem; i++)
		H_effective(i).reset();
	H_effective.reset();
	H_effective_RB.reset();
	U.reset();
	s.reset();
	V.reset();
	WB_pos.reset();
	rsp.reset();
	rsv_norm.reset();
	ts_CSIRS_signal.reset();
	temp_ts_signal.reset();
	UE_Beam.reset();
	rv_CSIRS_signal.reset();
	squeezed_ts_sig.reset();
	ch_est.reset();
	pos_freq.reset();
	temp_grid.reset();
	temp_rv_CSIRS_signal.reset();
	temp_CSIRS_signal.reset();
	ttmp_CSIRS_signal.reset();

}	/* end of method_beam_management() method definition */




/*
 * method_DMRS_allocation() - DMRS allocation.
 */
void ModuleBS_SUMIMO::method_DMRS_allocation(ModuleParameterMIMO &para, int &ind_slot) {

	/* local variable declarations */
	imat	w_f, w_t;
	uvec	delta, k_prime, l_prime;
	uword	l_0;
	uvec	n, k0, k1, k0_port, k1_port;
	uword	symb_num;
	uword	DMRS_length;
	int		c_init_DMRS;
	mat		c;
	mat		c1, c2;
	cx_mat	pseudo_sequence_DMRS;
	uword	DMRS_legnth;


	para.num_port_DMRS = zeros<irowvec>(this->RI + 1);
	for (uword i = 0; i < this->RI + 1; i++)
		para.num_port_DMRS(i) = para.port_set_DMRS(i);
	this->DMRS_position = zeros<ucube>(para.num_RB * para.num_sc, para.num_symb, this->RI);
	this->resource_grid_DMRS = zeros<cx_cube>(para.num_RB * para.num_sc, para.num_symb, this->RI);

	switch (para.DMRS_type) {
	case 1:
		w_f << 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr;
		w_t << 1 << 1 << endr
			<< 1 << 1 << endr
			<< 1 << 1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << -1 << endr
			<< 1 << -1 << endr
			<< 1 << -1 << endr;
		delta << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << endr;
		k_prime << 0 << 1 << endr;
		l_prime << 0 << 1 << endr;

		if (para.PDSCH_mapping_type == PDSCH_mapping_type::A)
			if (para.DL_DMRS_typeA_pos == 3)
				l_0 = 3;
			else
				l_0 = 2;
		else if (para.PDSCH_mapping_type == PDSCH_mapping_type::B)
			if (this->PDCCH_5G_position(0, 1))
				l_0 = 2;
			else
				l_0 = 1;
		else
			perror("PDSCH_mapping_type not supported");

		n = zeros<uvec>(3 * para.num_RB);
		k0 = zeros<uvec>(n.n_elem);
		k1 = zeros<uvec>(n.n_elem);
		for (uword i = 0; i < 3 * para.num_RB; i++)
			n(i) = i;
		k0 = 4 * n + 2 * k_prime(0);
		k1 = 4 * n + 2 * k_prime(1);

		if (this->RI < 5)
			symb_num = 1;
		else
			symb_num = 2;

		break;
	case 2:
		w_f << 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr;
		w_t << 1 << 1 << endr
			<< 1 << 1 << endr
			<< 1 << 1 << endr
			<< 1 << 1 << endr
			<< 1 << 1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << -1 << endr
			<< 1 << -1 << endr
			<< 1 << -1 << endr
			<< 1 << -1 << endr
			<< 1 << -1 << endr;
		delta << 0 << 0 << 2 << 2 << 4 << 4 << 0 << 0 << 2 << 2 << 4 << 4 << endr;
		k_prime << 0 << 1 << endr;
		l_prime << 0 << 1 << endr;

		if (para.PDSCH_mapping_type == PDSCH_mapping_type::A)
			if (para.DL_DMRS_typeA_pos == 3)
				l_0 = 3;
			else
				l_0 = 2;
		else if (para.PDSCH_mapping_type == PDSCH_mapping_type::B)
			if (this->PDCCH_5G_position(0, 1))
				l_0 = 2;
			else
				l_0 = 1;
		else
			perror("PDSCH_mapping_type not supported");

		n = zeros<uvec>(2 * para.num_RB);
		k0 = zeros<uvec>(n.n_elem);
		k1 = zeros<uvec>(n.n_elem);
		for (uword i = 0; i < 2 * para.num_RB - 1; i++)
			n(i) = i;
		k0 = 6 * n + 2 * k_prime(0);
		k1 = 6 * n + 2 * k_prime(1);

		if (this->RI < 7)
			symb_num = 1;
		else
			symb_num = 2;

	}	/* end of switch (para.DMRS_type) */


	if (para.cal_scenario) {	/* for calibration scenario */
		this->cell_ID = para.cell_ID;
		this->PDCCH_5G_position = zeros<umat>(para.num_RB * para.num_sc, para.num_symb);
		this->PDCCH_5G_position(span::all, span(0, 1)) = true;
		para.num_port_CSIRS = 1;
		this->CSIRS_position = zeros<ucube>(para.num_RB * para.num_sc, para.num_symb, 1);

		symb_num = 2;
		this->RI = 1;
		l_0 = 2;
		DMRS_length = 2 * n(n.n_elem - 1) + 2;
		this->DMRS_position_total = zeros<umat>(para.num_RB * para.num_sc, para.num_symb);

		for (uword ind_port = 0; ind_port < this->RI; ind_port++) {
			k0_port = k0 + delta(ind_port);
			k1_port = k1 + delta(ind_port);
			for (uword ind_symb = 0; ind_symb < symb_num; ind_symb++) {
				c_init_DMRS = imod( (pow(2, 17) * (14 * ind_slot + (l_0 + l_prime(ind_symb)) + 1) * (2 * this->cell_ID + 1) + 2 * this->cell_ID), pow(2, 31) );
				func_pseudo_sequence_generation(c_init_DMRS, 2 * DMRS_length, c);
				c1 = zeros<mat>(c.n_rows, c.n_cols / 2);
				c2 = zeros<mat>(c.n_rows, c.n_cols / 2);
				for (uword i = 0, idx = 0; i < c.n_cols; i += 2, idx++) {
					c1.col(idx) = c.col(i);
					c2.col(idx) = c.col(i + 1);
				}
				pseudo_sequence_DMRS = ((1 - 2 * c1) + cx_double(0, 1) * (1 - 2 * c2)) / sqrt(2);
				for (uword i = 0; i < k0_port.n_elem; i++) {
					this->DMRS_position(k0_port(i), l_0 + l_prime(ind_symb), ind_port) = true;
					this->DMRS_position(k1_port(i), l_0 + l_prime(ind_symb), ind_port) = true;
					for (uword j = 0; j < pseudo_sequence_DMRS.n_elem; j += 2) {
						this->resource_grid_DMRS(k0_port(i), l_0 + l_prime(ind_symb), ind_port) = (double)(w_f(ind_port, 0) * w_t(ind_port, ind_symb)) * pseudo_sequence_DMRS(j);
						this->resource_grid_DMRS(k1_port(i), l_0 + l_prime(ind_symb), ind_port) = (double)(w_f(ind_port, 1) * w_t(ind_port, ind_symb)) * pseudo_sequence_DMRS(j + 1);
					}
				}
			}
			this->DMRS_position_total = this->DMRS_position_total || this->DMRS_position.slice(ind_port);
		}
	}
	else {	/* Setting DMRS */
		DMRS_length = 2 * n(n.n_elem - 1) + 2;
		this->DMRS_position_total = zeros<umat>(para.num_RB * para.num_sc, para.num_symb);
		for (uword ind_port = 0; ind_port < this->RI; ind_port++) {
			k0_port = k0 + delta(ind_port);
			k1_port = k1 + delta(ind_port);
			for (uword ind_symb = 0; ind_symb < symb_num; ind_symb++) {
				c_init_DMRS = imod((pow(2, 17) * (14 * ind_slot + (l_0 + l_prime(ind_symb)) + 1) * (2 * this->cell_ID + 1) + 2 * this->cell_ID), pow(2, 31));
				func_pseudo_sequence_generation(c_init_DMRS, 2 * DMRS_length, c);
				c1 = zeros<mat>(c.n_rows, c.n_cols / 2);
				c2 = zeros<mat>(c.n_rows, c.n_cols / 2);
				for (uword i = 0, idx = 0; i < c.n_cols; i += 2, idx++) {
					c1.col(idx) = c.col(i);
					c2.col(idx) = c.col(i + 1);
				}
				pseudo_sequence_DMRS = ((1 - 2 * c1) + cx_double(0, 1) * (1 - 2 * c2)) / sqrt(2);
				for (uword i = 0, j = 0; i < k0_port.n_elem; i++, j += 2) {
					this->DMRS_position(k0_port(i), l_0 + l_prime(ind_symb), ind_port) = true;
					this->DMRS_position(k1_port(i), l_0 + l_prime(ind_symb), ind_port) = true;
					this->resource_grid_DMRS(k0_port(i), l_0 + l_prime(ind_symb), ind_port) = (double)(w_f(ind_port, 0) * w_t(ind_port, ind_symb)) * pseudo_sequence_DMRS(j);
					this->resource_grid_DMRS(k1_port(i), l_0 + l_prime(ind_symb), ind_port) = (double)(w_f(ind_port, 1) * w_t(ind_port, ind_symb)) * pseudo_sequence_DMRS(j + 1);
				}
			}
			this->DMRS_position_total = this->DMRS_position_total || this->DMRS_position.slice(ind_port);
		}
	}	/* end of if (para.cal_scenario) */

	this->DMRS_signal = zeros<cx_mat>(sum(sum(this->DMRS_position.slice(0))));
	for (uword ind_port = 0; ind_port < this->RI; ind_port++)
		for (uword j = 0, n = 0; j < this->DMRS_position.n_cols; j++)
			for (uword i = 0; i < this->DMRS_position.n_rows; i++)
				if (this->DMRS_position(i, j, ind_port))
					this->DMRS_signal(n++, ind_port) = this->resource_grid_DMRS(i, j, ind_port);


	/* local memory deallocation */
	w_f.reset();
	w_t.reset();
	delta.reset();
	k_prime.reset();
	l_prime.reset();
	n.reset();
	k0.reset();
	k1.reset();
	k0_port.reset();
	k1_port.reset();
	c.reset();
	c1.reset();
	c2.reset();
	pseudo_sequence_DMRS.reset();

}	/* end of method_DMRS_allocation() method definition */




/*
 * method_data_allocation() = Data allocation
 */
void ModuleBS_SUMIMO::method_data_allocation(ModuleParameterMIMO &para, int &ind_slot) {

	/* local variable declarations */
	umat	Data_position_temp;
	uword	data_position_total_for_LDPC;
	uword	RE_per_PRB;
	double	N_info;
	double	n, C;
	double	N_info_prime;
	uword	TBS_length;
	umat	table_for_LDPC_length;

	if (para.cal_scenario) {
		this->CQI = 1;
		this->RI = 1;
		this->num_layers = this->RI;

		this->CQI_table_num = para.MCS_table_index;
		if (this->CQI_table_num == 1)			/* MCS table 1 */
			for (uword i = 0; i < this->RI; i++)
				this->M_order(i) = CQI_Table_1.order(this->CQI[0] - 1);
		else if (this->CQI_table_num == 2)		/* MCS table 2 (including 256QAM) */
			for (uword i = 0; i < this->RI; i++)
				this->M_order(i) = CQI_Table_1.order(this->CQI[0] - 1);
		else
			perror("MCS table index not supported.");

		this->Data_position = this->Control_position || this->RS_position;
		this->Data_position = matlnot(this->Data_position);

		this->num_codewords = 1;
		para.num_coded_bits(0) = sum(sum(this->Data_position)) * 2;
		para.num_data_bits = 256;
	}
	else {
		this->CQI_table_num = para.MCS_table_index;
		this->M_order = zeros<irowvec>(this->RI);
		this->Coding_rate = zeros<rowvec>(this->RI);
		if (this->CQI_table_num == 1) {			/* MCS table 1 */
			for (uword i = 0; i < this->RI; i++) {
				this->M_order(i) = CQI_Table_1.order(this->CQI[0] - 1);
				this->Coding_rate(i) = CQI_Table_1.code_rate(this->CQI[0] - 1) / 1024;
			}
		}
		else if (this->CQI_table_num == 2) {	/* MCS table 2 (including 256QAM) */
			for (uword i = 0; i < this->RI; i++) {
				this->M_order(i) = CQI_Table_1.order(this->CQI[0] - 1);
				this->Coding_rate(i) = CQI_Table_1.code_rate(this->CQI[0] - 1) / 1024;
			}
		}
		else
			perror("MCS table index not supported.");


		/* set control, RS, and data position */
		this->Control_position = zeros<umat>(arma::size(this->PDCCH_5G_position));
		if (para.use_Sync)
			this->Control_position = this->PBCH_5G_position || this->PDCCH_5G_position || this->PSS_5G_position || this->SSS_5G_position;
		else
			this->Control_position = this->PDCCH_5G_position;

		this->RS_position = this->CSIRS_position_total || this->DMRS_position_total;
		this->Data_position = this->Control_position || this->RS_position;
		this->Data_position = matlnot(this->Data_position);

		/* # of data in a RB */
		this->num_Data_RB = zeros<umat>(para.num_RB, 1);
		Data_position_temp = zeros<umat>(para.num_sc, this->Data_position.n_cols);
		for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
			for (uword i = 0; i < para.num_sc; i++)
				Data_position_temp.row(i) = this->Data_position.row(ind_RB * para.num_sc + i);
			this->num_Data_RB(ind_RB) = sum(sum(Data_position_temp));
		}

		/* Determine the # of codewords */
		this->num_layers = this->RI;
		if (this->num_layers > 4)
			this->num_codewords = 2;
		else
			this->num_codewords = 1;

	}	/* end of if (para.cal_scenario) */

	/* Determine length of generating bits */
	if (para.cal_scenario) {
		para.num_data_bits(0) = 256;
		para.num_coded_bits(0) = 2112;
		this->Coding_rate(0) = 0.1354167;
	}
	else {		/* Setting length of TBS following TS 38.214 5.1.3.2 Transport block size determination */
		for (uword ind_codeword = 0; ind_codeword < para.num_codewords; ind_codeword++) {
			data_position_total_for_LDPC = sum(sum(this->Data_position));
			if (imod(data_position_total_for_LDPC, para.num_RB) != 0)
				perror("PRBs have different configuration");

			RE_per_PRB = floor(data_position_total_for_LDPC / para.num_RB);
			if (RE_per_PRB <= 9)
				RE_per_PRB = 6;
			else if (RE_per_PRB > 9 && RE_per_PRB <= 15)
				RE_per_PRB = 12;
			else if (RE_per_PRB > 15 && RE_per_PRB <= 30)
				RE_per_PRB = 18;
			else if (RE_per_PRB > 30 && RE_per_PRB <= 57)
				RE_per_PRB = 42;
			else if (RE_per_PRB > 57 && RE_per_PRB <= 90)
				RE_per_PRB = 72;
			else if (RE_per_PRB > 90 && RE_per_PRB <= 126)
				RE_per_PRB = 108;
			else if (RE_per_PRB > 126 && RE_per_PRB <= 150)
				RE_per_PRB = 144;
			else
				RE_per_PRB = 156;

			this->N_RE = RE_per_PRB * para.num_RB;
			N_info = (double)this->N_RE * (double)this->num_layers * this->Coding_rate(0) * (double)this->M_order(0);
			if (N_info > 3824.0) {
				n = floor(log2(N_info - 23.0)) - 5.0;
				N_info_prime = pow(2.0, n) * round((N_info - 24.0) / pow(2.0, n));
				if (this->Coding_rate(0) <= 1.0 / 4.0) {
					C = ceil((N_info_prime + 24.0) / 3816.0);
					TBS_length = 8 * (uword) C * (uword) ceil((N_info_prime + 24.0) / (8.0 * C)) - 24;
				}
				else {
					if (N_info_prime > 8424.0) {
						C = ceil((N_info_prime + 24.0) / 8424.0);
						TBS_length = 8 * (uword)C * (uword)ceil((N_info_prime + 24.0) / (8.0 * C)) - 24;
					}
					else
						TBS_length = 8 * (uword)ceil((N_info_prime + 24.0) / 8.0) - 24;
				}
			}
			else {
				n = max(3.0, floor(log2(N_info)) - 6.0);
				N_info_prime = max(24.0, pow(2.0, n) * floor(N_info / pow(2.0, n)));

				table_for_LDPC_length = zeros<umat>(1, 93);
				table_for_LDPC_length << 24 << 32 << 40 << 48 << 56 << 64 << 72 << 80 << 88 << 96
					<< 104 << 112 << 120 << 128 << 136 << 144 << 152 << 160 << 168 << 176
					<< 184 << 192 << 208 << 224 << 240 << 256 << 272 << 288 << 304 << 320
					<< 336 << 352 << 368 << 384 << 408 << 432 << 456 << 480 << 504 << 528
					<< 552 << 576 << 608 << 640 << 672 << 704 << 736 << 768 << 808 << 848
					<< 888 << 928 << 984 << 1032 << 1064 << 1128 << 1160 << 1192 << 1224 << 1256
					<< 1288 << 1320 << 1352 << 1416 << 1480 << 1544 << 1608 << 1672 << 1736 << 1800
					<< 1864 << 1928 << 2024 << 2088 << 2152 << 2216 << 2280 << 2408 << 2472 << 2536
					<< 2600 << 2664 << 2728 << 2792 << 2856 << 2976 << 3104 << 3240 << 3368 << 3496
					<< 3624 << 3752 << 3824 << endr;

				uword i = 0;
				while (N_info_prime > table_for_LDPC_length(0, i))
					i += 1;
				TBS_length = table_for_LDPC_length(i);

			}
			para.num_data_bits = zeros<irowvec>(para.num_codewords);
			para.num_data_bits(ind_codeword) = TBS_length;
		}
	}	/* end of if (para.cal_scenario)  */

	/* local memeory deallocation */
	Data_position_temp.reset();
	table_for_LDPC_length.reset();

}	/* end of method_data_allocation() method definition */




/*
 * method_sync_generation() - Synchronization signal generation
 */
void ModuleBS_SUMIMO::method_sync_generation(ModuleParameterMIMO &para, int &ind_slot) {

	/* local variable declaration */
	bool	empty = true;
	double	N_ID1, N_ID2;
	imat	xx, x0, x1;
	imat	mm;
	int		m0, m1;

	this->cell_ID = para.cell_ID;
	
	uvec Find_Result = find(para.slotset_SS_block == ind_slot + 1);
	if (!Find_Result.is_empty()) {
		N_ID1 = floor(this->cell_ID / 3.0);
		N_ID2 = fmod(this->cell_ID, 3.0);

		/* length of PSS/SSS : 127 */
		this->PSS_signal = zeros<imat>(127, 1);
		this->SSS_signal = zeros<imat>(127, 1);

		/* PSS generation */
		xx  << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr
			<< 0 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr
			<< 1 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr
			<< 0 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr
			<< 0 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr
			<< 0 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr
			<< 0 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr
			<< 1 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr << 1 << endr
			<< 0 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr
			<< 0 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr
			<< 1 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr
			<< 1 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr
			<< 1 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr;


		mm = zeros<imat>(1, 127);
		for (uword i = 0; i < 127; i++)
			mm(0, i) = mod(i + 43 * N_ID2, 127);

		for (uword i = 0; i < 127; i++)
			this->PSS_signal(i, 0) = 1 - 2 * xx(mm(i), 0);


		/* SSS generation */
		m0 = 15 * floor(N_ID1 / 112) + 5 * N_ID2;
		m1 = mod(N_ID1, 112);

		x0  << 0 << endr << 0 << endr << 0 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr
			<< 0 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr << 1 << endr
			<< 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 0 << endr
			<< 1 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr
			<< 1 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr
			<< 0 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr
			<< 1 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr
			<< 0 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr
			<< 0 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr
			<< 0 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr
			<< 0 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr
			<< 1 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr << 1 << endr
			<< 0 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr;

		x1  << 0 << endr << 0 << endr << 0 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr
			<< 0 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr
			<< 1 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr
			<< 1 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr
			<< 1 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr
			<< 1 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 0 << endr
			<< 1 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr
			<< 0 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr
			<< 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr
			<< 0 << endr << 0 << endr << 0 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr
			<< 0 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 1 << endr << 1 << endr << 1 << endr << 0 << endr << 0 << endr
			<< 1 << endr << 1 << endr << 0 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr << 1 << endr << 0 << endr
			<< 1 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr << 1 << endr;

		for (uword i = 0; i < 127; i++)
			this->SSS_signal(i) = (1 - 2 * x0(imod(i + m0, 127))) * (1 - 2 * x1(imod(i + m1, 127)));

	}

	/* memory deallocations */
	xx.reset();
	x0.reset();
	x1.reset();
	mm.reset();

}	/* end of method_sync_generation() method definition */




/*
 * method_HARQ_process() - HARQ process
 */
void ModuleBS_SUMIMO::method_HARQ_process(ModuleParameterMIMO &para, int &ind_subframe, feedback &feedback) {

	/* local variable declaration */
	uword	i;
	uvec	temp;
	uvec	HARQ_empty_process;

	if (this->HARQ_buffer.retrans(this->HARQ_process_index(0)).is_empty())
		this->HARQ_buffer.retrans(this->HARQ_process_index(0)) = zeros<uvec>(this->num_codewords);

	/* feedback.ACK == 1 ..originally */
	if (feedback.ACK(0) == 1 || (this->HARQ_buffer.retrans(this->HARQ_process_index(0))(0) == this->HARQ_buffer.max_retrans)) {
		this->HARQ_buffer.timing(this->HARQ_process_index(0)) = 0;
		this->HARQ_buffer.process(this->HARQ_process_index(0)) = false;
	}

	/* HARQ_buffer.retrans */
	for (uword ind_codeword = 0; ind_codeword < this->num_codewords; ind_codeword++) {
		if (feedback.ACK(ind_codeword) == 1 || feedback.ACK.is_empty())
			this->HARQ_buffer.retrans(this->HARQ_process_index(0))(ind_codeword) = 0;
		else if (this->HARQ_buffer.retrans(this->HARQ_process_index(0))(ind_codeword) < this->HARQ_buffer.max_retrans) {
			this->HARQ_buffer.timing(this->HARQ_process_index(0)) = (ind_subframe - 1) + 8;
			this->HARQ_buffer.retrans(this->HARQ_process_index(0))(ind_codeword) = this->HARQ_buffer.retrans(this->HARQ_process_index(0))(ind_codeword) + 1;
		}
		else
			this->HARQ_buffer.retrans(this->HARQ_process_index(0))(ind_codeword) = 0;
	}

	/* HARQ_buffer.process */
	if (this->HARQ_buffer.retrans(this->HARQ_process_index(0))(0) == 0) {
		this->HARQ_buffer.timing(this->HARQ_process_index(0)) = 0;
		this->HARQ_buffer.process(this->HARQ_process_index(0)) = false;
	}

	temp = find(this->HARQ_buffer.timing == ind_subframe + 1);
	if (temp.is_empty())
		this->HARQ_process_index.reset();
	else
		this->HARQ_process_index(0) = temp(0);

	if (this->HARQ_process_index.is_empty()) {
		HARQ_empty_process = find(this->HARQ_buffer.process == false);
		this->HARQ_process_index = zeros<uvec>(this->num_codewords);
		this->HARQ_process_index(0) = HARQ_empty_process(0);
		this->HARQ_newdata_indicator = zeros<uvec>(this->num_codewords);
		this->HARQ_newdata_indicator(span(0, this->num_codewords - 1)) = true;
		this->HARQ_buffer.process(this->HARQ_process_index(0)) = true;
	}
	else {
		for (uword ind_codeword = 0; ind_codeword < this->num_codewords; ind_codeword++)
			if (this->HARQ_buffer.retrans(this->HARQ_process_index(0))(ind_codeword) == 0)
				this->HARQ_newdata_indicator(ind_codeword) = false;
			else
				this->HARQ_newdata_indicator(ind_codeword) = true;
	}

	/* local memeory deallocation */
	temp.reset();
	HARQ_empty_process.reset();

}	/* end of method_HARQ_process() method definition */




/*
 * method_bit_generation() - Bit generation
 */
void ModuleBS_SUMIMO::method_bit_generation(ModuleParameterMIMO &para, int &ind_slot) {
	
	//cout << "BS.method_bit_generation(): START" << endl;

	/* local variable declaration */
	rowvec	temp_stream;
	irowvec length_bit(para.num_data_bits);

	this->bit_stream.reset();
	this->bit_stream = field<urowvec>(para.num_codewords);

	if (!fixRandom) {	/* Random number enabled */
		for (uword ind_codeword = 0; ind_codeword < para.num_codewords; ind_codeword++) {
			this->bit_stream(ind_codeword) = zeros<urowvec>(length_bit(ind_codeword));
			temp_stream = round(randu<rowvec>(length_bit(ind_codeword)));
			for (uword i = 0; i < temp_stream.n_elem; i++)
				this->bit_stream(ind_codeword)(i) = (unsigned int)temp_stream(i);	// type conversion to int from double

			temp_stream.reset();	/* memory deallocation */
		}
	}
	else {		/* Random number disabled. bit_stream value extracted from MATLAB running code that is random fixed */
		urowvec	fixed_bit_stream;
		switch (length_bit(0)) {
		case 1160:		// bit_stream length = 576
			fixed_bit_stream.load("MATLAB_fixed_bit_stream_1160.csv", csv_ascii);
			break;
		case 1256:		// bit_stream length = 888
			fixed_bit_stream.load("MATLAB_fixed_bit_stream_1256.csv", csv_ascii);
			break;
		default:		// There's no fixed bit stream. So use random bit stream.
			for (uword ind_codeword = 0; ind_codeword < para.num_codewords; ind_codeword++) {
				this->bit_stream(ind_codeword) = zeros<urowvec>(length_bit(ind_codeword));
				temp_stream = round(randu<rowvec>(length_bit(ind_codeword)));
				for (uword i = 0; i < temp_stream.n_elem; i++)
					this->bit_stream(ind_codeword)(i) = (unsigned int)temp_stream(i);	// type conversion to int from double

				temp_stream.reset();	/* memory deallocation */
			}
		}

		for (uword ind_codeword = 0; ind_codeword < para.num_codewords; ind_codeword++)
			this->bit_stream(ind_codeword) = fixed_bit_stream;

		fixed_bit_stream.reset();	/* memory deallocations */
	}

	/* local memory deallocation */
	temp_stream.reset();
	length_bit.reset();

	//cout << "BS.method_bit_generation(): END" << endl;

}	/* end of method_bit_generation() method definition */




/*
 * method_channel_coding() - Channel coding
 */
void ModuleBS_SUMIMO::method_channel_coding(ModuleParameterMIMO &para, NR_ChannelCoding &LDPC, int &ind_slot, feedback &feedback) {

	cout << "BS.method_channel_coding() START" << endl;

	/* local variable definition */
	uword ind_rep;

	this->ch_coding = field<ChCoding>(para.num_codewords);
	this->coded_bit_stream = field<urowvec>(para.num_codewords);
	para.num_coded_bits = zeros<irowvec>(para.num_codewords);

	for (uword ind_codeword = 0; ind_codeword < para.num_codewords; ind_codeword++) {
		if ((int)fix(ind_slot / (10 * pow(2, para.mu))) == 0) {
			if (LDPC_MODULE == LDPC_MATLAB) {
				if (this->N_RE < sum(sum(this->Data_position)))
					LDPC_initialization(para.num_data_bits(ind_codeword), this->M_order(0) * this->N_RE, this->M_order(0), this->Coding_rate(0), para.temp_ch_coding(ind_slot));
				else
					LDPC_initialization(para.num_data_bits(ind_codeword), this->M_order(0) * sum(sum(this->Data_position)), this->M_order(0), this->Coding_rate(0), para.temp_ch_coding(ind_slot));
			}
			else {
				if (this->N_RE < sum(sum(this->Data_position)))
					LDPC_C_initialization(LDPC, para.num_data_bits(ind_codeword), this->N_RE * this->M_order(0), this->M_order(0), this->Coding_rate(0), para.temp_ch_coding(ind_slot));
				else
					LDPC_C_initialization(LDPC, para.num_data_bits(ind_codeword), this->M_order(0) * sum(sum(this->Data_position)), this->M_order(0), this->Coding_rate(0), para.temp_ch_coding(ind_slot));
			}

			this->ch_coding(ind_codeword) = para.temp_ch_coding(ind_slot);
		}
		else {
			if (imod(ind_slot + 1, 10 * pow(2, para.mu)) == 0)
				ind_rep = 10 * pow(2, para.mu) - 1;
			else
				ind_rep = imod(ind_slot, 10 * pow(2, para.mu));

			this->ch_coding(ind_codeword) = para.temp_ch_coding(ind_rep);
		}

		cout << "BS.method_channel_coding(): STEP 00" << endl;

		/* LDPC decoding */
		this->ch_coding(ind_codeword).Decoding_type = 0;		// 0) Sum-product, 1) Min_sum
		this->ch_coding(ind_codeword).Max_LDPC_itr  = 50;
		this->ch_coding(ind_codeword).MS_factor     = 1;		// min sum scaling factor
		this->ch_coding(ind_codeword).LLR_Threshold = 10;

		if (ind_slot == 0) {
			this->ch_coding(ind_codeword).NACK_a = 0;			// NACK for CRC calculation
			this->ch_coding(ind_codeword).NACK_b = zeros<urowvec>(this->ch_coding(ind_codeword).C);	// NACK for code block segmentation
		}
		else {
			this->ch_coding(ind_codeword).NACK_a = feedback.NACK_a;
			this->ch_coding(ind_codeword).NACK_b = feedback.NACK_b;
		}

		this->ch_coding(ind_codeword).nr_codec.harq.HARQ_indicator = this->HARQ_newdata_indicator(ind_codeword);

		this->ch_coding(ind_codeword).A = para.num_data_bits(ind_codeword);
		this->ch_coding(ind_codeword).tx_a = this->bit_stream(ind_codeword);

		cout << "BS.method_channel_coding(): STEP 01" << endl;
		cout << "    this->ch_coding(ind_codeword).A = " << this->ch_coding(ind_codeword).A << endl;
		cout << "    size(this->ch_coding(ind_codeword).tx_a) = " << arma::size(this->ch_coding(ind_codeword).tx_a) << endl;

		if (LDPC_MODULE == LDPC_MATLAB)
			LDPC_tx(this->ch_coding(ind_codeword));
		else
			LDPC_C_tx(LDPC, this->ch_coding(ind_codeword));

		cout << "BS.method_channel_coding(): STEP 02" << endl;

		this->coded_bit_stream(ind_codeword) = this->ch_coding(ind_codeword).tx_g;
		para.num_coded_bits(ind_codeword) = this->coded_bit_stream(ind_codeword).n_elem;

		if (SLOT == 1100) {
			this->ch_coding(ind_codeword).tx_g.print("this->ch_coding(ind_codeword).tx_g:");
			//this->coded_bit_stream(ind_codeword).print("this->coded_bit_stream:");
			exit(1);
		}

		cout << "BS.method_channel_coding(): STEP 03" << endl;
	}

	cout << "BS.method_channel_coding() END" << endl;

}	/* end of method_channel_coding() method definition */




/*
 * method_scrambling() - Scrambling
 *                       TS 38.211, 7.3.1.1  (% TS 36.211, 6.3.1)
 */
void ModuleBS_SUMIMO::method_scrambling(ModuleParameterMIMO &para, int &ind_subframe) {
	
	//cout << "DEBUG: BS.method_scrambling() START" << endl;

	int			c_init;
	field<mat> scrambling_seq = field<mat>(1, this->num_codewords);
	this->scrambled_bit_stream = field<urowvec>(this->num_codewords);

	for (uword ind_codeword = 0; ind_codeword < this->num_codewords; ind_codeword++) {
		c_init = this->n_RNTI * pow(2, 14) + (ind_codeword + 1) * pow(2, 13) + ind_subframe * pow(2, 9), this->cell_ID;
		func_pseudo_sequence_generation(c_init, (double)this->coded_bit_stream(ind_codeword).n_elem, scrambling_seq(ind_codeword));
		this->scrambled_bit_stream(ind_codeword) = zeros<urowvec>(scrambling_seq(ind_codeword).n_elem);
		for (uword i = 0; i < this->scrambled_bit_stream(ind_codeword).n_elem; i++)
			this->scrambled_bit_stream(ind_codeword)(i) = mod(this->coded_bit_stream(ind_codeword)(i) + scrambling_seq(ind_codeword)(i), 2);
	}

	for (uword i = 0; i < scrambling_seq.n_elem; i++)
		if (!(scrambling_seq(0)(i) == 0 || scrambling_seq(0)(i) == 1))
			printf("BS.method_scrabling(): scrambling_seq(0)(%d) = %d\n", i, scrambling_seq(0)(i));

	for (uword i = 0; i < this->coded_bit_stream.n_elem; i++)
		if (!(this->coded_bit_stream(0)(i) == 0 || this->coded_bit_stream(0)(i) == 1))
			printf("BS.method_scrambling(): this->coded_bit_stream(0)(%d) = %d\n", i, this->coded_bit_stream(0)(i));

	for (uword i = 0; i < this->scrambled_bit_stream.n_elem; i++)
		if (!(this->scrambled_bit_stream(0)(i) == 0 || this->scrambled_bit_stream(0)(i) == 1))
			printf("BS.method_scrambling(): this->scrambled_bit_stream(0)(%d) = %d\n", i, this->scrambled_bit_stream(0)(i));

	if (SLOT == 1000) {
		scrambling_seq(0).print("scrambling_seq(0):");
		this->coded_bit_stream(0).print("this->coded_bit_stream:");	// 여기에서 에러가 파생됨.
		this->scrambled_bit_stream(0).print("BS.method_scrampling(): scrambled_bit_stream:");	// 여기서도 잘못됨.
		exit(1);
	}

	/* local memory deallocation */
	for (uword i = 0; i < scrambling_seq.n_elem; i++)
		scrambling_seq(i).reset();
	scrambling_seq.reset();

	//cout << "DEBUG: BS.method_scrambling() END" << endl;

}	/* end of method_scrambling() method definition */




/*
 * method_modulation_mapping() - Modulation mapping
 *								 CQI table: TS 38.214 V15.0.0, 5.2.2.1 (36.213 Table 7.2.3-1)
 *								 modulation : 38.211 7.3.1.2 and 5.1 (36.211 Table 7.1.1-1~7.1.5-1
 */
void ModuleBS_SUMIMO::method_modulation_mapping(void) {

	cout << "BS.method_modulation_mapping(): START" << endl;

	int order = this->M_order(0);

	cout << "    order = " << order << endl;

	this->modulated_signal_stream = field<cx_rowvec>(this->bit_stream.n_rows);
	umat temp_scrambled_bit_stream;

	for (int ind_codeword = 0; ind_codeword < this->bit_stream.n_rows; ind_codeword++) {
		this->modulated_signal_stream(ind_codeword) = zeros<cx_rowvec>(this->scrambled_bit_stream(ind_codeword).n_elem / order);
		temp_scrambled_bit_stream = reshape(this->scrambled_bit_stream(ind_codeword), order, this->modulated_signal_stream(ind_codeword).n_cols);

		//temp_scrambled_bit_stream = zeros<umat>(order, this->modulated_signal_stream(ind_codeword).n_cols);
		//for (uword i = 0; i < this->scrambled_bit_stream(ind_codeword).n_elem; i++)
		//	temp_scrambled_bit_stream(i) = this->scrambled_bit_stream(ind_codeword)(i);
		cout << "    size(this->modulated_signal_stream(" << ind_codeword << ") = " << arma::size(this->modulated_signal_stream(ind_codeword)) << endl;
		cout << "    size(temp_scrambled_bit_stream) =" << arma::size(temp_scrambled_bit_stream) << endl;
		cout << "    size(modulation_mapping) = " << arma::size(modulation_mapping) << ", size(modulation_mapping(order - 1)) = " << arma::size(modulation_mapping(order - 1)) << endl;
		if (SLOT > 10 && LDPC_MODULE == LDPC_CPP) {
			//temp_scrambled_bit_stream.print("temp_scrambled_bit_stream:");
			for (uword i = 0, j = 0; i < temp_scrambled_bit_stream.n_elem; i++)
				if (!(temp_scrambled_bit_stream(i) == 0 || temp_scrambled_bit_stream(i) == 1)) {
					printf("temp_scrambled_bit_stream(%d) = %u", i, temp_scrambled_bit_stream(i));
					printf(" == this->scrambled_bit_stream(%d)(%d) = %u(%d)\n", ind_codeword, i, this->scrambled_bit_stream(ind_codeword)(i), j++);
				}

			for (uword i = 0, j = 0; i < this->scrambled_bit_stream(0).n_elem; i++)
				if (!(this->scrambled_bit_stream(ind_codeword)(i) == 0 || this->scrambled_bit_stream(ind_codeword)(i) == 1))
					printf("    this->scrambled_bit_stream(%d)(%d) = %u[%u]\n", ind_codeword, i, this->scrambled_bit_stream(ind_codeword)(i), j++);
		}

		switch (order) {
		case 2:		/* QPSK */
			for (uword i = 0; i < temp_scrambled_bit_stream.n_cols; i++) {
				if (SLOT > 10 && LDPC_MODULE == LDPC_CPP) {
					cout << "modulated_signal_stream(" << ind_codeword << ")(" << i << ") = modulation_mapping(" << order - 1 << ")(" << 2 * temp_scrambled_bit_stream(1, i) + temp_scrambled_bit_stream(0, i) << ")";
					cout << ", temp_scrambled_bit_stream(1, " << i << ") = " << temp_scrambled_bit_stream(1, i) << ", temp_scrambled_bit_stream(0, " << i << ") = " << temp_scrambled_bit_stream(0, i) << endl;
				}
				this->modulated_signal_stream(ind_codeword)(i) = modulation_mapping(order - 1)(2 * temp_scrambled_bit_stream(1, i) + temp_scrambled_bit_stream(0, i));
			}
			break;
		case 4:		/* 16QAM */
			for (uword i = 0; i < temp_scrambled_bit_stream.n_cols; i++)
				this->modulated_signal_stream(ind_codeword)(i) = modulation_mapping(order - 1)(8 * temp_scrambled_bit_stream(3, i) + 4 * temp_scrambled_bit_stream(2, i) + 2 * temp_scrambled_bit_stream(1, i) + temp_scrambled_bit_stream(0, i));
			break;
		case 6:		/* 64QAM */
			for (uword i = 0; i < temp_scrambled_bit_stream.n_cols; i++)
				this->modulated_signal_stream(ind_codeword)(i) = modulation_mapping(order - 1)(32 * temp_scrambled_bit_stream(5, i) + 16 * temp_scrambled_bit_stream(4, i) + 8 * temp_scrambled_bit_stream(3, i) + 4 * temp_scrambled_bit_stream(2, i) + 2 * temp_scrambled_bit_stream(1, i) + temp_scrambled_bit_stream(0, i));
			break;
		case 8:		/* 256QAM */
			for (uword i = 0; i < temp_scrambled_bit_stream.n_cols; i++)
				this->modulated_signal_stream(ind_codeword)(i) = modulation_mapping(order - 1)(128 * temp_scrambled_bit_stream(7, i) + 64 * temp_scrambled_bit_stream(6, i) + 32 * temp_scrambled_bit_stream(5, i) + 16 * temp_scrambled_bit_stream(4, i) + 8 * temp_scrambled_bit_stream(3, i) + 4 * temp_scrambled_bit_stream(2, i) + 2 * temp_scrambled_bit_stream(1, i) + temp_scrambled_bit_stream(0, i));
			break;
		default:
			;
		}
	}

	temp_scrambled_bit_stream.reset();	/* memory deallocation */

	cout << "BS.method_modulation_mapping(): END" << endl;

}	/* end of method_modulation_mapping() method definition */




/*
 * method_layer_mapping() - Layer mapping
 *           Layer mapping, TS. 36.211, Table 6.3.3.2-1, V13.3.0.
 */
void ModuleBS_SUMIMO::method_layer_mapping(ModuleParameterMIMO &para) {

	/* local variable declaration */
	cx_mat		layered_PDCCH_stream_temp;
	cx_mat		layered_PBCH_stream_temp;
	cx_rowvec	modulated_signal_stream_for_grid;
	cx_mat		resized_modulated_signal_stream_for_grid;
	uvec		codeword1, codeword2;
	cx_rowvec	modulated_signal_stream_for_grid_1, modulated_signal_stream_for_grid_2;
	cx_rowvec	temp_modulated_signal_stream_for_grid_1, temp_modulated_signal_stream_for_grid_2;
	uword		add_elem;
	uword		layered_signal_len_temp;

	/* layered_PDCCH_stream */
	if (para.use_PDCCH) {
		layered_PDCCH_stream_temp = reshape(this->resource_grid(find(this->PDCCH_5G_position)), this->resource_grid.n_elem / min(4, para.num_Tx_antenna), min(4, para.num_Tx_antenna)).t();
		layered_PDCCH_stream_temp = repmat(layered_PDCCH_stream_temp, 2, 1);
		this->layered_PDCCH_stream.rows(0, min(4, para.num_port) - 1) = real(layered_PDCCH_stream_temp.rows(0, min(4, para.num_port) - 1));
		this->layered_PDCCH_stream.rows(min(4, para.num_port), this->layered_PDCCH_stream.n_rows - 1) = imag(layered_PDCCH_stream.rows(min(4, para.num_port), this->layered_PDCCH_stream.n_rows - 1));
	}

	/* layered_PBCH_stream */
	if (para.use_PBCH) {
		layered_PBCH_stream_temp = reshape(this->resource_grid(find(this->PBCH_5G_position)), this->resource_grid.n_elem / min(4, para.num_Tx_antenna), min(4, para.num_Tx_antenna)).t();
		layered_PBCH_stream_temp = repmat(layered_PBCH_stream_temp, 2, 1);
		this->layered_PBCH_stream.rows(0, min(4, para.num_port) - 1) = real(layered_PBCH_stream_temp.rows(0, min(4, para.num_port) - 1));
		this->layered_PBCH_stream.rows(min(4, para.num_port), this->layered_PBCH_stream.n_rows - 1) = imag(layered_PBCH_stream.rows(min(4, para.num_port), this->layered_PDCCH_stream.n_rows - 1));
	}

	/* layered signal */
	uword add_size = 0;
	
	switch (this->num_codewords) {
	case 1:
		modulated_signal_stream_for_grid = this->modulated_signal_stream(0);
		while (mod(modulated_signal_stream_for_grid.n_elem, this->num_layers) != 0)
			add_size++;
		resized_modulated_signal_stream_for_grid = zeros<cx_mat>(1, modulated_signal_stream_for_grid.n_elem + add_size);
		resized_modulated_signal_stream_for_grid.cols(0, modulated_signal_stream_for_grid.n_elem - 1) = modulated_signal_stream_for_grid;
		if (add_size > 0)
			resized_modulated_signal_stream_for_grid.cols(modulated_signal_stream_for_grid.n_elem, modulated_signal_stream_for_grid.n_elem + add_size) = 0;
		this->layered_signal = reshape(resized_modulated_signal_stream_for_grid, this->num_layers, resized_modulated_signal_stream_for_grid.n_elem / this->num_layers);
		break;
	case 2:
		codeword1 = zeros<uvec>(floor(this->num_layers / 2));
		for (uword i = 0; i < codeword1.n_elem; i++)
			codeword1(i) = i;
		codeword2 = zeros<uvec>(this->num_layers / 2);
		for (uword i = 0; i < codeword2.n_elem; i++)
			codeword2(i) = i + floor(this->num_layers / 2);
		add_elem = 0;
		while (mod(this->modulated_signal_stream(0).n_elem, codeword1.n_elem) != 0)
			add_elem++;
		temp_modulated_signal_stream_for_grid_1 = zeros<cx_rowvec>(this->modulated_signal_stream(0).n_elem + add_elem);
		temp_modulated_signal_stream_for_grid_1(span(0, this->modulated_signal_stream(0).n_elem - 1)) = this->modulated_signal_stream(0);
		temp_modulated_signal_stream_for_grid_1(span(this->modulated_signal_stream(0).n_elem, temp_modulated_signal_stream_for_grid_1.n_elem - 1)) = 0;
		add_elem = 0;
		while (mod(this->modulated_signal_stream(2).n_elem, codeword2.n_elem) != 0)
			add_elem++;
		temp_modulated_signal_stream_for_grid_2 = zeros<cx_rowvec>(this->modulated_signal_stream(1).n_elem + add_elem);
		temp_modulated_signal_stream_for_grid_2(span(0, this->modulated_signal_stream(1).n_elem - 1)) = this->modulated_signal_stream(1);
		temp_modulated_signal_stream_for_grid_2(span(this->modulated_signal_stream(1).n_elem, temp_modulated_signal_stream_for_grid_2.n_elem - 1)) = 0;

		if (temp_modulated_signal_stream_for_grid_1.n_elem > temp_modulated_signal_stream_for_grid_2.n_elem) {
			add_elem = 0;
			while (temp_modulated_signal_stream_for_grid_1.n_elem != temp_modulated_signal_stream_for_grid_2.n_elem)
				add_elem++;
			modulated_signal_stream_for_grid_2 = zeros<cx_rowvec>(temp_modulated_signal_stream_for_grid_2.n_elem + add_elem);
			modulated_signal_stream_for_grid_2(span(0, temp_modulated_signal_stream_for_grid_2.n_elem - 1)) = temp_modulated_signal_stream_for_grid_2;
			modulated_signal_stream_for_grid_2(span(temp_modulated_signal_stream_for_grid_2.n_elem, modulated_signal_stream_for_grid_2.n_elem - 1)) = 0;
		}
		else {
			add_elem = 0;
			while (temp_modulated_signal_stream_for_grid_1.n_elem != temp_modulated_signal_stream_for_grid_2.n_elem)
				add_elem++;
			modulated_signal_stream_for_grid_1 = zeros<cx_rowvec>(temp_modulated_signal_stream_for_grid_1.n_elem + add_elem);
			modulated_signal_stream_for_grid_1(span(0, temp_modulated_signal_stream_for_grid_1.n_elem - 1)) = temp_modulated_signal_stream_for_grid_1;
			modulated_signal_stream_for_grid_1(span(temp_modulated_signal_stream_for_grid_1.n_elem, modulated_signal_stream_for_grid_1.n_elem - 1)) = 0;
		}
		
		layered_signal_len_temp = max(modulated_signal_stream_for_grid_1.n_elem / codeword1.n_elem, modulated_signal_stream_for_grid_2.n_elem / codeword2.n_elem);
		this->layered_signal = zeros<cx_mat>(this->num_layers, layered_signal_len_temp);
		this->layered_signal.rows(codeword1) = reshape(modulated_signal_stream_for_grid_1, codeword1.n_elem, modulated_signal_stream_for_grid_1.n_elem / codeword1.n_elem);
		this->layered_signal.rows(codeword2) = reshape(modulated_signal_stream_for_grid_2, codeword2.n_elem, modulated_signal_stream_for_grid_2.n_elem / codeword2.n_elem);
	}

	/* local memory deallocation */
	layered_PDCCH_stream_temp.reset();
	layered_PBCH_stream_temp.reset();
	modulated_signal_stream_for_grid.reset();
	resized_modulated_signal_stream_for_grid.reset();
	modulated_signal_stream_for_grid_1.reset();
	modulated_signal_stream_for_grid_2.reset();
	temp_modulated_signal_stream_for_grid_1.reset();
	temp_modulated_signal_stream_for_grid_2.reset();

}	/* end of method_layer_mapping() method definition */




/*
 * method_precoding() - Precoding
 */
void ModuleBS_SUMIMO::method_precoding(ModuleParameterMIMO &para) {

	/* local variable definition */
	uword		length_temp;
	uvec		pos, len_data;
	umat		Data_position_temp;
	uvec		num_data_RB;
	uword		from, to;
	uword		Data_index;
	cx_mat		Data_temp;
	cx_mat      temp_PDCCH, temp_PBCH;


	/* Signal power allocation */
	if (para.use_Power_allocation) {
		this->power_DMRS  = para.power_DMRS;
		this->power_CSIRS = para.power_CSIRS;
		this->power_data  = para.power_data;
	}
	else {
		this->power_DMRS  = zeros<rowvec>(1);
		this->power_CSIRS = zeros<rowvec>(1);
		this->power_data  = zeros<rowvec>(1);
		this->power_DMRS(0)  = 1.0;
		this->power_CSIRS(0) = 1.0;
		this->power_data(0)  = 1.0;
	}

	/* Resource grid of DMRS and CSI-RS */
	this->resource_grid_DMRS = this->power_DMRS(0) * this->resource_grid_DMRS;
	this->resource_grid_CSIRS = this->power_CSIRS(0) * this->resource_grid_CSIRS;

	/* LDPC data fitting */
	for (uword ind_codeword = 0; ind_codeword < para.num_codewords; ind_codeword++) {
		length_temp = sum(sum(this->Data_position));
		if (this->layered_signal.n_elem <= length_temp) {
			pos = find(this->Data_position);
			len_data = pos(span(0, this->layered_signal.n_elem - 1));
			this->Data_position.zeros();
			for (uword i = 0; i < len_data.n_elem; i++)
				this->Data_position(len_data(i)) = true;
		}
		else
			perror("The number of generated data is not proper");
	}

	/* Calculate num_Data_RB */
	for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
		Data_position_temp = this->Data_position(span(0 + ind_RB * para.num_sc, para.num_sc - 1 + ind_RB * para.num_sc), span(0, para.num_symb - 1));
		this->num_Data_RB(ind_RB) = sum(sum(Data_position_temp));
	}

	/* Data precoding and antenna port mapping */
	this->precoded_signal_stream = zeros<cx_mat>(para.num_Tx_antenna, this->layered_signal.n_elem);
	num_data_RB = this->num_Data_RB;
	for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
		if (ind_RB == 0)
			from = 0;
		else
			from = sum(num_data_RB.rows(0, ind_RB - 1));
		to = sum(num_data_RB.rows(0, ind_RB)) - 1;

		this->precoded_signal_stream.cols(from, to) = this->Precoder(ind_RB).col(0) * this->power_data(0) * this->layered_signal.cols(from, to);
	}

	/* Resource grid of data */
	this->resource_grid = zeros<cx_cube>(para.num_RB * para.num_sc, para.num_symb, para.num_Tx_antenna);
	this->precoded_signal_stream = this->precoded_signal_stream.st();
	for (uword ind_port = 0; ind_port < para.num_Tx_antenna; ind_port++) {
		Data_index = 0;
		for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
			Data_temp = this->resource_grid(span(ind_RB * para.num_sc, (ind_RB + 1) * para.num_sc - 1), span(0, para.num_symb - 1), span(ind_port));
			Data_position_temp = this->Data_position(span(ind_RB * para.num_sc, (ind_RB + 1) * para.num_sc - 1), span(0, para.num_symb - 1));
			Data_temp(find(Data_position_temp)) = this->precoded_signal_stream(span(Data_index, Data_index + this->num_Data_RB(ind_RB) - 1), span(ind_port));
			Data_index += this->num_Data_RB(ind_RB);
			this->resource_grid(span(ind_RB * para.num_sc, (ind_RB + 1) * para.num_sc - 1), span(0, para.num_symb - 1), span(ind_port)) = Data_temp;
		}
	}

	/* PDCCH & PBCH precoding and antenna port mapping */
	if (para.use_PDCCH) {
		temp_PDCCH = this->Precoder_PDCCH * this->layered_PBCH_stream;
		this->precoded_PDCCH_stream = reshape(temp_PDCCH, 1, temp_PDCCH.n_elem);
		this->resource_grid(find(this->PDCCH_5G_position)) = this->precoded_PDCCH_stream;
	}

	if (para.use_PBCH) {
		temp_PBCH = this->Precoder_PDCCH * this->layered_PBCH_stream;
		this->precoded_PBCH_stream = reshape(temp_PBCH, 1, temp_PBCH.n_elem);
		this->resource_grid(find(this->PBCH_5G_position)) = this->precoded_PBCH_stream;
	}

	/* CSI-RS precoding */
	if (!para.cal_scenario) {
		for (uword ind_Tx = 0; ind_Tx < para.num_Tx_antenna; ind_Tx++)
			for (uword ind_port_CSIRS = 0; ind_port_CSIRS < para.num_port_CSIRS; ind_port_CSIRS++)
				this->resource_grid.slice(ind_Tx) = this->resource_grid.slice(ind_Tx) + this->Precoder_CSIRS(ind_Tx, ind_port_CSIRS) * this->resource_grid_CSIRS.slice(ind_port_CSIRS);
	}

	/* DMRS precoding */
	for (uword ind_Tx = 0; ind_Tx < para.num_Tx_antenna; ind_Tx++)
		for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++)
			for (uword ind_layer = 0; ind_layer < this->RI; ind_layer++)
				this->resource_grid(span(0 + ind_RB * para.num_sc, para.num_sc - 1 + ind_RB * para.num_sc), span::all, span(ind_Tx)) += this->Precoder(ind_RB)(ind_Tx, ind_layer) * this->resource_grid_DMRS(span(0 + ind_RB * para.num_sc, para.num_sc - 1 + ind_RB * para.num_sc), span::all, span(ind_layer));

	this->data_position_total = this->Data_position;


	/* local memory deallocation */
	pos.reset();
	len_data.reset();
	Data_position_temp.reset();
	num_data_RB.reset();
	Data_temp.reset();
	temp_PDCCH.reset();
	temp_PBCH.reset();

}	/* end of method_precoding() method definition */




/*
 * method_OFDM_generation() - OFDM generation
 *					TS 38.211 5.3
 *					For given mapped signal streams, generate OFDM signals
 */
void ModuleBS_SUMIMO::method_OFDM_generation(ModuleParameterMIMO &para) {

	/* local variable declaration */
	cx_mat		resource_grid_temp, resource_grid_2tmp, sq_resource_grid_1, sq_resource_grid_2;
	cx_mat		resource_grid_ifft;
	cx_mat		signal_stream_block_1, signal_stream_block_2, signal_stream_block_3, signal_stream_block_4, block_2_temp, block_4_temp;
	uword		size_block_1, size_block_2, size_block_3, size_block_4;

	/* calculate the size of this->OFDM_signal_stream and allocate memory to this->OFDM_signal_stream */
	size_block_1 = para.size_fft + para.length_CP(0);
	size_block_2 = (para.size_fft + para.length_CP(1)) * para.num_symb / 2;
	size_block_3 = size_block_1;
	size_block_4 = (para.size_fft + para.length_CP(1)) * (para.num_symb / 2 - 2);
	this->OFDM_signal_stream = zeros<cx_mat>(size_block_1 + size_block_2 + size_block_3 + size_block_4, para.num_Tx_antenna);

	this->mu = para.mu;
	for (uword ind_port = 0; ind_port < para.num_Tx_antenna; ind_port++) {
		sq_resource_grid_1 = this->resource_grid(span(para.num_RB * para.num_sc / 2, this->resource_grid.n_rows - 1), span::all, span(ind_port));
		resource_grid_2tmp = join_vert(sq_resource_grid_1, zeros<cx_mat>((para.size_fft - para.num_sc * para.num_RB), para.num_symb));
		sq_resource_grid_2 = this->resource_grid(span(0, para.num_RB * para.num_sc / 2 - 1), span::all, span(ind_port));
		resource_grid_temp = join_vert(resource_grid_2tmp, sq_resource_grid_2);

		resource_grid_ifft = para.size_fft / sqrt(para.num_RB * para.num_sc) * ifft(resource_grid_temp, para.size_fft);

		signal_stream_block_1 = join_vert(resource_grid_ifft(span(para.size_fft - para.length_CP(0), para.size_fft - 1), 0), resource_grid_ifft(span(0, para.size_fft - 1), 0));
		block_2_temp = join_vert(resource_grid_ifft(span(para.size_fft - para.length_CP(1), para.size_fft - 1), span(1, para.num_symb / 2 - 1)), resource_grid_ifft(span(0, para.size_fft - 1), span(1, para.num_symb / 2 - 1)));
		signal_stream_block_2 = reshape(block_2_temp, block_2_temp.n_elem, 1);
		signal_stream_block_3 = join_vert(resource_grid_ifft(span(para.size_fft - para.length_CP(0), para.size_fft - 1), para.num_symb / 2), resource_grid_ifft(span(0, para.size_fft -1), para.num_symb / 2));
		block_4_temp = join_vert(resource_grid_ifft(span(para.size_fft - para.length_CP(1), para.size_fft - 1), span(para.num_symb / 2 + 1, para.num_symb - 1)), resource_grid_ifft(span(0, para.size_fft - 1), span(para.num_symb / 2 + 1, para.num_symb - 1)));
		signal_stream_block_4 = reshape(block_4_temp, block_4_temp.n_elem, 1);
		this->OFDM_signal_stream.col(ind_port) = join_vert(join_vert(signal_stream_block_1, signal_stream_block_2), join_vert(signal_stream_block_3, signal_stream_block_4));
	}

	/* local memory deallocation */
	resource_grid_temp.reset();
	resource_grid_2tmp.reset();
	sq_resource_grid_1.reset();
	sq_resource_grid_2.reset();
	resource_grid_ifft.reset();
	signal_stream_block_1.reset();
	signal_stream_block_2.reset();
	signal_stream_block_3.reset();
	signal_stream_block_4.reset();
	block_2_temp, block_4_temp.reset();

}	/* end of method_OFDM_generation() method definition */




/*
 * method_genie() - Genie
 */
void ModuleBS_SUMIMO::method_genie(void) {

	this->genie.num_codewords = this->num_codewords;
	this->genie.PMI = this->PMI;
	this->genie.RI = this->RI;
	this->genie.CQI_table_num = this->CQI_table_num;
	this->genie.CQI = this->CQI;
	this->genie.m_order = this->M_order;

	/* position */
	this->genie.PDCCH_position = this->PDCCH_5G_position;
	this->genie.PSS_position = this->PSS_position;
	this->genie.SSS_position = this->SSS_position;
	this->genie.PCFICH_position = this->PCFICH_position;
	this->genie.PBCH_position = this->PBCH_5G_position;
	this->genie.DMRS_position = this->DMRS_position;
	this->genie.CSIRS_position = this->CSIRS_position;
	this->genie.CSIRS_position_total = this->CSIRS_position_total;
	this->genie.Data_position = this->Data_position;

	/* signal */
	this->genie.DMRS_signal = this->DMRS_signal;
	this->genie.CSIRS_signal = this->CSIRS_signal;

	/* resource grid */
	this->genie.resource_grid = this->resource_grid;
	this->genie.resource_grid_DMRS = this->resource_grid_DMRS;

	/* HARQ and channel coding */
	this->genie.HARQ_process_index = this->HARQ_process_index;
	this->genie.newdata_indicator = this->HARQ_newdata_indicator;
	this->genie.ch_coding = this->ch_coding;

	/* power allocation */
	this->genie.power_CSIRS = this->power_CSIRS;
	this->genie.power_data = this->power_data;
	this->genie.power_DMRS = this->power_DMRS;

	/* bit stream */
	this->genie.bit_stream = this->bit_stream;
	this->genie.coded_bit_stream = this->coded_bit_stream;
	this->genie.scrambled_bit_stream = this->scrambled_bit_stream;
	this->genie.OFDM_signal_stream = this->OFDM_signal_stream;

	/* MIMO */
	this->genie.num_Data_RB = this->num_Data_RB;
	this->genie.num_layers = this->num_layers;
	this->genie.Precoder = this->Precoder;

	this->genie.selected = true;
	
}	/* end of method_genie() method definition */




/* default constructor & destructor */
ModuleBS_SUMIMO::ModuleBS_SUMIMO()
{
}

ModuleBS_SUMIMO::~ModuleBS_SUMIMO()
{
}





/*
 * local function implementations
 *      - function name: func_CDM(), func_pseudo_sequence_generation()
 *		- MATLAB script file name: func_CDM.m, func_pseudo_sequence_generation.m
 */


/*
 * func_wide_beam_gen() - wide beam generation function
 */
void func_wide_beam_gen(ModuleParameterMIMO &para, field<cx_mat> &Widebeam, field<cx_mat> &Narrowbeam) {

	/* local variable declaration */
	rowvec		N1O1Temp, N1Temp, N2O2Temp, N2Temp;
	mat			Tx_theta_hor, Tx_theta_ver;
	cx_mat		TxB_hor, TxB_ver;
	cx_mat		Wood_TxB_hor, Wood_TxB_ver;
	uvec		beam_pos;
	umat		beam_pos_hor, beam_pos_ver;
	cx_mat		temp;
	cx_mat		polar1, polar2;
	cx_mat		temp_hor, temp_ver;


	int L  = para.num_SS_block;
	int N1 = para.N1;				// # of horizontal Tx antennas
	int N2 = para.N2;				// # of vertical Tx antennas
	int L1 = para.L1;				// # of horizontal wide beams
	int L2 = para.L2;				// # of vertical wide beams
	int O1 = para.O1;				// # of DFT oversampling of horizontal axis
	int O2 = para.O1;				// # of DFT oversampling of vertical axis
	int P1 = para.P1;				// # of horizontal narrow beams per a wide beam
	int P2 = para.P2;				// # of vertical narrow beams per a wide beam
	int S1 = para.S1;				// Narrow beam interval that composed of different wide beam (horizontal)
	int S2 = para.S2;				// Narrow beam interval that composed of different wide beam (vertical)

	Widebeam   = field<cx_mat>(L, 1);
	Narrowbeam = field<cx_mat>(L, 1);

	switch (N2) {
	case 1:
		N1O1Temp = zeros<rowvec>(N1 * O1);
		for (uword i = 0; i < N1 * O1; i++)
			N1O1Temp(i) = (double)i;
		N1Temp = zeros<rowvec>(N1);
		for (uword i = 0; i < N1; i++)
			N1Temp(i) = (double) i;

		Tx_theta_hor = ((N1O1Temp * O1 - 1.0) / N1 / O1).t() * N1Temp;
		TxB_hor = (exp(cx_double(0, -1) * datum::pi * Tx_theta_hor) / sqrt(N1)).t();
		Wood_TxB_hor = TxB_hor * diagmat(exp((cx_double(0, 1) * datum::pi * ((double)N1 - 1.0) / (double)N1 / (double)O1) / 2.0 * N1O1Temp));
		beam_pos = join_horiz(ones<urowvec>(P1), zeros<urowvec>(N1 * O1 - P1));
		for (uword ind_L1 = 0; ind_L1 < L1; ind_L1++) {
			temp = sum(Wood_TxB_hor(find(shift(beam_pos, S1 * ind_L1, 1))), 2);
			temp = temp / norm(temp, 2);
			Widebeam(ind_L1) = join_vert(temp, zeros<cx_mat>(N1 * N2, 1));
			polar1 = join_vert(TxB_hor.cols(find(shift(beam_pos, S1 * ind_L1, 1))), zeros<cx_mat>(N1 * N2, P1));
			polar2 = join_vert(zeros<cx_mat>(N1 * N2), TxB_hor.cols(find(shift(beam_pos, S1 * ind_L1, 1))));
			Narrowbeam(ind_L1) = join_horiz(polar1, polar2);
		}
		break;
	default:
		N1O1Temp = zeros<rowvec>(N1 * O1);
		for (uword i = 0; i <= N1 * O1 - 1; i++)
			N1O1Temp(i) = (double)i;
		N1Temp = zeros<rowvec>(N1);
		for (uword i = 0; i < N1; i++)
			N1Temp(i) = (double)i;
		Tx_theta_hor = (N1O1Temp / N1 / O1).st() * N1Temp;
		
		N2O2Temp = zeros<rowvec>(N2 * O2);
		for (uword i = 0; i <= N2 * O2 - 1; i++)
			N2O2Temp(i) = (double)i;
		N2Temp = zeros<rowvec>(N2);
		for (uword i = 0; i < N2; i++)
			N2Temp(i) = (double)i;
		Tx_theta_ver = (N2O2Temp / N2 / O2).st() * N2Temp;
		
		TxB_hor = (exp(cx_double(0, -1) * datum::pi * Tx_theta_hor) / sqrt(N1)).t();
		TxB_ver = (exp(cx_double(0, -1) * datum::pi * Tx_theta_ver) / sqrt(N2)).t();
		Wood_TxB_hor = TxB_hor * diagmat(exp((cx_double(0, 1) * datum::pi * ((double)N1 - 1.0) / (double)N1 / (double)O1) / 2.0 * N1O1Temp));
		Wood_TxB_ver = TxB_ver * diagmat(exp((cx_double(0, 1) * datum::pi * ((double)N2 - 1.0) / (double)N2 / (double)O2) / 2.0 * N2O2Temp));
		beam_pos_hor = join_horiz(ones<umat>(1, P1), zeros<umat>(1, N1 * O1 - P1));
		beam_pos_ver = join_horiz(ones<umat>(1, P2), zeros<umat>(1, N2 * O2 - P2));

		for (uword ind_L1 = 0; ind_L1 < L1; ind_L1++)
			for (uword ind_L2 = 0; ind_L2 < L2; ind_L2++) {
				temp_hor = sum(Wood_TxB_hor.cols(find(shift(beam_pos_hor, S1 * ind_L1, 1))), 1);
				temp_ver = sum(Wood_TxB_ver.cols(find(shift(beam_pos_ver, S2 * ind_L2, 1))), 1);
				temp = kron(temp_hor, temp_ver);
				temp = temp / norm(temp, 2);
				Widebeam(ind_L1 * L2 + ind_L2) = join_vert(temp, zeros<cx_mat>(N1 * N2, 1));
				polar1 = join_vert(kron(TxB_hor.cols(find(shift(beam_pos_hor, S1 * ind_L1, 1))), TxB_ver.cols(find(shift(beam_pos_ver, S2 * ind_L2, 1)))), zeros<cx_mat>(N1 * N2, P1 * P2));
				polar2 = join_vert(zeros<cx_mat>(N1 * N2, P1 * P2), kron(TxB_hor.cols(find(shift(beam_pos_hor, S1 * ind_L1, 1))), TxB_ver.cols(find(shift(beam_pos_ver, S2 * ind_L2, 1)))));
				Narrowbeam(ind_L1 * L2 + ind_L2) = join_horiz(polar1, polar2);
			}

		break;
	}

	uvec beamIDX;
	if (para.Tx_pol == 1) {
		for (uword ind_wb = 0; ind_wb < Widebeam.n_elem; ind_wb++) {
			for (uword i = 0; i < Widebeam(ind_wb).n_elem / 2; i++)
				beamIDX(i) = i;
			Widebeam(ind_wb) = Widebeam(ind_wb)(beamIDX);
		}
		for (uword ind_nb = 0; ind_nb < Narrowbeam.n_elem; ind_nb++) {
			for (uword i = 0; i < Widebeam(ind_nb).n_elem / 2; i++)
				beamIDX(i) = i;
			Widebeam(ind_nb) = Widebeam(ind_nb)(beamIDX);
		}
	}

	/* local memory dealocation */
	N1O1Temp.reset();
	N1Temp.reset();
	N2O2Temp.reset();
	N2Temp.reset();
	Tx_theta_hor.reset();
	Tx_theta_ver.reset();
	TxB_hor.reset();
	TxB_ver.reset();
	Wood_TxB_hor.reset();
	Wood_TxB_ver.reset();
	beam_pos.reset();
	beam_pos_hor.reset();
	beam_pos_ver.reset();
	temp.reset();
	polar1.reset();
	polar2.reset();
	temp_hor.reset();
	temp_ver.reset();
	beamIDX.reset();

}	/* end of func_wide_beam_gen() */




/*
 * func_CDM()
 */
void func_CDM(CSIRS_CDMType CDM_type, imat &w_f, imat &w_t, int &length_f, int &length_t) {

	switch (CDM_type) {
	case CSIRS_CDMType::NO_CDM:
		w_f << 1 << endr;
		w_t << 1 << endr;
		length_f = 1;
		length_t = 1;
		break;
	case CSIRS_CDMType::FD_CDM2:
		w_f << 1 <<  1 << endr
			<< 1 << -1 << endr;
		w_t << 1 << endr
			<< 1 << endr;
		length_f = 2;
		length_t = 1;
		break;
	case CSIRS_CDMType::CDM4:
		w_f << 1 <<  1 << endr
			<< 1 << -1 << endr
			<< 1 <<  1 << endr
			<< 1 << -1 << endr;
		w_t << 1 <<  1 << endr
			<< 1 <<  1 << endr
			<< 1 << -1 << endr
			<< 1 << -1 << endr;
		length_f = 2;
		length_t = 2;
		break;
	case CSIRS_CDMType::CDM8:
		w_f << 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr
			<< 1 << 1 << endr
			<< 1 << -1 << endr;
		w_t << 1 << 1 << 1 << 1 << endr
			<< 1 << 1 << 1 << 1 << endr
			<< 1 << -1 << 1 << -1 << endr
			<< 1 << -1 << 1 << -1 << endr
			<< 1 << 1 << -1 << -1 << endr
			<< 1 << 1 << -1 << -1 << endr
			<< 1 << -1 << -1 << 1 << endr
			<< 1 << -1 << -1 << 1 << endr;
		length_f = 2;
		length_t = 4;
		break;
	default:
		perror("CSIRS CDM Type not supported");
		break;
	}

}   /* end of func_CDM() */




/*
 * func_beam_selection() - Beam selection function 
 *		Calculate the received antenna signal of the per-layer transmit beam (the method defined in Type I CSI-RS codebook) for feedback
 *		Type of transmission beam: considering co-phase and Type I CSI-RS codebook
 *		In the case of co-phase feedback for digtal beam, feedback is considered for each RB.
 *
 *	  - INPUT: para, ch_est, Nbeam
 *    - OUTPUT: this->RI, this->Precoder
 */
void func_beam_selection(ModuleParameterMIMO &para, cx_cube &ch_est, cx_mat &Nbeam, int &RI, field<cx_mat> &Precoder) {

	Precoder = field<cx_mat>(para.num_RB, 1);

	/* local variable declaration */
	field<cx_mat>	rec_beam_signal;
	field<cx_mat>	beam_set, beam_set_1, beam_set_2;
	cube			MI_1, MI_2;
	uword			ana_ind_set;
	uword			dig_ind_set;
	cx_mat			phase;
	cx_mat			sq_ch_est_1, sq_ch_est_2, sq_ch_est_3, sq_ch_est_4;
	cx_mat			temp_beam, temp_beam1, temp_beam2;
	cx_mat			V;
	cx_mat			inv_temp;
	cx_mat			U;
	mat				sig_power;
	mat				interf_power(1, 1);
	mat				noise_enhancement;
	mat				SINR_temp(1, 1);
	
	uword			i3;
	umat			ind_set;
	umat			ac_ind_set;

	field<colvec>	Ana_sum;
	uword			max_rank;
	colvec			sq_MI;

	/* Mutual information calculation per beam per co-phase per RB */
	switch (para.Tx_pol) {
	case 1:
		/* Without Cross-Polarization at Tx antenna, only analog precoding is considered */
		for (uword ind_rank = 0; ind_rank < para.RI_max; ind_rank++) {
			switch (ind_rank) {
			case 0:
				rec_beam_signal = field<cx_mat>(Nbeam.n_cols, para.num_RB);
				for (uword i = 0; i < rec_beam_signal.n_elem; i++)
					rec_beam_signal(i) = zeros<cx_mat>(ch_est.n_rows, 1);
				beam_set_1 = field<cx_mat>(Nbeam.n_cols, 1, para.num_RB);
				for (uword i = 0; i < beam_set_2.n_elem; i++)
					beam_set_1(i) = zeros<cx_mat>(Nbeam.n_rows, 1);
				ana_ind_set = rec_beam_signal.n_rows;
				MI_1 = zeros<cube>(ana_ind_set, 1, para.num_RB);

				sq_ch_est_1 = zeros<cx_mat>(ch_est.n_rows, 1);
				V = zeros<cx_mat>(ch_est.n_rows, 1);

				for (uword ind_ana = 0; ind_ana < ana_ind_set; ind_ana++) {
					for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
						for (uword i = 0; i < ch_est.n_rows; i++)		// squeeze ch_est
							sq_ch_est_1(i, 0) = ch_est(i, ind_RB, ind_ana);
						rec_beam_signal(ind_ana, ind_RB) = sq_ch_est_1;
						beam_set_1(ind_ana, 0, ind_RB) = Nbeam.col(ind_ana);
						V = rec_beam_signal(ind_ana, ind_RB);

						switch (para.ue_detection_mode) {
						case ue_detection_mode::MMSE:
							inv_temp = solve((V.t() * V + para.sigma_n_freq * eye<cx_mat>(ind_rank + 1, ind_rank + 1)), V.t());
							U = inv_temp * V;
							break;
						case ue_detection_mode::ZF:
							inv_temp = pinv(V);
							U = inv_temp * V;
							break;
						}

						sig_power = pow(abs(diagvec(U)), 2);
						interf_power = sum(pow(abs(U), 2), 1) - sig_power;
						noise_enhancement = sum(pow(abs(inv_temp), 2), 1);
						SINR_temp = sig_power / (interf_power + para.sigma_n_freq * noise_enhancement);
						MI_1(ind_ana, 0, ind_RB) = sum(sum(log2(1.0 + SINR_temp)));
					}
				}
				break;
			case 1:
				if (para.N1 > para.N2 && para.N2 > 1) {
					i3 = 3;
					ind_set << para.O1     << 0		  << endr
							<< 0		   << para.O2 << endr
							<< 2 * para.O1 << 0		  << endr;
				}
				else if (para.N1 == para.N2) {
					i3 = 3;
					ind_set << para.O1 << 0		  << endr
							<< 0	   << para.O2 << endr
							<< para.O1 << para.O2 << endr;
				}
				else if (para.N1 == 2 && para.N2 == 1) {
					i3 = 1;
					ind_set << para.O1 << 0 << endr;
				}
				else if (para.N1 > 2 && para.N2 == 1) {
					i3 = 3;
					ind_set << para.O1 << 0 << endr
						<< 2 * para.O1 << 0 << endr
						<< 3 * para.O1 << 0 << endr;
				}

				rec_beam_signal = field<cx_mat>(i3 * Nbeam.n_cols, 1, para.num_RB);
				for (uword i = 0; i < rec_beam_signal.n_elem; i++)
					rec_beam_signal(i) = zeros<cx_mat>(ch_est.n_rows, 2);
				beam_set_2 = field<cx_mat>(i3 * Nbeam.n_cols, 1, para.num_RB);
				for (uword i = 0; i < beam_set_2.n_elem; i++)
					beam_set_2(i) = zeros<cx_mat>(Nbeam.n_rows, 2);
				ana_ind_set = rec_beam_signal.n_rows;
				MI_2 = zeros<cube>(ana_ind_set, 1, para.num_RB);

				sq_ch_est_1 = zeros<cx_mat>(ch_est.n_rows, 1);
				sq_ch_est_2 = zeros<cx_mat>(ch_est.n_rows, 1);
				temp_beam = zeros<cx_mat>(Nbeam.n_rows, 2);
				temp_beam1 = zeros<cx_mat>(Nbeam.n_rows, 1);
				temp_beam2 = zeros<cx_mat>(Nbeam.n_rows, 1);
				V = zeros<cx_mat>(ch_est.n_rows, 1);

				for (uword ind_ana = 0; ind_ana < ana_ind_set / i3; ind_ana++) {
					for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
						for (uword ind_temp = 0; ind_temp < i3; ind_temp++) {
							ac_ind_set << 1 + dmod(ind_ana + ind_set(ind_temp, 0) - 1, ana_ind_set / i3) << 1 + dmod(ind_ana + ind_set(ind_temp, 1) - 1, ana_ind_set / i3) << endr;
							for (uword i = 0; i < ch_est.n_rows; i++) {		// squeeze ch_est
								sq_ch_est_1(i, 0) = ch_est(i, ind_RB, ac_ind_set(0));
								sq_ch_est_2(i, 0) = ch_est(i, ind_RB, ac_ind_set(1));
							}
							rec_beam_signal(i3 * ind_ana + ind_temp, 1, ind_RB) = join_horiz(sq_ch_est_1, sq_ch_est_2) / 2;
							temp_beam1 = Nbeam.col(ac_ind_set(0));
							temp_beam2 = Nbeam.col(ac_ind_set(1));
							temp_beam = join_horiz(temp_beam1 / norm(temp_beam1, 2), temp_beam2 / norm(temp_beam2, 2)) / 2;
							beam_set_2(i3 * ind_ana + ind_temp, 0, ind_RB) = temp_beam;

							V = rec_beam_signal(i3 * ind_ana * ind_temp, 0, ind_RB);
							switch (para.ue_detection_mode) {
							case ue_detection_mode::MMSE:
								inv_temp = solve((V.t() * V + para.sigma_n_freq * eye<cx_mat>(ind_rank + 1, ind_rank + 1)), V.t());
								U = inv_temp * V;
								break;
							case ue_detection_mode::ZF:
								inv_temp = pinv(V);
								U = inv_temp * V;
								break;
							}

							sig_power = pow(abs(diagvec(U)), 2);
							interf_power = sum(pow(abs(U), 2), 1) - sig_power;
							noise_enhancement = sum(pow(abs(inv_temp), 2), 1);
							SINR_temp = sig_power / (interf_power + para.sigma_n_freq * noise_enhancement);
							MI_2(i3 * ind_ana + ind_temp, 0, ind_RB) = sum(sum(log2(1.0 + SINR_temp)));
						}
					}
				}
				break;
			}	/* end of inner switch(ind_rank) */

		}	/* end of for (uword ind_rank = 0; ind_rank < para.RI_max; ind_rank++) */
		break;
	case 2:
		for (uword ind_rank = 0; ind_rank < para.RI_max; ind_rank++) {
			switch (ind_rank) {
			case 0:
				rec_beam_signal = field<cx_mat>(Nbeam.n_cols / 2, 4, para.num_RB);
				for (uword i = 0; i < rec_beam_signal.n_elem; i++)
					rec_beam_signal(i) = zeros<cx_mat>(ch_est.n_rows, 1);
				beam_set_1 = field<cx_mat>(Nbeam.n_cols / 2, 4, para.num_RB);
				for (uword i = 0; i < beam_set_1.n_elem; i++)
					beam_set_1(i) = zeros<cx_mat>(Nbeam.n_rows, 1);
				ana_ind_set = rec_beam_signal.n_rows;
				dig_ind_set = rec_beam_signal.n_cols;
				phase << cx_double(1, 0) << exp(cx_double(0, 1) * datum::pi / 2.0) << exp(cx_double(0, 1) * datum::pi) << exp(cx_double(0, 1) * 3.0 * datum::pi / 2.0) << endr;
				MI_1 = zeros<cube>(ana_ind_set, dig_ind_set, para.num_RB);
				sq_ch_est_1 = zeros<cx_mat>(ch_est.n_rows, 1);
				sq_ch_est_2 = zeros<cx_mat>(ch_est.n_rows, 1);
				temp_beam = zeros<cx_mat>(Nbeam.n_rows, 1);
				V = zeros<cx_mat>(ch_est.n_rows, 1);

				for (uword ind_ana = 0; ind_ana < ana_ind_set; ind_ana++) {
					for (uword ind_dig = 0; ind_dig < dig_ind_set; ind_dig++) {
						for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
							for (uword i = 0; i < ch_est.n_rows; i++) {		// squeeze ch_est
								sq_ch_est_1(i, 0) = ch_est(i, ind_RB, ind_ana);
								sq_ch_est_2(i, 0) = ch_est(i, ind_RB, ind_ana + ana_ind_set);
							}
							rec_beam_signal(ind_ana, ind_dig, ind_RB) = (sq_ch_est_1 + phase(ind_dig) * sq_ch_est_2) / sqrt(2);
							temp_beam = Nbeam.col(ind_ana) + phase(ind_dig) * Nbeam.col(ind_ana + ana_ind_set);
							beam_set_1(ind_ana, ind_dig, ind_RB) = temp_beam / norm(temp_beam, 2);
							V = rec_beam_signal(ind_ana, ind_dig, ind_RB);

							switch (para.ue_detection_mode) {
							case ue_detection_mode::MMSE:
								inv_temp = solve((V.t() * V + para.sigma_n_freq * eye<cx_mat>(ind_rank + 1, ind_rank + 1)), V.t());
								U = inv_temp * V;
								break;
							case ue_detection_mode::ZF:
								inv_temp = pinv(V);
								U = inv_temp * V;
								break;
							}

							sig_power = pow(abs(diagvec(U)), 2);
							interf_power = sum(pow(abs(U), 2), 1) - sig_power;
							noise_enhancement = sum(pow(abs(inv_temp), 2), 1);
							SINR_temp = sig_power / (interf_power + para.sigma_n_freq * noise_enhancement);
							MI_1(ind_ana, ind_dig, ind_RB) = sum(sum(log2(1.0 + SINR_temp)));
						}
					}
				}
				break;
			case 1:
				if (para.N1 > para.N2 && para.N2 > 1) {
					i3 = 4;
					ind_set << 0	   << 0		  << endr
						    << para.O1 << 0		  << endr
							<< 0	   << para.O2 << endr
							<< 2 * para.O1 << 0   << endr;
				}
				else if (para.N1 == para.N2) {
					i3 = 4;
					ind_set << 0	   << 0		  << endr
							<< para.O1 << 0		  << endr
							<< 0	   << para.O2 << endr
							<< para.O1 << para.O2 << endr;
				}
				else if (para.N1 == 2 && para.N2 == 1) {
					i3 = 2;
					ind_set << 0	   << 0 << endr
							<< para.O1 << 0 << endr;
				}
				else if (para.N1 > 2 && para.N2 == 1) {
					i3 = 4;
					ind_set << 0		   << 0 << endr
							<< para.O1     << 0 << endr
							<< 2 * para.O1 << 0 << endr
							<< 3 * para.O1 << 0 << endr;
				}

				rec_beam_signal = field<cx_mat>(i3 * Nbeam.n_cols / 2, 2, para.num_RB);
				for (uword i = 0; i < rec_beam_signal.n_elem; i++)
					rec_beam_signal(i) = zeros<cx_mat>(ch_est.n_rows, 2);
				beam_set_2 = field<cx_mat>(i3 * Nbeam.n_cols / 2, 2, para.num_RB);
				for (uword i = 0; i < beam_set_2.n_elem; i++)
					beam_set_2(i) = zeros<cx_mat>(Nbeam.n_rows, 1);
				ana_ind_set = rec_beam_signal.n_rows;
				dig_ind_set = rec_beam_signal.n_cols;
				phase << cx_double(1, 0) << exp(cx_double(0, 1) * datum::pi / 2.0) << endr;
				MI_2 = zeros<cube>(ana_ind_set, dig_ind_set, para.num_RB);

				sq_ch_est_1 = zeros<cx_mat>(ch_est.n_rows, 1);
				sq_ch_est_2 = zeros<cx_mat>(ch_est.n_rows, 1);
				sq_ch_est_3 = zeros<cx_mat>(ch_est.n_rows, 1);
				sq_ch_est_4 = zeros<cx_mat>(ch_est.n_rows, 1);
				temp_beam  = zeros<cx_mat>(Nbeam.n_rows, 2);
				temp_beam1 = zeros<cx_mat>(Nbeam.n_rows, 1);
				temp_beam2 = zeros<cx_mat>(Nbeam.n_rows, 1);
				V = zeros<cx_mat>(ch_est.n_rows, 1);

				for (uword ind_ana = 0; ind_ana < ana_ind_set / i3; ind_ana++) {
					for (uword ind_dig = 0; ind_dig < dig_ind_set; ind_dig++) {
						for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
							for (uword ind_temp = 0; ind_temp < i3; ind_temp++) {
								ac_ind_set << imod(ind_ana + ind_set(ind_temp, 0), ana_ind_set / i3) << imod(ind_ana + ind_set(ind_temp, 1), ana_ind_set / i3) << endr;
								for (uword i = 0; i < ch_est.n_rows; i++) {		// squeeze ch_est
									sq_ch_est_1(i, 0) = ch_est(i, ind_RB, ind_ana);
									sq_ch_est_2(i, 0) = ch_est(i, ind_RB, ind_ana + ana_ind_set / i3);
									sq_ch_est_3(i, 0) = ch_est(i, ind_RB, ac_ind_set(0));
									sq_ch_est_4(i, 0) = ch_est(i, ind_RB, ac_ind_set(1) + ana_ind_set / i3);
								}
								rec_beam_signal(i3 * ind_ana + ind_temp, ind_dig, ind_RB) = join_horiz((sq_ch_est_1 + phase(ind_dig) * sq_ch_est_2), sq_ch_est_3 - phase(ind_dig) * sq_ch_est_4) / 2;
								temp_beam1 = Nbeam.col(ind_ana) + phase(ind_dig) * Nbeam.col(ac_ind_set(1) + ana_ind_set / i3);
								temp_beam2 = Nbeam.col(ac_ind_set(0)) - phase(ind_dig) * Nbeam.col(ac_ind_set(1) + ana_ind_set / i3);
								temp_beam = join_horiz(temp_beam1 / norm(temp_beam1, 2), temp_beam2 / norm(temp_beam2, 2)) / 2;
								beam_set_2(i3 * ind_ana + ind_temp, ind_dig, ind_RB) = temp_beam;

								V = rec_beam_signal(i3 * ind_ana + ind_temp, ind_dig, ind_RB);
								switch (para.ue_detection_mode) {
								case ue_detection_mode::MMSE:
									inv_temp = solve((V.t() * V + para.sigma_n_freq * eye<cx_mat>(ind_rank + 1, ind_rank + 1)), V.t());
									U = inv_temp * V;
									break;
								case ue_detection_mode::ZF:
									inv_temp = pinv(V);
									U = inv_temp * V;
									break;
								}

								sig_power = pow(abs(diagvec(U)), 2);
								interf_power = sum(pow(abs(U), 2), 1) - sig_power;
								noise_enhancement = sum(pow(abs(inv_temp), 2), 1);
								SINR_temp = sig_power / (interf_power + para.sigma_n_freq * noise_enhancement);
								MI_2(i3 * ind_ana + ind_temp, ind_dig, ind_RB) = sum(sum(log2(1.0 + SINR_temp)));
							}
						}
					}
				}
				break;
			}	/* end of inner switch(ind_rank) */

		}	/* end of for (uword ind_rank = 0; ind_rank < para.RI_max; ind_rank++) */
		break;
	}	/* end of switch(para.Tx_pol) */

	/* Analog beam selection */
	Ana_sum = field<colvec>(para.RI_max);
	for (uword ind_rank = 0; ind_rank < para.RI_max; ind_rank++) {
		switch (ind_rank) {
		case 0:
			Ana_sum(ind_rank) = zeros<colvec>(MI_1.n_rows);
			sq_MI = zeros<colvec>(MI_1.n_cols);
			for (uword ind_ana = 0; ind_ana < MI_1.n_rows; ind_ana++)
				for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
					for (uword i = 0; i < MI_1.n_cols; i++)
						sq_MI(i) = MI_1(ind_ana, i, ind_RB);
					Ana_sum(ind_rank)(ind_ana) += max(sq_MI);
				}
			break;
		case 1:
			Ana_sum(ind_rank) = zeros<mat>(MI_2.n_rows, 1);
			sq_MI = zeros<colvec>(MI_2.n_cols);
			for (uword ind_ana = 0; ind_ana < MI_2.n_rows; ind_ana++)
				for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
					for (uword i = 0; i < MI_2.n_cols; i++)
						sq_MI(i) = MI_2(ind_ana, i, ind_RB);
					Ana_sum(ind_rank)(ind_ana, 0) += max(sq_MI);
				}
		}
	}

	if (Ana_sum(0).max() > Ana_sum(1).max())
		max_rank = 1;
	else
		max_rank = 2;

	if (!para.fixed_RI.is_empty())
		max_rank = para.fixed_RI(0);

	uword max_ind_ana;
	max_ind_ana = Ana_sum(max_rank - 1).index_max();

	//cout << "BS.func_beam_management: START digital beam selection" << ", para.num_RB = " << para.num_RB << endl;
	//cout << "max_rank = " << max_rank << endl;
	//cout << "size(Precoder) = " << arma::size(Precoder) << ", size(beam_set_1) = " << arma::size(beam_set_1);
	//cout << ", size(MI_1) = " << arma::size(MI_1) << ", size(MI_2) = " << arma::size(MI_2) << endl;
	
	/* Digital beam selection */
	uvec	max_ind_dig = zeros<uvec>(para.num_RB);
	for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
		switch (max_rank) {
		case 1:
			sq_MI = zeros<colvec>(MI_1.n_cols);
			for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
				for (uword i = 0; i < MI_1.n_cols; i++)
					sq_MI(i) = MI_1(max_ind_ana, i, ind_RB);
				max_ind_dig(ind_RB) = sq_MI.index_max();
				Precoder(ind_RB) = beam_set_1(max_ind_ana, max_ind_dig(ind_RB), ind_RB);
			}
			break;
		case 2:
			sq_MI = zeros<colvec>(MI_2.n_cols);
			for (uword ind_RB = 0; ind_RB < para.num_RB; ind_RB++) {
				for (uword i = 0; i < MI_2.n_cols; i++)
					sq_MI(i) = MI_2(max_ind_ana, i, ind_RB);
				max_ind_dig(ind_RB) = sq_MI.index_max();
				Precoder(ind_RB) = beam_set_2(max_ind_ana, max_ind_dig(ind_RB), ind_RB);
			}
		}
	}
	
	RI = max_rank;

	//cout << "RI = " << RI << endl;

	/* local memory deallocation */
	for (uword i = 0; i < rec_beam_signal.n_elem; i++)
		rec_beam_signal(i).reset();
	rec_beam_signal.reset();
	for (uword i = 0; i < beam_set.n_elem; i++)
		beam_set(i).reset();
	beam_set.reset();
	for (uword i = 0; i < beam_set_1.n_elem; i++)
		beam_set_1(i).reset();
	beam_set_1.reset();
	for (uword i = 0; i < beam_set_2.n_elem; i++)
		beam_set_2(i).reset();
	beam_set_2.reset();
	MI_1.reset();
	MI_2.reset();
	phase.reset();
	sq_ch_est_1.reset();
	sq_ch_est_2.reset();
	sq_ch_est_3.reset();
	sq_ch_est_4.reset();
	temp_beam.reset();
	temp_beam1.reset();
	temp_beam2.reset();
	V.reset();
	inv_temp.reset();
	U.reset();
	sig_power.reset();
	interf_power.reset();
	noise_enhancement.reset();
	SINR_temp.reset();
	ind_set.reset();
	ac_ind_set.reset();
	for (uword i = 0; i < Ana_sum.n_elem; i++)
		Ana_sum(i).reset();
	Ana_sum.reset();
	sq_MI.reset();

}	/* end of func_beam_selection() */



/* end of module_BS_SUMIMO.cpp */