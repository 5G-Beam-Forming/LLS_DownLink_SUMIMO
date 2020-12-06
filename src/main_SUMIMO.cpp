/*
 * main_SUMIMO.cpp - A drive program for LLS Down Link SUMIMO simulation.
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by jwpark
 *      Last Updated: 2018-09-17
 *
 * Copyright(c) 2016-2018 KIAST. All rights reserved.
 */

#include <cstdio>
#include <ctime>
#include <armadillo>
#include "sumimo.h"
#include "module_Parameter_MIMO.h"
#include "module_BS_SUMIMO.h"
#include "module_UE_MIMO.h"
#include "module_Channel_5G.h"
#include "LDPC/method_ChannelCoding.h"


using namespace std;
using namespace arma;


int			SLOT;


/* local function declaration */
string prStr_simCase(Sim_Case);
string prStr_chType(Ch_type);
string prStr_chMode(Ch_mode);




/*
 * main() - main drive function for LLS DownLink SUMIMO simulation
 */
int main(int argc, char **argv)
{
	/* local variable declarations */
	int			num_slot_total;
	clock_t		main_start, main_end;
	execTime	exectm;

	exectm.start = clock();

	/* class instantiations */
	ParameterInit			init;
	ModuleParameterMIMO		para;
	ModuleChannel5G			Chan;
	ChannelOutput			Chan_Output;
	ModuleBS_SUMIMO			BS;
	ModuleUE_SUMIMO			UE;
	NR_ChannelCoding		LDPC;

	/* set init values */
	init.sim_Case = Sim_Case::SUMIMO;
	init.mu		  = 0;
	init.BW		  = 10e6;
	init.f		  = 4e9;
	init.num_subframe = 20;
	init.ch_type  = Ch_type::CDL_A;
	init.ch_mode  = Ch_mode::Block_Fading;
	init.Tx_pol   = 2;
	init.pol_tx   = Pol_tx::cross_pol;
	init.num_Tx_antenna			   = 4;
	init.num_Tx_antenna_horizontal = 2;
	init.num_Tx_antenna_vertical   = 2;
	init.Rx_pol	= 2;
	init.pol_rx = Pol_rx::cross_pol;
	init.num_Rx_antenna			   = 1;
	init.num_Rx_antenna_horizontal = 1;
	init.num_Rx_antenna_vertical   = 1;
	init.user_speed = 3.56;
	init.harq_switch = HARQ_switch::off;

	irowvec		CQI(CQI_LENGTH);

	/* for CQI test */
	CQI(0) = 1;
	unsigned int SNR_stepsize_num = 4;

	/* program banner printing */
	cout << endl << endl;
	cout << " *************************[ Simulation Starts ]*************************\n";
	cout << "               Simulation case           : " << prStr_simCase(init.sim_Case) << endl;
	cout << "               Bandwidth                 : "; printf("%5.2f [MHz]\n", init.BW / 1e6);
	cout << "               Channel type              : " << prStr_chType(init.ch_type) << endl;
	cout << "               Channel mode              : " << prStr_chMode(init.ch_mode) << endl;
	cout << "               Subframe                  : "; printf("%d\n", init.num_subframe);
	cout << "               Carrier frequency         : "; printf("%4.2f[GHz]\n", init.f / 1e9);
	cout << "               num_Tx_antenna_horizontal : " << init.num_Tx_antenna_horizontal << endl;
	cout << "               num_Tx_antenna_vertical   : " << init.num_Tx_antenna_vertical << endl;
	cout << "               num_Rx_antenna            : " << init.num_Rx_antenna << endl;
	cout << " ***********************************************************************" << endl << endl << endl;

	field<rowvec>	BLER;
	field<rowvec>	Throughput;
	rowvec			SNR_start(SNR_LENGTH);
	rowvec			SNR_end(SNR_LENGTH);
	rowvec			SNR_stepsize(SNR_LENGTH);
	mat				SNR_range;

	SNR_start << 16.0 << -14.0 << -13.0 << -10.0 << -7.0 << -5.0 << -3.0 << -1.0 <<  1.0 <<  4.0 <<  6.0 <<  8.0 << 10.0 << 12.0 << endr;
	SNR_end   << 26.0 <<  10.0 <<  11.0 <<  12.0 << 15.0 << 18.0 << 21.0 << 24.0 << 27.0 << 30.0 << 33.0 << 36.0 << 39.0 << 42.0 << endr;
	SNR_end += 30.0;

	SNR_stepsize = (SNR_end - SNR_start) / ((double)SNR_stepsize_num - 1.0);
	SNR_range = zeros<mat>(SNR_stepsize_num, SNR_LENGTH);
	SNR_range.row(0) = SNR_start;
	for (uword ind_cqi = 0; ind_cqi < SNR_start.n_elem; ind_cqi++) {
		if (SNR_stepsize(ind_cqi) != 0)
			for (uword row = 1; row < SNR_stepsize_num; row++)
				 SNR_range(row, ind_cqi) = SNR_range(row - 1, ind_cqi) + SNR_stepsize(ind_cqi);
	}

	para.module_Parameter_MIMO(init);
	para.Tx_pol = init.Tx_pol;
	para.Rx_pol = init.Rx_pol;
	para.fixed_RI = zeros<irowvec>(1);
	para.fixed_RI(0) = 1;

	num_slot_total = para.num_slot_in_subframe * init.num_subframe;

	BLER = field<rowvec>(CQI.n_elem);
	Throughput = field<rowvec>(CQI.n_elem);
	for (uword i = 0; i < CQI.n_elem; i++) {
		BLER(i) = zeros<rowvec>(SNR_stepsize_num);
		Throughput(i) = zeros<rowvec>(SNR_stepsize_num);
	}

	SLOT = 0;
	/* main loop */
	for (uword ind_cqi = 0; ind_cqi < CQI.n_elem; ind_cqi++) {
		printf("==========================[ cqi = %llu, (%llu of %llu) ]==========================\n", CQI(ind_cqi), ind_cqi + 1, CQI.size());
		init.SNR_range = SNR_range.col(CQI(ind_cqi) - 1);

		for (uword ind_SNR = 0; ind_SNR < SNR_stepsize_num; ind_SNR++) {
			para.Reset(init.SNR_range(ind_SNR), ind_SNR);
			printf("  ----------------[ cqi = %llu, SNR = %5.2f dB, (%2d of %lu) ]----------------\n", CQI(ind_cqi), para.SNR, para.ind_SNR + 1, init.SNR_range.n_elem);
			BS.module_BS_SUMIMO(para);
			UE.module_UE_SUMIMO(para);
			UE.feedBack.CQI = zeros<irowvec>(CQI_LENGTH);
			for (int ind_slot = 0; ind_slot < num_slot_total; ind_slot++) {
				SLOT = ind_slot + 1;
				printf("  ..............[ cqi = %llu, SNR = %2d (%5.2f dB), slot = %d ]..............\n", CQI(ind_cqi), para.ind_SNR + 1, para.SNR, ind_slot + 1);
				Chan.module_Channel_5G(para);
				Chan_Output = Chan.Output;
				UE.feedBack.CQI(0) = CQI(ind_cqi);
				BS.Process(para, LDPC, ind_slot, UE.feedBack, Chan_Output);
				UE.process(para, LDPC, BS.genie, Chan_Output, ind_slot);
			}

			/* Result */
			BLER(ind_cqi)(ind_SNR) = 1.0 - (double)sum(UE.result.ACK) / (double)UE.HARQ_buffer.num_TB;
			Throughput(ind_cqi)(ind_SNR) = (double)sum(UE.result.suc_data) / ((double)para.num_subframe * 0.001) / pow(10.0, 6.0);
		}

	}	/* end of for main loop */

	/* Test print the Result */
	cout << "\n\n*********************************[ PROGRAM RESULT ]*********************************" << endl << endl;
	BLER.print("BLER:");
	cout << endl;
	Throughput.print("Throughtput:");
	cout << "***********************************[ PROGRAM END ]***********************************" << endl << endl;

	/* write result to the CSV file */
	mat BLER_MAT = zeros<mat>(BLER.n_elem, BLER(0).n_elem);
	mat Throughput_MAT = zeros<mat>(Throughput.n_elem, Throughput(0).n_elem);
	for (uword i = 0; i < BLER.n_elem; i++)
		for (uword j = 0; j < BLER(0).n_elem; j++)
			BLER_MAT(i, j) = BLER(i)(j);

	for (uword i = 0; i < Throughput.n_elem; i++)
		for (uword j = 0; j < Throughput(0).n_elem; j++)
			Throughput_MAT(i, j) = Throughput(i)(j);

	BLER_MAT.save("OUTPUT_BLER.csv", csv_ascii);
	Throughput_MAT.save("OUTPUT_Throughput.csv", csv_ascii);

	/* Program Ending Display */
	printf("\n%s terminates normally. Program success! Execution Time is ", argv[0]);
	exectm.end = clock();
	exectm.exeTime = exectm.end - exectm.start;
	exectm.ms = exectm.exeTime % 60;
	exectm.sec = exectm.exeTime / CLOCKS_PER_SEC;
	if (exectm.sec >= 60 && exectm.sec < 3600) {
		exectm.min = exectm.sec / 60;
		exectm.sec = exectm.sec % 60;
		printf(" %d ms (%d min %d sec %d ms).\n\n", exectm.exeTime, exectm.min, exectm.sec, exectm.ms);
	}
	else if (exectm.sec >= 3600) {
		exectm.min = exectm.sec / 60;
		exectm.hr = exectm.min / 60;
		exectm.min = exectm.min % 60;
		exectm.sec = exectm.sec % 60;
		printf(" %d ms (%d hour %d min %d sec %d ms).\n\n", exectm.exeTime, exectm.hr, exectm.min, exectm.sec, exectm.ms);
	}
	else {
		printf(" %d ms (%d sec %d ms).\n\n", exectm.exeTime, exectm.sec, exectm.ms);
	}

}	/* end of main() */




/*
 * Misc local function definitions
 *****************************************************************************************/

string prStr_simCase(Sim_Case sim_case) {

	string strBox;

	switch (sim_case) {
	case Sim_Case::SUSISO:				strBox = "SUSISO";				break;
	case Sim_Case::SUMIMO:				strBox = "SUMIMO";				break;
	case Sim_Case::SUMIMO_CLSM:			strBox = "SUMIMO_CLSM";			break;
	case Sim_Case::SUMIMO_Noncodebook:	strBox = "SUMIMO_Noncodebook";	break;
	case Sim_Case::MUMIMO:				strBox = "MUMIMO";				break;
	case Sim_Case::MUMIMO_CLSM:			strBox = "MUMIMO_CLSM";			break;
	case Sim_Case::MUMIMO_Noncodebook:	strBox = "MUMIMO_Noncodebook";	break;
	default:							strBox = "UNSUPPORTED (Unknown)";
	}
	return strBox;

}   /* end of simCase() method */



string prStr_chType(Ch_type CType) {

	string strBox;

	switch (CType) {
	case Ch_type::PedA:				strBox = "PedA";			break;
	case Ch_type::flat_Rayleigh:	strBox = "flat_Rayleigh";	break;
	case Ch_type::Rayleigh2:		strBox = "Rayleigh2";		break;
	case Ch_type::PedB:				strBox = "PedB";			break;
	case Ch_type::PedBcorr:			strBox = "PedBcorr";		break;
	case Ch_type::VehA:				strBox = "VehA";			break;
	case Ch_type::VehB:				strBox = "VehB";			break;
	case Ch_type::EPA5Hz:			strBox = "EPA5Hz";			break;
	case Ch_type::EVA5Hz:			strBox = "EVA5Hz";			break;
	case Ch_type::EVA70Hz:			strBox = "EVA70Hz";			break;
	case Ch_type::ETU70Hz:			strBox = "ETU70Hz";			break;
	case Ch_type::ETU300Hz:			strBox = "ETU300Hz";		break;
	case Ch_type::AWGN:				strBox = "AWGN";			break;
	case Ch_type::TDL_A:			strBox = "TDL_A";			break;
	case Ch_type::TDL_B:			strBox = "TDL_B";			break;
	case Ch_type::TDL_C:			strBox = "TDL_C";			break;
	case Ch_type::TDL_D:			strBox = "TDL_D";			break;
	case Ch_type::TDL_E:			strBox = "TDL_E";			break;
	case Ch_type::CDL_A:			strBox = "CDL_A";			break;
	case Ch_type::CDL_B:			strBox = "CDL_B";			break;
	case Ch_type::CDL_C:			strBox = "CDL_C";			break;
	case Ch_type::CDL_D:			strBox = "CDL_D";			break;
	case Ch_type::CDL_E:			strBox = "CDL_E";			break;
	default:						strBox = "UNSUPPORTED";
	}
	return strBox;

}   /* end of chType() method */



string prStr_chMode(Ch_mode cMod) {

	string strBox;

	switch (cMod) {
	case Ch_mode::Block_Fading:		strBox = "Block_Fading";	break;
	case Ch_mode::Fast_Fading:		strBox = "Fast_Fading";		break;
	case Ch_mode::flat_Rayleigh:	strBox = "flat_Rayleigh";	break;
	default:						strBox = "UNSUPPORTED";
	}
	return strBox;

}   /* end of chMode() method */



/* end of main_SUMIMO.cpp */