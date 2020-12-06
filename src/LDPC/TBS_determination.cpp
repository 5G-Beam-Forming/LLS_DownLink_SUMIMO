/*--------------------------------------------------------------------------------------*/
// TBS_get.cpp
// - Set the TBS based on TS38.214
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "SISO_PDSCH.h"
#include "misc.h"


int *TBS_determination (int n_PRB, int N_symb_sh, int log2_Mod_order, NR_CODEC *nr_codec )
{
	int *TBS_N_RE;
	TBS_N_RE = new int[2];// [0]: TBS , [1]: N_RE

	int N_DMRS_PRB = 0;
	int N_oh_PRB = 0;
	//int N_RE ;
	int N_RE_p ;
//	int N_RE_p_bar ;		// removal coding
	int N_SC_RB = 12; 	
	
	N_RE_p = N_SC_RB * N_symb_sh - N_DMRS_PRB - N_oh_PRB ;

	//if (N_RE_p <= 9)						N_RE_p_bar = 6 ;		// removal coding
	//else if (N_RE_p>9 && N_RE_p<=15)		N_RE_p_bar = 12 ;		// removal coding
	//else if (N_RE_p>15 && N_RE_p<=30)		N_RE_p_bar = 18 ;		// removal coding
	//else if (N_RE_p>30 && N_RE_p<=57)		N_RE_p_bar = 42 ;		// removal coding
	//else if (N_RE_p>57 && N_RE_p<=90)		N_RE_p_bar = 72 ;		// removal coding
	//else if (N_RE_p>90 && N_RE_p<=126)		N_RE_p_bar = 108 ;	// removal coding
	//else if (N_RE_p>126 && N_RE_p<=150)		N_RE_p_bar = 144 ;	// removal coding
	//else									N_RE_p_bar = 156 ;		// removal coding

	//TBS_N_RE[1] = N_RE_p_bar * n_PRB ;							// removal coding


	if (nr_codec->ldpc.LBRM_indicator == 1){ // for (LBRM)	//	additional coding
		if (n_PRB < 33)			n_PRB = 33;					//	additional coding
		else if (n_PRB < 66)	n_PRB = 66;					//	additional coding
		else if (n_PRB < 107)	n_PRB = 107;				//	additional coding
		else if (n_PRB < 135)	n_PRB = 135;				//	additional coding
		else if (n_PRB < 162)	n_PRB = 162;				//	additional coding
		else if (n_PRB < 217)	n_PRB = 217;				//	additional coding
		else					n_PRB = 273;				//	additional coding

		TBS_N_RE[1] = 156 * n_PRB ;							//	additional coding
	}														//	additional coding
	else													//	additional coding
		TBS_N_RE[1] = min(156,N_RE_p) * n_PRB ;				//	additional coding
	

	int v_layer = 1;		
	
	double N_info =  TBS_N_RE[1] * nr_codec->code_rate * log2_Mod_order * v_layer ;
	int N_info_p ;
	int N_info_round ;
	double n ;
	int diff_temp ;
	double tempa ; 
	//int TBS ;
	int TBS_set[93] = {24,32,40,48,56,64,72,80,88,96,104,112,120,128,136,144,152,160,168,176,184,192,208,224,240,256,272,288,304,320,336,352,368,384,408,432,456,480,504,528,552,576,608,640,672,704,736,768,808,848,888,928,984,1032,1064,1128,1160,1192,1224,1256,1288,1320,1352,1416,1480,1544,1608,1672,1736,1800,1864,1928,2024,2088,2152,2216,2280,2408,2472,2536,2600,2664,2728,2792,2856,2976,3104,3240,3368,3496,3624,3752,3824} ;
	int k ;
	double CC ;
	if (N_info <=3824){
		n = max(3.0, floor( log( N_info )/log(2.0) ) - 6.0) ;
		N_info_p = (int) max(24.0, pow(2.0,n) * floor( N_info / pow(2.0,n)  ) ) ;
		
		for(k=0;k<93;k++){
			diff_temp = TBS_set[k] - N_info_p ;
			if (diff_temp >=0)
				break ;
		}
		TBS_N_RE[0] = TBS_set[k] ;
	}
	else{
		n = floor( log( N_info -24.0 )/log(2.0) ) - 5.0 ;
		tempa = (N_info-24.0) / (pow(2.0,n)) - (int)((N_info-24.0) / (pow(2.0,n))) ;
		
		if (tempa >= 0.5) {
			N_info_round = (int)((N_info-24.0) / (pow(2.0,n))) + 1 ;
		}
		else{
			N_info_round = (int)((N_info-24.0) / (pow(2.0,n))) ;
		}

		//N_info_p = (int) ( pow(2.0,n) * (double)N_info_round ) ;			//  removal coding
		N_info_p = (int)max(3840.0, pow(2.0,n) * (double)N_info_round ) ;	//	additional coding

		if(nr_codec->code_rate <= 0.25){
			CC =  ceil( (double)(N_info_p+24.0) / (3816.0) ) ;
			TBS_N_RE[0] = (int) (8.0 * CC * ceil( (double)(N_info_p+24)/(8.0*CC) ) -24.0) ;
		}
		else{
			if (N_info_p >8424){
				CC =  ceil( (double)(N_info_p+24.0) / (8424.0) ) ;
				TBS_N_RE[0] =  (int) (8.0 * CC * (int)ceil( (N_info_p+24.0)/(8.0 * CC)) - 24.0) ;
			}
			else {
				TBS_N_RE[0] = 8 * (int)ceil( (N_info_p+24.0)/8.0 ) - 24 ;
			}
		}
	}
	return TBS_N_RE  ;

}