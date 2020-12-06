/*--------------------------------------------------------------------------------------*/
// method_QAM_Modulation.cpp
// - Class part for QAM modulation
// - [REF] 3GPP TS 38.211
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/


#include "stdafx.h"
#include "method_QAM_modulation.h"


/*-------------- Construction & Destruction --------------*/

QAM_Modulation::QAM_Modulation(int M_ary, int TotalCodeLength, QAM *qam)
{
	/*----- General setting ------*/
	ModOrder = M_ary ;
	ModMapper = Func_QAM_set_gen (ModOrder) ;
	log2M = (int) ( log((double) ModOrder) / log(2.0) ) ;	
	//SymbolLength = (int) floor((double)TotalCodeLength / (double)log2M) ;
	SymbolLength = (int) ceil((double)TotalCodeLength / (double)log2M) ;
	//SymbolLength = TotalCodeLength / log2M ;


	/*----- Modulator ------*/
	qam->tx_real = new double[SymbolLength] ;
	qam->tx_imag = new double[SymbolLength] ;


	/*----- Demodulator (LLR calculator) ------*/
	qam->rx_real = new double[SymbolLength] ;
	qam->rx_imag = new double[SymbolLength] ;

	pdf = new double*[SymbolLength] ;
	for (int i=0; i<SymbolLength; i++ )
		pdf[i] = new double[ModOrder] ;

	LLR = new double[TotalCodeLength] ;

}


QAM_Modulation::~QAM_Modulation()
{

}





/*-------------- Public methods --------------*/

int QAM_Modulation::GetSymbolLength ()
{
	return SymbolLength ;
}


// Tx side
QAM QAM_Modulation::Tx_QAM_Modulator (int *codeword, QAM *qam)
{
	int idx ;
	int *log2M_bits = new int[log2M] ;


	for (int i=0; i<SymbolLength; i++) {

		for (int j=0; j<log2M; j++){
			log2M_bits[j] = codeword[log2M*i + j] ;
		}
		
		idx = Func_Bin2Dec (log2M_bits, log2M) ;
		
		qam->tx_real[i] = ModMapper[idx][0] ;	// Real
		qam->tx_imag[i] = ModMapper[idx][1] ;	// Imag
	}

	delete []log2M_bits ;

	return *qam ;
}


// Rx side
double* QAM_Modulation::Rx_QAM_Demodulator (double noise_var, QAM *qam)
{
	// Calculate Gaussian pdf's
	for (int i=0; i<SymbolLength; i++ ){
		for (int j=0; j<ModOrder; j++){
			pdf[i][j] = exp ( - ( pow(qam->rx_real[i] - ModMapper[j][0], 2.0)
						        + pow(qam->rx_imag[i] - ModMapper[j][1], 2.0) ) / noise_var ) ;
		}
	}
		

	// LLR calculation 

	if (ModOrder == 2)
		Func_rx_Demod_M2 () ;		// BPSK
	else if (ModOrder == 4)
		Func_rx_Demod_M4 () ;		// QPSK
	else if (ModOrder == 16)
		Func_rx_Demod_M16 () ;		// 16QAM
	else if (ModOrder == 64)
		Func_rx_Demod_M64 () ;		// 64QAM
	else if (ModOrder == 256)
		Func_rx_Demod_M256 () ;		// 256QAM

	// LLR limitting
	Func_rx_LLR_Limitter () ;

	return LLR ;
}





/*-------------- Private methods --------------*/

// Tx side
int QAM_Modulation::Func_Bin2Dec (int *BitSequence, int length)
{
	int Decimal_Num = 0 ;

	for (int i=0; i<length; i++)
		Decimal_Num += (int) pow(2.0, length-i-1) * BitSequence[i] ;

	return  Decimal_Num ;
}


double** QAM_Modulation::Func_QAM_set_gen (int M_ary) 
{
	double **mod_map = new double*[M_ary] ;
	for (int i=0; i<M_ary; i++)
		mod_map[i] = new double[2] ;

	if (M_ary == 2)	// BPSK
	{
		mod_map[0][0] = 1 ;			mod_map[0][1] = 1 ;
		mod_map[1][0] = -1 ;		mod_map[1][1] = -1 ;


		for (int i=0; i<M_ary; i++)
			for (int j=0; j<2; j++)
				mod_map[i][j] /= sqrt(2.0) ;
	}
	else if (M_ary == 4)		// QPSK
	{
		mod_map[0][0] = 1 ;			mod_map[0][1] = 1 ;
		mod_map[1][0] = 1 ;			mod_map[1][1] = -1 ;
		mod_map[2][0] = -1 ;		mod_map[2][1] = 1 ;
		mod_map[3][0] = -1 ;		mod_map[3][1] = -1 ;


		for (int i=0; i<M_ary; i++)
			for (int j=0; j<2; j++)
				mod_map[i][j] /= sqrt(2.0) ;
	}
	else if (M_ary == 16)		// 16QAM
	{
		mod_map[0][0] = 1 ;			mod_map[0][1] = 1 ;
		mod_map[1][0] = 1 ;			mod_map[1][1] = 3 ;
		mod_map[2][0] = 3 ;			mod_map[2][1] = 1 ;
		mod_map[3][0] = 3 ;			mod_map[3][1] = 3 ;
		mod_map[4][0] = 1 ;			mod_map[4][1] = -1 ;
		mod_map[5][0] = 1 ;			mod_map[5][1] = -3 ;
		mod_map[6][0] = 3 ;			mod_map[6][1] = -1 ;
		mod_map[7][0] = 3 ;			mod_map[7][1] = -3 ;
		mod_map[8][0] = -1 ;		mod_map[8][1] = 1 ;
		mod_map[9][0] = -1 ;		mod_map[9][1] = 3 ;
		mod_map[10][0] = -3 ;		mod_map[10][1] = 1 ;
		mod_map[11][0] = -3 ;		mod_map[11][1] = 3 ;
		mod_map[12][0] = -1 ;		mod_map[12][1] = -1 ;
		mod_map[13][0] = -1 ;		mod_map[13][1] = -3 ;
		mod_map[14][0] = -3 ;		mod_map[14][1] = -1 ;
		mod_map[15][0] = -3 ;		mod_map[15][1] = -3 ;


		for (int i=0; i<M_ary; i++)
			for (int j=0; j<2; j++)
				mod_map[i][j] /= sqrt(10.0) ;
	}
	else if (M_ary == 64)		// 64QAM
	{
		mod_map[0][0] = 3 ;			mod_map[0][1] = 3 ;
		mod_map[1][0] = 3 ;			mod_map[1][1] = 1 ;
		mod_map[2][0] = 1 ;			mod_map[2][1] = 3 ;
		mod_map[3][0] = 1 ;			mod_map[3][1] = 1 ;
		mod_map[4][0] = 3 ;			mod_map[4][1] = 5 ;
		mod_map[5][0] = 3 ;			mod_map[5][1] = 7 ;
		mod_map[6][0] = 1 ;			mod_map[6][1] = 5 ;
		mod_map[7][0] = 1 ;			mod_map[7][1] = 7 ;
		mod_map[8][0] = 5 ;			mod_map[8][1] = 3 ;
		mod_map[9][0] = 5 ;			mod_map[9][1] = 1 ;
		mod_map[10][0] = 7 ;		mod_map[10][1] = 3 ;
		mod_map[11][0] = 7 ;		mod_map[11][1] = 1 ;
		mod_map[12][0] = 5 ;		mod_map[12][1] = 5 ;
		mod_map[13][0] = 5 ;		mod_map[13][1] = 7 ;
		mod_map[14][0] = 7 ;		mod_map[14][1] = 5 ;
		mod_map[15][0] = 7 ;		mod_map[15][1] = 7 ;
		mod_map[16][0] = 3 ;		mod_map[16][1] = -3 ;
		mod_map[17][0] = 3 ;		mod_map[17][1] = -1 ;
		mod_map[18][0] = 1 ;		mod_map[18][1] = -3 ;
		mod_map[19][0] = 1 ;		mod_map[19][1] = -1 ;
		mod_map[20][0] = 3 ;		mod_map[20][1] = -5 ;
		mod_map[21][0] = 3 ;		mod_map[21][1] = -7 ;
		mod_map[22][0] = 1 ;		mod_map[22][1] = -5 ;
		mod_map[23][0] = 1 ;		mod_map[23][1] = -7 ;
		mod_map[24][0] = 5 ;		mod_map[24][1] = -3 ;
		mod_map[25][0] = 5 ;		mod_map[25][1] = -1 ;
		mod_map[26][0] = 7 ;		mod_map[26][1] = -3 ;
		mod_map[27][0] = 7 ;		mod_map[27][1] = -1 ;
		mod_map[28][0] = 5 ;		mod_map[28][1] = -5 ;
		mod_map[29][0] = 5 ;		mod_map[29][1] = -7 ;
		mod_map[30][0] = 7 ;		mod_map[30][1] = -5 ;
		mod_map[31][0] = 7 ;		mod_map[31][1] = -7 ;
		mod_map[32][0] = -3 ;		mod_map[32][1] = 3 ;
		mod_map[33][0] = -3 ;		mod_map[33][1] = 1 ;
		mod_map[34][0] = -1 ;		mod_map[34][1] = 3 ;
		mod_map[35][0] = -1 ;		mod_map[35][1] = 1 ;
		mod_map[36][0] = -3 ;		mod_map[36][1] = 5 ;
		mod_map[37][0] = -3 ;		mod_map[37][1] = 7 ;
		mod_map[38][0] = -1 ;		mod_map[38][1] = 5 ;
		mod_map[39][0] = -1 ;		mod_map[39][1] = 7 ;
		mod_map[40][0] = -5 ;		mod_map[40][1] = 3 ;
		mod_map[41][0] = -5 ;		mod_map[41][1] = 1 ;
		mod_map[42][0] = -7 ;		mod_map[42][1] = 3 ;
		mod_map[43][0] = -7 ;		mod_map[43][1] = 1 ;
		mod_map[44][0] = -5 ;		mod_map[44][1] = 5 ;
		mod_map[45][0] = -5 ;		mod_map[45][1] = 7 ;
		mod_map[46][0] = -7 ;		mod_map[46][1] = 5 ;
		mod_map[47][0] = -7 ;		mod_map[47][1] = 7 ;
		mod_map[48][0] = -3 ;		mod_map[48][1] = -3 ;
		mod_map[49][0] = -3 ;		mod_map[49][1] = -1 ;
		mod_map[50][0] = -1 ;		mod_map[50][1] = -3 ;
		mod_map[51][0] = -1 ;		mod_map[51][1] = -1 ;
		mod_map[52][0] = -3 ;		mod_map[52][1] = -5 ;
		mod_map[53][0] = -3 ;		mod_map[53][1] = -7 ;
		mod_map[54][0] = -1 ;		mod_map[54][1] = -5 ;
		mod_map[55][0] = -1 ;		mod_map[55][1] = -7 ;
		mod_map[56][0] = -5 ;		mod_map[56][1] = -3 ;
		mod_map[57][0] = -5 ;		mod_map[57][1] = -1 ;
		mod_map[58][0] = -7;		mod_map[58][1] = -3 ;
		mod_map[59][0] = -7 ;		mod_map[59][1] = -1 ;
		mod_map[60][0] = -5 ;		mod_map[60][1] = -5 ;
		mod_map[61][0] = -5 ;		mod_map[61][1] = -7 ;
		mod_map[62][0] = -7 ;		mod_map[62][1] = -5 ;
		mod_map[63][0] = -7 ;		mod_map[63][1] = -7 ;


		for (int i=0; i<M_ary; i++)
			for (int j=0; j<2; j++)
				mod_map[i][j] /= sqrt(42.0) ;
	}
	else if (M_ary == 256)		// 256QAM
	{
		mod_map[0][0] = 5 ;			mod_map[0][1] = 5 ;
		mod_map[1][0] = 5 ;			mod_map[1][1] = 7 ;
		mod_map[2][0] = 7 ;			mod_map[2][1] = 5 ;
		mod_map[3][0] = 7 ;			mod_map[3][1] = 7 ;
		mod_map[4][0] = 5 ;			mod_map[4][1] = 3 ;
		mod_map[5][0] = 5 ;			mod_map[5][1] = 1 ;
		mod_map[6][0] = 7 ;			mod_map[6][1] = 3 ;
		mod_map[7][0] = 7 ;			mod_map[7][1] = 1 ;
		mod_map[8][0] = 3 ;			mod_map[8][1] = 5 ;
		mod_map[9][0] = 3 ;			mod_map[9][1] = 7 ;
		mod_map[10][0] = 1 ;		mod_map[10][1] = 5 ;
		mod_map[11][0] = 1 ;		mod_map[11][1] = 7 ;
		mod_map[12][0] = 3 ;		mod_map[12][1] = 3 ;
		mod_map[13][0] = 3 ;		mod_map[13][1] = 1 ;
		mod_map[14][0] = 1 ;		mod_map[14][1] = 3 ;
		mod_map[15][0] = 1 ;		mod_map[15][1] = 1 ;
		mod_map[16][0] = 5 ;		mod_map[16][1] = 11 ;
		mod_map[17][0] = 5 ;		mod_map[17][1] = 9 ;
		mod_map[18][0] = 7 ;		mod_map[18][1] = 11 ;
		mod_map[19][0] = 7 ;		mod_map[19][1] = 9 ;
		mod_map[20][0] = 5 ;		mod_map[20][1] = 13 ;
		mod_map[21][0] = 5 ;		mod_map[21][1] = 15 ;
		mod_map[22][0] = 7 ;		mod_map[22][1] = 13 ;
		mod_map[23][0] = 7 ;		mod_map[23][1] = 15 ;
		mod_map[24][0] = 3 ;		mod_map[24][1] = 11 ;
		mod_map[25][0] = 3 ;		mod_map[25][1] = 9 ;
		mod_map[26][0] = 1 ;		mod_map[26][1] = 11 ;
		mod_map[27][0] = 1 ;		mod_map[27][1] = 9 ;
		mod_map[28][0] = 3 ;		mod_map[28][1] = 13 ;
		mod_map[29][0] = 3 ;		mod_map[29][1] = 15 ;
		mod_map[30][0] = 1 ;		mod_map[30][1] = 13 ;
		mod_map[31][0] = 1 ;		mod_map[31][1] = 15 ;
		mod_map[32][0] = 11 ;		mod_map[32][1] = 5 ;
		mod_map[33][0] = 11 ;		mod_map[33][1] = 7 ;
		mod_map[34][0] = 9 ;		mod_map[34][1] = 5 ;
		mod_map[35][0] = 9 ;		mod_map[35][1] = 7 ;
		mod_map[36][0] = 11 ;		mod_map[36][1] = 3 ;
		mod_map[37][0] = 11 ;		mod_map[37][1] = 1 ;
		mod_map[38][0] = 9 ;		mod_map[38][1] = 3 ;
		mod_map[39][0] = 9 ;		mod_map[39][1] = 1 ;
		mod_map[40][0] = 13 ;		mod_map[40][1] = 5 ;
		mod_map[41][0] = 13 ;		mod_map[41][1] = 7 ;
		mod_map[42][0] = 15 ;		mod_map[42][1] = 5 ;
		mod_map[43][0] = 15 ;		mod_map[43][1] = 7 ;
		mod_map[44][0] = 13 ;		mod_map[44][1] = 3 ;
		mod_map[45][0] = 13 ;		mod_map[45][1] = 1 ;
		mod_map[46][0] = 15 ;		mod_map[46][1] = 3 ;
		mod_map[47][0] = 15 ;		mod_map[47][1] = 1 ;
		mod_map[48][0] = 11 ;		mod_map[48][1] = 11 ;
		mod_map[49][0] = 11 ;		mod_map[49][1] = 9 ;
		mod_map[50][0] = 9 ;		mod_map[50][1] = 11 ;
		mod_map[51][0] = 9 ;		mod_map[51][1] = 9 ;
		mod_map[52][0] = 11 ;		mod_map[52][1] = 13 ;
		mod_map[53][0] = 11 ;		mod_map[53][1] = 15 ;
		mod_map[54][0] = 9 ;		mod_map[54][1] = 13 ;
		mod_map[55][0] = 9 ;		mod_map[55][1] = 15 ;
		mod_map[56][0] = 13 ;		mod_map[56][1] = 11 ;
		mod_map[57][0] = 13 ;		mod_map[57][1] = 9 ;
		mod_map[58][0] = 15 ;		mod_map[58][1] = 11 ;
		mod_map[59][0] = 15 ;		mod_map[59][1] = 9 ;
		mod_map[60][0] = 13 ;		mod_map[60][1] = 13 ;
		mod_map[61][0] = 13 ;		mod_map[61][1] = 15 ;
		mod_map[62][0] = 15 ;		mod_map[62][1] = 13 ;
		mod_map[63][0] = 15 ;		mod_map[63][1] = 15 ;
		mod_map[64][0] = 5 ;		mod_map[64][1] = -5 ;
		mod_map[65][0] = 5 ;		mod_map[65][1] = -7 ;
		mod_map[66][0] = 7 ;		mod_map[66][1] = -5 ;
		mod_map[67][0] = 7 ;		mod_map[67][1] = -7 ;
		mod_map[68][0] = 5 ;		mod_map[68][1] = -3 ;
		mod_map[69][0] = 5 ;		mod_map[69][1] = -1 ;
		mod_map[70][0] = 7 ;		mod_map[70][1] = -3 ;
		mod_map[71][0] = 7 ;		mod_map[71][1] = -1 ;
		mod_map[72][0] = 3 ;		mod_map[72][1] = -5 ;
		mod_map[73][0] = 3 ;		mod_map[73][1] = -7 ;
		mod_map[74][0] = 1 ;		mod_map[74][1] = -5 ;
		mod_map[75][0] = 1 ;		mod_map[75][1] = -7 ;
		mod_map[76][0] = 3 ;		mod_map[76][1] = -3 ;
		mod_map[77][0] = 3 ;		mod_map[77][1] = -1 ;
		mod_map[78][0] = 1 ;		mod_map[78][1] = -3 ;
		mod_map[79][0] = 1 ;		mod_map[79][1] = -1 ;
		mod_map[80][0] = 5 ;		mod_map[80][1] = -11 ;
		mod_map[81][0] = 5 ;		mod_map[81][1] = -9 ;
		mod_map[82][0] = 7 ;		mod_map[82][1] = -11 ;
		mod_map[83][0] = 7 ;		mod_map[83][1] = -9 ;
		mod_map[84][0] = 5 ;		mod_map[84][1] = -13 ;
		mod_map[85][0] = 5 ;		mod_map[85][1] = -15 ;
		mod_map[86][0] = 7 ;		mod_map[86][1] = -13 ;
		mod_map[87][0] = 7 ;		mod_map[87][1] = -15 ;
		mod_map[88][0] = 3 ;		mod_map[88][1] = -11 ;
		mod_map[89][0] = 3 ;		mod_map[89][1] = -9 ;
		mod_map[90][0] = 1 ;		mod_map[90][1] = -11 ;
		mod_map[91][0] = 1 ;		mod_map[91][1] = -9 ;
		mod_map[92][0] = 3 ;		mod_map[92][1] = -13 ;
		mod_map[93][0] = 3 ;		mod_map[93][1] = -15 ;
		mod_map[94][0] = 1 ;		mod_map[94][1] = -13 ;
		mod_map[95][0] = 1 ;		mod_map[95][1] = -15 ;
		mod_map[96][0] = 11 ;		mod_map[96][1] = -5 ;
		mod_map[97][0] = 11 ;		mod_map[97][1] = -7 ;
		mod_map[98][0] = 9 ;		mod_map[98][1] = -5 ;
		mod_map[99][0] = 9 ;		mod_map[99][1] = -7 ;
		mod_map[100][0] = 11 ;		mod_map[100][1] = -3 ;
		mod_map[101][0] = 11 ;		mod_map[101][1] = -1 ;
		mod_map[102][0] = 9 ;		mod_map[102][1] = -3 ;
		mod_map[103][0] = 9 ;		mod_map[103][1] = -1 ;
		mod_map[104][0] = 13 ;		mod_map[104][1] = -5 ;
		mod_map[105][0] = 13 ;		mod_map[105][1] = -7 ;
		mod_map[106][0] = 15 ;		mod_map[106][1] = -5 ;
		mod_map[107][0] = 15 ;		mod_map[107][1] = -7 ;
		mod_map[108][0] = 13 ;		mod_map[108][1] = -3 ;
		mod_map[109][0] = 13 ;		mod_map[109][1] = -1 ;
		mod_map[110][0] = 15;		mod_map[110][1] = -3 ;
		mod_map[111][0] = 15 ;		mod_map[111][1] = -1 ;
		mod_map[112][0] = 11 ;		mod_map[112][1] = -11 ;
		mod_map[113][0] = 11 ;		mod_map[113][1] = -9 ;
		mod_map[114][0] = 9 ;		mod_map[114][1] = -11 ;
		mod_map[115][0] = 9 ;		mod_map[115][1] = -9 ;
		mod_map[116][0] = 11 ;		mod_map[116][1] = -13 ;
		mod_map[117][0] = 11 ;		mod_map[117][1] = -15 ;
		mod_map[118][0] = 9 ;		mod_map[118][1] = -13 ;
		mod_map[119][0] = 9 ;		mod_map[119][1] = -15 ;
		mod_map[120][0] = 13 ;		mod_map[120][1] = -11 ;
		mod_map[121][0] = 13 ;		mod_map[121][1] = -9 ;
		mod_map[122][0] = 15 ;		mod_map[122][1] = -11 ;
		mod_map[123][0] = 15 ;		mod_map[123][1] = -9 ;
		mod_map[124][0] = 13 ;		mod_map[124][1] = -13 ;
		mod_map[125][0] = 13 ;		mod_map[125][1] = -15 ;
		mod_map[126][0] = 15 ;		mod_map[126][1] = -13 ;
		mod_map[127][0] = 15 ;		mod_map[127][1] = -15 ;
		mod_map[128][0] = -5 ;		mod_map[128][1] = 5 ;
		mod_map[129][0] = -5 ;		mod_map[129][1] = 7 ;
		mod_map[130][0] = -7 ;		mod_map[130][1] = 5 ;
		mod_map[131][0] = -7 ;		mod_map[131][1] = 7 ;
		mod_map[132][0] = -5 ;		mod_map[132][1] = 3 ;
		mod_map[133][0] = -5 ;		mod_map[133][1] = 1 ;
		mod_map[134][0] = -7 ;		mod_map[134][1] = 3 ;
		mod_map[135][0] = -7 ;		mod_map[135][1] = 1 ;
		mod_map[136][0] = -3 ;		mod_map[136][1] = 5 ;
		mod_map[137][0] = -3 ;		mod_map[137][1] = 7 ;
		mod_map[138][0] = -1 ;		mod_map[138][1] = 5 ;
		mod_map[139][0] = -1 ;		mod_map[139][1] = 7 ;
		mod_map[140][0] = -3 ;		mod_map[140][1] = 3 ;
		mod_map[141][0] = -3 ;		mod_map[141][1] = 1 ;
		mod_map[142][0] = -1 ;		mod_map[142][1] = 3 ;
		mod_map[143][0] = -1 ;		mod_map[143][1] = 1 ;
		mod_map[144][0] = -5 ;		mod_map[144][1] = 11 ;
		mod_map[145][0] = -5 ;		mod_map[145][1] = 9 ;
		mod_map[146][0] = -7 ;		mod_map[146][1] = 11 ;
		mod_map[147][0] = -7 ;		mod_map[147][1] = 9 ;
		mod_map[148][0] = -5 ;		mod_map[148][1] = 13 ;
		mod_map[149][0] = -5 ;		mod_map[149][1] = 15 ;
		mod_map[150][0] = -7 ;		mod_map[150][1] = 13 ;
		mod_map[151][0] = -7 ;		mod_map[151][1] = 15 ;
		mod_map[152][0] = -3 ;		mod_map[152][1] = 11 ;
		mod_map[153][0] = -3 ;		mod_map[153][1] = 9 ;
		mod_map[154][0] = -1 ;		mod_map[154][1] = 11 ;
		mod_map[155][0] = -1 ;		mod_map[155][1] = 9 ;
		mod_map[156][0] = -3 ;		mod_map[156][1] = 13 ;
		mod_map[157][0] = -3 ;		mod_map[157][1] = 15 ;
		mod_map[158][0] = -1 ;		mod_map[158][1] = 13 ;
		mod_map[159][0] = -1 ;		mod_map[159][1] = 15 ;
		mod_map[160][0] = -11 ;		mod_map[160][1] = 5 ;
		mod_map[161][0] = -11 ;		mod_map[161][1] = 7 ;
		mod_map[162][0] = -9 ;		mod_map[162][1] = 5 ;
		mod_map[163][0] = -9 ;		mod_map[163][1] = 7 ;
		mod_map[164][0] = -11 ;		mod_map[164][1] = 3 ;
		mod_map[165][0] = -11 ;		mod_map[165][1] = 1 ;
		mod_map[166][0] = -9 ;		mod_map[166][1] = 3 ;
		mod_map[167][0] = -9 ;		mod_map[167][1] = 1 ;
		mod_map[168][0] = -13 ;		mod_map[168][1] = 5 ;
		mod_map[169][0] = -13 ;		mod_map[169][1] = 7 ;
		mod_map[170][0] = -15 ;		mod_map[170][1] = 5 ;
		mod_map[171][0] = -15 ;		mod_map[171][1] = 7 ;
		mod_map[172][0] = -13 ;		mod_map[172][1] = 3 ;
		mod_map[173][0] = -13 ;		mod_map[173][1] = 1 ;
		mod_map[174][0] = -15 ;		mod_map[174][1] = 3 ;
		mod_map[175][0] = -15;		mod_map[175][1] = 1 ;
		mod_map[176][0] = -11 ;		mod_map[176][1] = 11 ;
		mod_map[177][0] = -11 ;		mod_map[177][1] = 9 ;
		mod_map[178][0] = -9 ;		mod_map[178][1] = 11 ;
		mod_map[179][0] = -9 ;		mod_map[179][1] = 9 ;
		mod_map[180][0] = -11 ;		mod_map[180][1] = 13 ;
		mod_map[181][0] = -11 ;		mod_map[181][1] = 15 ;
		mod_map[182][0] = -9 ;		mod_map[182][1] = 13 ;
		mod_map[183][0] = -9 ;		mod_map[183][1] = 15 ;
		mod_map[184][0] = -13 ;		mod_map[184][1] = 11 ;
		mod_map[185][0] = -13 ;		mod_map[185][1] = 9 ;
		mod_map[186][0] = -15 ;		mod_map[186][1] = 11 ;
		mod_map[187][0] = -15 ;		mod_map[187][1] = 9 ;
		mod_map[188][0] = -13 ;		mod_map[188][1] = 13 ;
		mod_map[189][0] = -13 ;		mod_map[189][1] = 15 ;
		mod_map[190][0] = -15 ;		mod_map[190][1] = 13 ;
		mod_map[191][0] = -15 ;		mod_map[191][1] = 15 ;
		mod_map[192][0] = -5 ;		mod_map[192][1] = -5 ;
		mod_map[193][0] = -5 ;		mod_map[193][1] = -7 ;
		mod_map[194][0] = -7 ;		mod_map[194][1] = -5 ;
		mod_map[195][0] = -7 ;		mod_map[195][1] = -7 ;
		mod_map[196][0] = -5 ;		mod_map[196][1] = -3 ;
		mod_map[197][0] = -5 ;		mod_map[197][1] = -1 ;
		mod_map[198][0] = -7 ;		mod_map[198][1] = -3 ;
		mod_map[199][0] = -7 ;		mod_map[199][1] = -1 ;
		mod_map[200][0] = -3 ;		mod_map[200][1] = -5 ;
		mod_map[201][0] = -3 ;		mod_map[201][1] = -7 ;
		mod_map[202][0] = -1 ;		mod_map[202][1] = -5 ;
		mod_map[203][0] = -1 ;		mod_map[203][1] = -7 ;
		mod_map[204][0] = -3 ;		mod_map[204][1] = -3 ;
		mod_map[205][0] = -3 ;		mod_map[205][1] = -1 ;
		mod_map[206][0] = -1 ;		mod_map[206][1] = -3 ;
		mod_map[207][0] = -1 ;		mod_map[207][1] = -1 ;
		mod_map[208][0] = -5 ;		mod_map[208][1] = -11 ;
		mod_map[209][0] = -5 ;		mod_map[209][1] = -9 ;
		mod_map[210][0] = -7 ;		mod_map[210][1] = -11 ;
		mod_map[211][0] = -7 ;		mod_map[211][1] = -9 ;
		mod_map[212][0] = -5 ;		mod_map[212][1] = -13 ;
		mod_map[213][0] = -5 ;		mod_map[213][1] = -15 ;
		mod_map[214][0] = -7 ;		mod_map[214][1] = -13 ;
		mod_map[215][0] = -7 ;		mod_map[215][1] = -15 ;
		mod_map[216][0] = -3 ;		mod_map[216][1] = -11 ;
		mod_map[217][0] = -3 ;		mod_map[217][1] = -9 ;
		mod_map[218][0] = -1 ;		mod_map[218][1] = -11 ;
		mod_map[219][0] = -1 ;		mod_map[219][1] = -9 ;
		mod_map[220][0] = -3 ;		mod_map[220][1] = -13 ;
		mod_map[221][0] = -3 ;		mod_map[221][1] = -15 ;
		mod_map[222][0] = -1 ;		mod_map[222][1] = -13 ;
		mod_map[223][0] = -1 ;		mod_map[223][1] = -15 ;
		mod_map[224][0] = -11 ;		mod_map[224][1] = -5 ;
		mod_map[225][0] = -11 ;		mod_map[225][1] = -7 ;
		mod_map[226][0] = -9 ;		mod_map[226][1] = -5 ;
		mod_map[227][0] = -9 ;		mod_map[227][1] = -7 ;
		mod_map[228][0] = -11 ;		mod_map[228][1] = -3 ;
		mod_map[229][0] = -11 ;		mod_map[229][1] = -1 ;
		mod_map[230][0] = -9 ;		mod_map[230][1] = -3 ;
		mod_map[231][0] = -9 ;		mod_map[231][1] = -1 ;
		mod_map[232][0] = -13 ;		mod_map[232][1] = -5 ;
		mod_map[233][0] = -13 ;		mod_map[233][1] = -7 ;
		mod_map[234][0] = -15 ;		mod_map[234][1] = -5 ;
		mod_map[235][0] = -15 ;		mod_map[235][1] = -7 ;
		mod_map[236][0] = -13 ;		mod_map[236][1] = -3 ;
		mod_map[237][0] = -13 ;		mod_map[237][1] = -1 ;
		mod_map[238][0] = -15 ;		mod_map[238][1] = -3 ;
		mod_map[239][0] = -15 ;		mod_map[239][1] = -1 ;
		mod_map[240][0] = -11 ;		mod_map[240][1] = -11 ;
		mod_map[241][0] = -11 ;		mod_map[241][1] = -9 ;
		mod_map[242][0] = -9 ;		mod_map[242][1] = -11 ;
		mod_map[243][0] = -9 ;		mod_map[243][1] = -9 ;
		mod_map[244][0] = -11 ;		mod_map[244][1] = -13 ;
		mod_map[245][0] = -11 ;		mod_map[245][1] = -15 ;
		mod_map[246][0] = -9 ;		mod_map[246][1] = -13 ;
		mod_map[247][0] = -9 ;		mod_map[247][1] = -15 ;
		mod_map[248][0] = -13 ;		mod_map[248][1] = -11 ;
		mod_map[249][0] = -13 ;		mod_map[249][1] = -9 ;
		mod_map[250][0] = -15 ;		mod_map[250][1] = -11 ;
		mod_map[251][0] = -15 ;		mod_map[251][1] = -9 ;
		mod_map[252][0] = -13 ;		mod_map[252][1] = -13 ;
		mod_map[253][0] = -13 ;		mod_map[253][1] = -15 ;
		mod_map[254][0] = -15 ;		mod_map[254][1] = -13 ;
		mod_map[255][0] = -15 ;		mod_map[255][1] = -15 ;


		for (int i=0; i<M_ary; i++)
			for (int j=0; j<2; j++)
				mod_map[i][j] /= sqrt(170.0) ;
	}


	
	return mod_map ;
}


// Rx side
void QAM_Modulation::Func_rx_LLR_Limitter ()
{
	for (int i=0; i<log2M*SymbolLength; i++) {
		if (LLR[i] > LLR_threshold)
			LLR[i] = LLR_threshold ;
		if (LLR[i] < -LLR_threshold)
			LLR[i] = -LLR_threshold ;
	}
}


void QAM_Modulation::Func_rx_Demod_M2 () 
{
	// LLR calculation for BPSK
	for (int i=0; i<SymbolLength; i++)
		LLR[i] = log ( pdf[i][1] / pdf[i][0] ) ;
}


void QAM_Modulation::Func_rx_Demod_M4 () 
{

	// LLR calculation for QPSK
	for (int i=0; i<SymbolLength; i++) {
		LLR[2*i]   = log ( (pdf[i][2]+pdf[i][3]) / (pdf[i][0]+pdf[i][1]) ) ;
		LLR[2*i+1] = log ( (pdf[i][1]+pdf[i][3]) / (pdf[i][0]+pdf[i][2]) ) ;
	}

}


void QAM_Modulation::Func_rx_Demod_M16 () 
{
	// LLR calculation for 16QAM
	for (int i=0; i<SymbolLength; i++) {
		LLR[4*i]   = log ( (pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]) ) ;

		LLR[4*i+1] = log ( (pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]) ) ;

		LLR[4*i+2] = log ( (pdf[i][2]+pdf[i][3]+pdf[i][6]+pdf[i][7]+pdf[i][10]+pdf[i][11]+pdf[i][14]+pdf[i][15]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][4]+pdf[i][5]+pdf[i][8]+pdf[i][9]+pdf[i][12]+pdf[i][13]) ) ;

		LLR[4*i+3] = log ( (pdf[i][1]+pdf[i][3]+pdf[i][5]+pdf[i][7]+pdf[i][9]+pdf[i][11]+pdf[i][13]+pdf[i][15]) 
						  / (pdf[i][0]+pdf[i][2]+pdf[i][4]+pdf[i][6]+pdf[i][8]+pdf[i][10]+pdf[i][12]+pdf[i][14]) ) ;
	}
}


void QAM_Modulation::Func_rx_Demod_M64 () 
{
	// LLR calculation for 64QAM
	for (int i=0; i<SymbolLength; i++) {
		LLR[6*i]   = log ( (pdf[i][32]+pdf[i][33]+pdf[i][34]+pdf[i][35]+pdf[i][36]+pdf[i][37]+pdf[i][38]+pdf[i][39]
							+pdf[i][40]+pdf[i][41]+pdf[i][42]+pdf[i][43]+pdf[i][44]+pdf[i][45]+pdf[i][46]+pdf[i][47]
							+pdf[i][48]+pdf[i][49]+pdf[i][50]+pdf[i][51]+pdf[i][52]+pdf[i][53]+pdf[i][54]+pdf[i][55]
							+pdf[i][56]+pdf[i][57]+pdf[i][58]+pdf[i][59]+pdf[i][60]+pdf[i][61]+pdf[i][62]+pdf[i][63]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]
							+pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]
							+pdf[i][16]+pdf[i][17]+pdf[i][18]+pdf[i][19]+pdf[i][20]+pdf[i][21]+pdf[i][22]+pdf[i][23]
							+pdf[i][24]+pdf[i][25]+pdf[i][26]+pdf[i][27]+pdf[i][28]+pdf[i][29]+pdf[i][30]+pdf[i][31]) ) ;

		LLR[6*i+1] = log ( (pdf[i][16]+pdf[i][17]+pdf[i][18]+pdf[i][19]+pdf[i][20]+pdf[i][21]+pdf[i][22]+pdf[i][23]
							+pdf[i][24]+pdf[i][25]+pdf[i][26]+pdf[i][27]+pdf[i][28]+pdf[i][29]+pdf[i][30]+pdf[i][31]
							+pdf[i][48]+pdf[i][49]+pdf[i][50]+pdf[i][51]+pdf[i][52]+pdf[i][53]+pdf[i][54]+pdf[i][55]
							+pdf[i][56]+pdf[i][57]+pdf[i][58]+pdf[i][59]+pdf[i][60]+pdf[i][61]+pdf[i][62]+pdf[i][63]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]
							+pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]
							+pdf[i][32]+pdf[i][33]+pdf[i][34]+pdf[i][35]+pdf[i][36]+pdf[i][37]+pdf[i][38]+pdf[i][39]
							+pdf[i][40]+pdf[i][41]+pdf[i][42]+pdf[i][43]+pdf[i][44]+pdf[i][45]+pdf[i][46]+pdf[i][47]) ) ;

		LLR[6*i+2] = log ( (pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]
							+pdf[i][24]+pdf[i][25]+pdf[i][26]+pdf[i][27]+pdf[i][28]+pdf[i][29]+pdf[i][30]+pdf[i][31]
							+pdf[i][40]+pdf[i][41]+pdf[i][42]+pdf[i][43]+pdf[i][44]+pdf[i][45]+pdf[i][46]+pdf[i][47]
							+pdf[i][56]+pdf[i][57]+pdf[i][58]+pdf[i][59]+pdf[i][60]+pdf[i][61]+pdf[i][62]+pdf[i][63]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]
							+pdf[i][16]+pdf[i][17]+pdf[i][18]+pdf[i][19]+pdf[i][20]+pdf[i][21]+pdf[i][22]+pdf[i][23]
							+pdf[i][32]+pdf[i][33]+pdf[i][34]+pdf[i][35]+pdf[i][36]+pdf[i][37]+pdf[i][38]+pdf[i][39]
							+pdf[i][48]+pdf[i][49]+pdf[i][50]+pdf[i][51]+pdf[i][52]+pdf[i][53]+pdf[i][54]+pdf[i][55]) ) ;

		LLR[6*i+3] = log ( (pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]
							+pdf[i][20]+pdf[i][21]+pdf[i][22]+pdf[i][23]+pdf[i][28]+pdf[i][29]+pdf[i][30]+pdf[i][31]
							+pdf[i][36]+pdf[i][37]+pdf[i][38]+pdf[i][39]+pdf[i][44]+pdf[i][45]+pdf[i][46]+pdf[i][47]
							+pdf[i][52]+pdf[i][53]+pdf[i][54]+pdf[i][55]+pdf[i][60]+pdf[i][61]+pdf[i][62]+pdf[i][63]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]
							+pdf[i][16]+pdf[i][17]+pdf[i][18]+pdf[i][19]+pdf[i][24]+pdf[i][25]+pdf[i][26]+pdf[i][27]
							+pdf[i][32]+pdf[i][33]+pdf[i][34]+pdf[i][35]+pdf[i][40]+pdf[i][41]+pdf[i][42]+pdf[i][43]
							+pdf[i][48]+pdf[i][49]+pdf[i][50]+pdf[i][51]+pdf[i][56]+pdf[i][57]+pdf[i][58]+pdf[i][59]) ) ;

		LLR[6*i+4] = log ( (pdf[i][2]+pdf[i][3]+pdf[i][6]+pdf[i][7]+pdf[i][10]+pdf[i][11]+pdf[i][14]+pdf[i][15]
							+pdf[i][18]+pdf[i][19]+pdf[i][22]+pdf[i][23]+pdf[i][26]+pdf[i][27]+pdf[i][30]+pdf[i][31]
							+pdf[i][34]+pdf[i][35]+pdf[i][38]+pdf[i][39]+pdf[i][42]+pdf[i][43]+pdf[i][46]+pdf[i][47]
							+pdf[i][50]+pdf[i][51]+pdf[i][54]+pdf[i][55]+pdf[i][58]+pdf[i][59]+pdf[i][62]+pdf[i][63]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][4]+pdf[i][5]+pdf[i][8]+pdf[i][9]+pdf[i][12]+pdf[i][13]
							+pdf[i][16]+pdf[i][17]+pdf[i][20]+pdf[i][21]+pdf[i][24]+pdf[i][25]+pdf[i][28]+pdf[i][29]
							+pdf[i][32]+pdf[i][33]+pdf[i][36]+pdf[i][37]+pdf[i][40]+pdf[i][41]+pdf[i][44]+pdf[i][45]
							+pdf[i][48]+pdf[i][49]+pdf[i][52]+pdf[i][53]+pdf[i][56]+pdf[i][57]+pdf[i][60]+pdf[i][61]) ) ;

		LLR[6*i+5] = log ( (pdf[i][1]+pdf[i][3]+pdf[i][5]+pdf[i][7]+pdf[i][9]+pdf[i][11]+pdf[i][13]+pdf[i][15]
							+pdf[i][17]+pdf[i][19]+pdf[i][21]+pdf[i][23]+pdf[i][25]+pdf[i][27]+pdf[i][29]+pdf[i][31]
							+pdf[i][33]+pdf[i][35]+pdf[i][37]+pdf[i][39]+pdf[i][41]+pdf[i][43]+pdf[i][45]+pdf[i][47]
							+pdf[i][49]+pdf[i][51]+pdf[i][53]+pdf[i][55]+pdf[i][57]+pdf[i][59]+pdf[i][61]+pdf[i][63]) 
						  / (pdf[i][0]+pdf[i][2]+pdf[i][4]+pdf[i][6]+pdf[i][8]+pdf[i][10]+pdf[i][12]+pdf[i][14]
							+pdf[i][16]+pdf[i][18]+pdf[i][20]+pdf[i][22]+pdf[i][24]+pdf[i][26]+pdf[i][28]+pdf[i][30]
							+pdf[i][32]+pdf[i][34]+pdf[i][36]+pdf[i][38]+pdf[i][40]+pdf[i][42]+pdf[i][44]+pdf[i][46]
							+pdf[i][48]+pdf[i][50]+pdf[i][52]+pdf[i][54]+pdf[i][56]+pdf[i][58]+pdf[i][60]+pdf[i][62]) ) ;
	//cout<< LLR[6*i+0]<<" "<<LLR[6*i+1]<<" "<<LLR[6*i+2]<<" "<<LLR[6*i+3]<<" "<<LLR[6*i+4]<<" "<<LLR[6*i+5]<<"\n";
	}
}


void QAM_Modulation::Func_rx_Demod_M256 () 
{
	// LLR calculation for 256QAM
	for (int i=0; i<SymbolLength; i++) {
		LLR[8*i]   = log ( (pdf[i][128]+pdf[i][129]+pdf[i][130]+pdf[i][131]+pdf[i][132]+pdf[i][133]+pdf[i][134]+pdf[i][135]
							+pdf[i][136]+pdf[i][137]+pdf[i][138]+pdf[i][139]+pdf[i][140]+pdf[i][141]+pdf[i][142]+pdf[i][143]
							+pdf[i][144]+pdf[i][145]+pdf[i][146]+pdf[i][147]+pdf[i][148]+pdf[i][149]+pdf[i][150]+pdf[i][151]
							+pdf[i][152]+pdf[i][153]+pdf[i][154]+pdf[i][155]+pdf[i][156]+pdf[i][157]+pdf[i][158]+pdf[i][159]
							+pdf[i][160]+pdf[i][161]+pdf[i][162]+pdf[i][163]+pdf[i][164]+pdf[i][165]+pdf[i][166]+pdf[i][167]
							+pdf[i][168]+pdf[i][169]+pdf[i][170]+pdf[i][171]+pdf[i][172]+pdf[i][173]+pdf[i][174]+pdf[i][175]
							+pdf[i][176]+pdf[i][177]+pdf[i][178]+pdf[i][179]+pdf[i][180]+pdf[i][181]+pdf[i][182]+pdf[i][183]
							+pdf[i][184]+pdf[i][185]+pdf[i][186]+pdf[i][187]+pdf[i][188]+pdf[i][189]+pdf[i][190]+pdf[i][191]
							+pdf[i][192]+pdf[i][193]+pdf[i][194]+pdf[i][195]+pdf[i][196]+pdf[i][197]+pdf[i][198]+pdf[i][199]
							+pdf[i][200]+pdf[i][201]+pdf[i][202]+pdf[i][203]+pdf[i][204]+pdf[i][205]+pdf[i][206]+pdf[i][207]
							+pdf[i][208]+pdf[i][209]+pdf[i][210]+pdf[i][211]+pdf[i][212]+pdf[i][213]+pdf[i][214]+pdf[i][215]
							+pdf[i][216]+pdf[i][217]+pdf[i][218]+pdf[i][219]+pdf[i][220]+pdf[i][221]+pdf[i][222]+pdf[i][223]
							+pdf[i][224]+pdf[i][225]+pdf[i][226]+pdf[i][227]+pdf[i][228]+pdf[i][229]+pdf[i][230]+pdf[i][231]
							+pdf[i][232]+pdf[i][233]+pdf[i][234]+pdf[i][235]+pdf[i][236]+pdf[i][237]+pdf[i][238]+pdf[i][239]
							+pdf[i][240]+pdf[i][241]+pdf[i][242]+pdf[i][243]+pdf[i][244]+pdf[i][245]+pdf[i][246]+pdf[i][247]
							+pdf[i][248]+pdf[i][249]+pdf[i][250]+pdf[i][251]+pdf[i][252]+pdf[i][253]+pdf[i][254]+pdf[i][255]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]
							+pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]
							+pdf[i][16]+pdf[i][17]+pdf[i][18]+pdf[i][19]+pdf[i][20]+pdf[i][21]+pdf[i][22]+pdf[i][23]
							+pdf[i][24]+pdf[i][25]+pdf[i][26]+pdf[i][27]+pdf[i][28]+pdf[i][29]+pdf[i][30]+pdf[i][31]
							+pdf[i][32]+pdf[i][33]+pdf[i][34]+pdf[i][35]+pdf[i][36]+pdf[i][37]+pdf[i][38]+pdf[i][39]
							+pdf[i][40]+pdf[i][41]+pdf[i][42]+pdf[i][43]+pdf[i][44]+pdf[i][45]+pdf[i][46]+pdf[i][47]
							+pdf[i][48]+pdf[i][49]+pdf[i][50]+pdf[i][51]+pdf[i][52]+pdf[i][53]+pdf[i][54]+pdf[i][55]
							+pdf[i][56]+pdf[i][57]+pdf[i][58]+pdf[i][59]+pdf[i][60]+pdf[i][61]+pdf[i][62]+pdf[i][63]
							+pdf[i][64]+pdf[i][65]+pdf[i][66]+pdf[i][67]+pdf[i][68]+pdf[i][69]+pdf[i][70]+pdf[i][71]
							+pdf[i][72]+pdf[i][73]+pdf[i][74]+pdf[i][75]+pdf[i][76]+pdf[i][77]+pdf[i][78]+pdf[i][79]
							+pdf[i][80]+pdf[i][81]+pdf[i][82]+pdf[i][83]+pdf[i][84]+pdf[i][85]+pdf[i][86]+pdf[i][87]
							+pdf[i][88]+pdf[i][89]+pdf[i][90]+pdf[i][91]+pdf[i][92]+pdf[i][93]+pdf[i][94]+pdf[i][95]
							+pdf[i][96]+pdf[i][97]+pdf[i][98]+pdf[i][99]+pdf[i][100]+pdf[i][101]+pdf[i][102]+pdf[i][103]
							+pdf[i][104]+pdf[i][105]+pdf[i][106]+pdf[i][107]+pdf[i][108]+pdf[i][109]+pdf[i][110]+pdf[i][111]
							+pdf[i][112]+pdf[i][113]+pdf[i][114]+pdf[i][115]+pdf[i][116]+pdf[i][117]+pdf[i][118]+pdf[i][119]
							+pdf[i][120]+pdf[i][121]+pdf[i][122]+pdf[i][123]+pdf[i][124]+pdf[i][125]+pdf[i][126]+pdf[i][127]) ) ;

		LLR[8*i+1] = log ( (pdf[i][64]+pdf[i][65]+pdf[i][66]+pdf[i][67]+pdf[i][68]+pdf[i][69]+pdf[i][70]+pdf[i][71]
							+pdf[i][72]+pdf[i][73]+pdf[i][74]+pdf[i][75]+pdf[i][76]+pdf[i][77]+pdf[i][78]+pdf[i][79]
							+pdf[i][80]+pdf[i][81]+pdf[i][82]+pdf[i][83]+pdf[i][84]+pdf[i][85]+pdf[i][86]+pdf[i][87]
							+pdf[i][88]+pdf[i][89]+pdf[i][90]+pdf[i][91]+pdf[i][92]+pdf[i][93]+pdf[i][94]+pdf[i][95]
							+pdf[i][96]+pdf[i][97]+pdf[i][98]+pdf[i][99]+pdf[i][100]+pdf[i][101]+pdf[i][102]+pdf[i][103]
							+pdf[i][104]+pdf[i][105]+pdf[i][106]+pdf[i][107]+pdf[i][108]+pdf[i][109]+pdf[i][110]+pdf[i][111]
							+pdf[i][112]+pdf[i][113]+pdf[i][114]+pdf[i][115]+pdf[i][116]+pdf[i][117]+pdf[i][118]+pdf[i][119]
							+pdf[i][120]+pdf[i][121]+pdf[i][122]+pdf[i][123]+pdf[i][124]+pdf[i][125]+pdf[i][126]+pdf[i][127]
							+pdf[i][192]+pdf[i][193]+pdf[i][194]+pdf[i][195]+pdf[i][196]+pdf[i][197]+pdf[i][198]+pdf[i][199]
							+pdf[i][200]+pdf[i][201]+pdf[i][202]+pdf[i][203]+pdf[i][204]+pdf[i][205]+pdf[i][206]+pdf[i][207]
							+pdf[i][208]+pdf[i][209]+pdf[i][210]+pdf[i][211]+pdf[i][212]+pdf[i][213]+pdf[i][214]+pdf[i][215]
							+pdf[i][216]+pdf[i][217]+pdf[i][218]+pdf[i][219]+pdf[i][220]+pdf[i][221]+pdf[i][222]+pdf[i][223]
							+pdf[i][224]+pdf[i][225]+pdf[i][226]+pdf[i][227]+pdf[i][228]+pdf[i][229]+pdf[i][230]+pdf[i][231]
							+pdf[i][232]+pdf[i][233]+pdf[i][234]+pdf[i][235]+pdf[i][236]+pdf[i][237]+pdf[i][238]+pdf[i][239]
							+pdf[i][240]+pdf[i][241]+pdf[i][242]+pdf[i][243]+pdf[i][244]+pdf[i][245]+pdf[i][246]+pdf[i][247]
							+pdf[i][248]+pdf[i][249]+pdf[i][250]+pdf[i][251]+pdf[i][252]+pdf[i][253]+pdf[i][254]+pdf[i][255]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]
							+pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]
							+pdf[i][16]+pdf[i][17]+pdf[i][18]+pdf[i][19]+pdf[i][20]+pdf[i][21]+pdf[i][22]+pdf[i][23]
							+pdf[i][24]+pdf[i][25]+pdf[i][26]+pdf[i][27]+pdf[i][28]+pdf[i][29]+pdf[i][30]+pdf[i][31]
							+pdf[i][32]+pdf[i][33]+pdf[i][34]+pdf[i][35]+pdf[i][36]+pdf[i][37]+pdf[i][38]+pdf[i][39]
							+pdf[i][40]+pdf[i][41]+pdf[i][42]+pdf[i][43]+pdf[i][44]+pdf[i][45]+pdf[i][46]+pdf[i][47]
							+pdf[i][48]+pdf[i][49]+pdf[i][50]+pdf[i][51]+pdf[i][52]+pdf[i][53]+pdf[i][54]+pdf[i][55]
							+pdf[i][56]+pdf[i][57]+pdf[i][58]+pdf[i][59]+pdf[i][60]+pdf[i][61]+pdf[i][62]+pdf[i][63]
							+pdf[i][128]+pdf[i][129]+pdf[i][130]+pdf[i][131]+pdf[i][132]+pdf[i][133]+pdf[i][134]+pdf[i][135]
							+pdf[i][136]+pdf[i][137]+pdf[i][138]+pdf[i][139]+pdf[i][140]+pdf[i][141]+pdf[i][142]+pdf[i][143]
							+pdf[i][144]+pdf[i][145]+pdf[i][146]+pdf[i][147]+pdf[i][148]+pdf[i][149]+pdf[i][150]+pdf[i][151]
							+pdf[i][152]+pdf[i][153]+pdf[i][154]+pdf[i][155]+pdf[i][156]+pdf[i][157]+pdf[i][158]+pdf[i][159]
							+pdf[i][160]+pdf[i][161]+pdf[i][162]+pdf[i][163]+pdf[i][164]+pdf[i][165]+pdf[i][166]+pdf[i][167]
							+pdf[i][168]+pdf[i][169]+pdf[i][170]+pdf[i][171]+pdf[i][172]+pdf[i][173]+pdf[i][174]+pdf[i][175]
							+pdf[i][176]+pdf[i][177]+pdf[i][178]+pdf[i][179]+pdf[i][180]+pdf[i][181]+pdf[i][182]+pdf[i][183]
							+pdf[i][184]+pdf[i][185]+pdf[i][186]+pdf[i][187]+pdf[i][188]+pdf[i][189]+pdf[i][190]+pdf[i][191]) ) ;

		LLR[8*i+2] = log ( (pdf[i][32]+pdf[i][33]+pdf[i][34]+pdf[i][35]+pdf[i][36]+pdf[i][37]+pdf[i][38]+pdf[i][39]
							+pdf[i][40]+pdf[i][41]+pdf[i][42]+pdf[i][43]+pdf[i][44]+pdf[i][45]+pdf[i][46]+pdf[i][47]
							+pdf[i][48]+pdf[i][49]+pdf[i][50]+pdf[i][51]+pdf[i][52]+pdf[i][53]+pdf[i][54]+pdf[i][55]
							+pdf[i][56]+pdf[i][57]+pdf[i][58]+pdf[i][59]+pdf[i][60]+pdf[i][61]+pdf[i][62]+pdf[i][63]
							+pdf[i][96]+pdf[i][97]+pdf[i][98]+pdf[i][99]+pdf[i][100]+pdf[i][101]+pdf[i][102]+pdf[i][103]
							+pdf[i][104]+pdf[i][105]+pdf[i][106]+pdf[i][107]+pdf[i][108]+pdf[i][109]+pdf[i][110]+pdf[i][111]
							+pdf[i][112]+pdf[i][113]+pdf[i][114]+pdf[i][115]+pdf[i][116]+pdf[i][117]+pdf[i][118]+pdf[i][119]
							+pdf[i][120]+pdf[i][121]+pdf[i][122]+pdf[i][123]+pdf[i][124]+pdf[i][125]+pdf[i][126]+pdf[i][127]
							+pdf[i][160]+pdf[i][161]+pdf[i][162]+pdf[i][163]+pdf[i][164]+pdf[i][165]+pdf[i][166]+pdf[i][167]
							+pdf[i][168]+pdf[i][169]+pdf[i][170]+pdf[i][171]+pdf[i][172]+pdf[i][173]+pdf[i][174]+pdf[i][175]
							+pdf[i][176]+pdf[i][177]+pdf[i][178]+pdf[i][179]+pdf[i][180]+pdf[i][181]+pdf[i][182]+pdf[i][183]
							+pdf[i][184]+pdf[i][185]+pdf[i][186]+pdf[i][187]+pdf[i][188]+pdf[i][189]+pdf[i][190]+pdf[i][191]
							+pdf[i][224]+pdf[i][225]+pdf[i][226]+pdf[i][227]+pdf[i][228]+pdf[i][229]+pdf[i][230]+pdf[i][231]
							+pdf[i][232]+pdf[i][233]+pdf[i][234]+pdf[i][235]+pdf[i][236]+pdf[i][237]+pdf[i][238]+pdf[i][239]
							+pdf[i][240]+pdf[i][241]+pdf[i][242]+pdf[i][243]+pdf[i][244]+pdf[i][245]+pdf[i][246]+pdf[i][247]
							+pdf[i][248]+pdf[i][249]+pdf[i][250]+pdf[i][251]+pdf[i][252]+pdf[i][253]+pdf[i][254]+pdf[i][255]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]
							+pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]
							+pdf[i][16]+pdf[i][17]+pdf[i][18]+pdf[i][19]+pdf[i][20]+pdf[i][21]+pdf[i][22]+pdf[i][23]
							+pdf[i][24]+pdf[i][25]+pdf[i][26]+pdf[i][27]+pdf[i][28]+pdf[i][29]+pdf[i][30]+pdf[i][31]
							+pdf[i][64]+pdf[i][65]+pdf[i][66]+pdf[i][67]+pdf[i][68]+pdf[i][69]+pdf[i][70]+pdf[i][71]
							+pdf[i][72]+pdf[i][73]+pdf[i][74]+pdf[i][75]+pdf[i][76]+pdf[i][77]+pdf[i][78]+pdf[i][79]
							+pdf[i][80]+pdf[i][81]+pdf[i][82]+pdf[i][83]+pdf[i][84]+pdf[i][85]+pdf[i][86]+pdf[i][87]
							+pdf[i][88]+pdf[i][89]+pdf[i][90]+pdf[i][91]+pdf[i][92]+pdf[i][93]+pdf[i][94]+pdf[i][95]
							+pdf[i][128]+pdf[i][129]+pdf[i][130]+pdf[i][131]+pdf[i][132]+pdf[i][133]+pdf[i][134]+pdf[i][135]
							+pdf[i][136]+pdf[i][137]+pdf[i][138]+pdf[i][139]+pdf[i][140]+pdf[i][141]+pdf[i][142]+pdf[i][143]
							+pdf[i][144]+pdf[i][145]+pdf[i][146]+pdf[i][147]+pdf[i][148]+pdf[i][149]+pdf[i][150]+pdf[i][151]
							+pdf[i][152]+pdf[i][153]+pdf[i][154]+pdf[i][155]+pdf[i][156]+pdf[i][157]+pdf[i][158]+pdf[i][159]
							+pdf[i][192]+pdf[i][193]+pdf[i][194]+pdf[i][195]+pdf[i][196]+pdf[i][197]+pdf[i][198]+pdf[i][199]
							+pdf[i][200]+pdf[i][201]+pdf[i][202]+pdf[i][203]+pdf[i][204]+pdf[i][205]+pdf[i][206]+pdf[i][207]
							+pdf[i][208]+pdf[i][209]+pdf[i][210]+pdf[i][211]+pdf[i][212]+pdf[i][213]+pdf[i][214]+pdf[i][215]
							+pdf[i][216]+pdf[i][217]+pdf[i][218]+pdf[i][219]+pdf[i][220]+pdf[i][221]+pdf[i][222]+pdf[i][223]) ) ;

		LLR[8*i+3] = log ( (pdf[i][16]+pdf[i][17]+pdf[i][18]+pdf[i][19]+pdf[i][20]+pdf[i][21]+pdf[i][22]+pdf[i][23]
							+pdf[i][24]+pdf[i][25]+pdf[i][26]+pdf[i][27]+pdf[i][28]+pdf[i][29]+pdf[i][30]+pdf[i][31]
							+pdf[i][48]+pdf[i][49]+pdf[i][50]+pdf[i][51]+pdf[i][52]+pdf[i][53]+pdf[i][54]+pdf[i][55]
							+pdf[i][56]+pdf[i][57]+pdf[i][58]+pdf[i][59]+pdf[i][60]+pdf[i][61]+pdf[i][62]+pdf[i][63]
							+pdf[i][80]+pdf[i][81]+pdf[i][82]+pdf[i][83]+pdf[i][84]+pdf[i][85]+pdf[i][86]+pdf[i][87]
							+pdf[i][88]+pdf[i][89]+pdf[i][90]+pdf[i][91]+pdf[i][92]+pdf[i][93]+pdf[i][94]+pdf[i][95]
							+pdf[i][112]+pdf[i][113]+pdf[i][114]+pdf[i][115]+pdf[i][116]+pdf[i][117]+pdf[i][118]+pdf[i][119]
							+pdf[i][120]+pdf[i][121]+pdf[i][122]+pdf[i][123]+pdf[i][124]+pdf[i][125]+pdf[i][126]+pdf[i][127]
							+pdf[i][144]+pdf[i][145]+pdf[i][146]+pdf[i][147]+pdf[i][148]+pdf[i][149]+pdf[i][150]+pdf[i][151]
							+pdf[i][152]+pdf[i][153]+pdf[i][154]+pdf[i][155]+pdf[i][156]+pdf[i][157]+pdf[i][158]+pdf[i][159]
							+pdf[i][176]+pdf[i][177]+pdf[i][178]+pdf[i][179]+pdf[i][180]+pdf[i][181]+pdf[i][182]+pdf[i][183]
							+pdf[i][184]+pdf[i][185]+pdf[i][186]+pdf[i][187]+pdf[i][188]+pdf[i][189]+pdf[i][190]+pdf[i][191]
							+pdf[i][208]+pdf[i][209]+pdf[i][210]+pdf[i][211]+pdf[i][212]+pdf[i][213]+pdf[i][214]+pdf[i][215]
							+pdf[i][216]+pdf[i][217]+pdf[i][218]+pdf[i][219]+pdf[i][220]+pdf[i][221]+pdf[i][222]+pdf[i][223]
							+pdf[i][240]+pdf[i][241]+pdf[i][242]+pdf[i][243]+pdf[i][244]+pdf[i][245]+pdf[i][246]+pdf[i][247]
							+pdf[i][248]+pdf[i][249]+pdf[i][250]+pdf[i][251]+pdf[i][252]+pdf[i][253]+pdf[i][254]+pdf[i][255]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]
							+pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]
							+pdf[i][32]+pdf[i][33]+pdf[i][34]+pdf[i][35]+pdf[i][36]+pdf[i][37]+pdf[i][38]+pdf[i][39]
							+pdf[i][40]+pdf[i][41]+pdf[i][42]+pdf[i][43]+pdf[i][44]+pdf[i][45]+pdf[i][46]+pdf[i][47]
							+pdf[i][64]+pdf[i][65]+pdf[i][66]+pdf[i][67]+pdf[i][68]+pdf[i][69]+pdf[i][70]+pdf[i][71]
							+pdf[i][72]+pdf[i][73]+pdf[i][74]+pdf[i][75]+pdf[i][76]+pdf[i][77]+pdf[i][78]+pdf[i][79]
							+pdf[i][96]+pdf[i][97]+pdf[i][98]+pdf[i][99]+pdf[i][100]+pdf[i][101]+pdf[i][102]+pdf[i][103]
							+pdf[i][104]+pdf[i][105]+pdf[i][106]+pdf[i][107]+pdf[i][108]+pdf[i][109]+pdf[i][110]+pdf[i][111]
							+pdf[i][128]+pdf[i][129]+pdf[i][130]+pdf[i][131]+pdf[i][132]+pdf[i][133]+pdf[i][134]+pdf[i][135]
							+pdf[i][136]+pdf[i][137]+pdf[i][138]+pdf[i][139]+pdf[i][140]+pdf[i][141]+pdf[i][142]+pdf[i][143]
							+pdf[i][160]+pdf[i][161]+pdf[i][162]+pdf[i][163]+pdf[i][164]+pdf[i][165]+pdf[i][166]+pdf[i][167]
							+pdf[i][168]+pdf[i][169]+pdf[i][170]+pdf[i][171]+pdf[i][172]+pdf[i][173]+pdf[i][174]+pdf[i][175]
							+pdf[i][192]+pdf[i][193]+pdf[i][194]+pdf[i][195]+pdf[i][196]+pdf[i][197]+pdf[i][198]+pdf[i][199]
							+pdf[i][200]+pdf[i][201]+pdf[i][202]+pdf[i][203]+pdf[i][204]+pdf[i][205]+pdf[i][206]+pdf[i][207]
							+pdf[i][224]+pdf[i][225]+pdf[i][226]+pdf[i][227]+pdf[i][228]+pdf[i][229]+pdf[i][230]+pdf[i][231]
							+pdf[i][232]+pdf[i][233]+pdf[i][234]+pdf[i][235]+pdf[i][236]+pdf[i][237]+pdf[i][238]+pdf[i][239]) ) ;

		LLR[8*i+4] = log ( (pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]
							+pdf[i][24]+pdf[i][25]+pdf[i][26]+pdf[i][27]+pdf[i][28]+pdf[i][29]+pdf[i][30]+pdf[i][31]
							+pdf[i][40]+pdf[i][41]+pdf[i][42]+pdf[i][43]+pdf[i][44]+pdf[i][45]+pdf[i][46]+pdf[i][47]
							+pdf[i][56]+pdf[i][57]+pdf[i][58]+pdf[i][59]+pdf[i][60]+pdf[i][61]+pdf[i][62]+pdf[i][63]
							+pdf[i][72]+pdf[i][73]+pdf[i][74]+pdf[i][75]+pdf[i][76]+pdf[i][77]+pdf[i][78]+pdf[i][79]
							+pdf[i][88]+pdf[i][89]+pdf[i][90]+pdf[i][91]+pdf[i][92]+pdf[i][93]+pdf[i][94]+pdf[i][95]
							+pdf[i][104]+pdf[i][105]+pdf[i][106]+pdf[i][107]+pdf[i][108]+pdf[i][109]+pdf[i][110]+pdf[i][111]
							+pdf[i][120]+pdf[i][121]+pdf[i][122]+pdf[i][123]+pdf[i][124]+pdf[i][125]+pdf[i][126]+pdf[i][127]
							+pdf[i][136]+pdf[i][137]+pdf[i][138]+pdf[i][139]+pdf[i][140]+pdf[i][141]+pdf[i][142]+pdf[i][143]
							+pdf[i][152]+pdf[i][153]+pdf[i][154]+pdf[i][155]+pdf[i][156]+pdf[i][157]+pdf[i][158]+pdf[i][159]
							+pdf[i][168]+pdf[i][169]+pdf[i][170]+pdf[i][171]+pdf[i][172]+pdf[i][173]+pdf[i][174]+pdf[i][175]
							+pdf[i][184]+pdf[i][185]+pdf[i][186]+pdf[i][187]+pdf[i][188]+pdf[i][189]+pdf[i][190]+pdf[i][191]
							+pdf[i][200]+pdf[i][201]+pdf[i][202]+pdf[i][203]+pdf[i][204]+pdf[i][205]+pdf[i][206]+pdf[i][207]
							+pdf[i][216]+pdf[i][217]+pdf[i][218]+pdf[i][219]+pdf[i][220]+pdf[i][221]+pdf[i][222]+pdf[i][223]
							+pdf[i][232]+pdf[i][233]+pdf[i][234]+pdf[i][235]+pdf[i][236]+pdf[i][237]+pdf[i][238]+pdf[i][239]
							+pdf[i][248]+pdf[i][249]+pdf[i][250]+pdf[i][251]+pdf[i][252]+pdf[i][253]+pdf[i][254]+pdf[i][255]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]
							+pdf[i][16]+pdf[i][17]+pdf[i][18]+pdf[i][19]+pdf[i][20]+pdf[i][21]+pdf[i][22]+pdf[i][23]
							+pdf[i][32]+pdf[i][33]+pdf[i][34]+pdf[i][35]+pdf[i][36]+pdf[i][37]+pdf[i][38]+pdf[i][39]
							+pdf[i][48]+pdf[i][49]+pdf[i][50]+pdf[i][51]+pdf[i][52]+pdf[i][53]+pdf[i][54]+pdf[i][55]
							+pdf[i][64]+pdf[i][65]+pdf[i][66]+pdf[i][67]+pdf[i][68]+pdf[i][69]+pdf[i][70]+pdf[i][71]
							+pdf[i][80]+pdf[i][81]+pdf[i][82]+pdf[i][83]+pdf[i][84]+pdf[i][85]+pdf[i][86]+pdf[i][87]
							+pdf[i][96]+pdf[i][97]+pdf[i][98]+pdf[i][99]+pdf[i][100]+pdf[i][101]+pdf[i][102]+pdf[i][103]
							+pdf[i][112]+pdf[i][113]+pdf[i][114]+pdf[i][115]+pdf[i][116]+pdf[i][117]+pdf[i][118]+pdf[i][119]
							+pdf[i][128]+pdf[i][129]+pdf[i][130]+pdf[i][131]+pdf[i][132]+pdf[i][133]+pdf[i][134]+pdf[i][135]
							+pdf[i][144]+pdf[i][145]+pdf[i][146]+pdf[i][147]+pdf[i][148]+pdf[i][149]+pdf[i][150]+pdf[i][151]
							+pdf[i][160]+pdf[i][161]+pdf[i][162]+pdf[i][163]+pdf[i][164]+pdf[i][165]+pdf[i][166]+pdf[i][167]
							+pdf[i][176]+pdf[i][177]+pdf[i][178]+pdf[i][179]+pdf[i][180]+pdf[i][181]+pdf[i][182]+pdf[i][183]
							+pdf[i][192]+pdf[i][193]+pdf[i][194]+pdf[i][195]+pdf[i][196]+pdf[i][197]+pdf[i][198]+pdf[i][199]
							+pdf[i][208]+pdf[i][209]+pdf[i][210]+pdf[i][211]+pdf[i][212]+pdf[i][213]+pdf[i][214]+pdf[i][215]
							+pdf[i][224]+pdf[i][225]+pdf[i][226]+pdf[i][227]+pdf[i][228]+pdf[i][229]+pdf[i][230]+pdf[i][231]
							+pdf[i][240]+pdf[i][241]+pdf[i][242]+pdf[i][243]+pdf[i][244]+pdf[i][245]+pdf[i][246]+pdf[i][247]) ) ;

		LLR[8*i+5] = log ( (pdf[i][4]+pdf[i][5]+pdf[i][6]+pdf[i][7]+pdf[i][12]+pdf[i][13]+pdf[i][14]+pdf[i][15]
							+pdf[i][20]+pdf[i][21]+pdf[i][22]+pdf[i][23]+pdf[i][28]+pdf[i][29]+pdf[i][30]+pdf[i][31]
							+pdf[i][36]+pdf[i][37]+pdf[i][38]+pdf[i][39]+pdf[i][44]+pdf[i][45]+pdf[i][46]+pdf[i][47]
							+pdf[i][52]+pdf[i][53]+pdf[i][54]+pdf[i][55]+pdf[i][60]+pdf[i][61]+pdf[i][62]+pdf[i][63]
							+pdf[i][68]+pdf[i][69]+pdf[i][70]+pdf[i][71]+pdf[i][76]+pdf[i][77]+pdf[i][78]+pdf[i][79]
							+pdf[i][84]+pdf[i][85]+pdf[i][86]+pdf[i][87]+pdf[i][92]+pdf[i][93]+pdf[i][94]+pdf[i][95]
							+pdf[i][100]+pdf[i][101]+pdf[i][102]+pdf[i][103]+pdf[i][108]+pdf[i][109]+pdf[i][110]+pdf[i][111]
							+pdf[i][116]+pdf[i][117]+pdf[i][118]+pdf[i][119]+pdf[i][124]+pdf[i][125]+pdf[i][126]+pdf[i][127]
							+pdf[i][132]+pdf[i][133]+pdf[i][134]+pdf[i][135]+pdf[i][140]+pdf[i][141]+pdf[i][142]+pdf[i][143]
							+pdf[i][148]+pdf[i][149]+pdf[i][150]+pdf[i][151]+pdf[i][156]+pdf[i][157]+pdf[i][158]+pdf[i][159]
							+pdf[i][164]+pdf[i][165]+pdf[i][166]+pdf[i][167]+pdf[i][172]+pdf[i][173]+pdf[i][174]+pdf[i][175]
							+pdf[i][180]+pdf[i][181]+pdf[i][182]+pdf[i][183]+pdf[i][188]+pdf[i][189]+pdf[i][190]+pdf[i][191]
							+pdf[i][196]+pdf[i][197]+pdf[i][198]+pdf[i][199]+pdf[i][204]+pdf[i][205]+pdf[i][206]+pdf[i][207]
							+pdf[i][212]+pdf[i][213]+pdf[i][214]+pdf[i][215]+pdf[i][220]+pdf[i][221]+pdf[i][222]+pdf[i][223]
							+pdf[i][228]+pdf[i][229]+pdf[i][230]+pdf[i][231]+pdf[i][236]+pdf[i][237]+pdf[i][238]+pdf[i][239]
							+pdf[i][244]+pdf[i][245]+pdf[i][246]+pdf[i][247]+pdf[i][252]+pdf[i][253]+pdf[i][254]+pdf[i][255]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][2]+pdf[i][3]+pdf[i][8]+pdf[i][9]+pdf[i][10]+pdf[i][11]
							+pdf[i][16]+pdf[i][17]+pdf[i][18]+pdf[i][19]+pdf[i][24]+pdf[i][25]+pdf[i][26]+pdf[i][27]
							+pdf[i][32]+pdf[i][33]+pdf[i][34]+pdf[i][35]+pdf[i][40]+pdf[i][41]+pdf[i][42]+pdf[i][43]
							+pdf[i][48]+pdf[i][49]+pdf[i][50]+pdf[i][51]+pdf[i][56]+pdf[i][57]+pdf[i][58]+pdf[i][59]
							+pdf[i][64]+pdf[i][65]+pdf[i][66]+pdf[i][67]+pdf[i][72]+pdf[i][73]+pdf[i][74]+pdf[i][75]
							+pdf[i][80]+pdf[i][81]+pdf[i][82]+pdf[i][83]+pdf[i][88]+pdf[i][89]+pdf[i][90]+pdf[i][91]
							+pdf[i][96]+pdf[i][97]+pdf[i][98]+pdf[i][99]+pdf[i][104]+pdf[i][105]+pdf[i][106]+pdf[i][107]
							+pdf[i][112]+pdf[i][113]+pdf[i][114]+pdf[i][115]+pdf[i][120]+pdf[i][121]+pdf[i][122]+pdf[i][123]
							+pdf[i][128]+pdf[i][129]+pdf[i][130]+pdf[i][131]+pdf[i][136]+pdf[i][137]+pdf[i][138]+pdf[i][139]
							+pdf[i][144]+pdf[i][145]+pdf[i][146]+pdf[i][147]+pdf[i][152]+pdf[i][153]+pdf[i][154]+pdf[i][155]
							+pdf[i][160]+pdf[i][161]+pdf[i][162]+pdf[i][163]+pdf[i][168]+pdf[i][169]+pdf[i][170]+pdf[i][171]
							+pdf[i][176]+pdf[i][177]+pdf[i][178]+pdf[i][179]+pdf[i][184]+pdf[i][185]+pdf[i][186]+pdf[i][187]
							+pdf[i][192]+pdf[i][193]+pdf[i][194]+pdf[i][195]+pdf[i][200]+pdf[i][201]+pdf[i][202]+pdf[i][203]
							+pdf[i][208]+pdf[i][209]+pdf[i][210]+pdf[i][211]+pdf[i][216]+pdf[i][217]+pdf[i][218]+pdf[i][219]
							+pdf[i][224]+pdf[i][225]+pdf[i][226]+pdf[i][227]+pdf[i][232]+pdf[i][233]+pdf[i][234]+pdf[i][235]
							+pdf[i][240]+pdf[i][241]+pdf[i][242]+pdf[i][243]+pdf[i][248]+pdf[i][249]+pdf[i][250]+pdf[i][251]) ) ;

		LLR[8*i+6] = log ( (pdf[i][2]+pdf[i][3]+pdf[i][6]+pdf[i][7]+pdf[i][10]+pdf[i][11]+pdf[i][14]+pdf[i][15]
							+pdf[i][18]+pdf[i][19]+pdf[i][22]+pdf[i][23]+pdf[i][26]+pdf[i][27]+pdf[i][30]+pdf[i][31]
							+pdf[i][34]+pdf[i][35]+pdf[i][38]+pdf[i][39]+pdf[i][42]+pdf[i][43]+pdf[i][46]+pdf[i][47]
							+pdf[i][50]+pdf[i][51]+pdf[i][54]+pdf[i][55]+pdf[i][58]+pdf[i][59]+pdf[i][62]+pdf[i][63]
							+pdf[i][66]+pdf[i][67]+pdf[i][70]+pdf[i][71]+pdf[i][74]+pdf[i][75]+pdf[i][78]+pdf[i][79]
							+pdf[i][82]+pdf[i][83]+pdf[i][86]+pdf[i][87]+pdf[i][90]+pdf[i][91]+pdf[i][94]+pdf[i][95]
							+pdf[i][98]+pdf[i][99]+pdf[i][102]+pdf[i][103]+pdf[i][106]+pdf[i][107]+pdf[i][110]+pdf[i][111]
							+pdf[i][114]+pdf[i][115]+pdf[i][118]+pdf[i][119]+pdf[i][122]+pdf[i][123]+pdf[i][126]+pdf[i][127]
							+pdf[i][130]+pdf[i][131]+pdf[i][134]+pdf[i][135]+pdf[i][138]+pdf[i][139]+pdf[i][142]+pdf[i][143]
							+pdf[i][146]+pdf[i][147]+pdf[i][150]+pdf[i][151]+pdf[i][154]+pdf[i][155]+pdf[i][158]+pdf[i][159]
							+pdf[i][162]+pdf[i][163]+pdf[i][166]+pdf[i][167]+pdf[i][170]+pdf[i][171]+pdf[i][174]+pdf[i][175]
							+pdf[i][178]+pdf[i][179]+pdf[i][182]+pdf[i][183]+pdf[i][186]+pdf[i][187]+pdf[i][190]+pdf[i][191]
							+pdf[i][194]+pdf[i][195]+pdf[i][198]+pdf[i][199]+pdf[i][202]+pdf[i][203]+pdf[i][206]+pdf[i][207]
							+pdf[i][210]+pdf[i][211]+pdf[i][214]+pdf[i][215]+pdf[i][218]+pdf[i][219]+pdf[i][222]+pdf[i][223]
							+pdf[i][226]+pdf[i][227]+pdf[i][230]+pdf[i][231]+pdf[i][234]+pdf[i][235]+pdf[i][238]+pdf[i][239]
							+pdf[i][242]+pdf[i][243]+pdf[i][246]+pdf[i][247]+pdf[i][250]+pdf[i][251]+pdf[i][254]+pdf[i][255]) 
						  / (pdf[i][0]+pdf[i][1]+pdf[i][4]+pdf[i][5]+pdf[i][8]+pdf[i][9]+pdf[i][12]+pdf[i][13]
							+pdf[i][16]+pdf[i][17]+pdf[i][20]+pdf[i][21]+pdf[i][24]+pdf[i][25]+pdf[i][28]+pdf[i][29]
							+pdf[i][32]+pdf[i][33]+pdf[i][36]+pdf[i][37]+pdf[i][40]+pdf[i][41]+pdf[i][44]+pdf[i][45]
							+pdf[i][48]+pdf[i][49]+pdf[i][52]+pdf[i][53]+pdf[i][56]+pdf[i][57]+pdf[i][60]+pdf[i][61]
							+pdf[i][64]+pdf[i][65]+pdf[i][68]+pdf[i][69]+pdf[i][72]+pdf[i][73]+pdf[i][76]+pdf[i][77]
							+pdf[i][80]+pdf[i][81]+pdf[i][84]+pdf[i][85]+pdf[i][88]+pdf[i][89]+pdf[i][92]+pdf[i][93]
							+pdf[i][96]+pdf[i][97]+pdf[i][100]+pdf[i][101]+pdf[i][104]+pdf[i][105]+pdf[i][108]+pdf[i][109]
							+pdf[i][112]+pdf[i][113]+pdf[i][116]+pdf[i][117]+pdf[i][120]+pdf[i][121]+pdf[i][124]+pdf[i][125]
							+pdf[i][128]+pdf[i][129]+pdf[i][132]+pdf[i][133]+pdf[i][136]+pdf[i][137]+pdf[i][140]+pdf[i][141]
							+pdf[i][144]+pdf[i][145]+pdf[i][148]+pdf[i][149]+pdf[i][152]+pdf[i][153]+pdf[i][156]+pdf[i][157]
							+pdf[i][160]+pdf[i][161]+pdf[i][164]+pdf[i][165]+pdf[i][168]+pdf[i][169]+pdf[i][172]+pdf[i][173]
							+pdf[i][176]+pdf[i][177]+pdf[i][180]+pdf[i][181]+pdf[i][184]+pdf[i][185]+pdf[i][188]+pdf[i][189]
							+pdf[i][192]+pdf[i][193]+pdf[i][196]+pdf[i][197]+pdf[i][200]+pdf[i][201]+pdf[i][204]+pdf[i][205]
							+pdf[i][208]+pdf[i][209]+pdf[i][212]+pdf[i][213]+pdf[i][216]+pdf[i][217]+pdf[i][220]+pdf[i][221]
							+pdf[i][224]+pdf[i][225]+pdf[i][228]+pdf[i][229]+pdf[i][232]+pdf[i][233]+pdf[i][236]+pdf[i][237]
							+pdf[i][240]+pdf[i][241]+pdf[i][244]+pdf[i][245]+pdf[i][248]+pdf[i][249]+pdf[i][252]+pdf[i][253]) ) ;

		LLR[8*i+7] = log ( (pdf[i][1]+pdf[i][3]+pdf[i][5]+pdf[i][7]+pdf[i][9]+pdf[i][11]+pdf[i][13]+pdf[i][15]
							+pdf[i][17]+pdf[i][19]+pdf[i][21]+pdf[i][23]+pdf[i][25]+pdf[i][27]+pdf[i][29]+pdf[i][31]
							+pdf[i][33]+pdf[i][35]+pdf[i][37]+pdf[i][39]+pdf[i][41]+pdf[i][43]+pdf[i][45]+pdf[i][47]
							+pdf[i][49]+pdf[i][51]+pdf[i][53]+pdf[i][55]+pdf[i][57]+pdf[i][59]+pdf[i][61]+pdf[i][63]
							+pdf[i][65]+pdf[i][67]+pdf[i][69]+pdf[i][71]+pdf[i][73]+pdf[i][75]+pdf[i][77]+pdf[i][79]
							+pdf[i][81]+pdf[i][83]+pdf[i][85]+pdf[i][87]+pdf[i][89]+pdf[i][91]+pdf[i][93]+pdf[i][95]
							+pdf[i][97]+pdf[i][99]+pdf[i][101]+pdf[i][103]+pdf[i][105]+pdf[i][107]+pdf[i][109]+pdf[i][111]
							+pdf[i][113]+pdf[i][115]+pdf[i][117]+pdf[i][119]+pdf[i][121]+pdf[i][123]+pdf[i][125]+pdf[i][127]
							+pdf[i][129]+pdf[i][131]+pdf[i][133]+pdf[i][135]+pdf[i][137]+pdf[i][139]+pdf[i][141]+pdf[i][143]
							+pdf[i][145]+pdf[i][147]+pdf[i][149]+pdf[i][151]+pdf[i][153]+pdf[i][155]+pdf[i][157]+pdf[i][159]
							+pdf[i][161]+pdf[i][163]+pdf[i][165]+pdf[i][167]+pdf[i][169]+pdf[i][171]+pdf[i][173]+pdf[i][175]
							+pdf[i][177]+pdf[i][179]+pdf[i][181]+pdf[i][183]+pdf[i][185]+pdf[i][187]+pdf[i][189]+pdf[i][191]
							+pdf[i][193]+pdf[i][195]+pdf[i][197]+pdf[i][199]+pdf[i][201]+pdf[i][203]+pdf[i][205]+pdf[i][207]
							+pdf[i][209]+pdf[i][211]+pdf[i][213]+pdf[i][215]+pdf[i][217]+pdf[i][219]+pdf[i][221]+pdf[i][223]
							+pdf[i][225]+pdf[i][227]+pdf[i][229]+pdf[i][231]+pdf[i][233]+pdf[i][235]+pdf[i][237]+pdf[i][239]
							+pdf[i][241]+pdf[i][243]+pdf[i][245]+pdf[i][247]+pdf[i][249]+pdf[i][251]+pdf[i][253]+pdf[i][255]) 
						  / (pdf[i][0]+pdf[i][2]+pdf[i][4]+pdf[i][6]+pdf[i][8]+pdf[i][10]+pdf[i][12]+pdf[i][14]
							+pdf[i][16]+pdf[i][18]+pdf[i][20]+pdf[i][22]+pdf[i][24]+pdf[i][26]+pdf[i][28]+pdf[i][30]
							+pdf[i][32]+pdf[i][34]+pdf[i][36]+pdf[i][38]+pdf[i][40]+pdf[i][42]+pdf[i][44]+pdf[i][46]
							+pdf[i][48]+pdf[i][50]+pdf[i][52]+pdf[i][54]+pdf[i][56]+pdf[i][58]+pdf[i][60]+pdf[i][62]
							+pdf[i][64]+pdf[i][66]+pdf[i][68]+pdf[i][70]+pdf[i][72]+pdf[i][74]+pdf[i][76]+pdf[i][78]
							+pdf[i][80]+pdf[i][82]+pdf[i][84]+pdf[i][86]+pdf[i][88]+pdf[i][90]+pdf[i][92]+pdf[i][94]
							+pdf[i][96]+pdf[i][98]+pdf[i][100]+pdf[i][102]+pdf[i][104]+pdf[i][106]+pdf[i][108]+pdf[i][110]
							+pdf[i][112]+pdf[i][114]+pdf[i][116]+pdf[i][118]+pdf[i][120]+pdf[i][122]+pdf[i][124]+pdf[i][126]
							+pdf[i][128]+pdf[i][130]+pdf[i][132]+pdf[i][134]+pdf[i][136]+pdf[i][138]+pdf[i][140]+pdf[i][142]
							+pdf[i][144]+pdf[i][146]+pdf[i][148]+pdf[i][150]+pdf[i][152]+pdf[i][154]+pdf[i][156]+pdf[i][158]
							+pdf[i][160]+pdf[i][162]+pdf[i][164]+pdf[i][166]+pdf[i][168]+pdf[i][170]+pdf[i][172]+pdf[i][174]
							+pdf[i][176]+pdf[i][178]+pdf[i][180]+pdf[i][182]+pdf[i][184]+pdf[i][186]+pdf[i][188]+pdf[i][190]
							+pdf[i][192]+pdf[i][194]+pdf[i][196]+pdf[i][198]+pdf[i][200]+pdf[i][202]+pdf[i][204]+pdf[i][206]
							+pdf[i][208]+pdf[i][210]+pdf[i][212]+pdf[i][214]+pdf[i][216]+pdf[i][218]+pdf[i][220]+pdf[i][222]
							+pdf[i][224]+pdf[i][226]+pdf[i][228]+pdf[i][230]+pdf[i][232]+pdf[i][234]+pdf[i][236]+pdf[i][238]
							+pdf[i][240]+pdf[i][242]+pdf[i][244]+pdf[i][246]+pdf[i][248]+pdf[i][250]+pdf[i][252]+pdf[i][254]) ) ;
	}
}



