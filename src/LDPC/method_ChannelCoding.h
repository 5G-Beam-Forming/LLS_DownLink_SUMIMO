/*--------------------------------------------------------------------------------------*/
// method_ChannelCoding.h
// - Interface for NR_ChannelCoding class
// - [REF] 3GPP TS 38.212
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/

#ifndef _METHOD_CHANNELCODING_H_
#define _METHOD_CHANNELCODING_H_

#include "stdafx.h"



class NR_ChannelCoding 
{
public:

	// NR_ChannelCoding(int Mod_order, NR_CODEC *nr_codec) ;
	NR_ChannelCoding();
	virtual ~NR_ChannelCoding() ;
	void Initialization(int Mod_order, NR_CODEC *nr_codec);

	/*----- Tx ------*/
	void Tx_data_memory (int *data, NR_CODEC *nr_codec ) ;	// additional coding
	void Tx_CRC_Attachment (int *data) ;	// CRC24a
	void Tx_Segmentation (NR_CODEC *nr_codec) ;				// CRC24b
	void Tx_LDPCEncoding ();
	void Tx_RateMatching (NR_CODEC *nr_codec);
	void Tx_CodeBlockConcatenation () ;
	int* Tx_GetCodeword () ;


	/*----- Rx ------*/
	void Rx_CodeBlockConcatenation (double *LLR) ;
	void Rx_RateMatching (NR_CODEC *nr_codec) ;				// additional coding	(NR_CODEC을 input으로 받음)
	void Rx_LDPCDecoding (NR_CODEC *nr_codec) ;				// 
	void Rx_Segmentation () ;
	void Rx_CRC_Detachment (NR_CODEC *nr_codec) ;

	int* Rx_GetDecodedData () ;


	/*----- Common ------*/
	int GetE () ;
	int GetN () ;
	int GetG () ;
	int GetK () ;
	int GetBG () ;

private:
	/*----- Required parameters ------*/
	int A ;				// The number of data bits
	double R ;	// Code rate


	/*----- General setting ------*/
	unsigned int CRC24A, CRC24B, CRC16 ;						// CRC polynomials defined in 3GPP TS 38.212 v.1.2.0 (2017-11)
	unsigned int crctable24a[256], crctable24b[256] ;	// CRC tables calculated from polynomials defined in 3GPP TS 38.212 v.1.2.0 (2017-11)
	unsigned int crctable16[256] ;										// CRC tables calculated from polynomials defined in 3GPP TS 38.212 v.1.2.0 (2017-11)


	int* Func_common_CRC_Calculation (int *BitData, int BitDataLength, unsigned int *crctable, int CRC_Length) ;
	void Func_common_GetCRCTable () ;				// Get CRC polynomials defined in 38.212
	
	short** Func_GetBaseGraphTable () ;
	void Func_GetBaseGraph() ;

	void Func_GetZ_table () ;
	void Func_GetSegmentPara () ;
	
	short* Func_GetRowIndex(short *row_vec, int BG_row_start,int BG_row_end, int Max_Row_Weight, int BG_col_start,int BG_col_end, int Max_Col_Weight ) ;
	short* Func_GetColIndex(short *col_vec, int BG_row_start,int BG_row_end, int Max_Row_Weight, int BG_col_start,int BG_col_end, int Max_Col_Weight ) ;

	void Func_GetWeight() ;

	
	short** Func_GetCyclePermMatrix (short** matrix, int mtx_size, int shift_idx) ;
	short** Func_GetInvMatrix (short** mtx) ;
	void Func_GetBBMatrix ();
	void Func_GetInvBB_AA () ;
	void Func_GetInvBB_AA_max_row_weight () ;
	void Func_GetIndexInvBB_AA()  ;
	
	/* ----- Data memory ------*/
	int *data_memory ;	// additional coding

		
	/*----- CRC attachment ------*/
	int B ;		// The number of bits after CRC attachment

	// For Tx
	int *tx_b ;

	// For Rx
	int *rx_a ;
	int *rx_b ;


	/*----- Code block segmentation ------*/
	
	int K_cb ;	// Maximum code block size
	int L ;		// The number of CRC bits
	int C ;		// The number of code blocks
	int B_p ;	// Size of (code blocks + CRC bits), B+C*L
	int K_p ;	
	int K_b ;
	

	int BG_index;		// LDPC base graph index, 1 - BG#1, 2 - BG#2
	int Z_c ;			// Lifting size
	int i_LS ;			// Lifting size set index
	int Z_table[8][8] ; // Table of sets of LDPC lifting size Z

	int K ;		// Code block size

	// For Tx
	int **tx_c ;

	// For Rx
	double **rx_c ;
	double **rx_c_LLR ;

	/*----- LDPC coding ------*/
	int N ;	// Codeword bit sequence legnth
	int M ; // Parity bit sequence length
	short** BaseGraph_info_mtx ; // Base graph - (1st col: row index), (2nd col: column index), (3rd V_i,j)
	short** BaseGraph ; // Base graph
	
	short ** H ;				// Parity check matrix
	short* col_index ;		// Concatenation of non-zero index of all rows
	short* row_index ;		// Concatenation of non-zero index of all columns

	//short* col_index_AA ;	// Concatenation of non-zero index of rows of AA
	short* row_index_AA ;	// Concatenation of non-zero index of columns of AA

	short* col_index_BB ;	// Concatenation of non-zero index of rows of BB
	//short* row_index_BB ;	// Concatenation of non-zero index of columns of BB

	short* col_index_DD ;	// Concatenation of non-zero index of rows of DD
	//short* row_index_DD ;	// Concatenation of non-zero index of columns of DD

	short* row_weight ;		// Concatenation of row weight of all rows
	short* col_weight ;		// Concatenation of column weight of all rows

	int max_row_weight ;	
	int max_col_weight ;

	int max_row_weight_AA ;	
	int max_col_weight_AA ;

	int max_row_weight_BB ;	
	int max_col_weight_BB ;

	int max_row_weight_DD ;	
	int max_col_weight_DD ;

	int N_a ;	// N / Z_c
	int M_a ;	// M / Z_c
	int K_a ;	// K / Z_c

	int N_H ;	// # of columns of H
	int M_H ;	// # of rows of H
	int N_h ;	// N_H / Z_c
	int M_h ;	// M_H / Z_c
	
	int N_H_AA ;	// # of columns of H
	//int M_H_AA ;	// # of rows of H
	int N_h_AA ;	// N_H / Z_c
	int M_h_AA ;	// M_H / Z_c

	//int N_H_BB ;	// # of columns of H
	int M_H_BB ;	// # of rows of H
	int N_h_BB ;	// N_H / Z_c
	int M_h_BB ;	// M_H / Z_c

	//int N_H_DD ;	// # of columns of H
	int M_H_DD ;	// # of rows of H
	int N_h_DD ;	// N_H / Z_c
	int M_h_DD ;	// M_H / Z_c

	short** cyc_perm_mtx ;
	//short** zero_mtx ;

	//short** AA;	// Element matrix of H
	short** BB;	// Element matrix of H
	//short** DD;	// Element matrix of H
	short** inv_BB;	// BB^-1
	short** inv_BB_AA;// BB^-1 * AA

	int* null_position ;	// Positions of filler bits of LDPC encoding input bit sequence
	int null_cnt;			// # of filler bits of LDPC encoding input bit sequence
	// For Tx
	int **tx_d ;

	// For Rx
	double **rx_d ;
	double **rx_d_memory ;	// additional coding

	void Func_rx_HardDecision (int *hard, double *soft, int length) ;

	// LDPC decoding
	//double LLR_limit ;
	double **Lr ; 
	double **Lq ; 
	double **LLR_out;
	int **X_hat ;
	int **HD_Lq ;
	int **HD_Y ;

	int MAX_ITER ;
	int decoding_type ;

	/*----- Rate matching ------*/
	short **circular_buffer ;
	int N_cb ;
	int I_LBRM ;
	int N_ref ; 
	int TBS_LBRM ;
	double R_LBRM ; // additional coding
	int n_PRB_LBRM ;
	int maximum_num_PRB ;
	int CBGTI ;
	int CBGTI_present ; 
	int* E ;
	int C_p ;
	int N_L ;
	int Q_m ;
	//int N_IR ;
	int rv_id ;
	int* k0 ;		// additional coding
	int Q ;

	// For Tx
	short **tx_e ;
	short **tx_f ;
	// For Rx
	double **rx_e ;
	double **rx_f ;


	/*----- Code block concatenation ------*/
	int G ;

	// For Tx
	
	int *tx_g ;
	// For Rx
	
	double *rx_g ;

};


#endif