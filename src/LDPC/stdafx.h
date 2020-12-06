/*--------------------------------------------------------------------------------------*/
// stdafx.h
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/

#ifndef _STDAFX_H_
#define _STDAFX_H_

/*---------- Header ----------*/
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
#include < windows.h >

using namespace std;
using std::ifstream;

/*---------- Define ----------*/
#define MaxNumFrameErr			100				// The maximum number of observing frame error samples
#define NumPrintFrame			2				// The number of observing frames for displaying

//#define NULLBitValue			3				// The value to set the filler bits (Before implementing 'Rate matcher', set the filler bits as 0)
#define CRC24_Length			24				// The length of CRC bits for CRC24A and CRC24B
#define CRC16_Length			16				// The length of CRC bits for CRC16

#define LLR_threshold			10				// Threshold of LLR magnitude


#define RR						16807.0         // multiplier
#define MM						2147483647.0    // modulus 
#define PI						3.141592		// $phi$



/*---------- Structures ----------*/

typedef struct LDPC {
	int LDPCMaxIter ;
	int LDPCDecOption ;
	double LDPC_MinSum_factor ;

	int LBRM_indicator ;// additional coding
	// 1:limited buffer rate matching, 0: othersize
	int TBS_LBRM ;// additional coding


} ldpc ;

//typedef struct RM {
//	int Indicator_LBRM ; // 1:limited buffer rate matching, 0: othersize
//
//} rm ;

typedef struct HARQ {		// additional coding
	int HARQ_indicator ;	// additional coding
	int HARQ_type ;			// additional coding
	int HARQ_iter ;			// additional coding
	int HARQ_max_iter ;		// additional coding
	int* HARQ_index ;		// additional coding
	
} harq ;// additional coding


typedef struct NR_CODEC {
	int data_length ;

	double code_rate ;
	int G ;
	int ACK ;
	int ACK_cnt ;			// additional coding
	int NACK ;
	int NACK_cnt ;			// additional coding

	LDPC ldpc ;
	HARQ harq ;				// additional coding

} nr_codec ;



// QAM modulation/demodulation
typedef struct QAM {
	double *tx_real ;
	double *tx_imag ;

	double *rx_real ;
	double *rx_imag ;
} qam ;



#endif
