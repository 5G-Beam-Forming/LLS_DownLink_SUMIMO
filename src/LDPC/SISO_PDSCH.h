/*--------------------------------------------------------------------------------------*/
// SISO_PDSCH.h
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/

#ifndef _SISO_PDSCH_H_
#define _SISO_PDSCH_H_

#include "method_ChannelCoding.h"
#include "method_QAM_Modulation.h"

/*---------- Tx ----------*/
int *TBS_determination(int n_PRB, int N_symb_sh, int log2_Mod_order, NR_CODEC *nr_codec );
int* Tx_ChannelEncoding (int *data, NR_ChannelCoding nr_channelcoding, NR_CODEC *nr_codec) ;

/*---------- Rx ----------*/
QAM Rx_AWGN_Channel (double std_dev, int SymbolLength, QAM *qam) ;
int* Rx_ChannelDecoding (double *LLR, NR_ChannelCoding nr_channelcoding, NR_CODEC *nr_codec) ;


#endif