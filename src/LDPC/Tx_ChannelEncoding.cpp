/*--------------------------------------------------------------------------------------*/
// Tx_ChannelEncoding.cpp 
// - Encoder for NR LDPC codes
// - 3GPP TS 38.212
//                                                                              by KAIST
/*--------------------------------------------------------------------------------------*/


#include "stdafx.h"
#include "method_ChannelCoding.h"


int* Tx_ChannelEncoding (int *data, NR_ChannelCoding nr_channelcoding, NR_CODEC *nr_codec)
{
	///*-------- Tx channel coding procedure  -------*/
	nr_channelcoding.Tx_data_memory(data, nr_codec) ;	// additional coding
	nr_channelcoding.Tx_CRC_Attachment (data) ;	// CRC24a
	nr_channelcoding.Tx_Segmentation (nr_codec) ;		// 
	nr_channelcoding.Tx_LDPCEncoding () ;		// 
	nr_channelcoding.Tx_RateMatching (nr_codec) ;
	nr_channelcoding.Tx_CodeBlockConcatenation () ;

	return nr_channelcoding.Tx_GetCodeword () ;
}

