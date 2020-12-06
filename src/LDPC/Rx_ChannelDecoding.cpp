/*--------------------------------------------------------------------------------------*/
// Rx_ChannelDecoding.cpp 
// - Encoder for NR LDPC codes
// - 3GPP TS 38.212
//                                                                              by KAIST
/*--------------------------------------------------------------------------------------*/


#include "stdafx.h"
#include "method_ChannelCoding.h"


int* Rx_ChannelDecoding (double *LLR, NR_ChannelCoding nr_channelcoding, NR_CODEC *nr_codec)
{
//	/*-------- Rx channel decoding procedure  -------*/
	nr_channelcoding.Rx_CodeBlockConcatenation (LLR) ;
	nr_channelcoding.Rx_RateMatching (nr_codec) ;		// additional coding	(NR_CODEC을 input으로 받음)
	nr_channelcoding.Rx_LDPCDecoding (nr_codec) ;
	nr_channelcoding.Rx_Segmentation () ;		// CRC24b
	nr_channelcoding.Rx_CRC_Detachment (nr_codec) ;	// CRC24a

	return nr_channelcoding.Rx_GetDecodedData () ;
}
