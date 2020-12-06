/*--------------------------------------------------------------------------------------*/
// AWGN_channel.cpp 
// - Add white Gaussian noise to the transmitted symbol
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "SISO_PDSCH.h"
#include "misc.h"


QAM Rx_AWGN_Channel (double std_dev, int SymbolLength, QAM *qam)
{

	for (int i=0; i<SymbolLength; i++) {
		qam->rx_real[i] = qam->tx_real[i] + std_dev*gasdis() ;
		qam->rx_imag[i] = qam->tx_imag[i] + std_dev*gasdis() ;

		
		//qam->rx_real[i] = qam->tx_real[i];
		//qam->rx_imag[i] = qam->tx_imag[i];
	}
	return *qam ;
}