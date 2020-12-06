/*--------------------------------------------------------------------------------------*/
// method_QAM_Modulation.h
// - Interface for QAM_Modulation class
// - [REF] 3GPP TS 38.211
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/

#ifndef _METHOD_QAM_MODULATION_H_
#define _METHOD_QAM_MODULATION_H_


class QAM_Modulation 
{
public:
	QAM_Modulation(int M_ary, int CodeLength, QAM *qam) ;
	virtual ~QAM_Modulation() ;

	/*----- Tx ------*/
	QAM Tx_QAM_Modulator (int *codeword, QAM *qam) ;

	/*----- Rx ------*/
	double *Rx_QAM_Demodulator (double noise_var, QAM *qam) ;

	/*----- Common ------*/
	int GetSymbolLength () ;

private:
	/*----- Modulator ------*/
	int Func_Bin2Dec (int *BitSequence, int length) ;


	/*----- Demodulator (LLR calculator) ------*/
	double **pdf ;
	double *LLR ;

	void Func_rx_LLR_Limitter () ;
	// Ver. 1
	void Func_rx_Demod_M2 () ; 		// BPSK
	void Func_rx_Demod_M4 () ; 		// QPSK
	void Func_rx_Demod_M16 () ; 	// 16QAM
	void Func_rx_Demod_M64 () ; 	// 64QAM
	void Func_rx_Demod_M256 () ;	// 256QAM
	// Ver. 2
	//void Func_rx_Demodulator () ;
	

	/*----- Common ------*/
	int ModOrder ;
	int log2M ;
	int SymbolLength ;
	double **ModMapper ;

	double **Func_QAM_set_gen (int M_ary) ;

};


#endif