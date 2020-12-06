/*--------------------------------------------------------------------------------------*/
// method_ChannelEnoding.cpp
// - Class part for NR CODEC
// - [REF] 3GPP TS 38.212
// - Includes methods for Transmitter
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "method_ChannelCoding.h"



/*-------------- Public methods --------------*/

int* NR_ChannelCoding::Tx_GetCodeword ()
{	
	return tx_g ;
}

// 아래의 Tx_data_memory code 추가		// additional coding
void NR_ChannelCoding::Tx_data_memory (int *data, NR_CODEC *nr_codec) 
{
	if (nr_codec->harq.HARQ_indicator == 1 && nr_codec->NACK == 0){
		for(int i=0;i<A;i++)
			data_memory[i] = data[i] ;
	}
	else if  (nr_codec->harq.HARQ_indicator == 1 && nr_codec->NACK == 1){
		for(int i=0;i<A;i++)
			data[i] = data_memory[i] ;
	}
}

void NR_ChannelCoding::Tx_CRC_Attachment (int *data) 
{
	int *parity ;
	if (A > 3824){
		parity = Func_common_CRC_Calculation (data, A, crctable24a, CRC24_Length) ;	// CRC calculation
	}
	else{
		parity = Func_common_CRC_Calculation (data, A, crctable16, CRC16_Length) ;	// CRC calculation
	}
	

	// CRC attachment
	for (int i=0; i<B; i++) {
		if (i < A)
			tx_b[i] = data[i] ;
		else
			tx_b[i] = parity[i-A] ;
	}
	delete []parity ;

}


void NR_ChannelCoding::Tx_Segmentation (NR_CODEC *nr_codec)
{
	int s, r, k ;
	int	*parity ;
	
	s = 0 ;
	for (r=0; r<C; r++){
		
		for(k=0; k<K_p-L; k++){
			tx_c[r][k] = tx_b[s];
			s++;
		}

		if (C > 1) {			
			parity = Func_common_CRC_Calculation (tx_c[r], K, crctable24b, CRC24_Length) ;	// CRC calculation
			for (k=K_p-L; k<K_p; k++)
				tx_c[r][k] = parity[k+L-K_p] ;

			delete []parity;	// moved by pjw
		}
		
		for(k=K_p; k<K; k++) 
			tx_c[r][k] = -1; //NULL bits

	}

}

void NR_ChannelCoding::Tx_LDPCEncoding ()
{
	int k, m, r;
	int *c_t ;
	int *w1, *w2, *w2_tmp;

	//null_cnt = 0;
	m = 0;

	c_t =		new int[K] ;
	w1 =		new int[4*Z_c]() ;
	w2 =		new int[N_H-K-4*Z_c]() ;
	w2_tmp =	new int[N_H-K-4*Z_c]() ;
	
	for (r=0; r<C; r++){
		
		for(k=0; k<4*Z_c; k++)
			w1[k] = 0;

		for(k=0; k<N_H-K-4*Z_c; k++){
			w2[k] = 0;
			w2_tmp[k] = 0;
		}

		for (k=2*Z_c; k<K; k++){
			if (tx_c[r][k] != -1){
				tx_d[r][k-2*Z_c] = tx_c[r][k] ;
				//tx_d[r][k-2*Z_c] = tx_c[0][k] ;
			}
			else{
				tx_c[r][k] = 0 ;	//for encoding of null bits
				tx_d[r][k-2*Z_c] = -1; // null bits to be removed in rate matching block
			}
		}


		for(k=0; k<K; k++){
			c_t[k] = tx_c[r][k] ;
		}


		for (k=0; k<4*Z_c; k++){// column
			for (m=0; m<max_row_weight_BB; m++){
				if (col_index_BB[m + max_row_weight_BB*k] != -1){
					w1[k] =  (w1[k] + tx_c[r][col_index_BB[m + max_row_weight_BB*k]])%2 ;
				}
				else
					break ;
			}
		}
		

		for (k=0; k<K + 4*Z_c; k++){
			if (k<K)
				w2_tmp[k] = c_t[k] ;
			else
				w2_tmp[k] = w1[k-K] ;
		}


		for (k=0; k<N_H-K-4*Z_c; k++){// column
			for (m=0; m<max_row_weight_DD; m++)
				if (col_index_DD[m + max_row_weight_DD*k] != -1)
					w2[k] =  (w2[k] + w2_tmp[col_index_DD[m + max_row_weight_DD*k]])%2 ;
		}

		
		for(k=K; k<N+2*Z_c; k++){
			if (k < K+4*Z_c)
				tx_d[r][k-2*Z_c] = w1[k-K] ;
			else
				tx_d[r][k-2*Z_c] = w2[k-K-4*Z_c] ;
		}

	}

	delete []c_t ;
	delete []w1 ;
	delete []w2 ;
	delete []w2_tmp ;

}



void NR_ChannelCoding::Tx_RateMatching (NR_CODEC *nr_codec) {

	int k, j, r, i;

	rv_id = nr_codec->harq.HARQ_index[nr_codec->harq.HARQ_iter] ;	// additional coding

	for (r=0; r<C; r++){
		k = 0 ;
		j = 0 ;

		while (k < E[r]){
			
			if(tx_d[r][ (k0[rv_id]+j)%N_cb ] != -1) {				// if (tx_d[r][ (k0+j)%N_cb ] != NULL)
				tx_e[r][k] = tx_d[r][ (k0[rv_id]+j)%N_cb ] ;
				k++ ;
			}	
			j++ ;
		}
	}

	for (r = 0; r<C; r++) {
		for (j = 0; j<E[r] / Q; j++) {
			for (i = 0; i<Q; i++) {
				tx_f[r][i + j*Q] = tx_e[r][i*E[r] / Q + j];
			}
		}
	}

}



void NR_ChannelCoding::Tx_CodeBlockConcatenation ()
{
	int k, r, j ;

	k = 0 ;
	r = 0 ;

	while ( r < C ) {
		j = 0 ;
		while ( j < E[r] ) {
			tx_g[k] = (int) tx_f[r][j] ;
			k++ ;
			j++ ;
		}
		r++ ;
	}

}




