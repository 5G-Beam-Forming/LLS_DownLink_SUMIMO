/*--------------------------------------------------------------------------------------*/
// method_ChannelDeoding.cpp
// - Class part for NR CODEC
// - [REF] 3GPP TS 38.212
// - Includes methods for Receiver
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/

#include "stdafx.h"
#include "method_ChannelCoding.h"



inline char Sign(double input)
{
    if( input >= 0 ) return(1);
    else return(-1);
}

inline double FLUT(double input){ 	
	double tmp;
	tmp = (double)fabs(input);
	if( tmp < (double)0.00001 ) 
		tmp = (double)0.00001;
	tmp = (double)log( fabs(tanh(tmp/2.)) );
	return(tmp);
}

inline char HD(double input){
    if( input <= 0. ) 
		return(1);
    else 
		return(0);
}

inline int xor2( int a, int b ) { 
	return( (a || b) && !(a && b) ); 
}

/*-------------- Public methods --------------*/
int* NR_ChannelCoding::Rx_GetDecodedData ()
{
	return rx_a ;
}


void NR_ChannelCoding::Rx_CodeBlockConcatenation (double *LLR)
{
	int cnt = 0 ;
	int k, r, j ;
	k = 0 ;

	for (k=0; k<G; k++)
		rx_g[k] = LLR[k] ;

	k = 0 ;
	for (r=0; r<C; r++) {
		for(j=0; j<E[r]; j++){
			rx_f[r][j] = rx_g[k];
			k++;
		}
	}

}


void NR_ChannelCoding::Rx_RateMatching (NR_CODEC *nr_codec)	// additional coding  (NR_CODEC�� input���� ����)
{

	int cnt = 0;
	int r, j, i, k, n ;	

	for (r=0; r<C; r++){
		for(j=0; j<E[r]/Q; j++){
			for(i=0; i<Q; i++){
				rx_e[r][j+i*E[r]/Q] = rx_f[r][i+j*Q] ;
			}
		}
	}

	for(r=0; r<C; r++){
		for (j = 0;j < N_H; j++)
			rx_d[r][j] = 0;
	}

	rv_id = nr_codec->harq.HARQ_index[nr_codec->harq.HARQ_iter] ;	// additional coding


	/* ���� �� code
	//for (r=0; r<C; r++){
	//	for(j=0; j<E[r]; j++){
	//		if(j<K_p-2*Z_c){
	//			rx_d[r][j+2*Z_c] = rx_e[r][j] ;
	//		}
	//		else if (j<N_cb-null_cnt){
	//			rx_d[r][j+2*Z_c+null_cnt] = rx_e[r][j] ;
	//		}
	//		// rx_e���� rx_d�� �ѹ��� �о� �ö��� �Ʒ� else if�� else�� �ּ�ó��
	//		else if ( (j>=N_cb-null_cnt) && (j%(N_cb-null_cnt)<K_p-2*Z_c) ){
	//			rx_d[r][j%(N_cb-null_cnt)+2*Z_c] = rx_d[r][j%(N_cb-null_cnt)+2*Z_c] + rx_e[r][j] ;
	//		}
	//		else{
	//			rx_d[r][j%(N_cb-null_cnt)+2*Z_c+null_cnt] = rx_d[r][j%(N_cb-null_cnt)+2*Z_c+null_cnt] + rx_e[r][j] ;
	//		}
	//	}
	//}
	*/
	

	/* ���� �� code */		// additional coding
	// �� code�� //from ~ //to ���� ����
	//from
	for (r=0; r<C; r++){		
		j = 0;					
		k = 0;					

		while(1){				
			if( ((k0[rv_id]+j)%N_cb+2*Z_c < null_position[0]) || ((k0[rv_id]+j)%N_cb+2*Z_c > null_position[K-K_p-1]) ){
				if (k==E[r]) break;
				else if (E[r]<N_cb){// shortening
					rx_d[r][(k0[rv_id]+j)%N_cb+2*Z_c] = rx_e[r][k] ;
					k++;
				}
				else if (E[r]>N_cb){//repetition
					if (j<N_cb)
						rx_d[r][(k0[rv_id]+j)%N_cb+2*Z_c] = rx_e[r][k] ;
					else
						rx_d[r][(k0[rv_id]+j)%N_cb+2*Z_c] = rx_d[r][(k0[rv_id]+j)%N_cb+2*Z_c] + rx_e[r][k] ;
					k++;
				}
			}
			j++ ;
		}

		// HARQ process
		if (nr_codec->harq.HARQ_indicator == 1){
			if (nr_codec->harq.HARQ_iter == 0){
				for (r=0; r<C; r++){
					for (n=0;n<N_H;n++)
						rx_d_memory[r][n] = rx_d[r][n] ;					
				}
			}
			else{
				for (r=0; r<C; r++){
					for (n=0;n<N_H;n++){
						rx_d[r][n] = (rx_d[r][n] + rx_d_memory[r][n]) ;	
						rx_d_memory[r][n] = rx_d[r][n] ;
					}
				}
			}
		}	
	}
	// to

}


void NR_ChannelCoding::Rx_LDPCDecoding (NR_CODEC *nr_codec) 
{

	int k, m, r ;
	
	//int decoding_type ;

	int real_col, real_row ;
	int array_col ;

	int iter_cnt ;

	double beta, tmp_beta ;
	int sgn, tmp_sgn ;

	double* sort_reg ;
	int QQ, tmp_QQ ;

	int *col_cnt ;
	double sum_Lr ;

	double factor = 0.75 ;

	int tmp_bit, check=0;

	int stop ;
	
	col_cnt = new int[M_H] ;
	
	
	// //////////////
	// Initialization
	// //////////////
	//nr_codec->LLR_limit;
	for (r=0; r<C; r++){
		for (k=0; k<2*Z_c; k++){
			rx_d[r][k] = 0 ;
		}
		for(k=0;k<null_cnt; k++){
			rx_d[r][null_position[k]] = -LLR_threshold ;
		}		
		for(k=0; k<N_H; k++){
			if(rx_d[r][k] != 0){
				rx_d[r][k] = -rx_d[r][k] ; // demodulator�� LLR �� log p(1)/p(0)�� ������
			}
		}
	}

	for(r=0; r<C; r++){
		for(k=0; k<N_H; k++){
			HD_Y[r][k] = HD(rx_d[r][k]) ;
		}
	}


	for (r=0; r<C; r++){
		for (k=0; k<M_H; k++){ //r
			for (m=0; m<row_weight[k]; m++){//c
				real_col = col_index[m+ max_row_weight*k] ;

				if(real_col > -1){
					Lq[r][m+ max_row_weight*k] = rx_d[r][real_col] ;
					HD_Lq[r][m+ max_row_weight*k] = HD_Y[r][real_col] ;
				}
			}
		}
	}


	for(r=0; r<C; r++){
		// //////////////////
		// Decoding iteration
		// //////////////////
		stop = 0;

		for (iter_cnt=0; iter_cnt<MAX_ITER; iter_cnt++){
			// Check node updates
			if (decoding_type == 0){ // Sum-product
				for (k=0; k<M_H; k++){
					beta = 0;
					sgn = 1;
					for(m=0; m<row_weight[k]; m++){
						sgn = sgn * Sign(Lq[r][m + max_row_weight*k]) ;
						beta = beta + FLUT(Lq[r][m + max_row_weight*k]) ; 
					}
					for(m=0; m<row_weight[k]; m++){
						tmp_sgn = sgn * Sign(Lq[r][m + max_row_weight*k]) ;
						tmp_beta = beta - FLUT(Lq[r][m + max_row_weight*k]) ; 
						Lr[r][m + max_row_weight*k] = -(double)tmp_sgn * FLUT(tmp_beta) ;
					}
				}			
			}
			else if (decoding_type == 1) {//Min-sum
				int i, j;
				double tmp ;
				for(m=0; m<M_H; m++){
					sort_reg = new double[row_weight[m]] ;
					QQ = 0 ;

					for(k=0; k<row_weight[m]; k++){
						QQ = xor2(QQ, HD_Lq[r][max_row_weight*m+k]) ;
						sort_reg[k] = fabs(Lq[r][max_row_weight*m+k]) ;
					}

					// sorting to find first and second minimum of Lq			
					for(j=0; j<row_weight[m]; j++){
						tmp = sort_reg[j] ;
						i = j ;
						while (i>0 && sort_reg[i-1]>tmp){
							sort_reg[i] = sort_reg[i-1] ;
							i-- ;
						}
						sort_reg[i] = tmp ;
					}
					for(k=0; k<row_weight[m]; k++){
						tmp_QQ    = xor2(QQ, HD_Lq[r][max_row_weight*m + k]) ;
						real_col  = col_index[max_row_weight*m + k] ;
						tmp_QQ    = xor2(tmp_QQ, HD_Y[r][real_col]) ;                   
                    
						if( fabs( Lq[r][max_row_weight*m + k]) == sort_reg[0] ){
							//Lr[r][max_row_weight*m + k] = (double)( xor2(tmp_QQ,1) - tmp_QQ) * sort_reg[1] ;
							Lr[r][max_row_weight*m + k] = (double)( xor2(tmp_QQ,1) - tmp_QQ) * sort_reg[1] *factor;
						}
						else{
							//Lr[r][max_row_weight*m + k] = (double)( xor2(tmp_QQ,1) - tmp_QQ) * sort_reg[0] ;
							Lr[r][max_row_weight*m + k] = (double)( xor2(tmp_QQ,1) - tmp_QQ) * sort_reg[0] *factor;
						}
					}
					delete sort_reg ;
				}
			}
		
			// Pseudopoteriori probability 
			for (k=0; k<M_H; k++)	
				col_cnt[k] = 0;


			for(m=0; m<N_H; m++)
			{			
				sum_Lr = 0;
				for (k=0; k<col_weight[m]; k++){
					real_row = row_index[k + max_col_weight*m] ; 
					array_col = col_cnt[real_row] ;
					sum_Lr = sum_Lr + Lr[r][array_col+max_row_weight*real_row] ;
					col_cnt[real_row] = array_col + 1 ;
				}

				if (decoding_type == 0){
					rx_c_LLR[r][m] = rx_d[r][m] + sum_Lr ;
					X_hat[r][m] = HD(rx_c_LLR[r][m]) ;
					LLR_out[r][m] = rx_c_LLR[r][m] ;
				}

				else if (decoding_type == 1){
					rx_c_LLR[r][m] = fabs(rx_d[r][m]) + sum_Lr ;
					X_hat[r][m] = xor2(HD(rx_c_LLR[r][m]), HD_Y[r][m]) ;
					LLR_out[r][m] = rx_c_LLR[r][m] * Sign(rx_d[r][m]) ;
				}
			}
			check = 0;
			// for parity bits
			for (m=0; m<M_H; m++){
				tmp_bit = 0;
				for(k=0; k<row_weight[m]; k++){
					real_col = col_index[k + max_row_weight*m] ;
					tmp_bit = xor2(tmp_bit, X_hat[r][real_col]) ;
				}
				check = check + tmp_bit ;
			}
			if( check == 0 )
			{
				stop = iter_cnt;
				iter_cnt  = MAX_ITER; // for early termination
			}

			// Variable node updates
			for (m=0; m<M_H; m++)
			{
				for (k=0; k<row_weight[m]; k++){
					real_col = col_index[k + max_row_weight*m] ;
					Lq[r][k + max_row_weight*m] = rx_c_LLR[r][real_col] - Lr[r][k + max_row_weight*m] ;

					if(Lq[r][k + max_row_weight*m] > LLR_threshold)
						Lq[r][k + max_row_weight*m] = LLR_threshold;
					else if (Lq[r][k + max_row_weight*m] < -LLR_threshold)
						Lq[r][k + max_row_weight*m] = -LLR_threshold;
					else{}


					if( Lq[r][max_row_weight*m + k] >= 0. )
						HD_Lq[r][max_row_weight*m + k] = HD_Y[r][real_col] ;
					else
						HD_Lq[r][max_row_weight*m + k] = xor2(HD_Y[r][real_col], 1);
				}
			}
		}//for (iter_cnt=0; iter_cnt<MAX_ITER; iter_cnt++){
	}//for(r=0; r<C; r++)
	delete []col_cnt;
}



void NR_ChannelCoding::Rx_Segmentation ()
{
	int r ;
	int k ;
	int s ;

	for (r=0; r<C; r++){
		for (k=0; k<K; k++){
			rx_c[r][k] = X_hat[r][k] ;
		}
	}

	s = 0;	
	for (r=0; r<C; r++){
		for (k=0; k<K_p-L; k++){
			rx_b[s] = (int) rx_c[r][k] ;
			s++;
		}
	}

}

// ���� Rx_CRC_Detachment code�� �Ʒ��� code�� ��ü // additional coding // additional coding // additional coding
void NR_ChannelCoding::Rx_CRC_Detachment (NR_CODEC *nr_codec)	
{
	
	int *parity ;

	int ack_cnt = 0; 

	// CRC detachment
	for (int i=0; i<A; i++)
		rx_a[i] = rx_b[i] ;		

	if (A>3824)
		parity = Func_common_CRC_Calculation (rx_a, A, crctable24a,  CRC24_Length) ;	// CRC calculation

	else
		parity = Func_common_CRC_Calculation (rx_a, A, crctable16,  CRC16_Length) ;	// CRC calculation

	// CRC chcking
	
	if (A>3824){
		for (int i=0; i<CRC24_Length; i++) {
			if (parity[i] != rx_b[A+i]) {
				nr_codec->NACK_cnt++ ;	

				nr_codec->NACK = 1 ;
				nr_codec->ACK = 0 ;

				nr_codec->harq.HARQ_iter++;
				break ;
			}
			else
				ack_cnt++ ;
			if(ack_cnt == CRC24_Length){
				ack_cnt = 0;
				nr_codec->ACK_cnt++ ;

				nr_codec->ACK = 1 ;
				nr_codec->NACK = 0 ;

				nr_codec->harq.HARQ_iter = 0 ;
			}
		}
	}
	else{
		for (int i=0; i<CRC16_Length; i++) {
			if (parity[i] != rx_b[A+i]) {
				nr_codec->NACK_cnt++ ;	

				nr_codec->NACK = 1 ;
				nr_codec->ACK = 0 ;

				nr_codec->harq.HARQ_iter++;
				break ;
			}
			else
				ack_cnt++ ;
		}
		if(ack_cnt == CRC16_Length){
			ack_cnt = 0;
			nr_codec->ACK_cnt++ ;

			nr_codec->ACK = 1 ;
			nr_codec->NACK = 0 ;

			nr_codec->harq.HARQ_iter = 0 ;
		}
	}


	// Delete memory

	delete []parity ;

}



void NR_ChannelCoding::Func_rx_HardDecision (int *hard, double *soft, int length) 
{
	for (int i=0; i<length; i++)
		hard[i] = (soft[i] > 0) ? 1 : 0 ;
}
