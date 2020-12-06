/*--------------------------------------------------------------------------------------*/
// method_ChannelCoding.cpp
// - Class part for NR CODEC
// - [REF] 3GPP TS 38.212
//                                                                           by LIT KAIST
/*--------------------------------------------------------------------------------------*/


#include "stdafx.h"
#include "method_ChannelCoding.h"



/*-------------- Construction & Destruction --------------*/

//NR_ChannelCoding::NR_ChannelCoding(int Mod_order, NR_CODEC *nr_codec)
void NR_ChannelCoding::Initialization(int Mod_order, NR_CODEC *nr_codec)
{
	int k, r, m, j, i ;
	m=0;

	/*----- Required parameters ------*/
	A = nr_codec->data_length ;	
	G = nr_codec->G ;
	R = nr_codec->code_rate ;
	

	/*----- General setting ------*/
	Func_common_GetCRCTable () ;			// Get CRC polynomials defined in 38.212 Sec. 5.1.1


	/*----- CRC attachment ------*/
	if(A>3824)
		B = A + CRC24_Length ;
	else//(A<=3824)
		B = A + CRC16_Length ;


	/*----- LDPC base graph selection ------*/
	if (A<=292 || (A<=3824 && R<=0.67) || R<=0.25)	// LDPC base graph index, 1 - BG#1, 2 - BG#2
		BG_index = 2 ;	
	else
		BG_index = 1 ;
	//BG_index = 1 ;


	/*----- Code block segmentation ------*/
	Func_GetZ_table () ;
	
	Func_GetSegmentPara() ;

	
	// Position and number of null bits in tx_c
	null_position = new int[K-K_p];
	null_cnt = K - K_p ;
	for (k=0; k<K - K_p;k++){
		null_position[k] = K_p+k;
	}


	//----- LDPC coding ------*/
		
	// Base graph
	if (BG_index == 1){
		BaseGraph_info_mtx = new short*[316];	
		for ( i=0; i<316; i++)
			BaseGraph_info_mtx[i] = new short[3] ;
	}
	else if (BG_index == 2){
		BaseGraph_info_mtx = new short*[197];	
		for ( i=0; i<197; i++)
			BaseGraph_info_mtx[i] = new short[3] ;
	}
	

	Func_GetBaseGraphTable() ;

	Func_GetBaseGraph ();
	

	// weight of row and column
	row_weight = new short[M_H]() ;
	col_weight = new short[N_H]() ;
	
	// index of row and column
	col_index = new short[max_row_weight * M_H] ;
	row_index = new short[max_col_weight * N_H] ;
	for (k=0;k<max_row_weight * M_H;k++)		col_index[k] = -1;
	for (k=0;k<max_col_weight * N_H;k++)		row_index[k] = -1;

	row_index_AA = new short[max_col_weight_AA * N_H_AA] ;
	for (k=0;k<max_col_weight_AA * N_H_AA;k++)		row_index_AA[k] = -1;

	col_index_DD = new short[max_row_weight_DD * M_H_DD] ;
	for (k=0;k<max_row_weight_DD * M_H_DD;k++)		col_index_DD[k] = -1;
	

// Parity check matrix structure 
	cyc_perm_mtx = new short*[Z_c] ;
	//zero_mtx = new short*[Z_c] ;
	for (k=0; k<Z_c; k++){
		cyc_perm_mtx[k] = new short[Z_c] ;
		//zero_mtx[k] = new short[Z_c]() ;
	}
	// get row and col index vectors
	row_index = Func_GetRowIndex(row_index, 0, M_h, max_row_weight, 0, N_h, max_col_weight)  ;
	col_index = Func_GetColIndex(col_index, 0, M_h, max_row_weight, 0, N_h, max_col_weight)  ;

	Func_GetWeight(); // row_weight, col_weight


	row_index_AA = Func_GetRowIndex(row_index_AA, 0, M_h_AA, max_row_weight_AA, 0, N_h_AA, max_col_weight_AA)  ;
	
	col_index_DD = Func_GetColIndex(col_index_DD, M_h_AA, M_h_AA+M_h_DD, max_row_weight_DD, 0, N_h_DD, max_col_weight_DD)  ;
	

	// inv(BB)*AA matrix structure
	BB = new short*[4*Z_c] ;		for (k=0; k<4*Z_c; k++)			BB[k] = new short[4*Z_c]() ;
	inv_BB = new short*[4*Z_c] ;	for (k=0; k<4*Z_c; k++)			inv_BB[k] = new short[4*Z_c]() ;
	inv_BB_AA = new short*[4*Z_c];	for (k=0; k<4*Z_c; k++)			inv_BB_AA[k] = new short[K]() ;

	Func_GetBBMatrix () ; // BB

	for(k=0; k<4*Z_c; k++){	for (m=0; m<4*Z_c; m++)		inv_BB[k][m] = BB[k][m] ;	}

	inv_BB = Func_GetInvMatrix(inv_BB);// inv(BB)
		
	Func_GetInvBB_AA ()  ;

	Func_GetInvBB_AA_max_row_weight () ;
		
	col_index_BB = new short[max_row_weight_BB * M_H_BB] ;

	for (k=0;k<max_row_weight_BB * M_H_BB;k++)		col_index_BB[k] = -1;

	Func_GetIndexInvBB_AA() ;
	

	for (k=0; k<4*Z_c; k++){
		delete []BB[k] ;
		delete []inv_BB[k] ;
		delete []inv_BB_AA[k];
	}
	delete []BB;
	delete []inv_BB;
	delete []inv_BB_AA;
	
	
	// Rate matching

	circular_buffer = new short*[C] ;
	for (k=0; k<C; k++)
		circular_buffer[k] = new short[N] ;
	
	CBGTI = 1; // 0 or 1, TBD
	CBGTI_present = 0 ; // 0 or 1, TBD
	N_L = 1;// 1 or 2,... TBD
	//Q_m = (int) ( log((double) Mod_order) / log(2.0) ) ;
	Q_m = Mod_order;
	if (CBGTI_present == 0)
		C_p = C;
	//else TBD
	E = new int[C] ;

	I_LBRM = nr_codec->ldpc.LBRM_indicator ;// 1: Limited buffer rate matching, 0: Otherwise	// additional coding
	R_LBRM = 2.0/3.0 ;		// additional coding
	TBS_LBRM = 1000 ;//Transport block size TBD, about ~273*12*7 or 273*12*14		// removal coding

	maximum_num_PRB = (int) ceil ((double)TBS_LBRM / 12.0) ;// TDB	


	if(I_LBRM == 0)
		N_cb = N ;
	else{
		TBS_LBRM = nr_codec->ldpc.TBS_LBRM ;						// additional coding
		N_ref = (int) floor(double(TBS_LBRM) / double(C*R_LBRM)) ;
		N_cb = min(N, N_ref) ; 
	}

	//if (maximum_num_PRB <= 32)								// removal coding
	//	n_PRB_LBRM = 32 ;										// removal coding
	//else if (maximum_num_PRB >= 33 && maximum_num_PRB <= 66)	// removal coding
	//	n_PRB_LBRM = 66 ;										// removal coding
	//else if (maximum_num_PRB >= 67 && maximum_num_PRB <= 107)	// removal coding
	//	n_PRB_LBRM = 107 ;										// removal coding
	//else if (maximum_num_PRB >= 108 && maximum_num_PRB <= 135)// removal coding	
	//	n_PRB_LBRM = 135 ;										// removal coding		
	//else if (maximum_num_PRB >= 136 && maximum_num_PRB <= 162)// removal coding
	//	n_PRB_LBRM = 162 ;										// removal coding
	//else if (maximum_num_PRB >= 162 && maximum_num_PRB <= 217)// removal coding
	//	n_PRB_LBRM = 217 ;										// removal coding
	//else														// removal coding
	//	n_PRB_LBRM = 273 ;										// removal coding

	j = 0 ;	

	for (r=0; r<C; r++){
		if (CBGTI == 0)
			E[r] = 0;
		else{			
			if(j <= C_p - (G/(N_L*Q_m))%C_p - 1)
				E[r] = N_L * Q_m * (int)floor(double(G)/double(N_L*Q_m*C_p))	;
			else
				E[r] = N_L * Q_m * (int)ceil(double(G)/double(N_L*Q_m*C_p))	;
			j++ ;
		}
	}

	//rv_id = 0;// 0,1,2,3		// removal coding
	
	k0 = new int[4];			// additional coding

	//if (rv_id == 0){			// removal coding
	//	if (BG_index == 1)		// removal coding
	//		k0 = 0 ;			// removal coding
	//	else					// removal coding
	//		k0 = 0 ;			// removal coding
	//}							// removal coding
	//else if (rv_id == 1){														// removal coding
	//	if (BG_index == 1)														// removal coding
	//		k0 = (int) floor((double) (17*N_cb) / (double) (66*Z_c) ) * Z_c ;	// removal coding
	//	else																	// removal coding
	//		k0 = (int) floor((double) (13*N_cb) / (double) (50*Z_c) ) * Z_c ;	// removal coding
	//}																			// removal coding
	//else if (rv_id == 2){{													// removal coding
	//	if (BG_index == 1){														// removal coding
	//		k0 = (int) floor((double) (33*N_cb) / (double) (66*Z_c) ) * Z_c ;	// removal coding
	//	else{																	// removal coding
	//		k0 = (int) floor((double) (25*N_cb) / (double) (50*Z_c) ) * Z_c ;	// removal coding
	//}																			// removal coding
	//else if (rv_id == 3){{													// removal coding
	//	if (BG_index == 1){														// removal coding
	//		k0 = (int) floor((double) (56*N_cb) / (double) (66*Z_c) ) * Z_c ;	// removal coding
	//	else{																	// removal coding
	//		k0 = (int) floor((double) (43*N_cb) / (double) (50*Z_c) ) * Z_c ;	// removal coding
	//}																			// removal coding

	if (BG_index == 1){		// additional coding
		k0[0] = 0 ;			// additional coding
		k0[1] = (int) floor((double) (17*N_cb) / (double) (66*Z_c) ) * Z_c ;	// additional coding
		k0[2] = (int) floor((double) (33*N_cb) / (double) (66*Z_c) ) * Z_c ;	// additional coding
		k0[3] = (int) floor((double) (56*N_cb) / (double) (66*Z_c) ) * Z_c ;	// additional coding
	}						// additional coding
	else if (BG_index == 2){// additional coding
		k0[0] = 0 ;			// additional coding
		k0[1] = (int) floor((double) (13*N_cb) / (double) (50*Z_c) ) * Z_c ;	// additional coding
		k0[2] = (int) floor((double) (25*N_cb) / (double) (50*Z_c) ) * Z_c ;	// additional coding
		k0[3] = (int) floor((double) (43*N_cb) / (double) (50*Z_c) ) * Z_c ;	// additional coding
	}

	//R1-1713231 test
		//k0[0] = 0 ;															// removal coding
		//k0[1] = (int) floor((double) (1*N_cb) / (double) (4*Z_c) ) * Z_c ;	// removal coding
		//k0[2] = (int) floor((double) (2*N_cb) / (double) (4*Z_c) ) * Z_c ;	// removal coding
		//k0[3] = (int) floor((double) (3*N_cb) / (double) (4*Z_c) ) * Z_c ;	// removal coding



	//// HARQ
	nr_codec->harq.HARQ_indicator = 1;//1: HARQ on, 0: HAQR off	// additional coding	
	nr_codec->harq.HARQ_type = 2 ;	 //1:CC-HARQ, 2:IR-HARQ		// additional coding	

	if (nr_codec->harq.HARQ_indicator == 1)			// additional coding
		nr_codec->harq.HARQ_max_iter = 1; // 1~16 ;	// additional coding
	else if (nr_codec->harq.HARQ_indicator == 0)	// additional coding
		nr_codec->harq.HARQ_max_iter = 1;			// additional coding
	nr_codec->harq.HARQ_iter = 0;	// additional coding
	nr_codec->ACK_cnt = 0;			// additional coding
	nr_codec->NACK_cnt = 0;			// additional coding

	nr_codec->harq.HARQ_index = new int[nr_codec->harq.HARQ_max_iter] ;	// additional coding
	int HARQ_index_temp[16] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3} ;	// can be manually modified		// additional coding

	for (int k=0; k<nr_codec->harq.HARQ_max_iter; k++){			// additional coding
		if (nr_codec->harq.HARQ_type == 1)						// additional coding
			nr_codec->harq.HARQ_index[k] = 0;					// additional coding
		else if (nr_codec->harq.HARQ_type == 2)					// additional coding
			nr_codec->harq.HARQ_index[k] = HARQ_index_temp[k] ;	// additional coding
	}



	Q = Q_m ;//# of coded bits

	// Data memory
	data_memory = new int[A] ;		// additional coding							

	// Code block segmenation
	tx_b = new int[B] ;

	rx_a = new int[A] ;
	rx_b = new int[B] ;
	//LDPC coding
	tx_c = new int*[C] ;
	rx_c = new double*[C] ;
	rx_c_LLR = new double*[C] ;
	for ( r=0; r<C; r++) {
		tx_c[r] = new int[K] ;
		rx_c[r] = new double[N_H] ;
		rx_c_LLR[r] = new double[N_H] ;
	}

	tx_d = new int*[C] ;
	rx_d = new double*[C] ;
	rx_d_memory = new double*[C] ;			// additional coding
	for ( r=0; r<C; r++) {
		tx_d[r] = new int[N] ;
		rx_d[r] = new double[N_H] ;
		rx_d_memory[r] = new double[N_H] ;	// additional coding	
	}
	for ( r=0; r<C; r++) {
		for (k=0; k<N; k++)
			tx_d[r][k] = -2 ;

		for (k=0; k<N_H; k++)
			rx_d[r][k] = 0 ;
	}
	// Rate matching 
	tx_e = new short*[C] ;
	rx_e = new double*[C] ;

	for ( r=0; r<C; r++) {
		tx_e[r] = new short[E[r]] ;
		rx_e[r] = new double[E[r]] ;
	}

	tx_f = new short*[C] ;
	rx_f = new double *[C] ;

	for (r=0; r<C; r++){
		tx_f[r] = new short[E[r]] ;
		rx_f[r] = new double[E[r]] ;
	}
	// Concatenation
	tx_g = new int[G] ;
	rx_g = new double[G] ;
	
	for (k=0;k<G;k++)	{
		tx_g[k] = 0 ;
		rx_g[k] = 0 ;
	}

	// LDPC decoding
	Lr = new double*[C];
	Lq = new double*[C];
	LLR_out = new double*[C];
	X_hat = new int*[C];
	HD_Lq = new int*[C] ;
	HD_Y = new int*[C] ;
	for(r=0;r<C; r++){
		Lr[r] = new double[M_H * max_row_weight];
		Lq[r] = new double[M_H * max_row_weight];
		LLR_out[r] = new double[N_H];
		X_hat[r] = new int[N_H];
		HD_Lq[r] = new int [M_H * max_row_weight] ;
		HD_Y[r] = new int [N_H] ;
	}

	MAX_ITER = nr_codec->ldpc.LDPCMaxIter;
	decoding_type = nr_codec->ldpc.LDPCDecOption;; // 0: sum-product,  1: min-sum

}


NR_ChannelCoding::NR_ChannelCoding()
{

}


NR_ChannelCoding::~NR_ChannelCoding()
{

}



/*-------------- Common public methods for Tx and Rx --------------*/

int NR_ChannelCoding::GetE ()
{
	return E[0] ;
}
int NR_ChannelCoding::GetN ()
{
	return N ;
}
int NR_ChannelCoding::GetG ()
{
	return G ;
}
int NR_ChannelCoding::GetK ()
{
	return K ;
}
int NR_ChannelCoding::GetBG ()
{
	return BG_index ;
}



/*-------------- Common private methods for Tx and Rx --------------*/

int* NR_ChannelCoding::Func_common_CRC_Calculation (int *BitData, int BitDataLength, unsigned int *crctable, int CRC_Length)
{
	int idx ;
	int ByteDataLength = BitDataLength/8 ;
	int *crc_bits = new int[CRC_Length] ;

	unsigned int byte ;
	unsigned int crc = 0 ;
	unsigned int *ByteData = new unsigned int[ByteDataLength] ;
	unsigned int *crc_bytes = new unsigned int[CRC_Length/8] ;
	unsigned int crc_ref ;
	
	switch (CRC_Length) {
	case 24:
		crc_ref = 0xFF0000L ;
		break ;
	case 16:
		crc_ref = 0xFF00L ;
		break ;
	case 8:
		crc_ref = 0xFFL ;
		break ;
	}


	// Bit-to-byte conversion
	for(int i=0; i<ByteDataLength; i++) {
		byte =   BitData[i*8  ] << 7  ;
		byte += (BitData[i*8+1] << 6) ;
		byte += (BitData[i*8+2] << 5) ;
		byte += (BitData[i*8+3] << 4) ;
		byte += (BitData[i*8+4] << 3) ;
		byte += (BitData[i*8+5] << 2) ;
		byte += (BitData[i*8+6] << 1) ;
		byte += (BitData[i*8+7]     ) ;

		ByteData[i]=byte ;
	}


	// CRC bytes calculation
	for(int j=0; j<ByteDataLength; j++) {
		idx = ((crc ^ (ByteData[j] << (CRC_Length-8))) >> (CRC_Length-8)) & 0xff ;

		switch (CRC_Length) {
		case 24:
			crc = (crctable[idx] ^ (crc<<8)) & 0xffffff ;
			break ;
		case 16:
			crc = (crctable[idx] ^ (crc<<8)) & 0xffff ;
			break ;
		case 8:
			crc = (crctable[idx] ^ (crc<<8)) & 0xff ;
			break ;
		}
	}

	for(int k=0; k<CRC_Length/8; k++) {
		crc_bytes[k] = (crc & crc_ref) >> (CRC_Length - 8*(k+1)) ;
		crc_ref = crc_ref >> 8 ;
	}


	// Byte-to-bit conversion
	for(int i=0; i<CRC_Length/8; i++) {
		crc_bits[8*i  ] = (crc_bytes[i] & 0x80) >> 7 ;
		crc_bits[8*i+1] = (crc_bytes[i] & 0x40) >> 6 ;
		crc_bits[8*i+2] = (crc_bytes[i] & 0x20) >> 5 ;
		crc_bits[8*i+3] = (crc_bytes[i] & 0x10) >> 4 ;
		crc_bits[8*i+4] = (crc_bytes[i] & 0x08) >> 3 ;
		crc_bits[8*i+5] = (crc_bytes[i] & 0x04) >> 2 ;
		crc_bits[8*i+6] = (crc_bytes[i] & 0x02) >> 1 ;
		crc_bits[8*i+7] = (crc_bytes[i] & 0x01) ;
	}


	// Delete memory
	delete []ByteData ;
	delete []crc_bytes ;

	return crc_bits ;
}



void NR_ChannelCoding::Func_common_GetCRCTable ()
{
	unsigned int crc_accum24a ;
	unsigned int crc_accum24b ;
	unsigned int crc_accum16 ;


	// CRC polynomials
	CRC24A = (0x864cfb) & 0xffffff ;
	CRC24B = (0x800063) & 0xffffff ;	
	CRC16  = (0x1021)   & 0xffff ;



	// Calculating CRC tables
	for (int i=0; i<256; i++) {
		crc_accum24a = ((unsigned int)i << 16) % 0xffffff ;
		crc_accum24b = ((unsigned int)i << 16) % 0xffffff ;
		crc_accum16  = ((unsigned int)i << 8)  % 0xffff ;
		
		for (int j=0; j<8; j++) {
			crc_accum24a = (crc_accum24a & 0x800000L)  ?  ((crc_accum24a << 1) ^ CRC24A) & 0xffffff  :  (crc_accum24a << 1) & 0xffffff ;
			crc_accum24b = (crc_accum24b & 0x800000L)  ?  ((crc_accum24b << 1) ^ CRC24B) & 0xffffff  :  (crc_accum24b << 1) & 0xffffff ;
			crc_accum16  = (crc_accum16  & 0x8000L)    ?  ((crc_accum16 << 1)  ^ CRC16)  & 0xffff    :  (crc_accum16  << 1) & 0xffff ;
		}
		
		crctable24a[i] = (crc_accum24a) & 0xffffff ;
		crctable24b[i] = (crc_accum24b) & 0xffffff ;
		crctable16[i]  = (crc_accum16)  & 0xffff ;
	}
}

void NR_ChannelCoding::Func_GetZ_table () {

	Z_table[0][0] = 2;	Z_table[0][1] = 4;	Z_table[0][2] = 8;	Z_table[0][3] = 16;	 Z_table[0][4] = 32;  Z_table[0][5] = 64 ;	Z_table[0][6] = 128 ;	Z_table[0][7] = 256 ;
	Z_table[1][0] = 3;	Z_table[1][1] = 6;	Z_table[1][2] = 12;	Z_table[1][3] = 24;	 Z_table[1][4] = 48;  Z_table[1][5] = 96;	Z_table[1][6] = 192;	Z_table[1][7] = 384;
	Z_table[2][0] = 5;	Z_table[2][1] = 10;	Z_table[2][2] = 20;	Z_table[2][3] = 40;	 Z_table[2][4] = 80;  Z_table[2][5] = 160;	Z_table[2][6] = 320;	Z_table[2][7] = 0;
	Z_table[3][0] = 7;	Z_table[3][1] = 14;	Z_table[3][2] = 28;	Z_table[3][3] = 56;	 Z_table[3][4] = 112; Z_table[3][5] = 224;	Z_table[3][6] = 0;		Z_table[3][7] = 0;
	Z_table[4][0] = 9;	Z_table[4][1] = 18;	Z_table[4][2] = 36;	Z_table[4][3] = 72;	 Z_table[4][4] = 144; Z_table[4][5] = 288;	Z_table[4][6] = 0;  	Z_table[4][7] = 0;
	Z_table[5][0] = 11;	Z_table[5][1] = 22;	Z_table[5][2] = 44;	Z_table[5][3] = 88;	 Z_table[5][4] = 176; Z_table[5][5] = 352;	Z_table[5][6] = 0;  	Z_table[5][7] = 0;
	Z_table[6][0] = 13;	Z_table[6][1] = 26;	Z_table[6][2] = 52;	Z_table[6][3] = 104; Z_table[6][4] = 208; Z_table[6][5] = 0;  	Z_table[6][6] = 0;  	Z_table[6][7] = 0;
	Z_table[7][0] = 15;	Z_table[7][1] = 30;	Z_table[7][2] = 60;	Z_table[7][3] = 120; Z_table[7][4] = 240; Z_table[7][5] = 0;	Z_table[7][6] = 0;  	Z_table[7][7] = 0;
}

void NR_ChannelCoding::Func_GetSegmentPara() {
	
	int rr; 
	int cc;

	int Z_c_temp = 1000 ; 
	int i_LS_temp = 0 ; 


	if (BG_index == 1)
		K_cb = 8448 ;
	else // (BG_index == 2)
		K_cb = 3840 ;

	// Total number of code block, C
	if (B <= K_cb) {
		L = 0 ;
		C = 1 ; 
		B_p = B ;
	}
	else {
		L = 24 ;
		C = (int) ceil((double) B / (double) (K_cb-L)) ;
		B_p = B + C*L ;
	}

	// The number of bits in each code block
	K_p = B_p/C ;

	if (BG_index == 1)
		K_b = 22 ;
	else if (BG_index == 2){
		if (B > 640)
			K_b = 10 ;
		else if (B > 560)
			K_b = 9 ;
		else if (B > 192)
			K_b = 8 ;
		else
			K_b = 6 ; 
	}
	
	// Determine the lifting size and the set index
	for ( rr=0; rr<8; rr++ )	{
		for( cc=0; cc<8; cc++) {
			if (K_b * Z_table[rr][cc] >= K_p){
				Z_c = Z_table[rr][cc] ;
				if (Z_c <= Z_c_temp){
					Z_c_temp = Z_c;
					i_LS_temp = rr + 1;
				}
			}
		}
	}

	Z_c = Z_c_temp;
	i_LS = i_LS_temp ;
	
	
	if (BG_index == 1){
		K_a = 22 ;
		K = K_a * Z_c ; 
	}
	else if (BG_index == 2){
		K_a = 10 ;
		K = 10 * Z_c ;
	}
}

void NR_ChannelCoding::Func_GetBaseGraph()
{
	int k, m;
	if (BG_index == 1){
		N_a = 66 ;
		M_a = N_a - K_a ;
		N = N_a*Z_c ;
		M = N - K ;

		N_h = 68;
		M_h = N_h - K_a ;
		N_H = N_h*Z_c ;
		M_H = N_H - K ;
		max_row_weight = 19 ;	max_col_weight = 30 ;

		// for submatrix of H
		N_h_AA = 22 ;
		M_h_AA = 4 ;
		N_H_AA = N_h_AA*Z_c ;
		//M_H_AA = M_h_AA*Z_c ;

		N_h_BB = 4 ;
		M_h_BB = 4 ;
		//N_H_BB = N_h_BB*Z_c ;
		M_H_BB = M_h_BB*Z_c ;

		N_h_DD = 26 ;
		M_h_DD = 42 ;
		//N_H_DD = N_h_DD*Z_c ;
		M_H_DD = M_h_DD*Z_c ;
		
		max_row_weight_AA = 17 ;	max_col_weight_AA = 4 ;
		//max_row_weight_BB = 3 ;	max_col_weight_BB = 3 ;
		max_row_weight_DD = 9 ;	max_col_weight_DD = 26 ;
	}
	else if (BG_index == 2){
		N_a = 50 ;
		M_a = N_a - K_a ;
		N = N_a*Z_c ;
		M = N - K ;

		N_h = 52;
		M_h = N_h - K_a ;
		N_H = N_h*Z_c ;
		M_H = N_H - K ;

		max_row_weight = 10 ;	max_col_weight = 23 ;	

		// for submatrix of H
		N_h_AA = 10 ;
		M_h_AA = 4 ;
		N_H_AA = N_h_AA*Z_c ;
		//M_H_AA = M_h_AA*Z_c ;

		N_h_BB = 4 ;
		M_h_BB = 4 ;
		//N_H_BB = N_h_BB*Z_c ;
		M_H_BB = M_h_BB*Z_c ;

		N_h_DD = 14 ;
		M_h_DD = 38 ;
		//N_H_DD = N_h_DD*Z_c ;
		M_H_DD = M_h_DD*Z_c ;
		
		max_row_weight_AA = 8 ;	max_col_weight_AA = 3 ;
		//max_row_weight_BB = 41 ;	max_col_weight_BB = 3 ;
		max_row_weight_DD = 5 ;	max_col_weight_DD = 20 ;
	}

	BaseGraph = new short*[M_h] ;
	for (k=0; k<M_h; k++)
		BaseGraph[k] = new short[N_h] ;
	
	for (k=0; k<M_h; k++){
		for(m=0; m<N_h; m++)
			BaseGraph[k][m] = -1 ;
	}
	k=0;
	if (BG_index == 1){
		while (k<316) {
			BaseGraph[ BaseGraph_info_mtx[k][0] ][ BaseGraph_info_mtx[k][1] ] = BaseGraph_info_mtx[k][2];
			k++;
		}		
	}
	else if (BG_index == 2){
		while (k<197) {
			BaseGraph[ BaseGraph_info_mtx[k][0] ][ BaseGraph_info_mtx[k][1] ] = BaseGraph_info_mtx[k][2];
			k++;
		}
	}
}

void NR_ChannelCoding::Func_GetWeight() {
	// Get weight
	int i, k;

	i = 0;
	for (k=0; k < max_row_weight * M_H; k++){
		if (col_index[k] != -1)						
			row_weight[i]++;
		if (k!=0 && (k+1)%max_row_weight == 0)		
			i++;				
	}
	i = 0;
	for (k=0; k < max_col_weight * N_H; k++){
		if (row_index[k] != -1)						
			col_weight[i]++;
		if (k!=0 && (k+1)%max_col_weight == 0)		
			i++;				
	}
}

short* NR_ChannelCoding::Func_GetRowIndex(short *row_index_vec, int BG_row_start,int BG_row_end, int Max_Row_Weight, int BG_col_start,int BG_col_end, int Max_Col_Weight ) {
	int i, j,  k, m, r;

	int cnt, tmp, *incc, **row_index_vec_temp ;
	
	incc = new int [Max_Col_Weight];

	row_index_vec_temp = new int*[(BG_col_end - BG_col_start) * Z_c] ;
	for(r=0;r<(BG_col_end - BG_col_start) * Z_c; r++)	
		row_index_vec_temp[r] = new int[Max_Col_Weight];
	
	for(r=0;r<(BG_col_end - BG_col_start) * Z_c; r++)	
		for(m=0; m<Max_Col_Weight; m++)
			row_index_vec_temp[r][m] = -1;
	
	for(j=BG_col_start;j<BG_col_end;j++){//?œìž‘ col index~??col index
        cnt=0;
        for(i=BG_row_start;i<BG_row_end;i++){//?œìž‘ row index~??row index
            if(BaseGraph[i][j]>-1){
                incc[cnt] = Z_c - BaseGraph[i][j] ;
				//cout<< cnt << " " ;
                for(k=0;k<Z_c;k++){
                    tmp=i*Z_c + incc[cnt]+k;
                    if(tmp < (i+1)*Z_c)		row_index_vec_temp[(j-BG_col_start)*Z_c+k][cnt] = tmp;
                    else					row_index_vec_temp[(j-BG_col_start)*Z_c+k][cnt] = tmp-Z_c;
                }
                cnt=cnt+1;
            }
        }
	}	

	for(k=0; k<(BG_col_end - BG_col_start) * Z_c; k++){
		for(m=0; m<Max_Col_Weight; m++) 	{
			row_index_vec[m+Max_Col_Weight*k]=row_index_vec_temp[k][m];	
		}
	}
		
	delete [] incc;
	for ( k=0; k<(BG_col_end - BG_col_start) * Z_c; k++) {
		delete []row_index_vec_temp[k] ;
	} 

	delete [] row_index_vec_temp;

	return row_index_vec ;

}


short* NR_ChannelCoding::Func_GetColIndex(short *col_index_vec, int BG_row_start,int BG_row_end, int Max_Row_Weight, int BG_col_start,int BG_col_end, int Max_Col_Weight ) {
	int i, j,  k, m, r;

	int cnt, tmp, *incr, **col_index_vec_temp ;
	
	incr = new int [Max_Row_Weight];

	col_index_vec_temp = new int*[(BG_row_end - BG_row_start) * Z_c] ;
	for(r=0;r<(BG_row_end - BG_row_start) * Z_c; r++)	
		col_index_vec_temp[r] = new int[Max_Row_Weight];
	
	for(r=0;r<(BG_row_end - BG_row_start) * Z_c; r++)	
		for(m=0; m<Max_Row_Weight; m++)
			col_index_vec_temp[r][m] = -1;

	
	// Get index of H
	for(i=BG_row_start;i<BG_row_end;i++){ //?œìž‘ row index~??row index
		cnt=0;
		for(j=BG_col_start;j<BG_col_end;j++){//?œìž‘ col index~??col index
			if(BaseGraph[i][j]>-1)	{
				incr[cnt] = BaseGraph[i][j] ;
				for(k=0;k<Z_c;k++)	{
					tmp=j*Z_c + incr[cnt]+k;
					if(tmp < (j+1)*Z_c)		
						col_index_vec_temp[(i - BG_row_start)*Z_c+k][cnt] = tmp;
					else					
						col_index_vec_temp[(i - BG_row_start)*Z_c+k][cnt] = tmp-Z_c;
				}
				cnt=cnt+1;
			}
		}
	}

	for(k=0; k< (BG_row_end - BG_row_start) * Z_c; k++){
		for(m=0; m<Max_Row_Weight; m++) 	col_index_vec[m+Max_Row_Weight*k]=col_index_vec_temp[k][m];		}

	
	delete [] incr;

	for ( k=0; k<(BG_row_end - BG_row_start) * Z_c; k++) {		delete []col_index_vec_temp[k] ;
	} 

	delete [] col_index_vec_temp;

	return col_index_vec ;


}

void NR_ChannelCoding::Func_GetIndexInvBB_AA() {
	
	int i, k, m ;

	// Get index
	i = 0;
	for (k=0; k<M_H_BB; k++){
		for (m=0; m<K; m++){
			if (inv_BB_AA[k][m] == 1){
				col_index_BB[i] = m ;
				i++; 
			}
		}
		i = (k+1)*max_row_weight_BB ;		
	}
}

void NR_ChannelCoding::Func_GetBBMatrix (){
	int k,m,i,j;

	for (k=0; k<M_h_BB; k++){
		for (m=N_h_AA; m<N_h_AA+N_h_BB; m++){
			if (BaseGraph[k][m] != -1){
				cyc_perm_mtx = Func_GetCyclePermMatrix(cyc_perm_mtx, Z_c, BaseGraph[k][m]) ;

				for (i=0; i<Z_c; i++){
					for (j=0; j<Z_c; j++)
						BB[k*Z_c+i][(m-N_h_AA)*Z_c+j] = cyc_perm_mtx[i][j] ;
				}
			}
		}
	}
}




short** NR_ChannelCoding::Func_GetCyclePermMatrix (short** matrix, int mtx_size, int shift_idx){

	//int** matrix;
	int aa ;
	int bb ;

	for (aa=0; aa<mtx_size; aa++){
		for(bb=0; bb<mtx_size; bb++){
			if (bb == (shift_idx+aa) % mtx_size)
				matrix[aa][bb] = 1 ;
			else
				matrix[aa][bb] = 0 ;
		}
	}
	
	return matrix;
}

void NR_ChannelCoding::Func_GetInvBB_AA () {

		for (int k=0;k<4*Z_c; k++){
		for (int m=0;m<K; m++){
			for (int i=0; i<max_col_weight_AA; i++){
				if (row_index_AA[i+max_col_weight_AA*m] != -1){
					inv_BB_AA[k][m] = (inv_BB_AA[k][m] + inv_BB[k][row_index_AA[i+max_col_weight_AA*m]]) % 2 ;
				}
				else
					break;
			}
		}
	}
}
void NR_ChannelCoding::Func_GetInvBB_AA_max_row_weight(){
	int max_row_weight_BB_temp = 0;

	max_row_weight_BB = 0;

	for (int k=0;k<4*Z_c; k++){
		for (int m=0;m<K; m++){
			max_row_weight_BB_temp = max_row_weight_BB_temp + inv_BB_AA[k][m];
		}
		if(max_row_weight_BB <= max_row_weight_BB_temp){
			max_row_weight_BB = max_row_weight_BB_temp ;
		}		
		max_row_weight_BB_temp = 0 ;
	}
}

short ** NR_ChannelCoding::Func_GetInvMatrix (short** mtx){

	int k ;
	int m ;
	int i ;
	int index ;
	short** eye_mtx ;
	short** mtx_1 ;

	eye_mtx = new short*[4*Z_c] ;
	for (k=0; k<4*Z_c; k++)
		eye_mtx[k] = new short[4*Z_c] ;

	mtx_1 = new short*[4*Z_c] ;
	for (k=0; k<4*Z_c; k++)
		mtx_1[k] = new short[8*Z_c] ;
	

	int aa ;
	int bb ;
	for (aa=0; aa<4*Z_c; aa++){
		for(bb=0; bb<4*Z_c; bb++){
			if (bb == aa % (4*Z_c))
				eye_mtx[aa][bb] = 1 ;
			else
				eye_mtx[aa][bb] = 0 ;
		}
	}

	for (k=0; k<4*Z_c; k++){
		for (m=0; m<8*Z_c; m++){
			if (m<4*Z_c)
				mtx_1[k][m] = mtx[k][m] ;
			else
				mtx_1[k][m] = eye_mtx[k][m-4*Z_c] ;
		}
	}

	for (i=0; i<4*Z_c; i++){
		
		for (k=i; k<4*Z_c; k++){
			if (mtx_1[k][i] == 1){
				index = k;
				break;
			}
		}

		if (index != i){
			for (m=0; m<8*Z_c; m++)
				mtx_1[i][m] = (mtx_1[i][m] + mtx_1[index][m]) % 2 ;
		}

		for (k=0; k<4*Z_c; k++){
			if (k != i && mtx_1[k][i] == 1){
				for (m=i; m<8*Z_c; m++)
					mtx_1[k][m] = (mtx_1[k][m] + mtx_1[i][m]) % 2 ;
			}
		}
	}


	for (k=0; k<4*Z_c; k++){
		for (m=0; m<4*Z_c; m++)
			mtx[k][m] = mtx_1[k][m+4*Z_c] ;
	}

	for (int k=0; k<4*Z_c; k++) {
		delete []eye_mtx[k] ;
		delete []mtx_1[k] ;
	} 
	delete eye_mtx ;
	delete mtx_1 ;

	return mtx;
}

short** NR_ChannelCoding::Func_GetBaseGraphTable()
{
	int BG_table_1[316][10] = {
		{	0	,	0	,	250	,	307	,	73	,	223	,	211	,	294	,	0	,	135	},
		{	0	,	1	,	69	,	19	,	15	,	16	,	198	,	118	,	0	,	227	},
		{	0	,	2	,	226	,	50	,	103	,	94	,	188	,	167	,	0	,	126	},
		{	0	,	3	,	159	,	369	,	49	,	91	,	186	,	330	,	0	,	134	},
		{	0	,	5	,	100	,	181	,	240	,	74	,	219	,	207	,	0	,	84	},
		{	0	,	6	,	10	,	216	,	39	,	10	,	4	,	165	,	0	,	83	},
		{	0	,	9	,	59	,	317	,	15	,	0	,	29	,	243	,	0	,	53	},
		{	0	,	10	,	229	,	288	,	162	,	205	,	144	,	250	,	0	,	225	},
		{	0	,	11	,	110	,	109	,	215	,	216	,	116	,	1	,	0	,	205	},
		{	0	,	12	,	191	,	17	,	164	,	21	,	216	,	339	,	0	,	128	},
		{	0	,	13	,	9	,	357	,	133	,	215	,	115	,	201	,	0	,	75	},
		{	0	,	15	,	195	,	215	,	298	,	14	,	233	,	53	,	0	,	135	},
		{	0	,	16	,	23	,	106	,	110	,	70	,	144	,	347	,	0	,	217	},
		{	0	,	18	,	190	,	242	,	113	,	141	,	95	,	304	,	0	,	220	},
		{	0	,	19	,	35	,	180	,	16	,	198	,	216	,	167	,	0	,	90	},
		{	0	,	20	,	239	,	330	,	189	,	104	,	73	,	47	,	0	,	105	},
		{	0	,	21	,	31	,	346	,	32	,	81	,	261	,	188	,	0	,	137	},
		{	0	,	22	,	1	,	1	,	1	,	1	,	1	,	1	,	0	,	1	},
		{	0	,	23	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	1	,	0	,	2	,	76	,	303	,	141	,	179	,	77	,	22	,	96	},
		{	1	,	2	,	239	,	76	,	294	,	45	,	162	,	225	,	11	,	236	},
		{	1	,	3	,	117	,	73	,	27	,	151	,	223	,	96	,	124	,	136	},
		{	1	,	4	,	124	,	288	,	261	,	46	,	256	,	338	,	0	,	221	},
		{	1	,	5	,	71	,	144	,	161	,	119	,	160	,	268	,	10	,	128	},
		{	1	,	7	,	222	,	331	,	133	,	157	,	76	,	112	,	0	,	92	},
		{	1	,	8	,	104	,	331	,	4	,	133	,	202	,	302	,	0	,	172	},
		{	1	,	9	,	173	,	178	,	80	,	87	,	117	,	50	,	2	,	56	},
		{	1	,	11	,	220	,	295	,	129	,	206	,	109	,	167	,	16	,	11	},
		{	1	,	12	,	102	,	342	,	300	,	93	,	15	,	253	,	60	,	189	},
		{	1	,	14	,	109	,	217	,	76	,	79	,	72	,	334	,	0	,	95	},
		{	1	,	15	,	132	,	99	,	266	,	9	,	152	,	242	,	6	,	85	},
		{	1	,	16	,	142	,	354	,	72	,	118	,	158	,	257	,	30	,	153	},
		{	1	,	17	,	155	,	114	,	83	,	194	,	147	,	133	,	0	,	87	},
		{	1	,	19	,	255	,	331	,	260	,	31	,	156	,	9	,	168	,	163	},
		{	1	,	21	,	28	,	112	,	301	,	187	,	119	,	302	,	31	,	216	},
		{	1	,	22	,	0	,	0	,	0	,	0	,	0	,	0	,	105	,	0	},
		{	1	,	23	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	1	,	24	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	2	,	0	,	106	,	205	,	68	,	207	,	258	,	226	,	132	,	189	},
		{	2	,	1	,	111	,	250	,	7	,	203	,	167	,	35	,	37	,	4	},
		{	2	,	2	,	185	,	328	,	80	,	31	,	220	,	213	,	21	,	225	},
		{	2	,	4	,	63	,	332	,	280	,	176	,	133	,	302	,	180	,	151	},
		{	2	,	5	,	117	,	256	,	38	,	180	,	243	,	111	,	4	,	236	},
		{	2	,	6	,	93	,	161	,	227	,	186	,	202	,	265	,	149	,	117	},
		{	2	,	7	,	229	,	267	,	202	,	95	,	218	,	128	,	48	,	179	},
		{	2	,	8	,	177	,	160	,	200	,	153	,	63	,	237	,	38	,	92	},
		{	2	,	9	,	95	,	63	,	71	,	177	,	0	,	294	,	122	,	24	},
		{	2	,	10	,	39	,	129	,	106	,	70	,	3	,	127	,	195	,	68	},
		{	2	,	13	,	142	,	200	,	295	,	77	,	74	,	110	,	155	,	6	},
		{	2	,	14	,	225	,	88	,	283	,	214	,	229	,	286	,	28	,	101	},
		{	2	,	15	,	225	,	53	,	301	,	77	,	0	,	125	,	85	,	33	},
		{	2	,	17	,	245	,	131	,	184	,	198	,	216	,	131	,	47	,	96	},
		{	2	,	18	,	205	,	240	,	246	,	117	,	269	,	163	,	179	,	125	},
		{	2	,	19	,	251	,	205	,	230	,	223	,	200	,	210	,	42	,	67	},
		{	2	,	20	,	117	,	13	,	276	,	90	,	234	,	7	,	66	,	230	},
		{	2	,	24	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	2	,	25	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	3	,	0	,	121	,	276	,	220	,	201	,	187	,	97	,	4	,	128	},
		{	3	,	1	,	89	,	87	,	208	,	18	,	145	,	94	,	6	,	23	},
		{	3	,	3	,	84	,	0	,	30	,	165	,	166	,	49	,	33	,	162	},
		{	3	,	4	,	20	,	275	,	197	,	5	,	108	,	279	,	113	,	220	},
		{	3	,	6	,	150	,	199	,	61	,	45	,	82	,	139	,	49	,	43	},
		{	3	,	7	,	131	,	153	,	175	,	142	,	132	,	166	,	21	,	186	},
		{	3	,	8	,	243	,	56	,	79	,	16	,	197	,	91	,	6	,	96	},
		{	3	,	10	,	136	,	132	,	281	,	34	,	41	,	106	,	151	,	1	},
		{	3	,	11	,	86	,	305	,	303	,	155	,	162	,	246	,	83	,	216	},
		{	3	,	12	,	246	,	231	,	253	,	213	,	57	,	345	,	154	,	22	},
		{	3	,	13	,	219	,	341	,	164	,	147	,	36	,	269	,	87	,	24	},
		{	3	,	14	,	211	,	212	,	53	,	69	,	115	,	185	,	5	,	167	},
		{	3	,	16	,	240	,	304	,	44	,	96	,	242	,	249	,	92	,	200	},
		{	3	,	17	,	76	,	300	,	28	,	74	,	165	,	215	,	173	,	32	},
		{	3	,	18	,	244	,	271	,	77	,	99	,	0	,	143	,	120	,	235	},
		{	3	,	20	,	144	,	39	,	319	,	30	,	113	,	121	,	2	,	172	},
		{	3	,	21	,	12	,	357	,	68	,	158	,	108	,	121	,	142	,	219	},
		{	3	,	22	,	1	,	1	,	1	,	1	,	1	,	1	,	0	,	1	},
		{	3	,	25	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	4	,	0	,	157	,	332	,	233	,	170	,	246	,	42	,	24	,	64	},
		{	4	,	1	,	102	,	181	,	205	,	10	,	235	,	256	,	204	,	211	},
		{	4	,	26	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	5	,	0	,	205	,	195	,	83	,	164	,	261	,	219	,	185	,	2	},
		{	5	,	1	,	236	,	14	,	292	,	59	,	181	,	130	,	100	,	171	},
		{	5	,	3	,	194	,	115	,	50	,	86	,	72	,	251	,	24	,	47	},
		{	5	,	12	,	231	,	166	,	318	,	80	,	283	,	322	,	65	,	143	},
		{	5	,	16	,	28	,	241	,	201	,	182	,	254	,	295	,	207	,	210	},
		{	5	,	21	,	123	,	51	,	267	,	130	,	79	,	258	,	161	,	180	},
		{	5	,	22	,	115	,	157	,	279	,	153	,	144	,	283	,	72	,	180	},
		{	5	,	27	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	6	,	0	,	183	,	278	,	289	,	158	,	80	,	294	,	6	,	199	},
		{	6	,	6	,	22	,	257	,	21	,	119	,	144	,	73	,	27	,	22	},
		{	6	,	10	,	28	,	1	,	293	,	113	,	169	,	330	,	163	,	23	},
		{	6	,	11	,	67	,	351	,	13	,	21	,	90	,	99	,	50	,	100	},
		{	6	,	13	,	244	,	92	,	232	,	63	,	59	,	172	,	48	,	92	},
		{	6	,	17	,	11	,	253	,	302	,	51	,	177	,	150	,	24	,	207	},
		{	6	,	18	,	157	,	18	,	138	,	136	,	151	,	284	,	38	,	52	},
		{	6	,	20	,	211	,	225	,	235	,	116	,	108	,	305	,	91	,	13	},
		{	6	,	28	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	7	,	0	,	220	,	9	,	12	,	17	,	169	,	3	,	145	,	77	},
		{	7	,	1	,	44	,	62	,	88	,	76	,	189	,	103	,	88	,	146	},
		{	7	,	4	,	159	,	316	,	207	,	104	,	154	,	224	,	112	,	209	},
		{	7	,	7	,	31	,	333	,	50	,	100	,	184	,	297	,	153	,	32	},
		{	7	,	8	,	167	,	290	,	25	,	150	,	104	,	215	,	159	,	166	},
		{	7	,	14	,	104	,	114	,	76	,	158	,	164	,	39	,	76	,	18	},
		{	7	,	29	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	8	,	0	,	112	,	307	,	295	,	33	,	54	,	348	,	172	,	181	},
		{	8	,	1	,	4	,	179	,	133	,	95	,	0	,	75	,	2	,	105	},
		{	8	,	3	,	7	,	165	,	130	,	4	,	252	,	22	,	131	,	141	},
		{	8	,	12	,	211	,	18	,	231	,	217	,	41	,	312	,	141	,	223	},
		{	8	,	16	,	102	,	39	,	296	,	204	,	98	,	224	,	96	,	177	},
		{	8	,	19	,	164	,	224	,	110	,	39	,	46	,	17	,	99	,	145	},
		{	8	,	21	,	109	,	368	,	269	,	58	,	15	,	59	,	101	,	199	},
		{	8	,	22	,	241	,	67	,	245	,	44	,	230	,	314	,	35	,	153	},
		{	8	,	24	,	90	,	170	,	154	,	201	,	54	,	244	,	116	,	38	},
		{	8	,	30	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	9	,	0	,	103	,	366	,	189	,	9	,	162	,	156	,	6	,	169	},
		{	9	,	1	,	182	,	232	,	244	,	37	,	159	,	88	,	10	,	12	},
		{	9	,	10	,	109	,	321	,	36	,	213	,	93	,	293	,	145	,	206	},
		{	9	,	11	,	21	,	133	,	286	,	105	,	134	,	111	,	53	,	221	},
		{	9	,	13	,	142	,	57	,	151	,	89	,	45	,	92	,	201	,	17	},
		{	9	,	17	,	14	,	303	,	267	,	185	,	132	,	152	,	4	,	212	},
		{	9	,	18	,	61	,	63	,	135	,	109	,	76	,	23	,	164	,	92	},
		{	9	,	20	,	216	,	82	,	209	,	218	,	209	,	337	,	173	,	205	},
		{	9	,	31	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	10	,	1	,	98	,	101	,	14	,	82	,	178	,	175	,	126	,	116	},
		{	10	,	2	,	149	,	339	,	80	,	165	,	1	,	253	,	77	,	151	},
		{	10	,	4	,	167	,	274	,	211	,	174	,	28	,	27	,	156	,	70	},
		{	10	,	7	,	160	,	111	,	75	,	19	,	267	,	231	,	16	,	230	},
		{	10	,	8	,	49	,	383	,	161	,	194	,	234	,	49	,	12	,	115	},
		{	10	,	14	,	58	,	354	,	311	,	103	,	201	,	267	,	70	,	84	},
		{	10	,	32	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	11	,	0	,	77	,	48	,	16	,	52	,	55	,	25	,	184	,	45	},
		{	11	,	1	,	41	,	102	,	147	,	11	,	23	,	322	,	194	,	115	},
		{	11	,	12	,	83	,	8	,	290	,	2	,	274	,	200	,	123	,	134	},
		{	11	,	16	,	182	,	47	,	289	,	35	,	181	,	351	,	16	,	1	},
		{	11	,	21	,	78	,	188	,	177	,	32	,	273	,	166	,	104	,	152	},
		{	11	,	22	,	252	,	334	,	43	,	84	,	39	,	338	,	109	,	165	},
		{	11	,	23	,	22	,	115	,	280	,	201	,	26	,	192	,	124	,	107	},
		{	11	,	33	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	12	,	0	,	160	,	77	,	229	,	142	,	225	,	123	,	6	,	186	},
		{	12	,	1	,	42	,	186	,	235	,	175	,	162	,	217	,	20	,	215	},
		{	12	,	10	,	21	,	174	,	169	,	136	,	244	,	142	,	203	,	124	},
		{	12	,	11	,	32	,	232	,	48	,	3	,	151	,	110	,	153	,	180	},
		{	12	,	13	,	234	,	50	,	105	,	28	,	238	,	176	,	104	,	98	},
		{	12	,	18	,	7	,	74	,	52	,	182	,	243	,	76	,	207	,	80	},
		{	12	,	34	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	13	,	0	,	177	,	313	,	39	,	81	,	231	,	311	,	52	,	220	},
		{	13	,	3	,	248	,	177	,	302	,	56	,	0	,	251	,	147	,	185	},
		{	13	,	7	,	151	,	266	,	303	,	72	,	216	,	265	,	1	,	154	},
		{	13	,	20	,	185	,	115	,	160	,	217	,	47	,	94	,	16	,	178	},
		{	13	,	23	,	62	,	370	,	37	,	78	,	36	,	81	,	46	,	150	},
		{	13	,	35	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	14	,	0	,	206	,	142	,	78	,	14	,	0	,	22	,	1	,	124	},
		{	14	,	12	,	55	,	248	,	299	,	175	,	186	,	322	,	202	,	144	},
		{	14	,	15	,	206	,	137	,	54	,	211	,	253	,	277	,	118	,	182	},
		{	14	,	16	,	127	,	89	,	61	,	191	,	16	,	156	,	130	,	95	},
		{	14	,	17	,	16	,	347	,	179	,	51	,	0	,	66	,	1	,	72	},
		{	14	,	21	,	229	,	12	,	258	,	43	,	79	,	78	,	2	,	76	},
		{	14	,	36	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	15	,	0	,	40	,	241	,	229	,	90	,	170	,	176	,	173	,	39	},
		{	15	,	1	,	96	,	2	,	290	,	120	,	0	,	348	,	6	,	138	},
		{	15	,	10	,	65	,	210	,	60	,	131	,	183	,	15	,	81	,	220	},
		{	15	,	13	,	63	,	318	,	130	,	209	,	108	,	81	,	182	,	173	},
		{	15	,	18	,	75	,	55	,	184	,	209	,	68	,	176	,	53	,	142	},
		{	15	,	25	,	179	,	269	,	51	,	81	,	64	,	113	,	46	,	49	},
		{	15	,	37	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	16	,	1	,	64	,	13	,	69	,	154	,	270	,	190	,	88	,	78	},
		{	16	,	3	,	49	,	338	,	140	,	164	,	13	,	293	,	198	,	152	},
		{	16	,	11	,	49	,	57	,	45	,	43	,	99	,	332	,	160	,	84	},
		{	16	,	20	,	51	,	289	,	115	,	189	,	54	,	331	,	122	,	5	},
		{	16	,	22	,	154	,	57	,	300	,	101	,	0	,	114	,	182	,	205	},
		{	16	,	38	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	17	,	0	,	7	,	260	,	257	,	56	,	153	,	110	,	91	,	183	},
		{	17	,	14	,	164	,	303	,	147	,	110	,	137	,	228	,	184	,	112	},
		{	17	,	16	,	59	,	81	,	128	,	200	,	0	,	247	,	30	,	106	},
		{	17	,	17	,	1	,	358	,	51	,	63	,	0	,	116	,	3	,	219	},
		{	17	,	21	,	144	,	375	,	228	,	4	,	162	,	190	,	155	,	129	},
		{	17	,	39	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	18	,	1	,	42	,	130	,	260	,	199	,	161	,	47	,	1	,	183	},
		{	18	,	12	,	233	,	163	,	294	,	110	,	151	,	286	,	41	,	215	},
		{	18	,	13	,	8	,	280	,	291	,	200	,	0	,	246	,	167	,	180	},
		{	18	,	18	,	155	,	132	,	141	,	143	,	241	,	181	,	68	,	143	},
		{	18	,	19	,	147	,	4	,	295	,	186	,	144	,	73	,	148	,	14	},
		{	18	,	40	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	19	,	0	,	60	,	145	,	64	,	8	,	0	,	87	,	12	,	179	},
		{	19	,	1	,	73	,	213	,	181	,	6	,	0	,	110	,	6	,	108	},
		{	19	,	7	,	72	,	344	,	101	,	103	,	118	,	147	,	166	,	159	},
		{	19	,	8	,	127	,	242	,	270	,	198	,	144	,	258	,	184	,	138	},
		{	19	,	10	,	224	,	197	,	41	,	8	,	0	,	204	,	191	,	196	},
		{	19	,	41	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	20	,	0	,	151	,	187	,	301	,	105	,	265	,	89	,	6	,	77	},
		{	20	,	3	,	186	,	206	,	162	,	210	,	81	,	65	,	12	,	187	},
		{	20	,	9	,	217	,	264	,	40	,	121	,	90	,	155	,	15	,	203	},
		{	20	,	11	,	47	,	341	,	130	,	214	,	144	,	244	,	5	,	167	},
		{	20	,	22	,	160	,	59	,	10	,	183	,	228	,	30	,	30	,	130	},
		{	20	,	42	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	21	,	1	,	249	,	205	,	79	,	192	,	64	,	162	,	6	,	197	},
		{	21	,	5	,	121	,	102	,	175	,	131	,	46	,	264	,	86	,	122	},
		{	21	,	16	,	109	,	328	,	132	,	220	,	266	,	346	,	96	,	215	},
		{	21	,	20	,	131	,	213	,	283	,	50	,	9	,	143	,	42	,	65	},
		{	21	,	21	,	171	,	97	,	103	,	106	,	18	,	109	,	199	,	216	},
		{	21	,	43	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	22	,	0	,	64	,	30	,	177	,	53	,	72	,	280	,	44	,	25	},
		{	22	,	12	,	142	,	11	,	20	,	0	,	189	,	157	,	58	,	47	},
		{	22	,	13	,	188	,	233	,	55	,	3	,	72	,	236	,	130	,	126	},
		{	22	,	17	,	158	,	22	,	316	,	148	,	257	,	113	,	131	,	178	},
		{	22	,	44	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	23	,	1	,	156	,	24	,	249	,	88	,	180	,	18	,	45	,	185	},
		{	23	,	2	,	147	,	89	,	50	,	203	,	0	,	6	,	18	,	127	},
		{	23	,	10	,	170	,	61	,	133	,	168	,	0	,	181	,	132	,	117	},
		{	23	,	18	,	152	,	27	,	105	,	122	,	165	,	304	,	100	,	199	},
		{	23	,	45	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	24	,	0	,	112	,	298	,	289	,	49	,	236	,	38	,	9	,	32	},
		{	24	,	3	,	86	,	158	,	280	,	157	,	199	,	170	,	125	,	178	},
		{	24	,	4	,	236	,	235	,	110	,	64	,	0	,	249	,	191	,	2	},
		{	24	,	11	,	116	,	339	,	187	,	193	,	266	,	288	,	28	,	156	},
		{	24	,	22	,	222	,	234	,	281	,	124	,	0	,	194	,	6	,	58	},
		{	24	,	46	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	25	,	1	,	23	,	72	,	172	,	1	,	205	,	279	,	4	,	27	},
		{	25	,	6	,	136	,	17	,	295	,	166	,	0	,	255	,	74	,	141	},
		{	25	,	7	,	116	,	383	,	96	,	65	,	0	,	111	,	16	,	11	},
		{	25	,	14	,	182	,	312	,	46	,	81	,	183	,	54	,	28	,	181	},
		{	25	,	47	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	26	,	0	,	195	,	71	,	270	,	107	,	0	,	325	,	21	,	163	},
		{	26	,	2	,	243	,	81	,	110	,	176	,	0	,	326	,	142	,	131	},
		{	26	,	4	,	215	,	76	,	318	,	212	,	0	,	226	,	192	,	169	},
		{	26	,	15	,	61	,	136	,	67	,	127	,	277	,	99	,	197	,	98	},
		{	26	,	48	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	27	,	1	,	25	,	194	,	210	,	208	,	45	,	91	,	98	,	165	},
		{	27	,	6	,	104	,	194	,	29	,	141	,	36	,	326	,	140	,	232	},
		{	27	,	8	,	194	,	101	,	304	,	174	,	72	,	268	,	22	,	9	},
		{	27	,	49	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	28	,	0	,	128	,	222	,	11	,	146	,	275	,	102	,	4	,	32	},
		{	28	,	4	,	165	,	19	,	293	,	153	,	0	,	1	,	1	,	43	},
		{	28	,	19	,	181	,	244	,	50	,	217	,	155	,	40	,	40	,	200	},
		{	28	,	21	,	63	,	274	,	234	,	114	,	62	,	167	,	93	,	205	},
		{	28	,	50	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	29	,	1	,	86	,	252	,	27	,	150	,	0	,	273	,	92	,	232	},
		{	29	,	14	,	236	,	5	,	308	,	11	,	180	,	104	,	136	,	32	},
		{	29	,	18	,	84	,	147	,	117	,	53	,	0	,	243	,	106	,	118	},
		{	29	,	25	,	6	,	78	,	29	,	68	,	42	,	107	,	6	,	103	},
		{	29	,	51	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	30	,	0	,	216	,	159	,	91	,	34	,	0	,	171	,	2	,	170	},
		{	30	,	10	,	73	,	229	,	23	,	130	,	90	,	16	,	88	,	199	},
		{	30	,	13	,	120	,	260	,	105	,	210	,	252	,	95	,	112	,	26	},
		{	30	,	24	,	9	,	90	,	135	,	123	,	173	,	212	,	20	,	105	},
		{	30	,	52	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	31	,	1	,	95	,	100	,	222	,	175	,	144	,	101	,	4	,	73	},
		{	31	,	7	,	177	,	215	,	308	,	49	,	144	,	297	,	49	,	149	},
		{	31	,	22	,	172	,	258	,	66	,	177	,	166	,	279	,	125	,	175	},
		{	31	,	25	,	61	,	256	,	162	,	128	,	19	,	222	,	194	,	108	},
		{	31	,	53	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	32	,	0	,	221	,	102	,	210	,	192	,	0	,	351	,	6	,	103	},
		{	32	,	12	,	112	,	201	,	22	,	209	,	211	,	265	,	126	,	110	},
		{	32	,	14	,	199	,	175	,	271	,	58	,	36	,	338	,	63	,	151	},
		{	32	,	24	,	121	,	287	,	217	,	30	,	162	,	83	,	20	,	211	},
		{	32	,	54	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	33	,	1	,	2	,	323	,	170	,	114	,	0	,	56	,	10	,	199	},
		{	33	,	2	,	187	,	8	,	20	,	49	,	0	,	304	,	30	,	132	},
		{	33	,	11	,	41	,	361	,	140	,	161	,	76	,	141	,	6	,	172	},
		{	33	,	21	,	211	,	105	,	33	,	137	,	18	,	101	,	92	,	65	},
		{	33	,	55	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	34	,	0	,	127	,	230	,	187	,	82	,	197	,	60	,	4	,	161	},
		{	34	,	7	,	167	,	148	,	296	,	186	,	0	,	320	,	153	,	237	},
		{	34	,	15	,	164	,	202	,	5	,	68	,	108	,	112	,	197	,	142	},
		{	34	,	17	,	159	,	312	,	44	,	150	,	0	,	54	,	155	,	180	},
		{	34	,	56	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	35	,	1	,	161	,	320	,	207	,	192	,	199	,	100	,	4	,	231	},
		{	35	,	6	,	197	,	335	,	158	,	173	,	278	,	210	,	45	,	174	},
		{	35	,	12	,	207	,	2	,	55	,	26	,	0	,	195	,	168	,	145	},
		{	35	,	22	,	103	,	266	,	285	,	187	,	205	,	268	,	185	,	100	},
		{	35	,	57	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	36	,	0	,	37	,	210	,	259	,	222	,	216	,	135	,	6	,	11	},
		{	36	,	14	,	105	,	313	,	179	,	157	,	16	,	15	,	200	,	207	},
		{	36	,	15	,	51	,	297	,	178	,	0	,	0	,	35	,	177	,	42	},
		{	36	,	18	,	120	,	21	,	160	,	6	,	0	,	188	,	43	,	100	},
		{	36	,	58	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	37	,	1	,	198	,	269	,	298	,	81	,	72	,	319	,	82	,	59	},
		{	37	,	13	,	220	,	82	,	15	,	195	,	144	,	236	,	2	,	204	},
		{	37	,	23	,	122	,	115	,	115	,	138	,	0	,	85	,	135	,	161	},
		{	37	,	59	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	38	,	0	,	167	,	185	,	151	,	123	,	190	,	164	,	91	,	121	},
		{	38	,	9	,	151	,	177	,	179	,	90	,	0	,	196	,	64	,	90	},
		{	38	,	10	,	157	,	289	,	64	,	73	,	0	,	209	,	198	,	26	},
		{	38	,	12	,	163	,	214	,	181	,	10	,	0	,	246	,	100	,	140	},
		{	38	,	60	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	39	,	1	,	173	,	258	,	102	,	12	,	153	,	236	,	4	,	115	},
		{	39	,	3	,	139	,	93	,	77	,	77	,	0	,	264	,	28	,	188	},
		{	39	,	7	,	149	,	346	,	192	,	49	,	165	,	37	,	109	,	168	},
		{	39	,	19	,	0	,	297	,	208	,	114	,	117	,	272	,	188	,	52	},
		{	39	,	61	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	40	,	0	,	157	,	175	,	32	,	67	,	216	,	304	,	10	,	4	},
		{	40	,	8	,	137	,	37	,	80	,	45	,	144	,	237	,	84	,	103	},
		{	40	,	17	,	149	,	312	,	197	,	96	,	2	,	135	,	12	,	30	},
		{	40	,	62	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	41	,	1	,	167	,	52	,	154	,	23	,	0	,	123	,	2	,	53	},
		{	41	,	3	,	173	,	314	,	47	,	215	,	0	,	77	,	75	,	189	},
		{	41	,	9	,	139	,	139	,	124	,	60	,	0	,	25	,	142	,	215	},
		{	41	,	18	,	151	,	288	,	207	,	167	,	183	,	272	,	128	,	24	},
		{	41	,	63	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	42	,	0	,	149	,	113	,	226	,	114	,	27	,	288	,	163	,	222	},
		{	42	,	4	,	157	,	14	,	65	,	91	,	0	,	83	,	10	,	170	},
		{	42	,	24	,	137	,	218	,	126	,	78	,	35	,	17	,	162	,	71	},
		{	42	,	64	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	43	,	1	,	151	,	113	,	228	,	206	,	52	,	210	,	1	,	22	},
		{	43	,	16	,	163	,	132	,	69	,	22	,	243	,	3	,	163	,	127	},
		{	43	,	18	,	173	,	114	,	176	,	134	,	0	,	53	,	99	,	49	},
		{	43	,	25	,	139	,	168	,	102	,	161	,	270	,	167	,	98	,	125	},
		{	43	,	65	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	44	,	0	,	139	,	80	,	234	,	84	,	18	,	79	,	4	,	191	},
		{	44	,	7	,	157	,	78	,	227	,	4	,	0	,	244	,	6	,	211	},
		{	44	,	9	,	163	,	163	,	259	,	9	,	0	,	293	,	142	,	187	},
		{	44	,	22	,	173	,	274	,	260	,	12	,	57	,	272	,	3	,	148	},
		{	44	,	66	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
		{	45	,	1	,	149	,	135	,	101	,	184	,	168	,	82	,	181	,	177	},
		{	45	,	6	,	151	,	149	,	228	,	121	,	0	,	67	,	45	,	114	},
		{	45	,	10	,	167	,	15	,	126	,	29	,	144	,	235	,	153	,	93	},
		{	45	,	67	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	}};



		int BG_table_2[197][10] = {
			{	0	,	0	,	9	,	174	,	0	,	72	,	3	,	156	,	143	,	145	},
			{	0	,	1	,	117	,	97	,	0	,	110	,	26	,	143	,	19	,	131	},
			{	0	,	2	,	204	,	166	,	0	,	23	,	53	,	14	,	176	,	71	},
			{	0	,	3	,	26	,	66	,	0	,	181	,	35	,	3	,	165	,	21	},
			{	0	,	6	,	189	,	71	,	0	,	95	,	115	,	40	,	196	,	23	},
			{	0	,	9	,	205	,	172	,	0	,	8	,	127	,	123	,	13	,	112	},
			{	0	,	10	,	0	,	0	,	0	,	1	,	0	,	0	,	0	,	1	},
			{	0	,	11	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	1	,	0	,	167	,	27	,	137	,	53	,	19	,	17	,	18	,	142	},
			{	1	,	3	,	166	,	36	,	124	,	156	,	94	,	65	,	27	,	174	},
			{	1	,	4	,	253	,	48	,	0	,	115	,	104	,	63	,	3	,	183	},
			{	1	,	5	,	125	,	92	,	0	,	156	,	66	,	1	,	102	,	27	},
			{	1	,	6	,	226	,	31	,	88	,	115	,	84	,	55	,	185	,	96	},
			{	1	,	7	,	156	,	187	,	0	,	200	,	98	,	37	,	17	,	23	},
			{	1	,	8	,	224	,	185	,	0	,	29	,	69	,	171	,	14	,	9	},
			{	1	,	9	,	252	,	3	,	55	,	31	,	50	,	133	,	180	,	167	},
			{	1	,	11	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	1	,	12	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	2	,	0	,	81	,	25	,	20	,	152	,	95	,	98	,	126	,	74	},
			{	2	,	1	,	114	,	114	,	94	,	131	,	106	,	168	,	163	,	31	},
			{	2	,	3	,	44	,	117	,	99	,	46	,	92	,	107	,	47	,	3	},
			{	2	,	4	,	52	,	110	,	9	,	191	,	110	,	82	,	183	,	53	},
			{	2	,	8	,	240	,	114	,	108	,	91	,	111	,	142	,	132	,	155	},
			{	2	,	10	,	1	,	1	,	1	,	0	,	1	,	1	,	1	,	0	},
			{	2	,	12	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	2	,	13	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	3	,	1	,	8	,	136	,	38	,	185	,	120	,	53	,	36	,	239	},
			{	3	,	2	,	58	,	175	,	15	,	6	,	121	,	174	,	48	,	171	},
			{	3	,	4	,	158	,	113	,	102	,	36	,	22	,	174	,	18	,	95	},
			{	3	,	5	,	104	,	72	,	146	,	124	,	4	,	127	,	111	,	110	},
			{	3	,	6	,	209	,	123	,	12	,	124	,	73	,	17	,	203	,	159	},
			{	3	,	7	,	54	,	118	,	57	,	110	,	49	,	89	,	3	,	199	},
			{	3	,	8	,	18	,	28	,	53	,	156	,	128	,	17	,	191	,	43	},
			{	3	,	9	,	128	,	186	,	46	,	133	,	79	,	105	,	160	,	75	},
			{	3	,	10	,	0	,	0	,	0	,	1	,	0	,	0	,	0	,	1	},
			{	3	,	13	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	4	,	0	,	179	,	72	,	0	,	200	,	42	,	86	,	43	,	29	},
			{	4	,	1	,	214	,	74	,	136	,	16	,	24	,	67	,	27	,	140	},
			{	4	,	11	,	71	,	29	,	157	,	101	,	51	,	83	,	117	,	180	},
			{	4	,	14	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	5	,	0	,	231	,	10	,	0	,	185	,	40	,	79	,	136	,	121	},
			{	5	,	1	,	41	,	44	,	131	,	138	,	140	,	84	,	49	,	41	},
			{	5	,	5	,	194	,	121	,	142	,	170	,	84	,	35	,	36	,	169	},
			{	5	,	7	,	159	,	80	,	141	,	219	,	137	,	103	,	132	,	88	},
			{	5	,	11	,	103	,	48	,	64	,	193	,	71	,	60	,	62	,	207	},
			{	5	,	15	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	6	,	0	,	155	,	129	,	0	,	123	,	109	,	47	,	7	,	137	},
			{	6	,	5	,	228	,	92	,	124	,	55	,	87	,	154	,	34	,	72	},
			{	6	,	7	,	45	,	100	,	99	,	31	,	107	,	10	,	198	,	172	},
			{	6	,	9	,	28	,	49	,	45	,	222	,	133	,	155	,	168	,	124	},
			{	6	,	11	,	158	,	184	,	148	,	209	,	139	,	29	,	12	,	56	},
			{	6	,	16	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	7	,	1	,	129	,	80	,	0	,	103	,	97	,	48	,	163	,	86	},
			{	7	,	5	,	147	,	186	,	45	,	13	,	135	,	125	,	78	,	186	},
			{	7	,	7	,	140	,	16	,	148	,	105	,	35	,	24	,	143	,	87	},
			{	7	,	11	,	3	,	102	,	96	,	150	,	108	,	47	,	107	,	172	},
			{	7	,	13	,	116	,	143	,	78	,	181	,	65	,	55	,	58	,	154	},
			{	7	,	17	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	8	,	0	,	142	,	118	,	0	,	147	,	70	,	53	,	101	,	176	},
			{	8	,	1	,	94	,	70	,	65	,	43	,	69	,	31	,	177	,	169	},
			{	8	,	12	,	230	,	152	,	87	,	152	,	88	,	161	,	22	,	225	},
			{	8	,	18	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	9	,	1	,	203	,	28	,	0	,	2	,	97	,	104	,	186	,	167	},
			{	9	,	8	,	205	,	132	,	97	,	30	,	40	,	142	,	27	,	238	},
			{	9	,	10	,	61	,	185	,	51	,	184	,	24	,	99	,	205	,	48	},
			{	9	,	11	,	247	,	178	,	85	,	83	,	49	,	64	,	81	,	68	},
			{	9	,	19	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	10	,	0	,	11	,	59	,	0	,	174	,	46	,	111	,	125	,	38	},
			{	10	,	1	,	185	,	104	,	17	,	150	,	41	,	25	,	60	,	217	},
			{	10	,	6	,	0	,	22	,	156	,	8	,	101	,	174	,	177	,	208	},
			{	10	,	7	,	117	,	52	,	20	,	56	,	96	,	23	,	51	,	232	},
			{	10	,	20	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	11	,	0	,	11	,	32	,	0	,	99	,	28	,	91	,	39	,	178	},
			{	11	,	7	,	236	,	92	,	7	,	138	,	30	,	175	,	29	,	214	},
			{	11	,	9	,	210	,	174	,	4	,	110	,	116	,	24	,	35	,	168	},
			{	11	,	13	,	56	,	154	,	2	,	99	,	64	,	141	,	8	,	51	},
			{	11	,	21	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	12	,	1	,	63	,	39	,	0	,	46	,	33	,	122	,	18	,	124	},
			{	12	,	3	,	111	,	93	,	113	,	217	,	122	,	11	,	155	,	122	},
			{	12	,	11	,	14	,	11	,	48	,	109	,	131	,	4	,	49	,	72	},
			{	12	,	22	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	13	,	0	,	83	,	49	,	0	,	37	,	76	,	29	,	32	,	48	},
			{	13	,	1	,	2	,	125	,	112	,	113	,	37	,	91	,	53	,	57	},
			{	13	,	8	,	38	,	35	,	102	,	143	,	62	,	27	,	95	,	167	},
			{	13	,	13	,	222	,	166	,	26	,	140	,	47	,	127	,	186	,	219	},
			{	13	,	23	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	14	,	1	,	115	,	19	,	0	,	36	,	143	,	11	,	91	,	82	},
			{	14	,	6	,	145	,	118	,	138	,	95	,	51	,	145	,	20	,	232	},
			{	14	,	11	,	3	,	21	,	57	,	40	,	130	,	8	,	52	,	204	},
			{	14	,	13	,	232	,	163	,	27	,	116	,	97	,	166	,	109	,	162	},
			{	14	,	24	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	15	,	0	,	51	,	68	,	0	,	116	,	139	,	137	,	174	,	38	},
			{	15	,	10	,	175	,	63	,	73	,	200	,	96	,	103	,	108	,	217	},
			{	15	,	11	,	213	,	81	,	99	,	110	,	128	,	40	,	102	,	157	},
			{	15	,	25	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	16	,	1	,	203	,	87	,	0	,	75	,	48	,	78	,	125	,	170	},
			{	16	,	9	,	142	,	177	,	79	,	158	,	9	,	158	,	31	,	23	},
			{	16	,	11	,	8	,	135	,	111	,	134	,	28	,	17	,	54	,	175	},
			{	16	,	12	,	242	,	64	,	143	,	97	,	8	,	165	,	176	,	202	},
			{	16	,	26	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	17	,	1	,	254	,	158	,	0	,	48	,	120	,	134	,	57	,	196	},
			{	17	,	5	,	124	,	23	,	24	,	132	,	43	,	23	,	201	,	173	},
			{	17	,	11	,	114	,	9	,	109	,	206	,	65	,	62	,	142	,	195	},
			{	17	,	12	,	64	,	6	,	18	,	2	,	42	,	163	,	35	,	218	},
			{	17	,	27	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	18	,	0	,	220	,	186	,	0	,	68	,	17	,	173	,	129	,	128	},
			{	18	,	6	,	194	,	6	,	18	,	16	,	106	,	31	,	203	,	211	},
			{	18	,	7	,	50	,	46	,	86	,	156	,	142	,	22	,	140	,	210	},
			{	18	,	28	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	19	,	0	,	87	,	58	,	0	,	35	,	79	,	13	,	110	,	39	},
			{	19	,	1	,	20	,	42	,	158	,	138	,	28	,	135	,	124	,	84	},
			{	19	,	10	,	185	,	156	,	154	,	86	,	41	,	145	,	52	,	88	},
			{	19	,	29	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	20	,	1	,	26	,	76	,	0	,	6	,	2	,	128	,	196	,	117	},
			{	20	,	4	,	105	,	61	,	148	,	20	,	103	,	52	,	35	,	227	},
			{	20	,	11	,	29	,	153	,	104	,	141	,	78	,	173	,	114	,	6	},
			{	20	,	30	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	21	,	0	,	76	,	157	,	0	,	80	,	91	,	156	,	10	,	238	},
			{	21	,	8	,	42	,	175	,	17	,	43	,	75	,	166	,	122	,	13	},
			{	21	,	13	,	210	,	67	,	33	,	81	,	81	,	40	,	23	,	11	},
			{	21	,	31	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	22	,	1	,	222	,	20	,	0	,	49	,	54	,	18	,	202	,	195	},
			{	22	,	2	,	63	,	52	,	4	,	1	,	132	,	163	,	126	,	44	},
			{	22	,	32	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	23	,	0	,	23	,	106	,	0	,	156	,	68	,	110	,	52	,	5	},
			{	23	,	3	,	235	,	86	,	75	,	54	,	115	,	132	,	170	,	94	},
			{	23	,	5	,	238	,	95	,	158	,	134	,	56	,	150	,	13	,	111	},
			{	23	,	33	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	24	,	1	,	46	,	182	,	0	,	153	,	30	,	113	,	113	,	81	},
			{	24	,	2	,	139	,	153	,	69	,	88	,	42	,	108	,	161	,	19	},
			{	24	,	9	,	8	,	64	,	87	,	63	,	101	,	61	,	88	,	130	},
			{	24	,	34	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	25	,	0	,	228	,	45	,	0	,	211	,	128	,	72	,	197	,	66	},
			{	25	,	5	,	156	,	21	,	65	,	94	,	63	,	136	,	194	,	95	},
			{	25	,	35	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	26	,	2	,	29	,	67	,	0	,	90	,	142	,	36	,	164	,	146	},
			{	26	,	7	,	143	,	137	,	100	,	6	,	28	,	38	,	172	,	66	},
			{	26	,	12	,	160	,	55	,	13	,	221	,	100	,	53	,	49	,	190	},
			{	26	,	13	,	122	,	85	,	7	,	6	,	133	,	145	,	161	,	86	},
			{	26	,	36	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	27	,	0	,	8	,	103	,	0	,	27	,	13	,	42	,	168	,	64	},
			{	27	,	6	,	151	,	50	,	32	,	118	,	10	,	104	,	193	,	181	},
			{	27	,	37	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	28	,	1	,	98	,	70	,	0	,	216	,	106	,	64	,	14	,	7	},
			{	28	,	2	,	101	,	111	,	126	,	212	,	77	,	24	,	186	,	144	},
			{	28	,	5	,	135	,	168	,	110	,	193	,	43	,	149	,	46	,	16	},
			{	28	,	38	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	29	,	0	,	18	,	110	,	0	,	108	,	133	,	139	,	50	,	25	},
			{	29	,	4	,	28	,	17	,	154	,	61	,	25	,	161	,	27	,	57	},
			{	29	,	39	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	30	,	2	,	71	,	120	,	0	,	106	,	87	,	84	,	70	,	37	},
			{	30	,	5	,	240	,	154	,	35	,	44	,	56	,	173	,	17	,	139	},
			{	30	,	7	,	9	,	52	,	51	,	185	,	104	,	93	,	50	,	221	},
			{	30	,	9	,	84	,	56	,	134	,	176	,	70	,	29	,	6	,	17	},
			{	30	,	40	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	31	,	1	,	106	,	3	,	0	,	147	,	80	,	117	,	115	,	201	},
			{	31	,	13	,	1	,	170	,	20	,	182	,	139	,	148	,	189	,	46	},
			{	31	,	41	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	32	,	0	,	242	,	84	,	0	,	108	,	32	,	116	,	110	,	179	},
			{	32	,	5	,	44	,	8	,	20	,	21	,	89	,	73	,	0	,	14	},
			{	32	,	12	,	166	,	17	,	122	,	110	,	71	,	142	,	163	,	116	},
			{	32	,	42	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	33	,	2	,	132	,	165	,	0	,	71	,	135	,	105	,	163	,	46	},
			{	33	,	7	,	164	,	179	,	88	,	12	,	6	,	137	,	173	,	2	},
			{	33	,	10	,	235	,	124	,	13	,	109	,	2	,	29	,	179	,	106	},
			{	33	,	43	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	34	,	0	,	147	,	173	,	0	,	29	,	37	,	11	,	197	,	184	},
			{	34	,	12	,	85	,	177	,	19	,	201	,	25	,	41	,	191	,	135	},
			{	34	,	13	,	36	,	12	,	78	,	69	,	114	,	162	,	193	,	141	},
			{	34	,	44	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	35	,	1	,	57	,	77	,	0	,	91	,	60	,	126	,	157	,	85	},
			{	35	,	5	,	40	,	184	,	157	,	165	,	137	,	152	,	167	,	225	},
			{	35	,	11	,	63	,	18	,	6	,	55	,	93	,	172	,	181	,	175	},
			{	35	,	45	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	36	,	0	,	140	,	25	,	0	,	1	,	121	,	73	,	197	,	178	},
			{	36	,	2	,	38	,	151	,	63	,	175	,	129	,	154	,	167	,	112	},
			{	36	,	7	,	154	,	170	,	82	,	83	,	26	,	129	,	179	,	106	},
			{	36	,	46	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	37	,	10	,	219	,	37	,	0	,	40	,	97	,	167	,	181	,	154	},
			{	37	,	13	,	151	,	31	,	144	,	12	,	56	,	38	,	193	,	114	},
			{	37	,	47	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	38	,	1	,	31	,	84	,	0	,	37	,	1	,	112	,	157	,	42	},
			{	38	,	5	,	66	,	151	,	93	,	97	,	70	,	7	,	173	,	41	},
			{	38	,	11	,	38	,	190	,	19	,	46	,	1	,	19	,	191	,	105	},
			{	38	,	48	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	39	,	0	,	239	,	93	,	0	,	106	,	119	,	109	,	181	,	167	},
			{	39	,	7	,	172	,	132	,	24	,	181	,	32	,	6	,	157	,	45	},
			{	39	,	12	,	34	,	57	,	138	,	154	,	142	,	105	,	173	,	189	},
			{	39	,	49	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	40	,	2	,	0	,	103	,	0	,	98	,	6	,	160	,	193	,	78	},
			{	40	,	10	,	75	,	107	,	36	,	35	,	73	,	156	,	163	,	67	},
			{	40	,	13	,	120	,	163	,	143	,	36	,	102	,	82	,	179	,	180	},
			{	40	,	50	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	},
			{	41	,	1	,	129	,	147	,	0	,	120	,	48	,	132	,	191	,	53	},
			{	41	,	5	,	229	,	7	,	2	,	101	,	47	,	6	,	197	,	215	},
			{	41	,	11	,	118	,	60	,	55	,	81	,	19	,	8	,	167	,	230	},
			{	41	,	51	,	0	,	0	,	0	,	0	,	0	,	0	,	0	,	0	}};


	if (BG_index == 1){
		for (int rr=0; rr<316; rr++){
			for (int cc=0; cc<2; cc++)
				BaseGraph_info_mtx[rr][cc] = BG_table_1[rr][cc] ;

		for (int rr=0; rr<316; rr++)
			BaseGraph_info_mtx[rr][2] = BG_table_1[rr][i_LS+1] % Z_c ;
		}

	}
	else if (BG_index == 2){
		for (int rr=0; rr<197; rr++){
			for (int cc=0; cc<2; cc++)
				BaseGraph_info_mtx[rr][cc] = BG_table_2[rr][cc] ;

		for (int rr=0; rr<197; rr++)
			BaseGraph_info_mtx[rr][2] = BG_table_2[rr][i_LS+1] % Z_c ;
		}	
	}

	return BaseGraph_info_mtx;
}

