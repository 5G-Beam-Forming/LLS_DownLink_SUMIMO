/*
 * CQITable.h - .
 *
 *		Programmed by PJW
 *		Last Updated: 2018-05-10
 *
 * Copyright(c) 2015-2018 KAIST. All rights reserved.
 */

#ifndef __CQITABLE_H__
#define __CQITABLE_H__

#include <cstdarg>
#include <vector>
#include <complex>
#include <armadillo>


using namespace std;
using namespace arma;


struct CQI_Table_n {
	rowvec CQI;
	rowvec order;
	rowvec code_rate;
	rowvec efficiency;
};


extern "C" CQI_Table_n			CQI_Table_1;
extern "C" CQI_Table_n			CQI_Table_2;
extern  field<cx_rowvec>		modulation_mapping;
extern  field<imat>				modulation_mapping_bittable;


void getCQITable(void);



#endif   /* end of CQITable.h */