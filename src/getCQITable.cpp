/*
 * loadCQITable.cpp - load CQI_Table from the CQI_Table_new.mat file.
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *		Last Updated: 2018-10-01
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */
#include <cstdio>
#include <cstring>
#include <cstdarg>
#include <cstdlib>
#include <complex>
#include <armadillo>
#include "getCQITable.h"


using namespace std;
using namespace arma;


CQI_Table_n			CQI_Table_1;
CQI_Table_n			CQI_Table_2;
field<cx_rowvec>	modulation_mapping(16);
field<imat>			modulation_mapping_bittable(16);


/*
 * getCQITable() - return CQI Table value
 */
void getCQITable() {

	uword numElem = 15;

	CQI_Table_1.CQI        << 1.0 << 2.0 << 3.0 << 4.0 << 5.0 << 6.0 << 7.0 << 8.0 << 9.0 << 10.0 << 11.0 << 12.0 << 13.0 << 14.0 << 15.0 << endr;
	CQI_Table_1.order      << 2.0 << 2.0 << 2.0 << 2.0 << 2.0 << 2.0 << 4.0 << 4.0 << 4.0 << 6.0 << 6.0 << 6.0 << 6.0 << 6.0 << 6.0 << endr;
	CQI_Table_1.code_rate  << 78.0 << 120.0 << 193.0 << 308.0 << 449.0 << 602.0 << 378.0 << 490.0 << 616.0 << 466.0 << 567.0 << 666.0 << 772.0 << 873.0 << 948.0 << endr;
	CQI_Table_1.efficiency << 0.1523 << 0.2344 << 0.3770 << 0.6016 << 0.8770 << 1.1758 << 1.4766 << 1.9141 << 2.4063 << 2.7305 << 3.3223 << 3.9023 << 4.5234 << 5.1152 << 5.5547 << endr;

	CQI_Table_2.CQI        << 1.0 << 2.0 << 3.0 << 4.0 << 5.0 << 6.0 << 7.0 << 8.0 << 9.0 << 10.0 << 11.0 << 12.0 << 13.0 << 14.0 << 15.0 << endr;
	CQI_Table_2.order      << 0.0 << 2.0 << 2.0 << 2.0 << 4.0 << 4.0 << 4.0 << 6.0 << 6.0 << 6.0 << 6.0 << 6.0 << 8.0 << 8.0 << 8.0 << 8.0 << endr;
	CQI_Table_2.code_rate  << 78.0 << 193.0 << 449.0 << 378.0 << 490.0 << 616.0 << 466.0 << 567.0 << 666.0 << 772.0 << 873.0 << 711.0 << 797.0 << 885.0 << 948.0 << endr;
	CQI_Table_2.efficiency << 0.1523 << 0.3770 << 0.8770 << 1.4766 << 1.9141 << 2.4063 << 2.7305 << 3.3223 << 3.9023 << 4.5234 << 5.1152 << 5.5547 << 6.2266 << 6.9141 << 7.4063 << endr;

	modulation_mapping(0) << cx_double(0.707106781186548, 0.707106781186548) << cx_double(-0.707106781186548, -0.707106781186548) << endr;
	modulation_mapping(1) << cx_double(0.707106781186548, 0.707106781186548) << cx_double(0.707106781186548, -0.707106781186548) << cx_double(-0.707106781186548, 0.707106781186548) << cx_double(-0.707106781186548, -0.707106781186548) << endr;
	modulation_mapping(2) << NULL << endr;
	modulation_mapping(3) << cx_double(0.316227766016838, 0.316227766016838) << cx_double(0.316227766016838, 0.948683298050514) << cx_double(0.948683298050514, 0.316227766016838) << cx_double(0.948683298050514, 0.948683298050514)
		<< cx_double(0.316227766016838, -0.316227766016838) << cx_double(0.316227766016838, -0.948683298050514) << cx_double(0.948683298050514, -0.316227766016838) << cx_double(0.948683298050514, -0.948683298050514)
		<< cx_double(-0.316227766016838, 0.316227766016838) << cx_double(-0.316227766016838, 0.948683298050514) << cx_double(-0.948683298050514, 0.316227766016838) << cx_double(-0.948683298050514, 0.948683298050514)
		<< cx_double(-0.316227766016838, -0.316227766016838) << cx_double(-0.316227766016838, -0.948683298050514) << cx_double(-0.948683298050514, -0.316227766016838) << cx_double(-0.948683298050514, -0.948683298050514) << endr;
	modulation_mapping(4) << NULL << endr;
	modulation_mapping(5) << cx_double(0.462910049886276, 0.462910049886276) << cx_double(0.462910049886276, 0.154303349962092) << cx_double(0.154303349962092, 0.462910049886276) << cx_double(0.154303349962092, 0.154303349962092)
		<< cx_double(0.462910049886276, 0.771516749810460) << cx_double(0.462910049886276, 1.08012344973464) << cx_double(0.154303349962092, 0.771516749810460) << cx_double(0.154303349962092, 1.08012344973464)
		<< cx_double(0.771516749810460, 0.462910049886276) << cx_double(0.771516749810460, 0.154303349962092) << cx_double(1.08012344973464, 0.462910049886276) << cx_double(1.08012344973464, 0.154303349962092)
		<< cx_double(0.771516749810460, 0.771516749810460) << cx_double(0.771516749810460, 1.08012344973464) << cx_double(1.08012344973464, 0.771516749810460) << cx_double(1.08012344973464, 1.08012344973464)
		<< cx_double(0.462910049886276, -0.462910049886276) << cx_double(0.462910049886276, -0.154303349962092) << cx_double(0.154303349962092, -0.462910049886276) << cx_double(0.154303349962092, -0.154303349962092)
		<< cx_double(0.462910049886276, -0.771516749810460) << cx_double(0.462910049886276, -1.08012344973464) << cx_double(0.154303349962092, -0.771516749810460) << cx_double(0.154303349962092, -1.08012344973464)
		<< cx_double(0.771516749810460, -0.462910049886276) << cx_double(0.771516749810460, -0.154303349962092) << cx_double(1.08012344973464, -0.462910049886276) << cx_double(1.08012344973464, -0.154303349962092)
		<< cx_double(0.771516749810460, -0.771516749810460) << cx_double(0.771516749810460, -1.08012344973464) << cx_double(1.08012344973464, -0.771516749810460) << cx_double(1.08012344973464, -1.08012344973464)
		<< cx_double(-0.462910049886276, 0.462910049886276) << cx_double(-0.462910049886276, 0.154303349962092) << cx_double(-0.154303349962092, 0.462910049886276) << cx_double(-0.154303349962092, 0.154303349962092)
		<< cx_double(-0.462910049886276, 0.771516749810460) << cx_double(-0.462910049886276, 1.08012344973464) << cx_double(-0.154303349962092, 0.771516749810460) << cx_double(-0.154303349962092, 1.08012344973464)
		<< cx_double(-0.771516749810460, 0.462910049886276) << cx_double(-0.771516749810460, 0.154303349962092) << cx_double(-1.08012344973464, 0.462910049886276) << cx_double(-1.08012344973464, 0.154303349962092)
		<< cx_double(-0.771516749810460, 0.771516749810460) << cx_double(-0.771516749810460, 1.08012344973464) << cx_double(-1.08012344973464, 0.771516749810460) << cx_double(-1.08012344973464, 1.08012344973464)
		<< cx_double(-0.462910049886276, -0.462910049886276) << cx_double(-0.462910049886276, -0.154303349962092) << cx_double(-0.154303349962092, -0.462910049886276) << cx_double(-0.154303349962092, -0.154303349962092)
		<< cx_double(-0.462910049886276, -0.771516749810460) << cx_double(-0.462910049886276, -1.08012344973464) << cx_double(-0.154303349962092, -0.771516749810460) << cx_double(-0.154303349962092, -1.08012344973464)
		<< cx_double(-0.771516749810460, -0.462910049886276) << cx_double(-0.771516749810460, -0.154303349962092) << cx_double(-1.08012344973464, -0.462910049886276) << cx_double(-1.08012344973464, -0.154303349962092)
		<< cx_double(-0.771516749810460, -0.771516749810460) << cx_double(-0.771516749810460, -1.08012344973464) << cx_double(-1.08012344973464, -0.771516749810460) << cx_double(-1.08012344973464, -1.08012344973464) << endr;
	modulation_mapping(6) << NULL << endr;
	modulation_mapping(7) << cx_double(0.383482494423685, 0.383482494423685) << cx_double(0.383482494423685, 0.536875492193159) << cx_double(0.536875492193159, 0.383482494423685) << cx_double(0.536875492193159, 0.536875492193159)
		<< cx_double(0.383482494423685, 0.230089496654211) << cx_double(0.383482494423685, 0.0766964988847370) << cx_double(0.536875492193159, 0.230089496654211) << cx_double(0.536875492193159, 0.0766964988847370)
		<< cx_double(0.230089496654211, 0.383482494423685) << cx_double(0.230089496654211, 0.536875492193159) << cx_double(0.0766964988847370, 0.383482494423685) << cx_double(0.0766964988847370, 0.536875492193159)
		<< cx_double(0.230089496654211, 0.230089496654211) << cx_double(0.230089496654211, 0.0766964988847370) << cx_double(0.0766964988847370, 0.230089496654211) << cx_double(0.0766964988847370, 0.0766964988847370)
		<< cx_double(0.383482494423685, 0.843661487732108) << cx_double(0.383482494423685, 0.690268489962633) << cx_double(0.536875492193159, 0.843661487732108) << cx_double(0.536875492193159, 0.690268489962633)
		<< cx_double(0.383482494423685, 0.997054485501582) << cx_double(0.383482494423685, 1.15044748327106) << cx_double(0.536875492193159, 0.997054485501582) << cx_double(0.536875492193159, 1.15044748327106)
		<< cx_double(0.230089496654211, 0.843661487732108) << cx_double(0.230089496654211, 0.690268489962633) << cx_double(0.0766964988847370, 0.843661487732108) << cx_double(0.0766964988847370, 0.690268489962633)
		<< cx_double(0.230089496654211, 0.997054485501582) << cx_double(0.230089496654211, 1.15044748327106) << cx_double(0.0766964988847370, 0.997054485501582) << cx_double(0.0766964988847370, 1.15044748327106)
		<< cx_double(0.843661487732108, 0.383482494423685) << cx_double(0.843661487732108, 0.536875492193159) << cx_double(0.690268489962633, 0.383482494423685) << cx_double(0.690268489962633, 0.536875492193159)
		<< cx_double(0.843661487732108, 0.230089496654211) << cx_double(0.843661487732108, 0.0766964988847370) << cx_double(0.690268489962633, 0.230089496654211) << cx_double(0.690268489962633, 0.0766964988847370)
		<< cx_double(0.997054485501582, 0.383482494423685) << cx_double(0.997054485501582, 0.536875492193159) << cx_double(1.15044748327106, 0.383482494423685) << cx_double(1.15044748327106, 0.536875492193159)
		<< cx_double(0.997054485501582, 0.230089496654211) << cx_double(0.997054485501582, 0.0766964988847370) << cx_double(1.15044748327106, 0.230089496654211) << cx_double(1.15044748327106, 0.0766964988847370)
		<< cx_double(0.843661487732108, 0.843661487732108) << cx_double(0.843661487732108, 0.690268489962633) << cx_double(0.690268489962633, 0.843661487732108) << cx_double(0.690268489962633, 0.690268489962633)
		<< cx_double(0.843661487732108, 0.997054485501582) << cx_double(0.843661487732108, 1.15044748327106) << cx_double(0.690268489962633, 0.997054485501582) << cx_double(0.690268489962633, 1.15044748327106)
		<< cx_double(0.997054485501582, 0.843661487732108) << cx_double(0.997054485501582, 0.690268489962633) << cx_double(1.15044748327106, 0.843661487732108) << cx_double(1.15044748327106, 0.690268489962633)
		<< cx_double(0.997054485501582, 0.997054485501582) << cx_double(0.997054485501582, 1.15044748327106) << cx_double(1.15044748327106, 0.997054485501582) << cx_double(1.15044748327106, 1.15044748327106)
		<< cx_double(0.383482494423685, 0.383482494423685) << cx_double(0.383482494423685, -0.536875492193159) << cx_double(0.536875492193159, -0.383482494423685) << cx_double(0.536875492193159, -0.536875492193159)
		<< cx_double(0.383482494423685, -0.230089496654211) << cx_double(0.383482494423685, -0.0766964988847370) << cx_double(0.536875492193159, -0.230089496654211) << cx_double(0.536875492193159, -0.0766964988847370)
		<< cx_double(0.230089496654211, -0.383482494423685) << cx_double(0.230089496654211, -0.536875492193159) << cx_double(0.0766964988847370, -0.383482494423685) << cx_double(0.0766964988847370, -0.536875492193159)
		<< cx_double(0.230089496654211, -0.230089496654211) << cx_double(0.230089496654211, -0.0766964988847370) << cx_double(0.0766964988847370, -0.230089496654211) << cx_double(0.0766964988847370, -0.0766964988847370)
		<< cx_double(0.383482494423685, -0.843661487732108) << cx_double(0.383482494423685, -0.690268489962633) << cx_double(0.536875492193159, -0.843661487732108) << cx_double(0.536875492193159, -0.690268489962633)
		<< cx_double(0.383482494423685, -0.997054485501582) << cx_double(0.383482494423685, -1.15044748327106) << cx_double(0.536875492193159, -0.997054485501582) << cx_double(0.536875492193159, -1.15044748327106)
		<< cx_double(0.230089496654211, -0.843661487732108) << cx_double(0.230089496654211, -0.690268489962633) << cx_double(0.0766964988847370, -0.843661487732108) << cx_double(0.0766964988847370, -0.690268489962633)
		<< cx_double(0.230089496654211, -0.997054485501582) << cx_double(0.230089496654211, -1.15044748327106) << cx_double(0.0766964988847370, -0.997054485501582) << cx_double(0.0766964988847370, -1.15044748327106)
		<< cx_double(0.843661487732108, -0.383482494423685) << cx_double(0.843661487732108, -0.536875492193159) << cx_double(0.690268489962633, -0.383482494423685) << cx_double(0.690268489962633, -0.536875492193159)
		<< cx_double(0.843661487732108, -0.230089496654211) << cx_double(0.843661487732108, -0.0766964988847370) << cx_double(0.690268489962633, -0.230089496654211) << cx_double(0.690268489962633, -0.0766964988847370)
		<< cx_double(0.997054485501582, -0.383482494423685) << cx_double(0.997054485501582, -0.536875492193159) << cx_double(1.15044748327106, -0.383482494423685) << cx_double(1.15044748327106, -0.536875492193159)
		<< cx_double(0.997054485501582, -0.230089496654211) << cx_double(0.997054485501582, -0.0766964988847370) << cx_double(1.15044748327106, -0.230089496654211) << cx_double(1.15044748327106, -0.0766964988847370)
		<< cx_double(0.843661487732108, -0.843661487732108) << cx_double(0.843661487732108, -0.690268489962633) << cx_double(0.690268489962633, -0.843661487732108) << cx_double(0.690268489962633, -0.690268489962633)
		<< cx_double(0.843661487732108, -0.997054485501582) << cx_double(0.843661487732108, -1.15044748327106) << cx_double(0.690268489962633, -0.997054485501582) << cx_double(0.690268489962633, -1.15044748327106)
		<< cx_double(0.997054485501582, -0.843661487732108) << cx_double(0.997054485501582, -0.690268489962633) << cx_double(1.15044748327106, -0.843661487732108) << cx_double(1.15044748327106, -0.690268489962633)
		<< cx_double(0.997054485501582, -0.997054485501582) << cx_double(0.997054485501582, -1.15044748327106) << cx_double(1.15044748327106, -0.997054485501582) << cx_double(1.15044748327106, -1.15044748327106)
		<< cx_double(-0.383482494423685, 0.383482494423685) << cx_double(-0.383482494423685, 0.536875492193159) << cx_double(-0.536875492193159, 0.383482494423685) << cx_double(-0.536875492193159, 0.536875492193159)
		<< cx_double(-0.383482494423685, 0.230089496654211) << cx_double(-0.383482494423685, 0.0766964988847370) << cx_double(-0.536875492193159, 0.230089496654211) << cx_double(-0.536875492193159, 0.0766964988847370)
		<< cx_double(-0.230089496654211, 0.383482494423685) << cx_double(-0.230089496654211, 0.536875492193159) << cx_double(-0.0766964988847370, 0.383482494423685) << cx_double(-0.0766964988847370, 0.536875492193159)
		<< cx_double(-0.230089496654211, 0.230089496654211) << cx_double(-0.230089496654211, 0.0766964988847370) << cx_double(-0.0766964988847370, 0.230089496654211) << cx_double(-0.0766964988847370, 0.0766964988847370)
		<< cx_double(-0.383482494423685, 0.843661487732108) << cx_double(-0.383482494423685, 0.690268489962633) << cx_double(-0.536875492193159, 0.843661487732108) << cx_double(-0.536875492193159, 0.690268489962633)
		<< cx_double(-0.383482494423685, 0.997054485501582) << cx_double(-0.383482494423685, 1.15044748327106) << cx_double(-0.536875492193159, 0.997054485501582) << cx_double(-0.536875492193159, 1.15044748327106)
		<< cx_double(-0.230089496654211, 0.843661487732108) << cx_double(-0.230089496654211, 0.690268489962633) << cx_double(-0.0766964988847370, 0.843661487732108) << cx_double(-0.0766964988847370, 0.690268489962633)
		<< cx_double(-0.230089496654211, 0.997054485501582) << cx_double(-0.230089496654211, 1.15044748327106) << cx_double(-0.0766964988847370, 0.997054485501582) << cx_double(-0.0766964988847370, 1.15044748327106)
		<< cx_double(-0.843661487732108, 0.383482494423685) << cx_double(-0.843661487732108, 0.536875492193159) << cx_double(-0.690268489962633, 0.383482494423685) << cx_double(-0.690268489962633, 0.536875492193159)
		<< cx_double(-0.843661487732108, 0.230089496654211) << cx_double(-0.843661487732108, 0.0766964988847370) << cx_double(-0.690268489962633, 0.230089496654211) << cx_double(-0.690268489962633, 0.0766964988847370)
		<< cx_double(-0.997054485501582, 0.383482494423685) << cx_double(-0.997054485501582, 0.536875492193159) << cx_double(-1.15044748327106, 0.383482494423685) << cx_double(-1.15044748327106, 0.536875492193159)
		<< cx_double(-0.997054485501582, 0.230089496654211) << cx_double(-0.997054485501582, 0.0766964988847370) << cx_double(-1.15044748327106, 0.230089496654211) << cx_double(-1.15044748327106, 0.0766964988847370)
		<< cx_double(-0.843661487732108, 0.843661487732108) << cx_double(-0.843661487732108, 0.690268489962633) << cx_double(-0.690268489962633, 0.843661487732108) << cx_double(-0.690268489962633, 0.690268489962633)
		<< cx_double(-0.843661487732108, 0.997054485501582) << cx_double(-0.843661487732108, 1.15044748327106) << cx_double(-0.690268489962633, 0.997054485501582) << cx_double(-0.690268489962633, 1.15044748327106)
		<< cx_double(-0.997054485501582, 0.843661487732108) << cx_double(-0.997054485501582, 0.690268489962633) << cx_double(-1.15044748327106, 0.843661487732108) << cx_double(-1.15044748327106, 0.690268489962633)
		<< cx_double(-0.997054485501582, 0.997054485501582) << cx_double(-0.997054485501582, 1.15044748327106) << cx_double(-1.15044748327106, 0.997054485501582) << cx_double(-1.15044748327106, 1.15044748327106)
		<< cx_double(-0.383482494423685, -0.383482494423685) << cx_double(-0.383482494423685, -0.536875492193159) << cx_double(-0.536875492193159, -0.383482494423685) << cx_double(-0.536875492193159, -0.536875492193159)
		<< cx_double(-0.383482494423685, -0.230089496654211) << cx_double(-0.383482494423685, -0.0766964988847370) << cx_double(-0.536875492193159, -0.230089496654211) << cx_double(-0.536875492193159, -0.0766964988847370)
		<< cx_double(-0.230089496654211, -0.383482494423685) << cx_double(-0.230089496654211, -0.536875492193159) << cx_double(-0.0766964988847370, -0.383482494423685) << cx_double(-0.0766964988847370, -0.536875492193159)
		<< cx_double(-0.230089496654211, -0.230089496654211) << cx_double(-0.230089496654211, -0.0766964988847370) << cx_double(-0.0766964988847370, -0.230089496654211) << cx_double(-0.0766964988847370, -0.0766964988847370)
		<< cx_double(-0.383482494423685, -0.843661487732108) << cx_double(-0.383482494423685, -0.690268489962633) << cx_double(-0.536875492193159, -0.843661487732108) << cx_double(-0.536875492193159, -0.690268489962633)
		<< cx_double(-0.383482494423685, -0.997054485501582) << cx_double(-0.383482494423685, -1.15044748327106) << cx_double(-0.536875492193159, -0.997054485501582) << cx_double(-0.536875492193159, -1.15044748327106)
		<< cx_double(-0.230089496654211, -0.843661487732108) << cx_double(-0.230089496654211, -0.690268489962633) << cx_double(-0.0766964988847370, -0.843661487732108) << cx_double(-0.0766964988847370, -0.690268489962633)
		<< cx_double(-0.230089496654211, -0.997054485501582) << cx_double(-0.230089496654211, -1.15044748327106) << cx_double(-0.0766964988847370, -0.997054485501582) << cx_double(-0.0766964988847370, -1.15044748327106)
		<< cx_double(-0.843661487732108, -0.383482494423685) << cx_double(-0.843661487732108, -0.536875492193159) << cx_double(-0.690268489962633, -0.383482494423685) << cx_double(-0.690268489962633, -0.536875492193159)
		<< cx_double(-0.843661487732108, -0.230089496654211) << cx_double(-0.843661487732108, -0.0766964988847370) << cx_double(-0.690268489962633, -0.230089496654211) << cx_double(-0.690268489962633, -0.0766964988847370)
		<< cx_double(-0.997054485501582, -0.383482494423685) << cx_double(-0.997054485501582, -0.536875492193159) << cx_double(-1.15044748327106, -0.383482494423685) << cx_double(-1.15044748327106, -0.536875492193159)
		<< cx_double(-0.997054485501582, -0.230089496654211) << cx_double(-0.997054485501582, -0.0766964988847370) << cx_double(-1.15044748327106, -0.230089496654211) << cx_double(-1.15044748327106, -0.0766964988847370)
		<< cx_double(-0.843661487732108, -0.843661487732108) << cx_double(-0.843661487732108, -0.690268489962633) << cx_double(-0.690268489962633, -0.843661487732108) << cx_double(-0.690268489962633, -0.690268489962633)
		<< cx_double(-0.843661487732108, -0.997054485501582) << cx_double(-0.843661487732108, -1.15044748327106) << cx_double(-0.690268489962633, -0.997054485501582) << cx_double(-0.690268489962633, -1.15044748327106)
		<< cx_double(-0.997054485501582, -0.843661487732108) << cx_double(-0.997054485501582, -0.690268489962633) << cx_double(-1.15044748327106, -0.843661487732108) << cx_double(-1.15044748327106, -0.690268489962633)
		<< cx_double(-0.997054485501582, -0.997054485501582) << cx_double(-0.997054485501582, -1.15044748327106) << cx_double(-1.15044748327106, -0.997054485501582) << cx_double(-1.15044748327106, -1.15044748327106) << endr;
	modulation_mapping(8) << NULL << endr;
	modulation_mapping(9) << NULL << endr;
	modulation_mapping(10) << NULL << endr;
	modulation_mapping(11) << NULL << endr;
	modulation_mapping(12) << NULL << endr;
	modulation_mapping(13) << NULL << endr;
	modulation_mapping(14) << NULL << endr;
	modulation_mapping(15) << NULL << endr;


	modulation_mapping_bittable(0) << 0 << 1 << endr;
	modulation_mapping_bittable(1) << 0 << 1 << 0 << 1 << endr << 0 << 0 << 1 << 1 << endr;
	modulation_mapping_bittable(2) << NULL << endr;
	modulation_mapping_bittable(3) << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << endr
		<< 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << endr;
	modulation_mapping_bittable(4) << NULL << endr;
	modulation_mapping_bittable(5) << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0
		<< 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << endr
		<< 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0
		<< 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0
		<< 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << endr;
	modulation_mapping_bittable(6) << NULL << endr;
	modulation_mapping_bittable(7) << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0
		<< 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0
		<< 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0
		<< 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0
		<< 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0
		<< 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0
		<< 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0
		<< 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << endr
		<< 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0
		<< 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1
		<< 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0
		<< 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1
		<< 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0
		<< 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1
		<< 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0
		<< 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << 0 << 0 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0
		<< 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0
		<< 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1
		<< 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1
		<< 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0
		<< 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0
		<< 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1
		<< 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1
		<< 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << 1 << endr;
	modulation_mapping_bittable(8) << NULL << endr;
	modulation_mapping_bittable(9) << NULL << endr;
	modulation_mapping_bittable(10) << NULL << endr;
	modulation_mapping_bittable(11) << NULL << endr;
	modulation_mapping_bittable(12) << NULL << endr;
	modulation_mapping_bittable(13) << NULL << endr;
	modulation_mapping_bittable(14) << NULL << endr;
	modulation_mapping_bittable(15) << NULL << endr;

}	/* end of getCQITable() */



/* end of getCQItable.cpp */