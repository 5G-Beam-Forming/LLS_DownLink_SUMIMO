/*
 * module_Channel_5G.h
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *		Last Updated: 2018-10-05
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */
#ifndef __MODULE_CHANNEL_5G_H__
#define __MODULE_CHANNEL_5G_H__

#include <armadillo>
#include "sumimo.h"
#include "module_Parameter_MIMO.h"



class ModuleChannel5G {
public:
	colvec			PDP_Delay;
	colvec			PDP_Power;
	double			New_Freq;
	ivec			New_PDP_Delay_samples;
	int				Faders_num;
	int				Fast_Fading_Chan_num;
	ChannelOutput	Output;
	field<cx_mat>	Output_time_tap;

	mat				Cov_Matrix_time;
	cx_mat			Cov_Matrix_freq;

	/* various variables */
	int				Subframe_Duration;
	int				FFT_pts_num;
	int				Subcarrier_Total;
	irowvec			CP_pts_num;


	/* class method defintions */
	void module_Channel_5G(ModuleParameterMIMO &);
	void Chan_Matrix(ModuleParameterMIMO &, int);
	void Chan_H_fft(ModuleParameterMIMO &);

};   /* end of moduleChannel5G class definition */



#endif   /* end of module_Channel_5G.h */