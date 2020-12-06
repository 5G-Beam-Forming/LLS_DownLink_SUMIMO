/*
 * module_Channel_5G.cpp
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *		Last Updated: 2018-05-03
 *
 * Copyright(c) 2015-2018 KAIST. All rights reserved.
 */

#include <cstdio>
#include <armadillo>
#include "module_Channel_5G.h"
#include "ChannelOutput_CDL.h"
#include "ChannelOutput_TDL.h"


using namespace arma;


/* local function declarations */
double sinc(const double);



/*
 * module_Channel_5G() - Gateway mtehod for module_Channel_5G class.
 */
void ModuleChannel5G::module_Channel_5G(ModuleParameterMIMO &para)
{

	this->Chan_Matrix(para, 1);
	this->Chan_H_fft(para);

}  /* end of module_Channel_5G() method definition */



/*
 * Chan_Matrix()
 */
void ModuleChannel5G::Chan_Matrix(ModuleParameterMIMO &para, int Subframe_index)
{
	/* local variables declarations */
	uword			i, j, k, l, m;										// index for loop and vector, matrix, cube
	cx_mat			time2D;

	icolvec			distinct_New_PDP_Delay_samples;
	icolvec			equal_New_PDP_Delay_samples;
	cx_double		equal_New_PDP_Delay_samples_sum = cx_double(0, 0);
	double			H_Norm;
	double			Number_of_Samples;
	icolvec			New_PDP_Delay_unique;

	double			temp_sum = 0.0;
	rowvec			Cau_Sinc;
	mat				Pre_Sincs;
	cx_cube			Chan_Sig;
	int				Chan_Sig_dims = 0;
	double			norm_coef;
	int				Time_temp;

	cube			Alpha;
	cube			Beta;
	cube			Gamma;
	mat				Beta_temp;
	mat				Gamma_temp;
	cube			Matrix_P;
	mat				Matrix_P_temp;
	mat				angle_temp;
	mat				Angle_of_incoming;
	mat				Rel_Part;
	mat				Img_Part;

	int				Subframe_Sym_num;
	int				Samp_Time;
	int				Time_symbol = 1;

	field<cx_mat>	Chan_Out_temp;
	field<cx_cube>	Chan_Out_temp5D;
	cx_cube			Chan_Out;
	field<cx_mat>	Chan_Out4D;

	cube			Rx_Correlation;
	cube			Tx_Correlation;

	/* code start */
	this->FFT_pts_num = para.size_fft;

	switch (para.channel.type) {
	case Ch_type::AWGN:
		this->Output.time = ones<cx_cube>(para.num_Rx_antenna, para.num_Tx_antenna, 1);
		break;
	case Ch_type::flat_Rayleigh:
		this->Output.time = zeros<cx_cube>(para.num_Rx_antenna, para.num_Tx_antenna, 1);
		time2D = zeros<cx_mat>(para.num_Rx_antenna, para.num_Tx_antenna);
		time2D = 1.0 / sqrt(2.0) * (randn<mat>(para.num_Rx_antenna, para.num_Tx_antenna) + cx_double(0, 1) * randn<mat>(para.num_Rx_antenna, para.num_Tx_antenna));
		this->Output.time.slice(0) = time2D;
		time2D.reset();		// local memory deallocation
		break;
	case Ch_type::CDL_A:
	case Ch_type::CDL_B:
	case Ch_type::CDL_C:
	case Ch_type::CDL_D:
	case Ch_type::CDL_E:
	case Ch_type::TDL_A:
	case Ch_type::TDL_B:
	case Ch_type::TDL_C:
	case Ch_type::TDL_D:
	case Ch_type::TDL_E:
		switch (para.channel.type) {
		case Ch_type::CDL_A:
		case Ch_type::CDL_B:
		case Ch_type::CDL_C:
		case Ch_type::CDL_D:
		case Ch_type::CDL_E:
			ChannelOutput_CDL(para, this->PDP_Delay, this->PDP_Power, this->Output_time_tap, this->Output.angle_first);
			break;
		case Ch_type::TDL_A:
		case Ch_type::TDL_B:
		case Ch_type::TDL_C:
		case Ch_type::TDL_D:
		case Ch_type::TDL_E:
			ChannelOutput_TDL(para, this->PDP_Delay, this->PDP_Power, this->Output_time_tap);
		}

		this->New_Freq = para.Fs;
		this->New_PDP_Delay_samples = zeros<icolvec>(this->PDP_Delay.n_elem);
		for (uword i = 0; i < this->PDP_Delay.n_elem; i++)
			this->New_PDP_Delay_samples(i) = round(this->PDP_Delay(i) * this->New_Freq);

		New_PDP_Delay_unique = unique(this->New_PDP_Delay_samples);
		this->Faders_num = New_PDP_Delay_unique.n_elem;

		switch (para.channel.time_fading) {
		case Ch_mode::Block_Fading:
			if (para.channel.time_correlation == Time_correlation::independent) {
				if (para.method_interpolation == Method_interpolation::nearest_neighbor) {
					this->Output.time = zeros<cx_cube>(para.num_Rx_antenna, para.num_Tx_antenna, this->New_PDP_Delay_samples(this->New_PDP_Delay_samples.n_elem - 1) + 1);
					distinct_New_PDP_Delay_samples = unique(this->New_PDP_Delay_samples);

					equal_New_PDP_Delay_samples = zeros<icolvec>(this->New_PDP_Delay_samples.n_elem);
					for (uword ind_t = 0; ind_t < para.t.n_elem; ind_t++) {
						for (uword Fader_Index = 0; Fader_Index < this->Faders_num; Fader_Index++) {
							for (uword i = 0; i < equal_New_PDP_Delay_samples.n_elem; i++)
								equal_New_PDP_Delay_samples(i) = distinct_New_PDP_Delay_samples(Fader_Index) == this->New_PDP_Delay_samples(i);
							for (uword i = 0; i < this->Output_time_tap.n_rows; i++)
								for (uword j = 0; j < this->Output_time_tap.n_cols; j++) {
									equal_New_PDP_Delay_samples_sum = cx_double(0, 0);
									for (uword k = 0; k < equal_New_PDP_Delay_samples.n_elem; k++)
										if (equal_New_PDP_Delay_samples(k) == 1)
											equal_New_PDP_Delay_samples_sum += this->Output_time_tap(i, j)(ind_t, k);
									this->Output.time(i, j, distinct_New_PDP_Delay_samples(Fader_Index)) = equal_New_PDP_Delay_samples_sum;
								}
						}
					}

					for (i = 0; i < this->PDP_Power.n_elem; i++)
						temp_sum += pow(10.0, this->PDP_Power(i) / 10.0);

					H_Norm = sqrt(temp_sum);
					this->Output.time = this->Output.time / H_Norm;
				}
				else if (para.method_interpolation == Method_interpolation::sinc_interpolation) {

					distinct_New_PDP_Delay_samples = unique(this->New_PDP_Delay_samples);

					Number_of_Samples = para.time_CP[0] * para.Fs;

					Cau_Sinc = zeros<rowvec>(Number_of_Samples);
					for (i = 0; i < Number_of_Samples; i++)
						Cau_Sinc(i) = (double)i / this->New_Freq;

					Pre_Sincs = zeros<mat>(this->Faders_num, Cau_Sinc.n_elem);
					Chan_Sig = zeros<cx_cube>(this->Output_time_tap.n_rows, this->Output_time_tap.n_cols, this->Output_time_tap(0, 0).n_cols);
					Chan_Sig_dims = 3;
					for (i = 0; i < Chan_Sig.n_rows; i++)
						for (j = 0; j < Chan_Sig.n_cols; j++)
							for (k = 0; k < Chan_Sig.n_slices; k++)
								Chan_Sig(i, j, k) = this->Output_time_tap(i, j)(0, k);

					for (uword Temp_Index = 0; Temp_Index < this->Faders_num; Temp_Index++) {
						for (i = 0; i < Number_of_Samples; i++)
							Pre_Sincs(Temp_Index, i) = sinc(2.0 * i - 2.0 * distinct_New_PDP_Delay_samples(Temp_Index))
							* sqrt(para.Cov_Matrix_time(distinct_New_PDP_Delay_samples(Temp_Index) + 1, distinct_New_PDP_Delay_samples(Temp_Index) + 1));
					}

					norm_coef = 0.0;
					for (i = 0; i < Pre_Sincs.n_elem; i++)
						norm_coef += pow(Pre_Sincs(i), 2.0);

					Chan_Out_temp = field<cx_mat>(para.num_Rx_antenna, para.num_Tx_antenna);
					for (i = 0; i < para.num_Rx_antenna; i++)
						for (j = 0; j < para.num_Tx_antenna; j++)
							Chan_Out_temp(i, j) = zeros<cx_mat>(this->Faders_num, Cau_Sinc.n_elem);
					Chan_Out = zeros<cx_cube>(para.num_Rx_antenna, para.num_Tx_antenna, Cau_Sinc.n_elem);

					for (uword Fader_Index = 0; Fader_Index < this->Faders_num; Fader_Index++)
						for (uword Ind_Rx = 0; Ind_Rx < para.num_Tx_antenna; Ind_Rx++)
							for (uword Ind_Tx = 0; Ind_Tx < para.num_Tx_antenna; Ind_Tx++)
								for (i = 0; i < Cau_Sinc.n_elem; i++)
									Chan_Out_temp(Ind_Rx, Ind_Tx)(Fader_Index, i) = Chan_Sig(Ind_Rx, Ind_Tx, Fader_Index) * Pre_Sincs(Fader_Index, i);

					for (i = 0; i < para.num_Tx_antenna; i++)
						for (j = 0; j < para.num_Rx_antenna; j++)
							for (k = 0; k < this->Faders_num; k++) {
								for (l = 0; l < Cau_Sinc.n_elem; l++)
									Chan_Out(i, j, k) += Chan_Out_temp(i, j)(k, l);
							}

					this->Output.time = Chan_Out / norm_coef;

				}
			}
			else if (para.channel.time_correlation == Time_correlation::correlated) {
				;	// empty statement.
			}
			break;
		case Ch_mode::Fast_Fading:
			;    // empty statment;
		}   /* end of switch (para.channel.time_fading) */
		break;

	default:
		switch (para.channel.time_fading) {
		case Ch_mode::Block_Fading:
			if (para.channel.time_correlation == Time_correlation::correlated) {

				Time_temp = (Subframe_index - 1) * this->Subframe_Duration;

				Alpha = (randu(1, 1, this->Faders_num) * 2.0 - 1.0) * datum::pi;
				Beta_temp = (randu(this->Faders_num, para.sin_num) * 2.0 - 1.0) * datum::pi;
				Gamma_temp = (randu(this->Faders_num, para.sin_num) * 2.0 - 1.0) * datum::pi;
				
				Matrix_P_temp = zeros<mat>(this->Faders_num, para.sin_num);
				for (uword i = 0; i < para.sin_num; i++)
					Matrix_P_temp(0, i) = (double)i * 2.0 * datum::pi;
				Matrix_P_temp = repmat(Matrix_P_temp(0, span::all), para.num_Tx_antenna, this->Faders_num);
				angle_temp = zeros<mat>(para.num_Tx_antenna, this->Faders_num);
				angle_temp = Alpha;
				angle_temp = repmat(angle_temp(span::all, span::all), 0, para.sin_num);
				Angle_of_incoming = (Matrix_P_temp - datum::pi + angle_temp) / (4.0 * para.sin_num);

				Rel_Part = cos(Gamma_temp) * cos(para.w_d * Time_temp * cos(Angle_of_incoming) + Beta_temp);
				Img_Part = sin(Gamma_temp) * cos(para.w_d * Time_temp * cos(Angle_of_incoming) + Beta_temp);
				Chan_Sig.slice(0) = 2.0 / sqrt(2.0 * para.sin_num) * sum(Rel_Part + cx_double(0, 1) * Img_Part, 2);

				printf("Module_channel_5G.cpp line 258: Code Incomplete!\n");
			}
			else if (para.channel.time_correlation == Time_correlation::independent) {
				Chan_Sig = (randn(para.num_Rx_antenna, para.num_Tx_antenna, this->Faders_num)
					+ cx_double(0, 1) * randn(para.num_Rx_antenna, para.num_Tx_antenna, this->Faders_num)) / sqrt(2.0);
			}
			break;
		case Ch_mode::Fast_Fading:		// coding incomplete
			Subframe_Sym_num = 2.0 * (this->CP_pts_num(0) + (para.num_symb - 1) * this->CP_pts_num(1) + this->FFT_pts_num * para.num_symb);
			Samp_Time = Time_symbol / this->FFT_pts_num;
			this->Fast_Fading_Chan_num = Subframe_Sym_num;
			
			printf("Module_channel_5G.cpp line 258: Code Incomplete!\n");
		}


		/* Interpolation Part */
		switch (para.method_interpolation) {
		case Method_interpolation::sinc_interpolation:
			Number_of_Samples = 144.0 / (15000.0 * 2048.0) * para.subcarrier_spacing * this->FFT_pts_num;

			for (i = 0; i <= Number_of_Samples; i++)
				Cau_Sinc(i) = (double)i / this->New_Freq;

			Pre_Sincs = zeros<mat>(this->Faders_num, Cau_Sinc.n_elem);

			for (uword Temp_Index = 0; Temp_Index < this->Faders_num; Temp_Index++)
				for (i = 0; i < Cau_Sinc.n_elem; i++)
					Pre_Sincs(Temp_Index, i) = sinc(2.0 * this->New_Freq * (Cau_Sinc(i) - this->PDP_Delay(Temp_Index))) * sqrt(pow(10.0, 0.1 * this->PDP_Power(Temp_Index)));

			norm_coef = 0;
			for (i = 0; i < Pre_Sincs.n_elem; i++)
				norm_coef += pow(Pre_Sincs(i), 2.0);

			if (Chan_Sig_dims == 3) {	/* Block_Fading */
				Chan_Out_temp = field<cx_mat>(para.num_Rx_antenna, para.num_Tx_antenna);
				for (i = 0; i < para.num_Rx_antenna; i++)
					for (j = 0; j < para.num_Tx_antenna; j++)
						Chan_Out_temp(i, j) = zeros<cx_mat>(this->Faders_num, Cau_Sinc.n_elem);
				Chan_Out = zeros<cx_cube>(para.num_Rx_antenna, para.num_Tx_antenna, Cau_Sinc.n_elem);

				for (uword Fader_Index = 0; Fader_Index < this->Faders_num; Fader_Index++)
					for (uword Ind_Rx = 0; Ind_Rx < Chan_Sig.n_rows; Ind_Rx++)
						for (uword Ind_Tx = 0; Ind_Tx < Chan_Sig.n_cols; Ind_Tx++)
							for (i = 0; i < Cau_Sinc.n_elem; i++)
								Chan_Out_temp(Ind_Rx, Ind_Tx)(Fader_Index, i) = Chan_Sig(Ind_Rx, Ind_Tx, Fader_Index) * Pre_Sincs(Fader_Index, i);

				for (i = 0; i < para.num_Tx_antenna; i++)
					for (j = 0; j < para.num_Rx_antenna; j++)
						for (k = 0; k < this->Faders_num; k++) {
							Chan_Out(i, j, k) = cx_double(0, 0);
							for (l = 0; l < Cau_Sinc.n_elem; l++)
								Chan_Out(i, j, k) += Chan_Out_temp(i, j)(k, l);
						}

				this->Output.time = Chan_Out / norm_coef;
			}
			else if (Chan_Sig_dims == 4) {	/* Fast Fading */
				Chan_Out_temp5D = field<cx_cube>(para.num_Rx_antenna, para.num_Tx_antenna);
				for (i = 0; i < para.num_Rx_antenna; i++)
					for (j = 0; j < para.num_Tx_antenna; j++)
						Chan_Out_temp5D(i, j) = zeros<cx_cube>(this->Fast_Fading_Chan_num, this->Faders_num, Cau_Sinc.n_elem);

				Chan_Out4D = field<cx_mat>(para.num_Rx_antenna, para.num_Tx_antenna);
				for (i = 0; i < para.num_Rx_antenna; i++)
					for (j = 0; j < para.num_Tx_antenna; j++)
						Chan_Out4D(i, j) = zeros<cx_mat>(this->Fast_Fading_Chan_num, Cau_Sinc.n_elem);

				for (uword Fader_Index = 0; Fader_Index < this->Faders_num; Fader_Index++)
					for (uword Ind_Rx = 0; Ind_Rx < Chan_Sig.n_cols; Ind_Rx++)
						for (uword Ind_Tx = 0; Ind_Tx < Chan_Sig.n_cols; Ind_Tx++)
							for (i = 0; i < this->Fast_Fading_Chan_num; i++)
								for (j = 0; j < Cau_Sinc.n_elem; j++)
									Chan_Out_temp5D(Ind_Rx, Ind_Tx)(i, Fader_Index, j) = Chan_Sig(Ind_Rx, Ind_Tx, i) * Pre_Sincs(Fader_Index, j);


				for (i = 0; i < para.num_Tx_antenna; i++)
					for (j = 0; j < para.num_Rx_antenna; j++)
						for (k = 0; k < this->Fast_Fading_Chan_num; k++)
							for (m = 0; m < this->Faders_num; m++) {
								Chan_Out4D(i, j)(k, m) = cx_double(0, 0);
								for (l = 0; l < Cau_Sinc.n_elem; l++)
									Chan_Out4D(i, j)(k, l) += Chan_Out_temp5D(i, j)(k, l, m);
							}

				for (uword i = 0; i < this->Output.time4D.n_elem; i++)
					this->Output.time4D(i) = Chan_Out4D(i) / norm_coef;
			}

			if (para.channel.type == Ch_type::Rayleigh2 && para.channel.time_fading == Ch_mode::Fast_Fading)
				this->Output.time = Chan_Sig;

			break;
		case Method_interpolation::nearest_neighbor:
			if (Chan_Sig_dims == 3) {	/* Block_Fading */

				switch (para.channel.type) {
				case Ch_type::PedA:
				case Ch_type::PedB:
				case Ch_type::VehA:
				case Ch_type::VehB:
				case Ch_type::Rayleigh2:
				case Ch_type::TDL_A:
				case Ch_type::TDL_B:
				case Ch_type::TDL_C:
				case Ch_type::TDL_D:
				case Ch_type::TDL_E:
					Rx_Correlation = ones<cube>(para.channel.PDP_dB.n_cols, para.num_Rx_antenna, para.num_Tx_antenna);
					Tx_Correlation = ones<cube>(para.channel.PDP_dB.n_cols, para.num_Rx_antenna, para.num_Tx_antenna);
					for (uword Index_corr_temp = 0; Index_corr_temp < para.channel.PDP_dB.n_cols; Index_corr_temp++) {
						Rx_Correlation.slice(Index_corr_temp) = eye(para.num_Rx_antenna, para.num_Rx_antenna);
						Tx_Correlation.slice(Index_corr_temp) = eye(para.num_Tx_antenna, para.num_Tx_antenna);
					}
					break;
				case Ch_type::PedBcorr:
				case Ch_type::EVA5Hz:
				case Ch_type::EVA70Hz:
				case Ch_type::ETU70Hz:
				case Ch_type::ETU300Hz:
					Rx_Correlation = ones(para.channel.PDP_dB.n_cols, para.num_Rx_antenna, para.num_Tx_antenna);
					Tx_Correlation = ones(para.channel.PDP_dB.n_cols, para.num_Rx_antenna, para.num_Tx_antenna);
					for (uword Index_corr_temp = 0; Index_corr_temp < para.channel.PDP_dB.n_cols; Index_corr_temp++) {
						Rx_Correlation.slice(Index_corr_temp) = eye(para.num_Rx_antenna, para.num_Rx_antenna) + para.Channel_coefRX * ones(para.num_Rx_antenna, para.num_Rx_antenna) - para.Channel_coefRX * eye(para.num_Rx_antenna, para.num_Rx_antenna);
						Tx_Correlation.slice(Index_corr_temp) = eye(para.num_Tx_antenna, para.num_Tx_antenna) + para.Channel_coefTx * ones(para.num_Tx_antenna, para.num_Tx_antenna) - para.Channel_coefTx * eye(para.num_Tx_antenna, para.num_Tx_antenna);
					}
					break;
				}

				this->Output.time = zeros<cx_cube>(para.num_Rx_antenna, para.num_Tx_antenna, this->New_PDP_Delay_samples(this->New_PDP_Delay_samples.n_elem - 1) + 1);
				distinct_New_PDP_Delay_samples = unique(this->New_PDP_Delay_samples);

				for (uword Fader_Index = 0; Fader_Index < this->Faders_num; Fader_Index++) {
					for (i = 0; i < this->New_PDP_Delay_samples.n_elem; i++)
						equal_New_PDP_Delay_samples(i) = distinct_New_PDP_Delay_samples(Fader_Index) == this->New_PDP_Delay_samples(i);

					// the following 4 lines are Not Impletmented yet.
					//this->Output.time.slice(distinct_New_PDP_Delay_samples(Fader_Index))
					//	= sqrt(sum(pow(10.0, (this->PDP_Power(find(equal_New_PDP_Delay_samples)) / 10.0) ) ) )
					//	* (sqrtmat(squeeze(Rx_Correlation(Fader_Index, :, : ))) * Chan_Sig(;,:,Fader_Index)
					//		* (sqrtmat(squeeze(Tx_Correlation(Fader_Index, :, : )))).t());
				}

			}
			else if (Chan_Sig_dims == 4) { /* Fast_Fading */
				this->Output.time4D = field<cx_mat>(para.num_Rx_antenna, para.num_Tx_antenna);
				for (i = 0; i < para.num_Rx_antenna; i++)
					for (j = 0; j < para.num_Tx_antenna; j++)
						this->Output.time4D(i, j) = zeros<cx_mat>(this->Fast_Fading_Chan_num, this->New_PDP_Delay_samples(this->New_PDP_Delay_samples.n_elem) + 1);

				distinct_New_PDP_Delay_samples = unique(this->New_PDP_Delay_samples);

				for (uword Fader_Index = 0; Fader_Index < this->Faders_num; Fader_Index++) {
					for (i = 0; i < this->New_PDP_Delay_samples.n_elem; i++)
						equal_New_PDP_Delay_samples(i) = distinct_New_PDP_Delay_samples(Fader_Index) == this->New_PDP_Delay_samples(i);

					// the following 4 lines are Not Impletmented yet.
					//this->Output.time3D.slice(distinct_New_PDP_Delay_samples(Fader_Index))
					//	= sqrt(sum(pow(10, (this->PDP_Power(equal_New_PDP_Delay_samples) / 10))))
					//	* (sqrtmat(squeeze(Rx_Correlation(Fader_Index, :, : ))) * Chan_Sig(;,:,Fader_Index)
					//		* (sqrtmat(squeeze(Tx_Correlation(Fader_Index, :, : )))).t());
				}
			}

		}   /* end of switch(para.channel.type) */

	}   /* end of switch (para.method_interpolation) */

	/* memory deallocations */
	distinct_New_PDP_Delay_samples.reset();
	equal_New_PDP_Delay_samples.reset();
	New_PDP_Delay_unique.reset();
	Cau_Sinc.reset();
	Pre_Sincs.reset();
	for (uword i = 0; i < Chan_Out_temp.n_elem; i++)
		Chan_Out_temp(i).reset();
	Chan_Out_temp.reset();
	for (uword i = 0; i < Chan_Out_temp5D.n_elem; i++)
		Chan_Out_temp5D(i).reset();
	Chan_Out_temp5D.reset();
	Chan_Out.reset();
	for (uword i = 0; i < Chan_Out4D.n_elem; i++)
		Chan_Out4D(i).reset();
	Chan_Out4D.reset();
	Rx_Correlation.reset();
	Tx_Correlation.reset();
	Chan_Sig.reset();

}   /* end of Chan_Matrix() method definition */





/*
 * Chan_H_fft()
 */
void ModuleChannel5G::Chan_H_fft(ModuleParameterMIMO &para) {

	cx_colvec		temp_fft, ttmp_fft;
	cx_colvec		sq_time, zr_time;
	cx_mat			rep_temp;
	field<cx_mat>	FFT_temp;

	this->Subcarrier_Total = (para.size_fft - (para.size_fft - (para.num_RB * para.num_sc) / 2)) + ((para.num_RB * para.num_sc) / 2);

	if (para.channel.type == Ch_type::AWGN) {
		field<cx_mat> FFT_Temp(1, 1);
		FFT_Temp(0, 0) = zeros<cx_mat>(para.num_Rx_antenna, para.num_Tx_antenna);
		FFT_Temp(0, 0) = this->Output.time;
		this->Output.fft1D = zeros<cx_colvec>(2 * para.num_symb);
		Output.fft1D.fill(FFT_Temp(0, 0)(0, 0).real());

		FFT_Temp(0, 0).reset();		/* memory deallocation */
		FFT_Temp.reset();
	}
	else if (para.channel.type == Ch_type::flat_Rayleigh) {
		field<cx_mat> H_temp(1, 1);
		H_temp(0, 0) = zeros<cx_mat>(para.num_Rx_antenna, para.num_Tx_antenna);
		H_temp(0, 0) = this->Output.time;
		this->Output.fft1D = zeros<cx_colvec>(2 * para.num_symb);
		Output.fft1D.fill(H_temp(0, 0)(0, 0));

		H_temp(0, 0).reset();	/* memory deallocation */
		H_temp.reset();
	}
	else {
		sq_time = zeros<cx_colvec>(this->Output.time.n_slices);
		zr_time = zeros<cx_colvec>(para.size_fft - sq_time.n_elem);
		temp_fft = zeros<cx_colvec>(para.size_fft);

		uword ttmp_fft_size = (para.size_fft - (para.size_fft - (para.num_RB * para.num_sc) / 2)) + ((para.num_RB * para.num_sc) / 2);
		ttmp_fft = zeros<cx_colvec>(ttmp_fft_size);

		this->Output.fft = field<cx_cube>(ttmp_fft_size);
		for (uword i = 0; i < ttmp_fft_size; i++)
			this->Output.fft(i) = zeros<cx_cube>(para.num_symb, para.num_Rx_antenna, para.num_Tx_antenna);

		if (para.channel.time_fading == Ch_mode::Block_Fading) {
			for (uword ind_t = 0; ind_t < para.t.n_elem; ind_t++)
				for (uword Ind_Rx = 0; Ind_Rx < para.num_Rx_antenna; Ind_Rx++)
					for (uword Ind_Tx = 0; Ind_Tx < para.num_Tx_antenna; Ind_Tx++) {
						for (uword i = 0; i < this->Output.time.n_slices; i++)
							sq_time(i) = this->Output.time(Ind_Rx, Ind_Tx, i);
						for (uword i = 0; i < para.size_fft - sq_time.n_elem; i++)
							zr_time(i) = cx_double(0, 0);
						temp_fft = fft(join_cols(sq_time, zr_time));

						uword i, sqi;
						for (i = 0, sqi = para.size_fft - (para.num_RB * para.num_sc) / 2; sqi < para.size_fft; i++, sqi++)
							ttmp_fft(i) = temp_fft(sqi);
						for (sqi = 0; sqi < (para.num_RB * para.num_sc) / 2; i++, sqi++)
							ttmp_fft(i) = temp_fft(sqi);
						rep_temp = repmat(ttmp_fft, 1, para.num_symb);

						this->Output.fft2D = rep_temp;
						for (uword i = 0; i < rep_temp.n_rows; i++)
							for (uword j = 0; j < rep_temp.n_cols; j++)
								this->Output.fft(i)(j, Ind_Rx, Ind_Tx) = rep_temp(i, j);

					}
		}
		else if (para.channel.time_fading == Ch_mode::Fast_Fading) {
			;	// empty statement
		}	/* end of if */

	}

	/* memory deallocations */
	temp_fft.reset();
	ttmp_fft.reset();
	sq_time.reset();
	zr_time.reset();
	rep_temp.reset();

}   /* end of Chan_H_fft() method definition */





/*
 * sinc()
 */
double sinc(const double x) {

	if (x == 0)
		return 1.0;

	return sin(x) / x;

}   /* end of sinc() */



/* end of module_Channel_5G.cpp */