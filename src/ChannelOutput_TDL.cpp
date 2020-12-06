/*
 * ChannelOutput_CDL.cpp
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *		Last Updated: 2018-10-01
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */

#include <complex>
#include <cmath>
#include <armadillo>
#include "sumimo.h"
#include "ChannelOutput_TDL.h"


using namespace arma;



/*
 * ChannelOutput_TDL() - TDL channel in TR 38.901.  7.7.1
 */
void ChannelOutput_TDL(ModuleParameterMIMO &para, colvec &pathDelays, colvec &pathPowers_dB, field<cx_mat> &H) {

	int			Nt = para.num_Tx_antenna;
	int			Nr = para.num_Rx_antenna;
	ch_type		channel_type = para.Channel.type;
	double		user_speed = para.user_speed;
	rowvec		t = para.t;

	int			Num_ray = 20;

	double		freq = para.f;
	double		DS = para.DS;						// Wanted delay spread in ns;
													// 10 * 10e-9; 30 * 10e-9; 300 * 10e-9; 1000 * 10e-9;
													// subject to the frequency f and the scenario type
													// sceType (TR 38.901 Table 7.7.3-1~2)
													// 100ns for calibration

	/* Antenna configuration */
	// Transmit array type (ULA or URA)
	TxArrayType TxArrayType = para.TxArrayType;	// for calibration
	// Receive array type (ULA or URA)
	RxArrayType RxArrayType = para.RxArrayType;	// for calibration
	// Transmit antenna spacing in wavelengths (0.1-100)
	double		TxAnt = para.Tx_d_lambda;
	// Receive antenna spacing in wavelengths (0.1-100)
	double		RxAnt = para.Rx_d_lambda;
	
	double		wl;
	double		omega;
	uword		T;

	cube		Rand_phase;
	mat			Rand_phase_LOS;
	double		N0;
	rowvec		ray;
	rowvec		alpha_n;
	rowvec		omega_n;
	rowvec		beta_n;

	field<cx_mat>	H_pol_each;
	cx_mat		H_temp;

	/* TDL channel type */
	mat TDL = para.TDL;

	uword nTap = size(TDL, 0);		// uword nTap = CDL.n_rows;
	pathDelays = TDL.col(0) * DS;
	pathPowers_dB = TDL.col(1);
	colvec pathPowers = zeros<colvec>(TDL.col(1).n_elem);
	for (uword i = 0; i < TDL.col(1).n_elem; i++)
		pathPowers(i) = pow(10.0, pathPowers_dB(i) / 10.0);

	wl = para.light_speed / freq;
	omega = 2.0 * datum::pi * user_speed / wl;
	T = t.n_elem;
	//  P.Dent, G.E.Beottomley and T.Croft "JAKES FADING MODEL REVISITED"

	H = field<cx_mat>(Nr, Nt);
	for (uword i = 0; i < (uword)Nr; i++)
		for (uword j = 0; j < (uword)Nt; j++)
			H(i, j) = zeros<cx_mat>(T, nTap);

	switch (channel_type) {
	case ch_type::TDL_A:
	case ch_type::TDL_B:
	case ch_type::TDL_C:
		for (uword ind_tap = 0; ind_tap < nTap; ind_tap++) {
			if (para.Tx_pol == 2 && para.Rx_pol == 2)
				Rand_phase = randu(Nr / 2, Nt / 2, Num_ray);
			else if (para.Tx_pol == 2 && para.Rx_pol == 1)
				Rand_phase = randu(Nr, Nt / 2, Num_ray);
			else if (para.Tx_pol == 1 && para.Rx_pol == 2)
				Rand_phase = randu(Nr / 2, Nt, Num_ray);
			else
				Rand_phase = randu(Nr, Nt, Num_ray);

			for (uword b = 0; b < t.n_elem; b++) {
				N0 = Num_ray / 4.0;
				ray = zeros<rowvec>(N0);
				for (uword i = 1; i < N0 + 1; i++)
					ray(i) = i;

				alpha_n = zeros<rowvec>(ray.n_elem);
				alpha_n = 2.0 * datum::pi * (ray - 0.5) / Num_ray;
				omega_n = zeros<rowvec>(alpha_n.n_elem);
				omega_n = omega * cos(alpha_n);
				beta_n = zeros<rowvec>(ray.n_elem);
				beta_n = datum::pi * ray / N0;

				H_pol_each = field<cx_mat>(para.Rx_pol, para.Tx_pol);
				if (para.Tx_pol == 2 && para.Rx_pol == 2) {
					for (uword ind_tx_pol = 0; ind_tx_pol < (uword)para.Tx_pol; ind_tx_pol++)
						for (uword ind_rx_pol = 0; ind_rx_pol < (uword)para.Rx_pol; ind_rx_pol++) {
							H_temp = zeros<cx_mat>(Nr / 2, Nt / 2);
							for (uword ind_ray = 0; ind_ray < N0; ind_ray++)
								H_temp = H_temp + (cos(beta_n(ind_ray)) + cx_double(0, 1) * sin(beta_n(ind_ray))) * ones<mat>(Nr / 2, Nt / 2) * cos(omega_n(ind_ray) * ones<mat>(Nr / 2, Nt / 2) *t(b) + 2.0 * datum::pi * Rand_phase.slice(ind_ray));

							H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
						}
					H(span::all, span::all)(b, ind_tap) = join_vert(join_horiz(H_pol_each(0, 0), H_pol_each(0, 1)), join_horiz(H_pol_each(1, 0), H_pol_each(1, 1)));
				}
				else if (para.Tx_pol == 2 && para.Rx_pol == 1) {
					for (uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
						for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
							H_temp = zeros<cx_mat>(Nr, Nt / 2);
							for (uword ind_ray = 0; ind_ray < N0; ind_ray++)
								H_temp = H_temp + (cos(beta_n(ind_ray)) + cx_double(0, 1) * sin(beta_n(ind_ray))) * ones<mat>(Nr, Nt / 2) * cos(omega_n(ind_ray) * ones<mat>(Nr, Nt / 2) * t(b) + 2.0 * datum::pi * Rand_phase.slice(ind_ray));
						
							H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
						}
					H(span::all, span::all)(b, ind_tap) = join_horiz(H_pol_each(0, 0), H_pol_each(0, 1));
				}
				else if (para.Tx_pol == 1 && para.Rx_pol == 2) {
					for (uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
						for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
							H_temp = zeros<cx_mat>(Nr / 2, Nt);
							for (uword ind_ray = 0; ind_ray < N0; ind_ray++)
								H_temp = H_temp + (cos(beta_n(ind_ray)) + cx_double(0, 1) * sin(beta_n(ind_ray))) * ones<mat>(Nr / 2, Nt) * cos(omega_n(ind_ray) * ones<mat>(Nr / 2, Nt) * t(b) + 2.0 * datum::pi * Rand_phase.slice(ind_ray));

							H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
						}
					H(span::all, span::all)(b, ind_tap) = join_vert(H_pol_each(0, 0), H_pol_each(1, 0));
				}
				else if (para.Tx_pol == 1 && para.Rx_pol == 1) {
					for (uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
						for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
							H_temp = zeros<cx_mat>(Nr, Nt);
							for (uword ind_ray = 0; ind_ray < N0; ind_ray++)
								H_temp = H_temp + (cos(beta_n(ind_ray)) + cx_double(0, 1) * sin(beta_n(ind_ray))) * ones<mat>(Nr, Nt) * cos(omega_n(ind_ray) * ones<mat>(Nr, Nt) * t(b) + 2.0 * datum::pi * Rand_phase.slice(ind_ray));

							H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
						}
					H(span::all, span::all)(b, ind_tap) = H_pol_each(0, 0);
				}
			}
		}
		break;
	case ch_type::TDL_D:
	case ch_type::TDL_E:
		for (uword ind_tap = 0; ind_tap < nTap; ind_tap++) {
			if (para.Tx_pol == 2 && para.Rx_pol == 2) {
				Rand_phase = randu(Nr / 2, Nt / 2, Num_ray);
				Rand_phase_LOS = randu(Nr / 2, Nt / 2);
			}
			else if (para.Tx_pol == 2 && para.Rx_pol == 1) {
				Rand_phase = randu(Nr, Nt / 2, Num_ray);
				Rand_phase_LOS = randu(Nr, Nt / 2);
			}
			else if (para.Tx_pol == 1 && para.Rx_pol == 2) {
				Rand_phase = randu(Nr / 2, Nt, Num_ray);
				Rand_phase_LOS = randu(Nr / 2, Nt);
			}
			else {
				Rand_phase = randu(Nr, Nt, Num_ray);
				Rand_phase_LOS = randu(Nr, Nt);
			}

			for (uword b = 0; b < t.n_elem; b++) {
				N0 = Num_ray / 4.0;
				ray = zeros<rowvec>(N0);
				for (uword i = 1; i < N0 + 1; i++)
					ray(i) = i;

				alpha_n = zeros<rowvec>(ray.n_elem);
				alpha_n = 2.0 * datum::pi * (ray - 0.5) / Num_ray;
				omega_n = zeros<rowvec>(alpha_n.n_elem);
				omega_n = omega * cos(alpha_n);
				beta_n = zeros<rowvec>(ray.n_elem);
				beta_n = datum::pi * ray / N0;

				H_pol_each = field<cx_mat>(para.Rx_pol, para.Tx_pol);
				if (para.Tx_pol == 2 && para.Rx_pol == 2) {
					if (ind_tap == 1)
						for (uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
							for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
								H_temp = (cos(0) + cx_double(0, 1) * sin(0)) * ones<mat>(Nr / 2, Nt / 2) * cos(omega * ones<mat>(Nr / 2, Nt / 2) * t(b) + 2.0 * datum::pi * Rand_phase_LOS);
								H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
							}
					else {
						for (uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
							for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
								H_temp = zeros<cx_mat>(Nr / 2, Nt / 2);
								for (uword ind_ray = 0; ind_ray < N0; ind_ray++)
									H_temp = H_temp + (cos(beta_n(ind_ray)) + cx_double(0, 1) * sin(beta_n(ind_ray))) * ones<mat>(Nr / 2, Nt / 2) * cos(omega_n(ind_ray) * ones<mat>(Nr / 2, Nt / 2) * t(b) + 2.0 * datum::pi * Rand_phase.slice(ind_ray));
								H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
							}
						H(span::all, span::all)(b, ind_tap) = join_vert(join_horiz(H_pol_each(0, 0), H_pol_each(0, 1)), join_horiz(H_pol_each(1, 0), H_pol_each(1, 1)));
					}
				}
				else if (para.Tx_pol == 2 && para.Rx_pol == 1) {
					if (ind_tap == 1)
						for (uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
							for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
								H_temp = (cos(0) + cx_double(0, 1)) * ones<mat>(Nr, Nt / 2) * cos(omega * ones<mat>(Nr / 2, Nt / 2) * t(b) + 2.0 * datum::pi * Rand_phase_LOS);
								H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
							}
					else {
						for (uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
							for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
								H_temp = zeros<cx_mat>(Nr, Nt / 2);
								for (uword ind_ray = 0; ind_ray < N0; ind_ray++)
									H_temp = H_temp + (cos(beta_n(ind_ray)) + cx_double(0, 1) * sin(beta_n(ind_ray))) * ones<mat>(Nr, Nt / 2) * cos(omega_n(ind_ray) * ones<mat>(Nr, Nt / 2) * t(b) + 2.0 * datum::pi * Rand_phase.slice(ind_ray));
								H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
							}
						H(span::all, span::all)(b, ind_tap) = join_horiz(H_pol_each(0, 0), H_pol_each(0, 1));
					}
				}
				else if (para.Tx_pol == 1 && para.Rx_pol == 2) {
					if (ind_tap == 1)
						for (uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
							for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
								H_temp = (cos(0) + cx_double(0, 1) * sin(0)) * ones<mat>(Nr / 2, Nt) * cos(omega * ones<mat>(Nr / 2, Nt) * t(b) + 2.0 * datum::pi * Rand_phase_LOS);
								H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
							}
					else {
						for(uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
							for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
								H_temp = zeros<cx_mat>(Nr / 2, Nt);
								for (uword ind_ray = 0; ind_ray < N0; ind_ray++)
									H_temp = H_temp + (cos(beta_n(ind_ray)) + cx_double(0, 1) * sin(beta_n(ind_ray))) * ones<mat>(Nr / 2, Nt) * cos(omega_n(ind_ray) * ones<mat>(Nr / 2, Nt) * t(b) + 2.0 * datum::pi * Rand_phase.slice(ind_ray));
								H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
							}
						H(span::all, span::all)(b, ind_tap) = join_vert(H_pol_each(0, 0), H_pol_each(1, 0));
					}
				}
				else if (para.Tx_pol == 1 && para.Rx_pol == 1) {
					if (ind_tap == 1)
						for (uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
							for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
								H_temp = (cos(0) + cx_double(0, 1) * sin(0)) * ones<mat>(Nr, Nt) * cos(omega * ones<mat>(Nr, Nt) * t(b) + 2.0 * datum::pi * Rand_phase_LOS);
								H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
							}
					else {
						for (uword ind_tx_pol = 0; ind_tx_pol < para.Tx_pol; ind_tx_pol++)
							for (uword ind_rx_pol = 0; ind_rx_pol < para.Rx_pol; ind_rx_pol++) {
								H_temp = zeros<cx_mat>(Nr, Nt);
								for (uword ind_ray = 0; ind_ray < N0; ind_ray++)
									H_temp = H_temp + (cos(beta_n(ind_ray)) * cx_double(0, 1) * sin(beta_n(ind_ray))) * ones<mat>(Nr, Nt) * cos(omega_n(ind_ray) * ones<mat>(Nr, Nt) * t(b) + 2.0 * datum::pi * Rand_phase.slice(ind_ray));
								H_pol_each(ind_rx_pol, ind_tx_pol) = sqrt(pathPowers(ind_tap) * 2.0 / N0) * H_temp;
							}
					}
					H(span::all, span::all)(b, ind_tap) = H_pol_each(0, 0);
				}
			}
		}
		break;
	}

	/* local memory deallocation */
	t.reset();
	Rand_phase.reset();
	Rand_phase_LOS.reset();
	ray.reset();
	alpha_n.reset();
	beta_n.reset();
	omega_n.reset();
	for (uword i = 0; i < H_pol_each.n_elem; i++)
		H_pol_each(i).reset();
	H_temp.reset();
	H_pol_each.reset();
	H_temp.reset();
	TDL.reset();
	pathPowers.reset();

}   /* end of ChannelOutput_TDL() */



/* end of ChannelOutput_CDL.cpp */