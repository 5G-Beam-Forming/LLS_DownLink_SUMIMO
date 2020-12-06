/*
 * ChannelOutput_CDL.cpp
 *
 *		Converted C++ code from MATLAB script code.
 *		Converted by PJW
 *		Last Updated: 2019-01-10
 *
 * Copyright(c) 2016-2018 KAIST. All rights reserved.
 */
#include <complex>
#include <cmath>
#include <armadillo>
#include "sumimo.h"
#include "ChannelOutput_CDL.h"


using namespace arma;



/* local function declaration */
void	GCS_to_LCS(double alpha, double beta, double gamma, double theta_GCS, double phi_GCS, double &theta_LCS, double &phi_LCS);
void	GCS_to_LCS2(double alpha, double beta, double gamma, colvec theta_GCS, colvec phi_GCS, colvec &theta_LCS, colvec &phi_LCS);
void	GCS_to_LCS3(double, double, double, colvec, colvec, vector<double> &, vector<double> &);
double  sind(double);
double  cosd(double);
double  acosd(double);
double  angle(complex<double>);

colvec  repelem(double, uword);
colvec  repvect(colvec, uword);
rowvec  randperm(uword);
mat     Conf_AzimuthAngles(mat);
mat     Conf_ZenithAngles(mat);
mat     sph_to_car(mat, mat);
rowvec  sph_to_carD(double, double);
mat		BS_anttena_patterns(mat, mat);
mat		UE_antenna_patterns(mat, mat);



/*
 * ChannelOutput_CDL()
 */
void ChannelOutput_CDL(ModuleParameterMIMO &para, colvec &pathDelays, colvec &pathPowers_dB, field<cx_mat> &H, angleFirst &Ang_strong_cluster) {

	uword		i, k;									// loop and vector index

	int			Nt = para.num_Tx_antenna;
	int			Nr = para.num_Rx_antenna;
	Ch_type		channel_type = para.channel.type;
	double		user_speed = para.user_speed;
	double		theta_v = para.theta_v;
	double		phi_v = para.phi_v;
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
	TxArrayType TxArrayType = para.txArrayType;	// for calibration
	// Receive array type (ULA or URA)
	RxArrayType RxArrayType = para.rxArrayType;	// for calibration
	// Transmit antenna spacing in wavelengths (0.1-100)
	double		TxAnt = para.Tx_d_lambda;
	// Receive antenna spacing in wavelengths (0.1-100)
	double		RxAnt = para.Rx_d_lambda;
	
	// Number of transmit antenna elements per row for URA i.e. # of columns
	// N in BS; M * Wt * p = Nt int the input
	int			Wt = para.N1;						// for calibration
	// Number of receive antenna elements per row for URA i.e. # of columns
	// N in UE; M * Wr * P = Nr in the input
	// Wr = para->num_Rx_antenna;
	int			Wr = para.M1;

	/* ration AS_desired/AS_model parameter */
	int			AS_ratio = para.AS_ratio;

	double		mu_model_AoD;
	double		mu_model_AoA;
	double		mu_model_ZoD;
	double		mu_model_ZoA;

	double		mu_desired_AoD;
	double		mu_desired_AoA;
	double		mu_desired_ZoD;
	double		mu_desired_ZoA;

	colvec		AODs_temp;
	colvec		ZODs_temp;
	colvec		AOAs_temp;
	colvec		ZOAs_temp;

	double		correctAzAOA;
	double		AzAOA;
	mat			azAOA_temp;
	double		correctElAOA;
	double		ElAOA;

	mat			rhat_rx;
	rowvec		vbar;

	double		wl;
	colvec		n_BS;
	colvec		n_MS;

	cx_mat		a_BS;
	cx_mat		a_MS;

	colvec		ZoDs_LCS;
	colvec		AoDs_LCS;
	colvec		ZoAs_LCS;
	colvec		AoAs_LCS;

	cx_colvec	dopplerTerm;
	uword		T;

	mat			Tx_patterns;
	mat			Rx_patterns;

	double		XPR;

	mat			Initial_Phase;
	cx_mat		H1, H2, H3, H4, Hmid;
	cx_rowvec	Hmid2(2);
	rowvec		Hrv1(2), Hrv2(2);
	cx_colvec	Hcv1(2);
	cx_mat		Hmat(2, 2);


	/*
	 * CDL Channel Type and related parameters
	 */
	/* TR 38.901 Table 7.5-3 */
	colvec offset_ang;
	offset_ang << 0.0447 << -0.0447 << 0.1413 << -0.1413 << 0.2492 << -0.2492 << 0.3715 << -0.3715 << 0.5129 << -0.5129
			   << 0.6797 << -0.6797 << 0.8844 << -0.8844 << 1.1481 << -1.1481 << 1.5195 << -1.5195 << 2.1551 << -2.1551 << endr;

	/* CDL channel type */
	mat CDL = para.CDL;

	uword nTap = size(CDL, 0);		// uword nTap = CDL.n_rows;
	pathDelays = CDL.col(0) * DS;
	pathPowers_dB = CDL.col(1);
	colvec pathPowers = zeros<colvec>(CDL.col(1).n_elem);
	for (uword i = 0; i < CDL.col(1).n_elem; i++)
		pathPowers(i) = pow(10.0, pathPowers_dB(i) / 10.0);

	/* R1-1700144 */
	switch (channel_type) {
	case Ch_type::CDL_A:	// R1-1700144
		mu_model_AoD =   -2.05;
		mu_model_AoA = -164.40;
		mu_model_ZoD =   97.27;
		mu_model_ZoA =   87.68;
		break;
	case Ch_type::CDL_B:	// R1-1700144
		mu_model_AoD =  -4.96;
		mu_model_AoA = 176.60;
		mu_model_ZoD = 105.64;
		mu_model_ZoA =  71.76;
		break;
	case Ch_type::CDL_C:	// calculated by using R1-1700144
		mu_model_AoD = -19.6808;
		mu_model_AoA = 149.8943;
		mu_model_ZoD =  99.3310;
		mu_model_ZoA =  73.4974;
		break;
	case Ch_type::CDL_D:	// calculated by using R1-1700144
		mu_model_AoD =   1.9930;
		mu_model_AoA = 178.1941;
		mu_model_ZoD =  98.0897;
		mu_model_ZoA =  81.5213;
		break;
	case Ch_type::CDL_E:	// calculated by using R1-1700144
		mu_model_AoD =   2.8396;
		mu_model_AoA = 179.2480;
		mu_model_ZoD =  99.7957;
		mu_model_ZoA =  80.5322;
		break;
	}

	/* Default : For calibration in R1-17-1823 */
	mu_desired_AoD = 120.0 * (randu<double>() - 0.5);
	mu_desired_AoA = 360.0 * (randu<double>() - 0.5);
	mu_desired_ZoD =  45.0 * (randu<double>() + 2.0);
	mu_desired_ZoA =  45.0 * (randu<double>() + 1.0);

	AODs_temp = CDL.col(2);
	ZODs_temp = CDL.col(4);
	AOAs_temp = CDL.col(3);
	ZOAs_temp = CDL.col(5);

	AODs_temp = AS_ratio * (AODs_temp - mu_model_AoD) + mu_desired_AoD;
	ZODs_temp = AS_ratio * (ZODs_temp - mu_model_ZoD) + mu_desired_ZoD;
	AOAs_temp = AS_ratio * (AOAs_temp - mu_model_AoA) + mu_desired_AoA;
	ZOAs_temp = AS_ratio * (ZOAs_temp - mu_model_ZoA) + mu_desired_ZoA;
	
	switch (channel_type) {
	case Ch_type::CDL_A:
		GCS_to_LCS(0, para.Tx_downtilt, 0, ZODs_temp(1), AODs_temp(1), Ang_strong_cluster.ZOD, Ang_strong_cluster.AOD);
		break;
	case Ch_type::CDL_C:
		GCS_to_LCS(0, para.Tx_downtilt, 0, ZODs_temp(5), AODs_temp(5), Ang_strong_cluster.ZOD, Ang_strong_cluster.AOD);
		break;
	case Ch_type::CDL_B:
	case Ch_type::CDL_D:
	case Ch_type::CDL_E:
		GCS_to_LCS(0, para.Tx_downtilt, 0, ZODs_temp(0), AODs_temp(0), Ang_strong_cluster.ZOD, Ang_strong_cluster.AOD);
		break;
	}

	double Pi_UT_alpha;
	Pi_UT_alpha = 360.0 * randu<double>();

	mat AODs = zeros<mat>(nTap * Num_ray, 1);
	mat ZODs = zeros<mat>(nTap * Num_ray, 1);
	mat AOAs = zeros<mat>(nTap * Num_ray, 1);
	mat ZOAs = zeros<mat>(nTap * Num_ray, 1);
	
	colvec offsetAngM = zeros<colvec>(Num_ray);
	uvec ridx;
	for (uword ind_tap = 0; ind_tap < nTap; ind_tap++) {
		ridx = zeros<uvec>(Num_ray * (ind_tap + 1) - (1 + Num_ray * ind_tap) + 1);
		for (uword i = 0, idx = Num_ray * ind_tap; idx < Num_ray * (ind_tap + 1); i++, idx++)
			ridx(i) = idx;

		AODs(ridx) = repelem(AODs_temp(ind_tap), Num_ray) + offset_ang(randi<uvec>(Num_ray, distr_param(0, Num_ray - 1))) * para.c_ASD;
		ZODs(ridx) = repelem(ZODs_temp(ind_tap), Num_ray) + offset_ang(randi<uvec>(Num_ray, distr_param(0, Num_ray - 1))) * para.c_ZSD;
		AOAs(ridx) = repelem(AOAs_temp(ind_tap), Num_ray) + offset_ang(randi<uvec>(Num_ray, distr_param(0, Num_ray - 1))) * para.c_ASA;
		ZOAs(ridx) = repelem(ZOAs_temp(ind_tap), Num_ray) + offset_ang(randi<uvec>(Num_ray, distr_param(0, Num_ray - 1))) * para.c_ZSA;
	}

	AODs = repvect(AODs_temp, Num_ray) + repmat(offset_ang * para.c_ASD, nTap, 1);
	ZODs = repvect(ZODs_temp, Num_ray) + repmat(offset_ang * para.c_ZSD, nTap, 1);
	AOAs = repvect(AOAs_temp, Num_ray) + repmat(offset_ang * para.c_ASA, nTap, 1);
	ZOAs = repvect(ZOAs_temp, Num_ray) + repmat(offset_ang * para.c_ZSA, nTap, 1);

	/* Confine azimuth angles in [-180,180] */
	AOAs = Conf_AzimuthAngles(AOAs);
	AODs = Conf_AzimuthAngles(AODs);

	/* Confine zenith angles in [0,360] and map [180,360] to [180,0] */
	ZOAs = Conf_ZenithAngles(ZOAs);
	ZODs = Conf_ZenithAngles(ZODs);

	/* For LOS environment, adjust subpath AoDs and AoAs such that the AoD
	 * and AoA of the LOS component are aligned in line */
	if (channel_type == Ch_type::CDL_D || channel_type == Ch_type::CDL_E) {
		// Calculate the correct azimuth AoA for LOS case, which differ from azimuth AoD by 180 degrees
		// Channel_Info(1, 4) denotes the azimuth AoD for LOS component
		correctAzAOA = 180.0;		// AoA = 0 for both of CDL_D and CDL_E.
		// Calculate the difference between the generated azimuth AoA and the correct azimuth AoA
		AzAOA = AOAs(0) - correctAzAOA;
		AOAs = AOAs - AzAOA;
		azAOA_temp = zeros<mat>(AOAs.size());
		azAOA_temp = AOAs;
		azAOA_temp(azAOA_temp < 0.0) = azAOA_temp(azAOA_temp < 0) + 360.0;
		AOAs = azAOA_temp;

		// Calculate the correct elevation AoA for LOS component, which is the additive inverse of the corresponding elevation AoD
		correctElAOA = -ZODs(0);
		// Calculate the difference between the generated elevation AoA and the correct elevation AoA
		ElAOA = ZOAs(0) - correctElAOA;
		ZOAs = ZOAs - ElAOA;
	}   /* end of if */

	rhat_rx = sph_to_car(AOAs, ZOAs);

	/* TR 38.901 Equation 7.5 - 22 / 7.5 - 28 / 7.5 - 29, exponentional term which is a
	 * function of the Doppler due to user velocity */

	vbar = user_speed * sph_to_carD(phi_v, theta_v);    // UT velocity vector

	wl = para.light_speed / freq;

	n_BS = zeros<colvec>(para.N1 * para.N2);
	n_MS = zeros<colvec>(para.M1 * para.M2);
	for (uword i = 0; i < para.N1 * para.N2; i++)
		n_BS(i) = i + 1.0;
	for (uword i = 0; i < para.M1 * para.M2; i++)
		n_MS(i) = i + 1.0;


	/* LCS of AoDs and ZoDs with BS downtilt */
	GCS_to_LCS2(0, para.Tx_downtilt, 0, ZODs, AODs, ZoDs_LCS, AoDs_LCS);
	GCS_to_LCS2(Pi_UT_alpha, 0, 0, ZOAs, AOAs, ZoAs_LCS, AoAs_LCS);

	a_BS = zeros<cx_mat>(para.N1 * para.N2, AODs.size());
	a_MS = zeros<cx_mat>(para.M1 * para.M2, AOAs.size());

	switch (TxArrayType) {
	case TxArrayType::ULA:
		for (i = 0; i < para.N1 * para.N2; i++)
			for (k = 0; k < AODs.size(); k++)
				a_BS(i, k) = exp(cx_double(0, 1) * (n_BS(i) - 1) * 2.0 * datum::pi * TxAnt * sind(AoDs_LCS(k)));
		break;
	case TxArrayType::URA:
		for (i = 0; i < para.N1 * para.N2; i++)
			for (k = 0; k < AODs.size(); k++)
				a_BS(i, k) = exp(cx_double(0, 1) * 2.0 * datum::pi * TxAnt * (dmod(i, Wt) * cosd(ZoDs_LCS(k)) + fix((n_BS(i) - 1) / Wt) * sind(ZoDs_LCS(k)) * sind(AoDs_LCS(k))));
	}

	switch (RxArrayType) {
	case RxArrayType::ULA:
		for (i = 0; i < para.M1 * para.M2; i++)
			for (k = 0; k < AOAs.size(); k++)
				a_MS(i, k) = exp(cx_double(0, 1) * (n_MS(i) - 1) * 2.0 * datum::pi * RxAnt * sind(AoAs_LCS(k)));
		break;
	case RxArrayType::URA:
		for (i = 0; i < para.M1 * para.M2; i++)
			for (k = 0; k < AOAs.size(); k++)
				a_MS(i, k) = exp(cx_double(0, 1) * 2.0 * datum::pi * RxAnt * (dmod(i, Wr) * cosd(ZoAs_LCS(k)) + fix((n_MS(i) - 1.0) / Wr) * sind(ZoAs_LCS(k)) * sind(AoAs_LCS(k))));
	}

	/* last term of 38.901. eq. 7.5-22 */
	dopplerTerm = zeros<cx_colvec>(rhat_rx.n_rows);
	dopplerTerm = exp(kron(t, (cx_double(0, 1) * 2.0 * datum::pi * rhat_rx * vbar.t() / wl)));

	T = t.size();


	switch (para.tx_pattern_type) {
	case Tx_pattern_type::Omni_directional:
		Tx_patterns = ones<mat>(1, Num_ray * nTap);
		break;
	case Tx_pattern_type::Pattern:
		Tx_patterns = BS_anttena_patterns(ZoDs_LCS, AoDs_LCS);
		break;
	}

	switch (para.rx_pattern_type) {
	case Rx_pattern_type::Omni_directional:
		Rx_patterns = ones<mat>(1, Num_ray * nTap);
		break;
	case Rx_pattern_type::Pattern:
		Rx_patterns = UE_antenna_patterns(ZoAs_LCS, AoAs_LCS);
		break;
	}


	/* TR 38.901 7.7.1. step 3 */
	XPR = pow(10.0, (para.XPR_dB / 10.0));


	/* TR 38.901 7.7.1. step 4 (7.5. step 10, 11) */
	H = field<cx_mat>(Nr, Nt);
	for (uword i = 0; i < Nr; i++)
		for (uword j = 0; j < Nt; j++)
			H(i, j) = zeros<cx_mat>(T, nTap);
	Initial_Phase = zeros<mat>(para.Tx_pol, para.Rx_pol);

	for (uword ind_tap = 0; ind_tap < nTap; ind_tap++) {
		for (uword b = 0; b < t.n_elem; b++) {
			for (uword ind_ray = 0; ind_ray < Num_ray; ind_ray++) {

				if (para.Tx_pol == 1 && para.Rx_pol == 1)
					Initial_Phase = 2.0 * datum::pi * (randu<double>() - 1.0 / 2.0);
				else
					Initial_Phase = 2.0 * datum::pi * (randu<mat>(para.Tx_pol, para.Rx_pol) - 1.0 / 2.0);

				if (para.Tx_pol == 2 && para.Rx_pol == 2) {
					Hrv1(0) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle(0));
					Hrv1(1) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle(0));
					Hmat(0, 0) = exp(cx_double(0, 1) * Initial_Phase(0));
					Hmat(0, 1) = sqrt(1 / XPR) * exp(cx_double(0, 1) * Initial_Phase(1));
					Hmat(1, 0) = sqrt(1 / XPR) * exp(cx_double(0, 1) * Initial_Phase(2));
					Hmat(1, 1) = exp(cx_double(0, 1) * Initial_Phase(3));
					Hrv2(0) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle(0));
					Hrv2(1) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle(0));

					H1 = Hrv1 * Hmat * Hrv2.t() * dopplerTerm(ind_ray, b) * a_MS.col(ind_ray + Num_ray * ind_tap) * a_BS.col(ind_ray + Num_ray * ind_tap).t();

					Hrv1(0) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle(1));
					Hrv1(1) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle(1));
					Hmat(0, 0) = exp(cx_double(0, 1) * Initial_Phase(0));
					Hmat(0, 1) = sqrt(1 / XPR) * exp(cx_double(0, 1) * Initial_Phase(1));
					Hmat(1, 0) = sqrt(1 / XPR) * exp(cx_double(0, 1) * Initial_Phase(2));
					Hmat(1, 1) = exp(cx_double(0, 1) * Initial_Phase(3));
					Hrv2(0) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle(0));
					Hrv2(1) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle(0));

					H2 = Hrv1 * Hmat * Hrv2.t() * dopplerTerm(ind_ray, b) * a_MS.col(ind_ray + Num_ray * ind_tap) * a_BS.col(ind_ray + Num_ray * ind_tap).t();

					Hrv1(0) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle(0));
					Hrv1(1) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle(0));
					Hmat(0, 0) = exp(cx_double(0, 1) * Initial_Phase(0));
					Hmat(0, 1) = sqrt(1 / XPR) * exp(cx_double(0, 1) * Initial_Phase(1));
					Hmat(1, 0) = sqrt(1 / XPR) * exp(cx_double(0, 1) * Initial_Phase(2));
					Hmat(1, 1) = exp(cx_double(0, 1) * Initial_Phase(3));
					Hrv2(0) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle(1));
					Hrv2(1) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle(1));

					H3 = Hrv1 * Hmat * Hrv2.t() * dopplerTerm(ind_ray, b) * a_MS.col(ind_ray + Num_ray * ind_tap) * a_BS.col(ind_ray + Num_ray * ind_tap).t();

					Hrv1(0) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle(1));
					Hrv1(1) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle(1));
					Hmat(0, 0) = exp(cx_double(0, 1) * Initial_Phase(0));
					Hmat(0, 1) = sqrt(1 / XPR) * exp(cx_double(0, 1) * Initial_Phase(1));
					Hmat(1, 0) = sqrt(1 / XPR) * exp(cx_double(0, 1) * Initial_Phase(2));
					Hmat(1, 1) = exp(cx_double(0, 1) * Initial_Phase(3));
					Hrv2(0) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle(1));
					Hrv2(1) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle(1));

					H4 = Hrv1 * Hmat * Hrv2.t() * dopplerTerm(ind_ray, b) * a_MS.col(ind_ray + Num_ray * ind_tap) * a_BS.col(ind_ray + Num_ray * ind_tap).t();

					Hmid = join_vert(join_horiz(H1, H3), join_horiz(H2, H4));

					for (uword i = 0; i < Nr; i++)
						for (uword j = 0; j < Nt; j++)
							H(i, j)(b, ind_tap) = H(i, j)(b, ind_tap) + Hmid(i, j);

					pow(norm(H1), 2);
					pow(norm(H3), 2);
					pow(norm(H2), 2);
					pow(norm(H4), 2);
				}
				else if (para.Tx_pol == 2 && para.Rx_pol == 1) {
					Hmat(0) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(0);
					Hrv1(0) = exp(imag(Initial_Phase(0)));
					Hrv1(1) = sqrt(1 / XPR) * exp(imag(Initial_Phase(1)));
					Hrv2(0) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle[0]);
					Hrv2(1) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle[0]);

					H1 = Hmat(0) * Hrv1 * Hrv2.t() * dopplerTerm(ind_ray, b) * a_MS.col(ind_ray + Num_ray * ind_tap) * a_BS.col(ind_ray + Num_ray * ind_tap).t();

					Hmat(0) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(0);
					Hrv1(0) = exp(imag(Initial_Phase(0)));
					Hrv1(1) = sqrt(1 / XPR) * exp(imag(Initial_Phase(1)));
					Hrv2(0) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle[1]);
					Hrv2(1) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle[1]);

					H3 = Hmat(0) * Hrv1 * Hrv2.t() * dopplerTerm(ind_ray, b) * a_MS.col(ind_ray + Num_ray * ind_tap) * a_BS.col(ind_ray + Num_ray * ind_tap).t();

					Hmid2(0) = H1(0, 0);
					Hmid2(1) = H3(0, 0);
					H(0, 0)(b, ind_tap) += H1(0, 0);
				}
				else if (para.Tx_pol == 1 && para.Rx_pol == 2) {
					Hrv1(0) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle[0]);
					Hrv1(1) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle[0]);
					Hcv1(0) = exp(imag(Initial_Phase(0)));
					Hcv1(1) = sqrt(1 / XPR) * exp(imag(Initial_Phase(1)));
					Hrv2(0) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(0);

					H1 = Hrv1 * Hcv1 * Hrv2(0) * dopplerTerm(ind_ray, b) * a_MS.col(ind_ray + Num_ray * ind_tap) * a_BS.col(ind_ray + Num_ray * ind_tap).t();

					Hrv1(0) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(para.pol_slant_angle[1]);
					Hrv1(1) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * sind(para.pol_slant_angle[1]);
					Hcv1(0) = exp(imag(Initial_Phase(0)));
					Hcv1(1) = sqrt(1 / XPR) * exp(imag(Initial_Phase(1)));
					Hrv2(0) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(0);

					H2 = Hrv1 * Hcv1 * Hrv2(0) * dopplerTerm(ind_ray, b) * a_MS.col(ind_ray + Num_ray * ind_tap) * a_BS.col(ind_ray + Num_ray * ind_tap).t();

					Hmid2(0) = H1(0, 0);
					Hmid2(1) = H2(0, 0);
					H(0, 0)(b, ind_tap) += H1(0, 0);
				}
				else if (para.Tx_pol == 1 && para.Rx_pol == 1) {
					Hrv1.fill(0); Hrv2.fill(0); Hcv1.fill(0);
					Hrv1(0) = sqrt(Rx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(0);
					Hcv1(0) = exp(cx_double(0, 1) * Initial_Phase(0));
					Hrv2(0) = sqrt(Tx_patterns(ind_ray + Num_ray * ind_tap)) * cosd(0);

					H(0, 0)(b, ind_tap) = H(0, 0)(b, ind_tap) + Hrv1(0) * Hcv1(0) * Hrv2(0) * dopplerTerm(ind_ray, b) * a_MS(ind_ray + Num_ray * ind_tap) * a_BS(ind_ray + Num_ray * ind_tap);
				}

			}   /* end of for (uword ind_ray = 0; ind_ray < Num_ray; ind_ray++) */

			for (uword i = 0; i < Nr; i++)
				for (uword j = 0; j < Nt; j++)
					H(i, j)(b, ind_tap) = sqrt(pathPowers(ind_tap) / Num_ray) * H(i, j)(b, ind_tap);

		}   /* end of for (uword b = 0; b < t.size(); b++) */

	}   /* end of for (uword ind_tap = 0; ind_tap < nTap; ind_tap++) */

	/* memory deallocations */
	t.reset();
	AODs_temp.reset();	ZODs_temp.reset();	AOAs_temp.reset();	ZOAs_temp.reset();
	azAOA_temp.reset();
	rhat_rx.reset();
	vbar.reset();
	n_BS.reset();	n_MS.reset();	a_BS.reset();	a_MS.reset();
	ZoDs_LCS.reset();	AoDs_LCS.reset();	ZoAs_LCS.reset();	AoAs_LCS.reset();
	dopplerTerm.reset();
	Tx_patterns.reset();	Rx_patterns.reset();
	Initial_Phase.reset();
	H1.reset(); H2.reset(); H3.reset(); H4.reset(); Hmid.reset(); Hmid2.reset();
	Hrv1.reset(); Hrv2.reset(); Hcv1.reset(); Hmat.reset();
	offset_ang.reset();
	CDL.reset();
	pathPowers.reset();
	offsetAngM.reset();
	ridx.reset();

}   /* end of ChannelOutput_CDL() */





/*
 * GCs_to_LCS() - double type operands(function parameters)
 *    - input parameter: double alpha, beta, gamma, colvec theta_GCS, phi_GCS
 *	  - output parameter: colvec theta_LCS, phi_LCS
 */
void GCS_to_LCS(double alpha, double beta, double gamma, double theta_GCS, double phi_GCS, double &theta_LCS, double &phi_LCS) {

	complex<double> phi;
	double angphi;
	double phir, phii;

	theta_LCS = acosd(cosd(beta) * cosd(gamma) * cosd(theta_GCS) + (sind(beta) * cosd(gamma) * cosd(phi_GCS - alpha) - sind(gamma) * sind(phi_GCS - alpha)) * sind(theta_GCS));
	phir = cosd(beta)*sind(theta_GCS)*cosd(phi_GCS - alpha) - sind(beta)*cosd(theta_GCS);
	phii = cosd(beta)*sind(gamma)*cosd(theta_GCS) + (sind(beta)*sind(gamma)*cosd(phi_GCS - alpha) + cosd(gamma)*sind(phi_GCS - alpha))*sind(theta_GCS);
	phi = complex<double>(phir, phii);
	angphi = angle(phi);
	phi_LCS = angphi / datum::pi * 180.0;

}   /* end of GCS_to_LCS() */





/*
 * GCS_to_LCS2 - arma::colvec type operands(function parameters)
 *    - input parameter: double alpha, beta, gamma, colvec theta_GCS, phi_GCS
 *	  - output parameter: colvec theta_LCS, phi_LCS
 */
void GCS_to_LCS2(double alpha, double beta, double gamma, colvec theta_GCS, colvec phi_GCS, colvec &theta_LCS, colvec &phi_LCS) {
	
	cx_double phi;
	double angphi;
	double phir, phii;
	colvec tLCS = zeros<colvec>(theta_GCS.n_elem);
	colvec pLCS = zeros<colvec>(theta_GCS.n_elem);

	for (uword i = 0; i < theta_GCS.n_elem; i++) {
		tLCS(i) = acosd(cosd(beta)*cosd(gamma)*cosd(theta_GCS(i)) + (sind(beta)*cosd(gamma)*cosd(phi_GCS(i) - alpha) - sind(gamma)*sind(phi_GCS(i) - alpha))*sind(theta_GCS(i)));
		phir = cosd(beta)*sind(theta_GCS(i))*cosd(phi_GCS(i) - alpha) - sind(beta)*cosd(theta_GCS(i));
		phii = cosd(beta)*sind(gamma)*cosd(theta_GCS(i)) + (sind(beta)*sind(gamma)*cosd(phi_GCS(i) - alpha) + cosd(gamma)*sind(phi_GCS(i) - alpha))*sind(theta_GCS(i));
		phi = cx_double(phir, phii);
		angphi = angle(phi);
		pLCS(i) = angphi / datum::pi * 180.0;
	}

	theta_LCS = tLCS;
	phi_LCS   = pLCS;

	tLCS.reset();
	pLCS.reset();

}   /* end of GCS_to_LCS() */





/*
 * GCS_to_LCS3 - std::vector type operands(function parameters)
 *    - input parameter: double alpha, beta, gamma, colvec theta_GCS, phi_GCS
 *	  - output parameter: colvec theta_LCS, phi_LCS
 */
void GCS_to_LCS3(double alpha, double beta, double gamma, colvec theta_GCS, colvec phi_GCS, vector<double> &theta_LCS, vector<double> &phi_LCS) {

	complex<double> phi;
	double angphi;
	double phir, phii;
	vector<double> tLCS(theta_GCS.n_elem);
	vector<double> pLCS(theta_GCS.n_elem);

	for (uword i = 0; i < theta_GCS.n_elem; i++) {
		tLCS[i] = acosd(cosd(beta)*cosd(gamma)*cosd(theta_GCS(i)) + (sind(beta)*cosd(gamma)*cosd(phi_GCS(i) - alpha) - sind(gamma)*sind(phi_GCS(i) - alpha))*sind(theta_GCS(i)));
		phir = cosd(beta)*sind(theta_GCS(i))*cosd(phi_GCS(i) - alpha) - sind(beta)*cosd(theta_GCS(i));
		phii = cosd(beta)*sind(gamma)*cosd(theta_GCS(i)) + (sind(beta)*sind(gamma)*cosd(phi_GCS(i) - alpha) + cosd(gamma)*sind(phi_GCS(i) - alpha))*sind(theta_GCS(i));
		phi = complex<double>(phir, phii);
		angphi = angle(phi);
		pLCS[i] = angphi / datum::pi * 180.0;
	}
	
	theta_LCS = tLCS;
	phi_LCS = pLCS;

	tLCS.clear();
	pLCS.clear();

}   /* end of GCS_to_LCS() */





/*
 * MATLAB built-in functions implemented with C++
 *****************************************************************/

/* sind(X) - Sine of argument in degrees. X is degree. */
double sind(double x) {

	return sin(x * datum::pi / 180.0);

}



/* cosd(x) - Cosine of argument in degrees. x is degree. */
double cosd(double x) {

	return cos(x * datum::pi / 180.0);

}



/* acsd(X) - Inverse consine, result in degrees. */
double acosd(double x) {

	return acos(x) * 180.0 / datum::pi;

}



/* angle(Z) - Phase angle */
double angle(cx_double z) {

	return atan2(z.imag(), z.real());

}



/*
 * repelems() - repeats double type scalar elements
 *		o input:
 *        - double data: scalar elements
 *		  - uword n: times to repeats
 *		o output: colvec res
 */
colvec repelem(double data, uword n) {
	
	colvec res = zeros<colvec>(n);
	
	res.fill(data);
	
	return res;
}




/*
 * repvect() - repeats double type column vector
 *		o input:
 *		  - colvec data: column vector
 *		  - uword n: tiomes to repeats
 *		o output: colvec res
 */
colvec repvect(colvec data, uword n) {

	uword vecSize = data.n_elem;
	colvec res = zeros<colvec>(vecSize * n);
	
	for (uword i = 0, k = 0; i < vecSize; i++)
		for (uword j = 0; j < n; j++, k++)
			res(k) = data(i);

	return res;
}




/*
 * Conf_AzimuthAngles()
 *		- input: mat phi
 */
mat Conf_AzimuthAngles(mat phi) {

	for (uword i = 0; i < phi.n_elem; i++)
		phi(i) = dmod((phi(i) + 180.0), 360.0) - 180.0;

	return phi;

}   /* end of Conf_AzimuthAngles() */




/*
 * Conf_ZenithAngles()
 */
mat Conf_ZenithAngles(mat theta) {

	double ta;

	for (uword i = 0; i < theta.n_elem; i++) {
		ta = dmod(theta(i), 360.0);

		if (ta > 180.0)
			theta(i) = 360.0 - ta;
		else
			theta(i) = ta;
	}

	return theta;

}   /* end of Conf_ZenithAngles() */




/*
 * sph_to_car() - matrix operand
 */
mat sph_to_car(mat phi, mat theta) {

	uword	i;

	mat sintheta = zeros<mat>(arma::size(theta));
	for (i = 0; i < theta.n_elem; i++)
		sintheta(i) = sind(theta(i));
	
	mat x = zeros<mat>(arma::size(phi));
	for (i = 0; i < phi.n_elem; i++)
		x(i) = sintheta(i) * cosd(phi(i));

	mat y = zeros<mat>(arma::size(phi));
	for (i = 0; i < phi.n_elem; i++)
		y(i) = sintheta(i) * sind(phi(i));

	mat z = zeros<mat>(arma::size(theta));
	for (i = 0; i < theta.n_elem; i++)
		z(i) = cosd(theta(i));

	mat out = zeros<mat>(phi.n_rows, 3);
	out = join_rows(x, y);
	out = join_rows(out, z);

	sintheta.reset();
	x.reset();
	y.reset();
	z.reset();

	return out;

}   /* end of sph_to_car() */




/*
 * sph_to_carD() - double typed scalar operand
 */
rowvec sph_to_carD(double phi, double theta) {

	double sintheta, x, y, z;
	rowvec out;
	
	sintheta = sind(theta);
	x = sintheta * cosd(phi);
	y = sintheta * sind(phi);
	z = cosd(theta);

	out << x << y << z << endr;

	return out;

}   /* end of sph_to_car() */




/*
 * BS_anttena_patterns() - based on antenna pattern inTR36.873
 */
mat BS_anttena_patterns(mat ZoDs, mat AoDs) {

	double SLA_V	 = 30.0;
	double theta_3dB = 65.0;
	double phi_3dB	 = 65.0;
	double A_max	 = 30.0;
	double G_max	 = 8.0;

	mat A_dB_theta(arma::size(ZoDs));
	mat A_dB_phi(arma::size(AoDs));
	mat A_dB(arma::size(A_dB_theta));
	mat A_pattern(arma::size(A_dB));

	uword i;

	for (i = 0; i < ZoDs.n_elem; i++)
		A_dB_theta(i) = -min(12.0 * pow(((ZoDs(i) - 90.0) / theta_3dB), 2.0), SLA_V);

	for (i = 0; i < AoDs.n_elem; i++)
		A_dB_phi(i) = -min(12.0 * pow((AoDs(i) / phi_3dB), 2.0), A_max);


	for (i = 0; i < A_dB_theta.n_elem; i++)
		A_dB(i) = -min(-(A_dB_theta(i) + A_dB_phi(i)), A_max);

	for (i = 0; i < A_dB.n_elem; i++)
		A_pattern(i) = pow(10.0, ((G_max + A_dB(i)) / 10.0));

	A_dB_theta.reset();
	A_dB_phi.reset();
	A_dB.reset();

	return A_pattern;

}   /* end of BS_anttena_patterns() */




mat UE_antenna_patterns(mat ZoAs, mat AoAs) {

	double SLA_V = 25.0;
	double theta_3dB = 90.0;
	double phi_3dB = 90.0;
	double A_max = 25.0;
	double G_max = 5.0;

	mat A_dB_theta(arma::size(ZoAs));
	mat A_dB_phi(arma::size(AoAs));
	mat A_dB(arma::size(A_dB_theta));
	mat A_pattern(arma::size(A_dB));

	uword i;

	for (i = 0; i < ZoAs.size(); i++)
		A_dB_theta(i) = -min(12.0 * pow(((ZoAs(i) - 90.0) / theta_3dB), 2.0), SLA_V);

	for (i = 0; i < AoAs.size(); i++)
		A_dB_phi(i) = -min(12.0 * pow((AoAs(i) / phi_3dB), 2.0), A_max);


	for (i = 0; i < A_dB_theta.size(); i++)
		A_dB(i) = -min(-(A_dB_theta(i) + A_dB_phi(i)), A_max);

	for (i = 0; i < A_dB.n_elem; i++)
		A_pattern(i) = sqrt(pow(pow(10.0, (G_max / 10.0) * 10.0), (A_dB(i) / 10.0)));

	A_dB_theta.reset();
	A_dB_phi.reset();
	A_dB.reset();

	return A_pattern;

}   /* end of UE_antenna_patterns() */



/* end of ChannelOutput_CDL.cpp */