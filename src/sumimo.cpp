/*
 * sumimo.cpp - misc function for Download SUMIMO simulator.
 *
 *			Programmed by PJW
 *
 * Copyright(c) 2016-2018 KIAST. All rights reserved.
 */

#include <cmath>
#include <armadillo>
#include "sumimo.h"



/*
* func_pseudo_sequence_generation() -  Generate pseudo random seqeunce according to
*										TS. 36.211, V13.3.0, 7.2
*										TS.38.211, V1.2.0, 5.2.1 (same)
*										c_init: initialization value
*										sequence_length: length of generated pseudo sequence
*/
void func_pseudo_sequence_generation(int c_init, double sequence_length, mat &c)
{
	uword			N_c = 1600;
	uword			l = 1;		// c_init.n_elem = 1
	uword			n;
	imat			x_1, x_2;
	vector<char>	temp;
	uword			k;

	n = N_c + (uword)sequence_length;
	x_1 = join_horiz(ones<imat>(l, 1), zeros<imat>(l, 30));
	x_1 = join_horiz(x_1, zeros<imat>(l, n - 31));
	x_2 = zeros<imat>(l, n);

	while (c_init > 1) {
		if (c_init % 2 == 1)
			temp.push_back(1);
		else
			temp.push_back(0);
		c_init /= 2;
	}
	if (c_init == 1)
		temp.push_back(1);


	k = temp.size();
	for (uword idx = 0; idx < k; idx++)
		x_2(0, idx) = (int)temp[idx];

	for (uword idx = 0; idx < n - 31; idx++) {
		for (uword i = 0; i < x_1.n_rows; i++) {
			x_1(i, idx + 31) = (x_1(i, idx + 3) + x_1(i, idx)) % 2;
			x_2(i, idx + 31) = (x_2(i, idx + 3) + x_2(i, idx + 2) + x_2(i, idx + 1) + x_1(i, idx)) % 2;
		}
	}

	c = zeros<mat>(l, n - (N_c));
	for (uword ridx = 0; ridx < x_1.n_rows; ridx++)
		for (uword cidx = 0, xidx = N_c; cidx < n - N_c; cidx++, xidx++)
			c(ridx, cidx) = (x_1(ridx, xidx) + x_2(ridx, xidx)) % 2;

	/* memory deallocations */
	x_1.reset();
	x_2.reset();
	temp.clear();


}   /* func_pseudo_sequence_generation() */





/*
 * mod() - integer modulo function.
 */
int mod(int a, int b) {

	int x;

	if (a < b)
		return a;
	else
		x = a / b;

	return a - b * x;

}	/* end of mod() */




/*
 * imod() - integer type MATLAB mod
 */
int imod(int a, int b) {

	int x;

	x = mod(a, b);

	if (x < 0)
		x = a + b;

	return x;

}	/* end of imod() */




/*
 * dmod() - double type MATLAB mod
 */
double dmod(double a, double b) {

	double x;

	x = fmod(a, b);

	if (x < 0.0)
		x = a + b;

	return x;

}	/* end of dmod() */




/*
 * matlnot() - element wide matrix logical not. Input matrix can be a 2-D matrix of the arma::umat type.
 */
umat matlnot(umat &X) {

	umat Y(arma::size(X));

	for (uword i = 0; i < X.n_elem; i++)
		if (X(i) == 0)
			Y(i) = 1;
		else
			Y(i) = 0;

	return Y;

}	/* end of matlnot() */




/*
 * fix() - Round Toward Zero with double type scalar value
 */
double fix(double x) {

	if (x > 0)
		return floor(x);
	else if (x == 0)
		return x;
	else
		return ceil(x);

}   /* end of fix() */



/* end of sumimo.cpp file */