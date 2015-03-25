/***************************************************************\
  APSampler. Looking for complex genetic interaction patterns 
 by a Metropolis-Hastings MCMC project. (c) A. Favorov 1999-2010
    $Id$
\***************************************************************/

#include "gauss.h"
#include <stdio.h>
#include <math.h>

double erfcc(double x);

/*
 Source reference:
 Foster, D.H. and Bischof, W.F. (1991) Thresholds from
 psychometric functions: superiority of bootstrap to
 incremental and probit variance estimators. Psychological
 Bulletin, 109, 152-159.
*/

double cgauss(double x)
//
//                   1            x           -t^2 
//  cgauss(x) is ----------- * [integral] (exp ----- * dt)
//               sqrt(2*pi)     -inf            2
//
/* Useful comments from CEPHES:
 *                             x
 *                             -
 *                   1        | |          2
 *    ndtr(x)  = ---------    |    exp( - t /2 ) dt
 *               sqrt(2pi)  | |
 *                           -
 *                          -inf.
 *
 *             =  ( 1 + erf(z) ) / 2
 *             =  erfc(z) / 2
 *
 * where z = x/sqrt(2). Computation is via the functions
 * erf and erfc.
 *
 *						erf.c
 *
 *	Error function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, erf();
 *
 * y = erf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * The integral is
 *
 *                           x 
 *                            -
 *                 2         | |          2
 *   erf(x)  =  --------     |    exp( - t  ) dt.
 *              sqrt(pi)   | |
 *                          -
 *                           0
 *
 * The magnitude of x is limited to 9.231948545 for DEC
 * arithmetic; 1 or -1 is returned outside this range.
 *
 * For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
 * erf(x) = 1 - erfc(x).
 *
 *
 *
 *
 *	Complementary error function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, erfc();
 *
 * y = erfc( x );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *  1 - erf(x) =
 *
 *                           inf. 
 *                             -
 *                  2         | |          2
 *   erfc(x)  =  --------     |    exp( - t  ) dt
 *               sqrt(pi)   | |
 *                           -
 *                            x
 *
 *
 * For small x, erfc(x) = 1 - erf(x); otherwise rational
 * approximations are computed.
 *
 */
{
/*
   Cumulative Gaussian.
   This routine takes a number between -infinity and infinity and returns
   the corresponding p-value of a cumulative Gaussian, the underlying
   Gaussian having zero mean and unit standard deviation.
*/

	double sqr2=sqrt(2.);
	if (x < 0.0) return((double)0.5*erfcc(-x/sqr2));
	else return((double)1.0 - (double)0.5*erfcc(x/sqr2));
}

double erfcc(double x)
{
/*
   Error function.
   Source: Cody, W. J. (1969). Rational Chebyshev approximations for the
   error function, Math. Comp., 23, 631-638. Described in Kennedy, W. J. &
   Gentle, J. E. (1980) Statistical Computing, New York: Dekker, pp. 90ff
// the last alternative is commented by A. Favorov in 2000 August
// because it does not work.
 
 */
	static double p1[4] = {
		2.4266795523053175e2,2.1979261618294152e1,6.9963834886191355,
		-3.5609843701815385e-2};
	static double q1[4] = {
		2.1505887586986120e2,9.1164905404514901e1,1.5082797630407787e1,
		1.0};
	static double p2[8] = {
		3.004592610201616005e2,4.519189537118729422e2,3.393208167343436870e2,
		1.529892850469404039e2,4.316222722205673530e1,7.211758250883093659,
		5.641955174789739711e-1,-1.368648573827167067e-7};
	static double q2[8] = {
		3.004592609569832933e2,7.909509253278980272e2,9.313540948506096211e2,
		6.389802644656311665e2,2.775854447439876434e2,7.700015293522947295e1,
		1.278272731962942351e1,1.0};
//	static double p3[5] = {
//		-2.99610707703542174e-3,-4.94730910623250734e-2,-2.26956593539686930e-1,
//		-2.78661308609647788e-1,-2.23192459734184686e-2};
//	static double q3[5] = {
//		1.06209230528467918e-2,1.91308926107829841e-1,1.05167510706793207,
//		1.98733201817135256,1.0};
	
	double z,zz,zpow,p,q,r;
	int j;
	
	z=fabs(x);
	if (z < 0.46875) {
		zz=z*z;
		zpow=zz;
		p=p1[0];
		q=q1[0];
		for (j=1;j<=3;j++) {
			p+=p1[j]*zpow;
			q+=q1[j]*zpow;
			zpow*=zz;
		}
		r=(double)1.-z*p/q;
	}
	else /* if (z < 4.0) */ {
		zpow=z;
		p=p2[0];
		q=q2[0];
		for (j=1;j<=7;j++) {
			p+=p2[j]*zpow;
			q+=q2[j]*zpow;
			zpow*=z;
		}
		r=exp((double)(-z*z))*p/q;
	}
/*	else {
		zz=1.0/(z*z);
		zpow=zz;
		p=p3[0];
		q=q3[0];
		for (j=1;j<=4;j++) {
			p+=p3[j]*zpow;
			q+=q3[j]*zpow;
			zpow*=zz;
		}
		r=exp((double)(-zz))/z*(0.564189583548+zz*p/q);
	} */
	return r;
}

