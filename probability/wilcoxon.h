/***************************************************************\
  APSampler. Looking for complex genetic interaction patterns 
 by a Metropolis-Hastings MCMC project. (c) A. Favorov 1999-2010
    $Id$
\***************************************************************/

#ifndef _WILCOXON_H
#define _WILCOXON_H

typedef enum {two-tail=0,lower-tail=1,upper-tail=2} hypothesis; 

double wilcoxon_p_value
				(unsigned long i,unsigned long m, unsigned long n, hypothesis hyp,
				double * frqncy, double * work);

/*
 * returns p-value of i inversions in set of m x's and n y's.
 * frqncy and work are double arrays prellocated as
 * desribed in comment to udist.
 * frqncy[0..m*n]
 * work is not less than (m*n + 2) / 2) + min(m,n) + 1
 * As far as we are going to use the Gaussian for m+n>=30, m>3,n>3
 * the maximal *work is (225+2)/2+15+1=130 (let it be 150) (it is for m,n>3) or (max(m,n)*3+2)/2 + 3 + 1
 * the maximal frqncy is: 226 ... let it be 250 (it is for m,n>3) or max(m,n)*3 + 1
*/

double wilcoxon_p_value_exp
				(unsigned long i,unsigned long m, unsigned long n, hypothesis hyp);

/*
 *  if m+n>=20 ; m>3, n>3 we can approximate p 
 *                                  /          i - mn/2             \
 *  as dustributed as Gauss[0,1]   | ------------------------------  |
 *                                  \ [ 1/12 mn (m + n + 1) ] ^ 1/2 /
 * 
 *  It is Gaussian approximaton for previous case.
 *  We use 30 !!!!
 *
 *  MATEMATISHE STATISTIK von
 *  Dr B.I. Van Der Varden
 *  Springer-Verlag BEGLIN-GOTTINGEN-HEIDELBERG, 1957
 *  chapter XII. Section 63.
 */

double wilcoxon_z_likelihood
				(unsigned long i,unsigned long m, unsigned long n,
				double * frqncy, double * work);
/*
 * returns likelihood of i inversions in set of m x's and n y's.
 * frqncy and work are double arrays prellocated as
 * desribed in comment to udist.
 * frqncy[0..m*n]
 * work is not less than (m*n + 2) / 2) + min(m,n) + 1
 * As far as we are going to use the Gaussian for m+n>=30, m>3,n>3
 * the maximal *work is (225+2)/2+15+1=130 (let it be 150) (it is for m,n>3) or (max(m,n)*3+2)/2 + 3 + 1
 * the maximal frqncy is: 226 ... let it be 250 (it is for m,n>3) or max(m,n)*3 + 1
*/
double wilcoxon_z_likelihood_exp
				(unsigned long i,unsigned long m, unsigned long n);

/*
 *  if m+n>=20 ; m>3, n>3 we can approximate p 
 *                                  /          i - mn/2             \
 *  as dustributed as Gauss[0,1]   | ------------------------------  |
 *                                  \ [ 1/12 mn (m + n + 1) ] ^ 1/2 /
 * 
 *  It is Gaussian approximaton for previous case.
 *  We use 30 !!!!
 *
 *  MATEMATISHE STATISTIK von
 *  Dr B.I. Van Der Varden
 *  Springer-Verlag BEGLIN-GOTTINGEN-HEIDELBERG, 1957
 *  chapter XII. Section 63.
 */

				
double * udist(unsigned long m, unsigned long n,
		double * frqncy, double * work);
/*
 * It will generate all frequencies for U (Wilcoxon) distribution and put it to
 * frqncy[0..n*m] and return frqncy once more.
 * frqncy is a double which is to be allocated not less than m*n+1
 * work is intermediate array. It is to be not less than 
 * (m*n + 2) / 2) + min(m,n) + 1
 * 
 * 
 * We are preallocate both arrays to save time
 * 
 */

unsigned long left_signed_cutoff (unsigned long m, unsigned long n);

/*
 * Returns the most-left number of inversions, which has Wilcoxon
 * probability that is not less than its flat one.
 *
 *
 *            U (Wilcoxon)  
 *             *
 *     ____*_______*_____    Flat (1..m*n)
 *       *           *
 * --------+-------+-------------------------------
 *       here
 *        it
 *        is
 *     
 *
 * Returns 0 if both are flat (one of m,n is 1 or 0)
 */ 



unsigned long inversions 
					(unsigned int xlen, double * x, 
					 unsigned int ylen, double * y);
/*
 *  x and y are 0-based arrays of size xlen and ylen 
 *  !!! They will be permuted inside themselves after 
 *  this operation
 *  
 *    inv = sum ( y_j<x_i )
 *          i,j
 *       i=0..xlen
 *       j=0..ylen    
 * 
 *  For any equality containing u x's and v y's cluster we do: 
 *  
 *  inv+=uv/2 (we use intger division, if uv is odd, 
 *  one equality is saved, we count of all this equalities and finally
 *  shift the inversions count to the middle)
 */
#endif
