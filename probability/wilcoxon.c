/***************************************************************\
  APSampler. Looking for complex genetic interaction patterns 
 by a Metropolis-Hastings MCMC project. (c) A. Favorov 1999-2010
    $Id$
\***************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "wilcoxon.h"
#include "gauss.h"

/* AS 62 generates the frequencies for the Mann-Whitney U-statistic. */
/* Users are much more likely to need the distribution function. */
/* Code to return the distribution function has been added at the end */
/* of AS 62 by Alan Miller.   Remove the C's in column 1 to activate it. */
/*     ALGORITHM AS 62  APPL. STATIST. (1973) VOL.22, NO.2 
       The distribution of the Mann-Whitney U-statistic is generated for 
       the two given sample sizes */

double * udist(unsigned long m, unsigned long n,
		double * frqncy, double * work)
/*
 * It will generate all frequencies for U distribution and put it to
 * frqncy[0..n*m] and return frqncy once more.
 * frqncy is a double which is to be allocated not less than m*n+1
 * work is intermediate array. It is to be not less than 
 * (m*n + 2) / 2) + min(m,n) + 1
 * 
 * We are preallocate both arrays to save time
 */
{
	unsigned long i, j, k, l, minmn, maxmn, n1, in, mn1;
	double sum;
	minmn = m>n?n:m;
	if (minmn < 1) 
		return frqncy;
	mn1 = m * n + 1;

/* Set up results for 1st cycle (one pf sets=1) and return if MINMN = 1 */

	maxmn = m>n?m:n;
	n1 = maxmn + 1;
	for (i = 0; i < n1; i++) 
		frqncy[i] = 1.;
	for (; i < mn1; ++i) 
		frqncy[i] = 0;
	if (minmn == 1) 
		return frqncy;

/* Generate successively higher order distributions */

	work[0] = 0.;
	in = maxmn;
	for (i = 2; i <= minmn; i++) //cycling over one of sets sizes
	{
		work[i-1] = 0.;
		in += maxmn;
		n1 = in + 1;
		l = in / 2 + 1;
		k = i-1;

/* Generate complete distribution from outside inwards */

		for (j = 0; j < l; j++) 
		{
			++k;
			--n1;
			sum = frqncy[j] + work[j];
			frqncy[j] = sum;
			work[k] = sum - frqncy[n1];
			frqncy[n1] = sum;
		}
	}
	return frqncy;	
} /* udist */

double wilcoxon_z_likelihood
				(unsigned long i,unsigned long m, unsigned long n,
				double * frqncy, double * work)
/*
 * returns posterior of i inversions in set of m x's an y y's.
 * frqncy and work are double arrays prellocated as
 * desribed in comment to udist.
 * frqncy[0..m*n]
 * work is not less than (m*n + 2) / 2) + min(m,n) + 1
 */
{
	double sum=0;
	unsigned long k,l;
	if (i>m*n) return 0;
	if (m==0 || n==0) return 1;
	//if one class if empty, inv. number == 0 whatever
	if (m==1)
		return 1./(double)(n+1);
	if (n==1)
		return 1./(double)(m+1);
	if (m>3&&n>3&&m+n>=30) return wilcoxon_z_likelihood_exp(i,m,n);
	udist(m,n,frqncy,work);
	l=m*n;
	for (k=0;k<=l;k++) sum+= frqncy[k];
	return frqncy[i]/sum;
}

double wilcoxon_norm_equvalent
				(double i,unsigned long m, unsigned long n)
{
	return 
			(
			 	(i - (double)m * (double)n * 0.5) / 
				sqrt
				(
					((double)1./(double)12.) * (double)m * (double)n * 
					((double)m + (double)n + 1.)
				)
			);
}

double wilcoxon_z_likelihood_exp
				(unsigned long i,unsigned long m, unsigned long n)
/*
 *  g>3,h>3
 * 
 *  if g+h>=20 we can approximate p 
 *                                  /          i - mn/2             \
 *  as dustributed as Gauss[0,1]   | ------------------------------  |
 *                                  \ [ 1/12 mn (m + n + 1) ] ^ 1/2 /
 * 
 *  It is Gaussian approximaton for previous case.
 * 
 *  MATEMATISHE STATISTIK von
 *  Dr B.I. Van Der Varden
 *  Springer-Verlag BEGLIN-GOTTINGEN-HEIDELBERG, 1957
 *  chapter XII. Section 63.
 */
{
	double doublei,c1,c2;
//	if (i<0) return 0;
	if (i>m*n) return 0;
//	if (i==0) return cgauss(wilcoxon_norm_equvalent(.5,m,n));
//	if (i==m*n) return (double)1.-cgauss(wilcoxon_norm_equvalent(m*n-.5,m,n));

// doublei=(double)i;
//So, we've replaced the line:
	doublei=(double)(2*i>n*m?m*n-i:i);
//Linux/Athlon gives a buggy result here for i close to mn.
	c1=	cgauss(wilcoxon_norm_equvalent(doublei+.5,m,n));
	c2= cgauss(wilcoxon_norm_equvalent(doublei-.5,m,n));
//	printf ("-> %30.20f   %30.20f\n",c1,c2);
	return (c1-c2); 
}


unsigned long left_signed_cutoff (unsigned long m, unsigned long n)

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
 */
{
	unsigned long current_point,left_point=0,right_point=m*n/2;
	double flat_posterior=1./(n*m+1),current_wilcoxon;
	if (m==1 || n ==1 || m==0||n==0) return 0;
	while ((right_point-left_point)>1)
	{
		current_point=(left_point+right_point)/2;
//		printf ("Currpo: %lu ***\n", current_point);
		current_wilcoxon=wilcoxon_z_likelihood_exp(current_point,m,n);
		if (current_wilcoxon>=flat_posterior)
			right_point=current_point;
		else left_point=current_point;  //we use the fact that Wilcoxon left
																		//half in monotonic
	}
	//if we are here, than right_point-left_point=1 or 0;
	//wilcoxon(right_point) >= flat_posterior;
	return right_point;	
}

unsigned long inversions
					(unsigned int xlen, double * x, 
					 unsigned int ylen, double * y)
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
{
	unsigned long i,j,inv=0;
	double x_value;
	unsigned long equal_x_count;
	unsigned int *x_used; //boolean vector of xlen
	unsigned long equal_y_count;
	unsigned long odd_equality_clusters=0;
	//initilisation
	if (!xlen || !ylen) return 0;
	x_used=(unsigned int *)calloc(xlen,sizeof(unsigned int));
	assert(x_used);
	//it is not a time leak - we are to clean it whatever
	//counting...
	for(i=0;i<xlen;i++)
	{
		if (x_used[i]) continue; 	//we have already used it
		x_used[i]=1;             	//we are doing it now
		x_value=x[i];							//for more speed
		equal_x_count=1;
		equal_y_count=0;
		for (j=i+1;j<xlen;j++)
			if(x[j]==x_value)
			{
				equal_x_count++;
				x_used[j]=1;         //and these, too
			};//the loop was to find all x-cluster
		for (j=0;j<ylen;j++)
		{
			if(y[j]<x_value)
			{
				inv+=equal_x_count;
				continue;
			}
			if(y[j]>x_value)
				continue;
			equal_y_count++;    //they appeared to be equal
		}
		if (equal_y_count)
		{
			inv+=(equal_y_count*equal_x_count)/2;
			odd_equality_clusters+=(equal_y_count*equal_x_count)%2;
		}
	}
	//finish
	free(x_used);
	if(inv*2<xlen*ylen)
	{
		inv+=odd_equality_clusters;
		if (inv*2>xlen*ylen) inv=(xlen*ylen)/2;
	}
	//if inv<mean(inv), we increase it
	//by odd_equality_clusters 
	//so, the odd cluster is a bit closer to mean than
	//just its half of invrsions.
	//if we jump over mean, inv is mean.
	return inv;
}

