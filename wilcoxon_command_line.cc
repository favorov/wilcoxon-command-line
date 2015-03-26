/***************************************************************\
  APSampler. Looking for complex genetic interaction patterns 
 by a Metropolis-Hastings MCMC project. (c) A. Favorov 1999-2010
    $Id$
\***************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <iomanip>


extern "C" {
	#include "wilcoxon.h"
}

int main(int argc,char ** argv)
{
	unsigned long tst0_len,tst1_len,inv,counter,expectation;
	double *tester0; //[gene_sets_number]
	double *tester1; //[gene_sets_number]
	double *frq; //[(gene_sets_number+1)/2) * ((gene_sets_number+1)/2)+1]
	double *work; //[gene_sets_number]

	tst0_len=3;
	tst1_len=10;

	double ar0[3]={9,11,12};
	double ar1[10]={1,2,3,4,5,6,7,8,9,10};

	tester0=ar0;
	tester1=ar1;
	
	unsigned long frq_size=31;
	unsigned long work_size=100;

	frq=(double*)calloc(frq_size,sizeof(double));
	assert(frq);
	if (!frq)
	{
		fprintf(stderr,"Allocation for sampler::frq failed.\n");
		exit(-500);
	}
	work=(double*)calloc(work_size,sizeof(double));
	assert(work);
	if (!work)
	{
		fprintf(stderr,"Allocation for sampler::work failed.\n");
		exit(-500);
	}

	inv=inversions(tst0_len, tester0, tst1_len, tester1);	
	std::cout<<wilcoxon_p_value( inv, tst0_len, tst1_len, twoTail, frq, work)
		<<std::endl;
}

