/***************************************************************\
  APSampler. Looking for complex genetic interaction patterns 
 by a Metropolis-Hastings MCMC project. (c) A. Favorov 1999-2010
    $Id$
\***************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern "C" {
	#include "wilcoxon.h"
}

#undef DEBUG_PRINT
//#define DEBUG_PRINT

int main(int argc,char ** argv)
{
	unsigned long tst0_len,tst1_len,inv,counter,expectation;
	double *tester0, //[gene_sets_number]
	double *tester1, //[gene_sets_number]
	double *frq, //[(gene_sets_number+1)/2) * ((gene_sets_number+1)/2)+1]
	double *work //[gene_sets_number]

	tst1_len=1;
	tst2_len=10;

	double ar1[1]={9};
	double ar2[10]={1,2,3,4,5,6,7,8,9,10};

	tester1=ar1;
	tester2=ar2
	
	frq_size=11;

	frq=(double*)calloc(frq_size,sizeof(double));
	assert(this->frq);
	if (!this->frq)
	{
		fprintf(stderr,"Allocation for sampler::frq failed.\n");
		exit(-500);
	}
	this->work=(double*)calloc(work_size,sizeof(double));
	assert(this->work);
	if (!this->work)
	{
		fprintf(stderr,"Allocation for sampler::work failed.\n");
		exit(-500);
	}

	inv=inversions(tst0_len, tester0, tst1_len, tester1);	
	wilcoxon_likelihood=
							wilcoxon_z_likelihood ( inv, tst0_len,tst1_len, frq,work);
}

extern unsigned int use_combinatorial_likelihood;
extern unsigned int use_exponential_likelihood;


unsigned long get_class_id (int * incmatrix_line, long patterns_number)
{
	int i;
	unsigned long result = 0;
	unsigned long mask=0x1;
	for (i=0;i<patterns_number;i++)
	{
		if (incmatrix_line[i]==1) 
		{
			result |= (mask << i);
			continue;
		}
		if (incmatrix_line[i]==0) continue;
		//if we are here, break.
		//it means that the patient's incmatrix line contains "-1"
		result=mask<<patterns_number;
		//put it to the absent class
		break;
	}

	return result;
}

//	p0=(this->zero_prior>0?this->zero_prior:
//			1./(this->dictionary[this->locus_offered_to_change-1]));
double patternset_zero_test 
			(
			 	double zero_hypothesis_prior,
				double uniformative_pattern_prior,
				unsigned int gene_sets_number, 
				unsigned int patterns_number, 
				double * levels, 
				int * incmatrix, 
				double * positive_hypotheses_posteriors,
				double * negative_hypotheses_posteriors,
				double *tester0, //[gene_sets_number]
				double *tester1, //[gene_sets_number]
				unsigned long *cls_ids, //[gene_sets_number]	
				double *frq, //[(gene_sets_number+1)/2) * ((gene_sets_number+1)/2)+1]
				double *work //[gene_sets_number]
				//the last five are work buffers
			)
/*
 * all double* are 0-based vectors, levels has gene_sets_number size,
 * roles and significances has patterns_number size.
 * incmatrix is 0-based-indexed matrix of 
 * gene_sets_number*patterns_number size
 * 
 * incmatrix[set#*patterns_number+pattern#]=1 means pattern pattern# 
 * is present in set no. set#, 0 value means it is absent there;
 * -1 means that we cannot test the pattern presence in the given genome
 * (the match is incorrect). So, we ascribe the case to an "absent class",
 * which clsid==2**patterns_number.
 *
 * we want to get the posteriors of zero-hypothesis about influence of 
 * every pattern on levels.
 * 
 * so, we split all the set to classes based on incidence matrix
 * a class for every different incmatrix string.
 * 
 * For every nonempty class pair we count Wilcoxon likelihood P(data | HO) as
 * 
 *       W(m,n)(i)
 *     ------------
 *        /m+n\
 *        \ n /      , where m and n are classes popuations and i is
 * 
 * the count of pairwise interclass comparisons which give y<x,
 *        
 *
 * We ascribe an exponential lokehood for both 
 *  
 * 
 * Using this, we can count 
 *
 *  
 * 
 *                                   P( data | H+ ) * P(H+) 
 * P( H+ | data ) = ----------------------------------------------------
 *   						     P(data|H+)*P(H+)+P(data|H-)*(P(H-))+P(data|H0)*P(H0)
 *
 * , etc.
 *
 * P(data|H+) = product ( P(i_e | H+) )
 *                 e
 *                                                   
 * P(data|H-) = product ( P(i_e | H-) )
 *                 e
 *
 * 
 * P(data|H-) = product ( P(i_e | H-) )   if (!use_combinatorial_likelihood)
 *                 e
 *
 *
 * P(data|H0) = product ( P(i_e | H-)+P(i_e | H0))+P(i_e | H+) ) 
 *                 e
 *
 *             -P(data | H+) - P(data | H-) if (use_combinatorial_likelihood)
 *
 *
 *
 *             
 * the closer is P(H0|data) posteior to 0, the more reliable is the pattern.
 * 
 * the returned value is multiplication of all (1-zero_posterior[pattern]),
 * which is an estimation for P(not h0|data) for all set.
 *
 * we allocate buffers outside the routine because we do not need
 * to clean them every time 
 * the size of work must be not less than (g*h+2)/(min(g,h)+2),
 * so we can say it is less tham gene_sets_number.
 */
{
	//our plan:
	//scan over patterns (i);
	//
	//for (pair_no = 0 .. 2**(p#-1) - 1
	//create class_p_id and class_a_id
	//put all points with class_a_id to test_0
	//put all points with class_p_id to test_1
	//count inversions 
	//count pair_sign
	//count wilcoxon likelihood
	//count flat likelihood
	//next pair_no
	//count h0-posterior
	//next pattern
	//count patternset posterior

	double patternset_z_h_rej_posterior;
	double wilcoxon_likelihood, pos_likelihood, neg_likelihood, 
		zero_hypothesis_posterior,posterior_norm;
	double H,L,w_0,w_mid;
	double power_2;
	double pos_likelihoods_product,
		neg_likelihoods_product,all_likelihoods_product,
		wilcoxon_likelihoods_product;
	unsigned long tst0_len,tst1_len,inv,min_class;
	unsigned long pair_no, max_pair_no, 
						pair_mask_low, pair_mask_high, class_0_id,class_1_id;
	unsigned long i,j;
	unsigned long left_cutoff,right_cutoff;
	unsigned int at_least_one_nonempty_pair;
	
	for (j=0;j<gene_sets_number;j++)
		cls_ids[j]=get_class_id(incmatrix+j*patterns_number,patterns_number);
#ifdef DEBUG_PRINT
	printf("CLS_IDS:\n");
	for (j=0;j<gene_sets_number;j++)
		printf("\t\t%lu : %lu  \n",j,cls_ids[j]);fflush(stdout);
	printf("\n");
	printf("[Combinatorial likelihood computation]=%i\n",use_combinatorial_likelihood);
	printf("[Exponential likelihood computation]=%i\n",use_exponential_likelihood);
#endif
	max_pair_no=(0x1<<(patterns_number-1)) - 1;
	patternset_z_h_rej_posterior=1.;
	for (i=0;i<patterns_number;i++)
	{
		tst0_len=tst1_len=0;
		for (j=0;j<gene_sets_number;j++)
		{
			if (*(incmatrix+j*patterns_number+i)==1)
				tst0_len++;
			else 
				tst1_len++;
		} 
		//we are counting pattern occurences in gene sets
		//to test if the pattern exist in some genomes and absent in other
		//if not, we omit it by prior(h0)=1
		if (tst0_len&&tst1_len)
		{
			pos_likelihoods_product=neg_likelihoods_product=
					all_likelihoods_product=wilcoxon_likelihoods_product=1.;
			at_least_one_nonempty_pair=0;
			for(pair_no=0;pair_no<=max_pair_no;pair_no++)
			{
				/*
				 *  pair_no : b01011   ; i=2 
				 *  we want class_0_id   =    b010011,
				 *                                *
				 *          class_1_id   =    b010111
				 *                           high | low
				 *                   
				 *  pair_mask_high = (b01011 >> 2) << 3 = (b010) << 3 = b010000
				 *  pair_mask_low = b01011 ^ (b010000>>1) = b01011 ^ b1000 = b11
				 *  We get what we want, now common part of the two
				 *  class ids is pair_mask_high|pair_mask_low
				 * 
				 */
				pair_mask_high=(pair_no>>i)<<(i+1);
				pair_mask_low=pair_no^(pair_mask_high>>1);
				class_0_id=pair_mask_high|pair_mask_low;
				class_1_id=class_0_id|(0x1<<i);
#ifdef DEBUG_PRINT
				printf("pair_no=%lu hi=%lu lo=%lu clsid0=%lu clsid2=%lu ",pair_no,pair_mask_high,pair_mask_low,class_0_id,class_1_id);
#endif
				tst0_len=tst1_len=0;
				for(j=0;j<gene_sets_number;j++)
				{
					if (cls_ids[j]==class_0_id)
						tester0[tst0_len++]=levels[j];
					if (cls_ids[j]==class_1_id)
						tester1[tst1_len++]=levels[j];
				}
#ifdef DEBUG_PRINT
				printf("len0=%lu ,len1=%lu ", tst0_len, tst1_len);
#endif
				if(tst0_len&&tst1_len)
				{
					at_least_one_nonempty_pair=1;
					inv=inversions(tst0_len, tester0, tst1_len, tester1);	
					wilcoxon_likelihood=
							wilcoxon_z_likelihood ( inv, tst0_len,tst1_len, frq,work);
					if(!use_exponential_likelihood)
					{
						left_cutoff=left_signed_cutoff (tst0_len,tst1_len);
						//the most-left number of inversions, which has Wilcoxon
						//probability that is not less than its flat one.
						//L lower wing for negetive pattern; 
						//R is right; if min_class==1,
						//they are equal and the distr is flat.
						right_cutoff=tst0_len*tst1_len-left_cutoff;
						min_class=tst1_len>tst0_len?tst0_len:tst1_len;
						power_2=pow(2,0.-(double)min_class);
						H=2*(1.-power_2)*((double)1./((double)tst0_len*(double)tst1_len + 1));
						L=2*power_2*((double)1./((double)tst0_len*(double)tst1_len + 1));
						//if it is flat, it is 1/max inv#
						//for every populator of min_class, we lower it twice.
						if (inv<=left_cutoff) 
						{
							pos_likelihood=H;
							neg_likelihood=L;
						}
						else if (inv>=right_cutoff)
						{
							pos_likelihood=L;
							neg_likelihood=H;
						}
						else //inv in ]left_cutoff,right_cutoff[
						{
							pos_likelihood=
									(
										L*((double)inv-(double)left_cutoff) +
										H*((double)right_cutoff-(double)inv)
									) /((double)right_cutoff-(double)left_cutoff);
							neg_likelihood=
									(
										H*((double)inv-(double)left_cutoff) +
										L*((double)right_cutoff-(double)inv)
									) /((double)right_cutoff-(double)left_cutoff);
						}
					}
					else //use_exponential_likelihood
					{
						w_0=
								wilcoxon_z_likelihood 
								(
									0,
									tst0_len,
									tst1_len,
									frq,
									work
								);
						w_mid=
								wilcoxon_z_likelihood 
								( 
									tst0_len*tst1_len/2, 
									tst0_len,
									tst1_len, 
									frq,
									work
								);
						pos_likelihood=
								exponential_likelihood
								(
									tst1_len,
									tst0_len,
									inv,
									&neg_likelihood,
									w_0,
									w_mid
								);
  				}
 
 					//the most probable is no inversions
				}
				else //one of the casses is empty
					//actually the code must work with 0,
					//but it do not make sense to test it.
				{
					inv=0;
					pos_likelihood=neg_likelihood=wilcoxon_likelihood=1.;
				}
				pos_likelihoods_product*=pos_likelihood;
				neg_likelihoods_product*=neg_likelihood;
				if (use_combinatorial_likelihood)
					all_likelihoods_product*=
							(wilcoxon_likelihood+pos_likelihood+neg_likelihood);
				else
					wilcoxon_likelihoods_product*=
							wilcoxon_likelihood;
							
#ifdef DEBUG_PRINT
				printf(": inv=%lu : p0=%g , p-=%g , p+=%g\n",
						inv, wilcoxon_likelihood, neg_likelihood, pos_likelihood);
#endif
			}
			if (use_combinatorial_likelihood || at_least_one_nonempty_pair)
			{
				positive_hypotheses_posteriors[i]=
						pos_likelihoods_product*.5*(1.-zero_hypothesis_prior);
				negative_hypotheses_posteriors[i]=
						neg_likelihoods_product*.5*(1.-zero_hypothesis_prior);
				if (use_combinatorial_likelihood)
					zero_hypothesis_posterior=
							(all_likelihoods_product-(pos_likelihoods_product+
							 neg_likelihoods_product))*zero_hypothesis_prior;
				else 
					zero_hypothesis_posterior=
							wilcoxon_likelihoods_product*zero_hypothesis_prior;

				posterior_norm=positive_hypotheses_posteriors[i]+
						negative_hypotheses_posteriors[i]+
						zero_hypothesis_posterior;

				positive_hypotheses_posteriors[i]/=posterior_norm;
				negative_hypotheses_posteriors[i]/=posterior_norm;
				zero_hypothesis_posterior/=posterior_norm;

#ifdef DEBUG_PRINT
				printf("whole pattern[%lu] : zero=%g pos/neg=%g/%g\n",
						i,zero_hypothesis_posterior,positive_hypotheses_posteriors[i],
						negative_hypotheses_posteriors[i]
						);
#endif
			}
			else
			{
				// we are here, ir means that there was no such a class pair,
				// that tst0 and tst1 are both nonempty, i.e. the pattern is 
				// uninformative
				//
				positive_hypotheses_posteriors[i]=negative_hypotheses_posteriors[i]=
					(1.-uniformative_pattern_prior)/2.;
#ifdef DEBUG_PRINT
				printf("pattern is uninformative, pos=neg=%g\n",positive_hypotheses_posteriors[i]);
#endif
			}
		}
		else
		{
			// we are here, ir means that tst0 or tst1 of that pattern alone 
			// scan is empty, i.e. it is uninformative
			//
			positive_hypotheses_posteriors[i]=negative_hypotheses_posteriors[i]=
				(1.-uniformative_pattern_prior)/2.;
#ifdef DEBUG_PRINT
			printf("pattern is uninformative, pos=neg=%g\n",positive_hypotheses_posteriors[i]);
#endif
		}
		patternset_z_h_rej_posterior*=
				(positive_hypotheses_posteriors[i]+negative_hypotheses_posteriors[i]);
	}
#ifdef DEBUG_PRINT
				printf("Patternset :  notzero=%g\n",
						patternset_z_h_rej_posterior
						);
#endif
	return (patternset_z_h_rej_posterior);
}
