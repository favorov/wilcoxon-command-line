/***************************************************************\
wilcoxon-command-line (c) A. Favorov 2015	
$Id$
\***************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <iomanip>
#include <iterator>
#include <string>

extern "C" {
	#include "wilcoxon.h"
}

/*
wilcoxon-command-line
Command-line utitlity to calculate Wilcoxon p-value

wilcoxon-command-line [swithes] < input 

If --taged switch is not given, the input is:
N_1 M_1 x_1_1.....x_N_1 y_1_1......y_M_1 
N_2 M_2 x_1_2.....x_N_2 y_1_2......y_M_2
.......
in cycle.

If --tagged (-t) switch is given, the input is:
tag1 N_1 M_1 x_1_1.....x_N_1 y_1_1......y_M_1 
tag2 N_2 M_2 x_1_2.....x_N_2 y_1_2......y_M_2
.......
in cycle. The tag has no spaces.

Each group N,M and M+N values after it is one test (one p-value)

the hypothesis parametesr is (only one is allowed in one run): 
-w (--two, --two-tail) for two-tail (default)
-l (--lower, --lower-tail) for alternative med(x) >= med(y) (lower tail)
-u (--upper, --upper-tail) for alternative med(x) <=med(y)  (upper tail)

The first two are integers, the other are doubles.
If the stream stops before all the N+M values are read
or the first two values 
	
*/
int main(int argc,char ** argv)
{
	unsigned long tst0_len,tst1_len,inv,counter=0;
	double *tester0; //[gene_sets_number]
	double *tester1; //[gene_sets_number]
	double *frq; //[(gene_sets_number+1)/2) * ((gene_sets_number+1)/2)+1]
	double *work; //[gene_sets_number]
	//we will test+alloc frq and size inside p_value code

	
	bool tagged;
	hypothesis hyp;













	frq=(double*)calloc(usual_frq_size,sizeof(double));
	assert(frq);
	if (!frq)
	{
		fprintf(stderr,"Allocation for main::frq failed.\n");
		exit(-500);
	}
	work=(double*)calloc(usual_work_size,sizeof(double));
	assert(work);
	if (!work)
	{
		fprintf(stderr,"Allocation for main::work failed.\n");
		exit(-500);
	}

 	std::istream_iterator<std::string> eos;              // end-of-stream iterator

	for(std::istream_iterator<std::string> input (std::cin);input!=eos;input++)
	{
		char * ptr; 
		const char * str;
		//std::cout<<*input<<std::endl;
		counter++;
		str=input->c_str();
		tst0_len=strtoul(str,&ptr,0);
		if (0 != str+input->length()-ptr)
		{
			std::cerr<<"Word "<<counter<<" (1-based) is not integer (population of list 0) as it was expected."<<std::endl;
			exit(-50);
		}
		
		if(0==tst0_len) break; //0 breaks

		tester0=(double*)calloc(tst0_len,sizeof(double));
		assert(tester0);
		if (!tester0)
		{
			std::cerr<<"Allocation for list 0  failed at word counter "<<counter<<" (1-based)."<<std::endl;
			exit(-500);
		}
		
		if (eos==++input)
		{
			std::cerr<<"Unexpected eos after word "<<counter<<" (1-based). We supposed to have population of list 1."<<std::endl;
			exit(-50);
		}
		counter++;
		str=input->c_str();
		tst1_len=strtoul(str,&ptr,0);
		if (0 != str+input->length()-ptr)
		{
			std::cerr<<"Word "<<counter<<" (1-based) is not integer (population of list 1) as it was expected."<<std::endl;
			exit(-50);
		}
		

		tester1=(double*)calloc(tst1_len,sizeof(double));
		assert(tester1);
		if (!tester1)
		{
			std::cerr<<"Allocation for list 1  failed at word counter "<<counter<<" (1-based)."<<std::endl;
			exit(-500);
		}
		
		for (unsigned int i=0;i<tst0_len;i++)
		{
			if (eos==++input)
			{
				std::cerr<<"Unexpected eos after word "<<counter<<" (1-based). We supposed to have element "<<i<<" (0-based) of list 0."<<std::endl;
				exit(-50);
			}
			counter++;
			str=input->c_str();
			tester0[i]=strtod(str,&ptr);
			if (0 != str+input->length()-ptr)
			{
				std::cerr<<"Word "<<counter<<" (1-based) is not double (supposed to be element "<<i<<" in list 0)."<<std::endl;
				exit(-50);
			}
		}

		for (unsigned int i=0;i<tst1_len;i++)
		{
			if (eos==++input)
			{
				std::cerr<<"Unexpected eos after word "<<counter<<" (1-based). We supposed to have element "<<i<<" (0-based) of list 1."<<std::endl;
				exit(-50);
			}
			counter++;
			str=input->c_str();
			tester1[i]=strtod(str,&ptr);
			if (0 != str+input->length()-ptr)
			{
				std::cerr<<"Word "<<counter<<" (1-based) is not double (supposed to be element "<<i<<" in list 1)."<<std::endl;
				exit(-50);
			}
		}

		inv=inversions(tst0_len, tester0, tst1_len, tester1);	
		std::cout<<wilcoxon_p_value( inv, tst0_len, tst1_len, twoTail, frq, work)
			<<std::endl;
	}

}

