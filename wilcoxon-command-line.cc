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
#include <vector>
#include <algorithm>
#include <fstream>

extern "C" {
	#include "wilcoxon.h"
}

/*
wilcoxon-command-line
Command-line utitlity to calculate Wilcoxon p-value

wilcoxon-command-line [switches] < input 

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

The first two numerics in cycle are to be integers, the other are doubles.
Stops on end-of-stream or if or the first given length (not tag) is zero.
	
*/

void help()
{
std::cout<<std::endl<<
"wilcoxon-command-line"<<std::endl<<
"Command-line utitlity to calculate Wilcoxon p-value"<<std::endl<<
""<<std::endl<<
"wilcoxon-command-line [switches] < input "<<std::endl<<
""<<std::endl<<
"If --taged switch is not given, the input is:"<<std::endl<<
"N_1 M_1 x_1_1.....x_N_1 y_1_1......y_M_1 "<<std::endl<<
"N_2 M_2 x_1_2.....x_N_2 y_1_2......y_M_2"<<std::endl<<
"......."<<std::endl<<
"in cycle."<<std::endl<<
""<<std::endl<<
"If --tagged (-t) switch is given, the input is:"<<std::endl<<
"tag1 N_1 M_1 x_1_1.....x_N_1 y_1_1......y_M_1 "<<std::endl<<
"tag2 N_2 M_2 x_1_2.....x_N_2 y_1_2......y_M_2"<<std::endl<<
"......."<<std::endl<<
"in cycle. The tag has no spaces."<<std::endl<<
""<<std::endl<<
"Each group N,M and M+N values after it is one test (one p-value)"<<std::endl<<
""<<std::endl<<
"the hypothesis parametesr is (only one is allowed in one run): "<<std::endl<<
"-w (--two, --two-tail) for two-tail (default)"<<std::endl<<
"-l (--lower, --lower-tail) for alternative med(x) >= med(y) (lower tail)"<<std::endl<<
"-u (--upper, --upper-tail) for alternative med(x) <=med(y)  (upper tail)"<<std::endl<<
""<<std::endl<<
"The first two are integers, the other are doubles."<<std::endl<<
"The first two numerics in cycle are to be integers, the other are doubles."<<std::endl<<
"Stops on end-of-stream or if or the first given length (not tag) is zero."<<std::endl;
}

int main(int argc,char ** argv)
{
	unsigned long tst0_len,tst1_len,inv,counter=0;
	double *tester0; //[gene_sets_number]
	double *tester1; //[gene_sets_number]
	double *frq; //[(gene_sets_number+1)/2) * ((gene_sets_number+1)/2)+1]
	double *work; //[gene_sets_number]
	//we will test+alloc frq and size inside p_value code

	std::istream * istream_ptr=&std::cin;
	
	bool tagged=false;
	hypothesis hyp=twoTail;

	unsigned int arguments_reflected=0;
	bool hyp_given=false;
	std::vector <std::string> args;
	for (int argn=1;argn<argc;argn++)
		args.push_back(argv[argn]);
	// command-line switches
	if
	(
		std::find(args.begin(),args.end(),
						"-?")
		!=args.end()
		||
		std::find(args.begin(),args.end(),
						"-h")
		!=args.end()
		||
		std::find(args.begin(),args.end(),
				"--help")
		!=args.end()
  )
	{
		help();
		return 0;
	}

	if
	(
		find(args.begin(),args.end(),
				"--tagged")
		!=args.end()
		||
		std::find(args.begin(),args.end(),
				"-t")
		!=args.end()
  )
	{
		tagged=true;
		arguments_reflected++;
	}

	if
	(
		find(args.begin(),args.end(),
				"-l")
		!=args.end()
		||
		std::find(args.begin(),args.end(),
				"--lower")
		!=args.end()
		||
		std::find(args.begin(),args.end(),
				"--lower-tail")
		!=args.end()
  )
	{
		if(hyp_given)
		{
			fprintf(stderr,"Hypothesis is given twice. Ask wilcoxon-command-line --help.\n");
			exit(-10);
		}
		hyp=lowerTail; 
		hyp_given=true;
		arguments_reflected++;
	}

	if
	(
		find(args.begin(),args.end(),
				"-u")
		!=args.end()
		||
		std::find(args.begin(),args.end(),
				"--upper")
		!=args.end()
		||
		std::find(args.begin(),args.end(),
				"--upper-tail")
		!=args.end()
  )
	{
		if(hyp_given)
		{
			fprintf(stderr,"Hypothesis is given twice. Ask wilcoxon-command-line --help.\n");
			exit(-10);
		}
		hyp=upperTail; 
		hyp_given=true;
		arguments_reflected++;
	}
	if
	(
		find(args.begin(),args.end(),
				"-w")
		!=args.end()
		||
		std::find(args.begin(),args.end(),
				"--two")
		!=args.end()
		||
		std::find(args.begin(),args.end(),
				"--two-tail")
		!=args.end()
  )
	{
		if(hyp_given)
		{
			fprintf(stderr,"Hypothesis is given twice. Ask wilcoxon-command-line --help.\n");
			exit(-10);
		}
		hyp=twoTail; 
		hyp_given=true;
		arguments_reflected++;
	}

	if (args.size()!=arguments_reflected)
	{
		bool lastarg_filename_ok=false;
		if (arguments_reflected+1==args.size())
		{
			std::string last_arg=args[arguments_reflected]; //last
			if ('-'!=last_arg.c_str()[0])
			{
				istream_ptr=new std::ifstream(last_arg.c_str());
				std::cout<<last_arg<<std::endl;
				if (*istream_ptr)
					lastarg_filename_ok=true;
				else
				{
					std::cerr<<"File "<<last_arg<<" is not readable.\n";
					exit(-20);
				}
			}
		};
		if (!lastarg_filename_ok)
		{
			std::cerr<<"It seemes to me that some arguments are unknown.  Ask wilcoxon-command-line --help.\n";
			exit(-10);
		}
	}


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

	for(std::istream_iterator<std::string> input(*istream_ptr);input!=eos;input++)
	{
		char * ptr; 
		const char * str;
		std::string tag;
		if (tagged)
		{
			counter++;
			tag=*input;
			if (eos==++input)
			{
				std::cerr<<"Unexpected eos after word "<<counter<<" (1-based). We supposed to have population of list 0 after tag."<<std::endl;
				exit(-50);
			}
		}
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
		if(tagged) std::cout<<tag<<" ";
		std::cout<<wilcoxon_p_value(inv, tst0_len, tst1_len, hyp, frq, work)
			<<std::endl;
	}
	delete istream_ptr;
}

