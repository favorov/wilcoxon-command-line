/***************************************************************\
wilcoxon-command-line (c) A. Favorov 2015	
$Id$
\***************************************************************/


#ifndef _GAUSS_H_
#define _GAUSS_H_

double cgauss(double x);

/*   Cumulative Gaussian. 


                   1            x           -t^2 
  cgauss(x) is ----------- * [integral] (exp ----- * dt)
               sqrt(2*pi)     -inf            2


   This routine takes a number between -infinity and infinity and returns
   the corresponding p-value of a cumulative Gaussian, the underlying
   Gaussian having zero mean and unit standard deviation.
*/
#endif
