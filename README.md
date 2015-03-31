# wilcoxon-command-line
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

The first two numerics in cycle are to be integers, the other are doubles.
Stops on end-of-stream or if or the first given length (not tag) is zero. 
