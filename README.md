# wilcoxon-command-line
Command-line utitlity to calculate Wilcoxon p-value

The input is:
N_1 M_1 x_1_1.....x_N_1 y_1_1......y_M_1 
N_2 M_2 x_1_2.....x_N_2 y_1_2......y_M_2
.......
in cycle.
Each group N,M and M+N values after it is one test (one p-value)

the only parameter is: 
-t for two-tail (default)
-l for alternative med(x) >= med(y) (lower tail)
-u for alternative med(x) <=med(y)  (upper tail)

The first two are integers, the other are doubles.
If the stream stops before all the N+M values are read
or the first two values 
