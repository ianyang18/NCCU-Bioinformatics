#Homework 2
##Goal
Sequence pairwise alignment by using dynamic programming.

##Files
###Data
* Score Table
	* pam100.txt
	* pam250.txt
* Test Files - *test_file* folder
* Result Files - *result_file* folder

###Code
* pro2_104753024.R

##Execution

**Previously, here are the arguments that you can execute.**

--input : read the file including two sequences, input file should follw the example format e.g. test.fasta.txt

--score (option) : to load score table into the program

--aln : global/local alignment choices

--gap_open (option) : additional gap penalty for the starting position, negative integer value

--gap_extend (option) : gap penalty, negative integer value

--output : output the sequence alignment result

**Instruction**

```
Rscript pro2_104753024.R --input <test_file> --score <score_table> --aln <global/local> --gap_open <negative_value> --gap_extend <negative_value> --output <result_file>
```


**Example (for the testing purpose)**

without --score, --gap option

```
Rscript pro2_104753024.R --input test_file/test1.txt --aln global --output result_file/result1_global.txt
```

```
Rscript pro2_104753024.R --input test_file/test1.txt --aln local --output result_file/result1_local.txt
```

```
Rscript pro2_104753024.R --input test_file/test2.txt --aln global --output result_file/result2_global.txt
```


with --gap option

```
Rscript pro2_104753024.R --input test_file/test2.txt --aln global --gap_open -4 --gap_extend -3 --output result_file/result2_gap_global.txt
```

```
Rscript pro2_104753024.R --input test_file/test2.txt --aln local --gap_open -4 --gap_extend -3 --output result_file/result2_gap_local.txt
```


##Reference
- [Global Alignment](http://www.csie.ntu.edu.tw/~kmchao/seq16fall/global_alignment.pdf)
- [Local Alignment](http://www.csie.ntu.edu.tw/~kmchao/seq16fall/local_alignment.pdf)
- [Various Scoring Schemes](http://www.csie.ntu.edu.tw/~kmchao/seq16fall/various_schemes.pdf)
- [An affine-gap-penalty example](http://www.csie.ntu.edu.tw/~kmchao/seq16fall/affine_gap_example.pdf)