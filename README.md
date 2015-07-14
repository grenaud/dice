# Readme for DICE

# Prerequisites

Python libraries:
- dadi
- numpy
- scipy

R libraries:
- bbmle

# Compiling and installing

cd libgab/

make

cd ..

cd src/

make

cd ..

# 2-Pop method: input data format

The input data for the 2-pop method is a file that should have at least four columns. Each row in the file denotes a particular configuration of ancestral reads, derived read and anchor/contaminant allele frequencies. An example row would be:

3	5	0.48	8

This means there are 8 sites in the genome where the ancient genome has 3 ancestral reads and 5 derived reads, and where the anchor/contaminant panel has a derived allele frequency of 0.48. If the file only has 4 columns, the program will assume the contaminant panel is the same as the anchor panel. The file should also have a header. We provide an example test file in the testData folder, under the name "test_twopop_4col.txt".

If the file has 5 columns, the program will assume that the contaminant and anchor panels are distinct. In this case, the third column refers to the anchor panel, while the fourth column refers to the contaminant panel. For example:

3	5	0.48	0.96	8

This row means there are 8 sites where the ancient genome has 3 ancestral reads and 5 derived reads, and where the anchor panel has a derived allele frequency of 0.48 and where the contaminant panel has a derived allele frequency of 0.96. We provide an example test file in the testData folder, under the name "test_twopop_5col.txt".

NOTE: when comparing multiple candidate contaminant populations, it is best to use only sites that are segregating in all populations that are tested as contaminants, so that the likelihood is composed of the same number of sites in each case.

# 2-Pop method: usage

To run type:

./dice -2p [options]  [input file]

Options:

-o     [output log]		Output log (default: stdout)

Computation options:

-s     [step]			MCMC interval space step (default: 1000)

-c     [#chains]		Max. number of Markov chains (default: 100000)

Starting values:

-e0     [error]			Error rate         (default: random)

-r0     [cont]			Contamination rate (default: random)

-tA0    [tauA]			Tau Ancient Genome        (default: random)

-tC0    [tauC]			Tau Anchor    (default: random)

Range for parameter values:

-e     el,eh			Error rate range          (default: 1e-05,0.1 )

-r     rl,rh			Contamination rate range  (default: 1e-05,0.5 )

-tA    tauAl,tauAh		Tau Ancient Genome range         (default: 1e-06,1   )

-tC    tauCl,tauCh		Tau Anchor range     (default: 1e-06,1   )

# 2-Pop method: BAM file option

[TO ADD]

# 2-Pop method: alternative error rate models

[TO ADD]

# 3-Pop method: input data format

The input data for the 3-pop method is a file that should have at least five columns. An example row would be:

7	0	0.33	0.23	10

This means there are 10 sites in the genome where the ancient genome has 7 ancestral reads and 0 derived reads, where the contaminant/first anchor (admixing) panel has a derived allele frequency of 0.33, and where the second anchor (non-admixing) panel has a derived allele frequency of 0.23. If the file only has 5 columns, the program will assume the contaminant panel is the same as the first anchor panel. The file should also have a header. We provide an example test file in the testData folder, under the name "test_threepop_5col.txt".

If the file has 6 columns, the program will assume that the contaminant and the two anchor panels are distinct. In this case, the third column refers to the first anchor (admixing) panel, the fourth column refers to the second anchor (non-admixing) panel, and the fifth column refers to the contaminant panel. For example:

8	1	0.40	0.55	0.83	4

This row means there are 4 sites where the ancient genome has 8 ancestral reads and 1 derived read, and where the admixing anchor panel has a derived allele frequency of 0.40, the non-admixing anchor panel has a derived allele frequency of 0.55, and where the contaminant panel has a derived allele frequency of 0.83. We provide an example test file in the testData folder, under the name "test_twopop_6col.txt".

NOTE: when comparing multiple candidate contaminant populations, it is best to use only sites that are segregating in all populations that are tested as contaminants, so that the likelihood is composed of the same number of sites in each case.

# 3-Pop method: usage

To run type:

./dice -3p [options]  [input file]

Options:

-o     [output log]		Output log (default: stdout)

Computation options:

-s     [step]			MCMC interval space step (default: 1000)

-c     [#chains]		Max. number of Markov chains (default: 100000)

Starting values:

-e0     [error]			Error rate         (default: random)

-r0     [cont]			Contamination rate (default: random)

-tA0    [tauA]			Tau Ancient Genome        (default: random)

-tC0    [tauC]			Tau Anchor    (default: random)

-aR0    [admR]			Admixture time     (default: random)

-aT0    [admT]			Admixture rate     (default: random)

Range for parameter values:

-e     el,eh			Error rate range          (default: 1e-05,0.1 )

-r     rl,rh			Contamination rate range  (default: 1e-05,0.5 )

-tA    tauAl,tauAh		Tau Ancient Genome range         (default: 1e-06,1   )

-tC    tauCl,tauCh		Tau Anchor range     (default: 1e-06,1   )

-aR    admRl,admRh		Admixture time range      (default: 1e-06,0.5 )

-aT    admTl,admTh		Admixture rate range      (default: 0.05,0.11 )

Population specific constants:

-idy     [drift]		Drift specific to admixing panel (default: 0.16)

-idz     [drift]		Drift specific to non-admixing panel (default: 0.16)

-nc      [num c]		Number of sampled chromosomes from the admixing panel (default: 20)

-nb      [num b]		Number of sampled chromosomes from the non-admixing panel (default: 20)

# 3-Pop method: calculating drifts specific to each anchor population

We provide an R script to calculate the drift times (inner drift Y and inner drift Z) specific to each anchor population, which should be inputted into the command line in the options -idy and -idz when running the 3-pop method. The input for this R script is a tab-separated file that should contain 5 rows, each describing a particular configuration of derived allele frequencies in the two populations. An example line would be:

11	100	78	100	21

This would mean that there are 21 sites where the panel from the first population has 11 derived alleles out of 100 sampled alleles, and the panel from the second population has 78 derived alleles out of 100 sampled alleles. We provide an example input file for the R script in the testData folder, under the name "test_calcdrifts_input.txt". To calculate drifts on this file, one would need to run the following script:

Rscript CalcDrifts.R test_calcdrifts_input.txt > test_calcdrifts_output.txt

# 3-Pop method: BAM file option

[TO ADD]

# 3-Pop method: alternative error rate models

[TO ADD]
