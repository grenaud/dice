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

# 2-Pop model: input data format

The input data for the 2-pop method is a file that should have at least four columns. Each row in the file denotes a particular configuration of ancestral reads, derived read and anchor/contaminant allele frequencies. An example row would be:

3	5	0.48	8

This means there are 8 sites in the genome where the ancient genome has 3 ancestral reads and 5 derived reads, and where the anchor/contaminant panel has a derived allele frequency of 0.48. If the file only has 4 columns, the program will assume the contaminant panel is the same as the anchor panel. The file should also have a header. We provide an example test file in the testData folder, under the name "test_twopop_4col.txt".

If the file has 5 columns, the program will assume that the contaminant and anchor panels are distinct. In this case, the third column refers to the anchor panel, while the fourth column refers to the contaminant panel. For example:

3	5	0.48	0.96	8

This row means there are 8 sites where the ancient genome has 3 ancestral reads and 5 derived reads, and where the anchor panel has a derived allele frequency of 0.48 and where the contaminant panel has a derived allele frequency of 0.96. We provide an example test file in the testData folder, under the name "test_twopop_5col.txt".


# 2-Pop model: usage

./dice -2p [options]  [input file]

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



# 2-Pop model: BAM file option

[TO ADD]

# 3-Pop model: input data format

The input data for the 3-pop method is a file that should have at least five columns. An example row would be:

7	0	0.33	0.23	10


This means there are 10 sites in the genome where the ancient genome has 7 ancestral reads and 0 derived reads, where the contaminant/first anchor (admixing) panel has a derived allele frequency of 0.33, and where the second anchor (non-admixing) panel has a derived allele frequency of 0.23. If the file only has 5 columns, the program will assume the contaminant panel is the same as the first anchor panel. The file should also have a header. We provide an example test file in the testData folder, under the name "test_threepop_5col.txt".

If the file has 6 columns, the program will assume that the contaminant and the two anchor panels are distinct. In this case, the third column refers to the first anchor (admixing) panel, the fourth column refers to the second anchor (non-admixing) panel, and the fifth column refers to the contaminant panel. For example:

8	1	0.40	0.55	0.83	4

This row means there are 4 sites where the ancient genome has 8 ancestral reads and 1 derived read, and where the admixing anchor panel has a derived allele frequency of 0.40, the non-admixing anchor panel has a derived allele frequency of 0.55, and where the contaminant panel has a derived allele frequency of 0.83. We provide an example test file in the testData folder, under the name "test_twopop_6col.txt".


# 3-Pop model: usage

./dice -3p [options]  [input file]

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

-idy     [drift]		Inner drift Y (default: 0.16)

-idz     [drift]		Inner drift Z (default: 0.16)

-nc      [num c]		Number nC (default: 20)

-nb      [num b]		Number nB (default: 20)


# 3-Pop model: calculating drifts specific to each anchor population

[TO ADD]

# 3-Pop model: BAM file option

[TO ADD]
