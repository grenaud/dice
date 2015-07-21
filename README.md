# Readme for DICE

DICE is a Bayesian method to jointly infer contamination from present-day humans in ancient DNA samples, error and drift parameters using MCMC. Our approach is applicable to nuclear data. DICE works by computing the likelihood of finding a certain derived allele as contaminant by using the derived allele frequency in a potentially contaminating population.
 
# Contact

Fernando Racimo fernandoracimo@gmail.com

Gabriel Renaud    gabriel.reno@gmail.com


# Prerequisites

C++ libraries:
- cmake
- zlib
- git 
- 
Python libraries:
- dadi
- numpy
- scipy

R libraries:
- bbmle

# Compiling and installing

NOTE: Make sure you are connected to the internet when you build the code. The compiler needs to be able to retrieve tabix from the samtools package.

cd bamtools/

cmake ..

make

cd ../..

cd libgab/

make

cd ..

cd src/

make

cd ..

# Running DICE

The starting data is raw aDNA fragments aligned to the nuclear genome in BAM format. We use the word "fragments" because, since aDNA molecules are small, we need the adapters trimmed and the overlapping portions of the reads to be merged (see http://grenaud.github.io/leehom for software to do this). This BAM file has to be sorted (wrt coordinates) and indexed. 

There are two ways to run DICE:
- Convert to native format (recommended)
- Run DICE directly on the BAM and use deamination profiles and quality scores to infer the error rate. This mode is a bit slower (see section below).

The native format is a simple text file that contains the derived/ancestral base counts and their frequencies in different panel populations (see examples below).

An intersection of the base count at each position in the BAM file and the derived allele frequency must be made. You can do this whichever way you want but we have created a small program to do this: src/BAM2DICE. This program takes the following arguments:

src/BAM2DICE [options] [fasta file] [bam file] [region or file with regions to use] [freq for pop1] [freq for pop2] ... 

Description:

[fasta file] : This is the fasta file you supplied the aligner

[bam]        : Sorted and indexed BAM file

[region]     : A list of regions that the program will produce data for. 
	       We recommend using regionswith a high
               mapability score. Try to aim to have a least 1M defined sites.
            

This file has the following format:
               
	       ------
	       
	       refID1:start-end
	       
	       refID2:start-end
	       
	       refID3:start-end
	       
	       ...

               ------
               

For example:
               
               -----
               
	       chr1:304012-419131
	       
	       chr1:518593-712340


[freq ..]    : A set of files containing allele frequencies from panel 
               population. Which will be used as contaminant, anchor or admixed 
               need to be specified as options. These frequencies use the same 
               used for a software package designed to import, store and 
               process allele frequencies (grenaud.github.io/mistartools).
               ex: 

#chr	coord	REF,ALT	root	anc	IndividualA

7	35190	G,T	0,1:0	0,1:0	122,1:0

DICE can handle gzipped text files so gzip whenever possible to save space.  By default, we discard CpG islands but they can be added back in using the -wcpg. Also, you can flag transitions and transversions using the "-t" option.

Example:
TODO

This will produce the files.. TODO. We combine sites with the same allele frequency base count in the BAM file to increase speed.

	

# Test data

TODO




# 2-Pop method: input data format

The input data for the 2-pop method is a file that should have at least four columns. Each row in the file denotes a particular configuration of ancestral reads, derived read and anchor/contaminant allele frequencies. An example row would be:

3	5	0.48	8

This means there are 8 sites in the genome where the ancient genome has 3 ancestral reads and 5 derived reads, and where the anchor/contaminant panel has a derived allele frequency of 0.48. If the file only has 4 columns, the program will assume the contaminant panel is the same as the anchor panel. The file should also have a header. We provide an example test file in the testData folder, under the name "test_twopop_4col.txt".

If the file has 5 columns, the program will assume that the contaminant and anchor panels are distinct. In this case, the third column refers to the anchor panel, while the fourth column refers to the contaminant panel. For example:

3	5	0.48	0.96	8

This row means there are 8 sites where the ancient genome has 3 ancestral reads and 5 derived reads, and where the anchor panel has a derived allele frequency of 0.48 and where the contaminant panel has a derived allele frequency of 0.96. We provide an example test file in the testData folder, under the name "test_twopop_5col.txt".

NOTES:
- We recommend using large panels for determining the anchor and contaminant frequencies, and rounding the frequencies to 2 decimals, to prevent the list of configurations from becoming too large.
- The allele frequencies for the anchor population should be larger than 0 and smaller than 1. This is not required for the contaminant population, if different from the anchor.
- When comparing multiple candidate contaminant populations, it is best to use only sites that are segregating in all populations that are tested as contaminants, so that the likelihood is composed of the same number of sites in each case.

# 2-Pop method: usage

To run type:

./src/dice -2p [options]  [input file]

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

[TODO]

# 3-Pop method: input data format

The input data for the 3-pop method is a file that should have at least five columns. An example row would be:

7	0	0.33	0.23	10

This means there are 10 sites in the genome where the ancient genome has 7 ancestral reads and 0 derived reads, where the contaminant/first anchor (admixing) panel has a derived allele frequency of 0.33, and where the second anchor (non-admixing) panel has a derived allele frequency of 0.23. If the file only has 5 columns, the program will assume the contaminant panel is the same as the first anchor panel. The file should also have a header. We provide an example test file in the testData folder, under the name "test_threepop_5col.txt".

If the file has 6 columns, the program will assume that the contaminant and the two anchor panels are distinct. In this case, the third column refers to the first anchor (admixing) panel, the fourth column refers to the second anchor (non-admixing) panel, and the fifth column refers to the contaminant panel. For example:

8	1	0.40	0.55	0.83	4

This row means there are 4 sites where the ancient genome has 8 ancestral reads and 1 derived read, and where the admixing anchor panel has a derived allele frequency of 0.40, the non-admixing anchor panel has a derived allele frequency of 0.55, and where the contaminant panel has a derived allele frequency of 0.83. We provide an example test file in the testData folder, under the name "test_twopop_6col.txt".

NOTES:
- We recommend using large panels for determining the anchor and contaminant frequencies, and rounding the frequencies to 2 decimals, to prevent the list of configurations from becoming too large.
- The allele frequencies for the anchor population should be larger than 0 and smaller than 1. This is not required for the contaminant population, if different from the anchor.
- When comparing multiple candidate contaminant populations, it is best to use only sites that are segregating in all populations that are tested as contaminants, so that the likelihood is composed of the same number of sites in each case.

# 3-Pop method: usage

To run type:

./src/dice -3p [options]  [input file]

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

-nc      [num c]		Number of sampled chromosomes from the admixing panel (default: 100)

-nb      [num b]		Number of sampled chromosomes from the non-admixing panel (default: 100)

# 3-Pop method: calculating drifts specific to each anchor population

We provide an R script to calculate the drift times (inner drift Y and inner drift Z) specific to each anchor population, which should be inputted into the command line in the options -idy and -idz when running the 3-pop method. The input for this R script is a tab-separated file that should contain 5 rows, each describing a particular configuration of derived allele frequencies in the two populations. An example line would be:

11	100	78	100	21

This would mean that there are 21 sites where the panel from the first population has 11 derived alleles out of 100 sampled alleles, and the panel from the second population has 78 derived alleles out of 100 sampled alleles. We provide an example input file for the R script in the testData folder, under the name "test_calcdrifts_input.txt". To calculate drifts on this file, one would need to run the following script:

Rscript CalcDrifts.R test_calcdrifts_input.txt > test_calcdrifts_output.txt

# BAM file option

You can also run DICE directly on the BAM file. This mode however is a bit slower than the normal mode since we cannot combine sites together and read fragment needs to be computed independently. The advantage is that, in this mode, the error rate parameter is not estimated genome-wide, but is computed directly at each site, using mapping quality, base quality and deamination rates. 

- First, you need to compute your deamination rates. The deamination profile is a simple substitution matrix with the following tab-delimited format:

-------------
A>C  A>G  A>T  C>A  C>G  C>T        G>A  G>C  G>T  T>A  T>C  T>G
0.0  0.0  0.0  0.0  0.0  0.0792496  0.0  0.0  0.0  0.0  0.0  0.0
0.0  0.0  0.0  0.0  0.0  0.0204847  0.0  0.0  0.0  0.0  0.0  0.0
0.0  0.0  0.0  0.0  0.0  0.0183053  0.0  0.0  0.0  0.0  0.0  0.0
0.0  0.0  0.0  0.0  0.0  0.0163882  0.0  0.0  0.0  0.0  0.0  0.0
0.0  0.0  0.0  0.0  0.0  0.0163684  0.0  0.0  0.0  0.0  0.0  0.0
0.0  0.0  0.0  0.0  0.0  0.0163688  0.0  0.0  0.0  0.0  0.0  0.0
-------------

Where the first base is the one next to the end. Ideally, you should have a deamination profile for the 5' and 3' end. You can use the simple "bam2prof" tool to generate those: https://github.com/grenaud/schmutzi/blob/master/bam2prof.cpp

- Second, run diceBAM

# Alternative error rate models

By default, DICE uses a single error parameter for the entire dataset. However, in ancient DNA datasets, transitions tend to have a greater error rate due to deamination. Hence we have the following error models:

- a single error parameter (default)
- a separate error parameter for transitions and transversion. This requires the data to have been flagged previously by BAM2DICE using the "-t" option. 
- Two different error parameters and a error probability parameter (pe) that will use the first error parameter with probability pe and the second one with probability 1-pe (see Fu et al. 2014). This mode can be triggered using the "-2e" option.
