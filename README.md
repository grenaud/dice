# Readme for DICE

DICE is a Bayesian method to jointly estimate error rates, demographic parameters and the rate of contamination from present-day humans in ancient DNA (aDNA) samples from modern or archaic humans. Our approach is based on an MCMC sampler, and is applicable to nuclear autosomal aDNA data. It can also serve to determine the most probable ancestry of the individual(s) that contaminated the sample. 
 
# Contact

Fernando Racimo: fernandoracimo@gmail.com

Gabriel Renaud: gabriel.reno@gmail.com

# Table of Contents

- Prerequisites
- Downloading
- Compiling and installing
- 2-Pop method: input data format
- 2-Pop method: usage
- Output
- 3-Pop method: input data format
- 3-Pop method: usage
- 3-Pop method: calculating drifts specific to each anchor population
- BAM to DICE format conversion
- BAM file option
- Alternative error models
- How do I get my allele frequencies for my anchor/contaminant/admixing populations?
- Example of running DICE from a raw BAM file

# Prerequisites

C++ libraries:
- cmake
- zlib
- git 

Python libraries:
- dadi
- numpy
- scipy

R libraries:
- bbmle

# Downloading

Make sure you have git installed and type:

    git clone  --depth=1 --recursive https://github.com/grenaud/dice.git


# Compiling and installing

Make sure you are connected to the internet when you build the code. The compiler needs to be able to retrieve tabix from the samtools package.

First, go to the folder where you downloaded the dice folder, and then type the following commands:


    cd bamtools/
    mkdir build
    cd build
    cmake ..
    make
    cd ../..
    cd libgab/
    make
    cd ..
    cd src/
    make
    cd ..


# 2-Pop method: input data format

DICE has two main demographic inference methods: the 2-pop method and the 3-pop method. The former is faster, but the latter allows for the inference of admixture between a present-day human population and the population to which the aDNA sample belongs.

The input data for the 2-pop method in DICE is a file that should have at least four columns. In the section "BAM to DICE format conversion", we provide a way to build this file starting from a BAM file, but the user can also do this with her/his own scripts. Each row in the file denotes a particular configuration of ancestral reads, derived reads and anchor/contaminant allele frequencies, obtained from a panel of present-day human individuals (like the 1000 Genomes Project population panels). An example row would be:

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

To run DICE with the 2-pop method, type:

    src/dice -2p [options]  [input file]

General options:

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

# Output

The output of running  DICE is a set of columns, each containing the steps in the MCMC for a particular parameter. The first column denotes the number of the step in the chain, the second is the log-posterior for that particular step, the third is the error rate, the fourth is the contamination rate, the fifth is the drift parameter for the anchor side of the tree ("tau C"), the sixth is the drift paramer for the ancient individual's side of the tree ("tau A"), the sixth and seventh columns are the admixture rate and time parameters (which are always zero in the two-pop method), and finally the eight column is the acceptance rate of the chain. Additional columns may appear after the last column if the user chose one of the alternative error models (see "Alternative error models" section). We recommend plotting all the parameters to ensure convergence has been reached. As in any other MCMC method, when estimating the posterior distribution, we recommend sampling after a certain burn-in period and every X number of steps, to ensure the chain is well-mixed.

# 3-Pop method: input data format

The input data for the 3-pop method in DICE is a file that should have at least five columns. An example row would be:

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

To run type DICE with the 3-pop method, type:

    src/dice -3p [options]  [input file]


General options:

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


# BAM to DICE format conversion

This section is for users whose starting data is raw aDNA fragments aligned to the nuclear genome in BAM format. We use the word "fragments" because, since aDNA molecules are small, we need the adapters trimmed and the overlapping portions of the reads to be merged (see http://grenaud.github.io/leehom for software to do this). This BAM file has to be sorted with respect to coordinates, and indexed. 

If the user has a BAM file, there are two ways to run DICE with it:
- Convert your BAM file to DICE's native format. This is recommended for first-time users.
- Run DICE directly on the BAM and use deamination profiles and quality scores to infer the error rates at each position. This mode is a bit slower (see "BAM file option" below).

DICE's native format is a simple text file that contains the derived/ancestral base counts and their frequencies in different panel populations (see "input data format" section). To convert your BAM file to DICE's native format, we have created a small program: src/BAM2DICE. This program takes the following arguments:

    src/BAM2DICE [options] [fasta file] [bam file] [region or file with regions to use] [freq for pop1] [freq for pop2] ... 

Description:

[fasta file] : This is the fasta file you supplied the aligner

[bam]        : Sorted and indexed BAM file

[region]     : A list of regions that the program will produce data for. 
	       We recommend using regions with a high
               mapability score. Try to aim to have a least 1M defined sites.
This file has the following format:
               
	       refID1:start-end
	       
	       refID2:start-end
	       
	       refID3:start-end
	       
	       ...

For example:
               
	       chr1:304012-419131
	       
	       chr1:518593-712340


[freq ..]    : A set of 1 or more tab-separated files containing derived allele frequencies from different panel populations. These files should be in mistar format (https://github.com/grenaud/mistar).

For example: 

	       chr	coord	REF,ALT	root	anc	IndividualA

	       7	35190	G,T	0,1:0	0,1:0	122,1:0

A pre-made list of panels for the 1000 Genomes data is available here: https://bioinf.eva.mpg.de/dice/. These can be downloaded automatically by typing  : 

    make allelefreq

Also, refer to section called "How do I get my allele frequencies for my anchor/contaminant/admixing populations?" in this README.

The user must specify which panel should be used as the contaminant panel and the anchor panel. If the user is planning to run the 3-pop method, she/he should also specify which is the admixing anchor panel. The name of the populations much be the same (case sensitive) as in the headers of the frequency files.

--anch                Comma-separated list of anchor populations         (default: all)

--cont                Comma-separated list of contaminant populations    (default: all)

--admx                Comma-separated list of admixing anchor populations      (default: all)

So, for example, if the user types -anch YRI,LWK -cont CEU,CHB, the program will generate 2x2 = 4 files for each possible combination of anchor and contaminant panels, and these files will be valid input for the 2-pop method. Instead, if the user types --anch YRI,LWK  --cont CEU,GBR --admix CEU,GBR as options, the program will produce 2x2x2 = 8 files for each possible combination of anchor, contaminant and admixing panels, and these files will be valid input for the 3-pop method.

DICE can handle gzipped text files so gzip whenever possible to save space.  By default, we discard CpG islands but they can be added back in using the -wcpg. Also, you can flag transitions and transversions using the "-t" option.


# BAM file option

You can also run DICE directly on the BAM file. This mode however is a bit slower than the normal mode. The advantage is that, in this mode, the error rate parameter is not estimated genome-wide, but is computed directly at each site, using mapping quality, base quality and predicted post-mortem damage rates. 

- First, you need to compute a post-mortem DNA damage profile. This is a simple substitution matrix with the following tab-delimited format:

	       A>C  A>G  A>T  C>A  C>G  C>T        G>A  G>C  G>T  T>A  T>C  T>G

	       0.0  0.0  0.0  0.0  0.0  0.0792496  0.0  0.0  0.0  0.0  0.0  0.0

	       0.0  0.0  0.0  0.0  0.0  0.0204847  0.0  0.0  0.0  0.0  0.0  0.0

	       0.0  0.0  0.0  0.0  0.0  0.0183053  0.0  0.0  0.0  0.0  0.0  0.0

	       0.0  0.0  0.0  0.0  0.0  0.0163882  0.0  0.0  0.0  0.0  0.0  0.0

	       0.0  0.0  0.0  0.0  0.0  0.0163684  0.0  0.0  0.0  0.0  0.0  0.0

	       0.0  0.0  0.0  0.0  0.0  0.0163688  0.0  0.0  0.0  0.0  0.0  0.0


The columns here denote the probability of post-mortem substitution for all 12 types of nucleotide changes. Starting from the top, each line represents the position with respect to the 5' end of the read. So, for example, the first row corresponds to the position next to the 5' end, the second row corresponds to 2 positions away from the 5' end, the third row to 3 positions away from the 5' end, etc.

Ideally, you should have two separate deamination profiles, one for the 5' end and one for the 3' end. You can use the simple "bam2prof" tool to generate those: https://github.com/grenaud/schmutzi/blob/master/bam2prof.cpp

- Second, run diceBAM

# Alternative error models

By default, DICE uses a single error parameter for the entire dataset. However, in ancient DNA datasets, transitions tend to have a greater error rate due to deamination, and there is also a chance the alleles were mis-polarized. Hence we have the following alternative error rate models:

- a single error parameter model (default)
- a model with two separate error parameters, one for transitions and one for transversions. This mode will be triggered automatically if the input data was  previously flagged by BAM2DICE using the "-t" option. 
- a probabilistic two-error rate model, with two different error parameters and a third parameter (pe) that determines what proportion of the genome is affected by the first error parameter, as opposed to the second (see Fu et al. 2014). This mode can be triggered using the "-2e" option when running DICE.
- a probabilistic ancestral state misidentification model, where we have an error rate parameter and a mispolarization parameter. This mode can be triggered using the "-pol" option when running DICE.
- a site-specific error rate model, based on base and mapping qualities, and a post-mortem DNA damage matrix. This mode is implemented in diceBAM (see "BAM file option").


# How do I get my allele frequencies for my anchor/contaminant/admixing populations?

You can either do it on your own using a large panel of individuals from the same population. You need to polarize which allele is the ancestral one.

We provide pre-parsed 1000 Genomes Phase III data in mistar format (https://github.com/grenaud/mistar) which counts the number of alleles as well as ancestral/chimp allele information from the Primate EPO alignments (http://www.ensembl.org/info/genome/compara/analyses.html). To obtain this data, make sure you are connected to the internet and type:

    cd src/   
    make allelefreq
    cd ..

This will use "wget" to retrieve pre-parsed allele information to be used by DICE2BAM (or diceBAM). The codes for the population are found here:  http://www.1000genomes.org/category/frequently-asked-questions/population



# Example of running DICE from a raw BAM file

This section provides an example of running DICE from a raw BAM file. First, make sure you are connected to the internet and download simulated test data using:

    cd src/
    make testdatasimulated
    cd ..

Then, transform the BAM file (all.bam) into native DICE format:
    
    src/BAM2DICE -o testData/simulated/input.dice  -2p testData/simulated/all.ref.fa testData/simulated/all.bam  testData/simulated/all.med.regions testData/simulated/all.mst.gz
    
the following file should be created:

    testData/simulated/input.dice_Cont_Anch_contHuman.dice

This data had about 5% contamination and tauA and tauC were both equal to 0.5. We will try to predict both parameters:

     src/dice -2p -o testData/simulated/input.dice_Cont_Anch_contHuman.dice.out testData/simulated/input.dice_Cont_Anch_contHuman.dice

This should produce a file called:
     
     testData/simulated/input.dice_Cont_Anch_contHuman.dice.out 

Then plot the posteriors using the following command: 

    src/log2plots.R  testData/simulated/input.dice_Cont_Anch_contHuman.dice.out testData/simulated/input.dice_Cont_Anch_contHuman.dice.out.pdf

This will create the following files:

    testData/simulated/input.dice_Cont_Anch_contHuman.dice.out.e.pdf	 for the error
    testData/simulated/input.dice_Cont_Anch_contHuman.dice.out.it.pdf	 for the log of the post. prob 
    testData/simulated/input.dice_Cont_Anch_contHuman.dice.out.r.pdf	 contamination rate
    testData/simulated/input.dice_Cont_Anch_contHuman.dice.out.tauA.pdf  for tauA
    testData/simulated/input.dice_Cont_Anch_contHuman.dice.out.tauC.pdf  for tauC


# Which genomic regions to take?

You should take some regions that are highly likely to have been mapped at the correct location to avoid wasting time with mismapped fragments. You can find examples of highly mappable regions in:
    https://bioinf.eva.mpg.de/dice/mapabilityTracks/


    
