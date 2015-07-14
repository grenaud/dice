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

The input data for the 2-pop method is a file that should have at least four columns. Each row in the file denotes a particular configuration of ancestral reads, derived read and anchor/contaminant allele frequencies. All columns should be separated by tabs. Here is an example row:

3	5	0.48	8

This means there are 8 sites in the genome where the ancient genome has 3 ancestral reads and 5 derived reads, and where the anchor/contaminant panel has a derived allele frequency of 0.48. If the file only has 4 columns, the program will assume the contaminant panel is the same as the anchor panel. The file should also have a header. We provide an example test file in the testData folder, under the name "test_twopop_4col.txt".

If the file has 5 columns, the program will assume that the contaminant and anchor panels are distinct. In this case, the third column refers to the anchor panel, while the fourth column refers to the contaminant panel. For example:

3	5	0.48	0.96	8

This row means there are 8 sites where the ancient genome has 3 ancestral reads and 5 derived reads, and where the anchor panel has a derived allele frequency of 0.48 and where the contaminant panel has a derived allele frequency of 0.96. We provide an example test file in the testData folder, under the name "test_twopop_5col.txt".


# 2-Pop model - command line:

The command line to run DICE on a file is as follows:

./src/dice -2p [input.txt] > [output.txt]

This will print the steps of the MCMC chain onto the output.txt file.

# 2-Pop model - BAM file option:
[TO ADD]

# 3-Pop model
[TO ADD]

# 3-Pop model - BAM file options
[TO ADD]

# 3-Pop model - Calculating drifts specific to each anchor population.
[TO ADD]
