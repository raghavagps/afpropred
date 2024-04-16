# **AfProPred**
A computational method to predict the Anti-Freezing proteins based on evolutionary profiles
## Introduction
AfProPred is a tool developed by Raghva-Lab in 2024. It is designed to predict whether a protein is Anti-Freezing or not. It utilizes both amino-acid compositions and PSSM as features to make predictions using an ExtraTrees Classifier. AfProPred is also available as web-server at https://webs.iiitd.edu.in/raghava/afpropred. Please read/cite the content about the AfProPred for complete information including algorithm behind the approach.

## PIP Installation
PIP version is also available for easy installation and usage of this tool. The following command is required to install the package 
```
pip install afpropred
```
To know about the available option for the pip package, type the following command:
```
afpropred -h
```
## Standalone
The Standalone version of AfProPred is written in python3 and following libraries are necessary for the successful run:

- scikit-learn==1.3.0
- argparse
- biopython
- numpy
- pandas

## Minimum USAGE
To know about the available option for the standalone, type the following command:
```
python afpropred.py -h
```
To run the example, type the following command:
```
python afpropred.py -i example_input.fa
```
This will predict the probability whether a submitted sequence will localize to the cytoplasm or nucleus. It will use other parameters by default. It will save the output in "outfile.csv" in CSV (comma separated variables).

## Full Usage
```
usage: afpropred.py [-h] -i INPUT [-o OUTPUT] [-t THRESHOLD] [-m {1,2}] [-d {1,2}]
                    [-wd WORKING]

```
```
Please provide following arguments.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence in FASTA format
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.48
  -m {1,2}, --model {1,2}
                        Model: 1: AAC feature based ExtraTrees Classifier , 2: AAC + PSSM
                        feature based ExtraTrees Classifier, by default 1
  -d {1,2}, --display {1,2}
                        Display: 1:AFP, 2: All peptides, by default 2
  -wd WORKING, --working WORKING
                        Working Directory: Temporary directory to write files
```

**Input File:** It allow users to provide input in the FASTA format.

**Output File:** Program will save the results in the CSV format, in case user does not provide output file name, it will be stored in "outfile.csv".

**Threshold:** User should provide threshold between 0 and 1, by default its 0.5.

**Display type:** This option allow users to display only Anti-Freezing proteins or all the input proteins.

**Working Directory:** Directory where intermediate files as well as final results will be saved

AfProPred Package Files
=======================
It contains the following files, brief description of these files given below

INSTALLATION	      : Installations instructions

LICENSE				      : License information

README.md			      : This file provide information about this package

model               : This folder contains two pickled models

afpropred.py        : Main python program

possum              : This folder contains the program POSSUM, that is used to calculate PSSM features

ncbi-blast-2.15.0+  : This folder contains the BLAST executables (not provided). Kindly download the BLAST executables from the following [link](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/) based on your OS. The blast directory should be in the same folder as afpropred.py

example_input.fasta : Example file containing protein sequences in FASTA format

example_output.csv	: Example output file for the program