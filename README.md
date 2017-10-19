## SAM2CONSENSUS

SAM2CONSENSUS is a C program written to get consensus nucleotide sequence from a sam file.

### Compilation
gcc sam2consensus.c -o SAM2CONSENSUS -lm 


## Usage
*SAM2CONSENSUS -i file.sam -o file.fa*

### For coverage, variation and shannon entropy, run the program with -e option (extended output)
*SAM2CONSENSUS -e -i file.sam -o file.fa*

## Please check your file before running
* It generates consensus from a sam file where a single genome used as a reference
* It work on a SAM file NOT on a BAM file (please convert your bam file to sam)
* It is tested on viral and small bacterial genome, large genomes should work, but not tested
* Use the program at your own risk

## Contact
Please send your comments, suggestions and bugs to Dr. Sreenu Vattipally (sreenu.vattipally@glasgow.ac.uk)
