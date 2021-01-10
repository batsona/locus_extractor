# The purpose of this script is to extract loci of interest from 1+ genome assemblies in a user-specified directory. This script will use the following steps to achieve this end:

1. Make blast databases for user-specified genome assemblies (genome assemblies must end in '.fa', '.fna', or '.fasta') if not already created.
   a. This step will only occur if a blastdatabse is not present in the working directory which the genomes/fastas of interest are present.

2. Use blastn to identify the coordinates of each fasta in the query file (>1 fasta sequences may be present in the query file).

3. **Only** the top hit for each query sequence is extracted from each genome assembly/fasta of interest (if that locus is identified by blastn).

4. All of the extracted loci which are derived from the same query sequence are combined into one fasta file. For example, if the input query.fasta file has five sequences, named sequentially (i.e., '>sequence_1', '>sequence_2',..., '>sequence_5'), the output files will be named 'sequence_1', 'sequence_2', etc. Each output file will contain the sequence of interest that was found in each input genome file. If the query file was search against 6 genome assemblies, and each assembly had one copy of sequence_1 to sequence_5, then each output file would contain 6 sequences. 

Required binaries:
* blastn (2.9.0)
* makeblastdb (2.9.0)

Non-standard Python libraries:
* biopython (Bio) (1.74)
* pandas (0.24.2)
* pyfaidx (0.5.7)

## Example usage
1. Use blastn to align DNA sequences in `my_path/to/query_sequences/my_queries.fasta` against all the genome assemblies in `my_path/to/assemblies/`:

`python extract_loci.py -a my_path/to/assemblies/ -q my_path/to/query_sequences/my_queries.fasta -o my_path/to/output_directory/`

2. Same as the above example, except include 500 bp of flanking sequence in the five prime direction:

`python extract_loci.py -a my_path/to/assemblies/ -q my_path/to/query_sequences/my_queries.fasta -o my_path/to/output_directory/ -l 500`
