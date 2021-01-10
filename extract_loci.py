import argparse, glob, os, time
from Bio import SeqIO
import pandas as pd
from pyfaidx import Fasta

#The main function of this script is extracting loci from a genome assembly
#First, a database will be made for the assembly (-a or --assembly)
#second, an input fasta file with n loci will be BLAST searched against the assembly (-q or --query)
#third, the loci will be extracted with a user specified buffer (-b or --buffer)
#finally, the loci will be split up into separate fasta files

#
help_message = 'The purpose of this program is to extract user specified loci from 1+ genome assemblies.'

#help arg
argument_parser = argparse.ArgumentParser(description = help_message)
#assemblies arg
argument_parser.add_argument("-a", "--assemblies", help = "Give the full path to the assemblies, do not give the assembly names. This argument expects that all of the assemblies of interest are in the same working directory.")
#query arg
argument_parser.add_argument("-q", "--query", help = "Give the full path and filename for the query. More than one locus may be used.")
#output directory
argument_parser.add_argument("-o", "--outdir", help = "Give the full path for the output directory of this script.")
#buffer argument
argument_parser.add_argument("-l", "--left_buffer", help = "Give the 5' distance that flanks the start and end positions. i.e., if a value of 30 is set, 30 bp in the 5' direction will be extracted with the locus of interest. Default = 0.")
#buffer argument
argument_parser.add_argument("-r", "--right_buffer", help = "Give the 3' distance that flanks the end positions. i.e., if a value of 30 is set, 30 bp in the 3' direction will be extracted with the locus of interest. Default = 0.")

args         = argument_parser.parse_args()

#set variables for the arguments passed to the assemblies, and the query file
assembly_dir = args.assemblies
query_dir    = args.query
output_dir   = args.outdir + "/extracted_loci_" + time.strftime("%Y%m%d_%H%M%S")

#this is just to check that the base pair buffer is equivalent to the user-specified buffer
#if not, the default buffer is '0'
if args.left_buffer:
	left_buffer = int(args.left_buffer)
else:
    left_buffer = 0

if args.right_buffer:
    right_buffer = int(args.right_buffer)
else:
    right_buffer = 0

####################################################################################################################
####################################################################################################################
##################################################                ##################################################
##################################################   FUNCTIONS    ##################################################
##################################################                ##################################################
####################################################################################################################
####################################################################################################################
def get_assemblies(full_path_to_assemblies):
	"""The purpose of this function is to gather the genomes assemblies in the user-specified path. This function 
	will only look for files that end with .fa, .fasta, or .fna extensions."""
	assemblies = []
	for extension in ['*.fasta', '*.fa', '*.fna']:
	    for file in glob.glob(full_path_to_assemblies + '/' + extension):
	        assemblies.append(file)
	return assemblies

def get_basename(path_to_assembly):
	"""The purpose of this function is to get the basename of any file with an expected extensions. It is a verbose command, but
	the function first splits off the extensions from the file path (os.path.splitext(path_to_assembly)[0]) and keeps the path only.
	Next, the file name is kept but the rest of the path is discarded with os.path.split(result)[1]. The returned string should 
	represent the name of the assembly."""
	return os.path.split(os.path.splitext(path_to_assembly)[0])[1]

def count_assemblies(path_to_assemblies):
	"""The purpose of this function is to count the number of assemblies that were identified in the function get_assemblies()"""
	return len(get_assemblies(path_to_assemblies))

def make_blast_database(full_path_to_assembly):
	"""The purpose of this function is to make a blast database for a specific genome assembly."""
	assembly_database = full_path_to_assembly + ".nin"
	if os.path.exists(assembly_database):
		print(assembly_database, "exists for", get_basename(assembly_database))
	elif os.path.exists(full_path_to_assembly):
		cmd = "makeblastdb -in " + full_path_to_assembly + " -dbtype nucl"
		os.system(cmd)
	else:
		print(f"{full_path_to_assembly} does not exist, exiting script.")
		exit()

def blastn(query_file, blast_database, output_path_name):
	"""This function calls blast with the user-specified query for each genome assembly."""
	cmd = "blastn -query " + query_file + " -db " + blast_database + " -outfmt 6 -out " + output_path_name
	os.system(cmd)

def fasta_iter(query_file_path):
	"""The purpose of this function is to collect all of the FASTA headers from the query file. These headers
	will be used for parsing the blastn output."""
	with open(query_file_path, "r", newline = '') as file:
		fasta_ids = []
		for fasta in SeqIO.parse(file, "fasta"):
			fasta_ids.append(fasta.id)
	return fasta_ids

def parse_outformat_six(blast_outformat_6_file, query_id):
	labels = ['qseqid', 'sseqid', 'pident', 'align_len', 'num_mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'ssend', 'evalue', 'bitscore']
	df = pd.read_csv(blast_outformat_6_file, sep = "\t", names = labels)
	df2 = df.filter(items=['qseqid', 'sseqid', 'sstart', 'ssend'])
	df3 = df2.loc[df2['qseqid'] == query_id]
	if len(df3) > 0:
		contig = df3.iloc[0]['sseqid']
		start_position = df3.iloc[0]['sstart']
		end_position = df3.iloc[0]['ssend']
		return contig, start_position, end_position #these poisitions are not necessarily phased correctly. This will be corrected in a different function.
	else:
		return 0

def extract_locus_with_pyfaidx(blast_outformat_6_file, query_id, full_path_to_assembly, full_path_to_output, left_buffer = left_buffer, right_buffer = right_buffer):
	"""The purpose of this function is to extract the regions from each contig of interest."""
	 #collects regions from blastout file
	result                 = parse_outformat_six(blast_outformat_6_file, query_id)
	if result != 0:
		output_header      = f"{get_basename(full_path_to_assembly)}_{query_id}"
		output_file_name   = f"{full_path_to_output}/{output_header}.tmp.fasta"
		output             = buffer_corrector(full_path_to_assembly, result[0], result[1], result[2], left_buffer, right_buffer) #result[0] = stored pyfaidx contig; result[1] = start_position; result[2] = end_position
		extracted_sequence = Fasta(full_path_to_assembly)[result[0]][output[0]:output[1]]
		if adjust_coordinates(result[1], result[2])[0] == "Reverse complement":	
			extracted_sequence = Fasta(full_path_to_assembly)[result[0]][output[0]:output[1]].reverse.complement
		else:
			pass		
		with open(output_file_name, "w") as file:
			file.write(f">{output_header}\n")
			file.write(f"{extracted_sequence.seq}\n")
			file.close()

def adjust_coordinates(start_position, end_position):
	"""The purpose of this function is to determine whether the start position of a locus is greater than the end position
	of the locus. If so, this will return False. This function will be helpful for determining whether coordinates need to be 
	swapped for extracting a locus and for reverse complementation of the resulting string."""
	if start_position < end_position:
		return "Don't reverse complement", int(start_position) - 1, int(end_position) #the minus one is used to correct the start position from 1-based (BLAST) to 0-based (Python) coordinates
	else:
		return "Reverse complement", int(end_position) - 1, int(start_position) #reshuffles start and end positions.

def buffer_corrector(full_path_to_assembly, contig, start_position, end_position, left_bufffer = left_buffer, right_buffer = right_buffer):
	"""The purpose of this function is to adjust the left and right (5' and 3') buffers. Sometimes, the buffer that is specified by the user will go out of bounds of the contig.
	In other words, the buffer exceeds the start point of the contig or the end point of the contig."""
	coordinates = adjust_coordinates(start_position, end_position)
	one_contig  = Fasta(full_path_to_assembly)[contig]
	new_start = coordinates[1] - left_buffer
	new_end   = coordinates[2] + right_buffer
	if new_start >= 0 and new_end <= len(one_contig):
		return new_start, new_end
	elif new_start >= 0 and new_end > len(one_contig):
		return new_start, len(one_contig)
	elif new_start < 0 and new_end <= len(one_contig):
		return 0, new_end
	else:
		return coordinates[1], coordinates[2]

def concatenate_tmp_fastas(full_path_to_outdir, query_id):
	"""This function concatenates the all the files with the same query_id and .tmp.fasta extension"""
	cmd = f"cat {full_path_to_outdir}/*{query_id}.tmp.fasta > {full_path_to_outdir}/{query_id}.fasta"
	os.system(cmd)

####################################################################################################################                                                                                                                                              
####################################################################################################################                                                                                                                                                 
##################################################                ##################################################                                                                                                                                                 
##################################################     SCRIPT     ##################################################                                                                                                                                                 
##################################################                ##################################################                                                                                                                                                 
####################################################################################################################                                                                                                                                                 
####################################################################################################################                                  
if __name__ == '__main__':
	print(f"creating output directory...")
	os.system(f"mkdir {output_dir}")

	print("blast searching all assemblies in", assembly_dir, "for", os.path.split(query_dir)[1], "with a buffer of", left_buffer, "and", right_buffer, "in the 5' and 3' direction, respectively.")
	print("resulting fasta files will be output into", output_dir)

	###first step - making blast directories and using BLAST to identify the locations of the query fastas
	print(f"Making BLAST directories and BLAST searching for {os.path.split(query_dir)[1]}")

	for assembly_path in get_assemblies(assembly_dir):
		output_path_name = output_dir + '/' + get_basename(assembly_path) + ".blastout"
		make_blast_database(assembly_path)
		blastn(query_dir, assembly_path, output_path_name)

	###second step - parse .blastout results to collect the top hit for each locus that is present in an assembly
	for fasta in fasta_iter(query_dir):
		print(f"blast searching for {fasta}...")
		for assembly_path in get_assemblies(assembly_dir):
			blast_outformat_file = output_dir + '/' + get_basename(assembly_path) + ".blastout"
			extract_locus_with_pyfaidx(blast_outformat_file, fasta, assembly_path, output_dir, left_buffer, right_buffer)
		concatenate_tmp_fastas(output_dir, fasta)

	####An inelegant solution to deleting all temporary fasta files
	remove_dir = output_dir + "/*.tmp.fasta*"
	print("cleaning up temporary fasta files...")
	for file in glob.glob(remove_dir):
		os.remove(file)









