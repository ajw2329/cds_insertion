#!/usr/bin/python

import sys
from operator import itemgetter
import argparse
from splice_lib import splice_lib
import copy
import re
import subprocess
import glob
import pickle as pkl
import time
from datetime import datetime


stop_codons = ["TAG", "TAA", "TGA"]
			

def add_sequence_to_transcript_dict(transcript_fasta, full_transcript_dict):

	transcript = ""

	with open(transcript_fasta) as file:

		for line in file:

			if line.startswith(">"):

				transcript = re.sub('>', '', line.split()[0].strip())
				full_transcript_dict[transcript]["sequence"] = ""

			else:

				full_transcript_dict[transcript]["sequence"] += line.strip()


def make_cds_dict(transcript_cds_dict):
	'''
		Returns CDS dictionary from CDS exon 
		input file (indexed by gene - this is OK 
		since CDS input file is generated using an 
		established annotation)

		As each CDS is collected, it is 
		matched via splice junctions to 
		potential matching transcripts 
		via the junction-indexed transcript dict.  

		The match is assessed by calling 
		cds_junction_search (which calls 
		assign cds_to_transcript).

	'''

	cds_dict = {}
		
	for transcript, transcript_entry in transcript_cds_dict.items():

		strand = transcript_entry["strand"]
		chrom = transcript_entry["chrom"]
		exon_list = copy.deepcopy(transcript_entry["exons"])

		gene = transcript_entry["gene"]

		start_end = splice_lib.get_feature_start_end(strand, 
													 transcript_entry["exons"])

		flat_exon_list = [str(i) for j in exon_list[:] for i in j]

		test_cds_key = chrom + "_" + "_".join(flat_exon_list) + "_" + strand



		if test_cds_key not in cds_dict:

			cds_dict[test_cds_key] = {
				"cds_start": start_end[0],
				"cds_end": start_end[1],
				"exons": exon_list,
				"junctions": splice_lib.get_junctions(exon_list),
				"chrom": chrom,
				"strand": strand,
				"gene": gene,
				"transcripts": [transcript]

									}

		else:

			cds_dict[test_cds_key]["transcripts"].append(transcript)

	return cds_dict


def generate_start_codon_bedfile(input_cds_dict, outdir):


	outfile = open(outdir + "/start_codons.bed", 'w')

	for cds in input_cds_dict:

		exons = input_cds_dict[cds]["exons"]

		strand = input_cds_dict[cds]["strand"]

		if strand == "+":

			start = str(exons[0][0] - 1)
			end = str(exons[0][0])

		elif strand == "-":

			start = str(exons[-1][1] - 1)
			end = str(exons[-1][1])

		chrom = input_cds_dict[cds]["chrom"]


		outfile.write("\t".join([chrom, 
								 start, 
								 end, 
								 cds, 
								 "1000", 
								 strand]) + "\n")

	outfile.close()



def generate_exon_bedfile(full_transcript_dict, 
						  outdir):

	outfile = open(outdir + "/exons.bed", 'w')
	outfile.truncate()

	for transcript in full_transcript_dict:

		strand = full_transcript_dict[transcript]["strand"]
		chrom = full_transcript_dict[transcript]["chrom"]

		for exon in full_transcript_dict[transcript]["exons"]:

			start = str(exon[0] - 1)
			end = str(exon[1])

			outfile.write("\t".join([chrom, 
									 start, 
									 end, 
									 transcript, 
									 "1000", 
									 strand]) + "\n")

	outfile.close()



def call_bedtools_intersect(outdir, 
							bedtools_exec = "bedtools"):

	bedtools_command = (bedtools_exec + 
					   " intersect -a " + 
					   outdir + 
					   "/start_codons.bed -b " + 
					   outdir + 
					   "/exons.bed -wa -wb -s > " + 
					   outdir + "/start_codon_overlaps.bed")

	try:
		subprocess.check_output(bedtools_command, 
							    shell = True, 
							    stderr = subprocess.STDOUT)

	except subprocess.CalledProcessError as e:
		print(e.output)
		sys.exit("bedtools intersect failed in " +
				 "'run_bedtools_intersect' of cds_insertion.py. Exiting . . . ")



def seek_cds_matches(outdir, full_transcript_dict, input_cds_dict):

	with open(outdir + "/start_codon_overlaps.bed", 'r') as file:

		for line in file:

			entry = line.strip().split()

			cds = entry[3].strip()
			transcript = entry[9].strip()

			cds_entry = copy.deepcopy(input_cds_dict[cds])

			assign_cds_to_transcript(cds_entry, 
									 transcript, 
									 full_transcript_dict, 
									 input_cds_dict)



def assign_cds_to_transcript(cds, 
							 transcript, 
							 full_transcript_dict, 
							 input_cds_dict):

	'''
		Takes as input an entry from cds exon dict, 
		a transcript id, and the full transcript 
		dict.  Checks to see if the CDS is a 
		full match for the transcript, and 
		(if it is) adds it to the full transcript dict.
	'''

	#full_verbose = False

	#if transcript == "MSTRG.11569.2":
#
	#	full_verbose = True

	transcript_entry = full_transcript_dict[transcript]

	annotated = False

	gene = cds["gene"]

	cds_exons = cds["exons"]
	transcript_exons = copy.deepcopy(transcript_entry["exons"])

	transcript_jl = copy.deepcopy(transcript_entry["junctions"])
	cds_jl = cds["junctions"]


	strand = transcript_entry["strand"]
	chrom = transcript_entry["chrom"]


	transcript_boundaries = splice_lib.get_feature_start_end(
		strand, 
		transcript_exons)

	transcript_start = int(transcript_boundaries[0])
	transcript_end = int(transcript_boundaries[1])
	
	cds_start = int(cds["cds_start"])
	cds_end = int(cds["cds_end"])
	transcript_seq = full_transcript_dict[transcript]["sequence"]
	transcript_len = len(transcript_seq)

	start_is_contained, start_containing_exon = splice_lib.position_contained(
		transcript_exons, 
		cds_start)

		## Above returns a list where the first element
		## is boolean and the second is either the 
		## containing exon (if True) or None (if False)

	if start_is_contained:

		valid_cds_start = cds_start

		valid_adj_cds_start = splice_lib.genome_to_transcript_coords(
			cds_start, 
			strand, 
			transcript_exons, 
			"GT")

		translation_attempt = splice_lib.translate_ORF(
			transcript_seq, 
			stop_codons, 
			valid_adj_cds_start)

		translation_attempt_is_successful = translation_attempt[0]

			###Note that there are two possible outcomes: 
			#1) a stop codon is found 
			#2) no stop codon is found, 
			#and we keep both of them 
			#(potentially we may find examples of non-stop decay substrates)

		if translation_attempt_is_successful:

			(translation_attempt_is_successful, 
			 valid_adj_cds_end, 
			 cds_seq, 
			 aa_seq) = translation_attempt

			valid_cds_end = splice_lib.genome_to_transcript_coords(
				valid_adj_cds_end, 
				strand, 
				transcript_exons, 
				"TG")
			
			stop_containing_exon = splice_lib.position_contained(
				transcript_exons, 
				valid_cds_end)[1]

			add_cds_to_transcript_dict(
				strand, 
				transcript_exons, 
				valid_cds_start, 
				valid_cds_end, 
				start_containing_exon, 
				stop_containing_exon, 
				chrom, 
				valid_adj_cds_start, 
				transcript_entry, 
				transcript, 
				gene, 
				annotated, 
				input_cds_dict, 
				cds_seq, 
				aa_seq)



		elif translation_attempt[1] is None:

			(translation_attempt_is_successful, 
			 valid_adj_cds_end, 
			 cds_seq, 
			 aa_seq) = translation_attempt

			if strand == "+":

				left = valid_cds_start
				right = transcript_end

			elif strand == "-":

				left = transcript_end
				right = valid_cds_start

			if valid_cds_start not in transcript_entry["nonstop"]:

				transcript_entry["nonstop"][valid_cds_start] = {
					"start": valid_cds_start,
					"exons": splice_lib.get_exon_subset(
						transcript_exons, 
						left, 
						right),
					"gene": gene
					}

				transcript_entry["gene"] = gene



def add_cds_to_transcript_dict(strand, 
							   transcript_exons, 
							   valid_cds_start, 
							   valid_cds_end,
							   start_containing_exon, 
							   stop_containing_exon, 
							   chrom, 
							   valid_adj_cds_start, 
							   transcript_entry, 
							   transcript, 
							   gene, 
							   annotated, 
							   input_cds_dict, 
							   cds_seq, 
							   aa_seq):


	cds_utr_exons = splice_lib.get_cds_utr_exons(
		strand, 
		transcript_exons, 
		valid_cds_start, 
		valid_cds_end, 
		start_containing_exon, 
		stop_containing_exon)

	five_utr_exons = cds_utr_exons[0]
	cds_exons = cds_utr_exons[1]
	three_utr_exons = cds_utr_exons[2]



	if strand == "+":

		junction_contained_three_utr_exons = copy.deepcopy(three_utr_exons[0:-1])

	elif strand == "-":

		junction_contained_three_utr_exons = copy.deepcopy(three_utr_exons[1:])

	junction_contained_three_utr_length = splice_lib.calc_length_exon_list(
		junction_contained_three_utr_exons)



	flat_cds_list = [i for j in 
					 cds_utr_exons[1] 
					 for i in j]

	flat_cds_list = list(map(str, flat_cds_list))

	cds_key = (chrom + 
			  "_" + 
			  "_".join(flat_cds_list) + 
			  "_" + 
			  strand)

	if not annotated:

		if cds_key in input_cds_dict:

			annotated = True



	if cds_key not in transcript_entry["CDS"]:


		valid_adj_cds_end = splice_lib.genome_to_transcript_coords(
			valid_cds_end, 
		    strand, 
		    transcript_exons, 
		    "GT")


		
		if strand == "+":


			stop_codon = splice_lib.get_exon_subset(
					transcript_exons, 
					genomic_start = splice_lib.genome_to_transcript_coords(
						valid_adj_cds_end + 1, 
					    strand, 
					    transcript_exons, 
					    "TG"), 
					genomic_end = splice_lib.genome_to_transcript_coords(
						valid_adj_cds_end + 3, 
						strand, 
						transcript_exons, 
						"TG"))

		elif strand == "-":


			stop_codon = splice_lib.get_exon_subset(
				transcript_exons, 
				genomic_start = splice_lib.genome_to_transcript_coords(
					valid_adj_cds_end + 3, 
					strand, 
					transcript_exons, 
					"TG"), 
				genomic_end = splice_lib.genome_to_transcript_coords(
					valid_adj_cds_end + 1, 
					strand, 
					transcript_exons, 
					"TG"))

		PTC_status = is_PTC(
			stop_codon[-1][-1], 
			transcript_exons, 
			strand)


		transcript_entry["CDS"][cds_key] = {
			"start": valid_cds_start,
			"end": valid_cds_end,
			"adj_start": valid_adj_cds_start,
			"adj_end": valid_adj_cds_end,
			"annotated": str(annotated),
			"exons": copy.deepcopy(cds_exons),
			"stop_codon": copy.deepcopy(stop_codon),
			"gene": gene,
			"PTC_exon": PTC_status[0],
			"PTC_distance": PTC_status[1],
			"CDS_length": splice_lib.calc_length_exon_list(cds_exons),
			"five_utr_exons": copy.deepcopy(five_utr_exons),
			"three_utr_exons": copy.deepcopy(three_utr_exons),
			"five_utr_seq": transcript_entry["sequence"][0:int(valid_adj_cds_start)],
			"three_utr_seq": transcript_entry["sequence"][int(valid_adj_cds_end) + 1:],
			"cds_seq": cds_seq,
			"aa_seq": aa_seq,
			"five_utr_length": splice_lib.calc_length_exon_list(five_utr_exons),
			"three_utr_length": splice_lib.calc_length_exon_list(three_utr_exons),
			"three_utr_junction_count": len(three_utr_exons) - 1,
			"junction_contained_three_utr_length": junction_contained_three_utr_length,
			"downstream_PTC_junction_distances": PTC_status[2]
			}

		transcript_entry["gene"] = gene

		if PTC_status[0] != "no_PTC" and PTC_status[0] != None:
			if len(junction_contained_three_utr_exons) == 0:

				print("WHOOPS!!!!!!!!!!!!!!!!!", three_utr_exons)

			if len(three_utr_exons) == 1:
				print("ACKACKACKACK", three_utr_exons)

	else: ##This is necessary in case translate ORF from another annotated (also from the same gene) resulted in addition of this same CDS (but without awareness of its annotated status)
		transcript_entry["CDS"][cds_key]["annotated"] = annotated

		#print "Marked existing CDS as annotated"




def is_PTC(
	stop_codon_end, 
	transcript_exons, 
	strand):

	'''
		Determines if stop codon is upstream of splice junction (PTC)

	'''

	
	if strand == "+":
		is_PTC_transcript_exons = copy.deepcopy(transcript_exons)
	elif strand == "-":
		is_PTC_transcript_exons = copy.deepcopy(transcript_exons)[::-1]

	PTC = False

	downstream_distances = []

	
	for i in is_PTC_transcript_exons:

		if stop_codon_end in range(int(i[0]),int(i[1])+1):

			if is_PTC_transcript_exons.index(i) == (len(is_PTC_transcript_exons) - 1):

				return ["no_PTC", "NA", "NA"]
			else:

				PTC = True

				if strand == "+":
					distance = int(i[1]) + 1 - stop_codon_end
				elif strand == "-":
					distance = stop_codon_end + 1 - int(i[0])

				
				exon_distance = [i, distance]
				downstream_distances.append(distance)

		elif PTC == True: #PTC-containing exon has been found - now collect remaining CDS -> downstream junction distances

			if is_PTC_transcript_exons.index(i) != (len(is_PTC_transcript_exons) - 1):

				if len(downstream_distances) == 0:
					downstream_distance = distance + int(i[1]) - int(i[0]) + 1
				else:
					downstream_distance = downstream_distances[-1] + int(i[1]) - int(i[0]) + 1
				downstream_distances.append(downstream_distance)

	else:

		if PTC == True:

			if len(downstream_distances) == 0:
				downstream_distances = "NA"

			exon_distance.append(downstream_distances)
			return exon_distance

		else:
			print("PTC Error ... printing transcript exons ...")
			print(stop_codon_end)
			print(is_PTC_transcript_exons)
			sys.exit("Error: stop codon not contained in spliced transcript")



def output_cds_inserted_gtf(full_transcript_dict, outdir):

	gtf_output = open(outdir + "/cds_inserted_no_ptc.gtf", 'w')
	nmd_gtf_output = open(outdir + "/cds_inserted_ptc.gtf", 'w')
	nonstop_gtf_output = open(outdir + "/cds_inserted_nonstop.gtf", 'w')

	for transcript, transcript_entry in full_transcript_dict.items():

		chrom = re.sub("&","_",transcript_entry["chrom"][:])
		source = "cds_insertion"
		tx_start = str(transcript_entry["exons"][0][0])
		tx_end = str(transcript_entry["exons"][-1][1])
		strand = transcript_entry["strand"][:]
		all_exons = copy.deepcopy(transcript_entry["exons"])

		for i in all_exons:
			i.insert(0, "exon")


		counter = 1
		
		for cds, cds_entry in transcript_entry["CDS"].items():

			nmd = False

			if cds_entry["PTC_exon"] != None and cds_entry["PTC_exon"] != "no_PTC":

				nmd = True


			transcript_id = transcript + "_" + str(counter)

			cds_entry["transcript_version"] = transcript_id


			field_nine = ('gene_id "' + 
						 transcript_entry["gene"] + 
						 '"; transcript_id "' + 
						 transcript_id + 
						 '"; gene_name "' + 
						 transcript_entry["gene"] + 
						 '";')

			if nmd:

				nmd_gtf_output.write("\t".join([chrom, 
												source, 
												"transcript", 
												tx_start, 
												tx_end, 
												".", 
												strand, 
												".", 
												field_nine]) + "\n")
			else:

				gtf_output.write("\t".join([chrom, 
											source, 
											"transcript", 
											tx_start, 
											tx_end, 
											".", 
											strand, 
											".", 
											field_nine]) + "\n")

			all_cds_exons = copy.deepcopy(cds_entry["exons"])

			for i in all_cds_exons:
				i.insert(0, "CDS")

			stop_codon = copy.deepcopy(cds_entry["stop_codon"])

			for i in stop_codon:
				i.insert(0, "stop_codon")


			elements = sorted(all_exons + 
							  all_cds_exons + 
							  stop_codon, key = lambda x: (x[1],x[2])) 

			#syntax help from user Padraic Cunningham stackoverflow https://stackoverflow.com/
			#questions/25046306/sort-a-2d-list-first-by-1st-column-and-then-by-2nd-column

			for element in elements:

				element_type = element[0]
				element_start = str(element[1])
				element_end = str(element[2])

				if nmd:
					nmd_gtf_output.write("\t".join([chrom, 
													source, 
													element_type, 
													element_start, 
													element_end, 
													".", 
													strand, 
													".", 
													field_nine]) + "\n")
				else:
					gtf_output.write("\t".join([chrom, 
												source, 
												element_type, 
												element_start, 
												element_end, 
												".", 
												strand, 
												".", 
												field_nine]) + "\n")

			counter += 1



		counter = 1


		for cds, cds_entry in transcript_entry["nonstop"].items():

			transcript_id = (transcript + 
							"_nonstop_" + 
							str(counter))

			field_nine = ('gene_id "' + 
						 transcript_entry["gene"] + 
						 '"; transcript_id "' + 
						 transcript_id + 
						 '"; gene_name "' + 
						 transcript_entry["gene"] + 
						 '";')

			nonstop_gtf_output.write("\t".join([chrom, 
												source, 
												"transcript", 
												tx_start, 
												tx_end, 
												".", 
												strand, 
												".", 
												field_nine]) + "\n")


			all_cds_exons = copy.deepcopy(cds_entry["exons"])

			for i in all_cds_exons:
				i.insert(0, "CDS")

			elements = sorted(all_exons + 
							  all_cds_exons, 
							  key = lambda x: (x[1],x[2])) 

			#syntax help for above from user Padraic Cunningham stackoverflow https://stackoverflow.com/questions/25046306/sort-a-2d-list-first-by-1st-column-and-then-by-2nd-column

			for element in elements:

				element_type = element[0]
				element_start = str(element[1])
				element_end = str(element[2])

				nonstop_gtf_output.write("\t".join([chrom, 
													source, 
													element_type, 
													element_start, 
													element_end, 
													".", 
													strand, 
													".", 
													field_nine]) + "\n")

			counter += 1

	gtf_output.close()
	nmd_gtf_output.close()
	nonstop_gtf_output.close()


def characterize_transcripts_as_nmd_nsd(full_transcript_dict, 
										ptc_distance_nmd_threshold):


	for transcript, transcript_entry in full_transcript_dict.items():

		try:
			normal_cds_count = len(transcript_entry["CDS"])
		except KeyError:
			normal_cds_count = 0
		try:
			nonstop_cds_count = len(transcript_entry["nonstop"])
		except KeyError:
			nonstop_cds_count = 0

		cds_count = normal_cds_count + nonstop_cds_count

		PTC_distances = []
		max_downstream_PTC_distances = []
		CDS_lengths = []
		five_utr_lengths = []
		three_utr_lengths = []
		three_utr_junction_counts = []
		junction_contained_three_utr_lengths = []

		always_nmd = "False"
		sometimes_nmd = "False"
		always_nonstop = "False"
		sometimes_nonstop = "False"



		for cds, cds_entry in transcript_entry["CDS"].items():

			if cds_entry["PTC_exon"] != "no_PTC":

				PTC_distances.append(cds_entry["PTC_distance"])

				if cds_entry["downstream_PTC_junction_distances"] != "NA":

					max_downstream_PTC_distances.append(
						max(cds_entry["downstream_PTC_junction_distances"]))


			CDS_lengths.append(
				cds_entry["CDS_length"])

			five_utr_lengths.append(
				cds_entry["five_utr_length"])

			three_utr_lengths.append(
				cds_entry["three_utr_length"])

			three_utr_junction_counts.append(
				cds_entry["three_utr_junction_count"])

			junction_contained_three_utr_lengths.append(
				cds_entry["junction_contained_three_utr_length"])


		if (all([
				 i >= ptc_distance_nmd_threshold 
				 for i in 
				 max_downstream_PTC_distances 
				 if i != "NA"
				 ]) and 
			len(max_downstream_PTC_distances) > 0 and
			len(max_downstream_PTC_distances) == normal_cds_count):

			always_nmd = "True"
			sometimes_nmd = "True"

		elif (any([
				   i >= ptc_distance_nmd_threshold 
				   for i in 
				   max_downstream_PTC_distances 
				   if i != "NA"]) and 
			  len(max_downstream_PTC_distances) > 0 and
			  normal_cds_count > 0):

			sometimes_nmd = "True"

		if nonstop_cds_count > 0:
			sometimes_nonstop = "True"
			if normal_cds_count == 0:
				always_nonstop = "True"


		transcript_entry["always_nmd"] = always_nmd
		transcript_entry["sometimes_nmd"] = sometimes_nmd
		transcript_entry["always_nonstop"] = always_nonstop
		transcript_entry["sometimes_nonstop"] = sometimes_nonstop


		if len(CDS_lengths) == 0:
			CDS_lengths = ["NA"]
		if len(five_utr_lengths) == 0:
			five_utr_lengths = ["NA"]
		if len(three_utr_lengths) == 0:
			three_utr_lengths = ["NA"]
		if len(three_utr_junction_counts) == 0:
			three_utr_junction_counts = ["NA"]
		if len(junction_contained_three_utr_lengths) == 0:
			junction_contained_three_utr_lengths = ["NA"]
		if len(PTC_distances) == 0:
			PTC_distances = ["NA"]
		if len(max_downstream_PTC_distances) == 0:
			max_downstream_PTC_distances = ["NA"]

		transcript_entry["CDS_lengths"] = ",".join(map(str, CDS_lengths))
		transcript_entry["five_utr_lengths"] = ",".join(map(str, five_utr_lengths))
		transcript_entry["three_utr_lengths"] = ",".join(map(str, three_utr_lengths))
		transcript_entry["three_utr_junction_counts"] = ",".join(map(str, three_utr_junction_counts))
		transcript_entry["junction_contained_three_utr_lengths"] = ",".join(map(str, junction_contained_three_utr_lengths))
		transcript_entry["max_downstream_PTC_distances"] = ",".join(map(str, max_downstream_PTC_distances))
		transcript_entry["PTC_distances"] = ",".join(map(str, PTC_distances))
		transcript_entry["normal_cds_count"] = normal_cds_count
		transcript_entry["nonstop_cds_count"] = nonstop_cds_count



def output_transcript_table(full_transcript_dict, 
						    outdir):

	output_table = open(outdir + "/transcript_characteristics.tsv", 'w')

	output_table.write("\t".join(["gene", 
								  "transcript_id", 
								  "normal_cds_count", 
								  "nonstop_cds_count", 
								  "always_nmd", 
								  "sometimes_nmd", 
								  "always_nonstop",
								  "sometimes_nonstop",
								  "PTC_distances",
								  "max_downstream_PTC_distances",
								  "CDS_lengths",
								  "five_utr_lengths",
								  "three_utr_lengths",
								  "three_utr_junction_counts"]) + "\n")


	for transcript, transcript_entry in full_transcript_dict.items():

		output_table.write("\t".join([transcript_entry["gene"], 
									  transcript, 
									  str(transcript_entry["normal_cds_count"]), 
									  str(transcript_entry["nonstop_cds_count"]), 
									  transcript_entry["always_nmd"], 
									  transcript_entry["sometimes_nmd"], 
									  transcript_entry["always_nonstop"], 
									  transcript_entry["sometimes_nonstop"], 
									  transcript_entry["PTC_distances"], 
									  transcript_entry["max_downstream_PTC_distances"], 
									  transcript_entry["CDS_lengths"], 
									  transcript_entry["five_utr_lengths"], 
									  transcript_entry["three_utr_lengths"], 
									  transcript_entry["three_utr_junction_counts"]]) + "\n")

	output_table.close()


def output_aa_sequence(full_transcript_dict, 
				       outdir):

	output_table = open(outdir + "/transcript_aa_seq.tsv", 'w')

	output_table.write("\t".join(["gene", 
								  "transcript_id",
								  "transcript_version", 
								  "cds_id", 
								  "aa_seq"]) + "\n")

	for transcript, transcript_entry in full_transcript_dict.items():

		for cds, cds_entry in transcript_entry["CDS"].items():

			output_table.write("\t".join([transcript_entry["gene"], 
										  transcript,
										  cds_entry["transcript_version"], 
										  cds, 
										  cds_entry["aa_seq"]]) + "\n")

	output_table.close()




def make_bb_file(prefix, 
				 bigGenePred_as, 
				 chrNameLength, 
				 outdir, 
				 gtfToGenePred_exe = "gtfToGenePred", 
				 genePredToBigGenePred_exe = "genePredToBigGenePred", 
				 bedToBigBed_exe = "bedToBigBed"):


	gtfToGenePred_command = (gtfToGenePred_exe + 
							" " + 
							outdir + 
							"/" + 
							prefix + 
							".gtf " +  
							outdir + 
							"/" + 
							prefix + 
							".genepred")

	genePredToBigGenePred_command = (genePredToBigGenePred_exe + 
									" " + 
									outdir + 
									"/" + 
									prefix + 
									".genepred " +  
									outdir + 
									"/" + 
									prefix + 
									".bgp")

	sort_command = ("sort -k1,1 -k2,2n " + 
				   outdir + 
				   "/" + 
				   prefix + 
				   ".bgp > " + 
				   outdir + 
				   "/" + 
				   prefix + 
				   ".sorted.bgp")

	bedToBigBed_command = (bedToBigBed_exe + 
						  " -type=bed12+8 -tab -as=" + 
						  bigGenePred_as + 
						  " " + 
						  outdir + 
						  "/" + 
						  prefix + 
						  ".sorted.bgp " + 
						  chrNameLength + 
						  " " + 
						  outdir + 
						  "/" + 
						  prefix + 
						  ".bb")

	#print gtfToGenePred_command


	try:
		subprocess.check_output(gtfToGenePred_command, 
								shell = True, 
								stderr = subprocess.STDOUT)

		subprocess.check_output(genePredToBigGenePred_command, 
								shell = True, 
								stderr = subprocess.STDOUT)

		subprocess.check_output(sort_command, 
								shell = True, 
								stderr = subprocess.STDOUT)

		subprocess.check_output(bedToBigBed_command, 
								shell = True, 
								stderr = subprocess.STDOUT)

	except subprocess.CalledProcessError as e:
		print(e.output)
		sys.exit("Kent tool failed in 'make_bb_file' " +
			     "of cds_insertion.py. Exiting . . . ")


def output_cds_utr_exons(full_transcript_dict, outdir):

	output_bed = open(outdir + "/all_cds_utr_exons.bed", 'w')

	for transcript, transcript_entry in full_transcript_dict.items():

		chrom = transcript_entry["chrom"]
		strand = transcript_entry["strand"]

		for cds, cds_entry in transcript_entry["CDS"].items():

			for exon in cds_entry["exons"]:

				entry = [chrom, 
						 str(exon[0] - 1), 
						 str(exon[1]), 
						 "CDS", 
						 "1000", 
						 strand]

				output_bed.write("\t".join(entry) + "\n")

			for exon in cds_entry["five_utr_exons"]:

				entry = [chrom, 
						 str(exon[0] - 1), 
						 str(exon[1]), 
						 "UTR5", 
						 "1000", 
						 strand]

				output_bed.write("\t".join(entry) + "\n")

			for exon in cds_entry["three_utr_exons"]:

				entry = [chrom, 
						 str(exon[0] - 1), 
						 str(exon[1]), 
						 "UTR3", 
						 "1000", 
						 strand]

				output_bed.write("\t".join(entry) + "\n")

	output_bed.close()


def main(args):

	start_time = time.time()

	print(("{0}: {1} seconds elapsed. Starting cds_insertion. " +
		   "Now parsing arguments."
		   ).format(str(datetime.now().replace(microsecond = 0)), 
		   			str(round(time.time() - start_time, 1))))

	parser = argparse.ArgumentParser()
	parser.add_argument("--transcript_gtf", 
						type = str, 
						help = "Input transcript gtf filename", 
						required = True)

	parser.add_argument("--transcript_fasta", 
						type = str, 
						help = "Input transcript fasta (assumes generation with gffread)", 
						required = True)

	parser.add_argument("--annotation_gtf", 
						type = str, 
						help = "Input gtf filename that contains CDS (preferably CCDS) features", 
						required = True)

	parser.add_argument("--CCDS", 
						action = "store_true", 
						help = "If set, CDS used will be restricted to CCDS.")

	parser.add_argument("--outdir", 
						type = str, 
						help = "Path for output files", 
						required = True)

	parser.add_argument("--bigGenePred_as_path", 
						type = str, 
						help = "Path to bigGenePred.as file")

	parser.add_argument("--gtfToGenePred_path", 
						type = str, 
						help = "Path to gtfToGenePred exe", 
						default = "gtfToGenePred")

	parser.add_argument("--genePredToBigGenePred_path", 
						type = str, 
						help = "Path to genePredToBigGenePred exe", 
						default = "genePredToBigGenePred")

	parser.add_argument("--bedToBigBed_path", 
						type = str, 
						help = "Path to bedToBigBed exe", 
						default = "bedToBigBed")

	parser.add_argument("--bedtools_path", 
						type = str, 
						help = "Path to bedtools exe", 
						default = "bedtools")

	parser.add_argument("--chrNameLength_path", 
						type = str, 
						help = "Path to chrNameLength file")

	parser.add_argument("--make_bigBed", 
						action = "store_true", 
						help = "If set, will attempt to use Kent " +
							   "utils to output bigbed file from gtfs. " +
							   "Requires chrNameLength_path and " + 
							   "bigGenePred_as_path to be set.  If " + 
							   "Kent utils are not in PATH, set the " + 
							   "paths to the executables as well.")

	parser.add_argument("--ptc_dist", 
						type = int, 
						help = "Set minimum PTC distance " + 
						"required for transcript to be " + 
						"considered putative NMD substrate. " + 
						"default = 55", 
						default = 55)



	args = parser.parse_args(args)


	transcript_gtf = args.transcript_gtf
	transcript_fasta = args.transcript_fasta
	annotation_gtf = args.annotation_gtf
	CCDS = args.CCDS
	output_directory = args.outdir
	gtfToGenePred_path = args.gtfToGenePred_path
	genePredToBigGenePred_path = args.genePredToBigGenePred_path
	bedToBigBed_path = args.bedToBigBed_path
	bigGenePred_as_path = args.bigGenePred_as_path
	chrNameLength_path = args.chrNameLength_path
	make_bigBed = args.make_bigBed
	ptc_dist = args.ptc_dist
	bedtools_path = args.bedtools_path

	subprocess.call("mkdir -p " + output_directory, shell = True)

	print(("{0}: {1} seconds elapsed. Arguments parsed. " + 
		   "Now importing (unannotated) transcript GTF file."
		   ).format(str(datetime.now().replace(microsecond = 0)), 
		            str(round(time.time() - start_time, 1))))

	full_transcript_dict = splice_lib.generate_standard_transcript_dict(transcript_gtf)
	splice_lib.sort_transcript_dict_exons(full_transcript_dict)
	splice_lib.add_junctions_to_transcript_dict(full_transcript_dict)


	print(("{0}: {1} seconds elapsed. Transcript GTF file imported. " + 
	   "Now importing annotated transcript GTF file."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))

	if CCDS == True:
		transcript_cds_dict = splice_lib.generate_standard_transcript_dict(annotation_gtf, ccds = CCDS, feature = "CDS")
		splice_lib.sort_transcript_dict_exons(transcript_cds_dict)
		splice_lib.add_junctions_to_transcript_dict(transcript_cds_dict)
	else:
		transcript_cds_dict = splice_lib.generate_standard_transcript_dict(annotation_gtf, feature = "CDS")
		splice_lib.sort_transcript_dict_exons(transcript_cds_dict)
		splice_lib.add_junctions_to_transcript_dict(transcript_cds_dict)


	print(("{0}: {1} seconds elapsed. Annotated transcript GTF file imported. " + 
	   "Now importing transcript sequences."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))

	add_sequence_to_transcript_dict(transcript_fasta, full_transcript_dict)


	print(("{0}: {1} seconds elapsed. Transcript sequences imported. " + 
	   "Indexing annotated transcripts by junction."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	junction_transcript_dict = splice_lib.index_transcripts_by_junctions(full_transcript_dict)



	print(("{0}: {1} seconds elapsed. Transcripts indexed by junction. " + 
	   "Now generating CDS dict."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	input_cds_dict = make_cds_dict(transcript_cds_dict)



	print(("{0}: {1} seconds elapsed. CDS dict generated. " + 
	   "Now generating start coding bed file."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	generate_start_codon_bedfile(input_cds_dict, output_directory)


	print(("{0}: {1} seconds elapsed. Start codon bedfile generated. " + 
	   "Now generating exon bed file."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	generate_exon_bedfile(full_transcript_dict, output_directory)


	print(("{0}: {1} seconds elapsed. Exon bed file generated. " + 
	   "Now running bedtools intersect to identify exon-overlapping start codons."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	call_bedtools_intersect(output_directory, bedtools_path)


	print(("{0}: {1} seconds elapsed. Start codon overlaps identified. " + 
	   "Now attempting in silico translation."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	seek_cds_matches(output_directory, full_transcript_dict, input_cds_dict)


	print(("{0}: {1} seconds elapsed. CDS identification complete. " + 
	   "Now assessing putative NMD, NSD status."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	characterize_transcripts_as_nmd_nsd(full_transcript_dict, ptc_distance_nmd_threshold = ptc_dist)


	print(("{0}: {1} seconds elapsed. NMD, NSD status assessment complete. " + 
	   "Now writing CDS-inserted GTF."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	output_cds_inserted_gtf(full_transcript_dict, output_directory)


	print(("{0}: {1} seconds elapsed. CDS-inserted GTF output complete. " + 
	   "Now writing transcript characteristics table."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	output_transcript_table(full_transcript_dict, output_directory)


	print(("{0}: {1} seconds elapsed. Transcript characteristics table output " + 
	   "complete. Now writing putative protein sequences."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	output_aa_sequence(full_transcript_dict, output_directory)


	print(("{0}: {1} seconds elapsed. Putative protein sequence output " + 
	   "complete. Now writing CDS and UTR exon bed files."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))


	output_cds_utr_exons(full_transcript_dict, output_directory)


	if make_bigBed:


		print(("{0}: {1} seconds elapsed. CDS and UTR bed files written. " + 
		   "Now creating bigGenePred files."
		   ).format(str(datetime.now().replace(microsecond = 0)), 
		            str(round(time.time() - start_time, 1))))


		if chrNameLength_path is None or bigGenePred_as_path is None:
			print(("Both --chrNameLength_path and " + 
				   "--bigGenePred_as_path must be set " + 
				   "in order to generate bigBed file.  Skipping . . . "))
		else:

			make_bb_file("cds_inserted_no_ptc", 
						 bigGenePred_as_path, 
						 chrNameLength_path, 
						 output_directory, 
						 gtfToGenePred_exe = gtfToGenePred_path, 
						 genePredToBigGenePred_exe = genePredToBigGenePred_path, 
						 bedToBigBed_exe = bedToBigBed_path)

			make_bb_file("cds_inserted_ptc", 
						 bigGenePred_as_path, 
						 chrNameLength_path, 
						 output_directory, 
						 gtfToGenePred_exe = gtfToGenePred_path, 
						 genePredToBigGenePred_exe = genePredToBigGenePred_path, 
						 bedToBigBed_exe = bedToBigBed_path)

			make_bb_file("cds_inserted_nonstop", 
						 bigGenePred_as_path, 
					 	 chrNameLength_path, 
					 	 output_directory, 
					 	 gtfToGenePred_exe = gtfToGenePred_path, 
					 	 genePredToBigGenePred_exe = genePredToBigGenePred_path, 
					 	 bedToBigBed_exe = bedToBigBed_path)


		print(("{0}: {1} seconds elapsed. bigGenePred files written. " + 
	   "Now writing transcript dict pkl."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))

	else:

		print(("{0}: {1} seconds elapsed. CDS and UTR bed files written. " + 
		   "Now writing transcript dict pkl."
		   ).format(str(datetime.now().replace(microsecond = 0)), 
		            str(round(time.time() - start_time, 1))))


	pkl.dump(full_transcript_dict, 
			 open(output_directory + "/cds_insertion_transcript_dict.pkl", "wb"))


	print(("{0}: {1} seconds elapsed. Transcript dict pkl written. " + 
	   "cds_insertion complete."
	   ).format(str(datetime.now().replace(microsecond = 0)), 
	            str(round(time.time() - start_time, 1))))

	print("Ceci n'est pas un algorithme bioinformatique.")


	return full_transcript_dict


if __name__ == '__main__':

	main(sys.argv[1:])


	




