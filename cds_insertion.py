#!/usr/bin/python

import sys
from operator import itemgetter

'''
Sample CDS input: 
chr1 874420 874509 + CDS SAMD11
chr1 874655 874671 + CDS SAMD11
chr1 874420 874509 + CDS SAMD11
chr1 874652 874671 + CDS SAMD11


Sample GENCODE transcript input:
chr1	11869	12227	+	ENST00000456328.2	DDX11L1
chr1	12613	12721	+	ENST00000456328.2	DDX11L1
chr1	13221	14409	+	ENST00000456328.2	DDX11L1
chr1	11872	12227	+	ENST00000515242.2	DDX11L1
chr1	12613	12721	+	ENST00000515242.2	DDX11L1
chr1	13225	14412	+	ENST00000515242.2	DDX11L1

'''

###TODO add sort method for each exon list as it is inducted
###TODO allow multiple input CCDS per gene
###TODO double-check numbering
###TODO output genePred file as well as file with tags
###

def add_exon_to_cds_dict(entry,exon_dict):
	##make field 4 CDS/stop_codon -> keep collecting entries until reaching new stop codon
	if entry[5] in exon_dict:

		if entry[4] == "stop_codon":
			exon_dict[entry[5]]["stop_codons"] += 1

		CDS_count = exon_dict[entry[5]]["stop_codons"] + 1

		if CDS_count in exon_dict[entry[5]]["cds"]:

			exon_dict[entry[5]]["cds"][CDS_count].append([entry[1],entry[2]])

		else:
			exon_dict[entry[5]]["cds"][CDS_count] = [[entry[1],entry[2]]]
		
	else:
		exon_dict[entry[5]] = {"cds":{1:[[entry[1],entry[2]]]}, "stop_codons": 0}
			
def add_exon_to_transcript_dict(entry,exon_dict):
	
	if entry[5] in exon_dict:
		if entry[4] in exon_dict[entry[5]]:		
			exon_dict[entry[5]][entry[4]]["exons"].append([entry[1],entry[2]])
		else:
			exon_dict[entry[5]][entry[4]] = {"chrom":entry[0],"strand":entry[3],"sequence": entry[6],"exons":[[entry[1],entry[2]]]}			
	else:
		exon_dict[entry[5]] = {entry[4]:{"chrom":entry[0],"strand":entry[3], "sequence":entry[6] ,"exons":[[entry[1],entry[2]]]}}

def get_junctions(exon_list):
	
	junction_list = []
	for i in range(0,len(exon_list) - 1):
		junction_list.append([exon_list[i][1],exon_list[i+1][0]])
	return junction_list

def junction_overlap(transcript_start, transcript_end, cds_start, cds_end, cds_jl, transcript_jl, transcript_seq, transcript_exons, strand):
	##TODO replace with subset operation

	
	cds_subset = [i for i in cds_jl if i in transcript_jl]

	if cds_subset == cds_jl:
		
		first_cds_junction_index = transcript_jl.index(cds_subset[0])
		last_cds_junction_index = transcript_jl.index(cds_subset[-1])
		transcript_junction_subset = transcript_jl[first_cds_junction_index:last_cds_junction_index+1]
		

		if transcript_junction_subset == cds_subset and cds_start in range(transcript_start, transcript_end + 1) and cds_end in range(transcript_start, transcript_end + 1):
			
			return [cds_start, cds_end]  ##If this is true, then we just output a line in genePred format, so that genePredToGtf will perform the UTR annotation for us

	##TODO fix this to go exon-by-exon
	else:
		for exon in transcript_exons:
			if cds_start in range(int(exon[0]), int(exon[1]) + 1): ##Really this should go exon by exon 
				return translate_ORF(transcript_start, cds_start, transcript_seq, transcript_exons, strand)
		

		else:
			return "uncontained_start"


def get_feature_start_end(strand, feature_exons):
	if strand == "+":
		return [feature_exons[0][0],feature_exons[-1][1]]

	elif strand == "-":
		return [feature_exons[-1][1], feature_exons[0][0]]


def genome_to_transcript_coords(position, transcript_exons, strand = "+", direction = "TG"): ##Where "TG" is transcript -> genome and "GT" is genome -> transcript
	
	position = int(position)
	
	if strand == "-":
		transcript_exons.reverse()

	exon_lengths = [0] ## zero is just to prevent list index out of range on first iteration of loop below

	for i in transcript_exons:
		#print i
		prev_length = exon_lengths[-1]
		exon_lengths += [prev_length + (int(i[1]) - int(i[0]))]

	#print exon_lengths
	#print position

	if direction == "TG":		
		for i in range(0,len(exon_lengths)-1):
			if position < exon_lengths[i+1]:				
				return int(transcript_exons[i][0]) + position - exon_lengths[i] 
		else:
			print "GARBAGE"
				
	elif direction == "GT":
		for i in transcript_exons:
			if position in range(int(i[0]),int(i[1]) + 1):
				return position - int(i[0]) + int(exon_lengths[transcript_exons.index(i)])		



def translate_ORF(transcript_start, cds_start, transcript_seq, transcript_exons, strand):
	
	##cds_start must be adjusted to reflect the fact that we are now dealing (through the sequence) with the spliced transcript (rather than genomic coords)

	adj_cds_start = genome_to_transcript_coords(cds_start, transcript_exons, strand, "GT")
	
	
	

	CDS_seq_left_bound = transcript_seq[adj_cds_start:].upper()
	codon_list = []

	for i in range(0,len(CDS_seq_left_bound) - 2,3):

		codon = CDS_seq_left_bound[i:i+3]
		if codon not in stop_codons:
			codon_list += [codon]
		else:
			codon_list += [codon]
			adj_cds_stop = adj_cds_start + 3*len(codon_list)
			print codon_list
			return [cds_start, genome_to_transcript_coords(adj_cds_stop, transcript_exons, strand, "TG")]
	else:
		return "uncontained_stop"


	
	##Goal is to return the position of the new stop codon

def is_NMD(cds_end, transcript_exons, strand):
	
	if strand == "-":
		transcript_exons.reverse()
	
	for i in transcript_exons:
		if cds_end in range(int(i[0]),int(i[1])+1):

			if transcript_exons.index(i) == (len(transcript_exons) - 1):
				return "no_PTC"			
			else:
				distance = int(i[1]) + 1 - cds_end
				return ["NMD_candidate", distance]

	else:
		return "Error: stop codon not contained in spliced transcript"


if __name__ == '__main__':

	stop_codons = ["TAG","TAA","TGA"]

	full_transcript_exon_file = open('exon_test.tsv','r')


	full_transcript_dict = {}
	cds_dict = {}
	all_events = []


	for line in full_transcript_exon_file:
		entry = line.split()
		add_exon_to_transcript_dict(entry,full_transcript_dict)
	
	full_transcript_exon_file.close()

	cds_file = open('cds_test2_missing_exon.tsv','r')

	for line in cds_file:
		entry = line.split()
	
		add_exon_to_cds_dict(entry,cds_dict)
	
	cds_file.close()

	for gene in cds_dict:
		for cds in cds_dict[gene]["cds"]:
			if gene in full_transcript_dict:
				for transcript in full_transcript_dict[gene]:
				
					cds_exons = sorted(cds_dict[gene]["cds"][cds], key=itemgetter(0))
					transcript_exons = sorted(full_transcript_dict[gene][transcript]["exons"], key=itemgetter(0))
					
					transcript_jl = get_junctions(transcript_exons)
					cds_jl = get_junctions(cds_exons)

					strand = full_transcript_dict[gene][transcript]["strand"]
			
		
					transcript_boundaries = get_feature_start_end(strand, transcript_exons)
					transcript_start = int(transcript_boundaries[0])
					transcript_end = int(transcript_boundaries[1])

					cds_boundaries = get_feature_start_end(strand, cds_exons)
					cds_start = int(cds_boundaries[0])
					cds_end = int(cds_boundaries[1])

			
					transcript_seq = full_transcript_dict[gene][transcript]["sequence"]

					start_stop_codons = junction_overlap(transcript_start, transcript_end, cds_start, cds_end, cds_jl, transcript_jl, transcript_seq, transcript_exons, strand)

					if len(start_stop_codons) == 2:
						cds_end = start_stop_codons[1]

					NMD_status = is_NMD(cds_end, transcript_exons, strand)

					all_events.append([transcript,transcript_start, transcript_end,start_stop_codons, NMD_status])


				
	for i in all_events:
		print i
				
		
			
						








