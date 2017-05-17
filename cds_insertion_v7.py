#!/usr/bin/python

import sys
from operator import itemgetter
import argparse
import splice_lib
import copy

'''
Sample CDS input: 
chr1 874420 874509 + CDS SAMD11
chr1 874655 874671 + CDS SAMD11
chr1 874420 874509 + CDS SAMD11
chr1 874652 874671 + CDS SAMD11


Sample GENCODE transcript input:
chr1	11869	12227	+	ENST00000456328.2	DDX11L1	sequence
chr1	12613	12721	+	ENST00000456328.2	DDX11L1	sequence
chr1	13221	14409	+	ENST00000456328.2	DDX11L1	sequence	
chr1	11872	12227	+	ENST00000515242.2	DDX11L1	sequence
chr1	12613	12721	+	ENST00000515242.2	DDX11L1	sequence
chr1	13225	14412	+	ENST00000515242.2	DDX11L1	sequence

'''

###TODO rewrite 



###TODO - fix 3'-UTR length (currently off by 3 due to including stop codon)
###TODO - fix minus strand PTC distance in is_NMD
###TODO - output identity of PTC exon
###TODO double-check numbering
###TODO output genePred file as well as file with tags


stop_codons = ["TAG", "TAA", "TGA"]
			
def generate_transcript_dict(transcript_exon_filename):

	'''
		Returns a dictionary indexed by transcript ID from transcript exon input file
	'''

	transcript_exon_file = open(transcript_exon_filename, 'r')

	full_transcript_dict = {}

	for line in transcript_exon_file:

		entry = line.split()
		chrom = entry[0].strip()
		start = entry[1].strip()
		end = entry[2].strip()
		strand = entry[3].strip()
		transcript_id = entry[4].strip()
		gene_id = entry[5].strip()
		transcript_seq  = entry[6].strip()
	
		if transcript_id in full_transcript_dict:		
			full_transcript_dict[transcript_id]["exons"].append([start,end])

		else:
			full_transcript_dict[transcript_id] = {"chrom":chrom,
										"strand":strand,
										"gene": gene_id,
										"sequence": transcript_seq,
										"exons":[[start,end]],
										"length": len(transcript_seq),
										"included_form_events": [],
										"excluded_form_events": [],
										"CDS": {},
										"nogo": {}
										}			

	transcript_exon_file.close()
	return full_transcript_dict


def index_transcripts_by_junctions(full_transcript_dict):
	'''
		Returns dictionary where splice junctions (specified by chrom, strand, and inner coordinates of flanking exons) are keys which point to lists of matching transcripts

		Takes the transcript dictionary as input.

		This facilitates the matching of CDS without needing prior gene information.

	'''

	junction_dict = {}

	for transcript in full_transcript_dict:

		transcript_exons = sorted(full_transcript_dict[transcript]["exons"][:], key=itemgetter(0))
		full_transcript_dict[transcript]["exons"] = transcript_exons

		transcript_jl = splice_lib.get_junctions(transcript_exons)
		full_transcript_dict[transcript]["junctions"] = transcript_jl

		for junction in transcript_jl:

			junction_key = full_transcript_dict[transcript]["chrom"] + "_" + "_".join(junction) + "_" + full_transcript_dict[transcript]["strand"]

			if junction_key not in junction_dict:

				junction_dict[junction_key] = [transcript]

			else:
				junction_dict[junction_key].append(transcript)

	return junction_dict


def add_exon_to_cds_dict(cds_filename, junction_transcript_dict, full_transcript_dict):
	'''
		Returns CDS dictionary from CDS exon input file (indexed by gene - this is OK since CDS input file is generated using an established annotation)

		As each CDS is collected, it is matched via splice junctions to potential matching transcripts via the junction-indexed transcript dict.  

		The match is assessed by calling cds_junction_search (which calls assign cds_to_transcript).

	'''
	##make field 4 CDS/stop_codon -> keep collecting entries until reaching new stop codon
	#chr1 874420 874509 + CDS SAMD11

	temp_exon_list = []
	cds_dict = {}

	cds_file = open(cds_filename, 'r')

	for line in cds_file:

		entry = line.strip().split()

		chrom = entry[0].strip()
		start = entry[1].strip()
		end = entry[2].strip()
		strand = entry[3].strip()
		feature = entry[4].strip()
		gene = entry[5].strip()

		if feature != "stop_codon":

			temp_exon_list.append([start, end])

		elif feature == "stop_codon" and len(temp_exon_list) != 0: ##Some genes have extraneous stop codons thrown in by gencode

			exon_list = sorted(temp_exon_list, key = itemgetter(0))


			start_end = splice_lib.get_feature_start_end(strand, exon_list)
			flat_temp_exon_list = [i for j in temp_exon_list[:] for i in j]

			test_cds_key = chrom + "_" + "_".join(flat_temp_exon_list) + "_" + strand



			if test_cds_key not in cds_dict:

				cds_dict[test_cds_key] = {
											"cds_start": start_end[0],
											"cds_end": start_end[1],
											"exons": sorted(exon_list, key=itemgetter(0)),
											"junctions": splice_lib.get_junctions(exon_list),
											"chrom": chrom,
											"strand": strand,
											"gene": gene,
											"transcripts": []

										}

				#####Call add cds to full transcript dict here
				print "Looking for matching transcripts . . . "
				cds_junction_search(cds_dict[test_cds_key], junction_transcript_dict, full_transcript_dict)

			temp_exon_list = []

	cds_file.close()

	return cds_dict


def cds_junction_search(cds, junction_transcript_dict, full_transcript_dict):

	'''
		Takes as input one cds exon dict entry, the junction-indexed transcript dict, and the original transcript id indexed transcript dict

		Tries to match each cds junction to transcripts via the junction-indexed transcript dict.  If a match is found, assign_cds_to_transcript is called,
		which then checks for a full CDS match.
	'''

	##Step one - cycle through junctions to find a match
	for junction in cds["junctions"]:

		print "Checking for matching junction . . . "

		test_key = cds["chrom"] + "_" + "_".join(junction) + "_" + cds["strand"]

		if test_key in junction_transcript_dict:

			print "Matching junction found . . . "

			for transcript in junction_transcript_dict[test_key]:

				print "Attempting to add cds to",transcript
				assign_cds_to_transcript(cds, transcript, full_transcript_dict)


def assign_cds_to_transcript(cds, transcript, full_transcript_dict):

	'''
		Takes as input an entry from cds exon dict, a transcript id, and the full transcript dict.  Checks to see if the CDS is a full match for the transcript, and (if it is) adds it to
		the full transcript dict.
	'''

	gene = cds["gene"]

	cds_exons = cds["exons"]
	transcript_exons = full_transcript_dict[transcript]["exons"]

	transcript_jl = full_transcript_dict[transcript]["junctions"]
	cds_jl = cds["junctions"]


	strand = full_transcript_dict[transcript]["strand"]
	chrom = full_transcript_dict[transcript]["chrom"]


	transcript_boundaries = splice_lib.get_feature_start_end(strand, transcript_exons)
	transcript_start = int(transcript_boundaries[0])
	transcript_end = int(transcript_boundaries[1])

	cds_start = int(cds["cds_start"])
	cds_end = int(cds["cds_end"])

	transcript_seq = full_transcript_dict[transcript]["sequence"]


	start_codon_contained = splice_lib.position_contained(transcript_exons, cds_start)
		##Above returns a list where the first element is boolean and the second is either the containing exon (if True) or None (if False)
	start_is_contained = start_codon_contained[0]
	start_containing_exon = start_codon_contained[1]

	stop_codon_contained = splice_lib.position_contained(transcript_exons, cds_end)
	stop_is_contained = stop_codon_contained[0]
	stop_containing_exon = stop_codon_contained[1]

	junctions_do_overlap = splice_lib.junction_overlap(cds_jl, transcript_jl)
		##Just returns a boolean

		##start_is_contained, stop_is_contained, and junction_overlap_status must all evaluate to True in order to add the original CCDS to the transcript
	if start_is_contained:

		valid_cds_start = cds_start

		valid_adj_cds_start = splice_lib.genome_to_transcript_coords(cds_start, strand, transcript_exons, "GT")

		if (stop_is_contained and junctions_do_overlap) or (stop_is_contained and (start_containing_exon == stop_containing_exon)):

			valid_cds_end = cds_end

			cds_utr_exons = splice_lib.get_cds_utr_exons(strand, transcript_exons, valid_cds_start, valid_cds_end, start_containing_exon, stop_containing_exon)

			flat_cds_list = [i for j in cds_utr_exons[1] for i in j]
			flat_cds_list = map(str, flat_cds_list)

			cds_key = chrom + "_" + "_".join(flat_cds_list) + "_" + strand

			if cds_key not in full_transcript_dict[transcript]["CDS"]:

				print "Adding CDS to transcript!"

				valid_adj_cds_end = splice_lib.genome_to_transcript_coords(valid_cds_end, strand, transcript_exons, "GT")

				NMD_status = is_NMD(valid_cds_end, transcript_exons, strand)



				full_transcript_dict[transcript]["CDS"][cds_key] = {
																							"start": valid_cds_start,
																							"end": valid_cds_end,
																							"adj_start": valid_adj_cds_start,
																							"adj_end": valid_adj_cds_end,
																							"CCDS": True,
																							"exons": cds_exons,
																							"gene": gene,
																							"PTC_exon": NMD_status[0],
																							"PTC_distance": NMD_status[1]
																							}

			else: ##This is necessary in case translate ORF from another CCDS (also from the same gene) resulted in addition of this same CDS (but without awareness of its CCDS status)
				full_transcript_dict[transcript]["CDS"][cds_key]["CCDS"] = True

				print "Marked existing CDS as CCDS"

		else:

			translation_attempt = splice_lib.translate_ORF(transcript_seq, stop_codons, valid_adj_cds_start)
			###Note that there are two possible outcomes: 1) a stop codon is found 2) no stop codon is found, and we keep both of them (potentially we may find examples of non-stop decay substrates)

			translation_attempt_is_successful = translation_attempt[0]

			if translation_attempt_is_successful:

				print "Translation attempt successful!"
				valid_adj_cds_end = translation_attempt[1]

				valid_cds_end = splice_lib.genome_to_transcript_coords(valid_adj_cds_end, strand, transcript_exons, "TG")

				
				stop_containing_exon = splice_lib.position_contained(transcript_exons, valid_cds_end)[1]

				
				cds_utr_exons = splice_lib.get_cds_utr_exons(strand, transcript_exons, valid_cds_start, valid_cds_end, start_containing_exon, stop_containing_exon)


				flat_cds_list = [i for j in cds_utr_exons[1] for i in j]
				flat_cds_list = map(str, flat_cds_list)

				cds_key = chrom + "_" + "_".join(flat_cds_list) + "_" + strand


				if cds_key not in full_transcript_dict[transcript]["CDS"]: 

					NMD_status = is_NMD(valid_cds_end, transcript_exons, strand)

					full_transcript_dict[transcript]["CDS"][cds_key] = {
																								"start": valid_cds_start,
																								"end": valid_cds_end,
																								"adj_start": valid_adj_cds_start,
																								"adj_end": valid_adj_cds_end,
																								"CCDS": True,
																								"exons": copy.deepcopy(cds_utr_exons[1]),
																								"gene": gene,
																								"PTC_exon": NMD_status[0],
																								"PTC_distance": NMD_status[1]
																								}
			else:

				if valid_cds_start not in full_transcript_dict[transcript]["nogo"]:
					full_transcript_dict[transcript]["nogo"][valid_cds_start] = {
																							"start": valid_cds_start
																							}


def is_NMD(cds_end, transcript_exons, strand):

	'''
		Determines if stop codon is upstream of splice junction (PTC)

	'''

	print "Checking NMD status . . . "
	
	
	if strand == "+":
		is_NMD_transcript_exons = transcript_exons
	elif strand == "-":
		is_NMD_transcript_exons = transcript_exons[::-1]
	
	for i in transcript_exons:

		if cds_end in range(int(i[0]),int(i[1])+1):

			if is_NMD_transcript_exons.index(i) == (len(is_NMD_transcript_exons) - 1):

				print "No evidence of NMD for transcript-CDS combination . . . continuing  . . . "
				return ["no_PTC", "NA"]
			else:
				if strand == "+":
					distance = int(i[1]) + 1 - cds_end
				elif strand == "-":
					distance = cds_end + 1 - int(i[0])

				print "This transcript may be an NMD substrate!"
					
				return [i, str(distance)]

	else:
		print "NMD Error ... printing transcript exons ..."
		print cds_end
		print is_NMD_transcript_exons
		return "Error: stop codon not contained in spliced transcript"


def import_events(event_ioe_filename, full_transcript_dict):

	'''
		Takes as input SUPPA ioe file (path to), full transcript dict

		Adds events to full transcript dict (separated by included vs excluded form)
	'''

	event_ioe_file = open(event_ioe_filename, 'r')

	event_ioe_file.readline() ##skip header

	for line in event_ioe_file:

		entry = line.strip().split()

		event_id = entry[2].strip().split(";")[1]
		included_form_transcripts = entry[3].strip().split(",")
		all_forms_transcripts = entry[4].strip().split(",")

		excluded_form_transcripts = list(set(all_forms_transcripts) - set(included_form_transcripts))

		for transcript in included_form_transcripts:

			if event_id not in full_transcript_dict[transcript]["included_form_events"]:
				full_transcript_dict[transcript]["included_form_events"].append(event_id)

		for transcript in excluded_form_transcripts:

			if event_id not in full_transcript_dict[transcript]["excluded_form_events"]:
				full_transcript_dict[transcript]["excluded_form_events"].append(event_id)

	event_ioe_file.close()



def make_cds_centric_dict(full_transcript_dict):

	'''
		Takes as input full transcript dict, returns a CDS (chrom, strand, full exon list all used to index) indexed dict containing transcripts and events that belong to CDS

	'''

	cds_centric_dict = {}

	for transcript in full_transcript_dict:

		for cds in full_transcript_dict[transcript]["CDS"]:

			if cds not in cds_centric_dict:

				cds_centric_dict[cds] = full_transcript_dict[transcript]["CDS"][cds].copy()
				cds_centric_dict[cds]["transcripts"] = [transcript]
				cds_centric_dict[cds]["included_form_events"] = full_transcript_dict[transcript]["included_form_events"][:]
				cds_centric_dict[cds]["excluded_form_events"] = full_transcript_dict[transcript]["excluded_form_events"][:]

			else:

				if transcript not in cds_centric_dict[cds]["transcripts"]:

					cds_centric_dict[cds]["transcripts"].append(transcript)


				for i in full_transcript_dict[transcript]["included_form_events"]:

					if i not in cds_centric_dict[cds]["included_form_events"]:

						cds_centric_dict[cds]["included_form_events"].append(i)

				for i in full_transcript_dict[transcript]["excluded_form_events"]:

					if i not in cds_centric_dict[cds]["excluded_form_events"]:

						cds_centric_dict[cds]["excluded_form_events"].append(i)

	return cds_centric_dict


def find_utr_events(cds_centric_dict, outdir):
	'''
		Takes as input cds centric dict, filename for output file that fill contain CDS-only events (i.e. both included and excluded events match the same CDS), returns gene-indexed
		dict containing UTR-only events

	'''

	print "Identifying UTR only events . . . "

	not_utr_only_events_at_least_sometimes = []

	utr_events_output_file = open(outdir + "/utr_events_output.tsv", 'w')
	utr_events_output_file.truncate()

	strict_utr_events_output_file = open(outdir + "/strict_utr_events_output.tsv", 'w')
	strict_utr_events_output_file.truncate()

	all_putative_utr_events = {}
	strict_putative_utr_events = {}

	for cds in cds_centric_dict:

		gene = cds_centric_dict[cds]["gene"]

		for i in cds_centric_dict[cds]["included_form_events"]:

			if i in cds_centric_dict[cds]["excluded_form_events"]:

				if gene not in all_putative_utr_events:

					all_putative_utr_events[gene] = [i]

				else:

					if i not in all_putative_utr_events[gene]:

						all_putative_utr_events[gene].append(i) 
			else:

				if i not in not_utr_only_events_at_least_sometimes:

					not_utr_only_events_at_least_sometimes.append(i)

		for i in cds_centric_dict[cds]["excluded_form_events"]:

			if i in cds_centric_dict[cds]["included_form_events"]:

				if gene not in all_putative_utr_events:

					all_putative_utr_events[gene] = [i]

				else:

					if i not in all_putative_utr_events[gene]:

						all_putative_utr_events[gene].append(i)
			else:

				if i not in not_utr_only_events_at_least_sometimes:

					not_utr_only_events_at_least_sometimes.append(i)



	for gene in all_putative_utr_events:

		strict_events = []

		for event in all_putative_utr_events[gene]:

			utr_events_output_file.write(gene + "\t" + event + "\n")

			if event not in not_utr_only_events_at_least_sometimes:

				strict_events.append(event)

		if len(strict_events) != 0:

			strict_putative_utr_events[gene] = strict_events[:]

	utr_events_output_file.close()

	for gene in strict_putative_utr_events:

		for event in strict_putative_utr_events[gene]:

			strict_utr_events_output_file.write(gene + "\t" + event + "\n")

	strict_utr_events_output_file.close()

	return all_putative_utr_events, strict_putative_utr_events

def get_genes_from_cds(cds_centric_dict, outdir):

	'''
		Takes as input cds centric dict, filename for all events file (contains gene, event id)

		Returns gene-indexed dictionary containing all events
		
	'''

	print "Constructing gene-event "

	all_events_output_file = open(outdir + "/all_events_output.tsv", 'w')
	all_events_output_file.truncate()

	all_events = {}

	for cds in cds_centric_dict:

		gene = cds_centric_dict[cds]["gene"]

		for i in cds_centric_dict[cds]["included_form_events"]:

			if gene not in all_events:

					all_events[gene] = [i]

			else:

				if i not in all_events[gene]:

					all_events[gene].append(i) 

		for i in cds_centric_dict[cds]["excluded_form_events"]:

			if gene not in all_events:

					all_events[gene] = [i]

			else:

				if i not in all_events[gene]:

					all_events[gene].append(i)

	for gene in all_events:

		for event in all_events[gene]:

			all_events_output_file.write(gene + "\t" + event + "\n")

	all_events_output_file.close()

	return all_events



def write_NMD_transcripts(full_transcript_dict, outdir):
	'''
		Identifies transcripts that always lead to PTC when translated from CCDS start codon, writes tab-sep output file "always_NMD_transcripts.tsv" in the following way:

		transcript_id	PTC_exons	PTC_distances	included_form_events	excluded_form_events

		Where fields 2-5 (1-base) may contain multiple items in a comma-separated list
	'''

	#sometimes_nmd_events = open(outdir + 'sometimes_NMD_events.tsv', 'w')
	#always_nmd_events = open(outdir + 'always_NMD_events.tsv', 'w')

	always_nmd_transcripts = open(outdir + '/always_NMD_transcripts.tsv', 'w')
	always_nmd_transcripts.truncate()

	never_nmd_transcripts = open(outdir + '/never_NMD_transcripts.tsv', 'w')
	never_nmd_transcripts.truncate()

	orphan_transcripts = open(outdir + '/orphan_transcripts.tsv', 'w')
	orphan_transcripts.truncate()

	always_nmd_switch_events = open(outdir + 'always_NMD_switch_events.tsv', 'w')
	always_nmd_switch_events.truncate()
	sometimes_nmd_switch_events = open(outdir + 'sometimes_NMD_switch_events.tsv', 'w')
	sometimes_nmd_switch_events.truncate()

	ambiguous_sometimes_nmd_events = open(outdir + 'ambiguous_sometimes_nmd_events.tsv', 'w')
	ambiguous_sometimes_nmd_events.truncate()

	included_sometimes_nmd_switch_event_dict = {}
	excluded_sometimes_nmd_switch_event_dict = {}

	included_always_nmd_switch_event_dict = {}
	excluded_always_nmd_switch_event_dict = {}
	included_always_mask = {}
	excluded_always_mask = {}

	#ambiguous_sometimes_nmd_event_dict = {}


	for transcript in full_transcript_dict:

		if len(full_transcript_dict[transcript]["CDS"]) == 0:

			orphan_transcripts.write(transcript + "\n")

			continue

		PTC_distances = []
		PTC_exons = []
		genes = []

		for CDS in full_transcript_dict[transcript]["CDS"]:

			if full_transcript_dict[transcript]["CDS"][CDS]["PTC_distance"] == "NA":

				for event in full_transcript_dict[transcript]["included_form_events"]:

					included_always_mask[event] = None

				for event in full_transcript_dict[transcript]["excluded_form_events"]:

					excluded_always_mask[event] = None

				break
			else:

				for event in full_transcript_dict[transcript]["included_form_events"]:
					included_always_nmd_switch_event_dict[event] = None

				for event in full_transcript_dict[transcript]["excluded_form_events"]:
					excluded_always_nmd_switch_event_dict[event] = None

				if full_transcript_dict[transcript]["CDS"][CDS]["PTC_distance"] not in PTC_distances:

					PTC_distances.append(full_transcript_dict[transcript]["CDS"][CDS]["PTC_distance"])

				if full_transcript_dict[transcript]["CDS"][CDS]["PTC_exon"] not in PTC_exons:

					PTC_exons.append("_".join(full_transcript_dict[transcript]["CDS"][CDS]["PTC_exon"]))

				if full_transcript_dict[transcript]["CDS"][CDS]["gene"] not in genes:

					genes.append(full_transcript_dict[transcript]["CDS"][CDS]["gene"])

		else:

			if len(PTC_distances) == 0 or len(PTC_exons) == 0:

				print "WHOOPS!!!!", transcript




			always_nmd_transcripts.write(transcript + "\t" + ",".join(genes) + "\t" + ",".join(PTC_exons) + "\t" + ",".join(PTC_distances) + "\t" + ",".join(full_transcript_dict[transcript]["included_form_events"]) + "\t" + ",".join(full_transcript_dict[transcript]["excluded_form_events"]) + "\n")



	for transcript in full_transcript_dict:

		if len(full_transcript_dict[transcript]["CDS"]) == 0:

			continue

		genes = []
		start_ends = []

		for CDS in full_transcript_dict[transcript]["CDS"]:

			if full_transcript_dict[transcript]["CDS"][CDS]["PTC_distance"] != "NA":

				break

			else:

				if full_transcript_dict[transcript]["CDS"][CDS]["gene"] not in genes:

					genes.append(full_transcript_dict[transcript]["CDS"][CDS]["gene"])

				start_end = str(full_transcript_dict[transcript]["CDS"][CDS]["start"]) + "_" + str(full_transcript_dict[transcript]["CDS"][CDS]["end"])

				if start_end not in start_ends:

					start_ends.append(start_end)

		else:

			never_nmd_transcripts.write(transcript + "\t" + ",".join(genes) + "\t" + ",".join(start_ends) + "\t" + str(len(start_ends)) + "\n")

	for event in included_always_nmd_switch_event_dict:

		ambiguous_sometimes_nmd_events.write(event + "\t" + "included" + "\n")

		if event not in excluded_always_nmd_switch_event_dict:

			sometimes_nmd_switch_events.write(event + "\t" + "included" + "\n")

			if event not in included_always_mask:
				always_nmd_switch_events.write(event + "\t" + "included" + "\n")

	for event in excluded_always_nmd_switch_event_dict:

		ambiguous_sometimes_nmd_events.write(event + "\t" + "excluded" + "\n")
		if event not in included_always_nmd_switch_event_dict:

			sometimes_nmd_switch_events.write(event + "\t" + "excluded" + "\n")

			if event not in excluded_always_mask:
				always_nmd_switch_events.write(event + "\t" + "excluded" + "\n")

	always_nmd_transcripts.close()
	never_nmd_transcripts.close()
	orphan_transcripts.close()
	always_nmd_switch_events.close()
	sometimes_nmd_switch_events.close()
	ambiguous_sometimes_nmd_events.close()


def write_CDS_gene_transcript_files(cds_centric_dict, outdir):
	'''
		Outputs a tsv file containing transcript ID, cds, and gene, and another with just gene and transcript ID

		Note that transcript IDs may appear multiple times if the transcript could conceivably be translated from multiple CCDS start codons
	'''

	cds_transcript_gene_output = open(outdir + "/cds_transcript_gene.tsv", 'w')
	cds_transcript_gene_output.truncate()

	transcript_gene_output = open(outdir + "/transcript_gene.tsv", 'w')
	transcript_gene_output.truncate()

	transcript_gene = []

	for cds in cds_centric_dict:

		for transcript in cds_centric_dict[cds]["transcripts"]:

			transcript_id = transcript
			gene = cds_centric_dict[cds]["gene"]

			gene_transcript_id = gene + "_" + transcript_id

			if gene_transcript_id not in transcript_gene:

				transcript_gene.append(gene_transcript_id)

			cds_transcript_gene_output.write(cds + "\t" + transcript + "\t" + cds_centric_dict[cds]["gene"] + "\n")

	cds_transcript_gene_output.close()

	for i in transcript_gene:

		transcript_gene_output.write("\t".join(i.split("_")) + "\n")

	
	transcript_gene_output.close()











def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--input_exons", type = str, help = "Input exon filename")
	parser.add_argument("--input_cds", type = str, help = "Input CDS filename")
	parser.add_argument("--events", type = str, help = "IOE file from SUPPA")
	parser.add_argument("--outdir", type = str, help = "Path for output files")

	args = parser.parse_args()


	input_exon_filename = args.input_exons
	input_cds_filename = args.input_cds
	suppa_event_ioe_filename = args.events
	output_directory = args.outdir

	full_transcript_dict = generate_transcript_dict(input_exon_filename)



	junction_transcript_dict = index_transcripts_by_junctions(full_transcript_dict)

	input_cds_dict = add_exon_to_cds_dict(input_cds_filename, junction_transcript_dict, full_transcript_dict) ##Note that cds_junction_search and assign_cds_to_transcript will also run as a part of this

	import_events(suppa_event_ioe_filename, full_transcript_dict)

	cds_centric_dict = make_cds_centric_dict(full_transcript_dict)

	find_utr_events(cds_centric_dict, output_directory)	

	get_genes_from_cds(cds_centric_dict, output_directory)

	write_NMD_transcripts(full_transcript_dict, output_directory)

	write_CDS_gene_transcript_files(cds_centric_dict, output_directory)


if __name__ == '__main__':

	main()


	




