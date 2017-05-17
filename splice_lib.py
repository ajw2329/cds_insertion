
import copy
import re
from operator import itemgetter

def get_junctions(exon_list):


	'''
		Takes a list of exons (list of lists with each element of the inner list containing the start, stop position of an exon) and returns a list of the junction positions (essentially the inner coordinates of any pair of exons).  This facilitates comparison of objects that may differ in start/stop positions but contain identical splice junctions.
	'''
	
	junction_list = []
	for i in range(0,len(exon_list) - 1):
		junction_list.append([exon_list[i][1],exon_list[i+1][0]])
	return junction_list


def junction_overlap(subset_jl, superset_jl):
	##TODO replace with subset operation
	#print "I am the junction_overlap function!"

	'''
		Takes as input two lists of junctions - expected subset and expected superset.  Checks whether subset is truly a subset of superset.  Returns boolean.
	'''

 	subset_subset = [i for i in subset_jl if i in superset_jl]


	if subset_subset == subset_jl and len(subset_subset) > 0:
		print "Subset subset contains all of the junctions - checking for connectivity . . . "

		first_subset_junction_index = superset_jl.index(subset_subset[0])
		last_subset_junction_index = superset_jl.index(subset_subset[-1])
		superset_junction_subset = superset_jl[first_subset_junction_index:last_subset_junction_index+1]
	
		if superset_junction_subset == subset_subset:
			print "Exact junction match!"
			return True
		else:
			print "Connectivity interrupted - no match"
			return False
	else:
		print "subset_subset missing junctions - no match"
		return False


def position_contained(exon_list, position):

	'''
		Takes as input list of exon coordinates (list of lists, each element is a list containing start, stop position of exon), position (single genomic coordinate)

		Checks whether position is contained within any of the exons, returns list containing boolean (True if the position was contained in one of the exons), the containing exon entry (or None of the position was not contain)
	'''

	for exon in exon_list:


		if ((int(position) >= int(exon[0])) and (int(position) <= int(exon[1]))):

			containing_exon = exon
			return [True, containing_exon]

	else:

		return [False, None]


def translate_ORF(transcript_seq, stop_codons, valid_adj_cds_start):
	'''
		Takes as input transcript sequence, list of stop codons (e.g. ["TAA", "TAG", "TGA"]), start codon.  Attempts to find in-frame stop-codon by translating the sequence from the start codon.

		Returns a list containing boolean (True if in-frame stop codon is found), stop codon coordinate (None if nothing is found) in transcript-centric coords

	'''

	CDS_seq_left_bound = transcript_seq[valid_adj_cds_start - 1:].upper()
	codon_list = []

	for i in range(0,len(CDS_seq_left_bound) - 2,3):

		codon = CDS_seq_left_bound[i:i+3]
		if codon not in stop_codons:
			codon_list += [codon]
		else:
			print "In-frame stop codon found!"
			#codon_list += [codon]
			#print "".join(codon_list)
			adj_cds_end = valid_adj_cds_start + 3*len(codon_list) - 1 ###Again check this
			return [True, adj_cds_end]

	else:
		print "No in-frame stop codon found"
		return [False, None]


def genome_to_transcript_coords(position, strand, transcript_exons,   direction = "TG"): ##Where "TG" is transcript -> genome and "GT" is genome -> transcript
	'''
		Interchanges transcript coords (distance along the spliced transcript with reference to the 5'-end of the transcript) with genomic coordinates
	'''

	position = int(position)
	
	if strand == "+":
		gtc_transcript_exons = transcript_exons[:]
	elif strand == "-":
		gtc_transcript_exons = transcript_exons[::-1]

	exon_lengths = [0] ## zero is just to prevent list index out of range on first iteration of loop below (cumulative lengths)

	individual_exon_lengths = []

	for i in gtc_transcript_exons:

		prev_length = exon_lengths[-1]
		exon_lengths += [prev_length + (int(i[1]) - int(i[0]) + 1)]
		individual_exon_lengths += [(int(i[1]) - int(i[0]) + 1)]

	#print individual_exon_lengths


	if direction == "TG":	
		
		print "Converting transcript to genomic coords"	
		for i in range(0,len(exon_lengths)-1):
			if position < exon_lengths[i+1]:
				if strand == "+":
					return int(gtc_transcript_exons[i][0]) + position - exon_lengths[i]   ###double check this!!!  ####MAKE IT TRIPLE!!!!!!!!!!!!!!!!!!!
				elif strand == "-":
					return int(gtc_transcript_exons[i][1]) -  position + exon_lengths[i] ###double check this!!!!

		else:
			print "GARBAGE"			
				
	elif direction == "GT":
		print "Converting genomic to transcript coords"

		for i in gtc_transcript_exons:
			if position in range(int(i[0]),int(i[1]) + 1):
				if strand == "+":
					return abs(position - int(i[0])) + int(exon_lengths[gtc_transcript_exons.index(i)]) + 1	
				elif strand == "-":
					return abs(position - int(i[1])) + int(exon_lengths[gtc_transcript_exons.index(i)]) + 1

#####QUICK EXAMPLE
#
#	genomic position: 236
#	start of exon 1: 231 (meaning that the 1st base in the exon is 231)
#	stop of exon 1: 270
#	start of CDS: 236
#	exon_lengths[transcript_exons.index([231,270])] is equal to exon_lengths[0], which is always going to be 0 the way it's written above
#	therefore, we return 236 - 231 + 0 = 5
#	However, the correct 1-base answer would be 6
#	We can fix this by supplanting 0 with 1 in the exon_lengths list after all the lenglth calculations are finished### NEVERMIND - just add 1 going GT and subtract 1 going TG
#
#	Going backards:
#	transcript coordinate = 6 + 231 - 1 = 236
#
#	More complicated: exon1: [231,270]; exon 2: [321,370]
#		exon_lengths = [0, 40, 90]
#		exon_lengths = [1, 40, 90] # after swapping out the 0 #####AGAIN DONT BOTHER WITH THIS - it's fixed by adding/subtracting 1 as below
#		genomic position = 326 ---> transcript position = 46
#		GT calculation: 326 - 321 + 40 + 1 = 5 + 40 + 1 = 46
#		TG calculation: 321 + 46 - 40 - 1 = 321 + 6 - 1 = 326
#
#	OK now what about example on the minus strand:
#		same exon positions (but the order is now swapped)	
#		exon_lengths = [0, 50, 90]
#		transcript position = 46
#		genomic position = 325
#		GT calculation: 370 - 325 + (length term which is 0 since it's in the first exon) + 1
#		TG calculation: 370 - 46 + (length term - again 0 here)
#
#
#

def get_feature_start_end(strand, feature_exons):
	'''
		Extracts the start, end positions of features based on the exons of the feature.  Start is the left-most position of the leftmost exon for features on the "+" strand, and 
		the right-most position of the right-most exon for features on the "-" strand.  End is the reverse for both cases.
	'''

	if strand == "+":
		return [feature_exons[0][0],feature_exons[-1][1]]

	elif strand == "-":
		return [feature_exons[-1][1], feature_exons[0][0]]


def get_cds_utr_exons(strand, feature_exons, cds_start, cds_end, start_container, end_container):
	'''

	'''

	feature_exons = copy.deepcopy(feature_exons)

	if strand == "+":

		five_utr_subset = copy.deepcopy(feature_exons[0:feature_exons.index(start_container) + 1])
		cds_subset = copy.deepcopy(feature_exons[feature_exons.index(start_container):feature_exons.index(end_container) + 1])
		three_utr_subset = copy.deepcopy(feature_exons[feature_exons.index(end_container):])

		five_utr_subset[-1][-1] = str(int(cds_start) - 1)

		cds_subset[0][0] = cds_start
		cds_subset[-1][-1] = cds_end

		three_utr_subset[0][0] = str(int(cds_end) + 1)

	if strand == "-":

		three_utr_subset = copy.deepcopy(feature_exons[0:feature_exons.index(start_container) + 1])
		cds_subset = copy.deepcopy(feature_exons[feature_exons.index(end_container):feature_exons.index(start_container) + 1])
		five_utr_subset = copy.deepcopy(feature_exons[feature_exons.index(end_container):])

		if len(cds_subset) == 0: 
			print feature_exons

		five_utr_subset[0][0] = str(int(cds_start) + 1)


		if len(cds_subset) == 0:
			print feature_exons
			print end_container
			print start_container
			print feature_exons[feature_exons.index(end_container)]

		cds_subset[0][0] = cds_end
		cds_subset[-1][-1] = cds_start

		three_utr_subset[-1][-1] = str(int(cds_end) - 1)

	return [five_utr_subset, cds_subset, three_utr_subset]



def generate_standard_event_dict(event_gtf_filename):
	'''
		Generates a dictionary of events indexed by event ID
		Presently limited to SUPPA formatting
	'''

	standard_event_dict = {}

	with open(event_gtf_filename) as gtf_file:

		for line in gtf_file:

			entry = line.split()

			if entry[2] == "exon":

				start = int(entry[3])
				end = int(entry[4])

				chrom = entry[0].strip()
				strand = entry[6].strip()

				event_id = re.sub('[;"]', '', entry[9])
				form = re.sub('[;"]', '', entry[11]).split(":")[-1]
				event_type = event_id.split(":")[0]

				if event_id not in standard_event_dict:

					standard_event_dict[event_id] = {
						"event_type": event_type,
						"included_exons": [],
						"excluded_exons": [],
						"chrom": chrom,
						"strand": strand,
						"included_count": 0,
						"excluded_count": 0,
						"included_form_transcripts": [],
						"excluded_form_transcripts": []
					}

				if form == "alternative1":

					standard_event_dict[event_id]["included_exons"].append([start, end])

				if form == "alternative2":

					standard_event_dict[event_id]["excluded_exons"].append([start, end])

	return standard_event_dict



def generate_standard_transcript_dict(transcript_gtf_filename):
	'''
		Generates a dictionary of transcripts indexed by transcript ID
		Presently limited to stringtie formatting
	'''

	standard_transcript_dict = {}

	with open(transcript_gtf_filename) as gtf_file:

		for line in gtf_file:

			entry = line.split()

			if entry[2] == "exon":

				start = int(entry[3])
				end = int(entry[4])

				chrom = entry[0].strip()
				strand = entry[6].strip()

				transcript_id = re.sub('[;"]', '', entry[11].strip())

				#print transcript_id


				if transcript_id not in standard_transcript_dict:

					standard_transcript_dict[transcript_id] = {
						"exons": [],
						"chrom": chrom,
						"strand": strand,
						"included_form_events": [],
						"excluded_form_events": []
					}

				standard_transcript_dict[transcript_id]["exons"].append([start, end])


	return standard_transcript_dict


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

		transcript_jl = get_junctions(transcript_exons)
		full_transcript_dict[transcript]["junctions"] = transcript_jl

		#print transcript_jl

		for junction in transcript_jl:

			junction_key = full_transcript_dict[transcript]["chrom"] + "_" + "_".join(map(str, junction)) + "_" + full_transcript_dict[transcript]["strand"]

			if junction_key not in junction_dict:

				junction_dict[junction_key] = [transcript]

			else:
				junction_dict[junction_key].append(transcript)

	return junction_dict

		