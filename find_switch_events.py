import argparse
import sys
import splice_lib


def string_bool(string):

	if string == "True":

		return True

	elif string == "False":

		return False

	else:

		sys.exit("Unsupported string passed to string_bool function.  Can only use 'True' or 'False'.  Exiting  . . . ")


def event_nmd_nsd_status(standard_event_dict, standard_transcript_dict):

	for event in standard_event_dict:

		included_always_nmd = []
		included_sometimes_nmd = []
		included_always_nonstop = []
		included_sometimes_nonstop = []

		excluded_always_nmd = []
		excluded_sometimes_nmd = []
		excluded_always_nonstop = []
		excluded_sometimes_nonstop = []

		for transcript in standard_event_dict[event]["included_form_transcripts"]:

			included_always_nmd.append(string_bool(standard_transcript_dict[transcript]["always_nmd"]))
			included_sometimes_nmd.append(string_bool(standard_transcript_dict[transcript]["sometimes_nmd"]))
			included_always_nonstop.append(string_bool(standard_transcript_dict[transcript]["always_nonstop"]))
			included_sometimes_nonstop.append(string_bool(standard_transcript_dict[transcript]["sometimes_nonstop"]))

		for transcript in standard_event_dict[event]["excluded_form_transcripts"]:

			excluded_always_nmd.append(string_bool(standard_transcript_dict[transcript]["always_nmd"]))
			excluded_sometimes_nmd.append(string_bool(standard_transcript_dict[transcript]["sometimes_nmd"]))
			excluded_always_nonstop.append(string_bool(standard_transcript_dict[transcript]["always_nonstop"]))
			excluded_sometimes_nonstop.append(string_bool(standard_transcript_dict[transcript]["sometimes_nonstop"]))

		always_nmd = False
		sometimes_nmd = False
		ambiguous_nmd = False
		always_nonstop = False
		sometimes_nonstop = False
		ambiguous_nonstop = False

		nmd_form = "NA"
		nonstop_form = "NA"


		if all(included_always_nmd) and not any(excluded_sometimes_nmd):

			always_nmd = True
			sometimes_nmd = True
			nmd_form = "included"

		elif all(excluded_always_nmd) and not any(included_sometimes_nmd):

			always_nmd = True
			sometimes_nmd = True
			nmd_form = "excluded"

		elif any(included_sometimes_nmd) and not any(excluded_sometimes_nmd):

			sometimes_nmd = True
			nmd_form = "included"

		elif any(excluded_sometimes_nmd) and not any(included_sometimes_nmd):

			sometimes_nmd = True
			nmd_form = "excluded"

		elif any(included_sometimes_nmd) and any(excluded_sometimes_nmd):

			ambiguous_nmd = True


		if all(included_always_nonstop) and not any(excluded_sometimes_nonstop):

			always_nonstop = True
			sometimes_nonstop = True
			nonstop_form = "included"

		elif all(excluded_always_nonstop) and not any(included_sometimes_nonstop):

			always_nonstop = True
			sometimes_nonstop = True
			nonstop_form = "excluded"

		elif any(included_sometimes_nonstop) and not any(excluded_sometimes_nonstop):

			sometimes_nonstop = True
			nonstop_form = "included"

		elif any(excluded_sometimes_nonstop) and not any(included_sometimes_nonstop):

			sometimes_nonstop = True
			nonstop_form = "excluded"

		elif any(included_sometimes_nonstop) and any(excluded_sometimes_nonstop):

			ambiguous_nonstop = True

		standard_event_dict[event]["always_nmd"] = always_nmd
		standard_event_dict[event]["sometimes_nmd"] = sometimes_nmd
		standard_event_dict[event]["ambiguous_nmd"] = ambiguous_nmd
		standard_event_dict[event]["always_nonstop"] = always_nonstop
		standard_event_dict[event]["sometimes_nonstop"] = sometimes_nonstop
		standard_event_dict[event]["ambiguous_nonstop"] = ambiguous_nonstop

		standard_event_dict[event]["nmd_form"] = nmd_form
		standard_event_dict[event]["nonstop_form"] = nonstop_form

		ptc_overlap = []

		if always_nmd:

			for transcript in standard_event_dict[event][nmd_form + "_form_transcripts"]:

				for cds in standard_transcript_dict[transcript]["CDS"]:

					if splice_lib.position_contained(standard_event_dict[event][nmd_form + "_alt_regions"], standard_transcript_dict[transcript]["CDS"][cds]["stop_codon"][0][0])[0] or splice_lib.position_contained(standard_event_dict[event][nmd_form + "_alt_regions"], standard_transcript_dict[transcript]["CDS"][cds]["stop_codon"][-1][-1])[0]:

						ptc_overlap.append(True)

					else:

						ptc_overlap.append(False)

		if len(ptc_overlap) > 0:

			if all(ptc_overlap):

				standard_event_dict[event]["ptc_overlap"] = True

			else:

				standard_event_dict[event]["ptc_overlap"] = False

		elif always_nmd:

			sys.exit("There is a problem in event_nmd_nsd_status(): ptc_overlap determination was not successful for an always_nmd event.")

		else:

			standard_event_dict[event]["ptc_overlap"] = "NA"

def output_table(standard_event_dict, outdir):

	output_table = open(outdir + "/event_nmd_nsd_status.tsv", 'w')

	for event in standard_event_dict:

		nmd_field = "NA"
		nmd_form = "NA"
		ptc_overlap = "NA"

		nonstop_field = "NA"
		nonstop_form = "NA"

		if standard_event_dict[event]["always_nmd"]:

			nmd_field = "always_nmd"
			nmd_form = standard_event_dict[event]["nmd_form"]

			if standard_event_dict[event]["ptc_overlap"]:

				ptc_overlap = "TRUE"

			else:

				ptc_overlap = "FALSE"

		elif standard_event_dict[event]["sometimes_nmd"]:

			nmd_field = "sometimes_nmd"
			nmd_form = standard_event_dict[event]["nmd_form"]

		elif standard_event_dict[event]["ambiguous_nmd"]:

			nmd_field = "ambiguous_nmd"


		if standard_event_dict[event]["always_nonstop"]:

			nonstop_field = "always_nonstop"
			nonstop_form = standard_event_dict[event]["nonstop_form"]

		elif standard_event_dict[event]["sometimes_nonstop"]:

			nonstop_field = "sometimes_nonstop"
			nonstop_form = standard_event_dict[event]["nonstop_form"]

		elif standard_event_dict[event]["ambiguous_nonstop"]:

			nonstop_field = "ambiguous_nonstop"

		output_table.write("\t".join([event, nmd_field, nmd_form, ptc_overlap, nonstop_field, nonstop_form]) + "\n")

	output_table.close()

def add_transcripts_to_event_dict(ioe_file, standard_event_dict):

	'''
		Creates pseudo event and transcript dicts to run event_nmd_nsd_status and output_table functions without having passed the full dicts as input to main.  Useful if running script as standalone.
	'''

	with open(ioe_file, 'r') as file:

		next(file)

		for line in file:

			entry = line.strip().split()
			event = entry[2].split(";")[1]

			included_form_transcripts = entry[3].split(",")
			all_transcripts = entry[4].split(",")
			excluded_form_transcripts = list(set(all_transcripts) - set(included_form_transcripts))

			standard_event_dict[event]["included_form_transcripts"] = list(included_form_transcripts) 
			standard_event_dict[event]["excluded_form_transcripts"] = list(excluded_form_transcripts)



def main(args, standard_transcript_dict = None, standard_event_dict = None):

	parser = argparse.ArgumentParser()

	parser.add_argument("--suppress_output", action = "store_true", help = "If set, not file will be output")
	parser.add_argument("--outdir", type = str, help = "Output directory")
	parser.add_argument("--ioe_file", type = str, help = "Event ioe file")
	#parser.add_argument("--transcript_table", type = str, help = "Table of transcript characteristics output by cds_insertion.py")
	parser.add_argument("--event_gtf", type = str, help = "Event gtf file - required for certain features to work e.g. PTC overlap determination if standard_event_dict not supplied as argument to main()")
	parser.add_argument("--transcript_dict_pkl", type = str, help = "Pickled transcript dict from cds_insertion.py.  Required for certain features to work e.g. PTC overlap determinatino if standard_transcript_dict not supplied as argument to main()")

	args = parser.parse_args(args)


	suppress_output = args.suppress_output
	outdir = args.outdir
	ioe_file = args.ioe_file
	#transcript_table = args.transcript_table
	event_gtf = args.event_gtf
	transcript_dict_pkl = args.transcript_dict_pkl 

	if standard_event_dict is None and event_gtf is not None:

		standard_event_dict = splice_lib.generate_standard_event_dict(event_gtf)
		add_transcripts_to_event_dict(ioe_file, standard_event_dict)

		splice_lib.find_exons_unique_to_form(standard_event_dict)

		for event in standard_event_dict:

			standard_event_dict[event]["included_alt_regions"], standard_event_dict[event]["excluded_alt_regions"] = splice_lib.get_alt_regions(standard_event_dict[event]["included_exons"], standard_event_dict[event]["excluded_exons"], standard_event_dict[event]["included_unique_exons"], standard_event_dict[event]["excluded_unique_exons"])




	else:

		sys.exit("Please supply either standard_event_dict as argument to find_switch_events main, or supply path to --event_gtf.  Supplying both is not an option (it is confusing!)")

	if standard_transcript_dict is None and transcript_dict_pkl is not None:

		import cPickle as pkl
		standard_transcript_dict = pkl.load(open(transcript_dict_pkl, "rb"))

	else:

		sys.exit("Please supply either standard_transcript_dict as arg to main() or supply path to --transcript_dict_pkl (not both)")


	if not suppress_output and outdir is None:

		sys.exit("Output directory required if suppress output not set in find_switch_events.py")

	event_nmd_nsd_status(standard_event_dict, standard_transcript_dict)
	output_table(standard_event_dict, outdir)
	

if __name__ == '__main__':

	main(sys.argv[1:])




