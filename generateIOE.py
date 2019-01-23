import splice_lib
import argparse
import copy
from operator import itemgetter
import subprocess


def add_events_to_transcripts(standard_transcript_dict, standard_event_dict, transcript_junction_dict):

	for event in standard_event_dict:


	 	inc_jl = standard_event_dict[event]["included_junctions"] = [map(int, i) for i in splice_lib.get_junctions(standard_event_dict[event]["included_exons"])]
		excl_jl = standard_event_dict[event]["excluded_junctions"] = [map(int, i) for i in splice_lib.get_junctions(standard_event_dict[event]["excluded_exons"])]


	 	for junction in excl_jl:

	 		junction_key = standard_event_dict[event]["chrom"] + "_" + "_".join(map(str, junction)) + "_" + standard_event_dict[event]["strand"]

	 		if junction_key in transcript_junction_dict:

	 			for transcript in transcript_junction_dict[junction_key]:

	 				if splice_lib.junction_overlap(excl_jl, standard_transcript_dict[transcript]["junctions"]):

	 					if event not in standard_transcript_dict[transcript]["excluded_form_events"]:
	 						standard_transcript_dict[transcript]["excluded_form_events"].append(event)

	 					if transcript not in standard_event_dict[event]["excluded_form_transcripts"]:
	 						standard_event_dict[event]["excluded_form_transcripts"].append(transcript)


	 	for junction in inc_jl:

	 		junction_key = standard_event_dict[event]["chrom"] + "_" + "_".join(map(str, junction)) + "_" + standard_event_dict[event]["strand"]

	 		if junction_key in transcript_junction_dict:

	 			for transcript in transcript_junction_dict[junction_key]:

	 				if splice_lib.junction_overlap(inc_jl, standard_transcript_dict[transcript]["junctions"]):

	 					if event not in standard_transcript_dict[transcript]["included_form_events"]:
	 						standard_transcript_dict[transcript]["included_form_events"].append(event)

	 					if transcript not in standard_event_dict[event]["included_form_transcripts"]:	
	 						standard_event_dict[event]["included_form_transcripts"].append(transcript)



def create_exon_indexed_transcript_dict(standard_transcript_dict):

	exon_indexed_transcript_dict = {}

	for transcript in standard_transcript_dict:

		for exon in standard_transcript_dict[transcript]["exons"]:

			exon_key = standard_transcript_dict[transcript]["chrom"] + "_" + str(exon[0]) + "_" + str(exon[1]) + "_" + standard_transcript_dict[transcript]["strand"]

			exon_indexed_transcript_dict.setdefault(exon_key, set()).add(transcript)

	return exon_indexed_transcript_dict


def add_included_RI_MR(standard_transcript_dict, standard_event_dict, exon_indexed_transcript_dict, outdir):

	for event in standard_event_dict:

		if standard_event_dict[event]["event_type"] in ["MR", "RI"]:

			query_key = standard_event_dict[event]["chrom"] + "_" + standard_event_dict[event]["included_exons"][0][0] + "_" + standard_event_dict[event]["included_exons"][-1][-1] + "_" + standard_event_dict[event]["strand"]

			if query_key in exon_indexed_transcript_dict:

				for transcript in exon_indexed_transcript_dict[query_key]:

					if event not in standard_transcript_dict[transcript]["included_form_events"]:

						standard_transcript_dict[transcript]["included_form_events"].append(event)

					if transcript not in standard_event_dict[event]["included_form_transcripts"]:

						standard_event_dict[event]["included_form_transcripts"].append(transcript)



def output_IOE(standard_event_dict, standard_transcript_dict, outdir, prefix):

	event_ioe_file = open(outdir + "/" + prefix + "_events.ioe", 'w')
	transcript_ioe_file = open(outdir + "/" + prefix + "_transcripts.ioe", 'w')


	event_ioe_file.write("\t".join(["seqname","gene_id","event_id", "alternative_transcripts", "total_transcripts"]) + "\n")

	for event in standard_event_dict:

		if len(standard_event_dict[event]["included_form_transcripts"]) > 0 and len(standard_event_dict[event]["excluded_form_transcripts"]) > 0 and set(standard_event_dict[event]["included_form_transcripts"]) != set(standard_event_dict[event]["excluded_form_transcripts"]):

			event_ioe_file.write(standard_event_dict[event]["chrom"] + "\t" + "gene_id" + "\t" + "gene_id;" + event + "\t" + ",".join(standard_event_dict[event]["included_form_transcripts"]) + "\t"  + ",".join(list(set(standard_event_dict[event]["included_form_transcripts"] + standard_event_dict[event]["excluded_form_transcripts"]))) + "\n")

	transcript_ioe_file.write("\t".join(["seqname","gene_id","transcript_id", "included_form_events", "excluded_form_events"]) + "\n")		

	for transcript in standard_transcript_dict:

		if len(standard_transcript_dict[transcript]["included_form_events"]) > 0 and len(standard_transcript_dict[transcript]["excluded_form_events"]) > 0:

			transcript_ioe_file.write(standard_transcript_dict[transcript]["chrom"] + "\t" + "gene_id" + "\t" + "gene_id;" + transcript + "\t" + ",".join(standard_transcript_dict[transcript]["included_form_events"]) + "\t" + ",".join(standard_transcript_dict[transcript]["excluded_form_events"]) + "\n")


def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--transcript_gtf", type = str, help = "Full transcript gtf file")
	parser.add_argument("--event_gtf", type = str, help = "Event gtf")
	parser.add_argument("--outdir", type = str, help = "Path to output directory")
	parser.add_argument("--prefix", type = str, help = "Prefix for output files", default = "generated")

	args = parser.parse_args()

	transcript_gtf = args.transcript_gtf
	event_gtf = args.event_gtf
	outdir = args.outdir
	prefix = args.prefix

	standard_transcript_dict = splice_lib.generate_standard_transcript_dict(transcript_gtf)

	splice_lib.sort_transcript_dict_exons(standard_transcript_dict)

	splice_lib.add_junctions_to_transcript_dict(standard_transcript_dict)

	exon_indexed_transcript_dict = create_exon_indexed_transcript_dict(standard_transcript_dict)

	standard_event_dict = splice_lib.generate_standard_event_dict(event_gtf)

	transcript_junction_dict = splice_lib.index_transcripts_by_junctions(standard_transcript_dict)

	add_events_to_transcripts(standard_transcript_dict, standard_event_dict, transcript_junction_dict)

	output_IOE(standard_event_dict, standard_transcript_dict, outdir, prefix)


if __name__ == '__main__':
	main()