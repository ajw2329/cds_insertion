import argparse
import sys
from splice_lib import splice_lib



###IMPORTANT NOTE: in python2.7, all([]) == True, while any([]) == False (WHY?!!)


def string_bool(string):

	if string == "True":

		return True

	elif string == "False":

		return False

	else:

		sys.exit("Unsupported string passed to " + 
			     "string_bool function.  Can only " + 
			     "use 'True' or 'False'.  Exiting  . . . ")

### write the following:

## event_statuses()
## cds_statuses()
## coding_switch()
## start_codon_switch()
## stop_codon_switch()
## mid_seq_only_switch()
## UTR lengths, CDS lengths, PTC lengths


def nmd_status(event_entry, standard_transcript_dict):


	property_dict = {"included": {"always_nmd": [],
								  "sometimes_nmd": []},
					 "excluded": {"always_nmd": [],
					 			  "sometimes_nmd": []}}

	for form in [ "included", "excluded" ]:

		for transcript in event_entry[form + "_form_transcripts"]:

			property_dict[form]["always_nmd"].append(
				string_bool(standard_transcript_dict[transcript]["always_nmd"]))
			property_dict[form]["sometimes_nmd"].append(
				string_bool(standard_transcript_dict[transcript]["sometimes_nmd"]))


	included_always_nmd = property_dict["included"]["always_nmd"]
	included_sometimes_nmd = property_dict["included"]["sometimes_nmd"]
	excluded_always_nmd = property_dict["excluded"]["always_nmd"]
	excluded_sometimes_nmd = property_dict["excluded"]["sometimes_nmd"]


	event_entry["nmd_status"] = "never"
	event_entry["nmd_form"] = "neither"

	if (all(included_always_nmd) and 
		len(included_always_nmd) > 0 and not 
		any(excluded_sometimes_nmd)):


		event_entry["nmd_status"] = "always"
		event_entry["nmd_form"] = "included"

	elif (all(excluded_always_nmd) and 
		  len(excluded_always_nmd) > 0 and not 
		  any(included_sometimes_nmd)):

		event_entry["nmd_status"] = "always"
		event_entry["nmd_form"] = "excluded"

	elif (any(included_sometimes_nmd) and not 
		  any(excluded_sometimes_nmd)):

		event_entry["nmd_status"] = "sometimes"
		event_entry["nmd_form"] = "included"

	elif (any(excluded_sometimes_nmd) and not 
		  any(included_sometimes_nmd)):

		event_entry["nmd_status"] = "sometimes"
		event_entry["nmd_form"] = "excluded"

	elif (any(included_sometimes_nmd) and 
		  any(excluded_sometimes_nmd)):

		event_entry["nmd_status"] = "sometimes"
		event_entry["nmd_form"] = "both"



	ptc_overlap = []

	if event_entry["nmd_form"] not in [ "both", "neither" ]:

		for transcript in event_entry[event_entry["nmd_form"] + "_form_transcripts"]:

			for cds, cds_entry in standard_transcript_dict[transcript]["CDS"].iteritems():

				if (splice_lib.position_contained(event_entry[event_entry["nmd_form"] + "_alt_regions"], 
												  cds_entry["stop_codon"][0][0])[0] or 
				    splice_lib.position_contained(event_entry[event_entry["nmd_form"] + "_alt_regions"], 
				    	                          cds_entry["stop_codon"][-1][-1])[0]):

					ptc_overlap.append(True)

				else:

					ptc_overlap.append(False)

	if len(ptc_overlap) > 0:

		if all(ptc_overlap):

			event_entry[ "ptc_overlap" ] = "always"

		elif any(ptc_overlap):

			event_entry[ "ptc_overlap" ] = "sometimes"

		else:

			event_entry[ "ptc_overlap" ] = "never"

	elif event_entry[ "nmd_form" ] not in [ "both", "neither" ]:

		sys.exit("There is a problem in event_nmd_nsd_status(): " +
			     "ptc_overlap determination was not successful for an always_nmd event.")

	else:

		event_entry["ptc_overlap"] = "NA"



def nsd_status(event_entry, standard_transcript_dict):


	property_dict = {"included": {"always_nonstop": [],
								  "sometimes_nonstop": []},
					 "excluded": {"always_nonstop": [],
					 			  "sometimes_nonstop": []}}

	for form in [ "included", "excluded" ]:

		for transcript in event_entry[form + "_form_transcripts"]:

			property_dict[form]["always_nonstop"].append(
				string_bool(standard_transcript_dict[transcript]["always_nonstop"]))
			property_dict[form]["sometimes_nonstop"].append(
				string_bool(standard_transcript_dict[transcript]["sometimes_nonstop"]))

	included_always_nonstop = property_dict["included"]["always_nonstop"]
	included_sometimes_nonstop = property_dict["included"]["sometimes_nonstop"]
	excluded_always_nonstop = property_dict["excluded"]["always_nonstop"]
	excluded_sometimes_nonstop = property_dict["excluded"]["sometimes_nonstop"]

	event_entry["nsd_status"] = "never"
	event_entry["nonstop_form"] = "neither"

	if (all(included_always_nonstop) and 
		len(included_always_nonstop) > 0 and not 
		any(excluded_sometimes_nonstop)):

		event_entry["nsd_status"] = "always"
		event_entry["nonstop_form"] = "included"

	elif (all(excluded_always_nonstop) and 
		  len(excluded_always_nonstop) > 0 and not 
		  any(included_sometimes_nonstop)):

		event_entry["nsd_status"] = "always"
		event_entry["nonstop_form"] = "excluded"

	elif (any(included_sometimes_nonstop) and not 
		  any(excluded_sometimes_nonstop)):

		event_entry["nsd_status"] = "sometimes"
		event_entry["nonstop_form"] = "included"

	elif (any(excluded_sometimes_nonstop) and not 
		  any(included_sometimes_nonstop)):

		event_entry["nsd_status"] = "sometimes"
		event_entry["nonstop_form"] = "excluded"

	elif (any(included_sometimes_nonstop) and 
		  any(excluded_sometimes_nonstop)):

		event_entry["nsd_status"] = "sometimes"
		event_entry["nonstop_form"] = "both"



def coding_noncoding_switch(event_entry, standard_transcript_dict):


	cds_counts = {"included": [],
				  "excluded": []}

	for form in ["included", "excluded"]:

		for transcript in event_entry[form + "_form_transcripts"]:

			normal_count = standard_transcript_dict[transcript]["normal_cds_count"]
			nonstop_count = standard_transcript_dict[transcript]["nonstop_cds_count"]

			cds_counts[form].append(normal_count + nonstop_count)


	if (all([i > 0 for i in cds_counts["excluded"]]) and
		all([i > 0 for i in cds_counts["included"]])):

		event_entry["coding_status"] = "always"
		event_entry["coding_switch"] = "never"
		event_entry["coding_form"] = "both"


	elif (any([i > 0 for i in cds_counts["excluded"]]) or
		  any([i > 0 for i in cds_counts["included"]])):

		event_entry["coding_status"] = "sometimes"	


		if (all([i > 0 for i in cds_counts["included"]]) and not
			 any([i > 0 for i in cds_counts["excluded"]])):

			event_entry["coding_switch"] = "always"
			event_entry["coding_form"] = "included"

		elif (all([i > 0 for i in cds_counts["excluded"]]) and not
			  any([i > 0 for i in cds_counts["included"]])):

			event_entry["coding_switch"] = "always"
			event_entry["coding_form"] = "excluded"

		elif (any([i > 0 for i in cds_counts["included"]]) and not
			  any([i > 0 for i in cds_counts["excluded"]])):

			event_entry["coding_switch"] = "sometimes"
			event_entry["coding_form"] = "included"

		elif (any([i > 0 for i in cds_counts["excluded"]]) and not
			  any([i > 0 for i in cds_counts["included"]])):

			event_entry["coding_switch"] = "sometimes"
			event_entry["coding_form"] = "excluded"

		else:

			event_entry["coding_switch"] = "sometimes"
			event_entry["coding_form"] = "both"


	else:

		event_entry["coding_status"] = "never"
		event_entry["coding_switch"] = "never"
		event_entry["coding_form"] = "neither"


def feature_overlap(event_entry, standard_transcript_dict):

	feature_overlap_dict = {"included": set(),
							"excluded": set()}

	feature_overlap_boolean_dict = {"included": {"cds": [],
												 "utr5": [],
												 "utr3": []},
									"excluded": {"cds": [],
												 "utr5": [],
												 "utr3": []}}


	for form in ["included", "excluded"]:

		for alt_region in event_entry[form + "_alt_regions"]:

			for transcript in event_entry[form + "_form_transcripts"]:

				for cds_entry in standard_transcript_dict[transcript]["CDS"].itervalues():

					if splice_lib.exon_has_overlap(cds_entry["exons"], alt_region):

						feature_overlap_dict[form].add("cds")
						feature_overlap_boolean_dict[form]["cds"].append(True)

					else:

						feature_overlap_boolean_dict[form]["cds"].append(False)


					if splice_lib.exon_has_overlap(cds_entry["five_utr_exons"], alt_region):

						feature_overlap_dict[form].add("utr5")
						feature_overlap_boolean_dict[form]["utr5"].append(True)

					else:

						feature_overlap_boolean_dict[form]["utr5"].append(False)

					if splice_lib.exon_has_overlap(cds_entry["three_utr_exons"], alt_region):

						feature_overlap_dict[form].add("utr3")
						feature_overlap_boolean_dict[form]["utr3"].append(True)

					else:
						feature_overlap_boolean_dict[form]["utr3"].append(False)


	### handle case for empty set

	for form, form_entry in feature_overlap_dict.iteritems():

		if len(form_entry) == 0:

			form_entry.add("NA")


	for form, form_entry in feature_overlap_boolean_dict.iteritems():

		for feature, feature_entry in form_entry.iteritems():

			if len(feature_entry) == 0: ## because any([]) evaluates to True

				feature_overlap_boolean_dict[form][feature] = "never"				

			elif all(feature_entry):

				feature_overlap_boolean_dict[form][feature] = "always"

			elif any(feature_entry):

				feature_overlap_boolean_dict[form][feature] = "sometimes"

			else:

				feature_overlap_boolean_dict[form][feature] = "never"				

	event_entry["feature_overlap"] = feature_overlap_dict
	event_entry["feature_overlap_boolean"] = feature_overlap_boolean_dict




def coding_seq_change(event_entry,
					  standard_transcript_dict):

	cds_set_dict = {"included": set(),
					"excluded": set()}



	for form in ["included", "excluded"]:

		for transcript in event_entry[form + "_form_transcripts"]:

			for cds in standard_transcript_dict[transcript]["CDS"]:

				cds_set_dict[form].add(cds)

	if (len(cds_set_dict["included"]) > 0 and 
		len(cds_set_dict["excluded"]) > 0):

		if (cds_set_dict["included"] ==
			   cds_set_dict["excluded"]):

			event_entry["cds_alteration"] = "never"
			event_entry["form_unique_cds"] = "neither"

		elif len(cds_set_dict["included"] &
			   cds_set_dict["excluded"]) == 0:

			event_entry["cds_alteration"] = "always"
			event_entry["form_unique_cds"] = "both"

		elif (len(cds_set_dict["included"] -
			   cds_set_dict["excluded"]) > 0 and 
			  len(cds_set_dict["excluded"] -
			   cds_set_dict["included"]) == 0):

			event_entry["cds_alteration"] = "sometimes"
			event_entry["form_unique_cds"] = "included"

		elif (len(cds_set_dict["excluded"] -
			   cds_set_dict["included"]) > 0 and 
			  len(cds_set_dict["included"] -
			   cds_set_dict["excluded"]) == 0):

			event_entry["cds_alteration"] = "sometimes"
			event_entry["form_unique_cds"] = "excluded"		

		else:

			event_entry["cds_alteration"] = "sometimes"
			event_entry["form_unique_cds"] = "both"

	else:

		event_entry["cds_alteration"] = "never"
		event_entry["form_unique_cds"] = "neither"




def check_process_property(property_source, 
						   property_target):

	if property_source != "NA":

		property_target.extend(
			map(int,
				property_source.split(",")))

	else:

		property_target.append("NA")



def run_checking_processing(event_entry, 
							standard_transcript_dict, 
							form_properties,
							properties_to_be_processed):

	for form in [ "included", "excluded" ]:

		for transcript in event_entry[form + "_form_transcripts"]:

			transcript_entry = standard_transcript_dict[transcript]

			form_properties[form]["form_transcripts"] += 1

			if transcript_entry["normal_cds_count"] + transcript_entry["nonstop_cds_count"] > 0:

				form_properties[form]["form_coding_transcripts"] += 1


			for property_type in properties_to_be_processed:

				check_process_property(transcript_entry[property_type],
								   form_properties[form][property_type])



def summarize_properties(form_properties, 
						 properties_to_be_processed):

	summarized_properties = {}

	for form, form_entry in form_properties.iteritems():

		summarized_properties[form] = {}

		for prop in properties_to_be_processed:

			prop_sorted = sorted(set(form_entry[prop])) # returns a list since sets are unsorted

			all_prop = ",".join(map(str, prop_sorted))
			no_na_prop = set([i for i in prop_sorted if i != "NA"])

			if len(no_na_prop) > 0:

				mean_prop = float(sum(set(no_na_prop)))/float(len(no_na_prop))

			else:

				mean_prop = "NA"

			summarized_properties[form]["mean_" + prop] = mean_prop
			summarized_properties[form]["all_" + prop] = all_prop

	return summarized_properties



def isoform_property_diffs(summarized_properties):

	isoform_property_diffs_dict = {}

	for prop, prop_entry in summarized_properties["included"].iteritems():

		if "mean" in prop:

			if (summarized_properties["included"][prop] != "NA" and
				summarized_properties["excluded"][prop] != "NA" and
				summarized_properties["excluded"][prop] != 0):

				isoform_property_diffs_dict[prop + "_inc_exc_ratio"] = (float(summarized_properties["included"][prop])/
																		float(summarized_properties["excluded"][prop]))

				isoform_property_diffs_dict[prop + "_inc_exc_diff"] = abs(float(summarized_properties["included"][prop])-
																		float(summarized_properties["excluded"][prop]))

			else:

				isoform_property_diffs_dict[prop + "_inc_exc_ratio"] = "NA"
				isoform_property_diffs_dict[prop + "_inc_exc_diff"] = "NA"

	return isoform_property_diffs_dict



def get_form_properties(event_entry, 
						standard_transcript_dict):

	form_properties = {"included": {"form_transcripts": 0,
									"form_coding_transcripts": 0,
									"five_utr_lengths": [],
									"three_utr_lengths": [],
									"CDS_lengths": [],
									"PTC_distances": [],
									"max_downstream_PTC_distances": [],
									"three_utr_junction_counts": []},
					   "excluded": {"form_transcripts": 0,
					   				"form_coding_transcripts": 0,
					   				"five_utr_lengths": [],
									"three_utr_lengths": [],
									"CDS_lengths": [],
									"PTC_distances": [],
									"max_downstream_PTC_distances": [],
									"three_utr_junction_counts": []}}


	properties_to_be_processed = ["five_utr_lengths",
								  "three_utr_lengths",
								  "CDS_lengths",
								  "PTC_distances",
								  "max_downstream_PTC_distances",
								  "three_utr_junction_counts"
								  ]

	run_checking_processing(event_entry, 
							standard_transcript_dict, 
							form_properties,
							properties_to_be_processed)


	summarized_properties = summarize_properties(form_properties,
												 properties_to_be_processed)

	summarized_properties["included"]["form_transcripts"] = form_properties["included"]["form_transcripts"]
	summarized_properties["excluded"]["form_transcripts"] = form_properties["excluded"]["form_transcripts"]

	summarized_properties["included"]["form_coding_transcripts"] = form_properties["included"]["form_coding_transcripts"]
	summarized_properties["excluded"]["form_coding_transcripts"] = form_properties["excluded"]["form_coding_transcripts"]

	isoform_property_diffs_dict = isoform_property_diffs(summarized_properties)

	event_entry["summarized_properties"] = summarized_properties
	event_entry["isoform_property_diffs"] = isoform_property_diffs_dict





def event_statuses(standard_event_dict, 
				   standard_transcript_dict):

	for event, event_entry in standard_event_dict.iteritems():

		if (len(event_entry["included_form_transcripts"]) == 0 or 
			len(event_entry["excluded_form_transcripts"]) == 0):

			event_entry["incomplete_transcript_information"] = True

		else:

			event_entry["incomplete_transcript_information"] = False

			nmd_status(event_entry, standard_transcript_dict)

			nsd_status(event_entry, standard_transcript_dict)

			coding_noncoding_switch(event_entry, standard_transcript_dict)

			feature_overlap(event_entry, standard_transcript_dict)

			coding_seq_change(event_entry, standard_transcript_dict)

			get_form_properties(event_entry, standard_transcript_dict)



def output_table(standard_event_dict, outdir):

	output_nmd_table = open(outdir + "/event_nmd_nsd_status.tsv", "w")

	nmd_header_content = ["event_id",
					      "event_type",
					      "nmd_status", 
					      "nmd_form", 
					      "ptc_overlap", 
					      "nonstop_status", 
					      "nonstop_form"]

	output_nmd_table.write("\t".join(nmd_header_content) + "\n")


	output_table = open(outdir + "/event_coding_status.tsv", 'w')


	header_content = ["event_id",
					  "event_type",
					  "inc_form_transcripts",
					  "exc_form_transcripts",
					  "inc_form_coding_transcripts",
					  "exc_form_coding_transcripts",
					  "inc_form_mean_cds_lengths",
					  "inc_form_all_cds_lengths",
					  "inc_form_mean_utr5_lengths",
					  "inc_form_all_utr5_lengths",
					  "inc_form_mean_utr3_lengths",
					  "inc_form_all_utr3_lengths",
					  "inc_form_mean_ptc_distances",
					  "inc_form_all_ptc_distances",
					  "inc_form_mean_max_downstream_ptc_distances",
					  "inc_form_all_max_downstream_ptc_distances",
					  "inc_form_mean_utr3_junction_counts",
					  "inc_form_all_utr3_junction_counts", 
					  "exc_form_mean_cds_lengths",
					  "exc_form_all_cds_lengths",
					  "exc_form_mean_utr5_lengths",
					  "exc_form_all_utr5_lengths",
					  "exc_form_mean_utr3_lengths",
					  "exc_form_all_utr3_lengths",
					  "exc_form_mean_ptc_distances",
					  "exc_form_all_ptc_distances",
					  "exc_form_mean_max_downstream_ptc_distances",
					  "exc_form_all_max_downstream_ptc_distances",
					  "exc_form_mean_utr3_junction_counts",
					  "exc_form_all_utr3_junction_counts",
					  "mean_cds_length_inc_exc_ratio",
					  "mean_cds_length_inc_exc_diff",
					  "mean_utr5_length_inc_exc_ratio",
					  "mean_utr5_length_inc_exc_diff",
					  "mean_utr3_length_inc_exc_ratio",
					  "mean_utr3_length_inc_exc_diff",
					  "mean_ptc_dist_inc_exc_ratio",
					  "mean_ptc_dist_inc_exc_diff",	
					  "mean_max_ds_ptc_dist_inc_exc_ratio",
					  "mean_max_ds_ptc_dist_inc_exc_diff",
					  "mean_utr3_junction_counts_inc_exc_ratio",						  				  					  
					  "mean_utr3_junction_counts_inc_exc_diff",
					  "inc_feature_overlap",
					  "exc_feature_overlap",
					  "inc_cds_overlap",
					  "exc_cds_overlap",
					  "inc_utr5_overlap",
					  "exc_utr5_overlap",
					  "inc_utr3_overlap",
					  "exc_utr3_overlap",
					  "nmd_status", 
					  "nmd_form", 
					  "ptc_overlap", 
					  "nonstop_status", 
					  "nonstop_form",
					  "coding_status",
					  "coding_switch",
					  "coding_form",
					  "cds_alteration",
					  "form_unique_cds"]

	output_table.write("\t".join(header_content) + "\n")

	for event, event_entry in standard_event_dict.iteritems():

		nmd_entry_content = [event,
						 	 event_entry["event_type"],
						 	 event_entry["nmd_status"],
							 event_entry["nmd_form"],
							 event_entry["ptc_overlap"],
							 event_entry["nsd_status"],
							 event_entry["nonstop_form"]]

		output_nmd_table.write("\t".join(nmd_entry_content) + "\n")

		feature_overlap = event_entry["feature_overlap"]
		feature_overlap_boolean = event_entry["feature_overlap_boolean"]
		summarized_properties = event_entry["summarized_properties"]
		isoform_property_diffs = event_entry["isoform_property_diffs"]


		entry_content = [event,
						 event_entry["event_type"],
						 str(summarized_properties["included"]["form_transcripts"]),
						 str(summarized_properties["excluded"]["form_transcripts"]),
						 str(summarized_properties["included"]["form_coding_transcripts"]),
						 str(summarized_properties["excluded"]["form_coding_transcripts"]),
						 str(summarized_properties["included"]["mean_CDS_lengths"]),
						 str(summarized_properties["included"]["all_CDS_lengths"]),
						 str(summarized_properties["included"]["mean_five_utr_lengths"]),
						 str(summarized_properties["included"]["all_five_utr_lengths"]),
						 str(summarized_properties["included"]["mean_three_utr_lengths"]),
						 str(summarized_properties["included"]["all_three_utr_lengths"]),
						 str(summarized_properties["included"]["mean_PTC_distances"]),
						 str(summarized_properties["included"]["all_PTC_distances"]),
						 str(summarized_properties["included"]["mean_max_downstream_PTC_distances"]),
						 str(summarized_properties["included"]["all_max_downstream_PTC_distances"]),
						 str(summarized_properties["included"]["mean_three_utr_junction_counts"]),
						 str(summarized_properties["included"]["all_three_utr_junction_counts"]),
						 str(summarized_properties["excluded"]["mean_CDS_lengths"]),
						 str(summarized_properties["excluded"]["all_CDS_lengths"]),
						 str(summarized_properties["excluded"]["mean_five_utr_lengths"]),
						 str(summarized_properties["excluded"]["all_five_utr_lengths"]),
						 str(summarized_properties["excluded"]["mean_three_utr_lengths"]),
						 str(summarized_properties["excluded"]["all_three_utr_lengths"]),
						 str(summarized_properties["excluded"]["mean_PTC_distances"]),
						 str(summarized_properties["excluded"]["all_PTC_distances"]),
						 str(summarized_properties["excluded"]["mean_max_downstream_PTC_distances"]),
						 str(summarized_properties["excluded"]["all_max_downstream_PTC_distances"]),
						 str(summarized_properties["excluded"]["mean_three_utr_junction_counts"]),
						 str(summarized_properties["excluded"]["all_three_utr_junction_counts"]),
						 str(isoform_property_diffs["mean_CDS_lengths_inc_exc_ratio"]),
						 str(isoform_property_diffs["mean_CDS_lengths_inc_exc_diff"]),
						 str(isoform_property_diffs["mean_five_utr_lengths_inc_exc_ratio"]),
						 str(isoform_property_diffs["mean_five_utr_lengths_inc_exc_diff"]),
						 str(isoform_property_diffs["mean_three_utr_lengths_inc_exc_ratio"]),
						 str(isoform_property_diffs["mean_three_utr_lengths_inc_exc_diff"]),	
						 str(isoform_property_diffs["mean_PTC_distances_inc_exc_ratio"]),
						 str(isoform_property_diffs["mean_PTC_distances_inc_exc_diff"]),	
						 str(isoform_property_diffs["mean_max_downstream_PTC_distances_inc_exc_ratio"]),				 						 						 
						 str(isoform_property_diffs["mean_max_downstream_PTC_distances_inc_exc_diff"]),
						 str(isoform_property_diffs["mean_three_utr_junction_counts_inc_exc_ratio"]),
						 str(isoform_property_diffs["mean_three_utr_junction_counts_inc_exc_diff"]),						 
					  	 ",".join(feature_overlap["included"]),
					  	 ",".join(feature_overlap["excluded"]),
					  	 str(feature_overlap_boolean["included"]["cds"]),
					  	 str(feature_overlap_boolean["excluded"]["cds"]),
					  	 str(feature_overlap_boolean["included"]["utr5"]),
					  	 str(feature_overlap_boolean["excluded"]["utr5"]),
					  	 str(feature_overlap_boolean["included"]["utr3"]),
					  	 str(feature_overlap_boolean["excluded"]["utr3"]),						 
						 event_entry["nmd_status"],
						 event_entry["nmd_form"],
						 event_entry["ptc_overlap"],
						 event_entry["nsd_status"],
						 event_entry["nonstop_form"],
						 str(event_entry["coding_status"]),
						 str(event_entry["coding_switch"]),
						 str(event_entry["coding_form"]),
						 str(event_entry["cds_alteration"]),
						 str(event_entry["form_unique_cds"])
						 ]

		output_table.write("\t".join(entry_content) + "\n")

	output_nmd_table.close()
	output_table.close()






def main(args, standard_transcript_dict = None, standard_event_dict = None):

	parser = argparse.ArgumentParser()

	parser.add_argument("--suppress_output", 
						action = "store_true", 
						help = "If set, not file will be output")

	parser.add_argument("--outdir", 
						type = str, 
						help = "Output directory")

	parser.add_argument("--ioe_files", 
						type = str, 
						help = "Event ioe file path")

	parser.add_argument("--event_gtf", 
						type = str, 
						help = "Event gtf file - required for certain " + 
							   "features to work e.g. PTC overlap " + 
							   "determination if standard_event_dict " + 
							   "not supplied as argument to main()")

	parser.add_argument("--transcript_dict_pkl", 
						type = str, 
						help = "Pickled transcript dict from " + 
							   "cds_insertion.py.  Required for " + 
							   "certain features to work e.g. PTC " + 
							   "overlap determination if " + 
							   "standard_transcript_dict not " + 
							   "supplied as argument to main(). ")

	args = parser.parse_args(args)


	suppress_output = args.suppress_output
	outdir = args.outdir
	ioe_files = args.ioe_files
	#transcript_table = args.transcript_table
	event_gtf = args.event_gtf
	transcript_dict_pkls = args.transcript_dict_pkls 

	if standard_event_dict is None and event_gtf is not None:

		print "Importing event data"

		standard_event_dict = splice_lib.generate_standard_event_dict(event_gtf)

		print "Event data imported.  Now adding transcripts to event dict via IOE file"

		splice_lib.add_transcripts_to_event_dict(ioe_files, standard_event_dict)

		print "Transcripts added to event dict.  Now finding exons unique to either form"

		splice_lib.find_exons_unique_to_form(standard_event_dict)

		print "Form-specific exons identified.  Now extracting alternative regions"

		for event, event_entry in standard_event_dict.iteritems():

			(
			 event_entry["included_alt_regions"], 
			 event_entry["excluded_alt_regions"]
			 ) = splice_lib.get_alt_regions(event_entry["included_exons"], 
										    event_entry["excluded_exons"], 
										    event_entry["included_unique_exons"], 
										    event_entry["excluded_unique_exons"])

	else:

		sys.exit("Please supply either standard_event_dict as " +
			     "argument to find_switch_events main, or supply " +
			     "path to --event_gtf.  Supplying both is not an " +
			     "option (it is confusing!)")


	if standard_transcript_dict is None and transcript_dict_pkls is not None:

		print "Importing transcript data"

		import cPickle as pkl

		transcript_dict_pkl_paths = transcript_dict_pkls.split(",")

		standard_transcript_dict = dict(chain(*map(dict.items, [pkl.load(open(i, "rb")) for i in
									transcript_dict_pkl_paths])))

	else:

		sys.exit("Please supply either standard_transcript_dict " +
			     "as arg to main() or supply path to --transcript_dict_pkls (not both)")


	if not suppress_output and outdir is None:

		sys.exit("Output directory required if suppress output not " +
			     "set in find_switch_events.py")

	print "Establishing event statuses"

	event_statuses(standard_event_dict, standard_transcript_dict)
	output_table(standard_event_dict, outdir)
	

if __name__ == '__main__':

	main(sys.argv[1:])




