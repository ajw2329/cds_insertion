topdir=/public/groups/sanfordlab/people/anjowall/projects/nmd_validation/
script_dir=$topdir/scripts/
util_dir=/public/groups/sanfordlab/people/anjowall/bin/

bash $script_dir/mapping.sh

bash $script_dir/pcr_dup_rm.sh 12

bash $script_dir/stringtie.sh

bash $script_dir/junctioncounts_dup_rm.sh

gffread $topdir/stringtie/all_merged.gtf -g /public/groups/sanfordlab/people/anjowall/genomes/GRCh38/STAR_genome_GRCh38/GRCh38.primary_assembly.genome.fa -w $topdir/stringtie/all_merged.fa

/public/home/anjowall/anaconda2/bin/python2.7 /public/home/anjowall/repos/cds_insertion/cds_insertion.py --transcript_gtf $topdir/stringtie/all_merged.gtf --transcript_fasta $topdir/stringtie/all_merged.fa --annotation_gtf $topdir/gencode.v29.basic.annotation.gtf --outdir $topdir/cds_insertion_out/ --CCDS --bigGenePred_as_path $util_dir/bigGenePred.as --gtfToGenePred_path $util_dir/gtfToGenePred --genePredToBigGenePred_path $util_dir/genePredToBigGenePred --bedToBigBed_path $util_dir/bedToBigBed --bedtools_path /public/home/anjowall/bedtools2/bin/bedtools --chrNameLength_path /public/groups/sanfordlab/people/anjowall/genomes/GRCh38/STAR_genome_GRCh38/chrNameLength.txt --make_bigBed

/public/home/anjowall/anaconda2/bin/python2.7 /public/home/anjowall/repos/cds_insertion/find_switch_events.py --outdir $topdir/cds_insertion_out/ --ioe_file $topdir/stringtie/splice_lib_events.ioe --event_gtf $topdir/stringtie/splice_lib_events.gtf --transcript_dict_pkl $topdir/cds_insertion_out/cds_insertion_transcript_dict.pkl

