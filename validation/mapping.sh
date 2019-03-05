source ~/.bash_profile

num_threads=20

export star_exe=/public/home/anjowall/STAR/bin/Linux_x86_64/STAR
export samtools_exe=/public/home/anjowall/samtools-1.8/bin/samtools

###These dirs must exist
export topdir=/public/groups/sanfordlab/people/anjowall/projects/nmd_validation/
export fastq_dir="$topdir"/fastq/
export genome_dir=/public/groups/sanfordlab/people/anjowall/genomes/GRCh38/STAR_genome_GRCh38/
export genome_fasta=/public/groups/sanfordlab/people/anjowall/genomes/GRCh38/STAR_genome_GRCh38/GRCh38.primary_assembly.genome.fa

###These dirs/files can be created
export bam_dir_pass_one="$topdir"/bam/

export fastq_prefix_list=$(ls "$fastq_dir"/*fastq.gz | awk -F\/ '{ print $NF}' | awk -F"_[12].fastq" '{ print $1 }' | sort | uniq)


mkdir -p "$bam_dir_pass_one"


#####STAR indexing
##Human

if [ ! -e "$genome_dir"/SA ]
        then
                echo "genome suffix array file not found . . . running STAR genome generate . . . "
                $star_exe --runMode genomeGenerate --genomeDir "$genome_dir" --genomeFastaFiles "$genome_fasta" --runThreadN "$num_threads" --genomeSAindexNbases 14 --limitGenomeGenerateRAM=140000000000 
else
        echo "Genome suffix array file found - assuming genome is indexed - attempting to proceed . . . "
fi


$star_exe --genomeDir $genome_dir --genomeLoad Remove
$star_exe --genomeDir $genome_dir --genomeLoad LoadAndExit


#####STAR mapping
echo "Attempting to map with STAR"

for prefix in ${fastq_prefix_list[*]}
        do
                if [ ! -e "$bam_dir_pass_one"/"$prefix""_pass_1_""Aligned.sortedByCoord.out.bam" ]
                        then
                                echo "Mapping sample  . . . "
                                $star_exe --runMode alignReads --readFilesIn "$fastq_dir""$prefix"_1.fastq.gz "$fastq_dir""$prefix"_2.fastq.gz --readFilesCommand zcat --genomeDir "$genome_dir" --outSAMtype BAM SortedByCoordinate --runThreadN "$num_threads" --alignEndsType Local --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20 --chimScoreMin 1 --outSAMattributes NH HI AS nM --outSAMattrIHstart 0 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix "$bam_dir_pass_one"/"$prefix""_pass_1_" --outFilterMultimapNmax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM=60688817693
                        echo "Sample already mapped . . .  proceeding . . . "
                fi
        done


$star_exe --genomeDir $genome_dir --genomeLoad Remove
