
bam_dir=/public/groups/sanfordlab/people/anjowall/projects/nmd_validation/bam/
ref_path=/public/groups/sanfordlab/people/anjowall/projects/nmd_validation/gencode.v29.basic.annotation.gtf
stringtie_out=/public/groups/sanfordlab/people/anjowall/projects/nmd_validation/stringtie/

mkdir -p $stringtie_out

for i in $bam_dir/*.bam
do
stringtie "$i" -G $ref_path -o $stringtie_out/$(echo "$i" | awk -F\/ '{ print $NF}').gtf & done

wait; stringtie --merge -i -G $ref_path -o $stringtie_out/all_merged.gtf $stringtie_out/*.gtf

