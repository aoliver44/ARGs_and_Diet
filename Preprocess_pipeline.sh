#!/bin/bash

# After checksum and unzipping downloaded files
raw_dir=
run_dir=

## Uncompress file
unzip_dir=$run_dir/unzipped
mkdir -p $unzip_dir

for file in $raw_dir/*.fastq.gz
do
    #get the base
    basename "$file"
    f="$(basename -- $file)"
    ##echo "$f"

    #without the gz
    file1="${f%.gz}"
    ##echo "$file1"

    if [ -f "$unzip_dir/$file1" ]
    then
	echo "$unzip_dir/$file1 already exist and will not be overwritten."
    else
	echo "$unzip_dir/$file1 does not exist. Decompressing now..."
	gunzip -c $file > $unzip_dir/$file1
    fi
done

echo "STEP 0 DONE AT: "; date


###################################################################
#
# STEP 1: REMOVING HUMAN READS USING BMTAGGER
# Note: paired-end files are usually named using R1 and R2 in the name.
# Note: if using single-end reads, only need to specify one input flag (-1)

## File paths
unzip_outdir=$run_dir/unzipped
mkdir -p $run_dir/step1_BMTagger/
bmt_outdir=$run_dir/step1_BMTagger
#
## Software paths
bmtagger_location=programs/bmtools/bmtagger
human_db=database/human_GRCh38_p13

PATH=$PATH:/share/lemaylab-backedup/milklab/programs/bmtools/bmtagger
PATH=$PATH:/share/lemaylab-backedup/milklab/programs/srprism/gnuac/app
module load blast
module load java bbmap
echo "NOW STARTING HUMAN READ REMOVAL STEP AT: "; date

for file in $unzip_outdir/*R1_001.fastq
do
  file1=$file
  file2=$(echo $file1 | sed 's/R1_001/R2_001/')
	filename=$(basename "$file1")
  basename=$(echo $filename | cut -f 1 -d "_")
	outname="$bmt_outdir/$basename"

  if [ -f $outname.human.txt ]
  then
    echo $outname.human.txt already exist and will not be overwritten.
  else
    echo $outname.human.txt does not exist. Running BMTagger now...
    $bmtagger_location/bmtagger.sh -b $human_db/GCF_000001405.39_GRCh38.p13_genomic.bitmask -x $human_db/GCF_000001405.39_GRCh38.p13_genomic.srprism -q 1 -1 $file1 -2 $file2 -o $outname.human.txt

    # filterbyname.sh is included in the bbmap module
    # bbmap module location: /software/bbmap/37.68/static/bbmap.sh
    # I included a local copy of this program: /share/lemaylab-backedup/milklab/programs/filterbyname_v37.68.sh
    # This script removes sequences in both R1 and R2 that matches the human reads
    # (sequence header is passed to the script in the outname.human.txt file)
    filterbyname.sh in=$file1 in2=$file2 out=$outname.R1_nohuman.fastq out2=$outname.R2_nohuman.fastq names=$outname.human.txt include=f
  fi
done

echo "STEP 1 DONE AT: "; date

####################################################################


###################################################################
#
# STEP 2: USE TRIMMOMATIC TO REMOVE LOW-QUALITY READS
# Take input from BMTagger (removal of the human DNA)
# Trimmomatic parameters followed this paper (Taft et al., 2018, mSphere) and using the paired end mode
# from link: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

echo "NOW STARTING READ CLEANING WITH TRIMMOMATIC AT: "; date

# File paths
mkdir -p $run_dir/step2_trim # only need to run once
trim_outdir=$run_dir/step2_trim

# Software paths
# Trimmomatic
trimmomatic_location=/software/trimmomatic/0.33/static/trimmomatic-0.33.jar

for file in $bmt_outdir/*.R1_nohuman.fastq
do
  STEM=$(basename "${file}" .R1_nohuman.fastq)

  file1=$file
  file2=$bmt_outdir/${STEM}.R2_nohuman.fastq

  # remove adapter and trimming at the same time (TrueSeq3-PE-2.fa, PE1_rc and PE2_rc)
  java -jar $trimmomatic_location PE -threads 15 $file1 $file2 $trim_outdir/${STEM}_R1_paired.fastq $trim_outdir/${STEM}_R1_unpaired.fastq.gz $trim_outdir/${STEM}_R2_paired.fastq $trim_outdir/${STEM}_R2_unpaired.fastq.gz ILLUMINACLIP:/share/lemaylab-backedup/milklab/programs/trimmomatic_adapter_input.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:99
done

echo "STEP 2 DONE AT: "; date

####################################################################
#
# STEP 3: Remove duplicated reads
# For FastUniq options: https://wiki.gacrc.uga.edu/wiki/FastUniq

echo "NOW STARTING REMOVING DUPLICATE READS AT: "; date

# File paths
mkdir -p $run_dir/step3_fastuniq
dup_outdir=$run_dir/step3_fastuniq

# Software paths
# FastUniq
fastuniq_location=/software/fastuniq/1.1/lssc0-linux/source/fastuniq

touch $dup_outdir/fastuniq_input_list.txt
fastuniq_input_list=$dup_outdir/fastuniq_input_list.txt

for file in $trim_outdir/*_R1_paired.fastq
do
  # clear the content of the list file so that FastUniq takes only 2 files every time it runs
  > $fastuniq_input_list

  STEM=$(basename "$file" _R1_paired.fastq)

  file1=$file
  file2=$trim_outdir/${STEM}_R2_paired.fastq

  echo $file1 >> $fastuniq_input_list
  echo $file2 >> $fastuniq_input_list

  # Althoug FastUniq takes list as input !!!!DO NOT SUPPLY ALL FILES IN ONE LIST!!!!
  # Because FastUniq can only write out 2 files for forward and reverse
  # Run FastUniq
  $fastuniq_location -i $fastuniq_input_list -t q -o $dup_outdir/${STEM}_R1_dup.fastq -p $dup_outdir/${STEM}_R2_dup.fastq
done

echo "STEP 3 DONE AT: "; date

###################################################################
#
# STEP 4: Merge paired-end reads with FLASH
# The FLASH manual link: http://ccb.jhu.edu/software/FLASH/MANUAL

echo "NOW STARTING PAIRED-END MERGING WITH FLASH AT: "; date

# File paths
mkdir -p $run_dir/step4_flash
flash_outdir=$run_dir/step4_flash

# Software paths
# FLASH
flash_location=programs/FLASH-1.2.11_2019/FLASH-1.2.11-Linux-x86_64/flash

for file in $dup_outdir/*_R1_dup.fastq
do
  STEM=$(basename "${file}" _R1_dup.fastq)

  file1=$file
  file2=$dup_outdir/${STEM}_R2_dup.fastq

  # -m: minium overlap length 10bp to be similar to pear
  # -M: max overlap length
  # -x: mismatch ratio, default is 0.25, which is quite high (e.g: 50bp overlap --> 12.5 mismatch by default)
  $flash_location $file1 $file2 -m 10 -M 100 -x 0.1 -o ${STEM} -d $flash_outdir
done

echo "STEP 4 DONE AT: "; date

####################################################################
#
# STEP 5: Normalize for sample average genome size and RPKG
# The MicrobeCensus help page: https://github.com/snayfach/MicrobeCensus

echo "NOW STARTING NORMALIZATION WITH MicrobeCensus AT: "; date

# File paths
mkdir -p $run_dir/step5_MicrobeCensus
mc_outdir=$run_dir/step5_MicrobeCensus

# Software paths
# MicrobeCensus location
microbecensus=/programs/MicrobeCensus-1.1.1/scripts/run_microbe_census_nomodule.py
# External RAPsearch2 v2.15 binary. I would like to have the binary NOT in the xzyao/.local home directory
# because the home directory gets wiped clean every time I log out
RAPSEARCH=/programs/MicrobeCensus-1.1.1/microbe_census/bin/rapsearch_Linux_2.15

for file in $dup_outdir/*_R1_dup.fastq
do
  STEM=$(basename "$file" _R1_dup.fastq)

  if [ -f $mc_outdir/${STEM}_mc.txt ]
  then
    echo "$file exist"
  else
    echo "Processing sample $file now"

    # Change the path to your home path on spitfire.
    # Actually no longer needed after specifying the RAPsearch2 location
    ##export PATH="/home/xzyao/miniconda3/bin:$PATH"
    ##export PATH="/home/xzyao/.local/bin:$PATH"

    # change dir for writing temporary files
    export TMPDIR=$mc_outdir

    file2=$dup_outdir/${STEM}_R2_dup.fastq

    # -h for help
    # -l read length to cut at, should be 150 for Novaseq paired end samples
    # -t thread number for rapsearch, microbecensus only uses 1 thread
    # -n number of reads to sample from seqfile and use for AGS estimation
    ## set at 100 million reads to use all reads (100 million reads should be more than the biggest lib size)
    $microbecensus $file,$file2 $mc_outdir/${STEM}_mc.txt \
    -r $RAPSEARCH \
    -l 150 -t 20 -n 100000000 #change per run depending library size (read numbers per sample)
  fi
done

echo "STEP 5 DONE AT: "; date

# STEP 5 for merged reads
# The caculated genome equivalents will be used to normalized KEGG, CAZy, beta-galactosidases
# and other count tables generated with merged reads (i.e. database is amino acid sequence)
# see /share/lemaylab-backedup/Zeya/scripts/ARG_metagenome/non_pipeline
# for scripts to get KEGG, CAZy and b-gal count tables

echo "NOW STARTING NORMALIZATION WITH MicrobeCensus FOR MERGED READS AT: "; date

## Set input and output file paths
mkdir -p $run_dir/step5_MicrobeCensus_merged
mc_outdir_merged=$run_dir/step5_MicrobeCensus_merged

for file in $flash_outdir/*.extendedFrags.fastq
do
  STEM=$(basename "$file" .extendedFrags.fastq)

  if [ -f $mc_outdir_merged/${STEM}_mc.txt ]
  then
    echo "$file exist"
  else
    echo "Processing sample $file now (merged reads)"

    # change dir for writing temporary files
    export TMPDIR=$mc_outdir_merged

    # -h for help
    # -l read length to cut at, should be 150 for Novaseq paired end samples
    # -t thread number for rapsearch, microbecensus only uses 1 thread
    # -n number of reads to sample from seqfile and use for AGS estimation
    $microbecensus $file $mc_outdir_merged/${STEM}_mc.txt \
    -r $RAPSEARCH \
    -l 150 -t 20 \
    -n 100000000 # Double check this n number to make sure it is more than the biggest library size for merged reads
  fi
done

echo "STEP 5 ON MERGED READS DONE AT: "; date

######################################################################################################################
#
# STEP 6: Align to MEGARes database with bwa

echo "NOW STARTING ALIGNMENT TO MEGARES_v2 WITH BWA AT: "; date

# File paths
mkdir -p $run_dir/step6_megares_bwa
megares_outdir=$run_dir/step6_megares_bwa

# Software paths
# megares db including both drug, metal and biocide resistance genes
megares_dir=/database/megares_v2/

# Make sure bash knows where to look for softwares
# Location and version of bwa module: /software/bwa/0.7.16a/lssc0-linux/bwa
module load bwa

for file in $dup_outdir/*_R1_dup.fastq
do
  STEM=$(basename "$file" _R1_dup.fastq)

  file1=$file
  file2=$dup_outdir/${STEM}_R2_dup.fastq

  # use the bwa mem method for fastest speed and accuracy
  # BWA aligner manual: http://bio-bwa.sourceforge.net/bwa.shtml
  # -t: thread
  bwa mem -t 20 $megares_dir/megares_modified_database_v2.00.fasta $file $file2 > $megares_outdir/${STEM}_align.sam
done

echo "STEP 6 DONE AT: "; date

#######################################################################################################################
#
# STEP 7. Count the MEGARes database alignment and normalize the count
# table with MicrobeCensus generated genome equivalents per sample

echo "NOW MAKING COUNT TABLE AT: "; date

# File paths
mkdir -p $run_dir/step7_norm_count_tab/
mkdir -p $run_dir/step7_norm_count_tab/resistomeanalyzer_output
ra_outdir=$run_dir/step7_norm_count_tab/resistomeanalyzer_output

# Software paths
ranalyzer=/programs/resistomeanalyzer/resistome

for file in $megares_outdir/*_align.sam
do
  STEM=$(basename "$file" _align.sam)

  # usage instructions: https://github.com/cdeanj/resistomeanalyzer
  $ranalyzer \
  -ref_fp $megares_dir/megares_modified_database_v2.00.fasta \
  -sam_fp $file \
  -annot_fp $megares_dir/megares_modified_annotations_v2.00.csv \
  -gene_fp $ra_outdir/${STEM}_gene.tsv \
  -group_fp $ra_outdir/${STEM}_group.tsv \
  -class_fp $ra_outdir/${STEM}_class.tsv \
  -mech_fp $ra_outdir/${STEM}_mechanism.tsv \
  -t 80 #Threshold to determine gene significance
done

echo "NOW STARTING COUNT TABLE NORMALIZATION AT: "; date

# File paths
mkdir -p $run_dir/step7_norm_count_tab/normalized_tab
norm_outdir=$run_dir/step7_norm_count_tab/normalized_tab

# Software paths
norm=/ARG_metagenome/make_RPKG_normtab.py
# Make sure Numpy and pandas are downloaded in the same directory as $norm
## See /share/lemaylab-backedup/Zeya/scripts/ARG_metagenome/init/install_python_packages

for file in $ra_outdir/*_gene.tsv
do
  STEM=$(basename "$file" _gene.tsv)

  python $norm \
  --mc $mc_outdir/${STEM}_mc.txt \
  --genelen /database/megares_v2/megares_modified_database_v2_GeneLen_org.tsv \
  --count $file \
  --out $norm_outdir/${STEM}_norm.csv \
  --mergein $file # This argument does not do anything here, but python2 requires *args to be not empty so I supply a dummy argument (if using python3, can omit this input)
done

# use the python script again to merge tables
python $norm --mergeout $norm_outdir/merge_norm_final.csv --mergein $norm_outdir/*_norm.csv #

echo "STEP 7 DONE AT: "; date



#####################################################################
##
## STEP 8. Assemble reads into contigs with megaHIT

echo "NOW STARTING ASEEMBLY WITH MEGAHIT AT: "; date

# File paths
mkdir -p $run_dir/step8_megahit_in_trimmomatic
megahit_outdir=$run_dir/step8_megahit_in_trimmomatic

# Software paths
megahit=/programs/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit

# assemble with paired end reads input (after trimmomatic)
for file in $trim_outdir/*_R1_paired.fastq
do
        STEM=$(basename "$file" _R1_paired.fastq)

        if [ -f $megahit_outdir/${STEM}_assembled/final.contigs.fa ]
        then
                echo "$megahit_outdir/${STEM}_assembled/final.contigs.fa exist"
        else
                # usage instruction at https://github.com/voutcn/megahit & http://www.metagenomics.wiki/tools/assembly/megahit
                # use 70% memory and 20 threads
                $megahit -1 $file -2 $trim_outdir/${STEM}_R2_paired.fastq \
                -m 0.7 -t 20 \
                -o $megahit_outdir/${STEM}_assembled
        fi
done

echo "STEP 8 DONE AT: "; date

########################################################################################################################
##
## STEP 9. ALIGN SHORT READS CONTAINING ARG TO CONTIGS

echo "NOW STARTING SHORT READ AND CONTIG ALIGNMENT AT: "; date

# File paths
mkdir -p $run_dir/step9_contig_bwa_nomerg
aln_outdir=$run_dir/step9_contig_bwa_nomerg
mkdir -p $aln_outdir/mapped_fastq

# Software paths
bbmap=/software/bbmap/37.68/static/
bwa=/software/bwa/0.7.16a/lssc0-linux/bwa
pyhd=/scripts/ARG_metagenome/prefix_to_compline.py

for file in $megares_outdir/*_align.sam
do
  STEM=$(basename "$file" _align.sam)

  if [ -f $aln_outdir/${STEM}_contig_aln.fasta ]
  then
    echo "$aln_outdir/${STEM}_contig_aln.fasta exist"
  else
    echo "Processing sample ${STEM} now......"

    # convert sam file to fastq and keep only reads that are mapped
    # usage: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/reformat-guide/
    # https://github.com/BioInfoTools/BBMap/blob/master/sh/reformat.sh
    $bbmap/reformat.sh in=$file out=$aln_outdir/mapped_fastq/${STEM}_megares_map.fastq mappedonly

    # use bwa to index and align (aln) contigs with short reads containing ARGs
    $bwa index $megahit_outdir/${STEM}_assembled/final.contigs.fa
    $bwa mem -t 15 $megahit_outdir/${STEM}_assembled/final.contigs.fa $aln_outdir/mapped_fastq/${STEM}_megares_map.fastq > $aln_outdir/${STEM}_contig_ARGread_aln.sam

    # get the contigs header containing ARG
    $bbmap/reformat.sh in=$aln_outdir/${STEM}_contig_ARGread_aln.sam out=$aln_outdir/${STEM}_contig_ARGread_aln_mappedonly.sam mappedonly

    # gather A0-based run header from the sam file & add a space after each line for the 1st field of the actual contig header (annoying reformating due to spaces in the contig header, e.g: ">k141_49608 flag=0 multi=1.0000 len=233")
    grep 'A0' $aln_outdir/${STEM}_contig_ARGread_aln_mappedonly.sam | cut -f 3 > $aln_outdir/${STEM}_header.txt
    python $pyhd --i $aln_outdir/${STEM}_header.txt --f $megahit_outdir/${STEM}_assembled/final.contigs.fa --o $aln_outdir/${STEM}_header.txt

    # filter use the above list to retain ARG read aligned contigs
    # http://seqanswers.com/forums/archive/index.php/t-75650.html
    $bbmap/filterbyname.sh in=$megahit_outdir/${STEM}_assembled/final.contigs.fa out=$aln_outdir/${STEM}_contig_ARGread_aln_mappedonly.fa names=$aln_outdir/${STEM}_header.txt include=t
  fi
done

echo "STEP 9 DONE AT: "; date

########################################################################################################################
##
## STEP 10. ID the taxonomy of contigs using CAT

echo "NOW STARTING TAXONOMY ID FOR ARG-containing contigs AT: "; date

# File paths
mkdir -p $run_dir/step10_CAT
CAT_outdir=$run_dir/step10_CAT
cd $CAT_outdir

# Software paths
CAT=/programs/CAT-5.0.3/CAT_pack/CAT # CAT version and database to be updated soon!
CATdb=/share/lemaylab-backedup/milklab/database/CAT_prepare_20190719
progdigal=/software/prodigal/2.6.3/x86_64-linux-ubuntu14.04/bin/prodigal
diamond=/programs/diamond


for file in $aln_outdir/*_contig_ARGread_aln_mappedonly.fa
do
        STEM=$(basename "$file" _contig_ARGread_aln_mappedonly.fa)

        echo "Processing sample $STEM now...... "
        # for help /share/lemaylab-backedup/milklab/programs/CAT-5.0.3/CAT_pack/CAT contigs -h
        # https://github.com/dutilh/CAT

        $CAT contigs -c $file -d $CATdb/2019-07-19_CAT_database -t $CATdb/2019-07-19_taxonomy --path_to_prodigal $progdigal --path_to_diamond $diamond -o $CAT_outdir/${STEM}_CAT -n 25

        $CAT add_names -i $CAT_outdir/${STEM}_CAT.contig2classification.txt -o $CAT_outdir/${STEM}.taxaid.txt -t $CATdb/2019-07-19_taxonomy --only_official

        $CAT summarise -c $file -i $CAT_outdir/${STEM}.taxaid.txt -o $CAT_outdir/${STEM}.taxaname.txt
done

echo "STEP 10 DONE AT: "; date