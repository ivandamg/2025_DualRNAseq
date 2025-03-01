# 2025_DualRNAseq


1. Download raw files with wget

         sbatch --partition=pshort_el8 --job-name=wget --time=0-02:00:00 --mem-per-cpu=12G --ntasks=1 --cpus-per-task=1 --output=wget.out --error=wget.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/10_dualRNAseq/01_RawData; wget https://lims.bioinformatics.unibe.ch/symlink/WT1_L1_R1_001_gGg9DyMWQgI6.fastq.gz; wget FILES"


2. Create  new folder and merge L1 and L2

         cat Argon_2_w1_L1_R1_001_kxBGHkWipUiF.fastq.gz > ../02_MergedData/ArgonW2_1_R1.fastq.gz
         cat Argon_2_w1_L1_R2_001_amSx0WBsRjBm.fastq.gz > ../02_MergedData/ArgonW2_1_R2.fastq.gz
         cat Argon_2_w1_L2_R1_001_YUpNROxlGVL0.fastq.gz >> ../02_MergedData/ArgonW2_1_R1.fastq.gz
         cat Argon_2_w1_L2_R2_001_sGU2ZcsylpCd.fastq.gz >> ../02_MergedData/ArgonW2_1_R2.fastq.gz


3. Clean reads with fastp

   
         cd /data/users/imateusgonzalez/2024_RNAseqTrial/01_RawData

         mkdir ../03_TrimmedData/


         for FILE in $(ls *R1.fastq.gz); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2)fastp --time=0-08:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=4 --output=$(echo $FILE | cut -d'_' -f1,2)_fastp.out --error=$(echo $FILE | cut -d'_' -f1,2)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/10_dualRNAseq/02_MergedData ; module load FastQC; module load fastp; fastp --in1 $FILE --in2 $(echo $FILE | cut -d'_' -f1,2)_R2.fastq.gz --out1 ../03_TrimmedData/$(echo $FILE | cut -d'_' -f1,2)_1_trimmed.fastq.gz --out2 ../03_TrimmedData/$(echo $FILE | cut -d'_' -f1,2)_2_trimmed.fastq.gz -h ../03_TrimmedData/$(echo $FILE | cut -d',' -f1,2)_fastp.html --thread 4; fastqc -t 4 ../03_TrimmedData/$(echo $FILE | cut -d'_' -f1,2)_1_trimmed.fastq.gz; fastqc -t 4 ../03_TrimmedData/$(echo $FILE | cut -d'_' -f1,2)_2_trimmed.fastq.gz"; sleep  1; done

2. Index reference Medicago

       sbatch --partition=pshort_el8 --job-name=StarIndex --time=0-01:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=StarIndex.out --error=StarIndex.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/10_dualRNAseq/00_References; module load STAR; STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /data/projects/p495_SinorhizobiumMeliloti/10_dualRNAseq/00_References --genomeFastaFiles GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna --sjdbGTFfile GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gff --sjdbOverhang 99 --genomeSAindexNbases 10"


3. Map reads to medicago

          for FILE in $(ls Argon_*_1_trimmed.fastq.gz ); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2)_1STAR --time=0-08:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_STAR.out --error=$(echo $FILE | cut -d'_' -f1,2)_STAR.error --mail-type=END,FAIL --wrap "module load STAR; cd /data/projects/p495_SinorhizobiumMeliloti/10_dualRNAseq/03_TrimmedData; STAR --runThreadN 8 --genomeDir /data/projects/p495_SinorhizobiumMeliloti/10_dualRNAseq/00_References --readFilesIn $FILE $(echo $FILE | cut -d'_' -f1,2)_2_trimmed.fastq.gz --readFilesCommand zcat --outFileNamePrefix $(echo $FILE | cut -d'_' -f1,2)_MappedMedicago --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --limitBAMsortRAM 1065539232; mv $(echo $FILE | cut -d'_' -f1,2)_MappedMedicagoUnmapped.out.mate1 ../04_MappedMedicago/$(echo $FILE | cut -d'_' -f1,2)_UnmappedMedicago_R1.fq; mv $(echo $FILE | cut -d'_' -f1,2)_MappedMedicagoUnmapped.out.mate2 ../04_MappedMedicago/$(echo $FILE | cut -d'_' -f1,2)_UnmappedMedicago_R2.fq; gzip ../04_MappedMedicago/$(echo $FILE | cut -d'_' -f1,2)_UnmappedMedicago_R1.fq; gzip ../04_MappedMedicago/$(echo $FILE | cut -d'_' -f1,2)_UnmappedMedicago_R2.fq"; sleep  1; done


3b. Map reads to Rhizobia directly. For Bacterial culture samples

          for FILE in $(ls trialNat*_1_trimmed.fastq.gz ); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2)_1STAR --time=0-04:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_STAR.out --error=$(echo $FILE | cut -d'_' -f1,2)_STAR.error --mail-type=END,FAIL --wrap "module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; cd /data/users/imateusgonzalez/2024_RNAseqTrial/02_TrimmedData; STAR --runThreadN 8 --limitOutSJcollapsed 2000000 --genomeDir /data/users/imateusgonzalez/2024_RNAseqTrial/00_ReferenceGenomes/01_Rhizobia --readFilesIn $FILE $(echo $FILE | cut -d'_' -f1,2)_2_trimmed.fastq.gz --readFilesCommand zcat --outFileNamePrefix $(echo $FILE | cut -d'_' -f1,2)_MappedRhizobia --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --limitBAMsortRAM 1065539232; mv $(echo $FILE | cut -d'_' -f1,2)_MappedRhizobiaUnmapped.out.mate1 ../03b_MapRhizobia_trialNat/$(echo $FILE | cut -d'_' -f1,2)_UnmappedRhizobia_R1.fastq; mv $(echo $FILE | cut -d'_' -f1,2)_MappedRhizobiaUnmapped.out.mate2 ../03b_MapRhizobia_trialNat/$(echo $FILE | cut -d'_' -f1,2)_UnmappedRhizobia_R2.fastq "; sleep  1; done


4. Index genome rhizobia

       sbatch --partition=pshort_el8 --job-name=StarIndex --time=0-01:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=StarIndex.out --error=StarIndex.error --mail-type=END,FAIL --wrap "cd /data/users/imateusgonzalez/2024_RNAseqTrial/00_ReferenceGenomes/01_Rhizobia; module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /data/users/imateusgonzalez/2024_RNAseqTrial/00_ReferenceGenomes/01_Rhizobia --genomeFastaFiles FribourgSMeliloti_Prokka.fna --sjdbGTFfile FribourgSMeliloti_Prokka_v2.gff --sjdbGTFfeatureExon CDS --sjdbOverhang 99 --genomeSAindexNbases 10"


5a. Noudles data: Map reads to Rhizobia

          for FILE in $(ls 150mg*_UnmappedMedicago_R1.fastq ); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2)_1STAR --time=0-04:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_STAR.out --error=$(echo $FILE | cut -d'_' -f1,2)_STAR.error --mail-type=END,FAIL --wrap "module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; cd /data/users/imateusgonzalez/2024_RNAseqTrial/03_UnmappedReads; STAR --runThreadN 8 --genomeDir /data/users/imateusgonzalez/2024_RNAseqTrial/00_ReferenceGenomes/01_Rhizobia --readFilesIn $FILE $(echo $FILE | cut -d'_' -f1,2)_UnmappedMedicago_R2.fastq --outFileNamePrefix $(echo $FILE | cut -d'_' -f1,2)_MappedRhizobia --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 1065539232"; sleep  1; done

5b. Bacterial liquid culture data: Map reads to Rhizobia. Changed  max number of collapsed junctions 1000000 -> 2000000
 
          for FILE in $(ls trial*_UnmappedMedicago_R1.fastq ); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2)_1STAR --time=0-04:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_STAR.out --error=$(echo $FILE | cut -d'_' -f1,2)_STAR.error --mail-type=END,FAIL --wrap "module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; cd /data/users/imateusgonzalez/2024_RNAseqTrial/03_UnmappedReads; STAR --runThreadN 8 --limitOutSJcollapsed 2000000 --genomeDir /data/users/imateusgonzalez/2024_RNAseqTrial/00_ReferenceGenomes/01_Rhizobia --readFilesIn $FILE $(echo $FILE | cut -d'_' -f1,2)_UnmappedMedicago_R2.fastq --outFileNamePrefix $(echo $FILE | cut -d'_' -f1,2)_MappedRhizobia --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 1065539232"; sleep  1; done



6. Count reads on features # need to install tool
https://subread.sourceforge.net/featureCounts.html



USE gff from Genome mapped reference
    
       for FILE in $(ls *.out.bam ); do echo $FILE; sbatch --partition=pshort_el8 --job-name=FC_$(echo $FILE | cut -d'_' -f1,2) --time=0-01:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=FC_$(echo $FILE | cut -d'_' -f1,2).out --error=FC_$(echo $FILE | cut -d'_' -f1,2).error --mail-type=END,FAIL --wrap "module load Subread; featureCounts -p --countReadPairs -t CDS -g ID -a /data/users/imateusgonzalez/2024_RNAseqTrial/00_ReferenceGenomes/01_Rhizobia/FribourgSMeliloti_Prokka_v2.gff  -o CountsTable_$(echo $FILE | cut -d'_' -f1,2).txt $FILE -T 8"; sleep  1; done


7. Import count tables in R to calculate the rpkm (reads per kilobase per million for each gene)

8. Calculate RNAseq metrics with CollectRnaSeqMetrics to disproportionate number of reads mapping to the 3’ ends of transcripts
Poly(A) selection is performed via oligo(dT) bead capture, which will lead to 5' regions not being captured for a large proportion of transcripts when degraded RNA is used as input. This ultimately affects accurate transcript abundance estimation, as transcript length is used as a basis for normalizing read counts.

a. transform annotation gff into refflat format with gffread

         /data/users/imateusgonzalez/SOFTS/gffread/gffread -E /data/users/imateusgonzalez/2024_RNAseqTrial/00_ReferenceGenomes/01_Rhizobia/FribourgSMeliloti_Prokka_v2.gff -T -o FribourgSMeliloti_Prokka_v2.refFlat

b. run piccard


         java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I /data/users/imateusgonzalez/2024_RNAseqTrial/03_UnmappedReads/150mg_L1_MappedRhizobiaAligned.sortedByCoord.out.bam -O /data/users/imateusgonzalez/2024_RNAseqTrial/03_UnmappedReads/150mg_L1_RNAseqMetrics.txt --REF_FLAT /data/users/imateusgonzalez/2024_RNAseqTrial/00_ReferenceGenomes/01_Rhizobia/FribourgSMeliloti_Prokka_v2.refFlat -STRAND FIRST_READ_TRANSCRIPTION_STRAND


      

 
