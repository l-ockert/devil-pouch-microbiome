#### TRAIN CLASSIFIER USING RESCRIPt#### 

# download SILVA files (ssu 138.1) #
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz
gunzip tax_slv_ssu_138.1.txt.gz

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz
gunzip taxmap_slv_ssu_ref_nr_138.1.txt.gz

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/tax_slv_ssu_138.1.tre.gz
gunzip tax_slv_ssu_138.1.tre.gz

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
gunzip SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz

# import SILVA files into QIIME2 #
qiime tools import \
--type 'FeatureData[SILVATaxonomy]' \
--input-path SILVA-base-files/tax_slv_ssu_138.1.txt \
--output-path taxranks-silva-138.1-ssu-nr99.qza \

qiime tools import \
--type 'FeatureData[SILVATaxidMap]' \
--input-path SILVA-base-files/taxmap_slv_ssu_ref_nr_138.1.txt \
--output-path taxmap-silva-138.1-ssu-nr99.qza \

qiime tools import \
--type 'Phylogeny[Rooted]' \
--input-path SILVA-base-files/tax_slv_ssu_138.1.tre \
--output-path taxtree-silva-138.1-nr99.qza \

qiime tools import \
--type 'FeatureData[RNASequence]' \
--input-path SILVA-base-files/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta \
--output-path silva-138.1-ssu-nr99-rna-seqs.qza

# reverse transcribe data #
qiime rescript reverse-transcribe \
--i-rna-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
--o-dna-sequences silva-138.1-ssu-nr99-seqs.qza

# combine files to creat reference database #
qiime rescript parse-silva-taxonomy \
--i-taxonomy-tree taxtree-silva-138.1-nr99.qza \
--i-taxonomy-map taxmap-silva-138.1-ssu-nr99.qza \
--i-taxonomy-ranks taxranks-silva-138.1-ssu-nr99.qza \
--o-taxonomy silva-138.1-ssu-nr99-tax.qza

# dereplicate sequences and taxonomy #
qiime rescript dereplicate \
--i-sequences silva-138.1-ssu-nr99-seqs.qza  \
--i-taxa silva-138.1-ssu-nr99-tax.qza \
--p-mode 'uniq' \
--o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
--o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza

# train classifier #
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
--i-reference-taxonomy silva-138.1-ssu-nr99-tax-derep-uniq.qza \
--o-classifier silva-138.1-ssu-nr99-classifier.qza

# make classifier amplicon specific (V3-V4) #
qiime feature-classifier extract-reads \
--i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
--p-f-primer CCTACGGGNGGCWGCAG \
--p-r-primer GACTACHVGGGTATCTAATCC \
--p-n-jobs 2 \
--p-read-orientation 'forward' \
--o-reads silva-138.1-ssu-nr99-seqs-341f-805r.qza

#### IMPORTING AND DENOISING DATA ####

# import paired-end demultiplexed casava 1.8 reads into QIIME2 #
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path raw-fastq-files \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path demultiplexed-files/demux-paired-end.qza

## visualise sequence count summary ##
qiime demux summarize \
--i-data demultiplexed-files/demux-paired-end.qza \
--o-visualization demultiplexed-files/demux-paired-end.qzv

# denoise and create feature table #
qiime dada2 denoise-paired \
--i-demultiplexed-seqs demultiplexed-files/demux-paired-end.qza \
--p-trunc-len-f 288 \
--p-trunc-len-r 256 \
--p-trim-left-f 21 \
--p-n-threads 0 \
--o-representative-sequences feature-table/rep-seqs.qza \
--o-table feature-table/unfiltered-table.qza \
--o-denoising-stats feature-table/stats.qza

## visualise feature table ##
qiime feature-table summarize \
--i-table feature-table/unfiltered-table.qza \
--o-visualization feature-table/unfiltered-table.qzv \
--m-sample-metadata-file devil-metadata.tsv

## visualise rep seqs ## 
qiime feature-table tabulate-seqs \
--i-data feature-table/rep-seqs.qza \
--o-visualization feature-table/rep-seqs.qzv

# generate phylogenetic tree #
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences feature-table/rep-seqs.qza \
--p-n-threads auto \
--o-alignment tree/aligned-rep-seqs.qza \
--o-masked-alignment tree/masked-aligned-rep-seqs.qza \
--o-tree tree/unrooted-tree.qza \
--o-rooted-tree tree/rooted-tree.qza

# assign taxonomy #
qiime feature-classifier classify-sklearn \
--i-classifier classification/silva-138.1-ssu-nr99-341f-805r-classifier.qza \
--i-reads feature-table/rep-seqs.qza \
--p-n-jobs -1 \
--o-classification classification/taxonomy.qza

# incorporate metadata #
qiime metadata tabulate \
--m-input-file classification/taxonomy.qza \
--o-visualization classification/taxonomy.qzv

# visualise taxonomy #
qiime taxa barplot \
--i-table feature-table/unfiltered-table.qza \
--i-taxonomy classification/taxonomy.qza \
--m-metadata-file devil-metadata.tsv \
--o-visualization classification/taxa-bar-plots.qzv

# export files for phyloseq #
qiime tools export \
--input-path feature-table/unfiltered-table.qza \
--output-path phyloseq

biom convert \
-i phyloseq/feature-table.biom \
-o phyloseq/otu_table.txt \
--to-tsv ## manually change column name "#OTU ID" to "OTUID" and remove first row comment

qiime tools export \
--input-path classification/taxonomy.qza \
--output-path phyloseq ## manually change "feature ID" to "OTUID"

qiime tools export \
--input-path tree/unrooted-tree.qza \
--output-path phyloseq

## NOW TO DECONTAM AND ANALYSE IN R...#


