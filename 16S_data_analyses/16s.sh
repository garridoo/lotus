#!/bin/bash

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# exits whenever a function returns 1
set -e

# get path to scripts
scripts_dir=$(dirname $0)

# load config file
config_file=$1
source $config_file

# load functions
source $scripts_dir/16s.functions.sh

# activate QIIME, etc.
source $scripts_dir/activate.sh

# cleanup
cleanup

# processing of the raw reads for all libraries:
# joining of pair reads, cuality control, demultiplexing

if [ "$l_list" != "" ]
then

    rm -f $working_dir/seqs.fasta

    for l in $l_list
    do    

        rm -f -r $working_dir/"$l"
        mkdir -p $working_dir/"$l"

        # merge paired reads
        log "["$l"] joining paired-end reads..."
        add_barcode_to_label.py $data_dir/"$l"_forward_reads.fastq.gz \
                                $data_dir/"$l"_barcodes.fastq.gz \
                                $working_dir/"$l"/forward_reads.fastq \
                                &>> $output
        add_barcode_to_label.py $data_dir/"$l"_reverse_reads.fastq.gz \
                                $data_dir/"$l"_barcodes.fastq.gz \
                                $working_dir/"$l"/reverse_reads.fastq \
                                &>> $output
        usearch -fastq_mergepairs $working_dir/"$l"/forward_reads.fastq \
                -reverse $working_dir/"$l"/reverse_reads.fastq \
                -fastqout $working_dir/"$l"/joined.fastq \
                &>> $output
        extract_barcodes.py -f $working_dir/"$l"/joined.fastq \
                            -c barcode_in_label \
                            --char_delineator 'BC=' \
                            --bc1_len 12 \
                            -o $working_dir/"$l"/ \
                            &>> $output

        # with fastq-join
        # join_paired_ends.py -f $data_dir/"$l"_forward_reads.fastq.gz \
        #              -r $data_dir/"$l"_reverse_reads.fastq.gz \
        #              -b $data_dir/"$l"_barcodes.fastq.gz \
        #              -o $working_dir/$l \
        #              &>> $output 
        # mv $working_dir/"$l"/fastq.join.fastq $working_dir/"$l"/joined.fastq
        # mv $working_dir/"$l"/fastq.join_barcodes.fastq $working_dir/"$l"/barcodes.fastq

        # demultiplex
        log "["$l"] demultiplexing and quality filtering..."
        split_libraries_fastq.py -i $working_dir/$l/joined.fastq \
                                 -b $working_dir/$l/barcodes.fastq \
                                 --rev_comp_mapping_barcodes \
                                 --max_barcode_errors $bc_err \
                                 -o $working_dir/$l \
                                 -m $data_dir/"$l"_mapping.txt \
                                 -q $phred \
                                 --phred_offset 33 \
                                 &>> $output

        # save and edit barcode label identifier for usearch compatibility
        rm -f $working_dir/$l/seqs.fasta
        cat $working_dir/$l/seqs.fna | \
            awk '/>/ {print $0, "barcodelabel="$1} !/>/ {print $0}' | \
            sed 's/=>/=/g;s/_[0-9]*$/;/g;s/ /;/g' >> \
            $working_dir/$l/seqs.fasta

        # concatenate all demultiplexed sequences
        cat $working_dir/$l/seqs.fasta >> $working_dir/seqs.fasta

        # generate table of sample sizes
        sampleSizes $working_dir/seqs.fasta \
                    $data_dir/"$l"_mapping.txt \
                    >> $working_dir/sample_sizes.txt

        mv $working_dir/"$l"/split_library_log.txt $working_dir/"$l"_split_library_log.txt

        rm -f -r $working_dir/"$l"

    done

fi

# NOTE: for projects with samples spread across several independent
# sequencing runs, repeat the above steps for each run and then concatenate
# all the labeled sequences into one seqs.fata file before continuing

# dereplication
log "dereplicating..."
usearch -derep_fulllength $working_dir/seqs.fasta \
        -fastaout $working_dir/seqs_unique.fasta \
        -sizeout \
        &>> $output

# abundance sort and discard singletons
log "sorting by abundance and discarding singeltons..."
usearch -sortbysize $working_dir/seqs_unique.fasta \
        -fastaout $working_dir/seqs_unique_sorted.fasta \
        -minsize $min_size \
        &>> $output

# OTU clustering
log "OTU clustering using UPARSE..."
usearch -cluster_otus $working_dir/seqs_unique_sorted.fasta \
        -otus $working_dir/otus.fasta \
        -id $id_threshold \
        &>> $output

# chimera detection
log "removing chimeras..."
usearch -uchime_ref $working_dir/otus.fasta \
        -db $gold_db \
        -strand plus \
        -nonchimeras $working_dir/otus_nc.fasta \
        -threads $n_cores \
        &>> $output

# align sequences to database using PyNAST and remove remaining
log "aligning OTU representative sequences to database..."
align_seqs.py -i $working_dir/otus_nc.fasta \
              -t $gg_core_aligned_db \
              -p $min_id_aln \
              -o $working_dir

sed -i 's/-//g' $working_dir/otus_nc_aligned.fasta

# rename OTUs and remove alignment gaps
log "renaming OTUs..."

cat $working_dir/otus_nc_aligned.fasta | \
    awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' \
    >> $working_dir/rep_seqs.fasta

# generate OTU table
log "generating OTU table..."
usearch -usearch_global $working_dir/seqs.fasta \
        -db $working_dir/rep_seqs.fasta \
        -strand plus \
        -id $id_threshold \
        -uc $working_dir/read_mapping.uc \
        &>> $output

# convert UC file to txt
log "converting UC OTU table file into text format..."
python $usearch_dir/uc2otutab.py $working_dir/read_mapping.uc \
    1> $working_dir/otu_table.txt \
    2>> $output

# taxonomy assignment
log "taxonomy assignment..."
assign_taxonomy.py -i $working_dir/rep_seqs.fasta \
                   -r $gg_core_db \
                   -t $gg_taxonomy \
                   -m $cls_method \
                   -o $working_dir/tax \
                   &>> $output

# cleanup
mv $working_dir/tax/rep_seqs_tax_assignments.txt $working_dir/taxonomy.txt
rm -f -r $working_dir/tax

# parse greengenes taxonomy table
sed -i 's/;/\t/g;s/ //g' $working_dir/taxonomy.txt

# convert OTU table to biom
log "converting OTU table to QIIME compatible biom format..."
biom convert -i $working_dir/otu_table.txt \
             -o $working_dir/otu_table.biom \
             --table-type="OTU table" \
             --to-json \
             &>> $output

# align the representative sequences
log "aligning representative sequences..."
clustalo --seqtype=DNA \
         --threads=$n_cores \
         -i $working_dir/rep_seqs.fasta \
         -o $working_dir/rep_seqs_aligned.fasta \
         --full \
         &>> $output

# filter the alignment
filter_alignment.py -i $working_dir/rep_seqs_aligned.fasta \
                    -o $working_dir \
                    &>> $output

# generate tree from alignment using FastTree
log "generating tree..."
make_phylogeny.py -i $working_dir/rep_seqs_aligned_pfiltered.fasta \
                  -o $working_dir/rep_seqs.tree \
                  &>> $output

### alpha-diversity

log "calculating indices of alpha-diversity..."

# perform rarefaction (only for alpha-diversity analysis)

single_rarefaction.py -i $working_dir/otu_table.biom \
                      -o $working_dir/otu_table_raref.biom \
                      -d $subsampling_depth \
                      &>> $output

# calculate alpha-diversity indices

alpha_diversity.py -i $working_dir/otu_table_raref.biom \
                   -o $working_dir/shannon.txt \
                   -m shannon \
                   &>> $output

alpha_diversity.py -i $working_dir/otu_table_raref.biom \
                   -o $working_dir/chao.txt \
                   -m chao1 \
                   &>> $output

alpha_diversity.py -i $working_dir/otu_table_raref.biom \
                   -o $working_dir/observed_otus.txt \
                   -m observed_species \
                   &>> $output

# extract rarefied OTU table to text format
biom convert -i $working_dir/otu_table_raref.biom \
             --table-type="OTU table" \
             --to-tsv \
             -o $working_dir/otu_table_raref.txt \
             &>> $output

sed -i '/# Const/d;s/#OTU ID.//g' $working_dir/otu_table_raref.txt \

# normalize OTU table
log "normalizing OTU table using the CSS method..."
$scripts_dir/normalize_otu_table.R $working_dir/otu_table.biom \
                        $working_dir/otu_table_norm.biom \
                        &>> $output

# extract normalized OTU table to text format
biom convert -i $working_dir/otu_table_norm.biom \
             --table-type="OTU table" \
             --to-tsv \
             -o $working_dir/otu_table_norm.txt \
             &>> $output

sed -i '/# Const/d;s/#OTU ID.//g' $working_dir/otu_table_norm.txt

### beta-diversity

log "calculating beta-diversity meassures..."

beta_diversity.py -i $working_dir/otu_table_norm.biom \
                  -m bray_curtis \
                  -o $working_dir \
                  &>> $output

mv $working_dir/bray_curtis_otu_table_norm.txt $working_dir/bray_curtis.txt
sed -i 's/^\t//g' $working_dir/bray_curtis.txt

beta_diversity.py -i $working_dir/otu_table_norm.biom \
                  -m weighted_unifrac \
                  -o $working_dir \
                  -t $working_dir/rep_seqs.tree \
                  &>> $output

mv $working_dir/weighted_unifrac_otu_table_norm.txt $working_dir/weighted_unifrac.txt
sed -i 's/^\t//g' $working_dir/weighted_unifrac.txt

beta_diversity.py -i $working_dir/otu_table_norm.biom \
                  -m unweighted_unifrac \
                  -o $working_dir \
                  -t $working_dir/rep_seqs.tree \
                  &>> $output

mv $working_dir/unweighted_unifrac_otu_table_norm.txt $working_dir/unweighted_unifrac.txt
sed -i 's/^\t//g' $working_dir/unweighted_unifrac.txt

# parse OTU TXT tables

sed -i 's/OTUId.//g' $working_dir/otu_table.txt
sed -i 's/OTUId.//g' $working_dir/otu_table_norm.txt
sed -i 's/.OTU.ID.//g' $working_dir/otu_table_raref.txt

log "DONE!"
 
