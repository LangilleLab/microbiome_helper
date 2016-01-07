Microbiome Helper
=================

An assortment of scripts to help process and automate various microbiome and metagenomic bioinformatic tools.

Including workflows (see below) for analyzing 16S and shotgun metagenomic data. 
Table of Contents
=================

  * [Microbiome Helper](#microbiome-helper)
  * [Table of Contents](#table-of-contents)
    * [Requirements](#requirements)
    * [Before getting started](#before-getting-started)
    * [Metagenomics Workflow](#metagenomics-workflow)
    * [16S Workflow](#16s-workflow)
    * [Additional QIIME analysis](#additional-qiime-analysis)
    * [PICRUSt workflow (for 16S data)](#picrust-workflow-for-16s-data)
    * [Brief description of scripts](#brief-description-of-scripts)
    * [Contact](#contact)


Requirements
------------
The following programs should be installed with commands accessible from the user's PATH, before trying to run any of the scripts included in this repository.

**Both pipelines**
* FastQC (v0.11.2) (optional): http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* PEAR: http://sco.h-its.org/exelixis/web/software/pear/doc.html
 
**Metagenomics**
* Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
* Human pre-indexed database: ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
* DIAMOND (> v.0.7.0): http://ab.inf.uni-tuebingen.de/software/diamond/
* MetaPhlAn2: https://bitbucket.org/biobakery/metaphlan2
* HUMAnN: http://huttenhower.sph.harvard.edu/humann

**16S**
* FASTX toolkit (v0.0.14): http://hannonlab.cshl.edu/fastx_toolkit/download.html
* BBMap (v35.59): http://sourceforge.net/projects/bbmap 
* USEARCH (v6.1.544): http://www.drive5.com/usearch/
* QIIME (v1.9): http://qiime.org
* SortMeRNA: http://bioinfo.lifl.fr/RNA/sortmerna/
* SUMACLUST: http://metabarcoding.org/sumatra
* PICRUSt: http://picrust.github.io/picrust/

**Visualization**
* STAMP: http://kiwi.cs.dal.ca/Software/STAMP

Before getting started
----------------------
* These workflows starts with raw paired-end MiSeq data in demultiplexed fastq format.
* Most of the commands allow parallel threads to be used with the `-p` option. Adjust according to your computer hardware abilities.
* These workflows can be completed in roughly 24 hours on a single quad-core desktop when starting with 25 million paired end reads.

Metagenomics Workflow
---------------------

1. (Optional) Run FastQC to allow manual inspection of the quality of sequences.

        mkdir fastqc_out
        fastqc -t 4 raw_data/* -o fastqc_out/

2. Stich paired end reads together (summary of stitching results are written to "pear_summary_log.txt").

        run_pear.pl -p 4 -o stitched_reads raw_data/*

3. Run Bowtie2 to screen out human sequences (Note: you can use run_deconseq.pl instead but it is much slower).
    
        run_human_filter.pl -p 4 -o screened_reads/ stitched_reads/*.assembled*

4. Run MetaPhlAn2 for taxonomic composition.

        run_metaphlan2.pl -p 4 -o metaphlan_taxonomy.txt screened_reads/*

5. Convert from MetaPhlAn to STAMP profile file.

        metaphlan_to_stamp.pl metaphlan_taxonomy.txt > metaphlan_taxonomy.spf

6. Run pre-HUMAnN (DIAMOND search).

        run_pre_humann.pl -p 4 -o pre_humann/ screened_reads/*

7. Run HUMAnN (link files to HUMAnN "input" directory and then run HUMAnN with scons command). Note that you can run this in parallel with `-j` option (e.g. scons -j 4), but I have found this often causes HUMAnN to unexpectedly error.

        ln -s $PWD/pre_humann/* ~/programs/humann-0.99/input/
        cd ~/programs/humann-0.99/
        scons

8. Convert HUMAnN output to STAMP format

        humann_to_stamp.pl 04b-hit-keg-mpm-cop-nul-nve-nve.txt > hummann_modules.spf
        humann_to_stamp.pl 04b-hit-keg-mpt-cop-nul-nve-nve.txt > hummann_pathways.spf
        humann_to_stamp.pl 01b-hit-keg-cat.txt > hummann_kos.spf

16S Workflow
------------

*Note that this workflow starts with raw paired-end MiSeq data in demultiplexed fastq format assumed to be located within a folder called `raw_data`*

1. (Optional) Run FastQC to allow manual inspection of the quality of sequences

        mkdir fastqc_out
        fastqc -t 4 raw_data/* -o fastqc_out/

2. Stich paired end reads together (summary of stitching results are written to "pear_summary_log.txt")

        run_pear.pl -p 4 -o stitched_reads raw_data/* 

3. Filter stitched reads by quality score, length and ensure forward and reverse primers match each read (summary written to "readFilter_log.txt" by default).

        readFilter.pl -q 30 -p 90 -l 400 stitched_reads/*.assembled.*
									
4. Convert FASTQ stitched files to FASTA AND remove any sequences that have an 'N' in them.

        run_fastq_to_fasta.pl -p -o fasta_files filtered_reads/*

5. Remove chimeric sequences with UCHIME (summary written to "chimeraFilter_log.txt" by default).

        chimeraFilter.pl -type 1 -db /usr/local/db/single_strand/Bacteria_RDP_trainset15_092015.udb fasta_files/*	

6. Create a QIIME "map.txt" file with the first column containing the sample names and another column called "FileInput" containing the filenames. This is a tab-delimited file and there must be columns named "BarcodeSequence" and "LinkerPrimerSequence" that are empy. This file can then contain other columns to group samples which will be used when figures are created later.

        create_qiime_map.pl non_chimeras/* > map.txt
		
7. Combine files into single QIIME "seqs.fna" file (~5 minutes).

        add_qiime_labels.py -i non_chimeras/ -m map.txt -c FileInput -o combined_fasta
		
8. Create OTU picking parameter file.

        echo "pick_otus:threads 4" >> clustering_params.txt
        echo "pick_otus:sortmerna_coverage 0.8" >> clustering_params.txt
        
9. Run the entire qiime open reference picking pipeline with the new sortmerna (for reference picking) and sumaclust (for de novo OTU picking). This does reference picking first, then subsamples failure sequences, de-novo OTU picks failures, ref picks against de novo OTUs, and de-novo picks again any left over failures. Note: You may want to change the subsampling percentage to a higher amount from the default -s 0.001 to -s 0.01 (e.g 1% of the failures) or -s 0.1 (e.g. 10% of the failures) (~24 hours).

        pick_open_reference_otus.py -i $PWD/combined_fasta/combined_seqs.fna -o $PWD/clustering/ -p $PWD/clustering_params.txt -m sortmerna_sumaclust -s 0.1 -v --min_otu_size 1 

10. Filter OTU table to remove singletons as well as low-confidence OTUs that are likely due to MiSeq bleed-through between runs (reported by Illumina to be 0.1% of reads). 

        remove_low_confidence_otus.py -i $PWD/clustering/otu_table_mc1_w_tax_no_pynast_failures.biom -o $PWD/clustering/otu_table_high_conf.biom

11. Summarize OTU table to determine number of sequences per sample.

        biom summarize-table -i clustering/otu_table_high_conf.biom -o clustering/otu_table_high_conf_summary.txt

12. Normalize OTU table to same sample depth - you will need to change the value of X shown below to match the read count of the sample with the lowest (acceptable) number of reads. Note: Don't like the idea of throwing away all that data? You may want to consider trying different normalization methods such as DESeq2 (see below).

        mkdir final_otu_tables
        single_rarefaction.py -i clustering/otu_table_high_conf.biom -o final_otu_tables/otu_table.biom -d X

13. Manually add column(s) to map.txt that contain information to group your samples (e.g. healthy vs disease).

14. Create UniFrac beta diversity plots.

        beta_diversity_through_plots.py -m map.txt -t clustering/rep_set.tre -i final_otu_tables/otu_table.biom -o plots/bdiv_otu

15. Create alpha diversity rarefaction plot (values min and max rare depth as well as number of steps should be based on the number of sequences within your OTU table).

        alpha_rarefaction.py -i final_otu_tables/otu_table.biom -o plots/alpha_rarefaction_plot -t clustering/rep_set.tre -m map.txt --min_rare_depth 1000 --max_rare_depth 35000 --num_steps 35

16. Convert BIOM OTU table to STAMP.
        
        biom_to_stamp.py -m taxonomy final_otu_tables/otu_table.biom >final_otu_tables/otu_table.spf

        #Note: there are a few OTUs where the genus Clostridium is within the wrong family. Filter these out here manually.
        grep -P -v "f__Erysipelotrichaceae\tg__Cl" otu_table.spf > tmp.spf
        mv tmp.spf otu_table.spf

17. Add sample metadata to BIOM file so that it can be used by other tools like phinch.org and phyloseq.

        biom add-metadata -i final_otu_tables/otu_table.biom -o final_otu_tables/otu_table_with_metadata.biom -m map.txt


Additional QIIME analysis
-------------------------

*Once you have created an OTU table using the above pipeline, there are several statistcal tests and visualization methods available within QIIME. Not all of these are useful for all data types or projects. This list is meant more as a catalog of what is possible.*

* Get statistical signfance of category groupings (be sure to read the help for this command to properly interpret the results). (Be sure to change the `-c` option to your particular field within your map.txt file).

        compare_categories.py -i final_otu_tables/otu_table.biom --method adonis -i plots/bdiv_otu/weighted_unifrac_dm.txt -m map.txt -c mouse_type -o adonis_ko_vs_wt
        compare_categories.py -i final_otu_tables/otu_table.biom --method anosim -i plots/bdiv_otu/weighted_unifrac_dm.txt -m map.txt -c mouse_type -o anosim_ko_vs_wt

* Compare within vs between b-diversity.

        make_distance_boxplots.py -m map.txt -d plots/bdiv_otu/weighted_unifrac_dm.txt -f mouse_type -o plots/bdiv_box_plots

* Make stacked bar charts.

        summarize_taxa_through_plots.py -i final_otu_tables/otu_table.biom -o plots/taxa_summary
		
		#Note: you can also collapse samples by a category
        summarize_taxa_through_plots.py -i final_otu_tables/otu_table.biom -o plots/taxa_summary -m map.txt -c mouse_type
		
* Compute MD index for IBD.

        compute_taxonomy_ratios.py -i final_otu_tables/otu_table.biom -e md -o map_with_md.txt -m map.txt

* Build and evaluate a classifier using Random Forests (this should only be done for large datasets).

        supervised_learning.py -i otu_table.biom -m map.txt -c mouse_type -o ml -v
		
		#Can also do 5-fold cross validation
		supervised_learning.py -i otu_table.biom -m map.txt -c mouse_type -o ml -v -e cv5

* Normalize OTU table using DESeq2 instead of rarefying (this method of normalization is still not mainstream. Check out the help of this script (normalize_table.py -h) to read more about these methods. 

        normalize_table.py -i clustering/otu_table_high_conf.biom -a DESeq2 -z -o final_otu_tables/otu_table_deseq2_norm_no_negatives.biom
		
PICRUSt workflow (for 16S data)
----------------

1. First reduce OTU table to just the reference OTUs.

        filter_otus_from_otu_table.py -i final_otu_tables/otu_table.biom -o final_otu_tables/closed_otus.biom --negate_ids_to_exclude -e /usr/local/lib/python2.7/dist-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta

2. Convert from BIOM encoding from HDF5 to JSON.

        biom convert -i final_otu_tables/closed_otus.biom -o final_otu_tables/closed_otus_json.biom --to-json --table-type "OTU table"

(Temporary fix) Install BIOM 1.3.1 and temporarily change your PYTHONPATH

        pip install --target ~/lib biom-format==1.3.1
        export PYTHONPATH=~/lib/

3. PICRUSt: normalize OTU table by predicted 16S copy numbers. NOTE: PICRUSt has not been updated yet for BIOM 2.1. Therefore you must change to a python environment with biom 1.3.1 for the next 3 commands).


        normalize_by_copy_number.py -i final_otu_tables/closed_otus_json.biom -o final_otu_tables/closed_otus_norm.biom

4. PICRUSt: Predict KOs .... this is magical :)

        predict_metagenomes.py -i final_otu_tables/closed_otus_norm.biom -o final_otu_tables/ko.biom

5. PICRUSt: Collapse KOs into KEGG Pathways.

        categorize_by_function.py -i ko.biom -l 3 -c KEGG_Pathways -o ko_L3.biom

(Temporary Fix) Revert back to newer BIOM version by exiting and logging back into the shell or:

        export PYTHONPATH=''

6. Convert BIOM to STAMP format.

        biom_to_stamp.py -m KEGG_Pathways ko_L3.biom > ko_L3.spf


Brief description of scripts
----------------------------

***Human Contamination***
        
* **run_deconseq.pl**: Wraps the deconseq program to filter out human reads from metagenomic data

* **run_human_filter.pl**: Wraps the Bowtie2 program to filter out human reads from metagenomic data (faster than deconseq) 

***PEAR (Paired-end stitching)***

* **run_pear.pl**: Makes running the PEAR program easier on many samples, by automatically identifying the paired-end files to merge together. 

***STAMP (Statistics and Visualization)***

I find STAMP to be a very useful tool for creating figures and doing statistical analyses. Here are scripts to covert tables from different tools into STAMP input profile files.

* **metaphlan_to_stamp.pl**: This is script that converts a merged MetaPhlAn output file that was created using all taxnomic ranks to a STAMP profile file. 

* **biom_to_stamp.py**: STAMP has built-in BIOM conversion, but depending on where the BIOM file comes from there can be slight format problems with STAMP. Specifically, this script handles PICRUSt BIOM output (both KOs and KEGG Pathways) and QIIME 16S OTU tables (with metadata 'taxonomy').

* **humann_to_stamp.pl**: Converts HUMAnN output to STAMP format (currently removes extra rows and renames samples ids so they are the same as the orginal file)

***MetaPhlan (Metagenome taxonomic annotation)***

* **run_metaphlan.pl**: Wraps the MetaPhlAn package and handles running multiple samples at once as well as handling paired end data a bit more cleaner. It also runs each sample in parallel and merges the results into a single output file. Also, easily allows gzipped or non-gzipped files.

***Humann (Metagenome functional annotation)***

* **run_pre_humann.pl**: Does similarity search against kegg database using search tool diamond using multiple threads. This output is then fed into HUMAnN. 

***Other assorted scripts***

* **run_fastq_to_fasta.pl**: Wraps the fastq_to_fasta command from the FASTX Toolkit to allow the use of multiple threads.

* **chimeraFilter.pl**: Wraps UCHIME (implemented in USEARCH) to filter out chimeric reads from a directory of reads in fasta format.  

* **readFilter.pl**: Wraps several read filtering commands together (using FASTX Toolkit and BBMap) to run on a directory of fastq files. 

* **remove_low_confidence_otus.py**: Filters a BIOM table to remove low confidence OTUs that result from MiSeq run-to-run bleed-through (based on 0.1% as reported by Illumina).   


Contact
-------

Questions or comments about this repository can be sent to morgan.g.i.langille@dal.ca


