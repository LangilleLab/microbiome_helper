Microbiome Helper
=================

An assortment of scripts to help process and automate various microbiome and metagenomic bioinformatic tools.

Including workflows (see below) for analyzing 16S and shotgun metagenomic data. 

Requirements
------------
The following programs should be installed with commands accessible from the user's PATH, before trying to run any of the scripts included in this repository.

**Both pipelines**
* FastQC (optional): http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* PEAR: http://sco.h-its.org/exelixis/web/software/pear/doc.html
 
**Metagenomics**
* Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
* Human pre-indexed database: ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
* MetaPhlAn2: https://bitbucket.org/biobakery/metaphlan2
* HUMAnN: http://huttenhower.sph.harvard.edu/humann

**16S**
* FASTX toolkit: http://hannonlab.cshl.edu/fastx_toolkit/download.html
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

1. (Optional) Run fastqc to allow manual inspection of the quality of sequences

        mkdir fastqc_out
        fastqc -t 4 raw_miseq_data/* -o fastqc_out/

2. Stich paired end reads together

        run_pear.pl -p 4 -o stitched_reads raw_miseq_data/* >pear.log

3. Run bowtie2 to screen out human sequences (Note: you can use run_deconseq.pl instead but it is much slower)
    
        run_human_filter.pl -p 4 -o screened_reads/ stitched_reads/*.assembled*

4. Run Metaphlan for taxanomic composition

        run_metaphlan2.pl -p 4 -o metaphlan_taxonomy.txt screened_reads/*

5. Convert from metaphlan to stamp profile file

        metaphlan2stamp.pl metaphlan_taxonomy.txt > metaphlan_taxonomy.spf

6. Run pre-humann (diamond search)

        run_pre_humann.pl -p 4 -o pre_humann/ screened_reads/*

7. Run humann (link files to humann "input" directory and then run humann with scons command). Note that you can run this in parallel with `-j` option (e.g. scons -j 4), but I have found this often causes humann to unexpectedly error.

        ln -s $PWD/pre_humann/* ~/programs/humann-0.99/input/
        cd ~/programs/human-0.99/
        scons

8. Convert humann output to stamp format

        humann_to_stamp.py 04b-hit-keg-mpm-cop-nul-nve-nve.txt > hummann_modules.spf
        humann_to_stamp.py 04b-hit-keg-mpt-cop-nul-nve-nve.txt > hummann_pathways.spf
        humann_to_stamp.py 01b-hit-keg-cat.txt > hummann_kos.spf

16S Workflow
------------

*Note that this workflow starts with raw paired-end MiSeq data in demultiplexed fastq format.*

1. Run fastqc to allow manual inspection of the quality of sequences (optional)

        mkdir fastqc_out
        fastqc -t 4 raw_miseq_data/* -o fastqc_out/

2. Stich paired end reads together (~3 hours)

        run_pear.pl -p 4 -o stitched_reads raw_miseq_data/* >pear.log
		
3. Convert FASTQ stitched files to FASTA AND remove any sequences that have an 'N' in them. (~20 minutes)

        run_fastq_to_fasta.pl -p -o fasta_files stitched_reads/*.assembled.*

4. Create a QIIME "map.txt" file with the first column containing the sample names and another column called "FileInput" containing the filenames. This is a tab-delimited file and there must be columns named "BarcodeSequence" and "LinkerPrimerSequence" that are empy. This file can then contain other columns to group samples which will be used when figures are created later.

        create_qiime_map.pl fasta_files/* > map.txt
		
5. Combine files into single QIIME "seqs.fna" file (~5 minutes)

        add_qiime_labels.py -i fasta_files/ -m map.txt -c FileInput -o combined_fasta
		
6. Create OTU picking parameter file. (Note: this is optional)

        echo "pick_otus:sortmerna_db /home/shared/sortmerna/97_otus" >> ucrss_smr_suma_params.txt
        echo "pick_otus:threads 4" >> ucrss_smr_suma_params.txt
		

7. Run the entire qiime open reference picking pipeline with the new sortmerna (for reference picking) and sumaclust (for de novo otu picking). This does reference picking first, then subsamples failure sequences, de-novo otu picks failures, ref picks against de novo otus, and de-novo picks again any left over failures. Note: the last de-novo picking step may have to be skipped if there are too many sequences by adding the option `--suppress_step4`. Note: You may want to change the subsampling percentage to a higher amount from the default -s 0.001 to -s 0.01 (e.g 1% of the failures) or -s 0.1 (e.g. 10% of the failures) (~24 hours)

        pick_open_reference_otus.py -i $PWD/combined_fasta/combined_seqs.fna -o $PWD/ucrss_sortmerna_sumaclust/ -p $PWD/ucrss_smr_suma_params.txt -m sortmerna_sumaclust -s 0.1 --suppress_step4 -v

8. Normalize OTU table to same sample depth (e.g. in this case 35566 sequences, but this value will depend on your OTU table)

        single_rarefaction.py -i ucrss_sortmerna_sumaclust/otu_table_mc2_w_tax_no_pynast_failures.biom -o final_otu_tables/otu_table.biom -d 35566


9. Create unifrac beta diversity plots

        beta_diversity_through_plots.py -m map.txt -t ucrss_sortmerna_sumaclust/rep_set.tre -i final_otu_tables/otu_table.biom -o plots/bdiv_otu

10. Create alpha diversity rarefaction plot

        alpha_rarefaction.py -i final_otu_tables/otu_table.biom -o plots/alpha_rarefaction_plot -t ucrss_sortmerna_sumaclust/rep_set.tre -m map.txt --min_rare_depth 1000 --max_rare_depth 35000 --num_steps 35

11. Convert BIOM otu table to STAMP

        biom convert -i final_otu_tables/otu_table.biom -o final_otu_tables/otu_table_json.biom --to-json --table-type "OTU table"
        
        #Note: must not be in virtualenv during this next command since it uses biom 1.3.1
        biom_to_stamp.py -m taxonomy final_otu_tables/otu_table_json.biom >final_otu_tables/otu_table.spf

        #Note: there is are a few OTUs where the genus Clostridium is within the wrong family. Filter these out here.
        grep -P -v "f__Erysipelotrichaceae\tg__Cl" otu_table.spf > tmp.spf
        mv tmp.spf otu_table.spf


PICRUSt workflow (for 16S data)
----------------

1. First reduce OTU table to just the reference OTUs

        filter_otus_from_otu_table.py -i final_otu_tables/otu_table.biom -o final_otu_tables/closed_otus.biom --negate_ids_to_exclude -e ~/.virtualenvs/qiime1.9/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta

2. Convert from hdf5 to json

        biom convert -i final_otu_tables/closed_otus.biom -o final_otu_tables/closed_otus_json.biom --to-json --table-type "OTU table"

3. PICRUSt normalize OTU table by predicted 16S copy numbers (must change to python evironment with biom 1.3.1)

        normalize_by_copy_number.py -i final_otu_tables/closed_otus_json.biom -o final_otu_tables/closed_otus_norm.biom

4. PICRUSt predict KOs .... this is magical :)

        predict_metagenomes.py -i final_otu_tables/closed_otus_norm.biom -o final_otu_tables/ko.biom

5. Collapse KOs into KEGG Pathways

        categorize_by_function.py -i ko.biom -l 3 -c KEGG_Pathways -o ko_L3.biom

6. Convert BIOM to STAMP format

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

* **metaphlan_to_stamp.pl**: This is script that converts a merged metaphlan output file that was created using all taxnomic ranks to a STAMP profile file. 

* **biom_to_stamp.py**: STAMP has built-in BIOM conversion, but depending on where the BIOM file comes from there can be slight format problems with STAMP. Specifically, this script handles PICRUSt BIOM output (both KOs and KEGG Pathways) and QIIME 16S OTU tables (with metadata 'taxonomy').

* **humann_to_stamp.pl**: Converts humann output to stamp format (currently removes extra rows and renames samples ids so they are the same as the orginal file)

***MetaPhlan (Metagenome taxonomic annotation)***

* **run_metaphlan.pl**: Wraps the Metaphlan package and handles running multiple samples at once as well as handling paired end data a bit more cleaner. It also runs each sample in parallel and merges the results into a single output file. Also, easily allows gzipped or non-gzipped files.

***Humann (Metagenome functional annotation)***

* **run_pre_humann.pl**: Does similarity search against kegg database using search tool diamond using multiple threads. This output is then fed into humann. 

***Other assorted scripts***

* **run_fastq_to_fasta.pl**: Wraps the fastq_to_fasta command from the FASTX Toolkit to allow the use of multiple threads.


Contact
-------

Questions or comments about this repository can be sent to morgan.g.i.langille@dal.ca

