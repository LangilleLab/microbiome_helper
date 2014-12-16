Microbiome Helper
=================

An assortment of scripts to help process and automate various microbiome and metagenomic bioinformatic tools.

As I analyze various microbiome datasets (16S, metagenomics, etc.), I tend to write scripts to automate processes or to convert file formats between tools. 

These are the collection of such tools. 

STAMP
-----

I find STAMP to be a very useful tool for creating figures and doing statistical analyses. Here are scripts to covert tables from different tools into STAMP input profile files.

metaphlan_to_stamp.pl -- This is a simple script that converts a merged metaphlan output file to a STAMP profile file.

biom_to_stamp.py -- STAMP has built-in BIOM conversion, but depending on where the BIOM file comes from there can be slight format problems with STAMP. Specifically, this script handles PICRUSt BIOM output (both KOs and KEGG Pathways) and QIIME 16S OTU tables (with metadata 'taxonomy').

MetaPhlan
---------

run_metaphlan.pl  -- Wraps the metaphlan.py package and handles running multiple samples at once as well as handling paired end data a bit more cleaner. It also runs each sample in parallel and merges the results into a single output file.


