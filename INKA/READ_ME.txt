Description for using the InKA pipe-line
========================================

The InKA pipeline consists of the following 
essential components or necessary code to
construct:

read_MQ_tables.R
inka.R
collect_output_tables.R

evquant - ELF executable (Source included)

complete-HGNC-21Apr2016.csv
Manning.Kinases_20Jan2017.Rdata
phomics_global_output.txt
PSP_human_03Jul2016.Rdata
PSP_KSR_human_03Jul2016.Rdata
TOP_NWK_nwk_proteome_unprot_ref_2014.Rdata
UniProt_data_08Jun2016.Rdata
valid_mapping.txt

make_intensity_pp.R

Running InKA 
============

The easiest way of running InKA is by copying the following 
required input files from a suitable MaxQuant search, into
the code directory:

evidence.txt
experimentalDesignTemplate.txt
modificationSpecificPeptides.txt
Phospho (STY)Sites.txt

Examples files are provided with the code
in the subdirectory example_data.

Execution of the scripts

./read_MQ_tables.R
./inka.R

This will create amongst others a file called "inka_out.pdf" 
containing the output of the analysis.

A utility script

./collect_output_tables.R

can be called to construct tables containing the InKA
scores for all kinases in all samples, and corresponding tables
for each of the submetrics that compose the InKA score.

For phosphoirylated tyrosine anti-body capture derived
data, we recommend using the scripts as shown above. For TiOx 
and other non-exclusively phospho tyrosine enriched samples, using 

./inka.R --tiox

might improve result presentation and also run much faster, as 
for serine/threonine phosphorylation both size and KS-relation volume 
can be considerably larger.

A second utility script

make_intensity_pp.R

is provided to be run in the following command sequence

./read_MQ_tables.R
./make_intensity_pp.R
mv pp-Peptide-report_intensities.txt  pp-Peptide-report.txt
./inka.R

resulting in an MS1 signal intensity-based InKA analysis.

Additional flags recognized by the inka.R script are

./inka.R --no.network 

which wiil supress the kinase substrate network generation
in the output.

Additional information
======================

Note that the above description assumes that the files

collect_output_tables.R
inka.R
make_intensity_pp.R
read_MQ_tables.R

are executable (Execute "chmod +x *.R" if this is not the case).

Note that the executable evquant is necessary and called from
./read_MQ_tables.R in order to speed up part of the R code
several orders of magnitude (>20K). It needs to be compiled from the 
following provided sources in the src directory:

bmslib_0.10.1.tar.gz
evquant_0.5.0.tar.gz

Note that the default library installation location is ~/lib. If this directory is
not present create it (execute "mkdir ~/lib"). Note that additional
libraries might be required for full use of the bmslib library.

----
Author: Alex Henneman
Date  : Mon Aug 27 12:04:37 CEST 2018
