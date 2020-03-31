# Functions to extract and plot fastqc data

fastqc is a frequently used program to evaluate the quality of
next-generation sequencing data. It creates a separate html page and
zip file per fastq file. This is useful when you have have small
number of sequence files, but is incovenient when you want to
compare the qualities from different files. 

The functions provided here makes it easy to incorporate the fastq
data into downstream R analyses and in addition to simply
pre-screening the data it is also possible to look back
retrospectively at the quality scores.
