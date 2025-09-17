# Bref4

**Bref4** is a program for compressing and decompressing VCF files with
phased non-missing genotypes. This is an alpha release.  The bref4 general 
release will include additional command line options and capabilities.

Version: 0.1 (alpha release) \
Last updated: September 17, 2025

## Contents

* [Download bref4](#download-bref4)
* [Bref4 Command Line](#bref4-command-line)
* [Input and Output File Parameters](#input-and-output-file-parameters)
* [General Parameters](#general-parameters)
* [License](#license)
* [References](#references)


## Download bref4

You can download **bref4.jar** with the command:

    wget https://faculty.washington.edu/browning/bref4.jar

or you can create **bref4.jar** with the commands:

    git clone https://github.com/browning-lab/bref4.git
    javac --release 11 -cp bref4/src/ bref4/src/bref4/Bref4Main.java
    jar cfe bref4.jar bref4/Bref4Main -C bref4/src/ ./

[Contents](#contents)

## Bref4 Command Line
The bref4 program requires Java version 1.11 or a later Java version.
The command to run bref4 is:

&nbsp; &nbsp; &nbsp; &nbsp; java -Xmx[_GB_]g -jar bref4.jar [_parameter_=_value_] [...]

where [_GB_] is the available GB of memory (e.g. "-Xmx16g"). Bref4 has two
required parameters: [**in**](#input-and-output-file-parameters)
and [**out**](#input-and-output-file-parameters). 
Entering "**java -jar bref4.jar**" at the command prompt will print a brief
description of the command line arguments. Executing the
[run.bref4.test](https://raw.githubusercontent.com/browning-lab/bref4/master/test/run.bref4.test)
shell script will download a small gzip-compressed VCF file, compress the 
VCF file, and decompress the resulting bref4 file.

The bref4 parameters are described in the following two sections.

## Input and Output File Parameters

The **in** and **out** parameters specify the input and output files.

* **in=**_input\_file_  &nbsp; &nbsp; &nbsp; An input file with phased genotypes.

* **out=**_output\_file_ &nbsp; &nbsp; &nbsp; The output file.

An input or output file can be a Variant Call Format (VCF) file[^1] or a bref4 file.
The input or output file type is indicated by its filename:
* Uncompressed VCF files must have names ending in **.vcf**
* Gzip-compressed VCF files must have names ending in **.vcf.gz** or **.vcf.bgz**
* Bref4-compressed VCF files must have names ending in **.bref4**

If _input\_file_ is a hyphen ("**in=-**"), bref4 will read an uncompressed
VCF file from standard input.  Similarly, if _output\_file_ is a hyphen
("**out=-**"), bref4 will write an uncompressed VCF file to standard output

All genotypes in an input VCF file must be nonmissing, phased, and have
the phased allele separator ('|'). You can phase genotypes and impute sporadic
missing genotypes with the
[Beagle](https://faculty.washington.edu/browning/beagle/beagle.html) program.
Gzip-compressed VCF output files can be indexed with the tabix program.[^2]

If _input\_file_ is a VCF file or a gzip-compressed VCF file, bref4
will ignore all non-GT FORMAT fields and their associated sample data. 
All other VCF data will be stored, including the meta-information lines, 
the sample identifiers, the CHROM, POS, ID, REF, ALT, QUAL, FILTER, and 
INFO fields for each VCF data line, and the phased genotypes for each sample. 
Statistics that summarize non-genotype sample data can be stored in the INFO 
field. For example, if each genotype has a GQ (conditional genotype quality) 
value, you can store the average GQ value in the INFO field.

[Contents](#contents)

## General Parameters

* **nthreads=**_integer_ (default: **number of CPU cores**)
&nbsp; &nbsp; &nbsp; Set the maximum number of computational threads to
_integer_ (a positive integer).

[Contents](#contents)

## License
The bref4 software is licensed under the Apache License, Version 2.0
(the License). You may obtain a copy of the License from
[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0).

[Contents](#contents)

## References

[^1]: Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE,
Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project
Analysis Group. The variant call format and VCFtools. Bioinformatics.
2011 Aug 1;27(15):2156-8. doi: 10.1093/bioinformatics/btr330.
PMID: 21653522; PMCID: PMC3137218.

[^2]: Li H. Tabix: fast retrieval of sequence features from generic TAB-delimited
files. Bioinformatics. 2011 Mar 1;27(5):718-9. doi: 10.1093/bioinformatics/btq671.
PMID: 21208982; PMCID: PMC3042176.
