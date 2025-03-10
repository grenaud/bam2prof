
# bam2prof

bam2prof is a tool designed to analyze BAM files and generate substitution profiles, particularly useful for assessing ancient DNA damage patterns. This standalone version is derived from the subprogram in schmutzi (https://github.com/grenaud/schmutzi) but has been rewritten to utilize htslib (https://github.com/samtools/htslib).

## Features

- Substitution Profiling: Computes substitution rates at both 5' and 3' ends of reads to assess DNA damage patterns.
- Customizable Parameters: Allows users to set minimum base quality scores, minimum read lengths, and specify the length of the profile to generate.

## Requirements

- htslib: Ensure that htslib (https://github.com/samtools/htslib) is installed on your system, as bam2prof depends on it for BAM file processing.

## Installation

To build bam2prof, follow these steps:

1. Clone the Repository:

git clone https://github.com/grenaud/bam2prof.git

2. Navigate to the Source Directory:

cd bam2prof/src

3. Compile the Program:

make

This will generate the bam2prof executable in the src directory.

## Usage

The general syntax for running bam2prof is:

bam2prof [options] <input.bam>

Key Options:

- -minq <int>: Set the minimum base quality score to consider. Bases with quality scores below this threshold will be ignored.
- -minl <int>: Define the minimum read length to process. Reads shorter than this length will be skipped.
- -length <int>: Specify the length of the profile to generate. This determines how many bases from the ends of reads are analyzed.
- -5p <file>: Output file for the 5' end substitution profile.
- -3p <file>: Output file for the 3' end substitution profile.

Example:

To analyze a BAM file with a minimum base quality of 30, minimum read length of 35, and generate profiles of length 10 for both 5' and 3' ends:

bam2prof -minq 30 -minl 35 -length 10 -5p output_5p.prof -3p output_3p.prof input.bam

This command will produce two files:

- output_5p.prof: Contains the substitution profile for the 5' end.
- output_3p.prof: Contains the substitution profile for the 3' end.

## Notes

- Ensure that your BAM files are properly indexed and that the reference genome used for alignment is accessible.
- For optimal results, consider subsampling your BAM file to a manageable size before running bam2prof. This can help in determining the ideal parameters for your specific dataset.

For more detailed information and updates, visit the bam2prof GitHub repository (https://github.com/grenaud/bam2prof).


# Developers 

- Gabriel Renaud
- Louis Kraft
- Thorfinn Korneliussen
