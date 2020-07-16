# Blast Filter

This script remove redundancy of blast output, based on qseqid, and a range of qstart and qend, the match with best bitscore is maintained, besides that, a column called sense its created to indicates the DNA sense of match.

## Dependencies

This script was build on python 3.6.5+ and have these dependencies:

- [pandas](https://pandas.pydata.org/);
- [numpy](https://numpy.org/);
- [argparse](https://docs.python.org/3/library/argparse.html);

## Usage

- python blast_filter.py -in blast_output

## Disclaimer

- I'm not a computer engineer or some related professional, I'm just write this script to study python and to automatize some bioinformatics tasks. So fell free to commit changes that makes the code more efficient or more clean.
- This script will continue to be developed to englobe others functions, such as filter by subject.