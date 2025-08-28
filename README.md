# FindBacksplice: a Tool for Locating Circular RNA Backsplice Coordinates

This script generates backsplice coordinates for a backsplice junction (BSJ) sequence, such as one generated for a probe/PCR primer, for use in Bioinformatic applications and pipelines. It requires a .fasta file of BSJ sequences, as well as an XML file of BLAST output ran on the same BSJ sequences. It will also require the .fasta file of the genome which was used to generate the BLAST database.

This was designed for use with nucleotide BLAST command-line implementations, such as BLAST+ (available at https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html).

It is ran with: 
    python blast_terminal_analysis.py -b <blast_xml_file> -i <probe_fasta> -g <genome_fasta>

Generates backsplice coordinates for use with circular RNA pipelines. Requires:
- the XML output of BLAST ran on a .fasta file of backsplice junction sequences/probes 
- the .fasta file of backsplice junction sequences/probes
- the .fasta file used to generate the BLAST database (i.e. a genome .fasta file)

It requires Python (version > 3.0), BioPython, and pandas.

The .xml can be generated with the -outfmt 5 tag for command-line nucleotide blast.

This was developed for the application note FindBacksplice. (there will be a preprint DOI here in time)
