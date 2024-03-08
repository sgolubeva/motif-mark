# motif-mark

The goal of this tool is to visualise a transcript and represent introns, exons and given protein binding motifs. 

The motif-mark-oop.py receives a FASTA file with one or more transcript sequences and a motif file with a list of motifs. The script will use pycairo to draw for each transcript:
1. Introns as straight black lines
2. Exons as black rectangles in correct positions on the sequence
3. Each motif represented by a different color in correct position. Overlapping motifs are represented by taller rectangles. 