#!/usr/bin/env python3
import cairo

class Motif():
    """Contains a list of motifs to be found in sequences"""

    def __init__(self, motif_file: str) -> None:
        '''
        Initialize class Motif
        '''
        self.motif_file = motif_file
        self.motif_list = []

    def read_mf(self) -> None:
        '''
        Reads motif file and adds motifs from it to a list
        '''
        with open(self.motif_file, 'r') as mf:
            for motif in mf:
                motif = motif.strip('\n')
                self.motif_list.append(motif)


class MotifMark():
    """Contains functionality for finding motifs in given fasta sequences and drawing a diagram
    representing introns, exons and different motifs"""


    def __init__(self, fasta_file: str, motif_obj) -> None:
        '''
        Initialize class MotifMark
        '''
        self.fasta_file = fasta_file
        self.canvas_height = 0
        self.canvas_width = 0
        self.seq_list = []
        self.motif_obj = motif_obj

    def process_fasta(self) -> None:
        '''
        Main method that controls the class operations going from reading a fasta file
        to outputting an image containing sequence diagrams and motifs
        '''
        self.parse_fasta()

    def parse_fasta(self) -> None:
        '''
        Parses fasta file containing sequences, generates one line sequences and
        saves then into a list
        '''
        header = ''
        seq = ''
        with open(self.fasta_file, 'r') as ff:
            for line in ff:
                line = line.strip('\n')
                if line.startswith('>'):
                    if seq == '':
                        header = line
                    else:
                        self.seq_list.append((header, seq))
                        seq = ''
                        header = line
                else:
                    seq += line.strip('\n')
            self.seq_list.append((header, seq))


    def create_canvas(self):
        '''
        Creates canvas to hold the digrams of sequences and legends
        '''
        
        surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, self.canvas_width, self.canvas_height)
        ctx = cairo.Context (surface)
        ctx.scale (self.canvas_width, self.canvas_height) # Normalizing the canvas

    def find_motif(self) -> None:
        '''
        Iterates over sequences and motifs finding all mitifs in each sequence
        '''

        for seq in self.seq_list:
            for motif in self.motif_obj.motif_list:
                pass    

if __name__ == "__main__":
    mf_file = 'Fig_1_motifs.txt'
    fasta_f = 'Figure_1.fasta'
    my_motif_ls = Motif(mf_file)
    my_motif_ls.read_mf()
    my_motif_mark = MotifMark(fasta_f, my_motif_ls)
    my_motif_mark.process_fasta()
    print(f'{my_motif_ls.motif_list=}')
    print(f'{my_motif_mark.motif_obj.motif_list=}')
    #print(f'{my_motif_mark.seq_list=}')