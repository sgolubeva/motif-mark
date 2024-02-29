#!/usr/bin/env python3
import cairo
from itertools import product
import re

class MotifList():
    """Contains a list of motifs to be found in sequences"""

    def __init__(self, motif_file: str) -> None:
        '''
        Initialize class Motif
        '''
        self.motif_file = motif_file
        self.motif_list = []
        self.reference = set()
        with open(self.motif_file, 'r') as mf:
            for motif in mf:
                motif = motif.strip('\n')
                self.motif_list.append(motif)

    def generate_reference(self) -> None:
        '''
        Creates a reference set by adding motifs without y character. If a motif contains y character,
        generates motifs with all possible combinations of t and c instead of y charachters and adds resulting strings 
        to the set
        '''
        rep = ['t', 'c']
        for or_motif in self.motif_list:
            or_motif = or_motif.lower()
            y_pos = [pos for pos, char in enumerate(or_motif) if char == 'y']
            if not y_pos:
                self.reference.add(or_motif)
        
            rep_comb = product(rep, repeat=len(y_pos))
        
            for comb in rep_comb:
                new_motif = list(or_motif)
                for j, repl in zip(y_pos, comb):
                    new_motif[j] = repl
                self.reference.add(''.join(new_motif))

class DNAList():
    """Processes fasta file, creates one line fasta, initializes instanses of DNA class and adds
    them to a list"""

    def __init__(self, fasta_f:str) -> None:
        
        self.fasta_f: str = fasta_f
        self.dna_list: list = []
        self.max_len: int = -1


    def __len__(self):
        """Returns the length of the dna list"""

        return len(self.dna_list)
    
    def __iter__(self):
        """Creates an iterator for the DNAList object to iterate over the dna_list"""

        return iter(self.dna_list)
    
    def parse_fasta(self) -> None:
        '''
        Parses fasta file containing sequences, generates one line sequences and
        saves then into a list
        '''
        header = ''
        seq = ''
        with open(self.fasta_f, 'r') as ff:
            for line in ff:
                line = line.strip('\n')
                if line.startswith('>'):
                    if seq == '':
                        header = line
                    else:
                        self.dna_list.append(DNA(seq, header))
                        #print(f'{header=} {seq=}')
                        seq = ''
                        header = line
                else:
                    seq += line.strip('\n')
            self.dna_list.append(DNA(seq, header))
            #print(f'{header=} {seq=}')
    
    def max_seq_len(self):
        """Iterates over the list of dna objects accessing their length, finds and returns
          the max seq length"""
        
        for dna in self.dna_list:
            if len(dna) > self.max_len:
                self.max_len = len(dna)
        return self.max_len


class DNA():
    """Contains DNA sequence information and draws the DNA backbone on canvas"""

    def __init__(self, dna_seq: str, dna_header: str ) -> None:
        """DNA class constructor. Takes DNA sequence and DNA header"""

        self.dna_seq = dna_seq
        self.dna_header = dna_header

    def __len__(self):
        """Returns lendth of the DNA sequence"""

        return len(self.dna_seq)
    
    def get_dna_seq(self):
        """Returns the DNA sequence"""

        return self.dna_seq
    
    def draw_dna(self, ctx: cairo.Context, startx: int, starty:int):
        """Draws the sequence backbone on canvas"""
        header = self.dna_header.split()
        
        fig_label = header[0][1::]
        
        # draw dna label
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_font_size(15)
        ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.move_to(startx, starty - 50)
        ctx.show_text(f'gene {fig_label}')
        ctx.stroke()
        # draw dna backbone
        ctx.set_line_width(3)
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(startx, starty)
        ctx.line_to(len(self.dna_seq), starty)
        ctx.stroke()


class Motif():
    """Contains a single motif sequence and functionality to draw it on provided canvas"""

    def __init__(self, motif:str, coordx:int, cordy:int, motif_height, color:tuple, ctx: cairo.Context) -> None:
        self.motif = motif
        self.coordx = coordx
        self.cordy = cordy 
        self.motif_height = motif_height
        self.rgb_red, self.rgb_green, self.rgb_blue = color
        self.canvas = ctx

    def draw_motif(self):
        """Draws the given motif on canvas"""

        for i in range(len(self.motif)):
            self.canvas.set_line_width(1)
            self.canvas.set_source_rgb(self.rgb_red, self.rgb_green, self.rgb_blue)
            self.canvas.move_to(self.coordx+i, self.cordy + motif_height)
            self.canvas.line_to(self.coordx+i, self.cordy - motif_height)
            self.canvas.stroke()
        
    # def draw_motif(self):
    #     """Draws motif on canvas"""

    #     self.canvas.set_source_rgb(self.rgb_red, self.rgb_green, self.rgb_blue)
    #     self.canvas.rectangle(self.coord1,self.motif_height,len(self.motif),self.motif_height*2)
    #     self.canvas.fill() 

class FindExon():
    """Finds exons in a given sequence and draws them on the canvas"""

    def __init__(self, seq: str) -> None:
        """Initialize find exone class"""

        self.seq = seq          
        self.result_exons = re.finditer('[A-Z]+', self.seq)

    def draw_exons(self, x_start: int, y_start: int, ctx: cairo.Context, exon_height):
        """Draws exons on a given sequence"""

        for exon in self.result_exons:
            ctx.set_source_rgb(0, 0, 0)
            ctx.canvas.rectangle(x_start + exon.start(), y_start - exon_height, len(exon.group()), exon_height*2)
            ctx.fill()

class MotifMark():
    """Contains functionality for finding motifs in given fasta sequences and drawing a diagram
    representing introns, exons and different motifs"""


    def __init__(self, fasta_file: str, motif_obj: MotifList) -> None:
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



    def longest_seq(self):
        '''
        Finds the longest sequence in the list of sequences generated by parse_fasta
        '''
        self.max_len = -1
        for seq in self.seq_list:
            if len(seq[1]) > self.max_len:
                self.max_len = len(seq[1])


    def create_canvas(self):
        '''
        Creates canvas to hold the digrams of sequences and legends
        '''
        self.canvas_height = len(self.seq_list)*50 + (len(self.motif_obj.motif_list))*20
        self.canvas_width = self.max_len + 50
        surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, self.canvas_width, self.canvas_height)
        ctx = cairo.Context (surface)
        return ctx

    def draw_seq(self):
        '''
        Draws sequences on the canvas representing introns as lines and exons as rectangles
        '''

    def find_motif(self) -> list:
        '''
        Iterates over sequences and motifs finding all mitifs in each sequence
        '''
        found_motifs = []
        for seq in self.seq_list:
            seq = seq[1].lower()
            for motif in self.motif_obj.motif_list:
                for i in range(len(seq)+1-len(motif)):
                    slice = seq[i:i + len(motif)]
                    print(f'{len(motif)} {motif} {len(slice)}')
                    if 'y' in motif:
                        is_match = True
                        for i  in range(len(motif)):
                            if motif[i] != slice[i] and slice[i] not in 'ct':
                                is_match = False
                                break
                        if is_match:
                            found_motifs.append(slice)
                    else:
                        if slice == motif:
                            found_motifs.append(slice)
        return found_motifs           

if __name__ == "__main__":
    mf_file = 'test_motif.txt'
    fasta_f = 'test.fasta'

    # my_motif_ls = MotifList(mf_file)
    # my_motif_ls.generate_reference()
    # #print(len(my_motif_ls.reference))
    # my_motif_mark = MotifMark(fasta_f, my_motif_ls)
    # my_motif_mark.process_fasta()
    # #print(my_motif_mark.find_motif())

    ### test motif class
    width = 200
    height = 200
    motif_height = 10
    surf = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context(surf)
    #set canvas color
    ctx.save()
    ctx.set_source_rgb(255, 255, 255)
    ctx.paint()
    ctx.restore()
    # draw line
    ctx.set_line_width(1)
    ctx.set_source_rgb(0.2, 0.23, 0.9)
    ctx.move_to(25,50)
    ctx.line_to(125,50)
    ctx.stroke()

    test_motif = Motif('aaaaaaaa', 50, 50, motif_height, (0.2, 0.23, 0.9), ctx)
    test_motif.draw_motif()
    surf.write_to_png('test_motif.png')
    ## end of test motif class

    ### test DNAList class
    fasta_f = 'Figure_1.fasta'
    dna_list = DNAList(fasta_f)
    dna_list.parse_fasta()
    #print(dna_list.max_seq_len())
    ### end of test DNAList class

    ### test DNA class specifically draw DNA funct
    dna_seq = DNA('tctgccttttgggtaactctttagtattttagcttctagttcctcctctctgccctgttctgctg', '>CLASP1 chr2:121444593-121445363')
    header = '>CLASP1 chr2:121444593-121445363'
    label = header.split()
    label = label[0][1::]
    
    surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
    context = cairo.Context(surface)
    context.save()
    context.set_source_rgb(255, 255, 255)
    context.paint()
    context.restore()
    dna_seq.draw_dna(context, 25, 25)
    surface.write_to_png('test_dna.png')
    ### end of test DNA class

    ### finally test if both DNAList, DNA, and FindExon are working together
    fasta_f = 'Figure_1.fasta'
    dna_list = DNAList(fasta_f)
    dna_list.parse_fasta()
    max_len = dna_list.max_seq_len()
    x_margin = 40
    y_margin = 150
    legend_y_margin = 30
    legend_square_size = 25
    motif_number = 4
    canv_width = max_len + (x_margin)
    canv_height = (y_margin * (len(dna_list)+2)) + (motif_number*legend_square_size + (legend_y_margin * 3))

    surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, canv_width, canv_height)
    context = cairo.Context(surface)
    context.save()
    context.set_source_rgb(255, 255, 255)
    context.paint()
    context.restore()
    update_by = y_margin
    for dna in dna_list:
        dna.draw_dna(context, x_margin, y_margin)
        y_margin = y_margin + update_by
    
    # legend
    print(f' before legend{y_margin=}')
    context.set_source_rgb(0, 0, 0)
    for i in range(motif_number):
        context.rectangle(x_margin, y_margin, legend_square_size, legend_square_size)
        context.fill()
        y_margin = y_margin + legend_y_margin
    

    # scale
    y_margin += 100
    # first vertical line
    context.set_line_width(3)
    context.set_source_rgb(0, 0, 0)
    context.move_to(x_margin, y_margin)
    context.line_to(x_margin, y_margin + 30)
    context.stroke()
    # 0 on the scale
    context.set_font_size(20)
    context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.move_to(x_margin - 5, y_margin + 50)
    context.show_text('0')
    context.stroke()
    # middle mark
    context.move_to(max_len/2, y_margin)
    context.line_to(max_len/2, y_margin + 30)
    context.stroke()
    # middle mark number
    context.set_font_size(20)
    context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.move_to(max_len/2 - 20, y_margin + 50)
    context.show_text(f'{round(max_len/2)}')
    context.stroke()
    # length of the scale
    context.move_to(x_margin, y_margin)
    context.line_to(max_len, y_margin)
    context.stroke()
    # second vertical line
    context.move_to(max_len, y_margin)
    context.line_to(max_len, y_margin + 30)
    context.stroke()
    # max number on the scale
    context.set_font_size(20)
    context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.move_to(max_len - 20, y_margin + 50)
    context.show_text(f'{max_len}')
    surface.write_to_png('test_sequences.png')

    
    ### end of testing DNAList, DNA, FindExon