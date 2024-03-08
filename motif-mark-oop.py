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
        with open(self.motif_file, 'r') as mf:
            for motif in mf:
                motif = motif.lower()
                motif = Motif(motif.strip('\n'))
                self.motif_list.append(motif)

    def __iter__(self):
        """
        Creates an iterator for the MotifList object to iterate over the motif_list
        """

        return iter(self.motif_list)
    
    def __len__(self):
        """
        Reurns length of the motif list
        """

        return len(self.motif_list)
    
    def __getitem__(self, index):
        """
        Returns an item from self.motif_list given an index
        """

        return self.motif_list[index]

    def get_motif_list(self):
        """
        Returns the list of motifs
        """

        return self.motif_list
    

class DNAList():
    """
    Processes fasta file, creates one line fasta, initializes instanses of DNA class and adds
    them to a list
    """

    def __init__(self, fasta_f:str) -> None:
        
        self.fasta_f: str = fasta_f
        self.dna_list: list = []
        self.max_len: int = -1


    def __len__(self):
        """
        Returns the length of the dna list
        """

        return len(self.dna_list)
    
    def __iter__(self):
        """
        Creates an iterator for the DNAList object to iterate over the dna_list
        """

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
                        print(f'{header}')
                        print(f'{line}')
                        
                        seq = ''
                        header = line
                else:
                    seq += line.strip('\n')
            self.dna_list.append(DNA(seq, header))
            print({header})
          
    
    def max_seq_len(self):
        """
        Iterates over the list of dna objects accessing their length, finds and returns
        the max seq length
        """
        
        for dna in self.dna_list:
            if len(dna) > self.max_len:
                self.max_len = len(dna)
        return self.max_len


class DNA():
    """Contains DNA sequence information and draws the DNA backbone on canvas"""

    def __init__(self, dna_seq: str, dna_header: str ) -> None:
        """
        DNA class constructor. Takes DNA sequence and DNA header
        """

        self.dna_seq = dna_seq
        self.dna_header = dna_header

    def __len__(self):
        """
        Returns lendth of the DNA sequence
        """

        return len(self.dna_seq)
    
    def get_dna_seq(self):
        """
        Returns the DNA sequence
        """

        return self.dna_seq
    
    def draw_dna(self, ctx: cairo.Context, startx: int, starty:int):
        """
        Draws the sequence backbone on canvas
        """

        header = self.dna_header.split()     
        fig_label = header[0][1::]
        
        # draw dna label
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_font_size(20)
        ctx.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.move_to(startx, starty - 100)
        ctx.show_text(f'gene {fig_label}')
        ctx.stroke()
        # draw dna backbone
        ctx.set_line_width(3)
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(startx, starty)
        ctx.line_to(len(self.dna_seq) + startx, starty)
        ctx.stroke()


class Motif():
    """Contains a single motif sequence and functionality to draw it on provided canvas"""

    def __init__(self, motif:str) -> None:
        self.motif = motif
 
    def draw_motif(self, coordx: int, coordy:int, motif_height:int, x_margin:int, color:tuple, ctx: cairo.Context):
        """
        Draws the given motif on canvas
        """

        rgb_red, rgb_green, rgb_blue = color
        ctx.set_source_rgb(rgb_red, rgb_green, rgb_blue)
        #coordx+=x_margin
        for i in range(len(self.motif)):
            ctx.set_line_width(1)         
            ctx.move_to(coordx + i + x_margin, coordy + motif_height)
            ctx.line_to(coordx + i + x_margin, coordy - motif_height)
            ctx.stroke()

    def get_motif(self):
        """
        returns motif sequence when called
        """

        return self.motif
        

class FindExons():
    """Finds exons in a given sequence and draws them on the canvas"""

    def __init__(self, seq: str) -> None:
        """
        Initialize find exone class
        """

        self.seq = seq          
        self.result_exons = re.finditer('[A-Z]+', self.seq)

    def draw_exons(self, x_start: int, y_start: int, ctx: cairo.Context, exon_height):
        """
        Draws exons on a given sequence
        """

        for exon in self.result_exons:
            ctx.set_source_rgb(0, 0, 0)
            ctx.rectangle(x_start + exon.start(), y_start - exon_height, len(exon.group()), exon_height*2)
            ctx.fill()


class MotifMark():
    """Contains functionality for finding motifs in given fasta sequences and drawing a diagram
    representing introns, exons and different motifs"""


    def __init__(self, DNAList: DNAList, MotifList: MotifList, x_margin: int, y_margin: int) -> None:
        '''
        Initialize class MotifMark
        '''
        self.dna_list = DNAList
        self.dna_list.parse_fasta()
        self.motif_list = MotifList
        self.x_margin = x_margin
        self.y_margin = y_margin
        self.max_len = self.dna_list.max_seq_len()
        self.legend_square_size = 25
        self.legend_margin = 30
        self.motif_height = 25
        self.motif_color_list = [(231/255, 76/255, 60/255),  (41/255, 128/255, 185/255),
                                 (39/255, 174/255, 96/255), (241/255, 196/255, 15/255), (211/255, 84/255, 0/255),
                                (0, 0, 255/255), (96/255, 96/255, 96/255), (0, 255/255, 255/255), (142/255, 68/255, 173/255),]
        self.motif_color_dict = {}
        for k, motif in enumerate(self.motif_list):
            motif = motif.get_motif()
            self.motif_color_dict[motif] = self.motif_color_list[k]
        

    def process_fasta(self) -> None:
        '''
        Main method that controls the class operations going from reading a fasta file
        to outputting an image containing sequence diagrams and motifs
        '''
        
        self.create_canvas()
        self.draw_dna_sequences()

        self.draw_legend()
        self.draw_scale()
        self.surface.write_to_png('lookah_sg_overlapping.png')


    def create_canvas(self):
        '''
        Creates canvas to hold the digrams of sequences and legends
        '''
        
        canvas_width = self.max_len + self.x_margin * 2
        canvas_height = (self.y_margin * (len(self.dna_list) + 2)) + (len(self.motif_list) * 
                                                                    self.legend_square_size + (self.legend_margin * len(self.motif_list)))
        #import ipdb; ipdb.set_trace()
        
        self.surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, canvas_width, canvas_height)
        self.context = cairo.Context(self.surface)
        self.context.save()
        self.context.set_source_rgb(255, 255, 255)
        self.context.paint()
        self.context.restore()
        
    
    def draw_dna_sequences(self):
        """
        Draws dna backbones and exons on the given canvas
        """

        update_by = self.y_margin
        for dna in self.dna_list:
            print(self.x_margin)
            dna.draw_dna(self.context, self.x_margin, self.y_margin)
            exons = FindExons(dna.get_dna_seq())
            exons.draw_exons(self.x_margin, self.y_margin, self.context, self.motif_height)
            self._find_motif(dna.get_dna_seq())
            self.y_margin = self.y_margin + update_by 

    def draw_legend(self):
        """
        Draws a legend on canvas
        """

        # draw exon label
        self.context.set_source_rgb(0, 0, 0)
        self.context.rectangle(self.x_margin, self.y_margin, self.legend_square_size, self.legend_square_size)
        self.context.fill()
        self.context.set_font_size(25)
        self.context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.context.move_to(self.x_margin * 2 + self.legend_square_size, self.y_margin + self.legend_square_size/2)
        self.context.show_text(f'exon')
        self.y_margin = self.y_margin + self.legend_margin

        # draw motif colors and motif sequence text
        for motif in self.motif_color_dict:   
            red, green, blue = self.motif_color_dict[motif]
            self.context.set_source_rgb(red, green, blue)
            self.context.rectangle(self.x_margin, self.y_margin, self.legend_square_size, self.legend_square_size)
            self.context.fill()

            self.context.set_font_size(25)
            self.context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            self.context.move_to(self.x_margin * 2 + self.legend_square_size, self.y_margin + self.legend_square_size/2)
            self.context.show_text(f'{motif}')
            self.y_margin = self.y_margin + self.legend_margin

        self.y_margin = self.y_margin + self.legend_square_size
        #draw example overlapping
        
        for i in range(10):
            self.context.set_source_rgb(231/255, 76/255, 60/255)
            self.context.set_line_width(1)
            self.context.move_to(self.x_margin + i, self.y_margin)
            self.context.line_to(self.x_margin + i, self.y_margin + self.motif_height * 2)
            self.context.stroke()

        for k in range(10):
            self.context.set_source_rgb(39/255, 174/255, 96/255)
            self.context.set_line_width(1)
            self.context.move_to(self.x_margin + k + 5, self.y_margin - 10)
            self.context.line_to(self.x_margin + k + 5, self.y_margin + self.motif_height * 2 + 10)
            self.context.stroke()

        self.context.set_font_size(25)
        self.context.set_source_rgb(0, 0, 0)
        self.context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.context.move_to(self.x_margin * 2 + self.legend_square_size, self.y_margin + self.legend_square_size)
        self.context.show_text(f'overlapping motifs')

    def draw_scale(self):
        """
        Draws a scale on canvas
        """

        self.y_margin += 100
        self.draw_vert_scale_line(self.x_margin, self.y_margin, add=30)
        self.draw_scale_text(self.x_margin, self.y_margin, txt_number='0', subtr=5, add=50)
        self.draw_vert_scale_line(self.max_len/2 + self.x_margin, self.y_margin, add=30)
        self.draw_scale_text(self.max_len/2 + self.x_margin, self.y_margin, txt_number=round(self.max_len/2), subtr=20, add=50)
        self.draw_horiz_scale_line(self.x_margin, self.y_margin, self.max_len + self.x_margin)
        self.draw_vert_scale_line(self.max_len + self.x_margin, self.y_margin, add=30)
        self.draw_scale_text(self.max_len + self.x_margin, self.y_margin, self.max_len, subtr=20, add=50 )

    def draw_vert_scale_line(self, x, y, add=0):
        """
        Draws a line for legend
        """

        self.context.set_line_width(3)
        self.context.set_source_rgb(0, 0, 0)
        self.context.move_to(x, y)
        self.context.line_to(x, y + add)

    def draw_horiz_scale_line(self, x, y, max_x):
        """
        Draws a horizontal line for scale
        """
        self.context.move_to(x, y)
        self.context.line_to(max_x, y)
        self.context.stroke()

    def draw_scale_text(self, pos_x, pos_y, txt_number, subtr=0, add=0):
        """
        Draws scale numbers
        """
        self.context.set_font_size(20)
        self.context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.context.move_to(pos_x - subtr, pos_y + add)
        self.context.show_text(f'{txt_number}')
        self.context.stroke()

   
    def _find_motif(self, dna) -> None:
        '''
        Iterates over sequences and motifs finding all motifs in each sequence
        '''

        amb_ncds = {'w': 'at',
                    's': 'cg',
                    'm': 'ac',
                    'k': 'gt',
                    'r': 'ag',
                    'y': 'ct',
                    'b': 'cgt',
                    'd': 'agt',
                    'h': 'act',
                    'v': 'acg',
                    'n': 'acgt',
                    'u': 't'    
                    }
        
        dna = dna.lower()   
        found_motifs = []
        for motif in self.motif_list:
            motif = motif.get_motif()
            regex = ''
            # for char in motif:
            #     if char in amb_ncds:
            #         regex+=f'[{amb_ncds[char]}]'
            #     else:
            #         regex+=f'[{char}]'
            #result = re.finditer(f'(?={regex})', dna)
            regex = ''.join([f'[{amb_ncds[char]}]' if char in amb_ncds else char for char in motif])
            regex = f'(?={regex})'
            result = re.finditer(regex, dna)
            for item in result:
                found_motifs.append((motif, item.start(), item.end()+len(motif)))
            

        sorted_found_motifs = sorted(found_motifs, key=lambda x: x[1])
        print(f'{sorted_found_motifs}')
        print(f'{dna=}')
        pos_track = []
        for k in range(len(sorted_found_motifs)):

            pos_track = list(filter(lambda x: x >= sorted_found_motifs[k][1], pos_track))
            pos_track.append(sorted_found_motifs[k][2])
            found_motif = Motif(sorted_found_motifs[k][0])
            found_motif.draw_motif(sorted_found_motifs[k][1], self.y_margin, self.motif_height + 10 * (len(pos_track) - 1),
                                   self.x_margin, self.motif_color_dict[sorted_found_motifs[k][0]], self.context)
                

if __name__ == "__main__":
    X_MARGIN = 25
    Y_MARGIN = 200
    fasta_f = 'Figure_1.fasta'
    motif_f = 'Fig_1_motifs.txt'

    dna_list = DNAList(fasta_f)
    motif_list = MotifList(motif_f)
    motif_mark = MotifMark(dna_list, motif_list, X_MARGIN, Y_MARGIN)
    motif_mark.process_fasta()

