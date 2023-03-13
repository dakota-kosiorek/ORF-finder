# Dakota Kosiorek
#   ORF and ORF_finder classes to find the open reading frames in a sequence

# DNA to RNA
def transcribe(dna_seq):
    rna_seq = dna_seq.replace("T", "U")
    return rna_seq

class ORF_finder:
    def __init__(self, rna_seq, seq_name, min_ORF_length):
        self.ORFs = list()
        self.seq_name = seq_name[1:]
        self.min_ORF_length = min_ORF_length
        
        # Dictionary for codons and their assosiated amino acids
        self.aa = {
            # UU (first charcter | second character)
            "UUU": "Phe",
            "UUC": "Phe",
            "UUA": "Leu",
            "UUG": "Leu",
            # CU
            "CUU": "Leu",
            "CUC": "Leu",
            "CUA": "Leu",
            "CUG": "Leu",
            # AU
            "AUU": "Ile",
            "AUC": "Ile",
            "AUA": "Ile",
            "AUG": "Met",   # Start
            # GU
            "GUU": "Val",
            "GUC": "Val",
            "GUA": "Val",
            "GUG": "Val",
            # UC
            "UCU": "Ser",
            "UCC": "Ser",
            "UCA": "Ser",
            "UCG": "Ser",
            # CC
            "CCU": "Pro",
            "CCC": "Pro",
            "CCA": "Pro",
            "CCG": "Pro",
            # AC
            "ACU": "Thr",
            "ACC": "Thr",
            "ACA": "Thr",
            "ACG": "Thr",
            # GC
            "GCU": "Ala",
            "GCC": "Ala",
            "GCA": "Ala",
            "GCG": "Ala",
            # UA
            "UAU": "Tyr",
            "UAC": "Tyr",
            "UAA": "STOP",  # Stop
            "UAG": "STOP",  # Stop
            # CA
            "CAU": "His",
            "CAC": "His",
            "CAA": "Gln",
            "CAG": "Gln",
            # AA
            "AAU": "Asn",
            "AAC": "Asn",
            "AAA": "Lys",
            "AAG": "Lys",
            # GA
            "GAU": "Asp",
            "GAC": "Asp",
            "GAA": "Glu",
            "GAG": "Glu",
            # UG
            "UGU": "Cys",
            "UGC": "Cys",
            "UGA": "STOP",  # Stop
            "UGG": "Trp",
            # CG
            "CGU": "Arg",
            "CGC": "Arg",
            "CGA": "Arg",
            "CGG": "Arg",
            # AG
            "AGU": "Ser",
            "AGC": "Ser",
            "AGA": "Arg",
            "AGG": "Arg",
            # GG
            "GGU": "Gly",
            "GGC": "Gly",
            "GGA": "Gly",
            "GGG": "Gly"
        }
        # Dictionary for the amino acid three letter code to the 1 letter code
        self.aa_single = {
            "Ala": "A",
            "Arg": "R",
            "Asn": "N",
            "Asp": "D",
            "Cys": "C",
            "Gln": "Q",
            "Glu": "E",
            "Gly": "G",
            "His": "H",
            "Ile": "I",
            "Leu": "L",
            "Lys": "K",
            "Met": "M",
            "Phe": "F",
            "Pro": "P",
            "Ser": "S",
            "Thr": "T",
            "Trp": "W",
            "Tyr": "Y",
            "Val": "V"
        }
        
        self.find(rna_seq)
    
    def find(self, rna_seq):
        self.orf_cntr = 1
        # Reverse complement the RNA sequence
        rev_comp_rna_seq = self.rev_comp(rna_seq)
        
        # Frames 1-3 (+ strand)
        self.orf_frame(rna_seq, "+")
        # Frames 1-3 (- strand)
        self.orf_frame(rev_comp_rna_seq, "-")
    
    # Gets the ORFs for a sequence on a strand
    def orf_frame(self, seq, strand):
        for shift in range(3):
            nt_total = shift
            orf_start = 0
            orf_stop = 0
            reading = False
            protein_seq = "M"
            
            # Loop through every codon in the ORF
            for i in range(shift, len(seq) - 2, 3):
                codon = seq[i:i+3]
                
                if len(codon) != 3:
                    break
                
                # If a stop codon was encountered after a start
                if self.aa[codon] == "STOP" and reading == True:
                    orf_stop = nt_total + 3
                    nt_count = orf_stop - orf_start
                    aa_count = int(nt_count / 3) - 1
                    
                    if strand == "-":
                        orf_start = len(seq) - orf_start - 1
                        orf_stop = len(seq) - nt_total - 2
                    
                    if nt_count > self.min_ORF_length:
                        myORF = ORF(label="ORF" + str(self.orf_cntr), 
                                    strand=strand, 
                                    frame=1+shift, 
                                    location=(orf_start + 1, orf_stop), 
                                    length=(nt_count, aa_count),
                                    min_ORF_length=self.min_ORF_length,
                                    protein_seq=protein_seq
                        )
                        # Add ORF to self.ORFs list
                        self.ORFs.append(myORF)
                        self.orf_cntr += 1
                    
                    protein_seq = "M"
                    reading = False
                
                # If a start codon was encountered for the first time
                elif self.aa[codon] == "Met" and reading == False:
                    reading = True
                    orf_start = nt_total
                
                elif reading == True:
                    # Get the single amino acid code from the current codon and add it to the ORF protein sequence
                    protein_seq += self.aa_single[self.aa[codon]]
                
                nt_total += 3

    # Reverse compliments the RNA sequence
    def rev_comp(self, seq):
        rev_seq = seq[::-1]
        rev_comp_seq = ""
        for nt in rev_seq:
            if nt == 'A':
                rev_comp_seq += 'U'
            elif nt == 'U':
                rev_comp_seq += 'A'
            elif nt == 'G':
                rev_comp_seq += 'C'
            elif nt == 'C':
                rev_comp_seq += 'G'
        
        return rev_comp_seq
    
    def display(self):
        for orf in self.ORFs:
            orf.display()

    def save_ORFs(self, file):
        file.write(f"{self.seq_name}\n")
        file.write("-------------------\n")
        for orf in self.ORFs:
            file.write(f"{orf.label}\n")
            file.write(f"\tFrame: \t\t\t\t{orf.frame}\n\tStrand: \t\t\t{orf.strand}\n") 
            file.write(f"\tStart: \t\t\t\t{orf.location[0]}\n\tStop: \t\t\t\t{orf.location[1]}\n")
            file.write(f"\tLength (nt): \t\t{orf.length[0]}\n\tLength (aa): \t\t{orf.length[1]}\n")
            file.write(f"\tProtein Sequence: \t{orf.protein_seq}\n\n") 
    
class ORF:
    def __init__(self, label, strand, frame, location, length, min_ORF_length, protein_seq):
        # Name string
        self.label = label
        # What strand the ORF is on (+/-)
        self.strand = strand
        # What frame the ORF is found in [1-3 (+) and 1-3 (-)]
        self.frame = frame
        # Start and stop location tuple for the ORF (start, stop)
        self.location = location
        # ORF length tuple (nucleotide length, amino acid length)
        self.length = length
        # The minimum length each open reading frame should be
        self.min_ORF_length = min_ORF_length
        # The protein sequence of the ORF
        self.protein_seq = protein_seq
        
    def display(self):
        print(f"{self.label}")
        print(f"\tFrame: \t\t{self.frame}\n\tStrand: \t{self.strand}") 
        print(f"\tStart: \t\t{self.location[0]}\n\tStop: \t\t{self.location[1]}")
        print(f"\tLength (nt): \t{self.length[0]}\n\tLength (aa): \t{self.length[1]}")
        print(f"\tProtein Sequence: {self.protein_seq}\n") 

