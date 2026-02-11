# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    species_list = [
        (gg_seq, gg_header, "Gallus gallus"),
        (mm_seq, mm_header, "Mus musculus"),
        (br_seq, br_header, "Balaeniceps rex"),
        (tt_seq, tt_header, "Tursiops truncatus")
    ]

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)

    alignments = []
    for seq, header, name in species_list:
        score = nw.align(hs_seq, seq)[0] # compare to human, only care about getting score
        alignments.append((score, name))
    alignments.sort(reverse=True)

    print("species by similarity to human BRD2, with scores:")
    for score, name in alignments:
        print(f"{name}: {score}")

if __name__ == "__main__":
    main()
