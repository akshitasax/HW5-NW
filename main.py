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

    # Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    sub_file = 'substitution_matrices/BLOSUM62.mat'
    sub_mat = NeedlemanWunsch(sub_file, gap_open=-10, gap_extend=-1)

    species_seqs = [gg_seq, mm_seq, br_seq, tt_seq]
    species_headers = [gg_header, mm_header, br_header, tt_header]

    alignments = []

    for seq, header in zip(species_seqs, species_headers):
        score, hs_alignment, species_alignment = sub_mat.align(hs_seq, seq)
        alignments.append((header, score, hs_alignment, species_alignment))
        print(f'{header} alignment with human BRD2: {score}')

    # Sort by score descending (most similar first)
    alignments_sorted = sorted(alignments, key=lambda x: x[1], reverse=True)
    print('Species similarity to human BRD2 (most to least similar):')
    for (header, score, _, _) in alignments_sorted:
        print(f"{header}: {score}")


        


    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix    

if __name__ == "__main__":
    main()
