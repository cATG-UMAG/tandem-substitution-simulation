#!/usr/bin/env python3
# Gets reference sequences from IG germline, and check the consensus sequence agains them
import re
from glob import glob

from Bio import SeqIO, pairwise2
from Bio.Align import AlignInfo, MultipleSeqAlignment


def main():
    # load germline sequences
    ref_seqs = {}
    for f in glob("fasta/full_ref/*.fasta"):
        for s in SeqIO.parse(f, "fasta"):
            s.id = s.name = s.description = s.id.split("|")[1].replace("*", "-")
            ref_seqs[s.id] = s

    # check and collect references
    references = []
    for f in glob("fasta/IG*.fasta"):
        name = re.sub("^.*/|[.].*$", "", f)

        ref = ref_seqs[name]
        seqs = [x for x in SeqIO.parse(f, "fasta")]

        # get consensus sequence from file
        m_align = MultipleSeqAlignment(seqs)
        align_summary = AlignInfo.SummaryInfo(m_align)
        consensus = align_summary.dumb_consensus()

        # cut reference seq to match other seqs length
        seq_len = len(consensus)
        ref = ref[:seq_len].upper()

        alignments = pairwise2.align.globalxs(ref.seq, consensus, -5, -2)
        if (
            alignments[0][-3] < seq_len - 1
        ):  # check score... left 1 base for errors with consensus
            print("Hey, check {}".format(f))
            print(pairwise2.format_alignment(*alignments[0]))
        else:
            references.append(ref)

    SeqIO.write(references, "fasta/references.fasta", "fasta-2line")


if __name__ == "__main__":
    main()
