import pandas as pd
from Bio import Align
import os
import sys
from tqdm import tqdm
import pysam
from pyfaidx import Fasta

aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 2
aligner.mismatch_score = -3
aligner.open_gap_score = -7
aligner.extend_gap_score = -2

query_spec = 'Angus'
reference = 'Brahman'

df_dict = {}

brahman = Fasta('/Users/callummacphillamy/PhD/Reference_Genomes/common_UOA-Brahman-Angus/Brahman_oriented2_ARS/Bos_indicus_hybrid.UOA_Brahman_1.ARS-UCD1.2.orientation.dna.toplevel.fa')
angus = Fasta('/Users/callummacphillamy/PhD/Reference_Genomes/common_UOA-Brahman-Angus/Angus/Angus.ARS-UCD1.2.orientation.fa')

stats = {}

for chromosome in range(1, 30):
    print(f"Working on chromosome {chromosome} of {query_spec} to {reference}")
    best_alignment = {}
    stats[f"chr{chromosome}"] = {}

    file = f'./{query_spec}_CpG/bams/{query_spec}.CpGs.chr{chromosome}.1kb.bam'
    bam = pysam.AlignmentFile(file, 'rb')
    for record in tqdm(bam.fetch(until_eof=True)):
        if record.qname not in best_alignment.keys():
            best_alignment[record.qname] = record
            continue
        elif record.qname in best_alignment.keys():
            # Get stats of saved record
            mapq = best_alignment[record.qname].mapq
            alignment_score = best_alignment[record.qname].get_tag("AS")

            # Check which one has better MapQ. We want to prioritise MapQ before alignment score.

            # If saved record mapq and alignment score is higher than the current bam record.
            # do nothing.
            if mapq >= record.mapq and alignment_score >= record.get_tag("AS"):
                continue

            # If saved record mapq is same or greater than but the alignment score
            elif mapq == record.mapq and alignment_score < record.get_tag("AS"):
                best_alignment[record.qname] = record
                continue

            # If saved record has greater mapq and lower alignment score
            # do nothing as we want to favour higher mapq.
            elif mapq > record.mapq and alignment_score < record.get_tag("AS"):
                continue

            # if saved record mapq is less than the current bam record. Update the saved record
            # to be the current bam record
            elif mapq < record.mapq:
                best_alignment[record.qname] = record
                continue
    #print(f'Number of unique CpGs in bam:\t{len(best_alignment.keys())}')
    stats[f"chr{chromosome}"]["num_CpGs"] = len(best_alignment.keys())
        #print(record)
    aligned_len_too_short = {}
    alen_above_900 = {}
    unmapped = {}
    secondary = {}
    supplementary = {}

    for k, record in tqdm(best_alignment.items()):
        if record.is_unmapped:
            unmapped[k] = None
            continue
        if record.is_secondary:
            secondary[k] = None
            continue
        if record.is_supplementary:
            supplementary[k] = None
            continue

        if record.is_mapped:
            if record.alen < 900:
                aligned_len_too_short[k] = None
                continue
            elif record.alen >= 900:
                alen_above_900[k] = None
                continue
    #print(f'Num unmapped: {len(unmapped.keys())}')
    #print(f'Num secondary: {len(secondary.keys())}')
    #print(f'Num supplementary: {len(supplementary.keys())}')
    #print(f'Too short alignment: {len(aligned_len_too_short.keys())}')
    #print(f'Alignment length above 900: {len(alen_above_900.keys())}')
    #print('Total',len(unmapped.keys()) + len(secondary.keys()) + len(supplementary.keys()) + len(aligned_len_too_short.keys()) +len(alen_above_900.keys()))

    stats[f"chr{chromosome}"]["num_unmapped"] = len(unmapped.keys())
    stats[f"chr{chromosome}"]["num_secondary"] = len(secondary.keys())
    stats[f"chr{chromosome}"]["num_supplementary"] = len(supplementary.keys())
    stats[f"chr{chromosome}"]["alen_too_short"] = len(aligned_len_too_short.keys())
    stats[f"chr{chromosome}"]["alen_>_900"] = len(alen_above_900.keys())
    stats[f"chr{chromosome}"]["num_CpGs"] = len(best_alignment.keys())

    too_many_sub_alignments = {}
    no_alignments = {}
    shared = {}
    snp_change = {}
    sub_alignment_too_short = {}

    for k, record in tqdm(best_alignment.items()):
        if k not in alen_above_900.keys():
            continue
        ref_seq = brahman[str(chromosome)][record.reference_start:record.reference_end]
        #ref_seq = angus[str(chromosome)][record.reference_start:record.reference_end]
        query_seq = record.query_sequence

        q_seq = str(record.query_sequence[450:552])
        r_seq = str(brahman[str(chromosome)][record.reference_start+450:record.reference_end-448])
        #r_seq = str(angus[str(chromosome)][record.reference_start + 450:record.reference_end - 448])
        alignments = aligner.align(q_seq, r_seq )
        if len(alignments) > 1:
            too_many_sub_alignments[k] = None
            continue
        elif len(alignments) == 0:
            no_alignments[k] = None
            continue
        elif len(alignments) == 1:
            alignment = alignments[0]
            q_a_start, q_a_stop = alignment.aligned[0][0]
            r_a_start, r_a_stop = alignment.aligned[1][0]
            q_a_seq = q_seq[q_a_start:q_a_stop]
            r_a_seq = r_seq[r_a_start:r_a_stop]
            if len(q_a_seq) != 102:
                sub_alignment_too_short[k] = None
                continue

            assert q_a_seq[50:52] == 'CG', f'{alignment}'
            if q_a_seq[50:52] == r_a_seq[50:52] == 'CG':
                shared[k] = f'{chromosome}\t{record.reference_start+450+r_a_start+50}\t{record.reference_start+450+r_a_start+50+2}\n'
            elif q_a_seq[50:52] == 'CG' and r_a_seq[50:52] != q_a_seq:
                snp_change[k] = f'{chromosome}\t{record.reference_start+450+r_a_start+50}\t{record.reference_start+450+r_a_start+50+2}\n'
            #rint(q_a_seq[50:52])
            #print(r_a_seq[50:52])
            #print(k)
            #print(record.reference_start+450+r_a_start+50,record.reference_start+450+r_a_start+50+2)
    #print(f'Shared CpGs: {len(shared.keys())}')
    #print(f'SNP change CpGs: {len(snp_change.keys())}')
    #print(f'No subalignments: {len(no_alignments.keys())}')
    #print(f'Subalignments too short: {len(sub_alignment_too_short.keys())}')
    #print(f'Too many subalignments: {len(too_many_sub_alignments.keys())}')

    stats[f"chr{chromosome}"]["shared_cpgs"] = len(shared.keys())
    stats[f"chr{chromosome}"]["snp_change"] = len(snp_change.keys())
    stats[f"chr{chromosome}"]["no_sub_alignments"] = len(no_alignments.keys())
    stats[f"chr{chromosome}"]["subalignment_too_short"] = len(sub_alignment_too_short.keys())
    stats[f"chr{chromosome}"]["too_many_subalignments"] = len(too_many_sub_alignments.keys())

    for k in shared.keys():
        if k in snp_change.keys():
            print(k)
            break
    bed = open(f'./{query_spec}_CpG/shared_CpGs/{query_spec}.CpGs.chr{chromosome}.to.{reference}_coords.SHARED.bed', 'w')
    for k, v in tqdm(shared.items()):
        chrom, start, stop = v.split()
        bed.write(f'{chrom}\t{start}\t{stop}\t{k}\n')
    bed.close()

    bed = open(f'./{query_spec}_CpG/shared_CpGs/{query_spec}.CpGs.chr{chromosome}.to.{reference}_coords.SNP_change.bed', 'w')
    for k, v in tqdm(snp_change.items()):
        chrom, start, stop = v.split()
        bed.write(f'{chrom}\t{start}\t{stop}\t{k}\n')
    bed.close()

print(f'Saving stats to {os.getcwd()}/{query_spec}.to.{reference}.CpG.stats.csv')
df = pd.DataFrame.from_dict(stats)
df.to_csv(f'./{query_spec}.to.{reference}.CpG.stats.csv')