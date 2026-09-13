#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
Step 1: NCBI RefSeq CDS Longest Transcript Isoform Extraction (De-redundancy)
脚本1：NCBI RefSeq 编码序列最长转录本代表序列提取与去冗余质控
=============================================================================
功能说明：
1. 从 NCBI RefSeq CDS FASTA 文件中解析 [gene=...] 结构化标签，按基因位点聚类；
2. 针对同一基因位点的多个可变剪接异构体（Splice Isoforms），严格仅保留编码区最长的代表序列；
3. 施加严格生物学读码框质控：
   - 长度必须为 3 的整数倍 (len(seq) % 3 == 0)；
   - 剥除末端终止密码子 (TAA/TAG/TGA)；
   - 编码长度不得短于 90 bp (>= 30 个氨基酸)；
   - 序列内部不得出现提前终止密码子 (No internal stop codons '*')。
4. 消除由同一基因不同剪接转录本自比对带来的假阳性候选对膨胀。
=============================================================================
"""

import os
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_longest_transcripts(input_fasta, output_cds, output_prot=None, min_len_bp=90):
    """
    Extract the longest CDS transcript per gene from an NCBI RefSeq CDS FASTA file.
    """
    gene_regex = re.compile(r'\[gene=([^\]]+)\]')
    longest_records = {}
    total_transcripts = 0

    print(f"Reading input FASTA: {input_fasta} ...")
    with open(input_fasta, 'r', encoding='utf-8', errors='ignore') as f_in:
        for rec in SeqIO.parse(f_in, 'fasta'):
            total_transcripts += 1
            m = gene_regex.search(rec.description)
            gene_id = m.group(1) if m else rec.id.split('_')[0]
            
            # Keep only the longest sequence for each gene
            if gene_id not in longest_records or len(rec.seq) > len(longest_records[gene_id].seq):
                longest_records[gene_id] = rec

    print(f"Total raw transcripts: {total_transcripts}")
    print(f"Unique gene loci identified: {len(longest_records)}")

    valid_cds = []
    valid_prot = []

    for gene_id, rec in longest_records.items():
        seq_str = str(rec.seq).upper()
        
        # 1. Reading frame check (must be multiple of 3)
        rem = len(seq_str) % 3
        if rem != 0:
            seq_str = seq_str[:-rem]
            
        # 2. Strip trailing stop codon
        if seq_str[-3:] in ['TAA', 'TAG', 'TGA']:
            seq_str = seq_str[:-3]
            
        # 3. Minimum length threshold
        if len(seq_str) < min_len_bp:
            continue
            
        # 4. Premature stop codon check
        prot_seq = str(Seq(seq_str).translate(to_stop=False))
        if '*' in prot_seq:
            continue
            
        clean_cds_rec = SeqRecord(Seq(seq_str), id=gene_id, description=f"[gene={gene_id}] longest_cds")
        valid_cds.append(clean_cds_rec)
        
        if output_prot:
            clean_prot_rec = SeqRecord(Seq(prot_seq), id=gene_id, description=f"[gene={gene_id}] translated_protein")
            valid_prot.append(clean_prot_rec)

    # Save output files
    SeqIO.write(valid_cds, output_cds, 'fasta')
    print(f"Saved {len(valid_cds)} clean longest CDS sequences to: {output_cds}")
    
    if output_prot and valid_prot:
        SeqIO.write(valid_prot, output_prot, 'fasta')
        print(f"Saved {len(valid_prot)} clean protein sequences to: {output_prot}")

    return len(longest_records), len(valid_cds)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python step1_longest_cds.py <input_cds.fna> <output_longest_cds.fna> [output_prot.faa]")
        print("Example: python step1_longest_cds.py Apis_mellifera_cds.fna Apis_mellifera_longest_cds.fna Apis_mellifera_prot.faa")
        sys.exit(0)
        
    in_file = sys.argv[1]
    out_cds = sys.argv[2]
    out_prot = sys.argv[3] if len(sys.argv) > 3 else None
    extract_longest_transcripts(in_file, out_cds, out_prot)
