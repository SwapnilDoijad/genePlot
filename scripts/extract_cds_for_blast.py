#!/usr/bin/env python3
import io
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_gff_for_cds(gff_path):
    cds_features = []
    reading_fasta = False
    fasta_lines = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                reading_fasta = True
                continue
            if reading_fasta:
                fasta_lines.append(line)
                continue
            if line.startswith("#"):
                continue
            cols = line.strip().split('\t')
            if len(cols) != 9:
                continue
            if cols[2] != "CDS":
                continue
            seqid = cols[0]
            start = int(cols[3]) - 1  # GFF is 1-based
            end = int(cols[4])
            strand = 1 if cols[6] == "+" else -1
            attrs = {key: value for key, value in (field.split("=") for field in cols[8].split(";") if "=" in field)}
            cds_features.append((seqid, start, end, strand, attrs))
    return cds_features, "".join(fasta_lines)

def extract_cds_sequences(cds_features, fasta_data):
    handle = io.StringIO(fasta_data) 
    contigs = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    cds_records = []
    for seqid, start, end, strand, attrs in cds_features:
        if seqid not in contigs:
            continue
        seq = contigs[seqid].seq[start:end]
        if strand == -1:
            seq = seq.reverse_complement()
        feature_id = attrs.get("ID", "unknown")
        gene = attrs.get("gene", "")
        product = attrs.get("product", "")
        header = f"{feature_id} {gene} {product}".strip()
        cds_records.append(SeqRecord(seq, id=header, description=""))
    return cds_records

def main():
    parser = argparse.ArgumentParser(description="Extract CDS features from GFF with embedded FASTA.")
    parser.add_argument("-i", "--input", required=True, help="Input GFF file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    cds_features, fasta_str = parse_gff_for_cds(args.input)
    cds_records = extract_cds_sequences(cds_features, fasta_str)
    SeqIO.write(cds_records, args.output, "fasta")
    # print(f"Extracted {len(cds_records)} CDS sequences to {args.output}")

if __name__ == "__main__":
    main()
