#!/usr/bin/env python3
"""
gff_to_fasta.py

Two use-cases:
  A) Extract the ##FASTA section from a GFF3 (no reference needed)
  B) Extract sequences for features (CDS/exon/gene/mRNA/transcript) from a reference FASTA + GFF3/GTF

Examples:
  # A) Just grab the FASTA section:
  python gff_to_fasta.py --gff annotation.gff3 --output out.fa

  # B) Extract CDS per transcript (default) from reference:
  python gff_to_fasta.py --gff annotation.gff3 --reference genome.fa --type CDS --output cds.fa

  # B) Extract exons per transcript:
  python gff_to_fasta.py --gff annotation.gff3 --reference genome.fa --type exon --output exons.fa

  # B) Extract genes (one record per gene):
  python gff_to_fasta.py --gff annotation.gff3 --reference genome.fa --type gene --output genes.fa

Notes:
  - Coordinates in GFF/GTF are 1-based inclusive; slicing below converts properly.
  - For CDS/exon, features belonging to the same parent/transcript are merged (sorted by coords).
  - Attribute parsing supports GFF3 (ID=,Parent=) and GTF (gene_id "x"; transcript_id "y"; ...).
"""

import sys
import argparse
from collections import defaultdict

def read_fasta_to_dict(path):
    seqs = {}
    current = None
    chunks = []
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith('>'):
                if current is not None:
                    seqs[current] = ''.join(chunks).replace('\n', '').replace('\r','')
                current = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            seqs[current] = ''.join(chunks).replace('\n', '').replace('\r','')
    return seqs

def write_fasta_record(outfh, header, seq, width=60):
    outfh.write(f'>{header}\n')
    for i in range(0, len(seq), width):
        outfh.write(seq[i:i+width] + '\n')

def revcomp(seq):
    comp = str.maketrans('ACGTRYMKBDHVNacgtrymkbdhvn',
                         'TGCAYRKMVHDBNtgcayrkmvhdbn')
    return seq.translate(comp)[::-1]

def detect_attr_style(attr_field):
    # Heuristic: GFF3 uses key=value; GTF uses key "value";
    if '=' in attr_field and ';' in attr_field:
        return 'gff3'
    if '"' in attr_field and ';' in attr_field:
        return 'gtf'
    # default guess: gff3
    return 'gff3'

def parse_attrs(attr_field, style):
    d = {}
    if style == 'gff3':
        # e.g., ID=abc;Parent=xyz;Name=foo
        parts = [p for p in attr_field.strip().split(';') if p]
        for p in parts:
            if '=' in p:
                k, v = p.split('=', 1)
                d[k.strip()] = v.strip()
    else:
        # GTF: key "value"; key2 "value2";
        parts = [p.strip() for p in attr_field.strip().split(';') if p.strip()]
        for p in parts:
            if ' ' in p:
                k, v = p.split(' ', 1)
                v = v.strip().strip('"')
                d[k.strip()] = v
    return d

def copy_fasta_section(gff_path, outfh):
    in_fasta = False
    wrote_any = False
    with open(gff_path, 'r', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('##FASTA'):
                in_fasta = True
                continue
            if in_fasta:
                if line.startswith('>'):
                    wrote_any = True
                outfh.write(line)
    if not wrote_any:
        raise RuntimeError("No ##FASTA section found in GFF. Provide --reference to extract feature sequences.")

def choose_defaults_for_type(feature_type, style):
    """
    Return (group_attr, id_attr) sensible defaults depending on feature type and attribute style.
    - group_attr: how to group multiple rows into a single sequence record (mainly for CDS/exon)
    - id_attr:    what becomes the FASTA header id (if available)
    """
    if feature_type in ('CDS', 'exon'):
        group_attr = 'Parent' if style == 'gff3' else 'transcript_id'
        id_attr = group_attr
    elif feature_type in ('mRNA', 'transcript'):
        group_attr = None
        id_attr = 'ID' if style == 'gff3' else 'transcript_id'
    elif feature_type == 'gene':
        group_attr = None
        id_attr = 'ID' if style == 'gff3' else 'gene_id'
    else:
        # generic defaults
        group_attr = 'Parent' if style == 'gff3' else 'transcript_id'
        id_attr = 'ID' if style == 'gff3' else 'transcript_id'
    return group_attr, id_attr

def main():
    ap = argparse.ArgumentParser(description="Convert GFF/GTF to FASTA.")
    ap.add_argument('--gff', required=True, help='GFF3/GTF file')
    ap.add_argument('--reference', help='Reference genome FASTA (required for feature extraction unless GFF has ##FASTA)')
    ap.add_argument('--type', default='CDS', help='Feature type to extract (e.g., CDS, exon, gene, mRNA, transcript). Default: CDS')
    ap.add_argument('--group-attr', help='Override grouping attribute (e.g., Parent or transcript_id). For CDS/exon, records with same group become one sequence.')
    ap.add_argument('--id-attr', help='Override attribute used for FASTA header id (e.g., ID, gene_id, transcript_id).')
    ap.add_argument('--name-attr', default=None, help='Optional attribute to append in header (e.g., Name, gene_name).')
    ap.add_argument('--respect-phase', action='store_true',
                    help='For CDS only: trim leading bases per phase (0/1/2) before concatenation. (Simple heuristic)')
    ap.add_argument('--output', default='-', help='Output FASTA file or "-" for stdout. Default: -')
    args = ap.parse_args()

    outfh = sys.stdout if args.output == '-' else open(args.output, 'w', encoding='utf-8')

    try:
        # If no reference is provided, try copying ##FASTA section (mode A)
        if not args.reference:
            try:
                copy_fasta_section(args.gff, outfh)
                return
            except RuntimeError:
                # Fall through to feature extraction error
                raise RuntimeError("No reference provided and no ##FASTA section found. Provide --reference to extract feature sequences.")

        # Mode B: feature extraction with reference
        seqs = read_fasta_to_dict(args.reference)

        feature_type = args.type
        groups = defaultdict(list)  # key -> list of (seqid, start, end, strand, phase, attrs)

        attr_style = None
        header_defaults_set = False
        group_attr = args.group_attr
        id_attr = args.id_attr

        with open(args.gff, 'r', encoding='utf-8') as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) != 9:
                    continue  # skip malformed
                seqid, source, ftype, start, end, score, strand, phase, attrs = parts
                if ftype != feature_type:
                    continue

                if attr_style is None:
                    attr_style = detect_attr_style(attrs)
                    if not header_defaults_set:
                        ga, ia = choose_defaults_for_type(feature_type, attr_style)
                        if group_attr is None:
                            group_attr = ga
                        if id_attr is None:
                            id_attr = ia
                        header_defaults_set = True

                ad = parse_attrs(attrs, attr_style)

                try:
                    start_i = int(start)
                    end_i = int(end)
                except ValueError:
                    continue

                phase_val = None
                if feature_type == 'CDS' and args.respect-phase:
                    # phase in GFF3 is 0/1/2; in GTF it's typically '.'
                    if phase in ('0', '1', '2'):
                        phase_val = int(phase)

                key = None
                if group_attr:
                    key = ad.get(group_attr)
                    if key is None:
                        # Fall back: try ID
                        key = ad.get('ID') or ad.get('transcript_id') or ad.get('gene_id')
                else:
                    # one record per line (gene/mRNA/transcript)
                    key = ad.get(id_attr) or f"{seqid}:{start}-{end}({strand})"

                groups[key].append((seqid, start_i, end_i, strand, phase_val, ad))

        if not groups:
            raise RuntimeError(f"No features of type '{feature_type}' found.")

        for key, feats in groups.items():
            # Validate seqid existence
            for (seqid, _, _, _, _, _) in feats:
                if seqid not in seqs:
                    raise RuntimeError(f"Sequence ID '{seqid}' from GFF not found in reference FASTA.")

            # Sort by genomic coordinate
            # For minus strand, we still collect in ascending coords and reverse-complement at the end
            feats_sorted = sorted(feats, key=lambda x: (x[0], x[1], x[2]))
            # Strand determination: assume all subfeatures share strand of first
            strand = feats_sorted[0][3]

            pieces = []
            spans = []
            for (seqid, s, e, strand_item, phase_val, ad) in feats_sorted:
                s0 = s - 1  # convert to 0-based
                e0 = e      # end is inclusive in GFF -> exclusive for slicing
                subseq = seqs[seqid][s0:e0]
                if feature_type == 'CDS' and args.respect-phase and phase_val:
                    # Trim leading bases per phase on the + strand pieces; heuristic for - handled after RC.
                    # This is a simplification; rigorous handling requires transcript-level frame tracking.
                    trim = phase_val if strand == '+' else 0
                    if trim:
                        subseq = subseq[trim:]
                pieces.append(subseq)
                spans.append(f"{seqid}:{s}-{e}")

            seq_concat = ''.join(pieces)
            if strand == '-':
                seq_concat = revcomp(seq_concat)

            # Header
            # Prefer id_attr if present on any feature row
            header_id = None
            if id_attr:
                for (_, _, _, _, _, ad) in feats_sorted:
                    if ad.get(id_attr):
                        header_id = ad[id_attr]
                        break
            if not header_id:
                header_id = key

            header = header_id
            if args.name-attr:
                # try grab a name-like attribute for readability
                for (_, _, _, _, _, ad) in feats_sorted:
                    if ad.get(args.name-attr):
                        header += f" {ad[args.name-attr]}"
                        break
            header += f" [{feature_type}] {','.join(spans)}({strand})"

            write_fasta_record(outfh, header, seq_concat)

    finally:
        if outfh is not sys.stdout:
            outfh.close()

if __name__ == '__main__':
    main()
