#!/usr/bin/env python3
import sys, argparse
import pysam

def parse_varid(v):
    chrom, pos = v.split(":")[:2]
    return chrom, int(pos)

def fetch_rec(vcf, chrom, pos):
    for rec in vcf.fetch(chrom, pos-1, pos):
        if rec.pos == pos:
            return rec
    return None

def get_gt(rec, sample):
    if (rec is None) or (sample not in rec.samples):
        return None
    return rec.samples[sample].get("GT")  # e.g. (0,1), (1,1), (None,None)

def is_missing(gt):
    return (gt is None) or any(a is None for a in gt)

def is_homozygote(gt, mode="alt"):
    """mode: 'alt' -> hom-nonref; 'ref' -> hom-ref; 'any' -> any homozygote (non-missing)"""
    if gt is None or any(a is None for a in gt):
        return False
    a, b = gt
    if a != b:
        return False
    if mode == "alt":
        return a > 0  # includes 1/1, 2/2, ...
    elif mode == "ref":
        return a == 0
    else:  # any
        return True

def load_pairs_or_vars(path):
    items = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            toks = line.strip().split()
            if len(toks) == 2:
                sample, varid = toks
                items.append(("pair", sample, varid))
            elif len(toks) == 1:
                items.append(("var", None, toks[0]))
            else:
                raise ValueError(f"Bad line (expect 1 or 2 columns): {line.strip()}")
    return items

def main():
    ap = argparse.ArgumentParser(description="Find cases where input GT is missing but output GT is homozygous.")
    ap.add_argument("-i","--input_vcf")
    ap.add_argument("-o",'--output_vcf')
    ap.add_argument("--focals", help="File with either 'sample chr:pos' or just 'chr:pos' per line")
    ap.add_argument("--hom", choices=["alt","ref","any"], default="alt", help="Homozygote type to flag (default: alt)")
    ap.add_argument("--warn-alleles", action="store_true", help="Warn if REF/ALT differ (debug aid)")
    args = ap.parse_args()

    in_vcf  = pysam.VariantFile(args.input_vcf)
    out_vcf = pysam.VariantFile(args.output_vcf)

    focals = load_pairs_or_vars(args.focals)

    total_checked = 0
    correct = 0
    flagged = 0
    num_other = 0

    # Precompute sample lists
    in_samples  = set(list(in_vcf.header.samples))
    out_samples = set(list(out_vcf.header.samples))

    for kind, sample, varid in focals:
        chrom, pos = parse_varid(varid)
        rec_in  = fetch_rec(in_vcf,  chrom, pos)
        rec_out = fetch_rec(out_vcf, chrom, pos)

        if rec_out is None:
            # Nothing to check if variant absent in output
            continue

        if args.warn_alleles and rec_in is not None:
            if rec_in.ref != rec_out.ref or rec_in.alts != rec_out.alts:
                sys.stderr.write(f"WARNING Allele mismatch at {varid}: "
                                 f"IN={rec_in.ref}>{rec_in.alts} OUT={rec_out.ref}>{rec_out.alts}\n")

        # Decide which samples to scan
        if kind == "pair":
            samples_to_check = [sample]
        else:  # kind == "var": scan all output samples (since we need OUT GTs)
            samples_to_check = list(out_samples)

        for s in samples_to_check:
            total_checked += 1
            gt_in  = get_gt(rec_in,  s) if s in in_samples else None
            gt_out = get_gt(rec_out, s)

            if is_missing(gt_in) and is_homozygote(gt_out, mode=args.hom):
                flagged += 1
                print(f"{s}\t{varid}\tIN=NA\tOUT={gt_out}")
            elif is_homozygote(gt_out, mode=args.hom) and is_homozygote(gt_in, mode=args.hom):
                correct += 1
            else:
                print(f"{s}\t{varid}\tIN=NA\tOUT={gt_out}")
                num_other += 1

    print(f"# Checked: {total_checked}\n# Flagged (IN missing & OUT homozygote [{args.hom}]): {flagged}")
    print(f"# Correct (IN & OUT homozygote [{args.hom}]): {correct}")
    print(f"# Other : {num_other}")


if __name__ == "__main__":
    main()
