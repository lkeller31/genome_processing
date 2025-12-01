#!/usr/bin/env python3
#Goal: filter rna gffs compared to a junctions file to only keep transcripts with introns supported by valid junctions
#Before running this, run gt gff3 -addintrons
import sys
import subprocess
import os

def parse_introns_from_gt(gff_with_introns, intron_bed):
    """
    Extract introns from a GFF that already contains intron features
    and convert to BED format.
    """
    with open(gff_with_introns) as gff, open(intron_bed, "w") as bed:
        for line in gff:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            ftype = fields[2]
            if ftype != "intron":
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            attrs = fields[8]

            # Extract Parent transcript
            parent = "NA"
            for item in attrs.split(";"):
                if item.startswith("Parent="):
                    parent = item.split("=")[1]

            bed.write(f"{chrom}\t{start}\t{end}\t{parent}\n")


def get_valid_transcripts(intersect_bed):
    """
    Collect transcripts whose introns (one or more but not explicitly all) are fully supported by junctions.
    """
    valid = set()
    with open(intersect_bed) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) >= 4:
                transcript = fields[3]
                valid.add(transcript)
    return valid


def filter_gff_by_transcripts(original_gff, valid_transcripts, outfile):
    """
    Output only features belonging to transcripts with all canonical introns.
    """
    with open(original_gff) as inp, open(outfile, "w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            attrs = fields[8]
            keep = False

            # Check transcript/exon/gene features
            for item in attrs.split(";"):
                if item.startswith("ID=") or item.startswith("Parent="):
                    tid = item.split("=")[1]
                    if tid in valid_transcripts:
                        keep = True

            if keep:
                out.write(line)


def main():
    if len(sys.argv) != 4:
        print("Usage: filter_rna_gt.py <rna_gff_with_introns.gff3> <mikado_junctions.bed> <output_prefix>")
        sys.exit(1)

    rna_gff = sys.argv[1]
    junctions_bed = sys.argv[2]
    prefix = sys.argv[3]

    intron_bed = prefix + ".bed"
    intersect_bed = prefix + "_supported.bed"
    filtered_gff = prefix + "_filtered.gff3"

    # 1. Extract introns → BED
    print("Extracting introns...")
    parse_introns_from_gt(rna_gff, intron_bed)

    # 2. Intersect introns with junctions
    print("Finding introns supported by junctions...")
    subprocess.run([
        "bedtools", "intersect", "-wa", "-u",
        "-a", intron_bed,
        "-b", junctions_bed
    ], stdout=open(intersect_bed, "w"), check=True)

    # 3. Get supported transcripts
    print("Collecting transcripts with supported introns...")
    valid_transcripts = get_valid_transcripts(intersect_bed)

    # 4. Filter the original GFF
    print("Filtering original GFF...")
    filter_gff_by_transcripts(rna_gff, valid_transcripts, filtered_gff)

    print("\nDone!")
    print(f"Output GFF: {filtered_gff}")


if __name__ == "__main__":
    main()
