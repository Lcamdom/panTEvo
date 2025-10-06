from Bio import SeqIO
import pandas as pd
import sys

def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def find_tsd_in_pair(seq1, seq2, min_len=5, max_mismatches=1):
    for length in range(6, min_len - 1, -1):  # Try 6 first, then 5
        for i in range(len(seq1) - length + 1):
            sub1 = seq1[i:i+length]
            for j in range(len(seq2) - length + 1):
                sub2 = seq2[j:j+length]
                mismatches = hamming_distance(sub1, sub2)
                if mismatches <= max_mismatches:
                    return sub1, mismatches
    return None, None

def extract_flanks_and_find_tsd(bedfile, fastafile):
    genome = SeqIO.to_dict(SeqIO.parse(fastafile, "fasta"))
    records = []

    with open(bedfile) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue

            fields = line.strip().split('\t')
            chrom, start, end, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
            seq = genome[chrom].seq

            # Get flanks
            pair1_left = seq[max(0, start - 6):start]                  # upstream of start
            pair1_right = seq[max(0, end - 6):end]                     # upstream of end

            pair2_left = seq[start:start + 6]                          # downstream of start
            pair2_right = seq[end:end + 6]                             # downstream of end

            # Find TSDs in each pair
            tsd1, mm1 = find_tsd_in_pair(str(pair1_left).upper(), str(pair1_right).upper())
            tsd2, mm2 = find_tsd_in_pair(str(pair2_left).upper(), str(pair2_right).upper())

            # Determine match source
            if tsd1 and tsd2:
                tsd_source = "both"
                tsd_final = tsd1 if len(tsd1) >= len(tsd2) else tsd2
                mismatches_final = mm1 if len(tsd1) >= len(tsd2) else mm2
            elif tsd1:
                tsd_source = "pair_1"
                tsd_final = tsd1
                mismatches_final = mm1
            elif tsd2:
                tsd_source = "pair_2"
                tsd_final = tsd2
                mismatches_final = mm2
            else:
                tsd_source = "none"
                tsd_final = "."
                mismatches_final = "."

            records.append({
                "TE_ID": name,
                "Chr": chrom,
                "Start": start,
                "End": end,
                "Pair1_UpstreamStart": str(pair1_left).upper(),
                "Pair1_UpstreamEnd": str(pair1_right).upper(),
                "Pair2_DownstreamStart": str(pair2_left).upper(),
                "Pair2_DownstreamEnd": str(pair2_right).upper(),
                "TSD_found": "Yes" if tsd_source != "none" else "No",
                "TSD_seq": tsd_final,
                "Mismatches": mismatches_final,
                "TSD_source": tsd_source
            })

    return pd.DataFrame(records)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python tsd_pipeline.py <input.bed> <genome.fa>")
        sys.exit(1)

    bed_file = sys.argv[1]
    genome_file = sys.argv[2]

    df = extract_flanks_and_find_tsd(bed_file, genome_file)
    df.to_csv("tsd_output.tsv", sep="\t", index=False)
    print("Done. Output saved to tsd_output.tsv")

