#!/usr/bin/env python3
from pathlib import Path
import click
import pandas as pd
import primer3

LHS = 'CCTTGGCACCCGAGAATTCCA'
RHS = 'A' * 30
MOD_RHS = '/5Phos/'

PRIMER3_ARGS = {
    "mv_conc": 10,
    "dv_conc": 20,
    "dntp_conc": 10,
    "dna_conc": 5000,
}

PREFERRED_JUNCTIONS = {"TA", "TT", "TC", "TG"}
NON_BAD_JUNCTIONS = {"TA", "TT", "TC", "TG", "AA", "AG", "CA", "CT", "GA", "GG", "GT", "TG"}

@click.command()
@click.option('--selection', required=True, help='TSV with probe selection and metadata')
@click.option('--probe_dir', required=True, help='Directory containing probe TSV files')
@click.option('--output_dir', required=True, help='Directory to write outputs')
@click.option('--min_tm', default=65, show_default=True)
@click.option('--max_tm', default=87, show_default=True)
@click.option('--homopolymer_threshold', default=8, show_default=True)
@click.option('--near_dup_threshold', default=5, show_default=True)
def main(selection, probe_dir, output_dir, min_tm, max_tm, homopolymer_threshold, near_dup_threshold):
    probe_dir = Path(probe_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    selection_df = pd.read_csv(selection, sep='\t')

    all_probes = {}
    for file in sorted(probe_dir.glob("*_probes.tsv")):
        df = pd.read_csv(file, sep='\t')
        for _, row in df.iterrows():
            all_probes[row['probe_id']] = row

    qc_results = []
    ordering_results = []

    for _, sel in selection_df.iterrows():
        original_id = sel['original_probe_id']
        new_id = sel['new_probe_id']
        gene_name = sel['gene_name']
        gene_symbol = sel['gene_symbol']
        composite_group = sel['composite_group']

        if original_id not in all_probes:
            print(f"Warning: {original_id} not found in probe TSVs.")
            continue

        probe = all_probes[original_id]
        left = probe['hyb_region_left']
        right = probe['hyb_region_right']
        bait = left + right
        ligation = left[-1] + right[0]

        left_tm = primer3.calc_tm(left, **PRIMER3_ARGS)
        right_tm = primer3.calc_tm(right, **PRIMER3_ARGS)
        tm_pass = min_tm <= left_tm <= max_tm and min_tm <= right_tm <= max_tm

        hp_left = primer3.calc_hairpin(left, **PRIMER3_ARGS).tm
        hp_right = primer3.calc_hairpin(right, **PRIMER3_ARGS).tm
        homodimer_left = primer3.calc_homodimer(left, **PRIMER3_ARGS).tm
        homodimer_right = primer3.calc_homodimer(right, **PRIMER3_ARGS).tm
        heterodimer_tm = primer3.calc_heterodimer(left, right, **PRIMER3_ARGS).tm

        flags = []
        if "N" in bait:
            flags.append("hard_masked_base")
        if has_homopolymer(bait, homopolymer_threshold):
            flags.append("homopolymer")
        if "GGGGG" in bait.upper():
            flags.append("5G_repeat")
        if not tm_pass:
            flags.append("bad_tm")
        if ligation not in NON_BAD_JUNCTIONS:
            flags.append("bad_junction")
        if hp_left > min_tm + 5 or hp_right > min_tm + 5:
            flags.append("hairpin_risk")
        if homodimer_left > (min_tm - 13) or homodimer_right > (min_tm - 13):
            flags.append("homodimer_risk")
        if heterodimer_tm > (min_tm - 20):
            flags.append("heterodimer_risk")

        is_preferred = ligation in PREFERRED_JUNCTIONS
        is_valid = ligation in NON_BAD_JUNCTIONS
        score = (1 if is_preferred else 0.5 if is_valid else 0.2) * (1 / (len(flags) + 1))

        qc_results.append({
            "probe_id": new_id,
            "original_probe_id": original_id,
            "gene_name": gene_name,
            "gene_symbol": gene_symbol,
            "composite_group": composite_group,
            "hyb_region_left": left,
            "hyb_region_right": right,
            "hyb_region_full": bait,
            "ligation_junction": ligation,
            "left_tm": round(left_tm, 1),
            "right_tm": round(right_tm, 1),
            "heterodimer_tm": round(heterodimer_tm, 1),
            "flags": ";".join(flags),
            "is_preferred_junction": is_preferred,
            "is_valid_junction": is_valid,
            "score": round(score, 3),
        })

        ordering_results.append({
            "probe_id_left": f"{new_id}_left",
            "sequence_left": LHS + left,
            "probe_id_right": f"{new_id}_right",
            "sequence_right": MOD_RHS + right + RHS,
        })

    qc_df = pd.DataFrame(qc_results)

    sequences = qc_df["hyb_region_full"].tolist()
    near_dup_flags = []
    for i, seq in enumerate(sequences):
        is_near_dup = any(
            hamming(seq, sequences[j]) <= near_dup_threshold
            for j in range(len(sequences)) if i != j
        )
        near_dup_flags.append(is_near_dup)
    qc_df["near_duplicate"] = near_dup_flags

    qc_df.to_csv(output_dir / "probe_qc.tsv", sep='\t', index=False)
    pd.DataFrame(ordering_results).to_csv(output_dir / "probe_ordering.tsv", sep='\t', index=False)


def has_homopolymer(seq, threshold=8):
    current = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            current += 1
            if current >= threshold:
                return True
        else:
            current = 1
    return False


def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


if __name__ == '__main__':
    main()