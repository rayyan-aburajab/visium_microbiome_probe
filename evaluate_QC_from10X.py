#!/usr/bin/env python3

import click
import pandas as pd
import primer3

# Constants
PRIMER3_ARGS = {
    "mv_conc": 10,
    "dv_conc": 20,
    "dntp_conc": 10,
    "dna_conc": 5000,
}

PREFERRED_JUNCTIONS = {"TA", "TT", "TC", "TG"}
NON_BAD_JUNCTIONS = {"TA", "TT", "TC", "TG", "AA", "AG", "CA", "CT", "GA", "GG", "GT", "TG"}

@click.command()
@click.option('--input_csv', required=True, help='CSV file with probe_id, hyb_region_lhs, hyb_region_rhs')
@click.option('--output_csv', required=True, help='Output CSV with QC annotations')
@click.option('--min_tm', default=65, show_default=True)
@click.option('--max_tm', default=87, show_default=True)
@click.option('--homopolymer_threshold', default=8, show_default=True)
def main(input_csv, output_csv, min_tm, max_tm, homopolymer_threshold):
    df = pd.read_csv(input_csv)

    required_cols = {'probe_id', 'hyb_region_lhs', 'hyb_region_rhs'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Input CSV must contain columns: {required_cols}")

    results = []
    for i, row in df.iterrows():
        probe_id = row['probe_id']
        lhs = row['hyb_region_lhs']
        rhs = row['hyb_region_rhs']
        lhs_id = f"{probe_id}_lhs"
        rhs_id = f"{probe_id}_rhs"

        bait = lhs + rhs
        ligation = lhs[-1] + rhs[0]

        # Tm checks
        lhs_tm = primer3.calc_tm(lhs, **PRIMER3_ARGS)
        rhs_tm = primer3.calc_tm(rhs, **PRIMER3_ARGS)
        tm_pass = min_tm <= lhs_tm <= max_tm and min_tm <= rhs_tm <= max_tm

        # Structure checks
        hp_lhs = primer3.calc_hairpin(lhs, **PRIMER3_ARGS).tm
        hp_rhs = primer3.calc_hairpin(rhs, **PRIMER3_ARGS).tm
        homodimer_lhs = primer3.calc_homodimer(lhs, **PRIMER3_ARGS).tm
        homodimer_rhs = primer3.calc_homodimer(rhs, **PRIMER3_ARGS).tm
        heterodimer_tm = primer3.calc_heterodimer(lhs, rhs, **PRIMER3_ARGS).tm

        # QC flags
        errors = []
        if "N" in bait:
            errors.append("hard_masked_base")
        if has_homopolymer(bait, homopolymer_threshold):
            errors.append("homopolymer")
        if "GGGGG" in bait.upper():
            errors.append("5G_repeat")
        if not tm_pass:
            errors.append("bad_tm")
        if ligation not in NON_BAD_JUNCTIONS:
            errors.append("bad_junction")
        if hp_lhs > min_tm + 5 or hp_rhs > min_tm + 5:
            errors.append("hairpin_risk")
        if homodimer_lhs > (min_tm - 13) or homodimer_rhs > (min_tm - 13):
            errors.append("homodimer_risk")
        if heterodimer_tm > (min_tm - 20):
            errors.append("heterodimer_risk")

        is_preferred = ligation in PREFERRED_JUNCTIONS
        is_valid = ligation in NON_BAD_JUNCTIONS
        score = (1 if is_preferred else 0.5 if is_valid else 0.2) * (1 / (len(errors) + 1))

        results.append([
            probe_id, lhs_id, lhs, rhs_id, rhs, ligation,
            round(lhs_tm, 1), round(rhs_tm, 1), round(heterodimer_tm, 1),
            ";".join(errors), is_preferred, is_valid, round(score, 3)
        ])

    out_df = pd.DataFrame(results, columns=[
        "probe_id", "lhs_id", "lhs_seq", "rhs_id", "rhs_seq", "ligation_junction",
        "lhs_tm", "rhs_tm", "heterodimer_tm",
        "flags", "is_preferred_junction", "is_valid_junction", "score"
    ])
    out_df.to_csv(output_csv, index=False)
    print(f"✅ Final QC complete: {len(out_df)} probe pairs written to {output_csv}")


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

if __name__ == "__main__":
    main()
