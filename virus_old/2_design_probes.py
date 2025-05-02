#!/usr/bin/env python3


from collections import Counter, defaultdict
import click
import rle
from Bio import SeqIO
from pathlib import Path

PROBE_LHS = 'CCTTGGCACCCGAGAATTCCA'
PROBE_RHS = 'A' * 30
MOD_RHS   = '/5Phos/'
IDT_SCALE = '25nm'
IDT_PURIF = 'STD'

PROBE_LEN = 25
GC_MAX = 0.73
GC_MIN = 0.44
HOMOP_MAX = 3
ALLOWED_NUCS = set("AGTC")

@click.command()
@click.option('--output_dir', required=True, help='Directory to write probe TSVs to')
@click.option('--output_fasta', is_flag=True, help='Write probe hybridization regions as FASTA')
@click.option('--idt', is_flag=True, help='Write probes formatted for IDT ordering')
@click.option('--species_name', default=None, help='Optional species tag for naming output')
@click.argument('target_fasta')
def main(target_fasta, output_dir, output_fasta, idt, species_name):
    target_path = Path(target_fasta)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tag = species_name or target_path.stem
    output_file = output_dir / f"{tag}_probes.tsv"

    cand_probes = defaultdict(list)
    for record in SeqIO.parse(target_path, "fasta"):
        rc_seq = record.seq.reverse_complement().upper()

        for pos, lhs, rhs in gen_probe_pairs(rc_seq):
            if not _check_lhs_seq(lhs):
                continue
            if not all(map(_check_gc, [lhs, rhs])):
                continue
            if not all(map(_check_homopolymer, [lhs, rhs])):
                continue

            cand_probes[record.id].append(ProbePair(lhs, rhs, pos))

    keep_probes = defaultdict(list)
    for record_id, probes in cand_probes.items():
        last_pos = None
        for probe in probes:
            if last_pos and abs(probe.lhs_start - last_pos) < (PROBE_LEN * 2):
                continue
            keep_probes[record_id].append(probe)
            last_pos = probe.lhs_start

    probes = [p for probe_list in keep_probes.values() for p in probe_list]

    if not probes:
        print(f"⚠️  {tag}: No valid probes generated.")
        return

    if output_fasta:
        with open(output_file.with_suffix('.fa'), 'w') as out:
            for record_id, probe_list in keep_probes.items():
                region_tag = record_id.split("|")[0]
                for i, p in enumerate(probe_list, 1):
                    out.write(f">{region_tag}_lhs_{i}\n{p.lhs}\n")
                    out.write(f">{region_tag}_rhs_{i}\n{p.rhs}\n")
    elif idt:
        with open(output_file.with_suffix('.idt.tsv'), 'w') as out:
            out.write("#id\tsequence\tscale\tpurification\n")
            for record_id, probe_list in keep_probes.items():
                region_tag = record_id.split("|")[0]
                for i, p in enumerate(probe_list, 1):
                    out.write(f"{region_tag}_lhs-{i}\t{PROBE_LHS + p.lhs}\t{IDT_SCALE}\t{IDT_PURIF}\n")
                    out.write(f"{region_tag}_rhs-{i}\t{MOD_RHS + p.rhs + PROBE_RHS}\t{IDT_SCALE}\t{IDT_PURIF}\n")
    else:
        with open(output_file, 'w') as out:
            out.write("#id\tpos\thyb_region\tprobe\n")
            for record_id, probe_list in keep_probes.items():
                region_tag = record_id.split("|")[0]
                for i, p in enumerate(probe_list, 1):
                    out.write(f"{region_tag}_lhs-{i}\t{p.lhs_start}\t{p.lhs}\t{PROBE_LHS + p.lhs}\n")
                    out.write(f"{region_tag}_rhs-{i}\t{p.lhs_start}\t{p.rhs}\t{MOD_RHS + p.rhs + PROBE_RHS}\n")

    print(f"✅ {tag}: {len(probes)} probes written to {output_file}")

class ProbePair:
    def __init__(self, lhs, rhs, pos):
        self.lhs = str(lhs)
        self.rhs = str(rhs)
        self.lhs_start = pos

def gen_probe_pairs(seq):
    for i in range(len(seq)):
        subseq = seq[i:i+(PROBE_LEN * 2)]
        if len(subseq) < PROBE_LEN * 2:
            break
        lhs = subseq[:PROBE_LEN]
        rhs = subseq[-PROBE_LEN:]
        yield i, lhs, rhs

def _check_gc(seq):
    ct = Counter(seq)
    if not set(ct).issubset(ALLOWED_NUCS):
        return False
    gc = (ct['G'] + ct['C']) / sum(ct.values())
    return GC_MIN <= gc <= GC_MAX

def _check_homopolymer(seq):
    hps = rle.encode(seq)
    return all(i <= HOMOP_MAX for i in hps[1])

def _check_lhs_seq(seq):
    return seq[-1] == 'T'

if __name__ == '__main__':
    main()
