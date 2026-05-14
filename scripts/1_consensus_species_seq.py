#!/usr/bin/env python3
import subprocess
from collections import Counter
from pathlib import Path
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import click

@click.command()
@click.option('--fasta_dir', required=True, help='Directory containing individual .fa files')
@click.option('--tsv_file', required=True, help='TSV mapping species to accessions')
@click.option('--output_dir', required=True, help='Base output directory')
def main(fasta_dir, tsv_file, output_dir):
    fasta_dir = Path(fasta_dir)
    output_dir = Path(output_dir)
    combined_dir = output_dir / "combined"
    aligned_dir = output_dir / "aligned"
    consensus_dir = output_dir / "consensus"

    for d in [combined_dir, aligned_dir, consensus_dir]:
        d.mkdir(parents=True, exist_ok=True)

    with open(tsv_file) as f:
        for line in f:
            if line.strip() == '':
                continue
            species, accessions_str = line.strip().split('\t')
            accessions = accessions_str.strip().split()
            species_tag = species.replace(' ', '_')

            combined_fasta = combined_dir / f"{species_tag}_combined.fa"
            aligned_fasta = aligned_dir / f"{species_tag}_aligned.fa"
            consensus_fasta = consensus_dir / f"{species_tag}_consensus.fa"

            with open(combined_fasta, 'w') as out:
                for acc in accessions:
                    fasta_file = fasta_dir / f"{acc}.fa"
                    if not fasta_file.exists():
                        print(f"Warning: skipping missing file: {fasta_file}")
                        continue
                    for record in SeqIO.parse(fasta_file, 'fasta'):
                        record.id = acc
                        record.description = ''
                        SeqIO.write(record, out, 'fasta')

            result = subprocess.run(['mafft', '--auto', combined_fasta], capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Error: MAFFT failed for {species_tag}\n{result.stderr}")
                continue
            with open(aligned_fasta, 'w') as out:
                out.write(result.stdout)

            try:
                alignment = AlignIO.read(aligned_fasta, 'fasta')
                consensus = ''
                for i in range(alignment.get_alignment_length()):
                    column = alignment[:, i]
                    counts = Counter(c for c in column if c not in '-N')
                    consensus += counts.most_common(1)[0][0] if counts else 'N'
                record = SeqRecord(Seq(consensus), id=f"{species_tag}_consensus", description="")
                SeqIO.write(record, consensus_fasta, "fasta")
            except Exception as e:
                print(f"Error: failed to generate consensus for {species_tag}: {e}")
                continue

if __name__ == '__main__':
    main()