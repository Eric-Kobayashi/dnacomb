#!/usr/bin/env python3
"""Run full pilot samples with anchored mate coverage and optional agreement."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time

root = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--index', type=int, required=True)
args = parser.parse_args()
samples = [s for s in json.loads((root / 'samples.json').read_text())
           if s['layout'] == 'pilot200' and not s['sample'].startswith('lib-GI')]
n_samples = len(samples)
assert n_samples == 14
control = args.index == 2 * n_samples
sample = samples[-1] if control else samples[args.index % n_samples]
assert not control or sample['sample'] == 'lib_S5_L001'
policy = 'baseline_control' if control else ('agreement' if args.index >= n_samples else 'coverage')
output = root / 'results_anchored' / policy
output.mkdir(parents=True, exist_ok=True)
prefix = output / sample['sample']
if prefix.with_suffix('.completed.json').exists():
    raise RuntimeError('Refusing to overwrite completed output')
spec = json.loads((root / 'burgold_pilot200.json').read_text())
spec['id'] += '_' + policy
if not control:
    guide = next(r for r in spec['regions'] if r['id'] == 'guide2')
    guide['require_forward_complete'] = True
    guide['require_reverse_complete'] = True
    guide['require_mate_agreement'] = policy == 'agreement'
spec_path = prefix.with_suffix('.design.json')
spec_path.write_text(json.dumps(spec, indent=2) + '\n')
project = root.parents[1] / 'software/dnacomb'
pixi = project.parent / 'pixi/bin/pixi'
command = [str(pixi), 'run', '--manifest-path', str(project / 'pixi.toml'), '--locked',
           'dnacomb', sample['r1'], sample['r2'], '--library-spec', str(spec_path),
           '--output', str(prefix), '--threads', os.environ.get('SLURM_CPUS_PER_TASK', '4'),
           '--no-cache', '--sort', '--skip-variants', '--mode',
           'pattern' if control else 'anchored-pattern', '--pattern-length', '8',
           '--pattern-tolerance', '1', '--distance-metric', 'exact',
           '--library', str(root / 'pilot_pairs.tsv')]
provenance = dict(sample=sample, policy=policy, command=command,
                  job_id=os.environ.get('SLURM_JOB_ID'), started_at=time.time())
prefix.with_suffix('.command.json').write_text(json.dumps(provenance, indent=2))
try:
    subprocess.run(command, cwd=project, check=True)
    with prefix.with_suffix('.summary.tsv').open() as stream:
        summary = list(csv.DictReader(stream, delimiter='\t'))
    total = next(int(r['count']) for r in summary if r['group'] == 'all' and r['metric'] == 'total')
    filtered = next(int(r['count']) for r in summary if r['group'] == 'filtered' and r['metric'] == 'total')
    metrics = dict(total=total, guide1=0, guide2=0, both=0, swap=0, designed_pair=0,
                   accepted_g1_missing_g2=0, accepted_g2_missing_g1=0, counted=0)
    with prefix.with_suffix('.counts.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            count = int(row['count'])
            a, b = row['guide1_distance'] == '0', row['guide2_distance'] == '0'
            metrics['counted'] += count
            for key, yes in [('guide1', a), ('guide2', b), ('both', a and b),
                             ('swap', row['combination_status'] == 'recombination'),
                             ('designed_pair', row['combination_status'] == 'match'),
                             ('accepted_g1_missing_g2', a and not row['guide2']),
                             ('accepted_g2_missing_g1', b and not row['guide1'])]:
                metrics[key] += count * yes
    assert metrics['counted'] + filtered == total
    with (root.parents[1] / 'published_datasets/burgold_2026/ena_runs.tsv').open() as stream:
        expected = next(int(r['read_count']) for r in csv.DictReader(stream, delimiter='\t')
                        if r['run_accession'] == sample['runs'][0])
    assert total == expected, (sample['sample'], total, expected)
    provenance.update(metrics=metrics, finished_at=time.time(), exit_code=0)
    prefix.with_suffix('.completed.json').write_text(json.dumps(provenance, indent=2))
    prefix.with_suffix('.exitcode').write_text('0\n')
    print(json.dumps(provenance), flush=True)
except BaseException:
    prefix.with_suffix('.exitcode').write_text('1\n')
    raise
