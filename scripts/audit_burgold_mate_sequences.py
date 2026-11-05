#!/usr/bin/env python3
"""Independent coordinate audit of the complete pilot plasmid read pairs."""
import collections
import csv
import gzip
import json
from pathlib import Path

root = Path(__file__).resolve().parent
sample = next(s for s in json.loads((root / 'samples.json').read_text()) if s['sample'] == 'lib_S5_L001')
with (root / 'pilot_pairs.tsv').open() as stream:
    catalog = list(csv.DictReader(stream, delimiter='\t'))
guide2_set = {r['guide2'].encode() for r in catalog}
complement = bytes.maketrans(b'ACGTN', b'TGCAN')
counts = collections.Counter()
examples = collections.defaultdict(list)

def records(path):
    with gzip.open(path, 'rb') as stream:
        while name := stream.readline():
            seq = stream.readline().strip()
            assert stream.readline().startswith(b'+')
            qual = stream.readline().strip()
            assert len(seq) == len(qual)
            yield name.strip(), seq, qual

for (id1, r1, q1), (id2, r2, q2) in zip(records(sample['r1']), records(sample['r2']), strict=True):
    # ENA assigns different accession prefixes to these two SINGLE deposits;
    # the original instrument identifier is the second header token.
    assert id1.split()[1].removesuffix(b'/1') == id2.split()[1].removesuffix(b'/2')
    # R1: 19 guide1 + 76 scaffold1 + 6 leader + 71 tRNA.
    forward = r1[172:192]
    reverse = r2[:20].translate(complement)[::-1]
    a, b = forward in guide2_set, reverse in guide2_set
    counts['total'] += 1
    counts['forward_exact'] += a
    counts['reverse_exact'] += b
    counts['both_exact'] += a and b
    counts['same_exact'] += a and b and forward == reverse
    counts['sequences_agree'] += forward == reverse
    if a and b and forward == reverse:
        category = 'both_exact_agree'
    elif a and not b:
        category = 'forward_exact_reverse_error'
    elif b and not a:
        category = 'reverse_exact_forward_error_or_layout'
    elif a and b:
        category = 'both_known_but_disagree'
    else:
        category = 'neither_exact_at_expected_coordinates'
    counts[category] += 1
    if len(examples[category]) < 4:
        examples[category].append(dict(read=id1.decode(), forward_guide2=forward.decode(),
            reverse_guide2=reverse.decode(), forward_exact=a, reverse_exact=b,
            forward_upstream_anchor=r1[164:172].decode(), forward_downstream_anchor=r1[192:200].decode(),
            forward_quality=q1[172:192].decode(), reverse_quality=q2[:20][::-1].decode()))

output = root / 'results_anchored' / 'coordinate_audit.json'
output.write_text(json.dumps(dict(sample=sample, counts=counts, examples=examples,
    definition='R1[172:192] and reverse-complement(R2[:20]); no flank correction or realignment'), indent=2))
print(json.dumps(counts))
