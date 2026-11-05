# Paired region extraction

Paired structured modes reverse-complement R2 into the construct orientation,
extract each mate, and merge the observations for each region. The default merge
uses a complete observation when available and chooses the observation with
higher mean quality when both are complete.

`--mode anchored-pattern` extends pattern extraction with known read boundaries.
For a fixed-length library region equal to `forward_start_region`, the start of
R1 can supply its opening boundary. For a region equal to
`reverse_start_region`, the end of reverse-complemented R2 can supply its closing
boundary. An observation is marked complete only when its extracted length is
exactly the configured fixed length and its inner flank was found. Shorter,
longer, and variable-length observations retain their original completeness.
This uses the declared read layout; reads should begin at the specified region
boundary, with adapters removed. Existing `pattern` behavior is preserved.

Library regions accept these optional booleans, each defaulting to `false`:

```json
{
  "id": "guide2",
  "seq_type": "Library",
  "min_length": 20,
  "max_length": 20,
  "max_distance": 0,
  "require_forward_complete": true,
  "require_reverse_complete": true,
  "require_mate_agreement": true
}
```

The first two options require complete extraction from the corresponding mate.
`require_mate_agreement` implies both coverage requirements and additionally
requires identical extracted sequences after R2 orientation normalization. With
`--distance-metric exact`, both mates therefore support the same exact library
sequence. These checks apply to paired alignment, positional and pattern modes.
Single-end and raw-mode inputs with enabled requirements are rejected.

When a requirement fails, that region is omitted while the read pair remains in
the count denominator. An omitted region means insufficient accepted evidence;
it does not establish a biological deletion or a single-guide construct.

For the Burgold pilot, R2 begins at guide 2 and ends within upstream constant
sequence. Its reverse complement contains the upstream anchor and the complete
20-base guide, followed by the read boundary. Ordinary pattern mode calls this
`Partial3Prime` because the downstream scaffold is outside the read.
`anchored-pattern` recognizes the configured endpoint and allows independent
mate evidence to be required. Guide 1 remains a 19-base partial observation
against the full 20-base reference because its leading G was not sequenced.
