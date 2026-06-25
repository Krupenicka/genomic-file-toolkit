# genomic-file-toolkit

A small, dependency-free **command-line toolkit for extracting metrics from
common genomics file formats**. A single object-oriented core handles four
formats through one shared interface, so the same commands work across files.

## Supported formats

| Format | GC content | Read quality | Counts | Region slice |
|--------|:----------:|:------------:|:------:|:------------:|
| FASTA  | ✅ | — | ✅ | — |
| FASTQ  | ✅ | ✅ | ✅ | — |
| SAM    | ✅ | ✅ | ✅ (per chromosome) | ✅ |
| VCF    | ✅ | ✅ | ✅ (per chromosome) | ✅ |

## Design

An abstract base class (`BIOFormat`) defines the metric interface; `FastaFile`,
`FastqFile`, `SamFile` and `VcfFile` implement it. Records are streamed with
generators, Phred quality is normalised to a 0–1 scale, and the right reader is
selected automatically from the file extension.

## Usage

```bash
python genomic_analysis.py <function> <input_file> [-o output] [-r chr:start-end]
```

| `<function>` | Description |
|--------------|-------------|
| `GC_content` | Mean GC content per record |
| `Quality` | Mean read quality (normalised 0–1) |
| `AlignmentsCount` | Record / per-chromosome counts |
| `SliceCount` | Number of alignments in `chr:start-end` |
| `SliceFile` | Alignments inside `chr:start-end` |

### Examples

```bash
python genomic_analysis.py GC_content reads.fasta
python genomic_analysis.py Quality reads.fastq -o quality.txt
python genomic_analysis.py SliceCount alignments.sam -r chr1:10000-20000
```

## Requirements

Python 3, standard library only (`argparse`, `re`, `collections`).

> Note: the FASTQ/SAM readers expect UTF‑16LE-encoded inputs, as used in the
> project's example files.

## License

MIT — see [LICENSE](LICENSE).
