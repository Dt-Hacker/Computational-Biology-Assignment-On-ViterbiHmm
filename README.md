# HMM Gene Finder

A Hidden Markov Model (HMM) implementation for identifying exon–intron boundaries in DNA sequences, using the Viterbi algorithm to find the most probable state path.

---

## Overview

This tool models a simplified gene structure with five hidden states and uses the Viterbi algorithm to decode the most likely sequence of states (exon, intron, splice donor, etc.) for a given DNA sequence.

It also supports direct log-probability scoring of hand-crafted state paths, making it useful for comparing candidate splice-site positions.

---

## Model

### Hidden States

| State | Symbol | Description |
|-------|--------|-------------|
| Start | `s` | Silent start state |
| Exon | `E` | Coding exon region |
| Splice Donor | `5` | 5′ splice-donor site (exon→intron boundary) |
| Intron | `I` | Non-coding intron region |
| End | `e` | Silent end state |

### Transitions

```
s  → E         (prob 1.0)
E  → E         (prob 0.9)
E  → 5         (prob 0.1)
5  → I         (prob 1.0)
I  → I         (prob 0.9)
I  → e         (prob 0.1)
```

### Emissions (nucleotide probabilities per state)

| State | A | C | G | T | Notes |
|-------|-----|-----|-----|-----|-------|
| `s` | 0.00 | 0.00 | 0.00 | 0.00 | Silent |
| `E` | 0.25 | 0.25 | 0.25 | 0.25 | Uniform (coding) |
| `5` | 0.05 | 0.00 | 0.95 | 0.00 | GU splice donor |
| `I` | 0.40 | 0.10 | 0.10 | 0.40 | AT-rich intron |
| `e` | 0.00 | 0.00 | 0.00 | 0.00 | Silent |

---

## Requirements

- Python 3.10+
- NumPy

Install dependencies:

```bash
pip install numpy
```

---

## Usage

Run directly with the built-in query sequence:

```bash
python hmm_gene_finder.py
```

The script will:
1. Score several hand-crafted candidate paths (varying splice-donor positions) and print their log-probabilities.
2. Run the Viterbi algorithm on the query sequence and print the globally optimal state path.

### Example Output

```
Log-probabilities for candidate paths:
  5-site at position  6    -43.2891
  5-site at position  8    -41.5104
  5-site at position 12    -39.8732
  ...

Query sequence :   CTTCATGTGAAAGCAGACGTAAGTCA
Best state path:   EEEEEEEEEE5IIIIIIIIIIIIIII
```

---

## Code Structure

| Function | Description |
|----------|-------------|
| `log_prob_of_path(state_path, sequence)` | Computes the log-probability of an explicit state path against the observed sequence |
| `initialise_viterbi(sequence)` | Allocates and seeds the Viterbi value and traceback matrices |
| `fill_viterbi(sequence, values, pointers)` | Fills the DP matrices column by column |
| `traceback(values, pointers)` | Recovers the best state path by walking backwards through the pointer matrix |

---

## Customisation

To run on your own sequence, change the `QUERY` variable:

```python
QUERY = "ATGCGTACGATCG..."
```

To adjust the model (e.g. allow intron→exon transitions for multi-exon genes), update the `TRANSITION` and `EMISSION` matrices directly.
