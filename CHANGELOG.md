# Changelog

## 1.1.0 (2026-11)

Version used for the revised manuscript.

* Fixed random seed (`RANDOM_SEED = 20261025`), so that repeated runs give the
  same reaction variants.
* The input file can be given as the first command-line argument
  (`python isodesmic_similarity.py my_molecules.txt`); without an argument the
  script reads `molecules.txt` as before.
* 2000 random cost vectors per molecule instead of 500, so that more distinct
  variants are found.
* `molecules.txt`: corrected SMILES of HMX (the ring had three carbon atoms
  instead of four) and NTO (tautomer of the triazolone ring).

## 1.0.0 (2026-09)

First release, accompanying the submitted manuscript.
