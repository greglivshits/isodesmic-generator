# Automated isodesmic reaction generator

Code accompanying:

> G. Livshits, M. Ben Noach, D. Hen, G. Nordon, J. Bernstein and D. T. Major,
> *Benchmarking DFT and Machine Learning Potentials for Heat of Formation and
> Combustion of High-Energy Materials*, Phys. Chem. Chem. Phys., submitted.

## What it does

For a target molecule the script enumerates isodesmic reactions against a curated
library of small, well characterised reference species, and ranks them.

1. The target and every library molecule are parsed into **bond signatures** that
   encode bond order together with the hybridisation and aromaticity of the two
   connected atoms.
2. Conservation of every bond type is written as a system of linear equations and
   solved for integer or fractional stoichiometric coefficients with a linear
   programming solver (`scipy.optimize.linprog`), under random cost-vector
   perturbations so that many distinct valid schemes are found.
3. Where strict bond conservation is impossible for fused ring systems, a
   secondary "ring-breaking" relaxation (Kekule form) is applied.
4. Candidate reactions are ranked by **IsoScore**:

       IsoScore = (10 - N) + SimScore - 0.0001 * C

   where *N* is the number of unique components, *SimScore* the mean Tanimoto
   similarity (Morgan fingerprints, radius 2) between the target and the library
   species, and *C* the total stoichiometric cost.

## Requirements

Python 3.9 or newer, with:

    numpy
    scipy
    rdkit
    pillow

Install with `pip install numpy scipy rdkit pillow`.

## Usage

    python isodesmic_similarity.py

The script reads `molecules.txt` from the working directory. Each line holds a
molecule label and its SMILES string, separated by whitespace:

    RDX O=[N+]([O-])N1CN([N+](=O)[O-])CN([N+](=O)[O-])C1

`molecules.txt` in this repository contains the twenty high-energy materials of
the benchmark set. With no input file present the script runs a single-molecule
demonstration on RDX.

## Output

* `isodesmic_structural_similarity.txt` - the ranked reaction schemes with their
  IsoScore, similarity, MCS score, component count and cost.
* `reaction_images_similarity/` - a drawn scheme for each variant.

The reference library is defined in `REF_LIBRARY_SMILES` at the top of the
script and can be extended.

## License

TO BE CHOSEN BY THE AUTHORS BEFORE DEPOSITION.
