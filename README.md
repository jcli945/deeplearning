basenji_motif.py

# Basenji motif state annotation using HMM

This repository contains a Python script for annotating motif-level regulatory states
based on Basenji-predicted scores using a two-stage Hidden Markov Model (HMM) approach.

The script is designed for downstream analysis of chromatin accessibility–related
sequence effects and was developed for cross-species regulatory element comparison.

---

## Overview

The script performs the following steps:

1. Load base-resolution Basenji prediction scores stored in HDF5 format
2. Compute mean sequence effect scores across nucleotide substitutions
3. Standardize (z-score) the mean scores
4. Apply a Gaussian HMM to estimate posterior state probabilities
5. Use these probabilities as emissions in a Multinomial HMM
6. Decode the final state sequence using the Viterbi algorithm
7. Output per-base regulatory states as a CSV file

---

## Requirements

The script was tested with Python 3 and requires the following packages:

- numpy
- pandas
- h5py
- scikit-learn
- scipy
- hmmlearn
- matplotlib
- seaborn

You can install the dependencies using:

```bash
pip install numpy pandas h5py scikit-learn scipy hmmlearn matplotlib seaborn
