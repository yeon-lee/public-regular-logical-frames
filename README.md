# Regular logical frames

Construction data and computational proof materials for **Designing Group-Valued
Codes with Full Regular Low-Weight Bases**, by **Jong Yeon Lee**.

This release contains only examples and claims in the paper and its appendices:
40 catalogue records and four specialized constructions. It includes their
construction inputs, logical seeds and witnesses, exact-search sources,
completion evidence, and finite proof certificates.

- [CATALOGUE.md](CATALOGUE.md): parameters and construction IDs.
- [catalogue.json](catalogue.json): coefficient/exponent arrays, canonical seeds,
  and physical witnesses for the 40 indexed examples.
- [CLAIMS.md](CLAIMS.md): proof packages supporting the paper's claims.
- [src/](src/): shared verification and exhaustive-search source, stored once.
- [proofs/](proofs/): final proof records grouped by claim or code family.

## Check the release

Use Python 3.10+ and a C++17 compiler:

```sh
python3 verify.py
```

This runs all 20 finite checks in a temporary workspace, including physical
reconstruction, canonical frames, witnesses, orbit-pool and coloring checks,
exclusion coverage, and small-instance engine tests. It does not repeat the
large exhaustive searches. For only catalogue reconstruction, use `--quick`;
`--code lp-f4-1088-128-22-w10-f22` checks one example.

## Inspect or reproduce a proof

```sh
python3 unpack.py /tmp/regular-frames-proof
```

Choose a new directory outside this release. The command restores the documented
verification layout. Read its `docs/VERIFICATION.md` and each proof package's
README for deep-replay commands. Repeated source and data files are stored once;
physical check matrices are reconstructed and required to match their recorded
hashes. [VERIFYING.md](VERIFYING.md) describes the format and evidence limits.

Discovery campaigns, unreported candidates, companion-only results, manuscript
files, Git history, and compiled executables are not included. The author retains
the full research material separately.

## How to cite

If you use or adapt the code constructions, logical-frame methods, or
distance-search and certification method, please cite the accompanying paper:

> Jong Yeon Lee, *Designing Group-Valued Codes with Full Regular Low-Weight Bases*
> (2026), [arXiv:2610.06820](https://arxiv.org/abs/2610.06820) [quant-ph].

```bibtex
@misc{lee2026designinggroupvaluedcodesregular,
  title={Designing Group-Valued Codes with Full Regular Low-Weight Bases},
  author={Jong Yeon Lee},
  year={2026},
  eprint={2610.06820},
  archivePrefix={arXiv},
  primaryClass={quant-ph},
  url={https://arxiv.org/abs/2610.06820},
}
```

### Efficient distance search and exact certification for LP codes

Please cite the paper when using or extending its method for efficient exhaustive
distance certification of lifted-product (LP) codes, including in independent
implementations. The method combines translation-symmetry reduction,
syndrome-directed branching, and stabilizer-overlap pruning to reduce the search
while preserving the lower-bound guarantee. The paper gives the completeness
arguments and demonstrates exact certification for LP examples with more than
a thousand physical qubits.

Searches on smaller constituent and quotient codes help reject weak candidates
early. Exact quantum distance requires completed full-code exclusions below the
weight of a verified logical witness, covering both CSS sectors or using a
verified sector-exchange symmetry. The method and its relation to prior work are
described in the paper's LP section and Appendix A.

Suggested attribution, when applicable:

> We use the exhaustive distance-certification method of Lee, combining
> translation-symmetry reduction with syndrome-directed branching and
> stabilizer-overlap pruning.

### Software and proof data

If you use the released implementation, constructions, or proof data, please
also cite the repository:

```bibtex
@misc{LeeRegularFramesSoftware,
  author = {Lee, Jong Yeon},
  title  = {regular-logical-frames: Constructions and Computational Proof Data},
  url    = {https://github.com/yeon-lee/public-regular-logical-frames},
  note   = {Software and computational proof data}
}
```

[CITATION.cff](CITATION.cff) provides machine-readable citation metadata.
