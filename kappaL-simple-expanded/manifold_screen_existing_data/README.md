# Structure–electronic manifold screening with existing data

This analysis performs no new DFT, BTE, phonon, or transport calculation.

## Question

Can known thermoelectric materials define reusable neighbourhoods in a joint
structure–electronic manifold, and can those neighbourhoods screen unreported
analogues?

## Construction

The complete-case core contains materials with all selected descriptors:

- Structure S1: composition statistics, density, volume/atom, atom count and
  scalar cell geometry.
- Electronic E2: gap, dielectric magnitude/anisotropy, MBJ gap, electron/hole
  effective masses, mass ratios and effective-mass anisotropy proxies.

Each block is robustly scaled and divided by the square root of its feature
count. The balanced joint metric is the Euclidean product metric of the two
equally weighted blocks. A 30-neighbour locally scaled graph is converted into
diffusion coordinates. PF, kL, and zT are excluded from graph construction.
The quantitative tests use the full 12-dimensional diffusion representation;
UMAP is applied only to those frozen coordinates to make the three maps easier
to read and does not determine any score or conclusion.

Known starrydata2 formulas with reported maximum `zT >= 1` are weak seeds after
the manifold is frozen. Their usefulness is tested by leaving out one complete
formula at a time and asking where it ranks by proximity to the remaining seed
formulas in the full 12-dimensional diffusion representation.

Unknown-candidate ranking is the geometric mean of:

1. proximity percentile to a known high-zT family in the joint manifold; and
2. the previously computed cross-validated dual-channel consistency score.

The second term prevents a close chemical analogue with an explicitly poor
electronic or low-kL score from ranking solely through seed proximity.

## Important boundary

The seed labels are reduced-formula matches. They do not resolve polymorph,
doping, microstructure, dimensionality, or measurement conditions. The output
is an analogue shortlist, not a zT prediction and not a list of confirmed good
thermoelectrics.

## Run

```bash
/home/wangchao/miniconda3/envs/te_manifold/bin/python \
  /home/wangchao/work_wc/2D_ZT/kappaL-simple-expanded/manifold_screen_existing_data/run_joint_manifold_screen.py
```

## Outputs

- `figures/joint_structure_electronic_manifold_screen.png`
- `figures/joint_structure_electronic_manifold_screen.pdf`
- `outputs/joint_manifold_points.csv`
- `outputs/manifold_candidate_ranking.csv`
- `outputs/manifold_formula_retrieval.csv`
- `outputs/joint_manifold_summary.json`

## Weight-free strict-AND audit

`run_strict_and_manifold.py` removes the equal-view coefficient entirely. For
every material pair it converts structure and electronic distances into local
neighbour-rank percentiles and defines

```text
r_AND(i,j) = max(r_structure(i,j), r_electronic(i,j)).
```

The worse view therefore controls the pair. A material that is structurally
close but electronically far cannot be rescued by the structure distance. The
main graph uses the 30 smallest worst-view ranks per node; `k=15,30,50` are all
reported as a locality sensitivity analysis. Purple points in the comparison
figure are selected only by the pre-existing dual score, not by manifold
proximity, which removes the visual selection circularity.

Additional outputs:

- `figures/strict_and_joint_manifold.png`
- `figures/strict_and_joint_manifold.pdf`
- `outputs/strict_and_manifold_points.csv`
- `outputs/strict_and_candidate_ranking.csv`
- `outputs/strict_and_sensitivity.csv`
- `outputs/strict_and_summary.json`

## Consensus audit: why purple points and formula stars disagree

`run_consensus_audit.py` separates three claims that were mixed in the first
joint-manifold figure:

1. the existing cross-validated low-kL and fixed-condition PF scores;
2. direct similarity to the same high-zT seed in the structure and electronic
   descriptor views; and
3. the two-dimensional UMAP layout, which is used only for display.

For every material, the audit chooses the high-zT seed that minimizes

```text
max(structure neighbour-rank percentile,
    electronic neighbour-rank percentile).
```

A purple consensus candidate must be in the top 10% of the existing dual
transport score and must also be within the top 10% neighbourhood of the same
seed in both descriptor views.  This is a screening definition, not an
independent validation.  Formula-matched high-zT rows are drawn as hollow stars
because the current StarryData2-to-JARVIS link does not resolve phase, doping,
dimensionality or microstructure.  In particular, the ten complete-case seed
formulas are displayed on 33 JARVIS structures; examples include `C` matched
from a `g=0.45` sample, `SnS2` from a `3L` sample and `Si` from a
`nano-bulk Si(model)` sample.

Run:

```bash
/home/wangchao/miniconda3/envs/te_manifold/bin/python \
  /home/wangchao/work_wc/2D_ZT/kappaL-simple-expanded/manifold_screen_existing_data/run_consensus_audit.py
```

Additional outputs:

- `figures/consensus_structure_electronic_audit.png`
- `figures/consensus_structure_electronic_audit.pdf`
- `outputs/consensus_audit_points.csv`
- `outputs/consensus_candidates.csv`
- `outputs/consensus_audit_summary.json`
