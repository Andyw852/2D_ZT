# Expanded simple-parameter dual space

This folder is independent of `kappaL-refactored/` and does not overwrite its
code or results.

It reproduces the visual logic of
`kappaL-refactored/mp_kappaL/figures/te_reference_dual_space_intersection.png`
with a larger JARVIS-DFT cohort and a deliberately small set of interpretable
scalar parameters.

## Data and cohort

- Source: local `jdft_3d-8-18-2021.json` snapshot, 55,723 records.
- Included: positive-gap entries with finite dielectric tensor and JARVIS n/p
  transport scalars.
- JARVIS transport condition: 600 K, carrier concentration `1e20 cm^-3`;
  conductivity uses the constant-relaxation-time approximation.
- The script records the resulting complete-case size in
  `outputs/expanded_simple_space_summary.json`.

## Simple parameter blocks

Structure/chemistry uses 13 scalar parameters: number of elements,
composition entropy, weighted atomic-number/mass/electronegativity statistics,
density, volume per atom, number of atoms, cell anisotropy, and angle
distortion.

Electronic/transport uses seven scalar parameters: band gap, dielectric
magnitude and anisotropy, n/p Seebeck coefficient, and n/p conductivity.
Power factors are retained in the output for diagnosis but are not used in the
distance calculation, avoiding double-counting `S` and conductivity.

The purple set means that a material lies in both top-5% reference
neighbourhoods. It is not an experimental-zT label.

## Run

```bash
/home/wangchao/miniconda3/envs/te_manifold/bin/python \
  /home/wangchao/work_wc/2D_ZT/kappaL-simple-expanded/run_expanded_dual_space.py
```

## Outputs

- `figures/te_reference_dual_space_intersection_expanded.png`
- `figures/te_reference_dual_space_intersection_expanded.pdf`
- `outputs/expanded_simple_space_membership.csv`
- `outputs/expanded_simple_space_features.csv`
- `outputs/expanded_simple_space_summary.json`
