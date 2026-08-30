# Cross-validated descriptor space without new first-principles calculations

This analysis uses only data already present in the workspace. It does not run
DFT, BTE, phonon, or new transport calculations.

## Targets

- Lattice channel: 137 existing starrydata2 experimental lattice-thermal-
  conductivity values near 300 K, formula-matched to JARVIS structures.
- Electronic channel: existing JARVIS n/p power factors at 600 K and
  `1e20 cm^-3`; the larger carrier-channel value is the regression target.
- External check: starrydata2 reduced-formula matches with maximum reported
  `zT >= 1`. These zT values are not used for fitting.

## Descriptor blocks

- `S0`: elemental count, composition entropy, atomic number, mass and
  electronegativity statistics.
- `S1`: `S0` plus density, volume per atom, atom count and cell-shape scalars.
- `S2`: `S1` plus bulk modulus, shear modulus and Poisson ratio.
- `E0`: band gap and dielectric magnitude/anisotropy.
- `E1`: `E0` plus MBJ gap and electron/hole effective masses.
- `E2`: `E1` plus effective-mass spectrum and mass-complexity proxies.

The mass-complexity columns are simple proxies, not an exact calculation of the
Fermi-surface complexity factor. The descriptor choices follow the logic of
the semi-empirical thermoelectric descriptor work
<https://doi.org/10.1039/C4EE03157A>, the electronic fitness function
<https://doi.org/10.1103/PhysRevMaterials.1.065405>, and the effective-mass /
Fermi-surface-complexity study <https://doi.org/10.1038/s41524-017-0013-3>.

## Validation

Both channels use five-fold cross-validation grouped by chemical system.
Electronic scores for all 9,029 materials are out-of-fold predictions. For the
lattice channel, the block/model is selected by grouped cross-validation on 137
labels and then fitted to those labels to score the unlabelled structures. The
137 labelled JIDs are replaced by their out-of-fold predictions before plotting;
therefore the full vertical axis should be read as a CV-calibrated prediction,
not as 9,029 measured kL values.

The zT external check is only formula-level. It does not prove that a JARVIS
polymorph is identical to the doped, nanostructured, or otherwise processed
experimental material.

## Run

```bash
/home/wangchao/miniconda3/envs/te_manifold/bin/python \
  /home/wangchao/work_wc/2D_ZT/kappaL-simple-expanded/descriptor_fit_no_new_dft/fit_descriptor_space.py
```

## Outputs

- `figures/cross_validated_structure_electronic_space.png`
- `figures/cross_validated_structure_electronic_space.pdf`
- `outputs/descriptor_model_comparison.csv`
- `outputs/cross_validated_descriptor_space.csv`
- `outputs/descriptor_space_summary.json`
