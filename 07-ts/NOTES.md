# 07-ts

ωB97X-D3/def2-TZVP, `DefGrid3 VeryTightSCF`, ORCA 6.1.1. Energies relative to trans (`03-scan`, torsion 180°,
−572.789576541 Eh), 1 Eh = 2625.4996 kJ/mol.

## Torsion 90°: closed-shell stability (`torsion90/`)

Geometry: point 90° of `03-scan/torsion_envelope.xyz` (from `torsion_from75`, RKS-optimized). 4 cores each.

- `rks_stab.inp`: RKS + `%scf STABPerform true STABRestartUHFifUnstable true end` (10–12 min)
- `uks_bs.inp`: UKS + `%scf BrokenSym 1,1 end` (triplet first, then the broken-symmetry singlet; 10 min)

```
nohup /usr/local/bin/orca_6_1_1/orca rks_stab.inp > rks_stab.out &
grep -A4 "eigenvalues of the stability matrix" rks_stab.out
grep -E "Expectation value of <S\*\*2>|Dispersion correction|FINAL SINGLE POINT ENERGY" rks_stab.out
grep -E "E\(High-Spin\)|E\(BrokenSym\)" uks_bs.out
```

- RKS is unstable: lowest eigenvalue of the stability matrix −0.0570 Eh; after the restart in UKS the
  solution is stable (lowest eigenvalue +0.0202 Eh).
- Both jobs converge to the same broken-symmetry solution (difference 1·10⁻⁶ Eh).
- `E(High-Spin)` and `E(BrokenSym)` in the BrokenSym block are SCF energies without D3; the D3 correction
  (−0.007792721 Eh, same geometry) is added below.

| Solution | ⟨S²⟩ | E / Eh | E − E(trans) / kJ mol⁻¹ |
|---|---|---|---|
| RKS (`03-scan`) | 0 | −572.704786235 | 222.62 |
| UKS broken-symmetry singlet | 1.071 | −572.731857262 | 151.54 |
| UKS triplet (SCF + D3) | 2.072 | −572.738988499 | 132.82 |

## Inversion TS (`inversion/`)

Start: point 180° of `03-scan/inversion_envelope.xyz` (`inversion_start.xyz`), 8 cores, 42 min.

```
! wB97X-D3 def2-TZVP DefGrid3 TightOpt VeryTightSCF OptTS Freq
%geom Calc_Hess true end
```

```
cp inversion_start.xyz inversion.xyz
nohup /usr/local/bin/orca_6_1_1/orca inversion.inp > inversion.out &
grep -A6 "Geometry convergence" inversion.out | tail -6
grep "imaginary mode" inversion.out | tail -1
grep -E "FINAL SINGLE POINT ENERGY|Zero point energy|Total enthalpy  |Final Gibbs free energy" inversion.out | tail -4
```

- Converged after 2 cycles (gradient and step criteria met; last energy change 9.8·10⁻⁶ Eh above TolE 1·10⁻⁶).
  The TS lies 0.02 kJ/mol above the scan point.
- Exactly one imaginary mode, −421.44 cm⁻¹. Displacing along it bends C2–N0=N1 (largest amplitudes on N0
  and C2) and leads to CNNC 180° on one side and 0° on the other.
- TS geometry: C2–N0=N1 179.4°, N0=N1–C13 116.4°, N=N 1.209 Å.
- The output contains two thermochemistry blocks (initial Hessian, final frequencies); the values below are
  from the last one. The imaginary mode is not included in the thermochemistry.

| | E / Eh | ZPE / Eh | H / Eh | G / Eh |
|---|---|---|---|---|
| TS | −572.725993661 | 0.19052383 | −572.52390013 | −572.57205784 |

Barriers (kJ/mol, gas phase, 298.15 K; minima from `02-freq`):

| | ΔE‡ | ΔE‡ + ZPE | ΔH‡ | ΔG‡ |
|---|---|---|---|---|
| TS − cis | 116.02 | 110.96 | 110.84 | 108.38 |
| TS − trans | 166.94 | 161.39 | 160.82 | 160.81 |

Ea = ΔH‡(TS − cis) + RT = 113.32 kJ/mol (unimolecular gas-phase reaction, 298.15 K).
