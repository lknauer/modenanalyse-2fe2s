# Sheet-Mapping Reference

This file documents the correspondence between sheet names in the German
and English versions of `modenanalyse_2fe2s`. Use this table when porting
Origin or Excel post-processing scripts between the two versions.

## Workbook 1: `*_analysis.xlsx`

| German (v1.0.0_de)         | English (v1.0.0_en)       | Content                                               |
|----------------------------|---------------------------|-------------------------------------------------------|
| `Modenanalyse`             | `Mode_analysis`           | Per-mode summary (freq, sym, type, OOP%, INP%, etc.)  |
| `Kern_Scores`              | `Core_scores`             | Cluster-core mode scores (in-plane + out-of-plane)    |
| `Abstaende_Gleichgewicht`  | `Equilibrium_distances`   | Equilibrium Fe-Fe, Fe-S, Fe-N, etc. distances         |
| `SCSD`                     | `SCSD`                    | D2h symmetry decomposition coefficients               |
| `Reorganisationsenergie`   | `Reorganization_energy`   | Per-mode lambda_X(i) and dr_X(i) for all channels     |
| `Reorg_Total`              | `Reorg_total`             | System totals Lambda_X = sum_i lambda_X(i)            |
| `Reorg_pro_Bindung`        | `Reorg_per_bond`          | Per-bond reorganization (Cys-S1, Cys-S2, His-N, etc.) |
| `Modulations_Spektren`     | `Modulation_spectra`      | Frequency-resolved M_X(omega)                         |
| `Lambda_kumulativ`         | `Lambda_cumulative`       | Cumulative Lambda_X(omega)                            |
| `B_Faktoren`               | `B_factors`               | Debye-Waller B-factors per atom                       |
| `Info`                     | `Info`                    | Run metadata                                          |

## Column header mapping (most common)

| German                      | English             | Unit / Notes               |
|-----------------------------|---------------------|----------------------------|
| `Modus#`                    | `Mode#`             | mode index                 |
| `Frequenz`                  | `Frequency`         | cm^-1                      |
| `Red.Masse`                 | `Red.mass`          | AMU                        |
| `Frc.Const.`                | `Frc.const.`        | mDyn/A                     |
| `Sym.`                      | `Sym.`              | irrep label                |
| `Typ`                       | `Type`              | In-plane / Out-of-plane / Torsional |
| `Typ (fein)`                | `Type (fine)`       | 7-level classification     |
| `Praez.`                    | `Prec.`             | precision: high / standard |
| `Kanal`                     | `Channel`           | Fe-Fe / Fe-N / Fe-S / N-H / H-acc |
| `Kern OOP%`                 | `Core OOP%`         | percentage                 |
| `Kern Lok%`                 | `Core loc%`         | percentage                 |
| `Kern-Modus`                | `Core mode`         | mode index                 |
| `Rang`                      | `Rank`              | integer                    |
| `Cluster-ID`                | `Cluster ID`        | integer                    |
| `Abstand`                   | `Distance`          | A                          |

## Classification value mapping

The `Typ (fein)` / `Type (fine)` and `Praez.` / `Prec.` columns contain
classification labels:

| German              | English           |
|---------------------|-------------------|
| `Pur INP`           | `Pure INP`        |
| `Pur OOP`           | `Pure OOP`        |
| `Stark INP`         | `Strong INP`      |
| `Stark OOP`         | `Strong OOP`      |
| `Mehrheitlich INP`  | `Majority INP`    |
| `Mehrheitlich OOP`  | `Majority OOP`    |
| `Gemischt`          | `Mixed`           |
| `hoch`              | `high`            |
| `signifikant`       | `significant`     |
| `trivial`           | `trivial`         |
| `standard`          | `standard`        |

## Output filename mapping

| German (v1.0.0_de)         | English (v1.0.0_en)        |
|----------------------------|----------------------------|
| `*_BEFUND.txt`             | `*_REPORT.txt`             |
| `*_analysis.xlsx`          | `*_analysis.xlsx`          |
| `*_analysis_Embeddings.xlsx` | `*_analysis_Embeddings.xlsx` |
| `*_analysis_NIS.xlsx`      | `*_analysis_NIS.xlsx`      |
| `*_analysis_interp*.xlsx`  | `*_analysis_interp*.xlsx`  |

## Notes

- The English version was prepared for international publication on
  Zenodo. Both versions produce **numerically identical** results — only
  user-facing strings differ.
- The German version includes the full 102-page manual
  (`docs/Handbuch.pdf`) with theory, validation, and full Excel sheet
  reference. The English version includes a 17-page English manual
  (`docs/Manual.pdf`) covering theory, configuration, output reference,
  validation, and troubleshooting. The original German manual is also
  shipped with the English edition as a supplementary reference
  (`docs/Manual_DE.pdf`).
- For Origin workflows: when porting between versions, only the sheet
  names and column headers change; the row order and numerical content
  are byte-for-byte identical for the same input data.
