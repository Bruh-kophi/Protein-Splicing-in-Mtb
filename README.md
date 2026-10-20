# Protein-Splicing-in-Mtb
Evidence of Protein splicing in M.tb
# Detecting intein protein splicing in *Mycobacterium tuberculosis* from proteomics data

Analysis code for the BMedSc (Honours) thesis *"Investigating Protein Splicing in
Mycobacterium tuberculosis Using Proteomics Data"*
(Alphonse Hounkpe, Computational Biology Division, University of Cape Town, 2026).

**Contents:** the scripts used to look for evidence that the three *M. tuberculosis*
inteins (in **RecA**, **SufB** and **DnaB**) are cut out of their host proteins
and that the flanking exteins are joined. The scripts work on peptide
identifications from three proteomic datasets.

---

## In plain words

Some proteins are made as one long chain with an extra piece in the middle, an
**intein**. The pieces either side are the **exteins**. After the chain is made,
the intein cuts itself out and joins the exteins together: **protein
splicing**. The genome only records the long chain, so splicing can only be
seen in the proteins themselves.

Proteins are measured by mass spectrometry after being cut into **peptides**
by the enzyme trypsin, which cuts after the amino acids K and R. Four kinds
of peptide can reveal splicing:

```
Unspliced chain:  [ N-extein ][ ......... INTEIN ......... ][ C-extein ]
Spliced protein:  [ N-extein ][ C-extein ]
                             ^ junction (the new bond)
```

| Type | Peptide | What it shows |
|---|---|---|
| 1 Boundary-spanning | runs from an extein into the intein | the unspliced chain is present |
| 2 Junction-spanning | crosses the new extein–extein bond | the exteins were joined (direct evidence) |
| 3 Boundary-start | starts at the first C-extein residue | indirect sign of splicing (RecA only) |
| 4 Intein-end | ends at the intein's last residue, where trypsin cannot cut | the intein's end bond was cut |

A fifth, indirect test asks whether the intein region has **fewer peptides than
expected** (the *depletion test*), and checks that test on hundreds of
ordinary proteins (the *pseudo-intein control*).

---

## The three inteins (1-based, H37Rv precursor numbering)

| Host | UniProt | Precursor | Intein | Residue before junction | Junction peptide |
|---|---|---|---|---|---|
| RecA | P9WHJ3 | 790 aa | 252–691 | K251 (trypsin cuts here) | needs a missed cleavage (`NKCSPPFK`) |
| SufB | P9WFP7 | 846 aa | 253–611 | G252 | `TLIIADEGSYVHYVEGCTAPIYK` |
| DnaB | P9WMR3 | 874 aa | 400–815 | G399 | `ESGSLEQDADVVILLHRPDAFDR` |

Because the residue before each junction differs, **each intein allows different
diagnostics**. A strategy built for one does not transfer to the others.

---

## Repository layout

```
intein-splicing-mtb/
├── README.md                  this file
├── scripts/                   all analysis scripts (keep them in ONE folder)
├── jobs/                      SLURM job and setup scripts for the UCT HPC
└── tests/                     tests on SYNTHETIC data (no real data in this repo)
```

All scripts sit in one folder on purpose: `intein_search.py` imports
`tagmatch.py` from its own folder, and `check_setup.py` looks for the other
scripts beside it.

---

## Scripts by stage

| Stage | Script | What it does | Dataset(s) | Thesis |
|---|---|---|---|---|
| 0 setup | `check_setup.py` | checks files, columns and modules before a long de novo run | 2 | 3.9 |
| 1 de novo I/O | `mztab_to_csv.py` | Casanovo mzTab → peptide CSV (no filtering) | 2 | 3.5 |
| 2 targeted | `check_known_inteins.py` | places peptides on RecA/SufB/DnaB; boundary motifs; trypsin at the junction; all four diagnostic types; 30-peptide panel | 1, 2 | 3.3, 3.4 |
| 2 targeted | `junction_scan.py` | any peptide touching a junction, no enzyme rule, no minimum overhang | 2 | 3.5 |
| 3 discovery | `intein_search.py` + `tagmatch.py` | hypothesis-free search for split peptides (unknown inteins), mass-validated, target–decoy FDR | 1, 2 (de novo) | 3.4, 3.5 |
| 3 discovery | `crosscheck_candidates.py` | discovery candidates vs the known boundaries (handles 0-/1-based coordinates) | 2 | 3.5 |
| 4 database search (MS-GF+) | `mzid_intein_peptides.py` | positions of identified peptides on the three proteins | 2 | 3.5 |
| 4 database search (MS-GF+) | `mzid_nontryptic.py` | tryptic vs non-tryptic peptide ends; background rate | 2 | 3.5 |
| 4 database search (MS-GF+) | `jh.py` | junction-crossing PSMs from the spliced-database search, at any q-value | 2 | 3.5 |
| 4 database search (FragPipe) | `fragpipe_intein_peptides.py` | the three above, for FragPipe `psm.tsv` + pepXML; per run and per group | 3 | 3.6 |
| 5 spectrum | `annot.py` | b/y fragment-ion check of one spectrum | 2 | 3.8 |
| 6 depletion | `depl.py` | digest + binomial depletion test, PXD009315 (record of that analysis) | 2 | 3.3, 3.7 |
| 6 depletion | `ctrl.py` | pseudo-intein control, PXD009315 (record of that analysis) | 2 | 3.7 |
| 6 depletion | `depletion_test.py` | general depletion test + control (reads counts / `psm.tsv`; DnaB included; proline rule as a switch) | 3 (and any) | 3.7 |

`depl.py` and `ctrl.py` are kept exactly as run for PXD009315: their
observed counts and file names are written into the code. For new data use
`depletion_test.py`, which is tested to give **identical** digest counts to
`depl.py` (see `tests/`).

---

## How each dataset was analysed

### Dataset 1 — Progenio, SAWC507 / SAWC5527 (protein-group export only)
The original spectra were not available. The analysis used the MaxQuant
export (repaired with `parse_proteingroups.py`, which is not in this
repository yet) and ran `check_known_inteins.py`, the split-match search and the
depletion test with a 135-protein control. The upstream search was fully
tryptic with a minimum peptide length of 7, so it could **not** report type 4
peptides or `CSPPFK`.

### Dataset 2 — PXD009315 (public; Q Exactive DDA)
```bash
# submitters' identifications
python3 mzid_intein_peptides.py 0686Mass_S7.mzid
python3 mzid_nontryptic.py 0686Mass_*.mzid
# repeat search with spliced sequences (MS-GF+ v2024.03.26), then
python3 jh.py results/*.mzid
# spectrum check of the key peptide
python3 annot.py 0686Mass_S7.mzML 34355 EAHEPVYNLDVENPDSYLAYGFAVHN
# depletion test and control (run in the deposit folder)
python3 depl.py
python3 ctrl.py
# de novo (Casanovo 5.2.1) on S7 and S10
python3 mztab_to_csv.py S7.mztab S7_peptides.csv
python3 check_setup.py --fasta UP000001584_83332.fasta --peptides S7_peptides.csv
python3 check_known_inteins.py UP000001584_83332.fasta S7_peptides.csv --out S7_known.txt
python3 intein_search.py --fasta UP000001584_83332.fasta --peptides S7_peptides.csv \
        --out S7_candidates.tsv --unspliced-out S7_unspliced.tsv
python3 junction_scan.py UP000001584_83332.fasta S7_peptides.csv --out S7_junction.txt
python3 crosscheck_candidates.py S7_candidates.tsv
```

### Dataset 3 — H37Rv, rifampicin 0.5× MIC (24 h) vs DMSO; timsTOF Pro 2 diaPASEF
```bash
# once: environment and database (see jobs/setup_fragpipe_env.sh, jobs/build_search_db.sh)
sbatch --cpus-per-task=2 --mem=8G jobs/fp_intein_job.sh dry   # setup check
sbatch jobs/fp_intein_job.sh                                   # search (semi-tryptic, min length 6)
bash jobs/run_intein_analysis.sh                               # intein analysis + depletion test
```
Run 10979 (empty injection) is not in the manifest. Run 10989 (Pool_2,
843 PSMs against 5,435–8,156 in the other runs) is excluded from pooled
counts with `--exclude-run 10989`.

---

## Requirements

| Use | Software (version used) |
|---|---|
| Most scripts | Python ≥ 3.8, standard library only |
| `intein_search.py` (via `tagmatch.py`) | pandas, numpy, biopython |
| Dataset 2 search | MS-GF+ v2024.03.26, OpenJDK 17.0.20.1 |
| Dataset 2 de novo | Casanovo 5.2.1, weights `casanovo_orbitrap_v5-2-0.ckpt` |
| Dataset 3 search | FragPipe 24.0; MSFragger 4.4.1; diaTracer 2.2.1; IonQuant 1.11.20; DIA-NN 1.8.2 beta 8; Philosopher 5.1.3-RC9; FragPipe-SpecLib 0.1.58 with pandas 2.3.3, scikit-learn 1.5.2 and lxml 6.1.3 (newer pandas/scikit-learn break it); Python 3.11.17 |
| Cluster | UCT ICTS HPC, SLURM, account `cbio`, partition `ada` |

---

## Conventions every script follows

- **Coordinates** are 1-based and inclusive, in precursor numbering.
  Exception: `intein_search.py` writes 0-based slice indices (`gap_start=251,
  gap_end=691` for RecA 252–691). `crosscheck_candidates.py` converts them.
- **I/L**: isoleucine and leucine have the same mass, so all matching treats
  them as identical. Reported residues come from the reference sequence.
- **Proline rule**: `check_known_inteins.py` and `depletion_test.py` (by
  default) apply it, while `depl.py`, `ctrl.py` and `mzid_nontryptic.py` do not.
  The effect on expected intein shares is about 2 percentage points. State
  the setting with every P value.
- **Evidence standard**: a junction or intein-end peptide is accepted only if
  its spectrum shows b and y ions bracketing the bond in question. Every
  script's hit is a lead, not a finding.
- **No numerical comparison across datasets.** Only within-dataset patterns
  are compared between datasets.

---

## Tests

Run from the repository root:
```bash
python3 tests/test_fragpipe_intein_peptides.py   # 20 checks on a synthetic fixture
python3 tests/test_depletion_test.py             # digest identical to depl.py (900 cases)
```
The fixtures are **entirely synthetic**: random sequences under the host
accessions, named `SYNTHETIC_*`, never mixed with real data.

---

## Data integrity

- No synthetic data are ever presented as real. Test inputs are labelled
  synthetic in file names and contents.
- `tagmatch.py` is **unmodified** (MD5 `efd6454395cfa710610c54753215f085`),
  so the statement in `intein_search.py` that it is imported unchanged can
  be checked.
- The comments added on 15 October 2026 changed no code. Each commented
  script was checked to have the same Python syntax tree as the version
  that produced the thesis results.
- Raw data are not in this repository. PXD009315 is available from
  ProteomeXchange/MassIVE (MSV000082209). The Progenio and Dataset 3 data
  belong to their producers.

---

## Provenance and licences — check before making the repository public

- `tagmatch.py` and the split-anchoring approach in `intein_search.py` come
  from the *proteomics-pipelines* code base (frameshift peptide detection).
- Choose a licence for the new code (for example MIT) and add a `LICENSE` file.
- Cite FragPipe, MSFragger, diaTracer, DIA-NN, MS-GF+ and Casanovo as listed
  in the thesis references.

## Citation

Hounkpe, A. (2026). *Investigating protein splicing in Mycobacterium
tuberculosis using proteomics data* [BMedSc (Hons) thesis].
University of Cape Town.

---

## Changes on 15 October 2026 (after the Dataset 3 search)

- `depletion_test.py`: `--exclude-run` now matches the run name only (text before the
  scan numbers). Before, a PSM whose scan number contained the excluded run ID could
  be dropped by mistake. Tests still pass.
- `annot.py`: accepts a whole spectrum id containing `=`, e.g. `frame=69390` for
  diaTracer pseudo-spectra (`python3 annot.py run_diatracer.mzML frame=69390 PEPTIDE`).
  The old `scan=` behaviour (Dataset 2) is unchanged.
- `figures/make_figures.py`: draws thesis Figures 1-6. All numbers in it are typed from
  the analysis outputs (`~/intein_results/TB_rif_DIA/intein_analysis/` and the DIA-NN
  `report.tsv`); nothing is simulated.

### Editable figures (15 October 2026)

* `figures/make_figures.py` now saves every figure twice: `.png` (picture) and
  `.svg` (editable). In the SVG the text is real text in Times New Roman, not
  drawn shapes, so any label can be retyped.
* `scripts/add_svg_figures.py` puts the SVG copies into the thesis Word file
  (the PNG stays as a fallback for old versions of Word):
  `python3 add_svg_figures.py thesis.docx thesis_editable.docx figures/`
* To edit a figure in Word (2016 or later, or Microsoft 365): right-click it >
  **Convert to Shape**, then click any label, arrow or bar to change it. The
  SVG files also open in Inkscape (free), Illustrator and PowerPoint.

---

## Update, 19 October 2026

**Figure numbering in the thesis changed.** The write-up now has five figures;
the SufB fragment-ion figure was dropped. The file names here keep their
original numbers, so they no longer match the thesis:

| File | Thesis figure |
|---|---|
| `figures/fig1_splicing_mechanism.*` | Figure 1 |
| `figures/fig2_inteins_to_scale.*` | Figure 2 |
| `figures/fig3_workflow.*` | Figure 3 |
| `figures/fig4_sufb_junction_ions.*` | not used in the current draft |
| `figures/fig5_depletion.*` | **Figure 4** |
| `figures/fig6_rif_vs_dmso.*` | **Figure 5** |

**Both data figures were redrawn on 19 October.**

- Depletion figure: Dataset 2 recomputed with the proline rule applied, so all
  three datasets are now calculated the same way (RecA expected .587 in each,
  SufB .450). Axis and legend reworded to "proportion" and "peptide-based".
- Rifampicin figure: reduced to the SufB junction-to-extein ratio alone. The
  DnaB panel was removed because, once non-detections were correctly excluded
  from the quantitative axis, it held one measured point and a count. Plotting
  three non-detections at zero and drawing a mean through them asserted a
  measurement that was never made. DnaB is reported in the Results text.

**New script.** `scripts/digest_mc.py` recomputes the expected intein
proportion at different missed-cleavage limits, with and without the proline
rule. Written to answer two questions:

1. Dataset 1's MaxQuant search allowed 3 missed cleavages while the simulation
   used 2. At 3, expectations move by less than 0.03 (RecA .587 -> .573,
   SufB .450 -> .475) and no conclusion changes.
2. Dataset 2's depletion recomputed with the proline rule, to match Datasets 1
   and 3: RecA p = .0049, SufB p = .051, DnaB p = .19.

Run: `python3 scripts/digest_mc.py UP000001584_83332.fasta`

**The seven diagnostic positions** used by the positional search in Dataset 1,
derived from the reference sequences rather than assumed: the intein start and
end boundaries of all three proteins (RecA 251/252 and 691/692, SufB 252/253
and 611/612, DnaB 399/400 and 815/816), plus RecA's first C-extein residue
C692. RecA alone contributes a seventh because K251 precedes its junction and
trypsin cuts after lysine; SufB and DnaB are preceded by glycine.
