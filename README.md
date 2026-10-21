# *Escherichia coli* S11 Genome Analysis & Reproducibility Toolkit

[![DOI: 10.6084/m9.figshare.34741515](https://img.shields.io/badge/DOI-10.6084%2Fm9.figshare.34741515-007EC6)](https://doi.org/10.6084/m9.figshare.34741515) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow)](LICENSE) [![Python: 3.13](https://img.shields.io/badge/Python-3.13-brightgreen)](https://www.python.org/downloads/) [![Galaxy: Europe](https://img.shields.io/badge/Galaxy-Europe-6F42C1)](https://usegalaxy.eu/) [![Checks](https://img.shields.io/github/actions/workflow/status/rk-nishan/ecoli-s11-analysis-scripts/validate.yml?branch=main&label=Checks)](https://github.com/rk-nishan/ecoli-s11-analysis-scripts/actions/workflows/validate.yml)

A reproducible Python toolkit for publication figures, local phylogroup typing, assembly and annotation checks, and Galaxy provenance capture for the *Escherichia coli* S11 project.

- **Publication figures** — regenerate Figures 2 and 3 as 600 dpi PNG/TIFF and vector PDF/SVG.
- **Standalone typing** — run the pinned ClermonTyping classifier and independently check primer evidence.
- **Galaxy provenance** — preserve tool versions, parameters, commands and available execution logs.
- **Result verification** — compare FASTA sequences and check existing QUAST and PGAP/ANI outputs.

**[Overview](#overview) · [Quick start](#quick-start) · [Figures](#reproduce-figures-2-and-3) · [ClermonTyping](#standalone-clermontyping) · [Galaxy](#galaxy-history--provenance) · [Script guide](#script-guide) · [Citation](#citation--data-access)**

---

## Overview

The repository and its [Figshare software archive](https://doi.org/10.6084/m9.figshare.34741515) contain **scripts and supporting documentation only**. Research inputs are supplied separately, and generated figures, results and logs remain in the local output folders.

| Analysis component | Execution environment |
|:---|:---|
| QUAST, fastp, Shovill, ABRicate and the other tools in **Manual Analysis S-11** | **[Galaxy Europe — usegalaxy.eu](https://usegalaxy.eu/)** |
| Figure generation and standalone ClermonTyping | Local Python environment; BLAST+ for typing |
| PGAP/ANI | Separate pinned NCBI workflow on a prepared Linux/Docker runtime |

The Galaxy export scripts retrieve existing history records. They do not execute the Galaxy tools locally.

## Quick start

### 1. Clone the repository

```bash
git clone https://github.com/rk-nishan/ecoli-s11-analysis-scripts.git
cd ecoli-s11-analysis-scripts
```

### 2. Create a Python 3.13 environment

<details>
<summary><strong>Windows · PowerShell</strong></summary>

```powershell
py -3.13 -m venv .venv
.venv\Scripts\Activate.ps1
```

</details>

<details>
<summary><strong>Linux / macOS</strong></summary>

```bash
python3.13 -m venv .venv
source .venv/bin/activate
```

</details>

### 3. Install the tested dependencies

```bash
python -m pip install -r requirements.txt
```

The complete environment is pinned in [requirements.lock.txt](requirements.lock.txt). BLAST+ is installed separately for ClermonTyping; the validated run used **BLAST+ 2.17.0+**. Select the workflow below and provide its external inputs.

---

## Reproduce Figures 2 and 3

**Required inputs:** the updated S11 **CARD and ResFinder ABRicate TSV reports** from Galaxy.

```bash
python scripts/prepare_figure_inputs.py --card /path/to/card_report.tsv --resfinder /path/to/resfinder_report.tsv --outdir inputs
python scripts/fetch_references.py --outdir inputs
python analyses/figures/scripts/reproduce_all.py --input-dir inputs --outdir outputs/figures
```

The input helper copies and validates the reports without modifying the source files. The reference helper retrieves the pinned public classification references separately. The orchestrator rebuilds both mapping tables, draws the figures and checks source integrity, ontology paths, counts, geometry and resolution.

**Local exports:** 600 dpi PNG/TIFF and vector PDF/SVG, using the manuscript colors and no overall title, caption or footer description. Use a new output folder for each run. Generated files are not uploaded automatically.

<details>
<summary><strong>Mapping and interpretation notes</strong></summary>

Figure 2 describes distinct gene-class associations. Figure 3 uses ontology-supported gene-locus/mechanism associations; a hit may belong to multiple mechanisms, so the categories are not exclusive. Coverage and identity filters are explicit in the scripts.

CARD ARO 4.0.2 classifies existing hits and does not replace the Galaxy screening database. A database build date in the screening logs must not be interpreted as a release version. ResFinder class metadata are pinned to commit `65c9ba3`.

The orchestrator sets `S11_FIGURE_DATA_DIR` for its child scripts. Developers calling the individual plotting or validation modules directly must set that variable to the prepared input folder and supply the generated mapping directory where applicable.

</details>

## Standalone ClermonTyping

**Required inputs:** the exact S11 assembly and the pinned upstream primers. Run `fetch_references.py` first to obtain the primers separately.

```bash
python analyses/clermon/scripts/run_local_clermont.py --assembly /path/to/S11_contigs.fasta --primers inputs/clermon/primers.fasta --outdir outputs/clermon
python analyses/clermon/scripts/audit_clermont.py --outdir outputs/clermon
```

Put `blastn` and `makeblastdb` on PATH, or pass `--blast-bin` to the run command. The wrapper checks assembly, classifier and primer hashes. The audit independently checks BLAST alignments and supported primer pairs against the assembly.

<details>
<summary><strong>Pinned implementation and platform notes</strong></summary>

The unmodified classifier is pinned to [ClermonTyping commit 69da4db](https://github.com/A-BN/ClermonTyping/tree/69da4db3cd98ab40528bacc289529323ff1407a0). This is its standalone in silico PCR mode without Mash or R.

The executable uses mismatch allowance 6 and a perfect 3′ seed of 3; its help text instead states 2 and 5. The audit records both settings. Reported genome intervals are not represented as experimentally amplified sequences.

Native Windows BLAST uses a temporary working directory without spaces where needed, and returned files are checksum-verified. Choose a fresh output directory for each run.

</details>

## Galaxy history & provenance

**Required access:** the history owner's Galaxy API key, supplied through the **`GALAXY_API_KEY` environment variable**. Keep credentials outside scripts and version control.

```bash
python scripts/download_galaxy_history.py --inventory-only --outdir outputs/history_inventory
python scripts/download_galaxy_history.py --outdir outputs/history_download
python scripts/export_galaxy_provenance.py --outdir outputs/tool_provenance
```

Downloads include hidden collection members and reject deleted or purged dataset downloads. Provenance captures actual settings, available commands and stdout/stderr, versions, metrics, dependencies and wrapper definitions. Missing logs are recorded as unavailable.

<details>
<summary><strong>Reading the exported records</strong></summary>

These commands write to the selected local output directory. They do not upload data to GitHub or Figshare.

Galaxy commands can contain server-specific paths and are execution evidence rather than portable local commands. Use the recorded settings and corresponding Galaxy tool versions for new upstream analyses. A wrapper's declared dependencies do not establish a complete runtime package inventory.

</details>

## Existing results & PGAP

### Compare assemblies and validate saved outputs

**Required inputs:** the existing FASTA files or the original S11 assembly, Galaxy QUAST report and preserved PGAP result directory.

```bash
python scripts/compare_fasta.py /path/to/original.fasta /path/to/galaxy_contigs.fasta --report outputs/comparison.json
python scripts/validate_results.py --assembly /path/to/S11_contigs.fasta --quast-report /path/to/report.tsv --pgap-dir /path/to/S11_PGAP --report outputs/validation.json
```

The checks calculate assembly statistics, verify saved QUAST metrics, compare annotated genomic sequences, check GFF coordinates and reconcile protein FASTA with GenBank translations. Input sequence files are not modified.

### Prepare a pinned PGAP command plan

Provide the exact assembly, an existing launcher from [NCBI commit fbadd450](https://github.com/ncbi/pgap/blob/fbadd450b385fa80b71ed11ecaa2a70df15f98ba/pgap.py) and the prepared reference directory:

```bash
python scripts/run_pgap.py --assembly /path/to/S11_contigs.fasta --launcher /path/to/pgap.py --reference-dir /path/to/PGAP_INPUT_DIR --outdir outputs/pgap_plan
```

The wrapper pins **PGAP 2026-06-18.build8602** and its container digest. By default it validates inputs and writes a command plan. `--execute` requires the matching references and an existing compatible Linux/Docker runtime.

<details>
<summary><strong>Validation scope and preserved configuration</strong></summary>

The result validator targets the original preserved PGAP output, which is distinct from the later identifier-corrected NCBI submission copy. Equivalent NCBI general-ID notation is normalized bijectively in memory.

The PGAP wrapper embeds the original analysis-use configuration, including the documented placeholder-author and temporary protein-ID validation exemptions. It is not a submission-ready configuration and does not install or download the large PGAP databases. A full PGAP rerun was not performed during software preparation.

</details>

---

## Script guide

| Task | Entry point |
|:---|:---|
| Prepare external figure reports | [prepare_figure_inputs.py](scripts/prepare_figure_inputs.py) |
| Retrieve pinned public references | [fetch_references.py](scripts/fetch_references.py) |
| Rebuild and check both figures | [reproduce_all.py](analyses/figures/scripts/reproduce_all.py) |
| Run standalone ClermonTyping | [run_local_clermont.py](analyses/clermon/scripts/run_local_clermont.py) |
| Independently audit primer evidence | [audit_clermont.py](analyses/clermon/scripts/audit_clermont.py) |
| Export non-deleted Galaxy datasets | [download_galaxy_history.py](scripts/download_galaxy_history.py) |
| Export Galaxy settings and logs | [export_galaxy_provenance.py](scripts/export_galaxy_provenance.py) |
| Compare FASTA files | [compare_fasta.py](scripts/compare_fasta.py) |
| Validate saved QUAST and PGAP/ANI outputs | [validate_results.py](scripts/validate_results.py) |
| Prepare or execute the pinned PGAP command | [run_pgap.py](scripts/run_pgap.py) |

## Reproducibility

- **Pinned dependencies and references:** the tested Python environment, upstream commits and input hashes are recorded in the code and dependency lock.
- **Local functional checks:** figure scripts were rerun with external inputs and compared with the saved exports; typing and its independent primer audit were repeated. Existing QUAST and PGAP/ANI outputs were checked separately.
- **Automated code checks:** [GitHub Actions](https://github.com/rk-nishan/ecoli-s11-analysis-scripts/actions/workflows/validate.yml) checks dependency consistency, Python syntax and command interfaces.

Released NumPy 2.2.0 replaces the original unavailable prerelease pin. Historical server setup, SSH, recovery and duplicate checkpoint helpers are excluded from the maintained code.

## Citation & data access

**Software:** [Figshare DOI · 10.6084/m9.figshare.34741515](https://doi.org/10.6084/m9.figshare.34741515) · [Versioned source release](https://github.com/rk-nishan/ecoli-s11-analysis-scripts/releases/tag/v1.0.0) · [CITATION.cff](CITATION.cff)

**Sequencing data:** [SRA · SRR39239575](https://www.ncbi.nlm.nih.gov/sra/SRR39239575) · [BioProject · PRJNA1478760](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1478760) · [BioSample · SAMN61018556](https://www.ncbi.nlm.nih.gov/biosample/SAMN61018556)

Cite the software DOI for these scripts and the NCBI accessions separately for sequencing data. All research inputs, references, figures, results and logs remain outside the repository and its software archive.

## License

First-party code is released under the [MIT License](LICENSE). The unmodified vendored ClermonTyping classifier retains GPL-3.0 and its original attribution. See [third-party notices](THIRD_PARTY_NOTICES.md) for details.
