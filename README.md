# E. coli S11 analysis scripts

Scripts for reproducing the S11 manuscript figures, running standalone ClermonTyping, checking existing assembly/annotation outputs and exporting Galaxy history provenance.

**Software DOI:** [10.6084/m9.figshare.34741515](https://doi.org/10.6084/m9.figshare.34741515).

This repository and its Figshare software archive contain **code and supporting documentation only**. Figures, input tables, reference files, assemblies, sequencing reads, annotation files, results and logs are excluded. Provide the required inputs separately. Sequencing data are associated with [SRA SRR39239575](https://www.ncbi.nlm.nih.gov/sra/SRR39239575), [BioProject PRJNA1478760](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1478760) and [BioSample SAMN61018556](https://www.ncbi.nlm.nih.gov/biosample/SAMN61018556).

QUAST, fastp, Shovill, ABRicate and the other tools in **Manual Analysis S-11** were run on **[usegalaxy.eu](https://usegalaxy.eu/)**. The export scripts preserve their existing job parameters and available logs; they do not run those tools locally. Figure generation and standalone ClermonTyping are local steps. The separate PGAP/ANI run retains its own pinned NCBI release.

## Install

Use Python 3.13. On Windows, `py -3` can create the environment before activation.

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

The complete tested dependency set is in `requirements.lock.txt`. BLAST+ is a separate dependency for ClermonTyping; the validated run used BLAST+ 2.17.0+.

## Reproduce Figures 2 and 3

Supply the **CARD and ResFinder ABRicate TSV reports** from the updated S11 Galaxy analysis. The input preparation script copies and validates them without modifying the source files. The reference helper downloads the pinned public CARD ontology, ResFinder class metadata and ClermonTyping primers into a local input folder; these reference files are not included in this repository or Figshare.

```bash
python scripts/prepare_figure_inputs.py --card /path/to/card_report.tsv --resfinder /path/to/resfinder_report.tsv --outdir inputs
python scripts/fetch_references.py --outdir inputs
python analyses/figures/scripts/reproduce_all.py --input-dir inputs --outdir outputs/figures
```

The orchestrator rebuilds the mapping tables, draws both figures and checks their sources, ontology paths, counts, dimensions and resolution. It refuses an existing output folder. The local exports are **600 dpi PNG/TIFF and vector PDF/SVG**, with the manuscript colors and no overall caption, title or footer description. None of the generated files is uploaded automatically.

Figure 2 describes distinct gene-class associations. Figure 3 uses ontology-supported gene-locus/mechanism associations; a hit may belong to multiple mechanisms, so its categories are not exclusive. Coverage and identity filtering is explicit in the scripts. CARD ARO 4.0.2 classifies existing hits and does not replace the Galaxy screening database. The screening logs' database build date must not be interpreted as a release version.

The entry point sets `S11_FIGURE_DATA_DIR` for its child scripts. Developers invoking the individual plotting or validation modules directly must set that environment variable to the prepared input folder and supply the generated mapping directory where applicable. Normal use requires only the orchestrator command above.

## Run standalone ClermonTyping

Prepare the pinned references with `fetch_references.py`, then supply the exact S11 assembly. The wrapper checks the assembly, upstream classifier and primer hashes.

```bash
python analyses/clermon/scripts/run_local_clermont.py --assembly /path/to/S11_contigs.fasta --primers inputs/clermon/primers.fasta --outdir outputs/clermon
python analyses/clermon/scripts/audit_clermont.py --outdir outputs/clermon
```

Put `blastn` and `makeblastdb` on PATH, or supply `--blast-bin` to the run command. Native Windows BLAST uses a temporary working directory without spaces where needed; returned files are checksum-verified. Keep a new output directory for each run.

The unmodified upstream classifier is pinned to [ClermonTyping commit 69da4db](https://github.com/A-BN/ClermonTyping/tree/69da4db3cd98ab40528bacc289529323ff1407a0). This is its standalone in silico PCR mode without Mash or R. The executable uses mismatch allowance 6 and a perfect 3′ seed of 3; its help text instead states 2 and 5. The audit records both settings and independently checks BLAST alignments and primer pairs against the supplied assembly. Genome intervals are not represented as experimentally amplified sequences.

## Export Galaxy records

The history owner must set **`GALAXY_API_KEY` in the environment**. Never put an API key in the scripts or commit it to Git.

```bash
python scripts/download_galaxy_history.py --inventory-only --outdir outputs/history_inventory
python scripts/download_galaxy_history.py --outdir outputs/history_download
python scripts/export_galaxy_provenance.py --outdir outputs/tool_provenance
```

Downloads include hidden collection members and reject deleted or purged dataset downloads. Provenance includes actual tool settings, available commands and stdout/stderr, versions, metrics, dependencies and wrapper definitions. Missing logs are recorded as unavailable. These commands write only to the chosen local output directory; they do not upload data to GitHub or Figshare.

Galaxy commands can contain server-specific paths and are execution evidence rather than portable local commands. Use the recorded settings and corresponding Galaxy tool versions for new upstream analyses. The wrapper's declared dependencies do not establish a complete runtime package inventory.

## Check existing results or prepare PGAP

```bash
python scripts/compare_fasta.py /path/to/original.fasta /path/to/galaxy_contigs.fasta --report outputs/comparison.json
python scripts/validate_results.py --assembly /path/to/S11_contigs.fasta --quast-report /path/to/report.tsv --pgap-dir /path/to/S11_PGAP --report outputs/validation.json
```

The checks independently calculate assembly statistics, verify existing QUAST metrics, compare annotated genomic sequences, check GFF coordinates and reconcile protein FASTA with GenBank translations. NCBI general-ID notation is normalized bijectively in memory; no input sequence file is modified. The validator targets the preserved original PGAP output, which is distinct from the later identifier-corrected NCBI submission copy.

For a PGAP command plan, provide the existing launcher from [NCBI commit fbadd450](https://github.com/ncbi/pgap/blob/fbadd450b385fa80b71ed11ecaa2a70df15f98ba/pgap.py), the exact assembly and a prepared reference directory:

```bash
python scripts/run_pgap.py --assembly /path/to/S11_contigs.fasta --launcher /path/to/pgap.py --reference-dir /path/to/PGAP_INPUT_DIR --outdir outputs/pgap_plan
```

The wrapper embeds the original tool configuration and pins PGAP release 2026-06-18.build8602 and its container digest. By default it validates inputs and writes a command plan. `--execute` requires an existing compatible Linux/Docker runtime and the matching references. It does not install or download the large PGAP databases. The embedded analysis-use YAML retains the documented placeholder-author and temporary protein-ID validation exemptions; it is not a submission-ready configuration. A full PGAP rerun was not performed during software preparation.

## Validation and license

The figure scripts were rerun with external inputs in a clean environment, and the exports were compared with the locally saved figures. The standalone typing and independent primer audit were rerun; saved QUAST and PGAP/ANI outputs were checked separately. The original unavailable NumPy prerelease pin was replaced with released NumPy 2.2.0. Historical server setup, SSH, recovery and duplicate checkpoint helpers are excluded from this maintained code.

First-party code is [MIT licensed](LICENSE). The vendored ClermonTyping classifier retains GPL-3.0 and its original attribution; see [third-party notices](THIRD_PARTY_NOTICES.md). Cite the software DOI using `CITATION.cff` and cite the NCBI accessions for sequencing data.
