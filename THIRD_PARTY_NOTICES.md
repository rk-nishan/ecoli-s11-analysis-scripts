# Third-party attribution

The unmodified `clermont.py` classifier is from [A-BN/ClermonTyping](https://github.com/A-BN/ClermonTyping), commit `69da4db3cd98ab40528bacc289529323ff1407a0`. Its GPL-3.0 license is retained alongside the source. Its optional Mash binary, R scripts, primer FASTA and reference databases are not bundled.

`fetch_references.py` obtains the required public references separately and checks their pinned SHA-256 hashes: CARD ARO 4.0.2, ResFinder antibiotic-class metadata at commit `65c9ba3` and the upstream ClermonTyping primers. These references retain their upstream terms, including CARD ARO's CC BY 4.0 attribution requirement and the ClermonTyping project's GPL-3.0 license. The helper does not redistribute the full CARD or ResFinder sequence databases in this repository.

NCBI PGAP, its container and large reference databases are supplied separately under NCBI's upstream terms. Their required source commit, release and container digest are recorded in the execution wrapper. Matplotlib and its DejaVu Sans font are installed as dependencies; no font package is redistributed.

The root MIT license applies to first-party code and does not override third-party licenses.
