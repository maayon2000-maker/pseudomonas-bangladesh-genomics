# Pseudomonas aeruginosa from Bangladesh, India and Pakistan: public-genome comparison

Analysis files for the manuscript on public *Pseudomonas aeruginosa* genomes from Bangladesh, India and Pakistan.
The study asks whether Bangladeshi genomes show a country-specific carbapenemase-gene profile or lineage once
clone and study structure are controlled.

Authors: Mohammad Motahar Hossain Ayon and Rupayan Parash Bhattacharjee (Department of Microbiology, BRAC University, Dhaka).

## Data source

All genomes are public GenBank assemblies selected from the NCBI Pathogen Detection isolates browser,
dataset PDG000000036.1843, downloaded on 22 November 2026. No new sequencing was done.
Accessions and reasons for exclusion are in `data/Supplementary_Tables_S1-S3.xlsx`.

| Table | Content |
|---|---|
| S1 | The 985 main-set genomes (60 Bangladesh, 674 India, 251 Pakistan) |
| S2 | The 122 genomes used in the tree |
| S3 | The 75 retrieved genomes outside the main set, with the reason |

`data/clusters_main_clean.csv` is the main-set table used by the scripts: assembly quality, sequence type,
country, BioProject, year, source, carbapenemase call and clone-group labels at three ANI thresholds.

## Analysis steps and software

| Step | Tool (version) | Purpose |
|---|---|---|
| Download | NCBI Datasets command-line tool | Retrieve assemblies by accession |
| Quality control | CheckM 1.3.2 (NCBI values for 849 genomes) | Completeness at least 95%, contamination below 5% |
| Sequence type | mlst 2.35.0 | PubMLST seven-locus scheme |
| Resistance genes | AMRFinderPlus 4.2.7, database 2026-10-29.1 | Acquired genes; intrinsic elements removed |
| Clone groups | skani 0.4.0 | ANI at least 99.9, 99.95 and 99.99%, alignment fraction at least 80%, single linkage |
| SNP distances | Snippy 4.6.0, snp-dists 0.8.2 | Pairs within shared clone groups |
| Tree | IQ-TREE 3.1.4, GTR+F+ASC+G4, 1,000 ultrafast bootstraps | Core-SNP tree of 122 genomes |
| Statistics and figures | R 4.5.1; Python with SciPy | Wilson intervals, Fisher's exact tests, Holm correction |

The carbapenemase definition is any gene symbol beginning blaNDM, blaVIM, blaIMP, blaKPC, blaSPM, blaGIM,
blaSIM, blaDIM or blaAIM. Full details are in the manuscript Methods.

## Files

- `data/`: supplementary tables and the main-set table.
- `scripts/`: the Python scripts used for the recount, sensitivity analyses, Holm correction, power calculation
  and the SNP, clone-sharing and tree figures. They were run from the project folder and read the tables named in
  their first lines; edit the paths if you run them elsewhere.

## Limits to keep in mind

Genes are predicted from public genomes that were not sampled at random. The analysis describes genotypes,
not resistance, and cannot show transmission between patients.

## Use of AI

Claude (Anthropic) helped write analysis code and draft and edit text. The authors carried out the analyses,
checked the outputs against the source data and take responsibility for the content.

## Citation

See `CITATION.cff`. The archived release has a Zenodo DOI: [add after the first release].

## Licence

Code and tables: MIT licence (`LICENSE`). Genome sequences remain under the terms of the public databases.
