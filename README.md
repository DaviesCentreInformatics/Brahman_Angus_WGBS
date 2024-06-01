# Brahman_Angus_WGBS
 Repo for code used in paper

Much of the code contained in the repo is based on tutorials and examples from the various tools used. Feel free to refer to those for additional information.

- Code to map WGBS reads to reference genome is [here](https://github.com/DaviesCentreInformatics/Brahman_Angus_WGBS/blob/main/code/main.nf)

- Code to convert `MethylDackel` output into `DNMTools` input is [here](./code/methylDackel2DNMtools.ipynb)

- Code to artificially mutate fasta file is [here](./code/mutateFasta.sh)

- Code to extract 1002bp CpG windows is [here](./code/createCpGFasta.ipynb).

- Code to map CpGs to reference is [here](./code/minimap2_align.sh)

- Code to filter alignments from CpG mapping is [here](./code/filterBam.py)

- Code to determine reference impact on methylation is [here](.code/methylation_ref_impact.ipynb)

- Code to identify DMRs between groups is [here](./code/DMR/aligned2Angus/) and [here](./code/DMR/aligned2Brahman/)
  
- Revised code to identify DEGs is [here](./code/revised_DEGs.R)

- Revised code to determine SNP and SV enrichment is [here](./code/enrichment.ipynb)

- Markdown file outlining steps taken to identify DMRs between references for the same group of samples is [here](./code/findComparableRegions.md)

- Bed files containing the shared and breed-specific CpGs are here:
>- [Brahman.breed-specific](CpGs/Brahman.Breed.Specific.CpGs.bed)
>- [Brahman.shared](./CpGs/Brahman_shared/)
>- [Angus.breed-specific](CpGs/Angus.Specific.CpGs.bed)
>- [Angus.shared](./CpGs/Angus_shared/)

- Bed files containing the DMRs identified using shared CpGs are [here](DMRs/):
