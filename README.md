# Brahman_Angus_WGBS
 Repo for code used in paper

Much of the code contained in the repo is based on tutorials and examples from the various tools used. Feel free to refer to those for additional information.

- Code to convert `MethylDackel` output into `DNMTools` input is [here](./code/methylDackel2DNMtools.ipynb)

- Code to extract 1002bp CpG windows is [here](./code/createCpGFasta.ipynb).

- Code to map CpGs to reference is [here](./code/minimap2_align.sh)

- Code to filter alignments from CpG mapping is [here](./code/filterBam.py)

- Code to identify DMRs between groups is [here](./code/DMR/aligned2Angus/) and [here](./code/DMR/aligned2Brahman/)
  
- Code to identify hypomethylated regions is [here](./code/idHMRs.sh)

- Code to identify enriched 3mers is [here](./code/enriched3mers.ipynb)

- Code to identify DEGs is [here](./code/idDEGs.r)

- Bed files containing the shared and breed-specific CpGs are here:
>- [Brahman.breed-specific](CpGs/Brahman.Breed.Specific.CpGs.bed)
>- [Brahman.shared](./CpGs/Brahman_shared/)
>- [Angus.breed-specific](CpGs/Angus.Specific.CpGs.bed)
>- [Angus.shared](./CpGs/Angus_shared/)