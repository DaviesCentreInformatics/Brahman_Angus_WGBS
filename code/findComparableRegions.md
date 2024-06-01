# Steps used to identify regions between Angus and Brahman reference genomes for DMR analysis

```shell
# Make 1000bp windows
bedtools makewindows -g /path/to/Angus.ARS-UCD1.2.orientation.fa.fai -w 100000 -i winnum > Angus_100kb_windows.bed
bedtools makewindows -g /path/to/Bos_indicus_hybrid.UOA_Brahman_1.ARS-UCD1.2.orientation.dna.toplevel.fa.fai -w 100000 -i winnum > Brahman_100kb_windows.bed

# Get the fasta sequences for the windows we just made
bedtools getfasta -fi /path/to/Angus.ARS-UCD1.2.orientation.fa -bed Angus_100kb_windows.bed -fo Angus_100kb_windows.fa
bedtools getfasta -fi /path/to/Bos_indicus_hybrid.UOA_Brahman_1.ARS-UCD1.2.orientation.dna.toplevel.fa -bed Brahman_100kb_windows.bed -fo Brahman_100kb_windows.fa

# Map them to the opposite reference
minimap2 -t 8 -c -x map-ont /path/to/Angus.ARS-UCD1.2.orientation.fa Brahman_100kb_windows.fa > Brahman_100kb_windows_to_Angus_ref.paf
minimap2 -t 8 -c -x map-ont /path/to/Bos_indicus_hybrid.UOA_Brahman_1.ARS-UCD1.2.orientation.dna.toplevel.fa Angus_100kb_windows.fa > Angus_100kb_windows_to_Brahman_ref.paf
```

Then use the code in `code/matchWindows.ipynb` to match 100kb windows between Brahman and Angus references

``` shell
# Get methylation values for the windows
bedtools intersect -a filtered_Angus_100kb_windows_w_Brahman_coords.bed -b angus_brahman_coord_methto_bed.bed -loj > angus_cpg_methylation_mapped_to_brahman.bed
bedtools intersect -a filtered_Angus_100kb_windows.bed -b angus_angus_coord_methto_bed.bed -loj > angus_cpg_methylation_mapped_to_angus.bed

samples=("AxA" "AxB" "BxA" "BxB")
for sample in ${samples[@]};
do
	bedtools intersect -a filtered_Angus_100kb_windows.bed -b ${sample}_angus_angus_coord_methto_bed_100kb.bed -loj > ${sample}_cpg_methylation_mapped_to_angus.bed
	bedtools intersect -a filtered_Angus_100kb_windows_w_Brahman_coords.bed -b ${sample}_angus_brahman_coord_methto_bed_100kb.bed -loj > ${sample}_cpg_methylation_mapped_to_brahman.bed
done
```

Finally, use the code in `code/performTest.ipynb` to determine significant DMRs.
