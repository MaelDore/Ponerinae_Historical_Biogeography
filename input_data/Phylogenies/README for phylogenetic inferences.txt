README File for Phylogenetic inferences sub-archive

Doré et al., 2025 - Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity

0_spades_bulk_contigs
- Description: Complete set of SPADes contig files used for analysis. Files labeled using {Genus}_{species}_{ExtractionID} format with all names updated to final manuscript names.

1_spruceup_files
- Description: Data files input into the spruceup trimming program. Taxon names not updated to the final publication names.
- ponerinae-792t-0p-prespruceup-supermatrix-bylocus-partition-file.nexus: By locus partitioning file for the initial concatenated, unfiltered supermatrix.
- ponerinae-792t-0p-prespruceup-supermatrix.phylip: Initial, unfiltered supermatrix file containing all UCE loci and taxa.
- ponerinae-792t-0p-spruceup-configuration-file.txt: Spruceup configuration file used to trim matrix.

2_concatenated_analyses
- Description: All input and log files for the concatenated analyses.
- ponerinae-792t-bylocus-iqtree-analysis-logfile.txt: IQ-Tree log file for the by locus + merging analysis.
- ponerinae-792t-bylocus-iqtree-result-names-updated.tre: Best ML tree for the by locus + merge analysis. Support values are UFBoot/SHaLRT. Terminal names updated to final publication names.
- ponerinae-792t-bylocus-merge-final-models.nexus: Final by locus + merge model file.
- ponerinae-792t-bylocus-partition-file.txt: Initial by locus partitioning file.
- ponerinae-792t-nopart-iqtree-analysis-logfile.txt: IQ-Tree log file for the unpartitioned analysis.
- ponerinae-792t-nopart-iqtree-result-names-updated.tre: Best ML tree for the unpartitioned analysis. Support values are UFBoot/SHaLRT. Terminal names updated to final publication names.
- ponerinae-792t-supermatrix.phylip: Complete concatenated data matrix. Filtered for 75 percent taxon completion. Used in all concatenated analyses.
- ponerinae-792t-swsc-iqtree-analysis-logfile.txt: IQ-Tree log file for the swsc + merging analysis.
- ponerinae-792t-swsc-iqtree-result-names-updated.tre: Best ML tree for the swsc + merge analysis. Support values are UFBoot/SHaLRT. Terminal names updated to final publication names.
- ponerinae-792t-swsc-merge-final-models.nexus: Final swsc + merge model file.
- ponerinae-792t-swsc-partition-file.nexus: Initial swsc partitioning file.

3_misc
- Description: Miscellaneous files.
- phylogenetics-command-list.txt: List of commands used to process reads and extract UCE loci and to analyze the resulting UCE data.
- taxonomy-updates.txt - File linking original taxon names to final publication names following the format: original-name,updated-name .