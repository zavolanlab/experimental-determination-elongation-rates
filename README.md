# Experimental determination of translation elongation rates
This is a repository that contains scripts for the determination of the translation elongation rates per gene from ribosome profiling runoff experiments. To date, these techniques have not been successfully applied on a per-gene level since ribosome profiling is agnostic to the situation of ribosomes on a single transcript.
Therefore, it is hard to extract information about ribosome stalling on an individual transcript from ribosome profiling data.
A starting point for our consideration was given by Ingolia *et. al.*: Ribosome Profiling of Mouse Embryonic Stem Cells Reveals the Complexity and Dynamics of Mammalian Proteomes, Cell (2011), where the so-called SL method is introduced.
This method is applied to single genes in Dana, Tuller: Determinants of Translation Elongation Speed and Ribosomal Profiling Biases in Mouse Embryonic Stem Cells, PLoS Comp. Biol. (2012).
Last but not least, we also endeavoured a modeling effort ourselves, which we will discuss in the last paragraph.

## Ribosome profiling runoff experiments
The general philisophy behind all determinations of translation elongation speed via ribosome profiling is the following: You block translation initiation [^1] with harringtonine and apply cycloheximide later on to block all ribosome movement. This can be done at different time differences. The crucial step is now to construct an observable from ribosome profiling data that changes at the speed of runoff. The following two paragraphs illustrate the SL method using the starting location as an observable that moves at the speed of translation elongation, and the CR method that presents the ratio of riboseq coverage of the first half of the transcript normalized by the coverage of the second half. The latter quantity should decay at the speed of translation elongation.

[^1]: There is one important caveat to this: Harringtonine actually blocks translation initiation by keeping the ribosome on the start codon to make its first elongation step. Therefore, the start codons are occupied, and no new ribosomes can enter the downstream parts of the CDS. However, that means that the first codon has to be discarded from all sorts of runoff analyses.

## SL method
This is the original method to compute global estimates of the translation elongation rates using riboseq runoff experiments from Ingolia *et. al.*. 2011. “Ribosome Profiling of Mouse Embryonic Stem Cells Reveals the Complexity and Dynamics of Mammalian Proteomes.” Cell 147 (4). Here, it is applied on a single-transcript-level (previously done by Dana, Tuller. 2012. "Determinants of Translation Elongation Speed and Ribosomal Profiling Biases in Mouse Embryonic Stem Cells". PLoS Comp. Biol.). The script that implements this method is called ```differential_translation_speed.py```.

### Input
Required input files are:
1. a .tsv file with the transcript ID, corresponding gene ID, and CDS coordinates on the transcript
2. a dictionary with time separations, .bam/.sam mapping files at the corresponding time separation, and the corresponding .json file with the p-site-offsets
3. column names indicating the method of reference location determination (```SL_location```in our case) and the corresponding error column
4. the size of bins for depletion curves (15nts in our case)
5. the number of jackknife bins (for error estimation, usually 10-20)
6. a cutoff value for chi^2_reduced (transcripts with fits that produce worse values are discarded)

### Procedure
The script processes the input as follows:
1. reading the CDS coordinates, transcript ID, and gene ID into internal dictionaries
2. creating profiles from the mapping (and offset) files
3. computing the (SL) locations for each transcript at each jackknife bin
4. compute the error (from jackknife bins) of for each transcript and time point
5. perform the fit of each transcript

### Output
The following output is created:
1. the SL curves in file ```p-sites.tsv```
2. the (SL) locations in file ```locations.tsv```
3. the locations with errors in ```locations_with_errors.tsv```
4. the elongation speeds in ```elongation_speeds.tsv```


## CR method
The goal of the CR method is to overcome some of the shortcomings of the SL method:
- The SL method uses one particular curve as reference (t = 0). This introduces additional noise.
- The recovery threshold of the SL method (0.5) is an arbitrary choice. Different choices have a high impact on the result.
- Genes that have a high ribosome flux (speed times coverage) dominate the global estimate. The average of single-gene SL estimates is ~ 3 times smaller than the global estimate as described above.
- It can only be applied to long enough genes (length of the CDS > 3000nt).
- It requires sufficient coverage of a transcript.

The alternative definition of a quantity that varies with the speed of translation elongation is the ratio of the coverage of the first half of the transcript divided by the coverage of the second half. This quantity should linearly fall with the translation elongation rate being the negative slope.

The method is implemented in ```compute_elong_speed_coverage_ratio.py```.

### Input
Required input files are:
1. a .tsv file with the transcript ID, corresponding gene ID, and CDS coordinates on the transcript
2. a dictionary with time separations, .bam/.sam mapping files at the corresponding time separation, and the corresponding .json file with the p-site-offsets
3. the number of jackknife bins (for error estimation, usually 10-20)
4. a cutoff value for chi^2_reduced (transcripts with fits that produce worse values are discarded)


### Procedure
The script processes the input as follows:
1. reading the CDS coordinates, transcript ID, and gene ID into internal dictionaries
2. computing the coverage ratios from .sam and .json files
3. compute the error (from jackknife bins) of for each transcript and time point
4. perform the fit of each transcript


### Output
The following output is created:
1. the coverage ratios in file ```p-sites.tsv```
2. the coverage ratios with errors in ```locations_with_errors.tsv```
3. the elongation speeds in ```elongation_speeds.tsv```


# License
The code is published under the MIT license.
