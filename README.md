# crispr_sgrna_library_design #

Several R scripts to make a whole genome library.

**ATTENTION**: The scripts are not well tested yet. And further
development is demanded.

## Requirements ##

### Packages ###

The scripts rely on several packages:
* [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
* [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
* [GenomicRanges](https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
* [DBI](https://dbi.r-dbi.org)
* [RSQLite](https://rsqlite.r-dbi.org)
* [ggplot2](https://ggplot2.tidyverse.org/)
* [dplyr](https://dplyr.tidyverse.org/)
* [tidyr](https://tidyr.tidyverse.org/)
* [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
* [sqlite3](https://sqlite.org/index.html)

### Data ###

To run the scripts also need reference data from NCBI or other
resources.
* reference genome in fasta format
* reference annotation in gff3 format

## Usage ##

1. set parameters in the `parameters.r` script.
2. run `Rscript spacer_db.r` to generate the database that stored
   informations for spacers.
3. (optional) run `Rscript select_library.r` to select spacers for a
   whole genome library.

## Notes ##

* Off-target finding uses bowtie.
* The calculation of spacer specific score is adopted from Perez et
  al. 2017 and Doench et al. 2016.
* The calculation of mismatch score is adopted from Stemmer et al. 2015

## Reference ##

* Perez, A.R., Pritykin, Y., Vidigal, J.A., Chhangawala, S., Zamparo, L., Leslie, C.S., and Ventura, A. (2017). GuideScan software for improved single and paired CRISPR guide RNA design. Nat Biotechnol 35, 347–349.
* Doench, J.G., Fusi, N., Sullender, M., Hegde, M., Vaimberg, E.W., Donovan, K.F., Smith, I., Tothova, Z., Wilen, C., Orchard, R., et al. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. Nat Biotechnol 34, 184–191.
* Stemmer, M., Thumberger, T., Del Sol Keyer, M., Wittbrodt, J., and Mateo, J.L. (2015). CCTop: An Intuitive, Flexible and Reliable CRISPR/Cas9 Target Prediction Tool. PLoS One 10, e0124633.

## LICENSE ##

See LICENSE file.
