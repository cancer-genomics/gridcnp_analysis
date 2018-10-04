# gridcnp_analysis

Currently, the script scripts/1-collectSampleInfo.R takes spreadsheets/info in the data directory and creates .csv files that will contain all the data we need for a particular analysis. Right now, we are planning on 3 analyses:

1. Ovarian cell lines (with both targeted + WGS)--no matched normals, done with Cp target set (need to find probe info). Unclear whether there is panel of normals that were done alongside this
2. CRC tissue for which we have matched plasma & normal--also have panel of normals, done with CpPa2 target set (we have this probe info)--should find out which genome data was aligned to
3. TCGA data--unclear which samples, most likely will verify CNV calls with SNP6 segmentation calls from Absolute. Will use WES--potentially can choose samples such that there are clear batches of cases & normals (either different sequencing dates or different tcga projects)

## TODO

* Find probes & comparable panel of normals for each analysis
* Build .csv files with sample info for all analyses
* Transfer BAMs into /dcl01/scharpf1/bams
* Do analyses :) 
