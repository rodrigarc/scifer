## Folder organization

### Inside this folder there are 3 different subfolders:

* script
    * contains a text file describing in detail how the raw data present in `extdata` was generated 

* rmd
    * contains Rmarkdown files that are used by `quality_report()` function to generate quality reports for all the sequences (HC_report.Rmd) and for individual reports for each sample (HC_individualized_report.Rmd)

* extdata
    * contains folders with raw data from single-cell sorted B cell receptor sanger sequencing (sorted_sangerseq) and flow cytometry data from when those cells were sorted (fcs_index_sorting)
        - fcs_index_sorting: `.fcs` files from one sample sorted in different plates. eg. E18_02 means sample E18, plate 02.
        
        - sorted_sangerseq: contains 2 folders for matching one of the samples in `fcs_index_sorting`.
            - E18_C1: contains raw `.ab1` files from 12 wells sequenced by sanger sequencing
            - E18_C1_R: contains raw `.ab1` files from 12 wells that were re-sequenced (repeated) to try to improve sequence quality
	

