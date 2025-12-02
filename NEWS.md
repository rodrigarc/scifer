# scifer 1.9.1

* fix igblast cdr3 bug with the contribution of Sebastian Ols (@sebbeols)
	- this affects igblast.R function, but mainly the igblast wrapper igblastwrap.py

# scifer 1.8.1

* Update citation, scifer is now published in ImmunoInformatics (DOI: 10.1016/j.immuno.2024.100046). 

# scifer 1.8.0

NEW FEATURES
* New igblast wrapper function to run igblast on fasta sequences
* igblast function returns the output from a local igblast run with any custom database
* New function requires additional steps for ARM processors, please read the README

# scifer 1.7.1 (2024-04-24)

BUG FIXES
* Fix compensation bug due to different compensation matrix positions in different `.fcs` file versions
* update minor nomenclature to AIRR-standards, column `file_ID` to `sequence_id`
* Fix error of `HC_individualized_report.Rmd` file not being found

NEW FEATURES
* Update quality thresholds based on a large B cell receptor test dataset

## scifer 0.99.4 (2022-10-31)

* Version bump to bioconductor, updated description


## scifer 0.99.3 (2022-09-20)

BUG FIXES

* Fix NOTE when building on windows computer

## scifer 0.99.2 (2022-09-08)

NEW FEATURES

* Update vignette according to Bioconductor
* Added unit tests for functions

SIGNIFICANT USER-VISIBLE CHANGES

* Split function fcs_processing() into two functions:
  * fcs_processing()
  * fcs_plot()

BUG FIXES

* Minor changes to avoid possible bugs



## scifer 0.99.1 (2022-08-26)

* Version bump to fix Bioconductor building error

## scifer 0.99.0 (2022-08-10)

* Adapting package for Bioconductor submission

## scifer 0.98.0 (2022-08-09)

* Major changes in functions
  * nomenclature and all documentation were added following R style guidelines
  * Processing flow cytometry data and adding to the quality_report() function

## scifer 0.1.0 (2022-04-15)

* Submitted to GitHub

## scifer 0.0.1 (2020-08-10)

* Primordial package created with simple functions and bugs
