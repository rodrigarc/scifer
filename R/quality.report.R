quality_report <- function(data_folder, outputfile, output_dir){

  rmarkdown::render("inst/rmd/HC_report.Rmd",
                    output_dir = output_dir,
                    params = list(
                      data_folder = data_folder
                    ),
                    output_file = outputfile)
}


quality_report(data_folder = "test/",
                outputfile = paste0(Sys.Date(),"_HC-QC-report.html"),
                output_dir = "test")

## Basic idea picked up from some random code as following

#render_report = function(region, year) {
#  rmarkdown::render(
#    "MyDocument.Rmd", params = list(
#      region = region,
#      year = year
#    ),
#    output_file = paste0("Report-", region, "-", year, ".pdf")
#  )
#}
