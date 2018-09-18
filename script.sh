cd evolutionUmea
R --no-save << FLAG1
   library(roxygen2)
   roxygenise()
FLAG1

cd -

R CMD build  evolutionUmea
R CMD INSTALL -l /home/xavier/DIR.myrlibrary evolutionUmea_0.97.tar.gz

#In R you run install.packages(path_to_file, repos = NULL, type="source")
