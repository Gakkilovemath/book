# book
Example programs for the book "Nonparametric Estimation under Shape Constraints" (2014),
Piet Groeneboom and Geurt Jongbloed. Cambridge University Press, 2014.

For running the R scripts, like bootstrap.R in the directory rbootstrap, one needs three
files: the R script, the input file (for bootstrap.R this is: dataThai.txt) and the C++
file (in this case icm.cpp). One needs to put the three files, R script, C++ file (which
has the extension .cpp) and input file (with extension .txt), into one directory, just like
it is done in the repository. The bootstrap procedure draws 1000 bootstrap samples and 
computes 95% confidence intervals for the hazards.

One needs recent versions of R and then the examples run in any case on Mac, Windows and
Linux (also on virtual versions of these sytems). One also needs the Rcpp package
for the command "sourceCpp", which makes the connection between the C++ code and R.
If one gets an error message of the type "sourceCpp cannot be found", the script will not work.
The directory "non-grouped" has the data file in a format where the observation times
are non-unique and where the ties at the observations are computed by the program.
The directory "grouped" has the data file in a format where the observation times
are unique and where the number of observations of a certain type are given at each
observation. Note the different structures of the input files in the grouped and non-grouped
directories! The rbootstrap directory has the input format in the "ungrouped style".
For more information on the present examples, see Rscripts.pdf.

In copying the files from the site, one should push the button "raw" above the file in GitHub,
to avoid getting html code instead of the file one needs.

