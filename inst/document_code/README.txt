################################################################################
##     
##  "Semi-Supervised: Scalable Semi-Supervised Routines for Real Data Problems"
##
## The instructions for executing the benchmark results from Section 2 
## and Figure results in Sections 1 and 2 are provided. 
##
################################################################################

Overview:
----------
All figures and tables can be regenerated from source. The way it works is that
the generate* R-scripts create the figures/tables using the csv files. The 
original csvs files are included from the execution on my linux machine. The 
csv files can be regenerated (takes 30 hours on my machine) on your 
machine using the run* R-scripts. 

Table 3 requires special consideration. You need to configure R to use 
the MKL or vector libraries depending on your computer. The results in the
csv directory now is from my execution under linux using MKL (column 4). Once 
R is set-up just run the script.

For best results, run these scripts on a computer that is using minimal 
functionality. 

If you want to quickly generate a row in my corresponding csv file, just
execute with nruns=1. It take about 20 minutes to cycle through all 10 data
sets for the performance table and then check whether the row matches up.

Note: Sometimes performances are different on different machines due to 
usage of more powerful BLAS libraries in R. My run was executed on 
Microsoft R Open 3.2.5 using MKL libraries with 8 threads.
 
Structure:
----------

Directory .: root directory. Open up a terminal and execute the commands from 
	     this directory, or set working directory in R to this directory
	     for execution of scripts.

R-Script Generate_Figures.R: The R-script generates all Figures in the paper
	 		     using the appropriate files in the csv directory.

R-Script Generate_Tables.R: The R-script generates the Tables from Section 2
	 		    using the appropriate files in the csv directory.

R-Script run_preformance.R: The R-script executes the s4pm, agraph, and glmnet
                            on each of the data sets given in Table 1. The
			    nruns variable in the R-script determines the 
			    number of runs. The csvs for all 100 executions on
			    my linux machine are provided. If you wish to verify
			    my results, it is easier to set nruns=1 and see that
			    1st runs are equal. The total 100 runs took about 
			    28 hours on my machine (the Song data took about 16). 
			    The execution of this script over-writes the csv 
			    files in the csv directory.

R-Script run_casp_CPU.R:   The R-script executes the s4pm, agraph, and glmnet
	 		   on the CASP data set for 10 runs with m=200. The 
			   script can be run on any R setup, i.e., change
			   the R BLAS library to execute the script. The 
			   execution of this script overrides the 
			   corresponding csv file.

R-Script run_casp_m.R      The R-script executes the s4pm, agraph, and glmnet
	 		   on the CASP data set for 10 runs with varying m.
			   The execution of this script overrides the 
			   corresponding csv file.

Directory csvs: The csvs resulting from the execution of the R-script. The 
	  	csvs from my computer using Linux are given. These same 
		results populated the Tables in Section 2. 

Directory source: The source files for an execution.

R-Script data_reader.R: This script gets the data files either from available 
			in R packages, or  downloaded from the internet. 
			The powerplant data could also have been downloaded 
			but the UCI repository's version is not in a 
			convenient form. As a result this data set was 
			added to the package.

R-Script base.R: The base R-script for executing the s4pm, agraph, and 
	 	 glmnet.

