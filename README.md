# mrclean
Two Mixed Integer Programs for cleaning a data file.

## To Use
Configure the Makefile with the locaion of IBM ILOG CPLEX libraries and binary

Compile with the Makefile by navigating to the root directory and entering: make

Update configuration file

Run the program. For example enter: ./mrclean data_file.tsv 0.05 NA 100 200 1 1

or                                  ./mrclean data_file.tsv 0.05 NA 100 200 1 1 greedy.sol

## Inputs
<data_file> - Tab seperated data file to clean.

<max_missing> - Maximum percent of data allowed in each row and column in the cleaned data file.

<row_lb> - Minimum number of rows allowed in a solution

<col_lb> - Minimum number of columns allowed in a solution

<miss_symbol> - String used to indicate missing data in <data_file>.

<num_header_rows> - Number of header rows in <data_file>.

<num_header_cols> - Number of header rows in <data_file>.

<incument_file> - Optional input used to provide starting integer solution to the MIPs. File should contian two binary vectors represented as lines of tab-seperated numbers. The first line represents the retained rows in the solutio and the second line represents the retained columns in the solution. 

## Configuration File
PRINT_SUMMARY - Boolean controlling if a summary of the cleaning results is printed.

WRITE_STATS - Boolean controlling if the statistic of the cleaning algorithms are written to a file.

RUN_ROW_COL - Boolean controlling if the RowCol IP is executed.

RUN_ELEMENT - Boolean controlling if the Element IP is executed.

## Outputs
RowCol_summary.csv - Statistics file for the RowCol IP containing <data_file>, <max_missing>, <miss_symbol>, run time, number of valid elements kept number of rows kept, and number of columns kept.

Element_summary.csv - Statistics file for the RowCol IP containing <data_file>, <max_missing>, <miss_symbol>, run time, number of valid elements kept number of rows kept, and number of columns kept.

<data_file>_cleaned.tsv - Tab-seperated cleaned data file.

## Notes
Recommend using [mrclean-greedy] to provide incumbent for MIPs

Requires IBM ILOG CPLEX

<data_file> should be tab sperated.
