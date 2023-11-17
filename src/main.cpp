#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>

#include "DataContainer.h"
#include "GreedySolver.h"
#include "RowColSolver.h"
#include "ElementSolver.h"
#include "Timer.h"
#include "ConfigParser.h"
#include "CleanSolution.h"

void summarize_results(const DataContainer &data,
                       const std::string &na_symbol,
                       const double elapsed_cpu_time,
                       const std::size_t num_rows_kept,
                       const std::size_t num_cols_kept,
                       const std::vector<bool> rows_to_keep,
                       const std::vector<bool> cols_to_keep);

void write_stats_to_file(const std::string &file_name,
                         const std::string &data_file,
                         const double max_perc_missing,
                         const double time,
                         const std::size_t num_valid_element,
                         const std::size_t num_rows_kept,
                         const std::size_t num_cols_kept);

int main(int argc, char *argv[]) {
  if (!((argc == 6) || (argc == 7))) {
    fprintf(stderr, "Usage: %s <data_file> <max_missing> <miss_symbol> <num_header_rows> <num_header_cols> (opt)<incument_file>", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Read in user inputs
  std::string data_file(argv[1]);
  double max_perc_missing = std::stod(argv[2]);
  std::string miss_symbol(argv[3]);
  std::size_t num_header_rows = std::stoul(argv[4]);
  std::size_t num_header_cols = std::stoul(argv[5]);
  std::string incumbent_file = "";
  bool SEEDING_MIP = false;
  if (argc == 7) {
    incumbent_file = argv[6];
    SEEDING_MIP = true;
  }

  // Check user inputs
  if ((max_perc_missing < 0) || (max_perc_missing > 1)) {
    throw std::runtime_error("max_perc_missing must be between [0,1].");
  }

  // Check config parser values
  ConfigParser parser("config.cfg");
  const bool PRINT_SUMMARY = parser.getBool("PRINT_SUMMARY");
  const bool WRITE_STATS = parser.getBool("WRITE_STATS");
  const bool RUN_ROW_COL = parser.getBool("RUN_ROW_COL");
  const bool RUN_ELEMENT = parser.getBool("RUN_ELEMENT");

  fprintf(stderr, "\n\nUser Inputs:\n");
  fprintf(stderr, "\t%s is MISSING SYMBOL\n", miss_symbol.c_str());
  fprintf(stderr, "\tAllowing %lf maximum missing data per row/col\n", max_perc_missing);
  fprintf(stderr, SEEDING_MIP ? "\tSeeding MIPs with solution\n" : "\tNOT seeding MIPS\n");

  Timer timer;

  // Read in the data
  DataContainer data(data_file, miss_symbol, num_header_rows, num_header_cols);
  fprintf(stderr, "\nSummary of raw data file\n");
  fprintf(stderr, "\tNum rows before cleaning: %lu\n", data.get_num_data_rows());
  fprintf(stderr, "\tNum cols before cleaning: %lu\n", data.get_num_data_cols());
  fprintf(stderr, "\tTotal number of valid entries: %lu\n", data.get_num_valid_data());
  fprintf(stderr, "\tMax perc missing - row: %lf\n", data.get_max_perc_miss_row());
  fprintf(stderr, "\tMax perc missing - col: %lf\n", data.get_max_perc_miss_col());
  fprintf(stderr, "\tMin perc missing - row: %lf\n", data.get_min_perc_miss_row());
  fprintf(stderr, "\tMin perc missing - col: %lf\n", data.get_min_perc_miss_col());
  fprintf(stderr, "\n");

  CleanSolution clean_sol(data.get_num_data_rows(), data.get_num_data_cols());
  if (!incumbent_file.empty()) {
    clean_sol.read_from_file(incumbent_file);
  }
  
  // Create and solve IP
  std::size_t rc_rows = 0, rc_cols = 0, rc_val_elements = 0;
  double rc_time = 0.0;
  std::vector<bool> rc_rows_to_keep(data.get_num_data_rows(), false), rc_cols_to_keep(data.get_num_data_cols(), false);
  if (RUN_ROW_COL) {
    fprintf(stderr, "Starting Row Col IP...\n");
    RowColSolver rc_solver(data, max_perc_missing);

    if (SEEDING_MIP) {      
      rc_solver.set_incumbent(clean_sol.get_rows_to_keep(), clean_sol.get_cols_to_keep());
    }
    
    timer.restart();
    rc_solver.solve();
    timer.stop();
  
    rc_rows_to_keep = rc_solver.get_rows_to_keep();
    rc_cols_to_keep = rc_solver.get_cols_to_keep();
    rc_rows = rc_solver.get_num_rows_to_keep();
    rc_cols = rc_solver.get_num_cols_to_keep();
    rc_time = timer.elapsed_cpu_time();
    rc_val_elements = data.get_num_valid_data_kept(rc_rows_to_keep, rc_cols_to_keep);

    if (PRINT_SUMMARY) {
      summarize_results(data, miss_symbol, rc_time, rc_rows, rc_cols, rc_rows_to_keep, rc_cols_to_keep);
    }
  }
  

  // Create and solve IP
  std::size_t element_rows = 0, element_cols = 0, element_val_elements = 0;
  double element_time = 0.0;
  std::vector<bool> element_rows_to_keep(data.get_num_data_rows(), false), element_cols_to_keep(data.get_num_data_cols(), false);
  if (RUN_ELEMENT) {
    fprintf(stderr, "Starting Element IP...\n");
    ElementSolver element_solver(data, max_perc_missing);

    if (SEEDING_MIP) {
      element_solver.set_incumbent(clean_sol.get_rows_to_keep(), clean_sol.get_cols_to_keep());
    }
    
    timer.restart();
    element_solver.solve();
    timer.stop();

    element_rows_to_keep = element_solver.get_rows_to_keep();
    element_cols_to_keep = element_solver.get_cols_to_keep();
    element_rows = element_solver.get_num_rows_to_keep();
    element_cols = element_solver.get_num_cols_to_keep();
    element_time = timer.elapsed_cpu_time();
    element_val_elements = data.get_num_valid_data_kept(element_rows_to_keep, element_cols_to_keep);

    if (PRINT_SUMMARY) {
      summarize_results(data, miss_symbol, element_time, element_rows, element_cols, element_rows_to_keep, element_cols_to_keep);
    }
  }  
  
  // Find the best result and use it to write the cleaned matrix
  std::size_t rc_num_valid = data.get_num_valid_data_kept(rc_rows_to_keep, rc_cols_to_keep);
  std::size_t element_num_valid = data.get_num_valid_data_kept(element_rows_to_keep, element_cols_to_keep);

  std::size_t lastindex = data_file.find_last_of("."); 
  std::string cleaned_file = data_file.substr(0, lastindex) + "_cleaned.tsv";

  // Determine best num_valid
  int best_alg = 1; // element
  std::size_t best_num_valid = element_num_valid;

  if (rc_num_valid > best_num_valid) {
    best_num_valid = rc_num_valid;
    best_alg = 0;
  }

  if (best_alg == 0) {
    data.write(cleaned_file, rc_rows_to_keep, rc_cols_to_keep);
  } else {
    data.write(cleaned_file, element_rows_to_keep, element_cols_to_keep);
  }
  
  // Wrtie statistics to file
  if (WRITE_STATS) {
    write_stats_to_file("RowCol_summary.csv", data_file, max_perc_missing, rc_time, rc_val_elements, rc_rows, rc_cols);
    write_stats_to_file("Element_summary.csv", data_file,max_perc_missing, element_time, element_val_elements, element_rows, element_cols);
  }
    
  return 0;
}

//------------------------------------------------------------------------------
// Print summary of algorithm to screen
//------------------------------------------------------------------------------
void summarize_results(const DataContainer &data,
                       const std::string &na_symbol,
                       const double elapsed_cpu_time,
                       const std::size_t num_rows_kept,
                       const std::size_t num_cols_kept,
                       const std::vector<bool> rows_to_keep,
                       const std::vector<bool> cols_to_keep)
{
  fprintf(stderr, "\tTook %lf seconds\n", elapsed_cpu_time);
  fprintf(stderr, "\tNum rows after cleaning: %lu\n", num_rows_kept);
  fprintf(stderr, "\tNum cols after cleaning: %lu\n", num_cols_kept);
  fprintf(stderr, "\tNumber of valid elements after cleaning: %lu\n", data.get_num_valid_data_kept(rows_to_keep, cols_to_keep));
  fprintf(stderr, "\n");
}

//------------------------------------------------------------------------------
// Write the statistics to a file
//------------------------------------------------------------------------------
void write_stats_to_file(const std::string &file_name,
                         const std::string &data_file,
                         const double max_perc_missing,
                         const double time,
                         const std::size_t num_valid_element,
                         const std::size_t num_rows_kept,
                         const std::size_t num_cols_kept) {
  FILE *summary;
  
  if((summary = fopen(file_name.c_str(), "a+")) == nullptr) {
    fprintf(stderr, "Could not open file (%s)", file_name.c_str());
    exit(EXIT_FAILURE);
  }

  fprintf(summary, "%s,%lf,%lf,%lu,%lu,%lu\n", data_file.c_str(), max_perc_missing, time, num_valid_element, num_rows_kept, num_cols_kept);

  fclose(summary);
}