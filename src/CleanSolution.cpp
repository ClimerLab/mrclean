#include "CleanSolution.h"

#include <fstream>
#include <sstream>
#include <algorithm>


CleanSolution::CleanSolution(const std::size_t num_rows,
                             const std::size_t num_cols) :  rows_to_keep(num_rows, 0),
                                                            cols_to_keep(num_cols, 0) {}

CleanSolution::CleanSolution(const std::vector<int> &_rows_to_keep,
                             const std::vector<int> &_cols_to_keep) : rows_to_keep(_rows_to_keep),
                                                                      cols_to_keep(_cols_to_keep) {}

CleanSolution::CleanSolution(const std::vector<bool> &_rows_to_keep,
                             const std::vector<bool> &_cols_to_keep) {
  rows_to_keep.resize(_rows_to_keep.size(), 0);
  cols_to_keep.resize(_cols_to_keep.size(), 0);

  for (std::size_t i = 0; i < rows_to_keep.size(); ++i) {
    if (_rows_to_keep[i]) {
      rows_to_keep[i] = 1;
    }
  }
  for (std::size_t j = 0; j < cols_to_keep.size(); ++j) {
    if (_cols_to_keep[j]) {
      cols_to_keep[j] = 1;
    }
  }
}

CleanSolution::~CleanSolution() {}

void CleanSolution::update(const std::vector<bool> &_rows_to_keep,
                           const std::vector<bool> &_cols_to_keep) {
  if (rows_to_keep.size() != _rows_to_keep.size()) {
    fprintf(stderr, "Error - Trying to update solution.\n");
    fprintf(stderr, "Old sol contains %lu rows and new sol contains %lu rows\n", rows_to_keep.size(), _rows_to_keep.size());
    exit(1);
  }
  if (cols_to_keep.size() != _cols_to_keep.size()) {
    fprintf(stderr, "Error - Trying to update solution.\n");
    fprintf(stderr, "Old sol contains %lu cols and new sol contains %lu cols\n", cols_to_keep.size(), _cols_to_keep.size());
    exit(1);
  }

  for (std::size_t i = 0; i < rows_to_keep.size(); ++i) {
    if (_rows_to_keep[i]) {
      rows_to_keep[i] = 1;
    }
  }
  for (std::size_t j = 0; j < cols_to_keep.size(); ++j) {
    if (_cols_to_keep[j]) {
      cols_to_keep[j] = 1;
    }
  }
}

void CleanSolution::update(const std::vector<int> &_rows_to_keep,
                           const std::vector<int> &_cols_to_keep) {
  if (rows_to_keep.size() != _rows_to_keep.size()) {
    fprintf(stderr, "Error - Trying to update solution.\n");
    fprintf(stderr, "Old sol contains %lu rows and new sol contains %lu rows\n", rows_to_keep.size(), _rows_to_keep.size());
    exit(1);
  }
  if (cols_to_keep.size() != _cols_to_keep.size()) {
    fprintf(stderr, "Error - Trying to update solution.\n");
    fprintf(stderr, "Old sol contains %lu cols and new sol contains %lu cols\n", cols_to_keep.size(), _cols_to_keep.size());
    exit(1);
  }

  rows_to_keep = _rows_to_keep;
  cols_to_keep = _cols_to_keep;
}

void CleanSolution::write_to_file(const std::string &file_name) {
  FILE *sol;

  if((sol = fopen(file_name.c_str(), "w+")) == nullptr) {
    fprintf(stderr, "Could not open file (%s)", file_name.c_str());
    exit(EXIT_FAILURE);
  }

  // Write rows_to_keep
  for (std::size_t i = 0; i < rows_to_keep.size() - 1; ++i) {
    fprintf(sol, rows_to_keep[i] ? "1\t" : "0\t");
  }
  fprintf(sol, rows_to_keep[rows_to_keep.size()-1] ? "1\n" : "0\n");
 
  // Write cols_to_keep
  for (std::size_t i = 0; i < cols_to_keep.size() - 1; ++i) {
    fprintf(sol, cols_to_keep[i] ? "1\t" : "0\t");
  }
  fprintf(sol, cols_to_keep[cols_to_keep.size()-1] ? "1\n" : "0\n");

  fclose(sol);
}

void CleanSolution::read_from_file(const std::string &file_name) {
  std::ifstream sol;
  std::string line, s;

  sol.open(file_name.c_str());
  if (!sol) {
    fprintf(stderr, "Could not open file (%s)", file_name.c_str());
    exit(EXIT_FAILURE);
  }

  std::getline(sol, line); // Read in first line  
  std::istringstream iss(line);
  std::vector<int> tmp_rows, tmp_cols;

  while (std::getline(iss, s, '\t')) {
    if (s.compare("1") == 0) {
      tmp_rows.push_back(1);
    } else if (s.compare("0") == 0) {
      tmp_rows.push_back(0);
    } else {
      fprintf(stderr, "ERROR - unkwown value in solutions file %s in %s\n", s.c_str(), file_name.c_str());
    }
  }
  if (tmp_rows.size() != rows_to_keep.size()) {
    fprintf(stderr, "ERROR - Read in %lu rows from solution, but expected %lu.\n", tmp_rows.size(), rows_to_keep.size());
  }

  std::getline(sol, line); // Read in first line 
  iss.clear();
  iss.str(line);
  while (std::getline(iss, s, '\t')) {
    if (s.compare("1") == 0) {
      tmp_cols.push_back(1);
    } else if (s.compare("0") == 0) {
      tmp_cols.push_back(0);
    } else {
      fprintf(stderr, "ERROR - unkwown value in solutions file %s in %s\n", s.c_str(), file_name.c_str());
    }
  }
  if (tmp_cols.size() != cols_to_keep.size()) {
    fprintf(stderr, "ERROR - Read in %lu cols from solution, but expected %lu.\n", tmp_cols.size(), cols_to_keep.size());
  }

  rows_to_keep = tmp_rows;
  cols_to_keep = tmp_cols;
}

std::vector<int> CleanSolution::get_rows_to_keep() const {
  return rows_to_keep;
}

std::vector<int> CleanSolution::get_cols_to_keep() const {
  return cols_to_keep;
}

std::size_t CleanSolution::get_num_rows_kept() const {
  std::size_t count = 0;
  for (auto r : rows_to_keep) {
    count += r;
  }
  return count;
}

std::size_t CleanSolution::get_num_cols_kept() const {
  std::size_t count = 0;
  for (auto c : cols_to_keep) {
    count += c;
  }
  return count;
}