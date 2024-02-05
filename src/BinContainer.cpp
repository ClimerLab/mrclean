#include "BinContainer.h"
#include <fstream>
#include <sstream>
#include <algorithm>

BinContainer::BinContainer(const std::string &_file_name,
                           const std::string &_na_symbol,
                           const std::size_t _num_header_rows,
                           const std::size_t _num_header_cols) :  file_name(_file_name),
                                                                  na_symbol(_na_symbol),
                                                                  num_header_rows(_num_header_rows),
                                                                  num_header_cols(_num_header_cols) {
  read();
  calc_num_valid();
}

BinContainer::~BinContainer() {}

void BinContainer::read() {
  std::string line;
  std::ifstream input;

  input.open(file_name.c_str());
  if (!input)
    throw std::runtime_error("Input file could not be opened.");
  
  // Determine the number of rows
  const std::size_t num_rows = std::count(std::istreambuf_iterator<char>(input),
                                          std::istreambuf_iterator<char>(), '\n');

  // Determine the number of columns
  input.seekg(std::ios_base::beg);  // Go to beginning of file
  std::getline(input, line); // Read in first line
  // const std::size_t num_cols = std::count(line.begin(), line.end(), '\t') + 1;
  std::size_t num_cols = 0;
    std::istringstream iss(line);
    std::string token;
    while (!iss.eof()) {
      std::getline(iss, token, '\t');
      if (!token.empty()) {
        ++num_cols;
      }
    }

  const std::size_t num_data_rows = num_rows - num_header_rows;
  const std::size_t num_data_cols = num_cols - num_header_cols;

  fprintf(stderr, "Num rows: %lu\n", num_rows);
  fprintf(stderr, "Num cols: %lu\n", num_cols);

  // Allocate memory
  std::vector<bool> data_row(num_data_cols);
  for (std::size_t i = 0; i < num_data_rows; ++i) {
    data.push_back(data_row);
  }
  fprintf(stderr, "Allocated memory\n");

  // Read in data
  input.seekg(std::ios_base::beg);  // Go to beginning of file
  for (std::size_t i = 0; i < num_header_rows; ++i) {
    std::getline(input, line);
  }

  for (std::size_t i = 0; i < num_data_rows; ++i) {
    std::getline(input, line);

    std::istringstream iss(line);
    std::string token;
    for (std::size_t j = 0; j < num_header_cols; ++j) {
      std::getline(iss, token, '\t');
    }
    for (std::size_t j = 0; j < num_data_cols; ++j) {
      std::getline(iss, token, '\t');
      token = trim(token);

      if (token.compare(na_symbol) == 0) {
        data[i][j] = false;
      } else {
        data[i][j] = true;
      }
    }
  }

  input.close();
}

void BinContainer::calc_num_valid() {
  num_valid_rows.resize(get_num_data_rows(), 0);
  num_valid_cols.resize(get_num_data_cols(), 0);

  for (std::size_t i = 0; i < get_num_data_rows(); ++i) {
    for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
      if (!is_data_na(i,j)) {
        ++num_valid_rows[i];
        ++num_valid_cols[j];
      }
    }
  }
}

std::string BinContainer::trim(std::string &str) const {
  size_t first = str.find_first_not_of(' ');
  if (std::string::npos == first) {
    return str;
  }
  size_t last = str.find_last_not_of(' ');
  return str.substr(first, (last - first + 1));
}

std::size_t BinContainer::get_num_header_rows() const {
  return num_header_rows;
}

std::size_t BinContainer::get_num_header_cols() const {
  return num_header_cols;
}

std::size_t BinContainer::get_num_data_rows() const {
  return data.size();
}

std::size_t BinContainer::get_num_data_cols() const {
  return data.empty() ? 0 : data[0].size();
}

std::size_t BinContainer::get_num_data() const {
  return get_num_data_rows() * get_num_data_cols();
}

std::size_t BinContainer::get_num_valid_data() const {
  std::size_t count = 0;
  for (std::size_t i = 0; i < get_num_data_rows(); ++i) {
    for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
      if (data[i][j]) {
        ++count;
      }
    }
  }
  return count;
}

std::size_t BinContainer::get_num_valid_data_kept(const std::vector<bool> &keep_row,
                                                  const std::vector<bool> &keep_col) const {
  if (keep_row.size() != get_num_data_rows()) {
    fprintf(stderr, "ERROR - BinContainer::get_num_valid_data_kept - The number of elements in 'keep_row' ");
    fprintf(stderr, "does not match the number of data rows\n (%lu vs. %lu)\n", keep_row.size(), get_num_data_rows());
    exit(EXIT_FAILURE);
  }
  if (keep_col.size() != get_num_data_cols()) {
    fprintf(stderr, "ERROR - BinContainer::get_num_valid_data_kept - The number of elements in 'keep_col' ");
    fprintf(stderr, "does not match the number of data columns\n (%lu vs. %lu)\n", keep_col.size(), get_num_data_cols());
    exit(EXIT_FAILURE);
  }

  std::size_t count = 0;
  for (std::size_t i = 0; i < get_num_data_rows(); ++i) {
    if (!keep_row[i]) { continue; }

    for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
      if (keep_col[j] && !is_data_na(i, j)) {
        ++count;
      }
    }
  }
  return count; 
}

std::size_t BinContainer::get_num_valid_data_kept(const std::vector<int> &keep_row,
                                                  const std::vector<int> &keep_col) const {
  if (keep_row.size() != get_num_data_rows()) {
    fprintf(stderr, "ERROR - BinContainer::get_num_valid_data_kept - The number of elements in 'keep_row' ");
    fprintf(stderr, "does not match the number of data rows\n (%lu vs. %lu)\n", keep_row.size(), get_num_data_rows());
    exit(EXIT_FAILURE);
  }
  if (keep_col.size() != get_num_data_cols()) {
    fprintf(stderr, "ERROR - BinContainer::get_num_valid_data_kept - The number of elements in 'keep_col' ");
    fprintf(stderr, "does not match the number of data columns\n (%lu vs. %lu)\n", keep_col.size(), get_num_data_cols());
    exit(EXIT_FAILURE);
  }

  std::size_t count = 0;
  for (std::size_t i = 0; i < get_num_data_rows(); ++i) {
    if (keep_row[i] == 0) { continue; }

    for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
      if ((keep_col[j]==1) && !is_data_na(i, j)) {
        ++count;
      }
    }
  }
  return count; 
}

std::size_t BinContainer::get_num_invalid_in_row(const std::size_t row) const {
  if (row >= get_num_data_rows()) {
    fprintf(stderr, "ERROR - BinContainer::get_num_invalid_in_row - Trying to access index out of bounds\n");
    exit(EXIT_FAILURE);
  }
  return get_num_data_rows() - num_valid_rows[row];
}

std::size_t BinContainer::get_num_invalid_in_col(const std::size_t col) const {
  if (col >= get_num_data_cols()) {
    fprintf(stderr, "ERROR - BinContainer::get_num_invalid_in_col - Trying to access index out of bounds\n");
    exit(EXIT_FAILURE);
  }
  return get_num_data_cols() - num_valid_cols[col];
}

std::size_t BinContainer::get_num_valid_in_row(const std::size_t row) const {
  if (row >= get_num_data_rows()) {
    fprintf(stderr, "ERROR - BinContainer::get_num_valid_in_row - Trying to access index out of bounds\n");
    exit(EXIT_FAILURE);
  }
  return num_valid_rows[row];
}

std::size_t BinContainer::get_num_valid_in_col(const std::size_t col) const {
  if (col >= get_num_data_cols()) {
    fprintf(stderr, "ERROR - BinContainer::get_num_valid_in_col - Trying to access index out of bounds\n");
    exit(EXIT_FAILURE);
  }
  return num_valid_cols[col];
}

double BinContainer::get_max_perc_miss_row() const {
  double max_perc_miss = 0.0;
  for (std::size_t i = 0; i < get_num_data_rows(); ++i) {
    std::size_t count_miss = 0;
    for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
      if (is_data_na(i,j)) {
        ++count_miss;
      }
    }
    double perc_miss = static_cast<double>(count_miss) / get_num_data_cols();
    if (perc_miss > max_perc_miss) {
      max_perc_miss = perc_miss;
    }
  }
  return max_perc_miss;
}

double BinContainer::get_max_perc_miss_col() const {
  double max_perc_miss = 0.0;
  for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
    std::size_t count_miss = 0;
    for (std::size_t i = 0; i < get_num_data_rows(); ++i) {    
      if (is_data_na(i,j)) {
        ++count_miss;
      }
    }
    double perc_miss = static_cast<double>(count_miss) / get_num_data_rows();
    if (perc_miss > max_perc_miss) {
      max_perc_miss = perc_miss;
    }
  }
  return max_perc_miss;
}

double BinContainer::get_min_perc_miss_row() const {
  double min_perc_miss = 1.0;
  for (std::size_t i = 0; i < get_num_data_rows(); ++i) {
    std::size_t count_miss = 0;
    for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
      if (is_data_na(i,j)) {
        ++count_miss;
      }
    }
    double perc_miss = static_cast<double>(count_miss) / get_num_data_cols();
    if (perc_miss < min_perc_miss) {
      min_perc_miss = perc_miss;
    }
  }
  return min_perc_miss;
}

double BinContainer::get_min_perc_miss_col() const {
  double min_perc_miss = 1.0;
  for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
    std::size_t count_miss = 0;
    for (std::size_t i = 0; i < get_num_data_rows(); ++i) {    
      if (is_data_na(i,j)) {
        ++count_miss;
      }
    }
    double perc_miss = static_cast<double>(count_miss) / get_num_data_rows();
    if (perc_miss < min_perc_miss) {
      min_perc_miss = perc_miss;
    }
  }
  return min_perc_miss;
}

bool BinContainer::is_data_na(const std::size_t i, const std::size_t j) const {
  return !data[i][j];
}

void BinContainer::write_orig(const std::string &out_file,
                              const std::vector<bool> &rows_to_keep,
                              const std::vector<bool> &cols_to_keep) const {
  if (get_num_data_rows() != rows_to_keep.size()) {
    fprintf(stderr, "ERROR - BinContainer::write_orig - Size of 'rows_to_keep' does not match the number of data rows\n");
    exit(EXIT_FAILURE);
  }
  if (get_num_data_cols() != cols_to_keep.size()) {
    fprintf(stderr, "ERROR - BinContainer::write_orig - Size of 'cols_to_keep' does not match the number of data cols\n");
    exit(EXIT_FAILURE);
  }

  const std::size_t num_header_rows = get_num_header_rows();
  const std::size_t num_header_cols = get_num_header_cols();
  const std::size_t num_data_rows = get_num_data_rows();
  const std::size_t num_data_cols = get_num_data_cols();

  fprintf(stderr, "Num header rows: %lu\n", num_header_rows);
  fprintf(stderr, "Num header cols: %lu\n", num_header_cols);
  fprintf(stderr, "Num data rows: %lu\n", num_data_rows);
  fprintf(stderr, "Num data cols: %lu\n", num_data_cols);

  std::string line;
  std::ifstream input;

  input.open(file_name.c_str());
  if (!input) {
    fprintf(stderr, "Input file could not be opened.\n");
    exit(EXIT_FAILURE);
  }

  FILE *output;
  if ((output = fopen(out_file.c_str(), "w")) == nullptr) {
    fprintf(stderr, "ERROR - Could not open file (%s).\n", out_file.c_str());
    exit(EXIT_FAILURE);
  }

  // Read in data
  for (std::size_t i = 0; i < num_header_rows; ++i) {
    std::getline(input, line);

    std::istringstream iss(line);
    std::string token;

    // Print first header col
    std::getline(iss, token, '\t');
    token = trim(token);
    fprintf(output, "%s", token.c_str());
    
    // Print remaining header columns
    for (std::size_t j = 1; j < num_header_cols; ++j) {
      std::getline(iss, token, '\t');
      token = trim(token);
      fprintf(output, "\t%s", token.c_str());
    }

    // Print data cols if kept
    for (std::size_t j = 0; j < num_data_cols; ++j) {
      std::getline(iss, token, '\t');
      if (cols_to_keep[j]) {
        token = trim(token);
        fprintf(output, "\t%s", token.c_str());
      }
    }
    fprintf(output, "\n");
  }

  for (std::size_t i = 0; i < num_data_rows; ++i) {
    std::getline(input, line);

    if (rows_to_keep[i]) {
      std::istringstream iss(line);
      std::string token;

      // Print first header col
      std::getline(iss, token, '\t');
      token = trim(token);
      fprintf(output, "%s", token.c_str());

      // Print remaining header columns
      for (std::size_t j = 1; j < num_header_cols; ++j) {
        std::getline(iss, token, '\t');
        token = trim(token);
        fprintf(output, "\t%s", token.c_str());
      }

      // Print data cols if kept
      for (std::size_t j = 0; j < num_data_cols; ++j) {
        std::getline(iss, token, '\t');
        if (cols_to_keep[j]) {
          token = trim(token);
          fprintf(output, "\t%s", token.c_str());
        }
      }
      fprintf(output, "\n");
    }
  }

  fclose(output);
  input.close();
}

void BinContainer::write_orig(const std::string &out_file,
                              const std::vector<int> &rows_to_keep,
                              const std::vector<int> &cols_to_keep) const {
  std::vector<bool> rows_to_keep_bool(rows_to_keep.size(), false);
  std::vector<bool> cols_to_keep_bool(cols_to_keep.size(), false);

  for (std::size_t i = 0; i < rows_to_keep.size(); ++i) {
    if (rows_to_keep[i] == 1) {
      rows_to_keep_bool[i] = true;
    }
  }
  for (std::size_t j = 0; j < cols_to_keep.size(); ++j) {
    if (cols_to_keep[j] == 1) {
      cols_to_keep_bool[j] = true;
    }
  }

  write_orig(out_file, rows_to_keep_bool, cols_to_keep_bool);
}

void BinContainer::write_orig_transpose(const std::string &_file_name) const {

}

void BinContainer::print_stats() const {
  const std::size_t M = get_num_data_rows();
  const std::size_t N = get_num_data_cols();

  std::vector<double> perc_miss_row(M, 0.0);
  std::vector<double> perc_miss_col(N, 0.0);
  double total_perc_miss = 0.0;

  for (std::size_t i = 0; i < M; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      if (!data[i][j]) {
        ++perc_miss_row[i];
        ++total_perc_miss;
      }
    }
    perc_miss_row[i] /= N;
  }

  for (std::size_t j = 0; j < N; ++j) {
    for (std::size_t i = 0; i < M; ++i) {    
      if (!data[i][j]) {
        ++perc_miss_col[j];
      }
    }
    perc_miss_col[j] /= M;
  }
  
  double min_row = 1.0, min_col = 1.0;
  double max_row = 0.0, max_col = 0.0;

  for (auto r : perc_miss_row) {
    if (r < min_row) {
      min_row = r;
    }

    if (r > max_row) {
      max_row = r;
    }
  }
  
  for (auto c : perc_miss_col) {
    if (c < min_col) {
      min_col = c;
    }

    if (c > max_col) {
      max_col = c;
    }
  }
  
  std::size_t total_elements = M * N;

  fprintf(stderr, "Stats:\n");
  fprintf(stderr, "Total perc missoing: %lf\n", total_perc_miss / total_elements *100);
  fprintf(stderr, "Min perc missing row: %lf\n", min_row * 100);
  fprintf(stderr, "Max perc missing row: %lf\n", max_row * 100);
  fprintf(stderr, "Min perc missing col: %lf\n", min_col * 100);
  fprintf(stderr, "Max perc missing col: %lf\n", max_col * 100);
}