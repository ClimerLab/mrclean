#include <fstream>
#include <sstream>
#include <algorithm>
#include <assert.h>

#include "DataContainer.h"

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
DataContainer::DataContainer(const std::string &file_name,
                             const std::string &_na_symbol,
                             const std::size_t num_header_rows,
                             const std::size_t num_header_cols ) :  na_symbol(_na_symbol) {
  read(file_name, num_header_rows, num_header_cols);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
DataContainer::~DataContainer() {}

//------------------------------------------------------------------------------
// Returns the number of header rows in the DataContainer object.
//------------------------------------------------------------------------------
std::size_t DataContainer::get_num_header_rows() const {
  return header_rows.size();
}

//------------------------------------------------------------------------------
// Returns the number of header columns in the DataContainer object.
//------------------------------------------------------------------------------
std::size_t DataContainer::get_num_header_cols() const {
  return header_cols.empty() ? 0 : header_cols[0].size();
}

//------------------------------------------------------------------------------
// Returns the number of data rows in the DataContaier object.
//------------------------------------------------------------------------------
std::size_t DataContainer::get_num_data_rows() const {
  return data.size();
}

//------------------------------------------------------------------------------
// Returns the number of data columns in the DataContainer object.
//------------------------------------------------------------------------------
std::size_t DataContainer::get_num_data_cols() const {
  return data.empty() ? 0 : data[0].size();
}

//------------------------------------------------------------------------------
// Returns the total number of data elements by taking the product of the 
// number of data rows and the number of data columns.
//------------------------------------------------------------------------------
std::size_t DataContainer::get_num_data() const {
  return get_num_data_rows() * get_num_data_cols();
}

//------------------------------------------------------------------------------
// Returns true if the (i,j) element in the data matrix matches the 'na_symbol'
// and false otherwise.
//------------------------------------------------------------------------------
bool DataContainer::is_data_na(const std::size_t i, const std::size_t j) const {
  assert(i < get_num_data_rows());
  assert(j < get_num_data_cols());
  return data[i][j].compare(na_symbol) == 0;
}

//------------------------------------------------------------------------------
// Counts the total number of valid elements in the data matrix. This equals
// the number of elements that do NOT match the na_symbol.
//------------------------------------------------------------------------------
std::size_t DataContainer::get_num_valid_data() const {
  std::size_t count = 0;
  for (std::size_t i = 0; i < get_num_data_rows(); ++i) {
    for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
      if (!is_data_na(i, j)) {
        ++count;
      }
    }
  }
  return count;
}

//------------------------------------------------------------------------------
// Reads in the data file given by 'file_name'. Throws an error if file does
// not exit. Expects elements in file to be seperated with a tab ('\t'). Number
// of header rows and columns defaults to 1 in the header file.
//------------------------------------------------------------------------------
void DataContainer::read(const std::string &file_name,
                         const std::size_t num_header_rows,
                         const std::size_t num_header_cols)
{
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
  const std::size_t num_cols = std::count(line.begin(), line.end(), '\t') + 1;

  const std::size_t num_data_rows = num_rows - num_header_rows;
  const std::size_t num_data_cols = num_cols - num_header_cols;


  // Allocate memory
  std::vector<std::string> header_rows_row(num_cols);
  for (std::size_t i = 0; i < num_header_rows; ++i) {
    header_rows.push_back(header_rows_row);
  }

  std::vector<std::string> header_cols_row(num_header_cols);
  for (std::size_t i = 0; i < num_data_rows; ++i) {
    header_cols.push_back(header_cols_row);
  }

  std::vector<std::string> data_row(num_data_cols);
  for (std::size_t i = 0; i < num_data_rows; ++i) {
    data.push_back(data_row);
  }

  // Read in data
  input.seekg(std::ios_base::beg);  // Go to beginning of file
  for (std::size_t i = 0; i < num_header_rows; ++i) {
    std::getline(input, line);

    std::istringstream iss(line);
    std::string token;
    for (std::size_t j = 0; j < num_cols; ++j) {
      std::getline(iss, token, '\t');
      header_rows[i][j] = token;
    }
  }

  for (std::size_t i = 0; i < num_data_rows; ++i) {
    std::getline(input, line);

    std::istringstream iss(line);
    std::string token;
    for (std::size_t j = 0; j < num_header_cols; ++j) {
      std::getline(iss, token, '\t');
      header_cols[i][j] = token;
    }
    for (std::size_t j = 0; j < num_data_cols; ++j) {
      std::getline(iss, token, '\t');
      data[i][j] = token;
    }
  }

  input.close();
}

//------------------------------------------------------------------------------
// Writes out the header rows, header columns, and data, based on the
// 'rows_to_keep' and 'cols_to_keep' vectors. For header rows, the header is
// printed if the 'rows_to_keep' element is true. For header cols, the header is
// printed if 'cols_to_keep' is true. The data elements are printed if the
// corresponding 'rows_to_keep' and 'cols_to_keep' are both true.
//------------------------------------------------------------------------------
void DataContainer::write(const std::string &file_name,
                          const std::vector<bool> &rows_to_keep,
                          const std::vector<bool> &cols_to_keep)
{
  assert(header_rows.size() > 0);
  assert(header_cols.size() > 0);
  assert(data.size() > 0);
  assert(data.size() == rows_to_keep.size());
  assert(data[0].size() == cols_to_keep.size());
  assert(header_cols.size() == data.size());
  assert(header_rows[0].size() == header_cols[0].size() + data[0].size());
  
  const std::size_t num_header_rows = get_num_header_rows();
  const std::size_t num_header_cols = get_num_header_cols();
  const std::size_t num_data_rows = get_num_data_rows();
  const std::size_t num_data_cols = get_num_data_cols();

  FILE *output;
  if ((output = fopen(file_name.c_str(), "w")) == nullptr) {
    fprintf(stderr, "ERROR - Could not open file (%s).\n", file_name.c_str());
    exit(EXIT_FAILURE);
  }
 
  // Write header rows
  for (std::size_t i = 0; i < num_header_rows; ++i) {
    fprintf(output, "%s", header_rows[i][0].c_str());

    for (std::size_t j = 1; j < num_header_cols; ++j) {
      fprintf(output, "\t%s", header_rows[i][j].c_str());
    }

    for (std::size_t j = 0; j < num_data_cols; ++j) {
      if (cols_to_keep[j]) {
        fprintf(output, "\t%s", header_rows[i][j + num_header_cols].c_str());
      }
    }
    fprintf(output, "\n");
  }

  // Write header cols and data
  for (std::size_t i = 0; i < num_data_rows; ++i) {
    if (rows_to_keep[i]) {
      fprintf(output, "%s", header_cols[i][0].c_str());

      for (std::size_t j = 1; j < num_header_cols; ++j) {
        fprintf(output, "\t%s", header_cols[i][j].c_str());
      }

      for (std::size_t j = 0; j < num_data_cols; ++j) {
        if (cols_to_keep[j]) {
          fprintf(output, "\t%s", data[i][j].c_str());
        }
      }
      fprintf(output, "\n");
    }
  }

  fclose(output);
}

//------------------------------------------------------------------------------
// Returns the total number of data elements kept. An element is only kept if
// the corresponding 'keep_row' and 'keep_col' are both true. Inpute vectors
// are boolean
//------------------------------------------------------------------------------
std::size_t DataContainer::get_num_valid_data_kept(const std::vector<bool> &keep_row,
                                                   const std::vector<bool> &keep_col) const
{
  assert(keep_row.size() == get_num_data_rows());
  assert(keep_col.size() == get_num_data_cols());

  std::size_t count = 0;
  for (std::size_t i = 0; i < get_num_data_rows(); ++i) {
    for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
      if (keep_row[i] && keep_col[j] && !is_data_na(i, j)) {
        ++count;
      }
    }
  }
  return count; 
}

//------------------------------------------------------------------------------
// Returns the total number of data elements kept. An element is only kept if
// the corresponding 'keep_row' and 'keep_col' are both true. Inpute vectors
// are integers
//------------------------------------------------------------------------------
std::size_t DataContainer::get_num_valid_data_kept(const std::vector<int> &keep_row,
                                                   const std::vector<int> &keep_col) const
{
  assert(keep_row.size() == get_num_data_rows());
  assert(keep_col.size() == get_num_data_cols());

  std::size_t count = 0;
  for (std::size_t i = 0; i < get_num_data_rows(); ++i) {
    for (std::size_t j = 0; j < get_num_data_cols(); ++j) {
      if ((keep_row[i] == 1) && (keep_col[j] == 1) && !is_data_na(i, j)) {
        ++count;
      }
    }
  }
  return count; 
}

//------------------------------------------------------------------------------
// Returns the maximum percent of missing data contained in any data row.
//------------------------------------------------------------------------------
double DataContainer::get_max_perc_miss_row() const
{
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

//------------------------------------------------------------------------------
// Returns the maximum percent of missing data contained in any data column.
//------------------------------------------------------------------------------
double DataContainer::get_max_perc_miss_col() const
{
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

//------------------------------------------------------------------------------
// Returns the minimum percent of missing data contained in any data row.
//------------------------------------------------------------------------------
double DataContainer::get_min_perc_miss_row() const
{
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

//------------------------------------------------------------------------------
// Returns the minimum percent of missing data contained in any data column.
//------------------------------------------------------------------------------
double DataContainer::get_min_perc_miss_col() const
{
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

//------------------------------------------------------------------------------
// Prints out a binary representation of the data matrix. The header rows and 
// columns are printed like the original data matrix. But the data elements are
// printed as 0 if they match the 'na_symbol' and 1 otherwise.
//------------------------------------------------------------------------------
void DataContainer::print_binary(const std::string &filename) const {
  assert(header_rows.size() > 0);
  assert(header_cols.size() > 0);
  assert(data.size() > 0);
  assert(header_cols.size() == data.size());
  assert(header_rows[0].size() == header_cols[0].size() + data[0].size());
  
  const std::size_t num_header_rows = get_num_header_rows();
  const std::size_t num_header_cols = get_num_header_cols();
  const std::size_t num_data_rows = get_num_data_rows();
  const std::size_t num_data_cols = get_num_data_cols();
  
  FILE *output;
  if ((output = fopen(filename.c_str(), "w")) == nullptr) {
    fprintf(stderr, "ERROR - Could not open file (%s).\n", filename.c_str());
    exit(EXIT_FAILURE);
  }

  // Write header rows
  for (std::size_t i = 0; i < num_header_rows; ++i) {
    fprintf(output, "%s", header_rows[i][0].c_str());

    for (std::size_t j = 1; j < num_header_cols + num_data_cols; ++j) {
      fprintf(output, "\t%s", header_rows[i][j].c_str());
    }

    fprintf(output, "\n");
  }

  // Write header cols and data
  for (std::size_t i = 0; i < num_data_rows; ++i) {
    fprintf(output, "%s", header_cols[i][0].c_str());

    for (std::size_t j = 1; j < num_header_cols; ++j) {
      fprintf(output, "\t%s", header_rows[i][j].c_str());
    }

    for (std::size_t j = 0; j < num_data_cols; ++j) {
      fprintf(output, is_data_na(i,j) ? "\t0" : "\t1");
    }
    fprintf(output, "\n");
  }

  fclose(output);
}