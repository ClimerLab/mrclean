#ifndef GREEDY_SOLVER_H
#define GREEDY_SOLVER_H

#include <vector>
#include "DataContainer.h"

class GreedySolver
{
private:
  const DataContainer *data;
  const std::size_t num_rows;
  const std::size_t num_cols;
  const double max_perc_miss;
  
  std::vector<std::size_t> alphas;
  std::vector<std::size_t> betas;
  std::vector<bool> keep_row;
  std::vector<bool> keep_col;
  std::size_t num_rows_kept;
  std::size_t num_cols_kept;
  
  void calc_alphas();
  void calc_betas();
  void remove_row(const std::size_t idx);
  void remove_col(const std::size_t idx);
  void update_rows(const std::size_t removed_col);
  void update_cols(const std::size_t removed_row);
  
  std::size_t get_num_missing_row(const std::size_t idx) const;
  std::size_t get_num_missing_col(const std::size_t idx) const;
  double get_perc_miss_row(const std::size_t idx) const;
  double get_perc_miss_col(const std::size_t idx) const;

  std::vector<std::size_t> get_missing_cols(const std::size_t rowIdx) const;
  std::vector<std::size_t> get_missing_rows(const std::size_t colIdx) const;

  std::size_t calc_num_rows_to_remove(const std::size_t idx) const;
  std::size_t calc_num_cols_to_remove(const std::size_t idx) const;
  
public:
  GreedySolver(const DataContainer &_data,
               const double max_perc_miss);
  ~GreedySolver();

  void solve();

  std::vector<bool> get_rows_kept_as_bool() const;
  std::vector<bool> get_cols_kept_as_bool() const; 
  std::size_t get_num_rows_kept() const;
  std::size_t get_num_cols_kept() const;
};

#endif