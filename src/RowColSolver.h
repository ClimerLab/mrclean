#ifndef ROW_COL_SOLVER_H
#define ROW_COL_SOLVER_H

#include <vector>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilomodel.h>

#include "BinContainer.h"

class RowColSolver
{
private:
  const BinContainer *data;
  const double TOL;
  const double max_perc_missing;
  const std::size_t row_lb;
  const std::size_t col_lb;

  std::size_t num_rows;
  std::size_t num_cols;
  std::size_t num_data;

  std::vector<double> r_var;
  std::vector<double> c_var;
  double obj_value;

  bool use_incumbent_obj;
  double limit;

  IloEnv env;
  IloCplex cplex;
  IloModel model;
  IloNumVarArray r;
  IloNumVarArray c;
  IloNumArray r_copy;
  IloNumArray c_copy;
  IloExpr obj;
  
  void build_model();
  void round_extreme_values();

public:
  RowColSolver(const BinContainer &_data,
               const double _min_perc_missing,
               const std::size_t _row_lb,
               const std::size_t _col_lb,
               const double _TOL = 0.00001);
  ~RowColSolver();

  void solve();
  double get_obj_value() const;
  std::vector<bool> get_rows_to_keep() const;
  std::vector<bool> get_cols_to_keep() const;
  std::size_t get_num_rows_to_keep() const;
  std::size_t get_num_cols_to_keep() const;

  void set_incumbent_obj(const double _obj);
};

#endif