#ifndef ROW_COL_SOLVER_H
#define ROW_COL_SOLVER_H

#include <vector>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilomodel.h>

#include "DataContainer.h"

class RowColSolver
{
private:
  const DataContainer *data;
  const double TOL;
  const double max_perc_missing;

  std::size_t num_rows;
  std::size_t num_cols;
  std::size_t num_data;

  std::vector<double> r_var;
  std::vector<double> c_var;
  double obj_value;
  bool use_incumbent;

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
  RowColSolver(const DataContainer &_data,
               const double _min_perc_missing,
               const double _TOL = 0.00001);
  ~RowColSolver();

  void solve();
  double get_obj_value() const;
  std::vector<bool> get_rows_to_keep() const;
  std::vector<bool> get_cols_to_keep() const;
  std::size_t get_num_rows_to_keep() const;
  std::size_t get_num_cols_to_keep() const;

  void set_incumbent(const std::vector<bool> &keep_rows, const std::vector<bool> &keep_cols);
  void set_incumbent(const std::vector<int> &keep_rows, const std::vector<int> &keep_cols);
};

#endif