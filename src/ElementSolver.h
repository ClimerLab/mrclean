#ifndef ELEMENT_SOLVER_H
#define ELEMENT_SOLVER_H

#include <vector>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilomodel.h>

#include "DataContainer.h"

class ElementSolver
{
private:
  const DataContainer *data;
  const double TOL;
  const double max_perc_missing;

  std::size_t num_rows;
  std::size_t num_cols;
  std::size_t num_data;
  std::vector<double> x_var;
  std::vector<double> r_var;
  std::vector<double> c_var;
  double obj_value;
  bool use_incumbent;

  IloEnv env;
  IloCplex cplex;
  IloModel model;
  IloNumVarArray x;
  IloNumVarArray r;
  IloNumVarArray c;
  IloNumArray x_copy;
  IloNumArray r_copy;
  IloNumArray c_copy;
  IloExpr obj;

  std::size_t get_data_idx(const std::size_t i, const std::size_t j) const;
  void build_model();
  void round_extreme_values();

public:
  ElementSolver(const DataContainer &_data,
                const double _min_perc_missing,
                const double _TOL = 0.000001);
  ~ElementSolver();

  void solve();
  double get_obj_value() const;
  std::vector<bool> get_rows_to_keep() const;
  std::vector<bool> get_cols_to_keep() const;
  std::size_t get_num_rows_to_keep() const;
  std::size_t get_num_cols_to_keep() const;

  void set_incumbent(const std::vector<bool> &keep_rows, const std::vector<bool> &keep_cols);
  void set_incumbent(const std::vector<int> &keep_rows, const std::vector<int> &keep_cols);
};
// class ElementSolver
// {
// private:
//   const DataContainer *data;
//   const std::string na_symbol;
//   const double TOL;
//   const double max_perc_missing;

//   std::size_t num_rows;
//   std::size_t num_cols;
//   std::size_t num_data;
//   std::vector<double> variables;
//   // std::vector<double> dVar;
//   // std::vector<double> rVar;
//   // std::vector<double> cVar;
//   double obj_value;
//   bool use_incumbent;

//   IloEnv env;
//   IloCplex cplex;
//   IloModel model;
//   IloNumVarArray x;
//   IloNumArray xCopy;
//   // IloNumVarArray d;
//   // IloNumVarArray r;
//   // IloNumVarArray c;
//   // IloNumArray dCopy;
//   // IloNumArray rCopy;
//   // IloNumArray cCopy;
//   IloExpr obj;

//   std::size_t get_data_idx(const std::size_t i, const std::size_t j) const;
//   void build_model();
//   void round_extreme_values();  

// public:
//   ElementSolver(const DataContainer &_data,
//             const std::string &_na_symbol,
//             const double _min_perc_missing,
//             const double _TOL = 0.000001);
//   ~ElementSolver();

//   void solve();
//   double get_obj_value() const;
//   std::vector<bool> get_rows_to_keep() const;
//   std::vector<bool> get_cols_to_keep() const;
//   std::size_t get_num_rows_to_keep() const;
//   std::size_t get_num_cols_to_keep() const;

//   void set_incumbent(const std::vector<bool> &keep_rows, const std::vector<bool> &keep_cols);
// };

#endif