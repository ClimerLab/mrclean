#ifndef ELEMENT_SOLVER_H
#define ELEMENT_SOLVER_H

#include <vector>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilomodel.h>

#include "BinContainer.h"

class ElementSolver {
  private:
    const BinContainer *data;
    const double TOL;
    const double max_perc_missing;

    std::size_t num_rows;
    std::size_t num_cols;
    std::size_t num_data;
    std::vector<double> x_var;
    std::vector<double> r_var;
    std::vector<double> c_var;
    
    double obj_value;

    bool use_incumbent_obj;
    double limit;

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
    ElementSolver(const BinContainer &_data,
                  const double _min_perc_missing,
                  const double _TOL = 0.000001);
    ~ElementSolver();

    void solve();
    double get_obj_value() const;
    std::vector<bool> get_rows_to_keep() const;
    std::vector<bool> get_cols_to_keep() const;
    std::size_t get_num_rows_to_keep() const;
    std::size_t get_num_cols_to_keep() const;

    void set_incumbent_obj(const double _obj);
};

#endif