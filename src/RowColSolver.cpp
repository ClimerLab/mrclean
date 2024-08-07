#include <assert.h>
#include "RowColSolver.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
RowColSolver::RowColSolver(const BinContainer &_data,
                           const double _min_perc_missing,
                           const std::size_t _row_lb,
                           const std::size_t _col_lb,
                           const double _TOL) : data(&_data),
                                                TOL(_TOL),
                                                max_perc_missing(_min_perc_missing),
                                                row_lb(_row_lb),
                                                col_lb(_col_lb),
                                                num_rows(data->get_num_data_rows()),
                                                num_cols(data->get_num_data_cols()),
                                                num_data(num_rows * num_cols),
                                                r_var(num_rows),
                                                c_var(num_cols),
                                                obj_value(0),
                                                use_incumbent_obj(false),
                                                limit(0.0),
                                                env(IloEnv()),
                                                cplex(IloCplex(env)),
                                                model(IloModel(env)),
                                                r(IloNumVarArray(env, num_rows, 0, 1, ILOINT)),
                                                c(IloNumVarArray(env, num_cols, 0, 1, ILOINT)),
                                                r_copy(IloNumArray(env, num_rows)),
                                                c_copy(IloNumArray(env, num_cols)),
                                                obj(IloExpr(env))
{
  build_model();
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
RowColSolver::~RowColSolver() {}


//------------------------------------------------------------------------------
// Builds the CPLEX model based on the BinContainer object.
//------------------------------------------------------------------------------
void RowColSolver::build_model() {
  // Add objective function
  for (std::size_t i = 0; i < num_rows; ++i) {
    double alpha = 0.0;
    for (std::size_t j = 0; j < num_cols; ++j) {
      if (!data->is_data_na(i,j)) {
        ++alpha;
      }
    }
    
    obj += static_cast<IloNum>(alpha) * r[i];
  }
  for (std::size_t j = 0; j < num_cols; ++j) {
    double beta = 0.0;
    for (std::size_t i = 0; i < num_rows; ++i) {
      if (!data->is_data_na(i,j)) {
        ++beta;
      }
    }
    obj += static_cast<IloNum>(beta) * c[j];
  }
  model.add(IloMaximize(env, obj, "Objective"));


  // Add contraints on rows
  for (std::size_t i = 0; i < num_rows; ++i) {    
    IloExpr col_sum(env);
    for (std::size_t j = 0; j < num_cols; ++j) {
      double b = data->is_data_na(i,j) ? 0.0 : 1.0;
      // col_sum -= x[num_rows + j] * static_cast<IloNum>((max_perc_missing + b - 1.0) / num_cols);
      col_sum -= c[j] * static_cast<IloNum>((max_perc_missing + b - 1.0) / num_cols);
    }
    model.add(r[i] + col_sum <= 1);
    col_sum.end();
  }
  
  // Add constrains on cols
  for (std::size_t j = 0; j < num_cols; ++j) {    
    IloExpr row_sum(env);
    for (std::size_t i = 0; i < num_rows; ++i) {
      double b = data->is_data_na(i,j) ? 0.0 : 1.0;
      // row_sum -= x[i] * static_cast<IloNum>((max_perc_missing + b - 1.0) / num_rows);
      row_sum -= r[i] * static_cast<IloNum>((max_perc_missing + b - 1.0) / num_rows);
    }
    // model.add(x[num_rows + j] + row_sum <= 1);
    model.add(c[j] + row_sum <= 1);
    row_sum.end();
  }

  // Add constraint for row_lb and col_lb
  IloExpr row_lb_sum(env);
  for (std::size_t i = 0; i < num_rows; ++i) {
    row_lb_sum += r[i];
  }
  model.add(row_lb_sum >= static_cast<double>(row_lb));
  row_lb_sum.end();

  IloExpr col_lb_sum(env);
  for (std::size_t j = 0; j < num_cols; ++j) {
    col_lb_sum += c[j];
  }
  model.add(col_lb_sum >= static_cast<double>(col_lb));
  col_lb_sum.end();
}

//------------------------------------------------------------------------------
// Rounds values close to 1 or 0 to 1 or 0 respectively. Values must be within
// 'TOL' of 1 or 0 to be rounded. Row decision variables and column decision 
// variables are both rounded.
//------------------------------------------------------------------------------
void RowColSolver::round_extreme_values() {
  for (std::size_t i = 0; i < r_var.size(); ++i) {
    if (r_var[i] > 1.0 - TOL) {
      r_var[i] = 1.0;
    } else if (r_var[i] < TOL) {
      r_var[i] = 0.0;
    }
  }
  for (std::size_t i = 0; i < c_var.size(); ++i) {
    if (c_var[i] > 1.0 - TOL) {
      c_var[i] = 1.0;
    } else if (c_var[i] < TOL) {
      c_var[i] = 0.0;
    }
  }
}


//------------------------------------------------------------------------------
// Calls CPLEX solver and rounds decision variables if a valid solution is
// found.
//------------------------------------------------------------------------------
void RowColSolver::solve() {
  cplex.extract(model);

  if (use_incumbent_obj) {
    cplex.setParam(IloCplex::Param::MIP::Tolerances::LowerCutoff, limit);
  }
  
  // cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, .04);
  cplex.setParam(IloCplex::Param::TimeLimit, 18000);
  cplex.setParam(IloCplex::Param::RandomSeed, 0);
  cplex.setParam(IloCplex::Param::Threads, 1);
  cplex.setOut(env.getNullStream());

  obj_value = 0.0;
  try {
    cplex.solve();

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
      printf("Infeasible Solution\n");
    } else if (cplex.getStatus() == IloAlgorithm::Optimal || cplex.getStatus() == IloAlgorithm::Feasible) {
      obj_value = cplex.getObjValue(); // Save objective value
      
      // Copy decision variables
      cplex.getValues(r, r_copy);
      cplex.getValues(c, c_copy);      
      for (std::size_t i = 0; i < r_var.size(); ++i) {
        r_var[i] = r_copy[i];
      }
      for (std::size_t i = 0; i < c_var.size(); ++i) {
        c_var[i] = c_copy[i];
      }
      
      round_extreme_values();
    }
  } catch (IloException& e) {
    std::cerr << "Concert exception caught: " << e << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception caught" << std::endl;
  }
}

//------------------------------------------------------------------------------
// Returns the objective value from CPLEX. If CPLEX has not been called, or did
// not find a solution, then 0 is returned.
//------------------------------------------------------------------------------
double RowColSolver::get_obj_value() const {
  return obj_value;
}

//------------------------------------------------------------------------------
// Returns the rows to keep as a boolean vector. In order for a row to be kept
// the corresponding decision variable must equal 1.0.
//------------------------------------------------------------------------------
std::vector<bool> RowColSolver::get_rows_to_keep() const {
  std::vector<bool> rows_to_keep(num_rows, false);

  for (std::size_t i = 0; i < num_rows; ++i) {
    if (r_var[i] == 1.0) {
      rows_to_keep[i] = true;
    }
  }

  return rows_to_keep;
}

//------------------------------------------------------------------------------
// Returns the colums to keep as a boolean vector. In order for a column to be
// kept the corresponding decision variable must equal 1.0.
//------------------------------------------------------------------------------
std::vector<bool> RowColSolver::get_cols_to_keep() const {
  std::vector<bool> cols_to_keep(num_cols, false);

  for (std::size_t j = 0; j < num_cols; ++j) {
    if (c_var[j] == 1.0) {
      cols_to_keep[j] = true;
    }
  }

  return cols_to_keep;
}

//------------------------------------------------------------------------------
// Return the number of rows keep. In order for a row to be kept the
// corresponding decision variable must equal 1.0.
//------------------------------------------------------------------------------
std::size_t RowColSolver::get_num_rows_to_keep() const {
  std::size_t count = 0;
  
  for (std::size_t i = 0; i < num_rows; ++i) {
    if (r_var[i] == 1.0) {
      ++count;
    }
  }

  return count;
}

//------------------------------------------------------------------------------
// Return the number of columns keep. In order for a column to be kept the
// corresponding decision variable must equal 1.0.
//------------------------------------------------------------------------------
std::size_t RowColSolver::get_num_cols_to_keep() const {
  std::size_t count = 0;
  
  for (std::size_t j = 0; j < num_cols; ++j) {
    // if (variables[num_rows + j] == 1.0) {
    if (c_var[j] == 1.0) {
      ++count;
    }
  }
  
  return count;
}

//------------------------------------------------------------------------------
// Takes a double representing a lower limit in the objective value.
//------------------------------------------------------------------------------
void RowColSolver::set_incumbent_obj(const double _obj) {
  use_incumbent_obj = true;
  limit = _obj;
}