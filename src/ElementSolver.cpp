#include <assert.h>
#include "ElementSolver.h"

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
ElementSolver::ElementSolver(const BinContainer &_data,      
                     const double _min_perc_missing,
                     const double _TOL) : data(&_data),
                                          TOL(_TOL),
                                          max_perc_missing(_min_perc_missing),
                                          num_rows(data->get_num_data_rows()),
                                          num_cols(data->get_num_data_cols()),
                                          num_data(num_rows * num_cols),
                                          x_var(num_data),
                                          r_var(num_rows),
                                          c_var(num_cols),
                                          obj_value(0),
                                          use_incumbent_obj(false),
                                          limit(0.0),
                                          env(IloEnv()),
                                          cplex(IloCplex(env)),
                                          model(IloModel(env)),
                                          x(IloNumVarArray(env, num_data, 0, 1, ILOINT)),
                                          r(IloNumVarArray(env, num_rows, 0, 1, ILOINT)),
                                          c(IloNumVarArray(env, num_cols, 0, 1, ILOINT)),
                                          x_copy(IloNumArray(env, num_data)),                                          
                                          r_copy(IloNumArray(env, num_rows)),
                                          c_copy(IloNumArray(env, num_cols)),
                                          obj(IloExpr(env))
{
  build_model();
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
ElementSolver::~ElementSolver() {}

std::size_t ElementSolver::get_data_idx(const std::size_t i, const std::size_t j) const {
  assert(i < num_rows);
  assert(j < num_cols);
  
  return j + (i * num_cols);
}

//------------------------------------------------------------------------------
// Builds the CPLEX model based on the BinContainer object.
//------------------------------------------------------------------------------
void ElementSolver::build_model() {
  // Add objective function
  std::size_t idx = 0;
  for (std::size_t i = 0; i < num_rows; ++i) {
    for (std::size_t j = 0; j < num_cols; ++j) {
      if (!data->is_data_na(i,j)) {
        obj += x[idx];
      }
      ++idx;
    }
  }  
  model.add(IloMaximize(env, obj, "Objective"));

  // Add constraints for each element based on the corresponding row and column
  idx = 0;
  for (std::size_t i = 0; i < num_rows; ++i) {
    for (std::size_t j = 0; j < num_cols; ++j) {
      model.add(x[idx] <= static_cast<IloNum>(0.5) * r[i] + static_cast<IloNum>(0.5) * c[j]);
      ++idx;
    }
  }

  // Add contraints on rows
  for (std::size_t i = 0; i < num_rows; ++i) {
    IloExpr col_sum(env);
    for (std::size_t j = 0; j < num_cols; ++j) {
      double b = data->is_data_na(i,j) ? 0.0 : 1.0;
      col_sum -= c[j] * static_cast<IloNum>((max_perc_missing + b - 1.0) / num_cols);
    }
    model.add(r[i] + col_sum <= 1);
    col_sum.end();
  }
  
  // Add constraints on cols
  for (std::size_t j = 0; j < num_cols; ++j) {
    IloExpr row_sum(env);
    for (std::size_t i = 0; i < num_rows; ++i) {
      double b = data->is_data_na(i,j) ? 0.0 : 1.0;
      row_sum -= r[i] * static_cast<IloNum>((max_perc_missing + b - 1.0) / num_rows);
    }
    model.add(c[j] + row_sum <= 1);
    row_sum.end();
  }
}

//------------------------------------------------------------------------------
// Rounds values close to 1 or 0 to 1 or 0 respectively. Values must be within
// 'TOL' of 1 or 0 to be rounded. Row decision variables, column decision 
// variables, and element decision are all rounded.
//------------------------------------------------------------------------------
void ElementSolver::round_extreme_values() {
  for (std::size_t i = 0; i < x_var.size(); ++i) {
    if (x_var[i] > 1.0 - TOL) {
      x_var[i] = 1.0;
    } else if (x_var[i] < TOL) {
      x_var[i] = 0.0;
    }
  }
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
void ElementSolver::solve() {
  // Extract CPLEX model
  cplex.extract(model);

  // Set incumbent solution if provided
  if (use_incumbent_obj) {
    cplex.setParam(IloCplex::Param::MIP::Tolerances::LowerCutoff, limit);
  }
  
  // Set CPLEX parameters
  // cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, .05);
  cplex.setParam(IloCplex::Param::TimeLimit, 18000);
  cplex.setParam(IloCplex::Param::RandomSeed, 0);
  cplex.setParam(IloCplex::Param::Threads, 1);
  // cplex.setOut(env.getNullStream());

  obj_value = 0.0;
  try {
    cplex.solve();

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
      printf("Infeasible Solution\n");
    } else if (cplex.getStatus() == IloAlgorithm::Optimal || cplex.getStatus() == IloAlgorithm::Feasible) {
      obj_value = cplex.getObjValue(); // Get objective value

      // Copy all decision variables
      cplex.getValues(x, x_copy);
      cplex.getValues(r, r_copy);
      cplex.getValues(c, c_copy);
      for (std::size_t i = 0; i < x_var.size(); ++i) {
        x_var[i] = x_copy[i];
      }
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
double ElementSolver::get_obj_value() const {
  return obj_value;
}

//------------------------------------------------------------------------------
// Returns the rows to keep as a boolean vector. In order for a row to be kept
// the corresponding decision variable must equal 1.0.
//------------------------------------------------------------------------------
std::vector<bool> ElementSolver::get_rows_to_keep() const {
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
std::vector<bool> ElementSolver::get_cols_to_keep() const {
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
std::size_t ElementSolver::get_num_rows_to_keep() const {
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
std::size_t ElementSolver::get_num_cols_to_keep() const {
  std::size_t count = 0;  
  for (std::size_t j = 0; j < num_cols; ++j) {
    if (c_var[j] == 1.0) {
      ++count;
    }
  }
  return count;
}

//------------------------------------------------------------------------------
// Takes a double representing a lower limit in the objective value.
//------------------------------------------------------------------------------
void ElementSolver::set_incumbent_obj(const double _obj) {
  use_incumbent_obj = true;
  limit = _obj;
}