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
                                          use_incumbent(false),
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
  if (use_incumbent) {
    cplex.addMIPStart(r, r_copy);
    cplex.addMIPStart(c, c_copy);
  }
  
  // Set CPLEX parameters
  // cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, .05);
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
// Takes two boolean vectors indicating which rows and columns to keep, and
// uses them to provide CPLEX an initial solution.
//------------------------------------------------------------------------------
void ElementSolver::set_incumbent(const std::vector<bool> &keep_rows, const std::vector<bool> &keep_cols) {
  assert(keep_rows.size() == num_rows);
  assert(keep_cols.size() == num_cols);
  assert(keep_rows.size() * keep_cols.size() == num_data);

  use_incumbent = true;

  for (std::size_t i = 0; i < num_rows; ++i) {
    if (keep_rows[i]) {
      r_copy[i] = 1.0;
    } else {
      r_copy[i] = 0.0;
    }
  }

  for (std::size_t j = 0; j < num_cols; ++j) {
    if (keep_cols[j]) {
      c_copy[j] = 1.0;
    } else {
      c_copy[j] = 0.0;
    }
  }

  // for (std::size_t i = 0; i < num_rows; ++i) {
  //   if (keep_rows[i]) {
  //     x_copy[num_data + i] = 1.0;
  //   } else {
  //     x_copy[num_data + i] = 0.0;
  //   }
  // }
  // for (std::size_t j = 0; j < num_cols; ++j) {
  //   if (keep_cols[j]) {
  //     x_copy[num_data + num_rows + j] = 1.0;
  //   } else {
  //     xCopy[num_data + num_rows + j] = 0.0;
  //   }
  // }
  // std::size_t idx = 0;
  // for (std::size_t i = 0; i < num_rows; ++i) {
  //   for (std::size_t j = 0; j < num_cols; ++j) {
  //     if (keep_rows[i] && keep_cols[j]) {
  //       xCopy[idx] = 1.0;
  //     } else {
  //       xCopy[idx] = 0.0;
  //     }
  //     ++idx;
  //   }
  // }
}

//------------------------------------------------------------------------------
// Takes two int vectors indicating which rows and columns to keep, and
// uses them to provide CPLEX an initial solution.
//------------------------------------------------------------------------------
void ElementSolver::set_incumbent(const std::vector<int> &keep_rows, const std::vector<int> &keep_cols) {
  assert(keep_rows.size() == num_rows);
  assert(keep_cols.size() == num_cols);
  assert(keep_rows.size() * keep_cols.size() == num_data);

  use_incumbent = true;

  for (std::size_t i = 0; i < num_rows; ++i) {
    if (keep_rows[i] == 1) {
      r_copy[i] = 1.0;
    } else if (keep_rows[i] == 0) {
      r_copy[i] = 0.0;
    }  else {
      fprintf(stderr, "ERROR - RowColSolver::set_incumbent - Unknown solution value in row %lu.\n", i);
      exit(EXIT_FAILURE);
    }
  }

  for (std::size_t j = 0; j < num_cols; ++j) {
    if (keep_cols[j] == 1) {
      c_copy[j] = 1.0;
    } else if (keep_cols[j] == 0) {
      c_copy[j] = 0.0;
    } else {
      fprintf(stderr, "ERROR - RowColSolver::set_incumbent - Unknown solution value in col %lu.\n", j);
      exit(EXIT_FAILURE);
    }
  }

  // for (std::size_t i = 0; i < num_rows; ++i) {
  //   if (keep_rows[i]) {
  //     x_copy[num_data + i] = 1.0;
  //   } else {
  //     x_copy[num_data + i] = 0.0;
  //   }
  // }
  // for (std::size_t j = 0; j < num_cols; ++j) {
  //   if (keep_cols[j]) {
  //     x_copy[num_data + num_rows + j] = 1.0;
  //   } else {
  //     xCopy[num_data + num_rows + j] = 0.0;
  //   }
  // }
  // std::size_t idx = 0;
  // for (std::size_t i = 0; i < num_rows; ++i) {
  //   for (std::size_t j = 0; j < num_cols; ++j) {
  //     if (keep_rows[i] && keep_cols[j]) {
  //       xCopy[idx] = 1.0;
  //     } else {
  //       xCopy[idx] = 0.0;
  //     }
  //     ++idx;
  //   }
  // }
}

// //------------------------------------------------------------------------------
// // Constructor.
// //------------------------------------------------------------------------------
// ElementSolver::ElementSolver(const BinContainer &_data,
//                      const std::string &_na_symbol,                     
//                      const double _min_perc_missing,
//                      const double _TOL) : data(&_data),
//                                           na_symbol(_na_symbol),
//                                           TOL(_TOL),
//                                           max_perc_missing(_min_perc_missing),
//                                           num_rows(data->get_num_data_rows()),
//                                           num_cols(data->get_num_data_cols()),
//                                           num_data(num_rows * num_cols),
//                                           variables(num_data + num_rows + num_cols),
//                                           // dVar(num_data),
//                                           // rVar(num_rows),
//                                           // cVar(num_cols),
//                                           obj_value(0),
//                                           use_incumbent(false),
//                                           env(IloEnv()),
//                                           cplex(IloCplex(env)),
//                                           model(IloModel(env)),
//                                           x(IloNumVarArray(env, variables.size(), 0, 1, ILOINT)),
//                                           xCopy(IloNumArray(env, variables.size())),
//                                           // d(IloNumVarArray(env, dVar.size(), 0, 1, ILOINT)),
//                                           // r(IloNumVarArray(env, rVar.size(), 0, 1, ILOINT)),
//                                           // c(IloNumVarArray(env, cVar.size(), 0, 1, ILOINT)),
//                                           // dCopy(IloNumArray(env, dVar.size())),
//                                           // rCopy(IloNumArray(env, rVar.size())),                                          
//                                           // cCopy(IloNumArray(env, cVar.size())),
//                                           obj(IloExpr(env))
// {
//   build_model();
// }

// //------------------------------------------------------------------------------
// // Destructor.
// //------------------------------------------------------------------------------
// ElementSolver::~ElementSolver() {}

// std::size_t ElementSolver::get_data_idx(const std::size_t i, const std::size_t j) const {
//   assert(i < num_rows);
//   assert(j < num_cols);
  
//   return j + (i * num_cols);
// }

// //------------------------------------------------------------------------------
// // Builds the CPLEX model based on the BinContainer object.
// //------------------------------------------------------------------------------
// void ElementSolver::build_model() {
//   // Add objective function
//   std::size_t idx = 0;
//   for (std::size_t i = 0; i < num_rows; ++i) {
//     for (std::size_t j = 0; j < num_cols; ++j) {
//       if (!data->is_data_na(i,j)) {
//         // obj += d[idx];
//         obj += x[idx];
//       }
//       ++idx;
//     }
//   }  
//   model.add(IloMaximize(env, obj, "Objective"));

//   idx = 0;
//   for (std::size_t i = 0; i < num_rows; ++i) {
//     for (std::size_t j = 0; j < num_cols; ++j) {
//       // model.add(d[idx] <= r[i]);
//       // model.add(d[idx] <= c[j]);
//       model.add(x[idx] <= x[num_data + i]);
//       model.add(x[idx] <= x[num_data + num_rows + j]);
//       // model.add(x[idx] <= static_cast<IloNum>(0.5) * x[num_data + i] + static_cast<IloNum>(0.5) * x[num_data + num_rows + j]);
//       ++idx;
//     }
//   }

//   // Add contraints on rows
//   for (std::size_t i = 0; i < num_rows; ++i) {
//     IloExpr col_sum(env);
//     for (std::size_t j = 0; j < num_cols; ++j) {
//       double b = data->is_data_na(i,j) ? 0.0 : 1.0;
//       idx = get_data_idx(i, j);
//       col_sum -= x[num_data + num_rows + j] * static_cast<IloNum>((max_perc_missing + b - 1.0) / num_cols);
//       // col_sum -= c[j] * static_cast<IloNum>((max_perc_missing + b - 1.0) / num_cols);
//     }
//     model.add(x[num_data + i] + col_sum <= 1);
//     // model.add(r[i] + col_sum <= 1);
//     col_sum.end();
//   }
  
//   // Add constraints on cols
//   for (std::size_t j = 0; j < num_cols; ++j) {
//     IloExpr row_sum(env);
//     for (std::size_t i = 0; i < num_rows; ++i) {
//       double b = data->is_data_na(i,j) ? 0.0 : 1.0;
//       idx = get_data_idx(i, j);
//       row_sum -= x[num_data + i] * static_cast<IloNum>((max_perc_missing + b - 1.0) / num_rows);
//       // row_sum -= r[i] * static_cast<IloNum>((max_perc_missing + b - 1.0) / num_rows);
//     }
//     model.add(x[num_data + num_rows + j] + row_sum <= 1);
//     // model.add(c[j] + row_sum <= 1);
//     row_sum.end();
//   }

//   // for (std::size_t i = 0; i < num_rows; ++i) {
//   //   for (std::size_t j = 0; j < num_cols; ++j) {
//   //     if (data->is_data_na(i,j)) {
//   //       model.add(x[num_data + i] + x[num_data + num_rows + j] <= static_cast<IloNum>(1));
//   //     }
//   //   }
//   // }
// }

// //------------------------------------------------------------------------------
// // Rounds values close to 1 or 0 to 1 or 0 respectively. Values must be within
// // 'TOL' of 1 or 0 to be rounded. Row decision variables, column decision 
// // variables, and element decision are all rounded.
// //------------------------------------------------------------------------------
// void ElementSolver::round_extreme_values() {
//   for (std::size_t i = 0; i < variables.size(); ++i) {
//     if (variables[i] > 1.0 - TOL)
//       variables[i] = 1.0;
//     else if (variables[i] < TOL)
//       variables[i] = 0.0;
//   }
//   // for (std::size_t i = 0; i < dVar.size(); ++i) {
//   //   if (dVar[i] > 1.0 - TOL)
//   //     dVar[i] = 1.0;
//   //   else if (dVar[i] < TOL)
//   //     dVar[i] = 0.0;
//   // }
//   // for (std::size_t i = 0; i < rVar.size(); ++i) {
//   //   if (rVar[i] > 1.0 - TOL)
//   //     rVar[i] = 1.0;
//   //   else if (rVar[i] < TOL)
//   //     rVar[i] = 0.0;
//   // }
//   // for (std::size_t i = 0; i < cVar.size(); ++i) {
//   //   if (cVar[i] > 1.0 - TOL)
//   //     cVar[i] = 1.0;
//   //   else if (cVar[i] < TOL)
//   //     cVar[i] = 0.0;
//   // }
// }

// //------------------------------------------------------------------------------
// // Calls CPLEX solver and rounds decision variables if a valid solution is
// // found.
// //------------------------------------------------------------------------------
// void ElementSolver::solve() {
//   cplex.extract(model);
//   if (use_incumbent) {
//     cplex.addMIPStart(x, xCopy);
//   }
  
//   // cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, .05);
//   cplex.setParam(IloCplex::Param::RandomSeed, 0);
//   cplex.setParam(IloCplex::Param::Threads, 1);
//   cplex.setOut(env.getNullStream());

//   obj_value = 0.0;
//   try {
//     cplex.solve();

//     if (cplex.getStatus() == IloAlgorithm::Infeasible) {
//       printf("Infeasible Solution\n");
//     } else if (cplex.getStatus() == IloAlgorithm::Optimal) {
//       obj_value = cplex.getObjValue();
//       cplex.getValues(x, xCopy);
//       // cplex.getValues(d, dCopy);
//       // cplex.getValues(r, rCopy);
//       // cplex.getValues(c, cCopy);

//       for (std::size_t i = 0; i < variables.size(); ++i) {
//         variables[i] = xCopy[i];
//       }
//       // for (std::size_t i = 0; i < dVar.size(); ++i) {
//       //   dVar[i] = dCopy[i];
//       // }
//       // for (std::size_t i = 0; i < rVar.size(); ++i) {
//       //   rVar[i] = rCopy[i];
//       // }
//       // for (std::size_t i = 0; i < cVar.size(); ++i) {
//       //   cVar[i] = cCopy[i];
//       // }
//       round_extreme_values();
//     }
//   } catch (IloException& e) {
//     std::cerr << "Concert exception caught: " << e << std::endl;
//   } catch (...) {
//     std::cerr << "Unknown exception caught" << std::endl;
//   }
// }

// //------------------------------------------------------------------------------
// // Returns the objective value from CPLEX. If CPLEX has not been called, or did
// // not find a solution, then 0 is returned.
// //------------------------------------------------------------------------------
// double ElementSolver::get_obj_value() const {
//   return obj_value;
// }

// //------------------------------------------------------------------------------
// // Returns the rows to keep as a boolean vector. In order for a row to be kept
// // the corresponding decision variable must equal 1.0.
// //------------------------------------------------------------------------------
// std::vector<bool> ElementSolver::get_rows_to_keep() const {
//   std::vector<bool> rows_to_keep(num_rows, false);
//   for (std::size_t i = 0; i < num_rows; ++i) {
//     if (variables[num_data + i] == 1.0) {
//     // if (rVar[i] == 1.0) {
//       rows_to_keep[i] = true;
//     }
//   }
//   return rows_to_keep;
// }

// //------------------------------------------------------------------------------
// // Returns the colums to keep as a boolean vector. In order for a column to be
// // kept the corresponding decision variable must equal 1.0.
// //------------------------------------------------------------------------------
// std::vector<bool> ElementSolver::get_cols_to_keep() const {
//   std::vector<bool> cols_to_keep(num_cols, false);
//   for (std::size_t j = 0; j < num_cols; ++j) {
//     if (variables[num_data + num_rows + j] == 1.0) {
//     // if (cVar[j] == 1.0) {
//       cols_to_keep[j] = true;
//     }
//   }
//   return cols_to_keep;
// }

// //------------------------------------------------------------------------------
// // Return the number of rows keep. In order for a row to be kept the
// // corresponding decision variable must equal 1.0.
// //------------------------------------------------------------------------------
// std::size_t ElementSolver::get_num_rows_to_keep() const {
//   std::size_t count = 0;
//   for (std::size_t i = 0; i < num_rows; ++i) {
//     if (variables[num_data + i] == 1.0) {
//     // if (rVar[i] == 1.0) {
//       ++count;
//     }
//   }
//   return count;
// }

// //------------------------------------------------------------------------------
// // Return the number of columns keep. In order for a column to be kept the
// // corresponding decision variable must equal 1.0.
// //------------------------------------------------------------------------------
// std::size_t ElementSolver::get_num_cols_to_keep() const {
//   std::size_t count = 0;  
//   for (std::size_t j = 0; j < num_cols; ++j) {
//     if (variables[num_data + num_rows + j] == 1.0) {
//     // if (cVar[j] == 1.0) {
//       ++count;
//     }
//   }
//   return count;
// }

// //------------------------------------------------------------------------------
// // Takes two boolean vectors indicating which rows and columns to keep, and
// // uses them to provide CPLEX an initial solution.
// //------------------------------------------------------------------------------
// void ElementSolver::set_incumbent(const std::vector<bool> &keep_rows, const std::vector<bool> &keep_cols) {
//   assert(keep_rows.size() == num_rows);
//   assert(keep_cols.size() == num_cols);
//   assert(keep_rows.size() * keep_cols.size() == num_data);

//   use_incumbent = true;

//   for (std::size_t i = 0; i < num_rows; ++i) {
//     if (keep_rows[i]) {
//       xCopy[num_data + i] = 1.0;
//     } else {
//       xCopy[num_data + i] = 0.0;
//     }
//   }

//   for (std::size_t j = 0; j < num_cols; ++j) {
//     if (keep_cols[j]) {
//       xCopy[num_data + num_rows + j] = 1.0;
//     } else {
//       xCopy[num_data + num_rows + j] = 0.0;
//     }
//   }

//   std::size_t idx = 0;
//   for (std::size_t i = 0; i < num_rows; ++i) {
//     for (std::size_t j = 0; j < num_cols; ++j) {
//       if (keep_rows[i] && keep_cols[j]) {
//         xCopy[idx] = 1.0;
//       } else {
//         xCopy[idx] = 0.0;
//       }
//       ++idx;
//     }
//   }
// }