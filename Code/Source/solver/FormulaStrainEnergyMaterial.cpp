// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "FormulaStrainEnergyMaterial.h"

#include <exprtk.hpp>

#include <array>
#include <cstring>
#include <stdexcept>

namespace formula_strain_energy {

namespace {

struct ThreadState {
  std::string formula;
  std::array<double, 9> I{};
  exprtk::symbol_table<double> sym_table;
  exprtk::expression<double> expression;
  exprtk::parser<double> parser;

  void compile(const std::string& expr_string)
  {
    formula = expr_string;
    sym_table.clear();
    sym_table.add_constants();
    for (int k = 0; k < 9; ++k) {
      I[k] = 0.0;
      const std::string name = "I" + std::to_string(k + 1);
      sym_table.add_variable(name, I[k]);
    }

    expression.register_symbol_table(sym_table);

    if (!parser.compile(expr_string, expression)) {
      throw std::runtime_error(
          "Formula strain energy: ExprTk parse error: " + parser.error() +
          "\nExpression: " + expr_string);
    }
  }
};

ThreadState& tls_state()
{
  thread_local ThreadState state;
  return state;
}

}  // namespace

void ensure_runtime(const std::string& formula)
{
  if (formula.empty()) {
    throw std::runtime_error("Formula strain energy: empty expression string.");
  }
  ThreadState& st = tls_state();
  if (st.formula != formula) {
    st.compile(formula);
  }
}

double eval_psi(const double inv[9])
{
  ThreadState& st = tls_state();
  for (int k = 0; k < 9; ++k) {
    st.I[k] = inv[k];
  }
  return st.expression.value();
}

void gradient_hessian_fd(const double inv[9], double h, double grad[9], double hess[9][9])
{
  if (h <= 0.0) {
    throw std::runtime_error("Formula strain energy: finite_difference_step must be positive.");
  }

  ThreadState& st = tls_state();
  std::array<double, 9> work{};

  for (int a = 0; a < 9; ++a) {
    for (int b = 0; b < 9; ++b) {
      hess[a][b] = 0.0;
    }
  }

  const double h2 = h * h;
  std::memcpy(work.data(), inv, 9 * sizeof(double));

  const double psi0 = eval_psi(inv);

  for (int k = 0; k < 9; ++k) {
    std::memcpy(work.data(), inv, 9 * sizeof(double));
    work[k] = inv[k] + h;
    const double pp = eval_psi(work.data());

    work[k] = inv[k] - h;
    const double pm = eval_psi(work.data());

    grad[k] = (pp - pm) / (2.0 * h);
    hess[k][k] = (pp - 2.0 * psi0 + pm) / h2;
  }

  for (int a = 0; a < 9; ++a) {
    for (int b = a + 1; b < 9; ++b) {
      std::memcpy(work.data(), inv, 9 * sizeof(double));
      work[a] = inv[a] + h;
      work[b] = inv[b] + h;
      const double f_pp = eval_psi(work.data());

      work[b] = inv[b] - h;
      const double f_pm = eval_psi(work.data());

      work[a] = inv[a] - h;
      const double f_mm = eval_psi(work.data());

      work[b] = inv[b] + h;
      const double f_mp = eval_psi(work.data());

      const double v = (f_pp - f_pm - f_mp + f_mm) / (4.0 * h * h);
      hess[a][b] = v;
      hess[b][a] = v;
    }
  }
}

} 
