// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FORMULA_STRAIN_ENERGY_MATERIAL_H
#define FORMULA_STRAIN_ENERGY_MATERIAL_H

#include <string>

/// Provides a thread-local compiled evaluator for isochoric strain energy
/// density \f$\Psi(I_1,\ldots,I_9)\f$ in terms of the same invariants used
/// by the CANN model (see ArtificialNeuralNetMaterial::computeInvariantsAndDerivatives).
///
/// The expression is parsed with ExprTk. In XML use variables `I1` .. `I9`.
/// First and second derivatives of \f$\Psi\f$ with respect to invariants are
/// obtained by central finite differences (step \p h).
namespace formula_strain_energy {

void ensure_runtime(const std::string& formula);

double eval_psi(const double inv[9]);

/// @param[out] grad  dPsi/dI_k, k = 0..8
/// @param[out] hess  d2Psi/(dI_a dI_b), symmetric
void gradient_hessian_fd(const double inv[9], double h, double grad[9], double hess[9][9]);

}  // namespace formula_strain_energy

#endif
