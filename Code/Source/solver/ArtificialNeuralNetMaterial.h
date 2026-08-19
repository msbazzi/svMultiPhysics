// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/* This material model implementation is based on the following paper: 
Peirlinck, M., Hurtado, J.A., Rausch, M.K. et al. A universal material model subroutine 
for soft matter systems. Engineering with Computers 41, 905–927 (2025). 
https://doi.org/10.1007/s00366-024-02031-w */

#ifndef ArtificialNeuralNet_model_H
#define ArtificialNeuralNet_model_H

#include "mat_fun.h"
#include "utils.h"
#include "Parameters.h"
#include <vector>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/unsupported/Eigen/CXX11/Tensor"

using namespace mat_fun;

// Class for parameter table for material models discovered by constitutive artificial neural network (CANN)

/* This material model implementation is based on the following paper: 
Peirlinck, M., Hurtado, J.A., Rausch, M.K. et al. A universal material model subroutine 
for soft matter systems. Engineering with Computers 41, 905–927 (2025). 
https://doi.org/10.1007/s00366-024-02031-w */

class ArtificialNeuralNetMaterial
{
  public:
    static constexpr int NUM_INVARIANTS = 18;

    // Invariant indices
    Vector<int> invariant_indices;

    // Activation functions
    Array<int> activation_functions;

    // Weights
    Array<double> weights;

    // Optional GOH dispersion kappa per legacy row. Zero preserves aligned rows.
    Vector<double> row_dispersion;

    // Optional beta-distributed fiber recruitment parameters per legacy row.
    Vector<int> row_recruitment_enabled;
    Vector<double> row_recruitment_lower_stretch;
    Vector<double> row_recruitment_upper_stretch;
    Vector<double> row_recruitment_tau;
    Vector<double> row_recruitment_alpha;
    Vector<double> row_recruitment_beta;
    Vector<int> row_recruitment_quadrature_points;

    // Number of rows in parameter table
    int num_rows;

    // Outputs from each layer
    void uCANN_h0(const double x, const int kf, double &f, double &df, double &ddf) const;
    void uCANN_h1(const double x, const int kf, const double W, double &f, double &df, double &ddf) const;
    void uCANN_h2(const double x, const int kf, const double W, double &f, double &df, double &ddf) const;
    void uCANN_scalar(const double xInv, const int kf0, const int kf1, const int kf2,
           const double W0, const double W1, const double W2,
           double &psi, double &dpsi, double &ddpsi) const;
    bool rowHasDispersion(const int row) const;
    bool rowHasRecruitment(const int row) const;

    // Strain energy and derivatives
    void uCANN(const double xInv, const int kInv,
           const int kf0, const int kf1, const int kf2,
           const double W0, const double W1, const double W2,
           double &psi, double (&dpsi)[NUM_INVARIANTS], double (&ddpsi)[NUM_INVARIANTS]) const;


    void evaluate(const double aInv[NUM_INVARIANTS], double &psi,
           double (&dpsi)[NUM_INVARIANTS], double (&ddpsi)[NUM_INVARIANTS]) const;

    // Helper for compute_pk2cc
    template<size_t nsd>
    void computeInvariantsAndDerivatives(
    const Matrix<nsd>& C, const Eigen::Matrix<double, nsd, Eigen::Dynamic>& fl, int nfd, double J2d, double J4d, const Matrix<nsd>& Ci,
    const Matrix<nsd>& Idm, const double Tfa, const double kap, Matrix<nsd>& N1, double& psi, double (&Inv)[NUM_INVARIANTS],
    std::array<Matrix<nsd>,NUM_INVARIANTS>& dInv, std::array<Tensor<nsd>,NUM_INVARIANTS>& ddInv) const;

    template<size_t nsd>
    void recruitedFiberInput(const int row, const double i4_eff_m1, const Matrix<nsd>& di4_eff_m1,
        const Tensor<nsd>& ddi4_eff_m1, double& x_rec, Matrix<nsd>& dx_rec, Tensor<nsd>& ddx_rec) const;
    
};

#endif // ArtificialNeuralNet_model_H
