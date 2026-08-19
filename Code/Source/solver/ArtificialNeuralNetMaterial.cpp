// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/* This material model implementation is based on the following paper: 
Peirlinck, M., Hurtado, J.A., Rausch, M.K. et al. A universal material model subroutine 
for soft matter systems. Engineering with Computers 41, 905–927 (2025). 
https://doi.org/10.1007/s00366-024-02031-w */

// The functions required for CANN material model implementations are defined here

#include "ArtificialNeuralNetMaterial.h"
#include "ComMod.h"
#include "mat_fun.h"
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>
using namespace mat_fun;

namespace {

double stable_softplus(const double x)
{
    if (x > 40.0) {
        return x;
    }
    if (x < -40.0) {
        return std::exp(x);
    }
    return std::log1p(std::exp(x));
}

double stable_sigmoid(const double x)
{
    if (x >= 0.0) {
        const double exp_neg_x = std::exp(-x);
        return 1.0 / (1.0 + exp_neg_x);
    }
    const double exp_x = std::exp(x);
    return exp_x / (1.0 + exp_x);
}

std::pair<std::vector<double>, std::vector<double>> gauss_legendre_rule(const int n)
{
    const double eps = 1.0e-14;
    const double pi = std::acos(-1.0);
    std::vector<double> x(n);
    std::vector<double> w(n);
    const int m = (n + 1) / 2;
    for (int i = 0; i < m; ++i) {
        double z = std::cos(pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n) + 0.5));
        double z_prev = 0.0;
        double p1 = 0.0;
        double p2 = 0.0;
        double pp = 0.0;
        do {
            p1 = 1.0;
            p2 = 0.0;
            for (int j = 1; j <= n; ++j) {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z_prev = z;
            z = z_prev - p1 / pp;
        } while (std::abs(z - z_prev) > eps);

        x[i] = -z;
        x[n - 1 - i] = z;
        const double weight = 2.0 / ((1.0 - z * z) * pp * pp);
        w[i] = weight;
        w[n - 1 - i] = weight;
    }
    for (int i = 0; i < n; ++i) {
        x[i] = 0.5 * (x[i] + 1.0);
        w[i] *= 0.5;
    }
    return {x, w};
}

}

/// @brief 0th layer output of CANN for activation func kf, input x
void ArtificialNeuralNetMaterial::uCANN_h0(const double x, const int kf, double &f, double &df, double &ddf) const {
    if (kf == 1) {
        f = x;
        df = 1;
        ddf = 0;
    } else if (kf == 2) {
        if (x == 0) {
            f = (std::abs(x) + x) / 2;
            df = 0;
            ddf = 0;
        }
        else {
            f = (std::abs(x) + x) / 2;
            df = 0.5 * (std::abs(x) / x + 1);
            ddf = 0;
        }
    } else if (kf == 3) {
        f = std::abs(x);
        df = std::abs(x) / x;
        ddf = 0;
    }
}

/// @brief 1st layer output of CANN for activation func kf, input x, weight W
void ArtificialNeuralNetMaterial::uCANN_h1(const double x, const int kf, const double W, double &f, double &df, double &ddf) const {
    if (kf == 1) {
        f = W * x;
        df = W;
        ddf = 0;
    } else if (kf == 2) {
        f = W * W * x * x;
        df = 2 * W * W * x;
        ddf = 2 * W * W;
    } else if (kf == 3) {
        f = W * W * W * x * x * x;
        df = 3 * W * W * W * x * x;
        ddf = 6 * W * W * W * x;
    }
}

/// @brief 2nd layer output of CANN for activation func kf, input x, weight W
void ArtificialNeuralNetMaterial::uCANN_h2(const double x, const int kf, const double W, double &f, double &df, double &ddf) const {
    if (kf == 1) {
        f = W * x;
        df = W;
        ddf = 0;
    } else if (kf == 2) {
        f = std::exp(W * x) - 1;
        df = W * std::exp(W * x);
        ddf = W * W * std::exp(W * x);
    } else if (kf == 3) {
        f = -std::log(1 - W * x);
        df = W / (1 - W * x);
        ddf = -W * W / ((1 - W * x) * (1 - W * x));
    }
}

void ArtificialNeuralNetMaterial::uCANN_scalar(
    const double xInv, const int kf0, const int kf1, const int kf2,
    const double W0, const double W1, const double W2,
    double &psi, double &dpsi, double &ddpsi) const
{
    double f0, df0, ddf0;
    uCANN_h0(xInv, kf0, f0, df0, ddf0);
    double f1, df1, ddf1;
    uCANN_h1(f0, kf1, W0, f1, df1, ddf1);
    double f2, df2, ddf2;
    uCANN_h2(f1, kf2, W1, f2, df2, ddf2);

    psi = W2 * f2;
    dpsi = W2 * df2 * df1 * df0;
    ddpsi = W2 * ((ddf2 * df1 * df1 + df2 * ddf1) * df0 * df0 + df2 * df1 * ddf0);
}

bool ArtificialNeuralNetMaterial::rowHasDispersion(const int row) const
{
    return row >= 0 && row < row_dispersion.size() && row_dispersion(row) != 0.0;
}

bool ArtificialNeuralNetMaterial::rowHasRecruitment(const int row) const
{
    return row >= 0 && row < row_recruitment_enabled.size() && row_recruitment_enabled(row) != 0;
}

template<size_t nsd>
void ArtificialNeuralNetMaterial::recruitedFiberInput(
    const int row, const double i4_eff_m1, const Matrix<nsd>& di4_eff_m1,
    const Tensor<nsd>& ddi4_eff_m1, double& x_rec, Matrix<nsd>& dx_rec, Tensor<nsd>& ddx_rec) const
{
    const double lower = row_recruitment_lower_stretch(row);
    const double upper = row_recruitment_upper_stretch(row);
    const double tau = row_recruitment_tau(row);
    const double alpha = row_recruitment_alpha(row);
    const double beta = row_recruitment_beta(row);
    const int nqp = row_recruitment_quadrature_points(row);

    const double i4_eff = i4_eff_m1 + 1.0;
    const double i4_safe = std::max(i4_eff, 1.0e-12);
    const double lambda = std::sqrt(i4_safe);
    const double dlambda_di4 = 0.5 / lambda;
    const double ddlambda_di4 = -0.25 / (lambda * lambda * lambda);

    const auto quadrature = gauss_legendre_rule(nqp);
    const auto& points = quadrature.first;
    const auto& weights = quadrature.second;
    const double log_beta_norm = std::lgamma(alpha) + std::lgamma(beta) - std::lgamma(alpha + beta);

    double integral = 0.0;
    double dintegral_dlambda = 0.0;
    double ddintegral_dlambda = 0.0;
    for (int q = 0; q < nqp; ++q) {
        const double u = points[q];
        if (u <= 0.0 || u >= 1.0) {
            continue;
        }
        const double lambda_s = lower + (upper - lower) * u;
        const double log_density = (alpha - 1.0) * std::log(u) + (beta - 1.0) * std::log1p(-u) - log_beta_norm;
        const double density = std::exp(log_density);
        const double y = lambda / lambda_s - 1.0;
        const double z = y / tau;
        const double sig = stable_sigmoid(z);
        const double r = tau * stable_softplus(z);
        const double dr_dlambda = sig / lambda_s;
        const double ddr_dlambda = sig * (1.0 - sig) / (tau * lambda_s * lambda_s);
        const double weighted_density = weights[q] * density;
        integral += weighted_density * r;
        dintegral_dlambda += weighted_density * dr_dlambda;
        ddintegral_dlambda += weighted_density * ddr_dlambda;
    }

    const double lambda_rec = 1.0 + integral;
    const double dlambda_rec_di4 = dintegral_dlambda * dlambda_di4;
    const double ddlambda_rec_di4 =
        ddintegral_dlambda * dlambda_di4 * dlambda_di4 + dintegral_dlambda * ddlambda_di4;

    x_rec = lambda_rec * lambda_rec - 1.0;
    const double dx_di4 = 2.0 * lambda_rec * dlambda_rec_di4;
    const double ddx_di4 = 2.0 * (dlambda_rec_di4 * dlambda_rec_di4 + lambda_rec * ddlambda_rec_di4);
    const Matrix<nsd> input_derivative = di4_eff_m1;
    const Tensor<nsd> input_second_derivative = ddi4_eff_m1;
    dx_rec = dx_di4 * input_derivative;
    ddx_rec = dx_di4 * input_second_derivative +
        ddx_di4 * dyadic_product<nsd>(input_derivative, input_derivative);
}

/// @brief Updates psi and its derivatives
void ArtificialNeuralNetMaterial::uCANN(
    const double xInv, const int kInv,
    const int kf0, const int kf1, const int kf2,
    const double W0, const double W1, const double W2,
    double &psi, double (&dpsi)[NUM_INVARIANTS], double (&ddpsi)[NUM_INVARIANTS]
) const {
    double row_psi = 0.0, row_dpsi = 0.0, row_ddpsi = 0.0;
    uCANN_scalar(xInv, kf0, kf1, kf2, W0, W1, W2, row_psi, row_dpsi, row_ddpsi);
    psi += row_psi;
    dpsi[kInv - 1] += row_dpsi;
    ddpsi[kInv - 1] += row_ddpsi;
}

/// @brief function to build psi and dpsidI1 to supported invariants
void ArtificialNeuralNetMaterial::evaluate(const double aInv[NUM_INVARIANTS], double &psi,
    double (&dpsi)[NUM_INVARIANTS], double (&ddpsi)[NUM_INVARIANTS]) const {
    // Initializing
    psi = 0;
    for (int i = 0; i < NUM_INVARIANTS; ++i) {
        dpsi[i] = 0;
        ddpsi[i] = 0;
    }

    double ref[NUM_INVARIANTS] = {3, 3, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (int i = 0; i < num_rows; ++i) {
        int kInv = this->invariant_indices(i);
        int kf0 = this->activation_functions(i, 0);
        int kf1 = this->activation_functions(i, 1);
        int kf2 = this->activation_functions(i, 2);
        double W0 = this->weights(i, 0);
        double W1 = this->weights(i, 1);
        double W2 = this->weights(i, 2);

        double xInv = aInv[kInv - 1] - ref[kInv - 1];
        const double kappa = rowHasDispersion(i) ? row_dispersion(i) : 0.0;
        if (kappa != 0.0) {
            xInv = kappa * (aInv[0] - ref[0]) + (1.0 - 3.0 * kappa) * xInv;
        }
        uCANN(xInv, kInv, kf0, kf1, kf2, W0, W1, W2, psi, dpsi, ddpsi);
    }
}

template<size_t nsd>
void ArtificialNeuralNetMaterial::computeInvariantsAndDerivatives(
const Matrix<nsd>& C, const Eigen::Matrix<double, nsd, Eigen::Dynamic>& fl, int nfd, double J2d, double J4d, const Matrix<nsd>& Ci,
const Matrix<nsd>& Idm, const double Tfa, const double kap, Matrix<nsd>& N1, double& psi, double (&Inv)[NUM_INVARIANTS],
std::array<Matrix<nsd>,NUM_INVARIANTS>& dInv, std::array<Tensor<nsd>,NUM_INVARIANTS>& ddInv) const {

    Matrix<nsd> C2 = C * C;
    const double fiber_scale = 1.0 - 3.0*kap;

    Inv[0] = J2d * C.trace();
    Inv[1] = 0.50 * (Inv[0]*Inv[0] - J4d * (C*C).trace());
    Inv[2] = C.determinant();
    N1 = fl.col(0)*fl.col(0).transpose();
    Matrix<nsd> N1_disp = kap*Idm + fiber_scale*N1;
    Inv[3] = kap*Inv[0] + fiber_scale * J2d * (fl.col(0).dot(C * fl.col(0)));
    Inv[4] = kap*J4d*C2.trace() + fiber_scale * J4d * (fl.col(0).dot(C2 * fl.col(0)));

    bool uses_theta_or_z_strain = false;
    for (int row = 0; row < num_rows; ++row) {
        const int invariant = invariant_indices(row);
        if (invariant == 10 || invariant == 11 || invariant == 18) {
            uses_theta_or_z_strain = true;
            break;
        }
    }
    if (uses_theta_or_z_strain && !(nfd == 2 || nfd >= 4)) {
        throw std::runtime_error(
            "[CANN] Etheta/Ez invariant rows require either two D1/D2 fiber directions "
            "or four f/s/el/smc fiber directions.");
    }

    Matrix<nsd> dInv1 = -Inv[0]/3 * Ci + J2d * Idm;
    Matrix<nsd> dInv2 = (C2.trace()/3)*Ci + Inv[0]*dInv1 + J4d*C;
    Matrix<nsd> dInv3 = Inv[2]*Ci;
    Matrix<nsd> dInv4 = -Inv[3]/3*Ci + J2d*N1_disp;
    Matrix<nsd> dInv5 = J4d*(N1_disp*C + C*N1_disp) - Inv[4]/3*Ci;

    Matrix<nsd> dInv6, dInv7, dInv8, dInv9, dInv10, dInv11;
    Tensor<nsd> ddInv6, ddInv7, ddInv8, ddInv9, ddInv10, ddInv11;
    Matrix<nsd> dE11, dE22, dE33, dE12, dE23, dE13;
    Tensor<nsd> ddE11, ddE22, ddE33, ddE12, ddE23, ddE13;
    Matrix<nsd> dEthetaEz;
    Tensor<nsd> ddEthetaEz;
    dInv6.setZero();
    dInv7.setZero();
    dInv8.setZero();
    dInv9.setZero();
    dInv10.setZero();
    dInv11.setZero();
    ddInv6.setZero();
    ddInv7.setZero();
    ddInv8.setZero();
    ddInv9.setZero();
    ddInv10.setZero();
    ddInv11.setZero();
    dE11.setZero();
    dE22.setZero();
    dE33.setZero();
    dE12.setZero();
    dE23.setZero();
    dE13.setZero();
    ddE11.setZero();
    ddE22.setZero();
    ddE33.setZero();
    ddE12.setZero();
    ddE23.setZero();
    ddE13.setZero();
    dEthetaEz.setZero();
    ddEthetaEz.setZero();

    Inv[11] = 0.5 * (C(0,0) - 1.0);
    Inv[12] = 0.5 * (C(1,1) - 1.0);
    dE11(0,0) = 0.5;
    dE22(1,1) = 0.5;
    if constexpr (nsd == 3) {
        Inv[13] = 0.5 * (C(2,2) - 1.0);
        Inv[14] = 0.5 * C(0,1);
        Inv[15] = 0.5 * C(1,2);
        Inv[16] = 0.5 * C(0,2);
        dE33(2,2) = 0.5;
        dE12(0,1) = 0.25;
        dE12(1,0) = 0.25;
        dE23(1,2) = 0.25;
        dE23(2,1) = 0.25;
        dE13(0,2) = 0.25;
        dE13(2,0) = 0.25;
    } else {
        Inv[13] = 0.0;
        Inv[14] = 0.5 * C(0,1);
        Inv[15] = 0.0;
        Inv[16] = 0.0;
        dE12(0,1) = 0.25;
        dE12(1,0) = 0.25;
    }

    Tensor<nsd> dCidC = -symmetric_dyadic_product<nsd>(Ci, Ci);
    Matrix<nsd> dJ4ddC = -2.0/3.0 * J4d * Ci;

    Tensor<nsd> ddInv1 = (-1.0/3.0)*(dyadic_product<nsd>(dInv1,Ci) + Inv[0]*dCidC + J2d*dyadic_product(Ci,Idm));
    Tensor<nsd> ddInv2 = dyadic_product<nsd>(dInv1,dInv1) + Inv[0]*ddInv1 + (1.0/3.0)*C2.trace()*dCidC 
                        + (1.0/3.0)*dyadic_product<nsd>((C2.trace()*dJ4ddC + 2*J4d*C),Ci) 
                        + dyadic_product<nsd>(dJ4ddC,C) - J4d*fourth_order_identity<nsd>();
    Tensor<nsd> ddInv3 = dyadic_product<nsd>(dInv3,Ci) + Inv[2]*dCidC;
    Tensor<nsd> ddInv4 = (-1.0/3.0)*(dyadic_product<nsd>(dInv4,Ci) + J2d*dyadic_product<nsd>(Ci,N1_disp) + Inv[3]*dCidC);
    Matrix<nsd> sum1 = (N1_disp*C + C*N1_disp);
    Tensor<nsd> ddInv5 = (-1.0/3.0)*(dyadic_product<nsd>(dInv5,Ci) + Inv[4]*dCidC + 2*J4d*dyadic_product(Ci,sum1))
                        + J4d*(2*symmetric_dyadic_product(N1_disp,Idm) - dyadic_product(N1_disp,Idm)
                        + 2*symmetric_dyadic_product(Idm,N1_disp) - dyadic_product(Idm,N1_disp));

    if (nfd >= 2) {
        Matrix<nsd> N2 = fl.col(1)*fl.col(1).transpose();
        Matrix<nsd> N12 = 0.5*(fl.col(0)*fl.col(1).transpose() + fl.col(1)*fl.col(0).transpose());
        Matrix<nsd> N2_disp = kap*Idm + fiber_scale*N2;
        Matrix<nsd> N12_disp = fiber_scale*N12;

        Inv[5] = fiber_scale * J2d * (fl.col(0).dot(C * fl.col(1)));
        Inv[6] = fiber_scale * J4d * (fl.col(0).dot(C2 * fl.col(1)));
        Inv[7] = kap*Inv[0] + fiber_scale * J2d * (fl.col(1).dot(C * fl.col(1)));
        Inv[8] = kap*J4d*C2.trace() + fiber_scale * J4d * (fl.col(1).dot(C2 * fl.col(1)));

        dInv6 = -Inv[5]/3*Ci + J2d*N12_disp;
        dInv7 = J4d*(N12_disp*C + C*N12_disp) - Inv[6]/3*Ci;
        dInv8 = -Inv[7]/3*Ci + J2d*N2_disp;
        dInv9 = J4d*(N2_disp*C + C*N2_disp) - Inv[8]/3*Ci;

        ddInv6 = -1.0/3.0*(dyadic_product(dInv6,Ci) + J2d*dyadic_product(Ci,N12_disp) + Inv[5]*dCidC);
        Matrix<nsd> sum12 = (N12_disp*C + C*N12_disp);
        ddInv7 = -1.0/3.0*(dyadic_product(dInv7,Ci) + Inv[6]*dCidC + 2*J4d*dyadic_product(Ci,sum12))
                + J4d*(2*symmetric_dyadic_product(N12_disp,Idm) - dyadic_product(N12_disp,Idm)
                + 2*symmetric_dyadic_product(Idm,N12_disp) - dyadic_product(Idm,N12_disp));
        ddInv8 = -1.0/3.0*(dyadic_product(dInv8,Ci) + J2d*dyadic_product(Ci,N2_disp) + Inv[7]*dCidC);
        Matrix<nsd> sum2 = (N2_disp*C + C*N2_disp);
        ddInv9 = -1.0/3.0*(dyadic_product(dInv9,Ci) + Inv[8]*dCidC + 2*J4d*dyadic_product(Ci,sum2))
                + J4d*(2*symmetric_dyadic_product(N2_disp,Idm) - dyadic_product(N2_disp,Idm)
                + 2*symmetric_dyadic_product(Idm,N2_disp) - dyadic_product(Idm,N2_disp));

        if (uses_theta_or_z_strain) {
            Eigen::Matrix<double, nsd, 1> theta_vec;
            Eigen::Matrix<double, nsd, 1> z_vec;
            if (nfd >= 4) {
                z_vec = fl.col(2);
                theta_vec = fl.col(3);
            } else {
                theta_vec = fl.col(0) + fl.col(1);
                z_vec = fl.col(0) - fl.col(1);
            }
            const double theta_norm = theta_vec.norm();
            const double z_norm = z_vec.norm();
            if (theta_norm < 1.0e-12 || z_norm < 1.0e-12) {
                throw std::runtime_error("[CANN] Unable to derive theta/z directions from the supplied fiber directions.");
            }
            theta_vec /= theta_norm;
            z_vec /= z_norm;

            Matrix<nsd> Htheta = theta_vec * theta_vec.transpose();
            Matrix<nsd> Hz = z_vec * z_vec.transpose();
            Inv[9] = 0.5 * (theta_vec.dot(C * theta_vec) - 1.0);
            Inv[10] = 0.5 * (z_vec.dot(C * z_vec) - 1.0);
            dInv10 = 0.5 * Htheta;
            dInv11 = 0.5 * Hz;
            Inv[17] = Inv[9] * Inv[10];
            dEthetaEz = Inv[10] * dInv10 + Inv[9] * dInv11;
            ddEthetaEz = dyadic_product<nsd>(dInv10, dInv11) + dyadic_product<nsd>(dInv11, dInv10);
        }
    }

    dInv = {dInv1, dInv2, dInv3, dInv4, dInv5, dInv6, dInv7, dInv8, dInv9, dInv10, dInv11,
        dE11, dE22, dE33, dE12, dE23, dE13, dEthetaEz};
    ddInv = {ddInv1, ddInv2, ddInv3, ddInv4, ddInv5, ddInv6, ddInv7, ddInv8, ddInv9, ddInv10, ddInv11,
        ddE11, ddE22, ddE33, ddE12, ddE23, ddE13, ddEthetaEz};
}


// Template instantiation
template void ArtificialNeuralNetMaterial::computeInvariantsAndDerivatives<2>(
const Matrix<2>& C, const Eigen::Matrix<double, 2, Eigen::Dynamic>& fl, int nfd, double J2d, double J4d, const Matrix<2>& Ci,
const Matrix<2>& Idm, const double Tfa, const double kap, Matrix<2>& N1, double& psi, double (&Inv)[NUM_INVARIANTS],
std::array<Matrix<2>,NUM_INVARIANTS>& dInv, std::array<Tensor<2>,NUM_INVARIANTS>& ddInv) const;

template void ArtificialNeuralNetMaterial::computeInvariantsAndDerivatives<3>(
const Matrix<3>& C, const Eigen::Matrix<double, 3, Eigen::Dynamic>& fl, int nfd, double J2d, double J4d, const Matrix<3>& Ci,
const Matrix<3>& Idm, const double Tfa, const double kap, Matrix<3>& N1, double& psi, double (&Inv)[NUM_INVARIANTS],
std::array<Matrix<3>,NUM_INVARIANTS>& dInv, std::array<Tensor<3>,NUM_INVARIANTS>& ddInv) const;

template void ArtificialNeuralNetMaterial::recruitedFiberInput<2>(
const int row, const double i4_eff_m1, const Matrix<2>& di4_eff_m1, const Tensor<2>& ddi4_eff_m1,
double& x_rec, Matrix<2>& dx_rec, Tensor<2>& ddx_rec) const;

template void ArtificialNeuralNetMaterial::recruitedFiberInput<3>(
const int row, const double i4_eff_m1, const Matrix<3>& di4_eff_m1, const Tensor<3>& ddi4_eff_m1,
double& x_rec, Matrix<3>& dx_rec, Tensor<3>& ddx_rec) const;
