// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MATERIAL_MODEL_H
#define MATERIAL_MODEL_H

#include "consts.h"
#include "mat_fun.h"
#include "ComMod.h"
#include "Parameters.h"

#include "eigen3/Eigen/Core"
#include "eigen3/unsupported/Eigen/CXX11/Tensor"

#include <memory>
#include <string>

/// @brief Abstract base class for object-oriented material models
/// Provides a clean alternative to the existing switch-based approach
template<size_t nsd>
class MaterialModel {
public:
    using Matrix = Eigen::Matrix<double, nsd, nsd>;
    using Tensor = Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>;

    virtual ~MaterialModel() = default;
    
    /// @brief Compute 2nd Piola-Kirchhoff stress and material tangent
    /// This follows the same computational pattern as existing mat_models.cpp functions
    /// @param[in] C Right Cauchy-Green deformation tensor
    /// @param[in] J Jacobian of deformation gradient  
    /// @param[in] J2d Isochoric scaling factor J^(-2/3)
    /// @param[in] Ci Isochoric right Cauchy-Green tensor  
    /// @param[in] Idm Identity matrix
    /// @param[in] nfd Number of fiber directions
    /// @param[in] fl Fiber directions
    /// @param[in] Tfa Active stress contribution
    /// @param[out] S 2nd Piola-Kirchhoff stress tensor
    /// @param[out] CC Material tangent modulus tensor
    virtual void compute_pk2cc(const Matrix& C, double J, double J2d, const Matrix& Ci,
                              const Matrix& Idm, int nfd, 
                              const Eigen::Matrix<double, nsd, Eigen::Dynamic>& fl,
                              const Matrix& Tfa, Matrix& S, Tensor& CC) = 0;
    
    /// @brief Get material model type identifier
    virtual consts::ConstitutiveModelType get_type() const = 0;
    
    /// @brief Get human-readable material name
    virtual std::string get_name() const = 0;
    
    /// @brief Get parameter value by name (for introspection and mixtures)
    virtual double get_parameter(const std::string& name) const { return 0.0; }
    
    /// @brief Set parameter value by name (for dynamic parameter updates)
    virtual void set_parameter(const std::string& name, double value) {}
    
    /// @brief Factory method to create material model instances
    /// This will allow gradual migration from existing parameter structures
    static std::unique_ptr<MaterialModel<nsd>> create_from_stM(const stModelType& stM);
    
    /// @brief Factory method to create from domain parameters (for new code)
    static std::unique_ptr<MaterialModel<nsd>> create_from_domain_params(
        consts::ConstitutiveModelType type,
        const DomainParameters* domain_params,
        double mu, double kap, double lam);
};

/// @brief Neo-Hookean material model - OOP implementation
/// Replicates the behavior of the existing Neo-Hookean switch case
template<size_t nsd>
class NeoHookeanOOP : public MaterialModel<nsd> {
private:
    double C10_;  ///< Neo-Hookean material parameter (μ/2)
    
public:
    using Matrix = typename MaterialModel<nsd>::Matrix;
    using Tensor = typename MaterialModel<nsd>::Tensor;
    
    /// @brief Constructor from C10 parameter
    explicit NeoHookeanOOP(double C10) : C10_(C10) {}
    
    /// @brief Compute Neo-Hookean stress and tangent
    /// Replicates existing mat_models.cpp logic exactly
    void compute_pk2cc(const Matrix& C, double J, double J2d, const Matrix& Ci,
                      const Matrix& Idm, int nfd, 
                      const Eigen::Matrix<double, nsd, Eigen::Dynamic>& fl,
                      const Matrix& Tfa, Matrix& S, Tensor& CC) override {
        
        using namespace mat_fun;
        
        // Neo-Hookean fictitious stress (matches existing implementation)
        Matrix S_bar = 2.0 * C10_ * Idm;
        
        // Add fiber reinforcement/active stress
        S_bar += Tfa;
        
        // Neo-Hookean has zero isochoric tangent modulus  
        Tensor CC_bar;
        CC_bar.setZero();
        
        // Transform to isochoric form using existing utility
        auto [S_iso, CC_iso] = bar_to_iso<nsd>(S_bar, CC_bar, J2d, C, Ci);
        S = S_iso;
        CC = CC_iso;
    }
    
    consts::ConstitutiveModelType get_type() const override {
        return consts::ConstitutiveModelType::stIso_nHook;
    }
    
    std::string get_name() const override {
        return "Neo-Hookean (OOP)";
    }
    
    double get_parameter(const std::string& name) const override {
        if (name == "C10") return C10_;
        return 0.0;
    }
    
    void set_parameter(const std::string& name, double value) override {
        if (name == "C10") C10_ = value;
    }
};

/// @brief Mooney-Rivlin material model - OOP implementation
/// Replicates the behavior of the existing Mooney-Rivlin switch case
template<size_t nsd>
class MooneyRivlinOOP : public MaterialModel<nsd> {
private:
    double C10_;  ///< First Mooney-Rivlin parameter
    double C01_;  ///< Second Mooney-Rivlin parameter
    
public:
    using Matrix = typename MaterialModel<nsd>::Matrix;
    using Tensor = typename MaterialModel<nsd>::Tensor;
    
    /// @brief Constructor from Mooney-Rivlin parameters
    MooneyRivlinOOP(double C10, double C01) : C10_(C10), C01_(C01) {}
    
    /// @brief Compute Mooney-Rivlin stress and tangent
    /// Replicates existing mat_models.cpp logic exactly
    void compute_pk2cc(const Matrix& C, double J, double J2d, const Matrix& Ci,
                      const Matrix& Idm, int nfd, 
                      const Eigen::Matrix<double, nsd, Eigen::Dynamic>& fl,
                      const Matrix& Tfa, Matrix& S, Tensor& CC) override {
        
        using namespace mat_fun;
        
        // Compute isochoric invariants (matches existing implementation)
        double Inv1 = Ci.trace();
        double J4d = std::pow(J, -4.0/3.0);
        
        // Mooney-Rivlin fictitious stress (matches existing implementation)
        Matrix S_bar = 2.0 * (C10_ + Inv1 * C01_) * Idm - 2.0 * C01_ * J2d * C;
        
        // Add fiber reinforcement/active stress
        S_bar += Tfa;
        
        // Mooney-Rivlin tangent modulus (matches existing implementation)
        Tensor CC_bar = 4.0 * J4d * C01_ * (dyadic_product<nsd>(Idm, Idm) - fourth_order_identity<nsd>());
        
        // Transform to isochoric form using existing utility
        auto [S_iso, CC_iso] = bar_to_iso<nsd>(S_bar, CC_bar, J2d, C, Ci);
        S = S_iso;
        CC = CC_iso;
    }
    
    consts::ConstitutiveModelType get_type() const override {
        return consts::ConstitutiveModelType::stIso_MR;
    }
    
    std::string get_name() const override {
        return "Mooney-Rivlin (OOP)";
    }
    
    double get_parameter(const std::string& name) const override {
        if (name == "C10" || name == "c1") return C10_;
        if (name == "C01" || name == "c2") return C01_;
        return 0.0;
    }
    
    void set_parameter(const std::string& name, double value) override {
        if (name == "C10" || name == "c1") C10_ = value;
        if (name == "C01" || name == "c2") C01_ = value;
    }
};

/// @brief Linear elastic material model - OOP implementation  
/// Replicates the behavior of the existing linear switch case
template<size_t nsd>
class LinearOOP : public MaterialModel<nsd> {
private:
    double mu_;  ///< Shear modulus (stored in C10 field in existing code)
    
public:
    using Matrix = typename MaterialModel<nsd>::Matrix;
    using Tensor = typename MaterialModel<nsd>::Tensor;
    
    /// @brief Constructor from shear modulus
    explicit LinearOOP(double mu) : mu_(mu) {}
    
    /// @brief Compute linear elastic stress and tangent
    /// Replicates existing mat_models.cpp logic exactly
    void compute_pk2cc(const Matrix& C, double J, double J2d, const Matrix& Ci,
                      const Matrix& Idm, int nfd, 
                      const Eigen::Matrix<double, nsd, Eigen::Dynamic>& fl,
                      const Matrix& Tfa, Matrix& S, Tensor& CC) override {
        
        using namespace mat_fun;
        
        // Linear model: S = μ * I (matches existing implementation)
        Matrix S_bar = mu_ * Idm;
        
        // Add fiber reinforcement/active stress
        S_bar += Tfa;
        
        // Linear model has zero tangent modulus
        Tensor CC_bar;
        CC_bar.setZero();
        
        // Transform to isochoric form using existing utility
        auto [S_iso, CC_iso] = bar_to_iso<nsd>(S_bar, CC_bar, J2d, C, Ci);
        S = S_iso;
        CC = CC_iso;
    }
    
    consts::ConstitutiveModelType get_type() const override {
        return consts::ConstitutiveModelType::stIso_lin;
    }
    
    std::string get_name() const override {
        return "Linear (OOP)";
    }
    
    double get_parameter(const std::string& name) const override {
        if (name == "C10" || name == "mu") return mu_;
        return 0.0;
    }
    
    void set_parameter(const std::string& name, double value) override {
        if (name == "C10" || name == "mu") mu_ = value;
    }
};

/// @brief Factory implementation - creates OOP materials from existing stModelType
/// This enables gradual migration from existing code
template<size_t nsd>
std::unique_ptr<MaterialModel<nsd>> MaterialModel<nsd>::create_from_stM(const stModelType& stM) {
    
    switch (stM.isoType) {
        case consts::ConstitutiveModelType::stIso_nHook: {
            return std::make_unique<NeoHookeanOOP<nsd>>(stM.C10);
        }
        
        case consts::ConstitutiveModelType::stIso_lin: {
            return std::make_unique<LinearOOP<nsd>>(stM.C10);
        }
        
        case consts::ConstitutiveModelType::stIso_MR: {
            return std::make_unique<MooneyRivlinOOP<nsd>>(stM.C10, stM.C01);
        }
        
        default:
            // Return nullptr for unsupported models (allows fallback to existing switch)
            return nullptr;
    }
}

/// @brief Factory implementation - creates OOP materials from domain parameters
/// This is for new code that can use OOP materials from the start
template<size_t nsd>
std::unique_ptr<MaterialModel<nsd>> MaterialModel<nsd>::create_from_domain_params(
    consts::ConstitutiveModelType type,
    const DomainParameters* domain_params,
    double mu, double kap, double lam) {
    
    switch (type) {
        case consts::ConstitutiveModelType::stIso_nHook: {
            double C10 = 0.5 * mu;  // Standard Neo-Hookean relationship
            return std::make_unique<NeoHookeanOOP<nsd>>(C10);
        }
        
        case consts::ConstitutiveModelType::stIso_lin: {
            return std::make_unique<LinearOOP<nsd>>(mu);
        }
        
        case consts::ConstitutiveModelType::stIso_MR: {
            if (domain_params && domain_params->constitutive_model.mooney_rivlin.c1.defined()) {
                double c1 = domain_params->constitutive_model.mooney_rivlin.c1.value();
                double c2 = domain_params->constitutive_model.mooney_rivlin.c2.value();
                return std::make_unique<MooneyRivlinOOP<nsd>>(c1, c2);
            } else {
                // Fallback to default parameters
                double c1 = 0.3 * mu;
                double c2 = 0.2 * mu;
                return std::make_unique<MooneyRivlinOOP<nsd>>(c1, c2);
            }
        }
        
        default:
            // Return nullptr for unsupported models
            return nullptr;
    }
}

#endif // MATERIAL_MODEL_H