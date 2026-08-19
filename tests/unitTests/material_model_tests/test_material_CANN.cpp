/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 */

#include "test_material_common.h"
#include "ArtificialNeuralNetMaterial.h"

#include <array>
#include <cmath>

namespace {

struct TestRecruitmentParams {
    bool enabled = false;
    double lower = 1.0;
    double upper = 1.3;
    double tau = 0.02;
    double alpha = 2.0;
    double beta = 5.0;
    int quadrature_points = 32;
};

class TestCANNLegacySingleRow : public TestMaterialModel {
public:
    TestCANNLegacySingleRow(const int invariant_index, const std::array<int,3>& activation_functions,
                            const std::array<double,3>& weights,
                            const std::array<double,3>& d1, const std::array<double,3>& d2,
                            const double dispersion = 0.0,
                            const TestRecruitmentParams recruitment = TestRecruitmentParams{})
      : TestMaterialModel(consts::ConstitutiveModelType::stArtificialNeuralNet,
                          consts::ConstitutiveModelType::stVol_ST91)
    {
        auto& dmn = com_mod.mockEq.mockDmn;
        dmn.stM.paramTable.num_rows = 1;
        dmn.stM.paramTable.invariant_indices.resize(1);
        dmn.stM.paramTable.activation_functions.resize(1, 3);
        dmn.stM.paramTable.weights.resize(1, 3);
        dmn.stM.paramTable.row_dispersion.resize(1);
        dmn.stM.paramTable.row_recruitment_enabled.resize(1);
        dmn.stM.paramTable.row_recruitment_lower_stretch.resize(1);
        dmn.stM.paramTable.row_recruitment_upper_stretch.resize(1);
        dmn.stM.paramTable.row_recruitment_tau.resize(1);
        dmn.stM.paramTable.row_recruitment_alpha.resize(1);
        dmn.stM.paramTable.row_recruitment_beta.resize(1);
        dmn.stM.paramTable.row_recruitment_quadrature_points.resize(1);

        dmn.stM.paramTable.invariant_indices[0] = invariant_index;
        dmn.stM.paramTable.row_dispersion[0] = dispersion;
        dmn.stM.paramTable.row_recruitment_enabled[0] = recruitment.enabled ? 1 : 0;
        dmn.stM.paramTable.row_recruitment_lower_stretch[0] = recruitment.lower;
        dmn.stM.paramTable.row_recruitment_upper_stretch[0] = recruitment.upper;
        dmn.stM.paramTable.row_recruitment_tau[0] = recruitment.tau;
        dmn.stM.paramTable.row_recruitment_alpha[0] = recruitment.alpha;
        dmn.stM.paramTable.row_recruitment_beta[0] = recruitment.beta;
        dmn.stM.paramTable.row_recruitment_quadrature_points[0] = recruitment.quadrature_points;
        for (int i = 0; i < 3; ++i) {
            dmn.stM.paramTable.activation_functions(0, i) = activation_functions[i];
            dmn.stM.paramTable.weights(0, i) = weights[i];
        }
        dmn.stM.Kpen = 0.0;

        fN(0,0) = d1[0]; fN(1,0) = d1[1]; fN(2,0) = d1[2];
        fN(0,1) = d2[0]; fN(1,1) = d2[1]; fN(2,1) = d2[2];
    }

    void printMaterialParameters() override {}

    double computeStrainEnergy(const Array<double>&) override
    {
        return 0.0;
    }
};

Array<double> legacy_test_F()
{
    return {{1.35, 0.02, 0.00},
            {0.01, 1.15, 0.00},
            {0.00, 0.00, 1.0 / (1.35 * 1.15)}};
}

std::array<double,3> test_d1()
{
    const double angle = 35.0 * 3.14159265358979323846 / 180.0;
    return {std::cos(angle), std::sin(angle), 0.0};
}

std::array<double,3> test_d2()
{
    const double angle = 35.0 * 3.14159265358979323846 / 180.0;
    return {std::cos(angle), -std::sin(angle), 0.0};
}

std::array<double,3> cylindrical_family_direction(const double angle_degrees)
{
    const double angle = angle_degrees * 3.14159265358979323846 / 180.0;
    return {std::cos(angle), 0.0, std::sin(angle)};
}

}

TEST(CANNLegacyMaterialTest, LegacyCubicH1StressAndTangentMatchFiniteDifference)
{
    const auto F = legacy_test_F();
    ArtificialNeuralNetMaterial cann;
    double psi = 0.0;
    double dpsi = 0.0;
    double ddpsi = 0.0;
    cann.uCANN_scalar(1.2, 1, 3, 1, 0.75, 1.0, 1250.0, psi, dpsi, ddpsi);
    EXPECT_NEAR(psi, 1250.0 * std::pow(0.75 * 1.2, 3), 1.0e-12);
    EXPECT_NEAR(dpsi, 1250.0 * 3.0 * std::pow(0.75, 3) * std::pow(1.2, 2), 1.0e-12);
    EXPECT_NEAR(ddpsi, 1250.0 * 6.0 * std::pow(0.75, 3) * 1.2, 1.0e-12);

    TestCANNLegacySingleRow material(
        1,                    // invariant 1: I1
        {1, 3, 1},             // identity input, cubic hidden layer, linear output
        {0.75, 1.0, 1250.0},   // Psi = W2 * (W0 * (I1 - 3))^3
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0});
    material.ustruct = false;

    material.testMaterialElasticityConsistentWithPK2Stress(F, 5e-3, 1e-6, 1e-6, false);
}

TEST(CANNLegacyMaterialTest, LegacyGreenStrainInputsProduceExpectedStress)
{
    const auto F = legacy_test_F();
    const double weight = 42.0;
    TestCANNLegacySingleRow material(
        15,                  // invariant 15: E12
        {1, 1, 1},            // linear row
        {1.0, 1.0, weight},   // Psi = weight * E12
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0});
    material.ustruct = false;

    Array<double> S(3, 3), Dm(6, 6);
    material.compute_pk2cc(F, S, Dm);

    EXPECT_NEAR(S(0,0), 0.0, 1.0e-12);
    EXPECT_NEAR(S(1,1), 0.0, 1.0e-12);
    EXPECT_NEAR(S(2,2), 0.0, 1.0e-12);
    EXPECT_NEAR(S(0,1), 0.5 * weight, 1.0e-12);
    EXPECT_NEAR(S(1,0), 0.5 * weight, 1.0e-12);
    material.testMaterialElasticityConsistentWithPK2Stress(F, 5e-3, 1e-6, 1e-6, false);
}

TEST(CANNLegacyMaterialTest, LegacyRecruitmentDisabledMatchesPlainFiberRow)
{
    const auto F = legacy_test_F();
    const double weight = 2.5e4;

    TestCANNLegacySingleRow plain(
        4,
        {2, 1, 1},
        {1.0, 1.0, weight},
        test_d1(),
        test_d2());
    plain.ustruct = false;

    TestCANNLegacySingleRow disabled(
        4,
        {2, 1, 1},
        {1.0, 1.0, weight},
        test_d1(),
        test_d2(),
        0.0,
        TestRecruitmentParams{});
    disabled.ustruct = false;

    Array<double> S_plain(3, 3), Dm_plain(6, 6);
    Array<double> S_disabled(3, 3), Dm_disabled(6, 6);
    plain.compute_pk2cc(F, S_plain, Dm_plain);
    disabled.compute_pk2cc(F, S_disabled, Dm_disabled);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(S_plain(i,j), S_disabled(i,j), 1.0e-8 * std::max(1.0, std::abs(S_disabled(i,j))));
        }
    }
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            EXPECT_NEAR(Dm_plain(i,j), Dm_disabled(i,j), 1.0e-8 * std::max(1.0, std::abs(Dm_disabled(i,j))));
        }
    }
}

TEST(CANNLegacyMaterialTest, LegacyRecruitedI4D1StressAndTangentMatchFiniteDifference)
{
    const auto F = legacy_test_F();
    const TestRecruitmentParams recruitment{true, 1.0, 1.35, 0.05, 2.0, 5.0, 24};
    TestCANNLegacySingleRow material(
        4,
        {1, 2, 2},
        {1.0, 4.0, 2.0e3},
        test_d1(),
        test_d2(),
        0.0,
        recruitment);
    material.ustruct = false;

    material.testMaterialElasticityConsistentWithPK2Stress(F, 5e-3, 1e-6, 1e-6, false);
}

TEST(CANNLegacyMaterialTest, LegacyRecruitedI4D2StressAndTangentMatchFiniteDifference)
{
    const auto F = legacy_test_F();
    const TestRecruitmentParams recruitment{true, 1.0, 1.35, 0.05, 2.0, 5.0, 24};
    TestCANNLegacySingleRow material(
        8,
        {1, 2, 2},
        {1.0, 4.0, 2.0e3},
        test_d1(),
        test_d2(),
        0.0,
        recruitment);
    material.ustruct = false;

    material.testMaterialElasticityConsistentWithPK2Stress(F, 5e-3, 1e-6, 1e-6, false);
}

TEST(CANNLegacyMaterialTest, LegacyDispersedRecruitedI4StressAndTangentMatchFiniteDifference)
{
    const auto F = legacy_test_F();
    const TestRecruitmentParams recruitment{true, 1.0, 1.35, 0.05, 2.0, 5.0, 24};
    TestCANNLegacySingleRow material(
        4,
        {1, 2, 2},
        {1.0, 4.0, 2.0e3},
        test_d1(),
        test_d2(),
        0.0231,
        recruitment);
    material.ustruct = false;

    material.testMaterialElasticityConsistentWithPK2Stress(F, 5e-3, 1e-6, 1e-6, false);
}

TEST(CANNLegacyMaterialTest, LegacyDispersedRecruitedThetaFiberStressAndTangentMatchFiniteDifference)
{
    const auto F = legacy_test_F();
    const TestRecruitmentParams recruitment{true, 1.0, 1.35, 0.05, 2.0, 5.0, 24};
    TestCANNLegacySingleRow material(
        10,                  // invariant 10: Etheta, transformed to I4theta - 1 for fiber recruitment
        {1, 3, 1},
        {1.0, 1.0, 2.0e3},   // Psi = 2000 * (I4theta_rec - 1)^3
        test_d1(),
        test_d2(),
        0.0231,
        recruitment);
    material.ustruct = false;

    material.testMaterialElasticityConsistentWithPK2Stress(F, 5e-3, 1e-6, 1e-6, false);
}

TEST(CANNLegacyMaterialTest, LegacyRecruitedZFiberStressAndTangentMatchFiniteDifference)
{
    const auto F = legacy_test_F();
    const TestRecruitmentParams recruitment{true, 0.95, 1.8, 0.05, 2.0, 5.0, 24};
    TestCANNLegacySingleRow material(
        11,                  // invariant 11: Ez, transformed to I4z - 1 for fiber recruitment
        {1, 1, 2},
        {1.0, 1.0, 2.0e-3},  // Psi = 0.002 * (exp(I4z_rec - 1) - 1)
        test_d1(),
        test_d2(),
        0.0,
        recruitment);
    material.ustruct = false;

    material.testMaterialElasticityConsistentWithPK2Stress(F, 5e-3, 1e-6, 1e-6, false);
}

TEST(CANNLegacyMaterialTest, LegacyFourFamilyFiberRowsUseDirectThetaAndZDirections)
{
    const auto F = legacy_test_F();
    const TestRecruitmentParams recruitment{true, 0.95, 1.8, 0.05, 2.0, 5.0, 24};
    TestCANNLegacySingleRow material(
        10,
        {1, 3, 1},
        {1.0, 1.0, 2.0e3},
        cylindrical_family_direction(28.807),
        cylindrical_family_direction(64.160),
        0.0,
        recruitment);
    material.nFn = 4;
    material.fN = Array<double>(3, 4);
    const auto f = cylindrical_family_direction(28.807);
    const auto s = cylindrical_family_direction(64.160);
    const auto el = cylindrical_family_direction(91.362);
    const auto smc = cylindrical_family_direction(-1.068);
    for (int i = 0; i < 3; ++i) {
        material.fN(i,0) = f[i];
        material.fN(i,1) = s[i];
        material.fN(i,2) = el[i];
        material.fN(i,3) = smc[i];
    }
    material.ustruct = false;

    material.testMaterialElasticityConsistentWithPK2Stress(F, 5e-3, 1e-6, 1e-6, false);
}

TEST(CANNLegacyMaterialTest, LegacyEthetaEzCubicH1StressAndTangentMatchFiniteDifference)
{
    const auto F = legacy_test_F();
    TestCANNLegacySingleRow material(
        18,                  // invariant 18: Etheta * Ez
        {1, 3, 1},            // cubic hidden layer on the cross term
        {1.0, 1.0, 1250.0},   // Psi = 1250 * (Etheta * Ez)^3
        {0.7071067811865475, 0.7071067811865475, 0.0},
        {0.7071067811865475, -0.7071067811865475, 0.0});
    material.ustruct = false;

    material.testMaterialElasticityConsistentWithPK2Stress(F, 5e-3, 1e-6, 1e-6, false);
}
