/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 */

#include "../test_common.h"
#include "Core/Exception.h"
#include "Parameters.h"

#include <string>

namespace {

tinyxml2::XMLElement* parse_cann_xml(tinyxml2::XMLDocument& document, const std::string& body)
{
    const std::string xml = "<Constitutive_model type=\"CANN\">" + body + "</Constitutive_model>";
    EXPECT_EQ(document.Parse(xml.c_str()), tinyxml2::XML_SUCCESS);
    return document.FirstChildElement("Constitutive_model");
}

void expect_cann_parse_error(const std::string& body)
{
    tinyxml2::XMLDocument document;
    auto* element = parse_cann_xml(document, body);
    ASSERT_NE(element, nullptr);

    CANNParameters params;
    EXPECT_THROW(params.set_values(element), svmp::ParseException);
}

}

TEST(CANNLegacyParserTest, ParsesLegacyRowDispersion)
{
    tinyxml2::XMLDocument document;
    auto* element = parse_cann_xml(document,
        "<Add_row row_name=\"hgo_D1\">"
        "  <Invariant_num> 4 </Invariant_num>"
        "  <Activation_functions> (1, 2, 2) </Activation_functions>"
        "  <Weights> (1.0, 10.0, 2500.0) </Weights>"
        "  <Dispersion> 0.0231 </Dispersion>"
        "</Add_row>");
    ASSERT_NE(element, nullptr);

    CANNParameters params;
    params.set_values(element);

    ASSERT_EQ(params.rows.size(), 1);
    EXPECT_DOUBLE_EQ(params.rows[0]->row.dispersion.value(), 0.0231);
}

TEST(CANNLegacyParserTest, ParsesLegacyRowBetaRecruitment)
{
    tinyxml2::XMLDocument document;
    auto* element = parse_cann_xml(document,
        "<Add_row row_name=\"recruited_D1\">"
        "  <Invariant_num> 4 </Invariant_num>"
        "  <Activation_functions> (1, 2, 2) </Activation_functions>"
        "  <Weights> (1.0, 10.0, 2500.0) </Weights>"
        "  <Recruitment distribution=\"beta\">"
        "    <Lower_stretch> 1.0 </Lower_stretch>"
        "    <Upper_stretch> 1.3 </Upper_stretch>"
        "    <Tau> 0.02 </Tau>"
        "    <Alpha> 2.0 </Alpha>"
        "    <Beta> 5.0 </Beta>"
        "    <Quadrature_points> 32 </Quadrature_points>"
        "  </Recruitment>"
        "</Add_row>");
    ASSERT_NE(element, nullptr);

    CANNParameters params;
    params.set_values(element);

    ASSERT_EQ(params.rows.size(), 1);
    EXPECT_TRUE(params.rows[0]->row.recruitment_enabled);
    EXPECT_DOUBLE_EQ(params.rows[0]->row.recruitment_lower_stretch, 1.0);
    EXPECT_DOUBLE_EQ(params.rows[0]->row.recruitment_upper_stretch, 1.3);
    EXPECT_DOUBLE_EQ(params.rows[0]->row.recruitment_tau, 0.02);
    EXPECT_DOUBLE_EQ(params.rows[0]->row.recruitment_alpha, 2.0);
    EXPECT_DOUBLE_EQ(params.rows[0]->row.recruitment_beta, 5.0);
    EXPECT_EQ(params.rows[0]->row.recruitment_quadrature_points, 32);
}

TEST(CANNLegacyParserTest, ParsesLegacyRowBetaRecruitmentDefaultQuadrature)
{
    tinyxml2::XMLDocument document;
    auto* element = parse_cann_xml(document,
        "<Add_row row_name=\"recruited_D2\">"
        "  <Invariant_num> 8 </Invariant_num>"
        "  <Activation_functions> (1, 2, 2) </Activation_functions>"
        "  <Weights> (1.0, 10.0, 2500.0) </Weights>"
        "  <Recruitment distribution=\"beta\">"
        "    <Lower_stretch> 1.0 </Lower_stretch>"
        "    <Upper_stretch> 1.3 </Upper_stretch>"
        "    <Tau> 0.02 </Tau>"
        "    <Alpha> 2.0 </Alpha>"
        "    <Beta> 5.0 </Beta>"
        "  </Recruitment>"
        "</Add_row>");
    ASSERT_NE(element, nullptr);

    CANNParameters params;
    params.set_values(element);

    ASSERT_EQ(params.rows.size(), 1);
    EXPECT_TRUE(params.rows[0]->row.recruitment_enabled);
    EXPECT_EQ(params.rows[0]->row.recruitment_quadrature_points, 32);
}

TEST(CANNLegacyParserTest, ParsesLegacyThetaFiberDispersionAndRecruitment)
{
    tinyxml2::XMLDocument document;
    auto* element = parse_cann_xml(document,
        "<Add_row row_name=\"theta_recruited\">"
        "  <Invariant_num> 10 </Invariant_num>"
        "  <Activation_functions> (1, 3, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 2500.0) </Weights>"
        "  <Dispersion> 0.0231 </Dispersion>"
        "  <Recruitment distribution=\"beta\">"
        "    <Lower_stretch> 1.0 </Lower_stretch>"
        "    <Upper_stretch> 1.3 </Upper_stretch>"
        "    <Tau> 0.02 </Tau>"
        "    <Alpha> 2.0 </Alpha>"
        "    <Beta> 5.0 </Beta>"
        "  </Recruitment>"
        "</Add_row>");
    ASSERT_NE(element, nullptr);

    CANNParameters params;
    params.set_values(element);

    ASSERT_EQ(params.rows.size(), 1);
    EXPECT_EQ(params.rows[0]->row.invariant_index.value(), 10);
    EXPECT_DOUBLE_EQ(params.rows[0]->row.dispersion.value(), 0.0231);
    EXPECT_TRUE(params.rows[0]->row.recruitment_enabled);
    EXPECT_EQ(params.rows[0]->row.recruitment_quadrature_points, 32);
}

TEST(CANNLegacyParserTest, ParsesLegacyZFiberRecruitment)
{
    tinyxml2::XMLDocument document;
    auto* element = parse_cann_xml(document,
        "<Add_row row_name=\"z_recruited\">"
        "  <Invariant_num> 11 </Invariant_num>"
        "  <Activation_functions> (1, 1, 2) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 0.002) </Weights>"
        "  <Recruitment distribution=\"beta\">"
        "    <Lower_stretch> 0.95 </Lower_stretch>"
        "    <Upper_stretch> 1.8 </Upper_stretch>"
        "    <Tau> 0.02 </Tau>"
        "    <Alpha> 2.0 </Alpha>"
        "    <Beta> 5.0 </Beta>"
        "  </Recruitment>"
        "</Add_row>");
    ASSERT_NE(element, nullptr);

    CANNParameters params;
    params.set_values(element);

    ASSERT_EQ(params.rows.size(), 1);
    EXPECT_EQ(params.rows[0]->row.invariant_index.value(), 11);
    EXPECT_TRUE(params.rows[0]->row.recruitment_enabled);
}

TEST(CANNLegacyParserTest, ParsesLegacyCubicH1Activation)
{
    tinyxml2::XMLDocument document;
    auto* element = parse_cann_xml(document,
        "<Add_row row_name=\"cubic_I1\">"
        "  <Invariant_num> 1 </Invariant_num>"
        "  <Activation_functions> (1, 3, 1) </Activation_functions>"
        "  <Weights> (0.75, 1.0, 1250.0) </Weights>"
        "</Add_row>");
    ASSERT_NE(element, nullptr);

    CANNParameters params;
    params.set_values(element);

    ASSERT_EQ(params.rows.size(), 1);
    EXPECT_EQ(params.rows[0]->row.activation_functions.value()[1], 3);
}

TEST(CANNLegacyParserTest, ParsesLegacyGreenStrainRows)
{
    tinyxml2::XMLDocument document;
    auto* element = parse_cann_xml(document,
        "<Add_row row_name=\"E11\">"
        "  <Invariant_num> 12 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 2.0) </Weights>"
        "</Add_row>"
        "<Add_row row_name=\"E12\">"
        "  <Invariant_num> 15 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 3.0) </Weights>"
        "</Add_row>"
        "<Add_row row_name=\"EthetaEz\">"
        "  <Invariant_num> 18 </Invariant_num>"
        "  <Activation_functions> (1, 3, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 4.0) </Weights>"
        "</Add_row>");
    ASSERT_NE(element, nullptr);

    CANNParameters params;
    params.set_values(element);

    ASSERT_EQ(params.rows.size(), 3);
    EXPECT_EQ(params.rows[0]->row.invariant_index.value(), 12);
    EXPECT_EQ(params.rows[1]->row.invariant_index.value(), 15);
    EXPECT_EQ(params.rows[2]->row.invariant_index.value(), 18);
    EXPECT_EQ(params.rows[2]->row.activation_functions.value()[1], 3);
}

TEST(CANNLegacyParserTest, ParsesLegacyEthetaAndEzRows)
{
    tinyxml2::XMLDocument document;
    auto* element = parse_cann_xml(document,
        "<Add_row row_name=\"Etheta\">"
        "  <Invariant_num> 10 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 1.0) </Weights>"
        "</Add_row>"
        "<Add_row row_name=\"Ez\">"
        "  <Invariant_num> 11 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 1.0) </Weights>"
        "</Add_row>");
    ASSERT_NE(element, nullptr);

    CANNParameters params;
    params.set_values(element);

    ASSERT_EQ(params.rows.size(), 2);
    EXPECT_EQ(params.rows[0]->row.invariant_index.value(), 10);
    EXPECT_EQ(params.rows[1]->row.invariant_index.value(), 11);
}

TEST(CANNLegacyParserTest, RejectsLegacyInvariantOutsideSupportedRange)
{
    expect_cann_parse_error(
        "<Add_row row_name=\"bad\">"
        "  <Invariant_num> 19 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 1.0) </Weights>"
        "</Add_row>");
}

TEST(CANNLegacyParserTest, RejectsLegacyRowDispersionOutsideGohRange)
{
    expect_cann_parse_error(
        "<Add_row row_name=\"hgo_D1\">"
        "  <Invariant_num> 4 </Invariant_num>"
        "  <Activation_functions> (1, 2, 2) </Activation_functions>"
        "  <Weights> (1.0, 10.0, 2500.0) </Weights>"
        "  <Dispersion> 0.5 </Dispersion>"
        "</Add_row>");
}

TEST(CANNLegacyParserTest, RejectsLegacyRowDispersionOnNonI4Invariant)
{
    expect_cann_parse_error(
        "<Add_row row_name=\"iso_I1\">"
        "  <Invariant_num> 1 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 5000.0) </Weights>"
        "  <Dispersion> 0.0231 </Dispersion>"
        "</Add_row>");
}

TEST(CANNLegacyParserTest, RejectsLegacyRowRecruitmentUnsupportedDistribution)
{
    expect_cann_parse_error(
        "<Add_row row_name=\"recruited_D1\">"
        "  <Invariant_num> 4 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 1.0) </Weights>"
        "  <Recruitment distribution=\"lognormal\">"
        "    <Lower_stretch> 1.0 </Lower_stretch>"
        "    <Upper_stretch> 1.3 </Upper_stretch>"
        "    <Tau> 0.02 </Tau>"
        "    <Alpha> 2.0 </Alpha>"
        "    <Beta> 5.0 </Beta>"
        "  </Recruitment>"
        "</Add_row>");
}

TEST(CANNLegacyParserTest, RejectsLegacyRowRecruitmentInvalidParameters)
{
    expect_cann_parse_error(
        "<Add_row row_name=\"bad_bounds\">"
        "  <Invariant_num> 4 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 1.0) </Weights>"
        "  <Recruitment distribution=\"beta\">"
        "    <Lower_stretch> 1.3 </Lower_stretch>"
        "    <Upper_stretch> 1.0 </Upper_stretch>"
        "    <Tau> 0.02 </Tau>"
        "    <Alpha> 2.0 </Alpha>"
        "    <Beta> 5.0 </Beta>"
        "  </Recruitment>"
        "</Add_row>");
    expect_cann_parse_error(
        "<Add_row row_name=\"bad_tau\">"
        "  <Invariant_num> 4 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 1.0) </Weights>"
        "  <Recruitment distribution=\"beta\">"
        "    <Lower_stretch> 1.0 </Lower_stretch>"
        "    <Upper_stretch> 1.3 </Upper_stretch>"
        "    <Tau> 0.0 </Tau>"
        "    <Alpha> 2.0 </Alpha>"
        "    <Beta> 5.0 </Beta>"
        "  </Recruitment>"
        "</Add_row>");
    expect_cann_parse_error(
        "<Add_row row_name=\"bad_shape\">"
        "  <Invariant_num> 4 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 1.0) </Weights>"
        "  <Recruitment distribution=\"beta\">"
        "    <Lower_stretch> 1.0 </Lower_stretch>"
        "    <Upper_stretch> 1.3 </Upper_stretch>"
        "    <Tau> 0.02 </Tau>"
        "    <Alpha> 0.0 </Alpha>"
        "    <Beta> 5.0 </Beta>"
        "  </Recruitment>"
        "</Add_row>");
    expect_cann_parse_error(
        "<Add_row row_name=\"bad_quad\">"
        "  <Invariant_num> 4 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 1.0) </Weights>"
        "  <Recruitment distribution=\"beta\">"
        "    <Lower_stretch> 1.0 </Lower_stretch>"
        "    <Upper_stretch> 1.3 </Upper_stretch>"
        "    <Tau> 0.02 </Tau>"
        "    <Alpha> 2.0 </Alpha>"
        "    <Beta> 5.0 </Beta>"
        "    <Quadrature_points> 3 </Quadrature_points>"
        "  </Recruitment>"
        "</Add_row>");
}

TEST(CANNLegacyParserTest, RejectsLegacyRowRecruitmentOnNonI4Invariant)
{
    expect_cann_parse_error(
        "<Add_row row_name=\"iso_I1\">"
        "  <Invariant_num> 1 </Invariant_num>"
        "  <Activation_functions> (1, 1, 1) </Activation_functions>"
        "  <Weights> (1.0, 1.0, 1.0) </Weights>"
        "  <Recruitment distribution=\"beta\">"
        "    <Lower_stretch> 1.0 </Lower_stretch>"
        "    <Upper_stretch> 1.3 </Upper_stretch>"
        "    <Tau> 0.02 </Tau>"
        "    <Alpha> 2.0 </Alpha>"
        "    <Beta> 5.0 </Beta>"
        "  </Recruitment>"
        "</Add_row>");
}
