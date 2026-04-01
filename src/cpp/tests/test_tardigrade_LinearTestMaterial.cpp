/**
 * \file test_tardigrade_LinearTestMaterial.cpp
 *
 * Tests for tardigrade_cmml
 */

#include <fstream>
#include <iostream>
#include <sstream>

#include "tardigrade_LinearTestMaterial.h"

#define BOOST_TEST_MODULE test_tardigrade_LinearTestMaterial
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element()

struct cout_redirect {
    cout_redirect(std::streambuf *new_buffer) : old(std::cout.rdbuf(new_buffer)) {}

    ~cout_redirect() { std::cout.rdbuf(old); }

   private:
    std::streambuf *old;
};

struct cerr_redirect {
    cerr_redirect(std::streambuf *new_buffer) : old(std::cerr.rdbuf(new_buffer)) {}

    ~cerr_redirect() { std::cerr.rdbuf(old); }

   private:
    std::streambuf *old;
};

BOOST_AUTO_TEST_CASE(test_LinearTestMaterial, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the linear hydra test model Jacobian
     */

    double              current_time = 1.23;
    double              dt           = 0.234;
    std::vector<double> parameters   = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

    std::vector<double> dof          = {0.1, 0.2, 0.3, 0.4};
    std::vector<double> previous_dof = {-0.1, -0.2, -0.3, -0.4};

    std::vector<double> answer = {3, 7, 11};

    std::vector<double> sdvs(0, 0);

    std::vector<double> result   = { 1, 1, 1 };
    std::vector<double> result_J = { 1, 1, 1 };

    std::vector<double> jacobian(3 * 4, 0);
    std::vector<double> additional(0, 0);
    std::string         output_message;

    tardigradeCMML::LinearTestMaterial::LinearTestMaterial model;
    tardigradeCMML::LinearTestMaterial::LinearTestMaterial model_J;

    int error_code =
        model.evaluate_model(current_time, dt, dof.data(), previous_dof.data(), dof.size(), parameters.data(),
                             parameters.size(), sdvs.data(), sdvs.size(), result.data(), result.size(), output_message);

    if (error_code != 0) {
        std::cerr << output_message << "\n";
        BOOST_CHECK(false);
    }

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    error_code =
        model_J.evaluate_model(current_time, dt, dof.data(), previous_dof.data(), dof.size(), parameters.data(),
                               parameters.size(), sdvs.data(), sdvs.size(), result_J.data(), result_J.size(),
                               jacobian.data(), additional.data(), additional.size(), output_message);

    if (error_code != 0) {
        std::cerr << output_message << "\n";
        BOOST_CHECK(false);
    }

    BOOST_TEST(answer == result_J, CHECK_PER_ELEMENT);

    // Check the Jacobian
    double eps = 1e-6;
    {
        constexpr unsigned int vardim = 4;
        constexpr unsigned int outdim = 3;
        std::vector<double>    x      = dof;
        std::vector<double>    jacobian_answer(vardim * outdim, 0);

        for (unsigned int i = 0; i < vardim; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            tardigradeCMML::LinearTestMaterial::LinearTestMaterial modelp;
            tardigradeCMML::LinearTestMaterial::LinearTestMaterial modelm;

            std::vector<double> rp(outdim, 0);
            std::vector<double> rm(outdim, 0);

            int ep =
                model.evaluate_model(current_time, dt, xp.data(), previous_dof.data(), dof.size(), parameters.data(),
                                     parameters.size(), sdvs.data(), sdvs.size(), rp.data(), rp.size(), output_message);
            int em =
                model.evaluate_model(current_time, dt, xm.data(), previous_dof.data(), dof.size(), parameters.data(),
                                     parameters.size(), sdvs.data(), sdvs.size(), rm.data(), rm.size(), output_message);

            if (ep != 0) {
                BOOST_CHECK(false);
            }
            if (em != 0) {
                BOOST_CHECK(false);
            }

            for (unsigned int j = 0; j < outdim; ++j) {
                jacobian_answer[vardim * j + i] = (rp[j] - rm[j]) / (2 * delta);
            }
        }
        BOOST_TEST(jacobian_answer == jacobian, CHECK_PER_ELEMENT);
    }
}
