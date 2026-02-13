/**
 * \file test_tardigrade_BasicSolid.cpp
 *
 * Tests for tardigrade_cmml
 */

#include <fstream>
#include <iostream>
#include <sstream>

#include "tardigrade_BasicSolid.h"

#define BOOST_TEST_MODULE test_tardigrade_BasicSolid
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

BOOST_AUTO_TEST_CASE(test_basic_functionality, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the basic functionality of the BasicSolid class
     */

    class modelMock : public tardigradeCMML::BasicSolid::BasicSolid {
       public:
        modelMock() : tardigradeCMML::BasicSolid::BasicSolid() {}

        void public_extract_parameters(double *parameters_begin, const unsigned int parameters_size) {
            extract_parameters(parameters_begin, parameters_size);
        }
    };

    modelMock model;

    std::vector<double> parameters = {1, 2, 3, 4, 5, 6, 7};

    model.public_extract_parameters(parameters.data(), 7);

    BOOST_TEST(model.getEvaluateModelResultSize() == 23);

    BOOST_TEST(model.getIsCurrent());

    BOOST_TEST(model.getDisplacementGradientIndex() == 1);

    BOOST_TEST(model.getTemperatureIndex() == 2);
}

BOOST_AUTO_TEST_CASE(test_evaluate_model, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Test the evaluate_model function of the BasicSolid::BasicSolid class
     */

    class modelMock : public tardigradeCMML::BasicSolid::BasicSolid {
       public:
        modelMock() : tardigradeCMML::BasicSolid::BasicSolid() {}

        std::vector<double> F, Fp;

        void public_extract_parameters(double *parameters_begin, const unsigned int parameters_size) {
            extract_parameters(parameters_begin, parameters_size);
        }

       protected:
        virtual void formDeformationGradients(const double *current_dof_begin, const double *previous_dof_begin,
                                              std::vector<double> &deformationGradient,
                                              std::vector<double> &previousDeformationGradient) override {
            tardigradeCMML::BasicSolid::BasicSolid::formDeformationGradients(current_dof_begin, previous_dof_begin,
                                                                             deformationGradient,
                                                                             previousDeformationGradient);

            F  = deformationGradient;
            Fp = previousDeformationGradient;
        }

        virtual void formDeformationGradients(const double *current_dof_begin, const double *previous_dof_begin,
                                              std::vector<double> &deformationGradient,
                                              std::vector<double> &previousDeformationGradient,
                                              std::vector<double> &dFdGradU) override {
            tardigradeCMML::BasicSolid::BasicSolid::formDeformationGradients(current_dof_begin, previous_dof_begin,
                                                                             deformationGradient,
                                                                             previousDeformationGradient, dFdGradU);

            F  = deformationGradient;
            Fp = previousDeformationGradient;
        }
    };

    modelMock model;

    double current_time = 1.23;
    double dt           = 0.345;

    std::vector<double> current_dof(26, 0);
    std::vector<double> previous_dof(26, 0);
    for (unsigned int i = 0; i < 26; ++i) {
        current_dof[i] += 0.01 * (i + 1);
        previous_dof[i] -= 0.01 * (i + 1);
    }

    std::vector<double> parameters = {2, 3, 4, 5, 6, 7, 8};

    std::vector<double> sdvs(4, 0);

    std::vector<double> result(23, -1);

    std::string output_message;

    std::vector<double> F_answer = {1.03996448, 0.05150977, 0.06305506, 0.07688404, 1.08957117,
                                    0.10225831, 0.1138036,  0.12763258, 1.14146156};

    std::vector<double> Fp_answer = {0.97641119, -0.03261049, -0.04163218, -0.04916405, 0.94255918,
                                     -0.0657176, -0.07473928, -0.08227115, 0.91019699};

    std::vector<double> answer = {2.0707740, 1.2922023, 1.7429554, 1.2922023, 3.0650159, 2.37918, 1.7429554, 2.37918,
                                  4.24473,   0.28,      0,         0,         0,         0,       0,         0,
                                  0,         -0.4,      -0.48,     -0.56,     0,         0,       0};

    int error_code = model.evaluate_model(current_time, dt, current_dof.data(), previous_dof.data(), 26,
                                          parameters.data(), 7, sdvs.data(), 4, result.data(), 23, output_message);

    BOOST_TEST(error_code == 0);

    if (error_code != 0) {
        std::cout << "ERROR MESSAGE:\n" << output_message << "\n";
    }

    BOOST_TEST(F_answer == model.F, CHECK_PER_ELEMENT);

    BOOST_TEST(Fp_answer == model.Fp, CHECK_PER_ELEMENT);

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    // Test the jacobians
    sdvs   = std::vector<double>(4, 0);
    result = std::vector<double>(23, -1);
    std::vector<double> jacobian(23 * 26, 0);
    std::vector<double> additional(0, 0);

    modelMock modelJ;
    error_code =
        modelJ.evaluate_model(current_time, dt, current_dof.data(), previous_dof.data(), 26, parameters.data(), 7,
                              sdvs.data(), 4, result.data(), 23, jacobian.data(), additional.data(), 0, output_message);

    BOOST_TEST(error_code == 0);

    if (error_code != 0) {
        std::cout << "ERROR MESSAGE:\n" << output_message << "\n";
    }

    BOOST_TEST(F_answer == modelJ.F, CHECK_PER_ELEMENT);

    BOOST_TEST(Fp_answer == modelJ.Fp, CHECK_PER_ELEMENT);

    BOOST_TEST(answer == result, CHECK_PER_ELEMENT);

    double eps = 1e-5;

    {
        constexpr unsigned int NUM_VAR = 26;
        constexpr unsigned int NUM_OUT = 23;

        std::vector<double> x = current_dof;

        for (unsigned int i = 0; i < NUM_VAR; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            std::vector<double> rp(23, 0);
            std::vector<double> rm(23, 0);

            modelMock mp, mm;

            int ep = mp.evaluate_model(current_time, dt, xp.data(), previous_dof.data(), 26, parameters.data(), 7,
                                       sdvs.data(), 4, rp.data(), 23, output_message);

            BOOST_TEST(ep == 0);

            int em = mm.evaluate_model(current_time, dt, xm.data(), previous_dof.data(), 26, parameters.data(), 7,
                                       sdvs.data(), 4, rm.data(), 23, output_message);

            BOOST_TEST(em == 0);

            for (unsigned int j = 0; j < NUM_OUT; ++j) {
                BOOST_TEST(jacobian[NUM_VAR * j + i] == (rp[j] - rm[j]) / (2 * delta));
            }
        }
    }
}
