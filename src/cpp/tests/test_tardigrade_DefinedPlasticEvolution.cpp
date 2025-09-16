/**
  * \file test_tardigrade_DefinedPlasticEvolution.cpp
  *
  * Tests for tardigrade_cmml
  */

#include "tardigrade_DefinedPlasticEvolution.h"
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_tardigrade_DefinedPlasticEvolution
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer )
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

struct cerr_redirect{
    cerr_redirect( std::streambuf * new_buffer )
        : old( std::cerr.rdbuf( new_buffer ) )
    { }

    ~cerr_redirect( ) {
        std::cerr.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

BOOST_AUTO_TEST_CASE( test_basic_functionality, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the basic functionality of the DefinedPlasticEvolution class
     */

    class modelMock : public tardigradeCMML::DefinedPlasticEvolution::DefinedPlasticEvolution{

        public:

            modelMock( ) : tardigradeCMML::DefinedPlasticEvolution::DefinedPlasticEvolution( ){ }

            void public_extract_parameters( double * parameters_begin, const unsigned int parameters_size ){

                extract_parameters( parameters_begin, parameters_size );

            }

    };

    modelMock model;

    std::vector< double > parameters = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    try{
    model.public_extract_parameters( parameters.data( ), parameters.size( ) );
    }catch(std::exception &e){tardigradeErrorTools::printNestedExceptions(e); throw;}

    BOOST_TEST( model.getEvaluateModelResultSize( ) == 26 );

    BOOST_TEST( model.getIsCurrent( ) );

    BOOST_TEST( model.getDensityIndex( ) == 1 );

    BOOST_TEST( model.getDOFInternalEnergyIndex( ) == 2 );

    BOOST_TEST( model.getDefinedVelocityGradientIndex( ) == 3 );

    BOOST_TEST( model.getInternalEnergyScaledByDensity( ) );

    BOOST_TEST( model.getDensityGradientIndex( ) == 5 );

    BOOST_TEST( model.getDisplacementGradientIndex( ) == 6 );

    BOOST_TEST( model.getTemperatureIndex( ) == 7 );

}

BOOST_AUTO_TEST_CASE( test_basic_functionality_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the basic functionality of the DefinedPlasticEvolution class
     */

    class modelMock : public tardigradeCMML::DefinedPlasticEvolution::DefinedPlasticEvolution{

        public:

            modelMock( ) : tardigradeCMML::DefinedPlasticEvolution::DefinedPlasticEvolution( ){ }

            void public_extract_parameters( double * parameters_begin, const unsigned int parameters_size ){

                extract_parameters( parameters_begin, parameters_size );

            }

    };

    modelMock model;

    std::vector< double > parameters = { 1, 2, 3, 0.2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    model.public_extract_parameters( parameters.data( ), parameters.size( ) );

    BOOST_TEST( model.getEvaluateModelResultSize( ) == 26 );

    BOOST_TEST( model.getIsCurrent( ) );

    BOOST_TEST( model.getDensityIndex( ) == 1 );

    BOOST_TEST( model.getDOFInternalEnergyIndex( ) == 2 );

    BOOST_TEST( model.getDefinedVelocityGradientIndex( ) == 3 );

    BOOST_TEST( !model.getInternalEnergyScaledByDensity( ) );

    BOOST_TEST( model.getDensityGradientIndex( ) == 5 );

    BOOST_TEST( model.getDisplacementGradientIndex( ) == 6 );

    BOOST_TEST( model.getTemperatureIndex( ) == 7 );

}

BOOST_AUTO_TEST_CASE( test_evaluate_model, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the evaluate_model function of the DefinedPlasticEvolution class
     */

    class modelMock : public tardigradeCMML::DefinedPlasticEvolution::DefinedPlasticEvolution{

        public:

            modelMock( ) : tardigradeCMML::DefinedPlasticEvolution::DefinedPlasticEvolution( ){ }

            std::vector< double > F, Fp;

            void public_extract_parameters( double * parameters_begin, const unsigned int parameters_size ){

                extract_parameters( parameters_begin, parameters_size );

            }

        protected:

            virtual void formDeformationGradients(
                const double *current_dof_begin, const double *previous_dof_begin,
                std::vector< double > &deformationGradient, std::vector< double > &previousDeformationGradient
            ) override{

                tardigradeCMML::BasicSolid::BasicSolid::formDeformationGradients(
                    current_dof_begin,   previous_dof_begin,
                    deformationGradient, previousDeformationGradient
                );

                F  = deformationGradient;
                Fp = previousDeformationGradient;

            }

            virtual void formDeformationGradients(
                const double *current_dof_begin, const double *previous_dof_begin,
                std::vector< double > &deformationGradient, std::vector< double > &previousDeformationGradient,
                std::vector< double > &dFdGradU
            ) override{

                tardigradeCMML::BasicSolid::BasicSolid::formDeformationGradients(
                    current_dof_begin,   previous_dof_begin,
                    deformationGradient, previousDeformationGradient,
                    dFdGradU
                );

                F  = deformationGradient;
                Fp = previousDeformationGradient;

            }

    };

    modelMock model;

    double current_time = 1.23;
    double dt = 0.345;

    std::vector< double >  current_dof( 26, 0 );
    std::vector< double > previous_dof( 26, 0 );
    for ( unsigned int i = 0; i < 26; ++i ){
        current_dof[ i ]  += 0.01 * ( i + 1 );
        previous_dof[ i ] -= 0.01 * ( i + 1 );
    }

    std::vector< double > parameters = { 1, 25, 16, 1, 14, 2, 3, 4, 5, 6, 7, 8, 0.78, 0.89, 9 };

    std::vector< double > sdvs_base( 18, 0 );
    std::vector< double > sdvs = sdvs_base;
    std::vector< double > sdvsJ = sdvs_base;

    std::vector< double > result( 26, -1 );

    std::string output_message;

    std::vector< double > F_answer = {
        1.03996448, 0.05150977, 0.06305506,
        0.07688404, 1.08957117, 0.10225831,
        0.1138036 , 0.12763258, 1.14146156
    };

    std::vector< double > Fp_answer = {
        0.97641119, -0.03261049, -0.04163218,
       -0.04916405,  0.94255918, -0.0657176 ,
       -0.07473928, -0.08227115,  0.91019699
    };

    std::vector< double > answer = {
        2.0100243 ,  1.2556396 ,  1.73524007,  1.2556396 ,  3.0772629 ,
        2.4309158 ,  1.73524007,  2.4309158 ,  4.360577  ,  0.28      ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        , -0.4       , -0.48      , -0.56      ,
        0.        ,  0.        ,  0.        , -1.35      , -1.44      ,
       -1.53
    };

    std::vector< double > sdvs_answer = {
        3.76430005e-03,  5.15758065e-03,  6.55086124e-03, -1.26713080e-03,
        4.24258814e-05,  1.35198255e-03, -6.29856164e-03, -5.07272888e-03,
       -3.84689614e-03,  2.80000000e-01, -4.00000000e-01, -4.80000000e-01,
       -5.60000000e-01, -1.35000000e+00, -1.44000000e+00, -1.53000000e+00,
        0, 0
    };

    std::cerr << "entering evaluate_model\n";
    int error_code = model.evaluate_model(
        current_time, dt,
        current_dof.data( ),
        previous_dof.data( ), 26,
        parameters.data( ),   15,
        sdvs.data( ),         18,
        result.data( ),       26,
        output_message
    );
    std::cerr << "exiting evaluate_model\n";

    BOOST_TEST( error_code == 0 );

    if ( error_code != 0 ){
        std::cout << "ERROR MESSAGE:\n" << output_message << "\n";
    }

    BOOST_TEST(   F_answer == model.F, CHECK_PER_ELEMENT );

    BOOST_TEST( Fp_answer == model.Fp, CHECK_PER_ELEMENT );

    BOOST_TEST(   sdvs == sdvs_answer, CHECK_PER_ELEMENT );

    BOOST_TEST(      answer == result, CHECK_PER_ELEMENT );

    // Test the jacobians
    result = std::vector< double >( 26, -1 );
    std::vector< double > jacobian( 26 * 26, 0 );
    std::vector< double > additional( 0, 0 );

    modelMock modelJ;
    error_code = modelJ.evaluate_model(
        current_time, dt,
        current_dof.data( ),
        previous_dof.data( ), 26,
        parameters.data( ),   15,
        sdvsJ.data( ),        18,
        result.data( ),       26,
        jacobian.data( ),
        additional.data( ),    0,
        output_message
    );

    BOOST_TEST( error_code == 0 );

    if ( error_code != 0 ){
        std::cout << "ERROR MESSAGE:\n" << output_message << "\n";
    }

    BOOST_TEST(   F_answer == modelJ.F, CHECK_PER_ELEMENT );

    BOOST_TEST( Fp_answer == modelJ.Fp, CHECK_PER_ELEMENT );

    BOOST_TEST(   sdvsJ == sdvs_answer, CHECK_PER_ELEMENT );

    BOOST_TEST(       answer == result, CHECK_PER_ELEMENT );

    double eps = 1e-5;

    {

        constexpr unsigned int NUM_VAR = 26;
        constexpr unsigned int NUM_OUT = 26;

        std::vector< double > x = current_dof;

        for ( unsigned int i = 0; i < NUM_VAR; ++i ){

            double delta = eps * std::fabs( x[ i ] ) + eps;

            std::vector< double > xp = x;
            std::vector< double > xm = x;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::vector< double > rp( NUM_OUT, 0 );
            std::vector< double > rm( NUM_OUT, 0 );

            modelMock mp, mm;

            sdvs = sdvs_base;
            int ep = mp.evaluate_model(
                current_time, dt,
                xp.data( ),
                previous_dof.data( ), 26,
                parameters.data( ),   15,
                sdvs.data( ),         18,
                rp.data( ),           26,
                output_message
            );

            BOOST_TEST( ep == 0 );

            sdvs = sdvs_base;
            int em = mm.evaluate_model(
                current_time, dt,
                xm.data( ),
                previous_dof.data( ), 26,
                parameters.data( ),   15,
                sdvs.data( ),         18,
                rm.data( ),           26,
                output_message
            );

            BOOST_TEST( em == 0 );

            for ( unsigned int j = 0; j < NUM_OUT; ++j ){

                BOOST_TEST( jacobian[ NUM_VAR * j + i ] == ( rp[ j ] - rm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}
