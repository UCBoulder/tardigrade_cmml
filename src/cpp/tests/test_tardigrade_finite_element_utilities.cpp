/**
  * \file test_tardigrade_finite_element_utilities.cpp
  *
  * Tests for tardigrade_finite_element_utilities
  */

#include<tardigrade_finite_element_utilities.h>
#include<sstream>
#include<fstream>
#include<iostream>
#include<array>

#define BOOST_TEST_MODULE test_tardigrade_finite_element_utilities
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

typedef tardigradeBalanceEquations::finiteElementUtilities::floatType floatType; //!< Define the float type to be the same as in the balance of mass

typedef tardigradeBalanceEquations::finiteElementUtilities::floatVector floatVector; //!< Define the float vector type to be the same as in the balance of mass

typedef tardigradeBalanceEquations::finiteElementUtilities::secondOrderTensor secondOrderTensor; //!< Define the second order tensor type to be the same as in the balance of mass

BOOST_AUTO_TEST_CASE( test_computeGradientSpatialJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the computation of the jacobian of the gradient of a quantity w.r.t. the spatial displacement
     */

    std::array< floatType, 12 > grad_a = { 0,  1,  2,
                                           3,  4,  5,
                                           6,  7,  8,
                                           9, 10, 11 };

    floatVector grad_test = { 0.1, 0.2, 0.3 };

    std::array< floatType, 12 > answer = { -0.2, -0.4, -0.6,
                                           -0.5, -1. , -1.5,
                                           -0.8, -1.6, -2.4,
                                           -1.1, -2.2, -3.3 };

    std::array< floatType, 12 > result;
    std::fill( std::begin( result ), std::end( result ), 0 );

    tardigradeBalanceEquations::finiteElementUtilities::computeGradientSpatialJacobian( std::begin( grad_a ), 12, grad_test, 2, std::begin( result ) );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

}

