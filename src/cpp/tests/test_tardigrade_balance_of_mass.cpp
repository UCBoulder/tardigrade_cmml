/**
  * \file test_tardigrade_balance_equations_balance_of_mass.cpp
  *
  * Tests for tardigrade_balance_equations_balance_of_mass
  */

#include<tardigrade_balance_of_mass.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_tardigrade_balance_equations_balance_of_mass
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

typedef tardigradeBalanceEquations::balanceOfMass::floatType floatType; //!< Define the float type to be the same as in the balance of mass

typedef tardigradeBalanceEquations::balanceOfMass::floatVector floatVector; //!< Define the float vector type to be the same as in the balance of mass

typedef tardigradeBalanceEquations::balanceOfMass::secondOrderTensor secondOrderTensor; //!< Define the second order tensor type to be the same as in the balance of mass

BOOST_AUTO_TEST_CASE( test_computeBalanceOfMass, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test computing the balance of mass
     */

    floatType density = 1.4;

    floatType density_dot = 2.8;

    floatVector density_gradient = { 1., 2., 3. };

    floatVector velocity = { 0.1, 0.2, 0.3 };

    secondOrderTensor velocity_gradient = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

    floatType answer = 6.3;

    floatType result;
    
    tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, density_gradient, velocity, velocity_gradient, result );

    BOOST_TEST( answer == result );

}

