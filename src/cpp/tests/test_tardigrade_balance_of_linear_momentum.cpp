/**
  * \file test_tardigrade_balance_equations_balance_of_linear_momentum.cpp
  *
  * Tests for tardigrade_balance_equations_balance_of_linear_momentum
  */

#include<tardigrade_balance_of_linear_momentum.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_tardigrade_balance_equations_balance_of_linear_momentum
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

typedef tardigradeBalanceEquations::balanceOfLinearMomentum::floatType floatType; //!< Define the float type to be the same as in the balance of linear momentum

typedef tardigradeBalanceEquations::balanceOfLinearMomentum::floatVector floatVector; //!< Define the float vector type to be the same as in the balance of linear momentum

typedef tardigradeBalanceEquations::balanceOfLinearMomentum::secondOrderTensor secondOrderTensor; //!< Define the second order tensor type to be the same as in the balance of linear momentum

BOOST_AUTO_TEST_CASE( test_computeBalanceOfLinearMomentum, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

}
