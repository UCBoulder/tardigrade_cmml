/**
  * \file test_tardigrade_balance_equations_constraint_equations.cpp
  *
  * Tests for tardigrade_balance_equations_constraint_equations
  */

#include<tardigrade_constraint_equations.h>
#include<tardigrade_finite_element_utilities.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_hydraLinearTestMaterial.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define USE_EIGEN
#include<tardigrade_vector_tools.h>

#define BOOST_TEST_MODULE test_tardigrade_balance_equations_constraint_equations
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

BOOST_AUTO_TEST_CASE( test_computeInternalEnergyConstraint, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

}
