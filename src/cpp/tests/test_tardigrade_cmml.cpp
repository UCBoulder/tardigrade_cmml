/**
  * \file test_tardigrade_cmml.cpp
  *
  * Tests for tardigrade_cmml
  */

#include<tardigrade_cmml.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_tardigrade_cmml
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

#ifdef TARDIGRADE_CMML_BUILD_AS_SHARED
BOOST_AUTO_TEST_CASE( existent_material, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test if a good material can be extracted
     */

    auto &factory = tardigradeCMML::MaterialFactory::Instance( );
    auto material = factory.GetMaterial("tardigradeCMML::BasicSolid");

    double current_time = 1.23;
    double dt = 0.345;

    std::vector< double >  current_dof( 26, 0 );
    std::vector< double > previous_dof( 26, 0 );
    for ( unsigned int i = 0; i < 26; ++i ){
        current_dof[ i ]  += 0.01 * ( i + 1 );
        previous_dof[ i ] -= 0.01 * ( i + 1 );
    }

    std::vector< double > parameters = { 2, 3, 4, 5, 6, 7, 8 };

    std::vector< double > sdvs( 4, 0 );

    std::vector< double > result( 23, -1 );

    std::string output_message;

    std::vector< double > answer = {
        2.0707740, 1.2922023, 1.7429554, 1.2922023, 3.0650159,
        2.37918  , 1.7429554, 2.37918  , 4.24473  , 0.28     ,
        0        , 0        ,  0       ,  0       ,  0       ,
        0        , 0        , -0.4     , -0.48    , -0.56    ,
        0        , 0        ,  0
    };

    int error_code = material->evaluate_model(
        current_time, dt,
        current_dof.data( ),
        previous_dof.data( ), 26,
        parameters.data( ),    7,
        sdvs.data( ),          4,
        result.data( ),       23,
        output_message
    );

    BOOST_TEST( error_code == 0 );

    if ( error_code != 0 ){
        std::cout << "ERROR MESSAGE:\n" << output_message << "\n";
    }

    BOOST_TEST(    answer == result,   CHECK_PER_ELEMENT );
    

}
#endif

BOOST_AUTO_TEST_CASE( non_existent_material_error, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test if an error is thrown for a bad material
     */

    auto &factory = tardigradeCMML::MaterialFactory::Instance( );
    BOOST_CHECK_THROW( factory.GetMaterial("not_a_material"), std::runtime_error );

}

