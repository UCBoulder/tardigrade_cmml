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

BOOST_AUTO_TEST_CASE( test_LinearHex, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    std::array< floatType, 24 > xi = { -1, -1, -1,
                                        1, -1, -1,
                                        1,  1, -1,
                                       -1,  1, -1,
                                       -1, -1,  1,
                                        1, -1,  1,
                                        1,  1,  1,
                                       -1,  1,  1 };

    std::array< floatType, 24 > X = { 0, 0, 0,
                                      1, 0, 0,
                                      1, 1, 0,
                                      0, 1, 0,
                                      0, 0, 1,
                                      1, 0, 1,
                                      1, 1, 1,
                                      0, 1, 1 };

    std::array< floatType, 9 > A = { 1.01964692, -0.02138607, -0.02731485,
                                     0.00513148,  1.0219469 , -0.00768935,
                                     0.04807642,  0.01848297,  0.99809319 };

    std::array< floatType, 3 > b = { -0.21576496, -0.31364397,  0.45809941 };

    std::array< floatType, 24 > x;

    std::fill( std::begin( x ), std::end( x ), 0 );

    for ( unsigned int i = 0; i < 8; ++i ){
        for ( unsigned int j = 0; j < 3; ++j ){
            for ( unsigned int k = 0; k < 3; ++k ){
                x[ 3 * i + j ] += X[ 3 * i + k ] * A[ 3 * j + k ];
            }
            x[ 3 * i + j ] += b[ j ];
        }
    }

    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<floatType, typename std::array< floatType, 24 >::const_iterator, typename std::array< floatType, 3>::const_iterator, typename std::array< floatType, 8>::iterator, typename std::array< floatType, 24>::iterator> e( std::cbegin( x ), std::cend( x ), std::cbegin( X ), std::cend( X ) );

    std::array< floatType, 8 > result;

    std::array< floatType, 8 > answer;

    for ( unsigned int i = 0; i < 8; i++ ){

        e.GetShapeFunctions( std::cbegin( xi ) + 3 * i, std::cbegin( xi ) + 3 * ( i + 1 ), std::begin( result ), std::end( result ) );

        std::fill( std::begin( answer ), std::end( answer ), 0 );

        answer[ i ] = 1;

        BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    }

    std::array< floatType, 3 > point = { 0, 0, 0 };

    e.GetShapeFunctions( std::cbegin( point ), std::cend( point ), std::begin( result ), std::end( result ) );

    std::fill( std::begin( answer ), std::end( answer ), 0.125 );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    std::array< floatType, 8 * 3 > dNdxi;

    e.GetLocalShapeFunctionGradients( std::cbegin( point ), std::cend( point ), std::begin( dNdxi ), std::end( dNdxi ) );

    floatType eps = 1e-6;

    std::array< floatType, 8 > vp, vm;

    for ( unsigned int i = 0; i < 3; ++i ){

        floatType delta = eps * std::fabs( point[ i ] ) + eps;

        std::array< floatType, 3 > pp = point;
        std::array< floatType, 3 > pm = point;

        pp[ i ] += delta;
        pm[ i ] -= delta;

        e.GetShapeFunctions( std::cbegin( pp ), std::cend( pp ), std::begin( vp ), std::end( vp ) );
        e.GetShapeFunctions( std::cbegin( pm ), std::cend( pm ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < 8; j++ ){

            BOOST_TEST( ( ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) ) == dNdxi[ 3 * j + i ] );

        }

    }

}

BOOST_AUTO_TEST_CASE( test_LinearHex2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE) ) {

    // Test the interpolation
    std::array< floatType, 40 > quantities = { 0.03096734, 0.25428207, 0.91240044, 0.9701742 , 0.24661978,
                                               0.69880741, 0.86642932, 0.63952185, 0.45684365, 0.49879475,
                                               0.72047501, 0.27125319, 0.80061903, 0.50925268, 0.81992102,
                                               0.30771324, 0.65927822, 0.31350889, 0.54465191, 0.31984146,
                                               0.45915279, 0.81309162, 0.95755166, 0.0039935 , 0.98393864,
                                               0.98961971, 0.09200859, 0.46202094, 0.87754638, 0.37499295,
                                               0.23868804, 0.37254533, 0.54245204, 0.73809604, 0.87751833,
                                               0.53958266, 0.07906357, 0.74954525, 0.70097189, 0.83657464 };

    std::array< floatType, 24 > X = { 0, 0, 0,
                                      1, 0, 0,
                                      1, 1, 0,
                                      0, 1, 0,
                                      0, 0, 1,
                                      1, 0, 1,
                                      1, 1, 1,
                                      0, 1, 1 };

    std::array< floatType, 9 > A = { 1.01964692, -0.02138607, -0.02731485,
                                     0.00513148,  1.0219469 , -0.00768935,
                                     0.04807642,  0.01848297,  0.99809319 };

    std::array< floatType, 3 > b = { -0.21576496, -0.31364397,  0.45809941 };

    std::array< floatType, 24 > x;

    std::fill( std::begin( x ), std::end( x ), 0 );

    for ( unsigned int i = 0; i < 8; ++i ){
        for ( unsigned int j = 0; j < 3; ++j ){
            for ( unsigned int k = 0; k < 3; ++k ){
                x[ 3 * i + j ] += X[ 3 * i + k ] * A[ 3 * j + k ];
            }
            x[ 3 * i + j ] += b[ j ];
        }
    }

    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<floatType, typename std::array< floatType, 24 >::const_iterator, typename std::array< floatType, 3>::const_iterator, typename std::array< floatType, 8>::iterator, typename std::array< floatType, 24>::iterator> e( std::cbegin( x ), std::cend( x ), std::cbegin( X ), std::cend( X ) );

    std::array< floatType, 3 > point = { 0.1626388 , 0.45020513, 0.22368613 };

    std::array< floatType, 5 > answer = { 0.48974732, 0.3641928 , 0.64009587, 0.62582472, 0.70145625 };

    std::array< floatType, 5 > result;

    e.InterpolateQuantity( std::cbegin( point ), std::cend( point ), std::cbegin( quantities ), std::cend( quantities ), std::begin( result ), std::end( result ) );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    std::array< floatType, 15 > dQuantitydxi;

    e.GetLocalQuantityGradient( std::cbegin( point ), std::cend( point ), std::cbegin( quantities ), std::cend( quantities ), std::begin( dQuantitydxi ), std::end( dQuantitydxi ) );

    floatType eps = 1e-6;

    std::array< floatType, 5 > vp, vm;

    for ( unsigned int i = 0; i < 3; ++i ){

        floatType delta = eps * std::fabs( point[ i ] ) + eps;

        std::array< floatType, 3 > pp = point;
        std::array< floatType, 3 > pm = point;

        pp[ i ] += delta;
        pm[ i ] -= delta;

        e.InterpolateQuantity( std::cbegin( pp ), std::cend( pp ), std::begin( quantities ), std::end( quantities ), std::begin( vp ), std::end( vp ) );
        e.InterpolateQuantity( std::cbegin( pm ), std::cend( pm ), std::begin( quantities ), std::end( quantities ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < 5; j++ ){

            BOOST_TEST( ( ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) ) == dQuantitydxi[ 3 * j + i ] );

        }

    }

}

BOOST_AUTO_TEST_CASE( test_LinearHex3, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE) ) {

    // Test the global shape function gradients
    std::array< floatType, 24 > X = { 0, 0, 0,
                                      1, 0, 0,
                                      1, 1, 0,
                                      0, 1, 0,
                                      0, 0, 1,
                                      1, 0, 1,
                                      1, 1, 1,
                                      0, 1, 1 };

    std::array< floatType, 9 > A = { 1.01964692, -0.02138607, -0.02731485,
                                     0.00513148,  1.0219469 , -0.00768935,
                                     0.04807642,  0.01848297,  0.99809319 };

    std::array< floatType, 3 > b = { -0.21576496, -0.31364397,  0.45809941 };

    std::array< floatType, 24 > x;

    std::fill( std::begin( x ), std::end( x ), 0 );

    for ( unsigned int i = 0; i < 8; ++i ){
        for ( unsigned int j = 0; j < 3; ++j ){
            for ( unsigned int k = 0; k < 3; ++k ){
                x[ 3 * i + j ] += X[ 3 * i + k ] * A[ 3 * j + k ];
            }
            x[ 3 * i + j ] += b[ j ];
        }
    }

    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<floatType, typename std::array< floatType, 24 >::const_iterator, typename std::array< floatType, 3>::const_iterator, typename std::array< floatType, 8>::iterator, typename std::array< floatType, 24>::iterator> e( std::cbegin( x ), std::cend( x ), std::cbegin( X ), std::cend( X ) );

    std::array< floatType, 3 > point = { 0.1626388 , 0.45020513, 0.22368613 };

    std::array< floatType, 24 > answer = { -0.09822601, -0.15892292, -0.11922661,  0.11321351, -0.21555842,
                                           -0.15867086,  0.2942983 ,  0.23441585, -0.41246145, -0.26220936,
                                            0.15914528, -0.31011619, -0.16879031, -0.25616407,  0.10872132,
                                            0.15907466, -0.347635  ,  0.16178372,  0.41277349,  0.34878489,
                                            0.43630492, -0.45013428,  0.23593439,  0.29366517 };

    std::array< floatType, 24 > result;

    e.GetGlobalShapeFunctionGradients( std::cbegin( point ), std::cend( point ), std::cbegin( x ), std::cend( x ), std::begin( result ), std::end( result ) );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    answer = { -0.10670335, -0.16251378, -0.11509422,  0.10670335, -0.22564316,
               -0.15980321,  0.28145359,  0.22564316, -0.42151619, -0.28145359,
                0.16251378, -0.30358638, -0.16819409, -0.25616682,  0.11509422,
                0.16819409, -0.35567624,  0.15980321,  0.44364897,  0.35567624,
                0.42151619, -0.44364897,  0.25616682,  0.30358638 };

    e.GetGlobalShapeFunctionGradients( std::cbegin( point ), std::cend( point ), std::cbegin( X ), std::cend( X ), std::begin( result ), std::end( result ) );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_LinearHex4, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE) ) {

    // Test the global gradient of a quantity
    std::array< floatType, 24 > X = { 0, 0, 0,
                                      1, 0, 0,
                                      1, 1, 0,
                                      0, 1, 0,
                                      0, 0, 1,
                                      1, 0, 1,
                                      1, 1, 1,
                                      0, 1, 1 };

    std::array< floatType, 9 > A = { 1.01964692, -0.02138607, -0.02731485,
                                     0.00513148,  1.0219469 , -0.00768935,
                                     0.04807642,  0.01848297,  0.99809319 };

    std::array< floatType, 3 > b = { -0.21576496, -0.31364397,  0.45809941 };

    std::array< floatType, 24 > x;

    std::fill( std::begin( x ), std::end( x ), 0 );

    for ( unsigned int i = 0; i < 8; ++i ){
        for ( unsigned int j = 0; j < 3; ++j ){
            for ( unsigned int k = 0; k < 3; ++k ){
                x[ 3 * i + j ] += X[ 3 * i + k ] * A[ 3 * j + k ];
            }
            x[ 3 * i + j ] += b[ j ];
        }
    }

    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<floatType, typename std::array< floatType, 24 >::const_iterator, typename std::array< floatType, 3>::const_iterator, typename std::array< floatType, 8>::iterator, typename std::array< floatType, 24>::iterator> e( std::cbegin( x ), std::cend( x ), std::cbegin( X ), std::cend( X ) );

    std::array< floatType, 3 > point = { 0.1626388 , 0.45020513, 0.22368613 };

    std::array< floatType, 15 > answer = { 1,  2,  3,
                                           4,  5,  6,
                                           7,  8,  9,
                                          -6, -5, -4,
                                          -1, -2, -3 };

    std::array< floatType, 5 > b2 = { -0.1, -0.2, -0.3, -0.4, -0.5 };

    std::array< floatType, 40 > quantity;
    std::fill( std::begin( quantity ), std::end( quantity ), 0 );

    for ( unsigned int i = 0; i < 8; ++i ){
        for ( unsigned int j = 0; j < 5; ++j ){
            for ( unsigned int k = 0; k < 3; ++k ){
                quantity[ 5 * i + j ] += answer[ 3 * j + k ] * x[ 3 * i + k ];
            }
            quantity[ 5 * i + j ] += b2[ j ];
        }
    }

    std::array< floatType, 15 > result;

    e.GetGlobalQuantityGradient( std::cbegin( point ), std::cend( point ), std::cbegin( quantity ), std::cend( quantity ), std::begin( result ), std::end( result ) );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::fill( std::begin( quantity ), std::end( quantity ), 0 );

    for ( unsigned int i = 0; i < 8; ++i ){
        for ( unsigned int j = 0; j < 5; ++j ){
            for ( unsigned int k = 0; k < 3; ++k ){
                quantity[ 5 * i + j ] += answer[ 3 * j + k ] * X[ 3 * i + k ];
            }
            quantity[ 5 * i + j ] += b2[ j ];
        }
    }

    e.GetGlobalQuantityGradient( std::cbegin( point ), std::cend( point ), std::cbegin( quantity ), std::cend( quantity ), std::begin( result ), std::end( result ), 0 );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    std::array< floatType, 24 > dNdx;
    e.GetGlobalShapeFunctionGradients( std::cbegin( point ), std::cend( point ), std::cbegin( x ), std::cend( x ), std::begin( dNdx ), std::end( dNdx ) );

    for ( unsigned int i = 0; i < 40; ++i ){

        floatType delta = eps * std::fabs( quantity[ i ] ) + eps;

        std::array< floatType, 40 > qp, qm;

        qp = quantity;
        qm = quantity;

        qp[ i ] += delta;
        qm[ i ] -= delta;

        std::array< floatType, 15 > vp, vm;

        e.GetGlobalQuantityGradient(
            std::cbegin( point ), std::cend( point ),
            std::cbegin( qp ), std::cend( qp ),
            std::begin( vp ), std::end( vp )
        );

        e.GetGlobalQuantityGradient(
            std::cbegin( point ), std::cend( point ),
            std::cbegin( qm ), std::cend( qm ),
            std::begin( vm ), std::end( vm )
        );

        unsigned int node = i / 5;
        unsigned int component = i - node * 5;

        for ( unsigned int j = 0; j < 15; ++j ){

            unsigned int row = j / 3;
            unsigned int col = j - 3 * row;

            if ( component == row ){

                BOOST_TEST( dNdx[ 3 * node + col ] == ( ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) ) );

            }
            else{

                BOOST_TEST( ( ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) ) == 0 );

            }
            
        }

    }

}

BOOST_AUTO_TEST_CASE( test_LinearHex5, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE) ) {

    // Test the calculation of the volume integral Jacobian of transformation
    std::array< floatType, 24 > X = { 0, 0, 0,
                                      1, 0, 0,
                                      1, 1, 0,
                                      0, 1, 0,
                                      0, 0, 1,
                                      1, 0, 1,
                                      1, 1, 1,
                                      0, 1, 1 };

    std::array< floatType, 9 > A = { 1.01964692, -0.02138607, -0.02731485,
                                     0.00513148,  1.0219469 , -0.00768935,
                                     0.04807642,  0.01848297,  0.99809319 };

    std::array< floatType, 3 > b = { -0.21576496, -0.31364397,  0.45809941 };

    std::array< floatType, 24 > x;

    std::fill( std::begin( x ), std::end( x ), 0 );

    for ( unsigned int i = 0; i < 8; ++i ){
        for ( unsigned int j = 0; j < 3; ++j ){
            for ( unsigned int k = 0; k < 3; ++k ){
                x[ 3 * i + j ] += X[ 3 * i + k ] * A[ 3 * j + k ];
            }
            x[ 3 * i + j ] += b[ j ];
        }
    }

    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<floatType, typename std::array< floatType, 24 >::const_iterator, typename std::array< floatType, 3>::const_iterator, typename std::array< floatType, 8>::iterator, typename std::array< floatType, 24>::iterator> e( std::cbegin( x ), std::cend( x ), std::cbegin( X ), std::cend( X ) );

    std::array< floatType, 3 > point = { 0.1626388 , 0.45020513, 0.22368613 };

    floatType answer = 0.125;

    floatType result;

    e.GetVolumeIntegralJacobianOfTransformation( std::cbegin( point ), std::cend( point ), result, 0 );

    BOOST_TEST( result == answer );

    answer = 0.125 * 1.0416398444780686;

    e.GetVolumeIntegralJacobianOfTransformation( std::cbegin( point ), std::cend( point ), result );

    BOOST_TEST( result == answer );

}
