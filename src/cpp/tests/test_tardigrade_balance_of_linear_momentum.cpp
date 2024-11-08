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

    constexpr unsigned int dim = 3;

    floatType density = 0.69646919;

    floatType density_dot = 0.28613933;

    floatVector density_gradient = { 0.22685145, 0.55131477, 0.71946897 };

    floatVector velocity = { 0.42310646, 0.9807642,  0.68482974 };

    floatVector velocity_dot = { 0.4809319,  0.39211752, 0.34317802 };

    secondOrderTensor velocity_gradient = { 0.72904971, 0.43857224, 0.0596779,  0.39804426, 0.73799541, 0.18249173, 0.17545176, 0.53155137, 0.53182759 };

    floatVector body_force = { 0.63440096, 0.84943179, 0.72445532 };

    floatVector answer = { -1.62394626, -3.14362652, -2.32569958 };

    floatVector result;

    tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                      std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                      std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ), std::end( body_force ),
                                                                                                      std::begin( result ), std::end( result ) );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0. );

    floatVector dRdRho, dRdRhoDot;

    secondOrderTensor dRdGradRho, dRdV, dRdVDot, dRdB;

    std::array< floatType, dim * dim * dim > dRdGradV;

    tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                      std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                      std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ), std::end( body_force ),
                                                                                                      std::begin( result ),     std::end( result ),
                                                                                                      std::begin( dRdRho ),     std::end( dRdRho ),
                                                                                                      std::begin( dRdRhoDot ),  std::end( dRdRhoDot ),
                                                                                                      std::begin( dRdGradRho ), std::end( dRdGradRho ),
                                                                                                      std::begin( dRdV ),       std::end( dRdV ),
                                                                                                      std::begin( dRdVDot ),    std::end( dRdVDot ),
                                                                                                      std::begin( dRdGradV ),   std::end( dRdGradV ),
                                                                                                      std::begin( dRdB ),       std::end( dRdB ) );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( density ) + eps;

        floatType xp = density;
        floatType xm = density;

        xp += delta;
        xm -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( xp, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( xm, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdRho[ j ] == grad );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( density_dot ) + eps;

        floatType xp = density_dot;
        floatType xm = density_dot;

        xp += delta;
        xm -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, xp, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, xm, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdRhoDot[ j ] == grad );

        }

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( density_gradient[ i ] ) + eps;

        floatVector xp = density_gradient;
        floatVector xm = density_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( xp ), std::end( xp ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( xm ), std::end( xm ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdGradRho[ dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( velocity[ i ] ) + eps;

        floatVector xp = velocity;
        floatVector xm = velocity;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( xp ), std::end( xp ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( xm ), std::end( xm ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdV[ dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( velocity_dot[ i ] ) + eps;

        floatVector xp = velocity_dot;
        floatVector xm = velocity_dot;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( xp ), std::end( xp ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( xm ), std::end( xm ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdVDot[ dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( velocity_gradient[ i ] ) + eps;

        secondOrderTensor xp = velocity_gradient;
        secondOrderTensor xm = velocity_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( xp ), std::end( xp ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( xm ), std::end( xm ), std::begin( body_force ),
                                                                                                          std::end( body_force ),
                                                                                                          std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdGradV[ dim * dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( body_force[ i ] ) + eps;

        floatVector xp = body_force;
        floatVector xm = body_force;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( xp ),
                                                                                                          std::end( xp ),
                                                                                                          std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
                                                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( xm ),
                                                                                                          std::end( xm ),
                                                                                                          std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdB[ dim * j + i ] == grad );

        }

    }

}
