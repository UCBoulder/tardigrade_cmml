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

    floatType dCdRho, dCdRhoDot;

    floatVector dCdGradRho, dCdV;

    secondOrderTensor dCdGradV;

    tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, density_gradient, velocity, velocity_gradient, result,
                                                                     dCdRho,  dCdRhoDot,   dCdGradRho      , dCdV,     dCdGradV );

    BOOST_TEST( answer == result );

    floatType eps = 1e-6;

    // Derivative w.r.t. density
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( density ) + eps;

        floatType xp = density + delta;

        floatType xm = density - delta;

        floatType vp;

        floatType vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( xp, density_dot, density_gradient, velocity, velocity_gradient, vp );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( xm, density_dot, density_gradient, velocity, velocity_gradient, vm );

        BOOST_TEST( dCdRho == ( vp - vm ) / ( 2 * delta ) );

    }

    // Derivative w.r.t. density dot
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( density_dot ) + eps;

        floatType xp = density_dot + delta;

        floatType xm = density_dot - delta;

        floatType vp;

        floatType vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, xp, density_gradient, velocity, velocity_gradient, vp );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, xm, density_gradient, velocity, velocity_gradient, vm );

        BOOST_TEST( dCdRhoDot == ( vp - vm ) / ( 2 * delta ) );

    }

    // Derivative w.r.t. density gradient
    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( density_gradient[ i ] ) + eps;

        floatVector xp = density_gradient;

        floatVector xm = density_gradient;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        floatType vp;

        floatType vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, xp, velocity, velocity_gradient, vp );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, xm, velocity, velocity_gradient, vm );

        BOOST_TEST( dCdGradRho[ i ] == ( vp - vm ) / ( 2 * delta ) );

    }

    // Derivative w.r.t. velocity
    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( velocity[ i ] ) + eps;

        floatVector xp = velocity;

        floatVector xm = velocity;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        floatType vp;

        floatType vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, density_gradient, xp, velocity_gradient, vp );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, density_gradient, xm, velocity_gradient, vm );

        BOOST_TEST( dCdV[ i ] == ( vp - vm ) / ( 2 * delta ) );

    }

    // Derivative w.r.t. velocity gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( velocity_gradient[ i ] ) + eps;

        secondOrderTensor xp = velocity_gradient;

        secondOrderTensor xm = velocity_gradient;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        floatType vp;

        floatType vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, density_gradient, velocity, xp, vp );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, density_gradient, velocity, xm, vm );

        BOOST_TEST( dCdGradV[ i ] == ( vp - vm ) / ( 2 * delta ) );

    }

}

BOOST_AUTO_TEST_CASE( test_multiphase_computeBalanceOfMass, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test computing the balance of mass
     */

    constexpr unsigned int nphases = 5;
    constexpr unsigned int dim     = 3;

    std::array<floatType,nphases> density = { 0.44775226, 0.10257244, 0.02078009, 0.34655081, 0.63241094 };

    std::array<floatType,nphases> density_dot = { 0.21140475, 0.84172361, 0.90589625, 0.65667824, 0.12870103 };

    std::array<floatType,nphases*dim> density_gradient = { 0.84973079, 0.69994479, 0.19507225, 0.71997378, 0.40400888,
                                                           0.10209729, 0.45538381, 0.78085324, 0.38934166, 0.80650313,
                                                           0.37377701, 0.52589524, 0.21972558, 0.6619266 , 0.71868273 };

    std::array<floatType,nphases*dim> velocity = { 1.78845045e-01, 3.43184239e-04, 3.03570624e-01, 2.61149246e-01,
                                                   6.22089112e-01, 4.77309376e-01, 2.83884563e-01, 8.54215703e-01,
                                                   6.62676921e-01, 1.26475940e-01, 7.91005742e-01, 6.58327605e-01,
                                                   9.64615247e-01, 5.82284188e-01, 8.29503047e-01 };

    std::array<floatType,nphases*dim*dim> velocity_gradient = { 0.08057272, 0.30451766, 0.56826345, 0.23438303, 0.61955306,
                                                                0.44571224, 0.15461077, 0.98556686, 0.01420748, 0.41255824,
                                                                0.1579279 , 0.98449085, 0.57013379, 0.15921173, 0.96501503,
                                                                0.15371508, 0.64557532, 0.08269366, 0.98004204, 0.59460498,
                                                                0.68023366, 0.88135227, 0.68412015, 0.51093195, 0.35998666,
                                                                0.62551418, 0.90176422, 0.15462472, 0.90648864, 0.51429834,
                                                                0.65717033, 0.29536049, 0.57881457, 0.58403957, 0.12585209,
                                                                0.99548878, 0.1307144 , 0.92037728, 0.54106702, 0.43851027,
                                                                0.73620888, 0.66502109, 0.64664504, 0.8294774 , 0.9394895 };

    std::array<floatType,nphases> answer = { 0.74267764, 1.39693567, 2.01351769, 1.90148278, 2.46462579 };

    std::array<floatType,nphases> result;

    tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                     std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                     std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( result ), std::end( result ) );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    std::array<floatType,nphases> dCdRho, dCdRhoDot;

    std::array<floatType,nphases*dim> dCdGradRho, dCdV;

    std::array<floatType,nphases*dim*dim> dCdGradV;

    tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                     std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                     std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( result ), std::end( result ),
                                                                     std::begin( dCdRho ),     std::end( dCdRho ),
                                                                     std::begin( dCdRhoDot ),  std::end( dCdRhoDot ),
                                                                     std::begin( dCdGradRho ), std::end( dCdGradRho ),
                                                                     std::begin( dCdV ),       std::end( dCdV ),
                                                                     std::begin( dCdGradV ),   std::end( dCdGradV ) );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    // Derivative w.r.t. density
    for ( unsigned int i = 0; i < nphases; i++ ){

        floatType delta = eps * std::fabs( density[ i ] ) + eps;

        std::array<floatType,nphases> xp = density;

        std::array<floatType,nphases> xm = density;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( xp ), std::end( xp ), std::begin( density_dot ), std::end( density_dot ),
                                                                         std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                         std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( xm ), std::end( xm ), std::begin( density_dot ), std::end( density_dot ),
                                                                         std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                         std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( j == i ){
                BOOST_TEST( dCdRho[ i ] == grad );
            }
            else{
                BOOST_TEST( grad == 0 );
            }

        }

    }

    // Derivative w.r.t. density dot
    for ( unsigned int i = 0; i < nphases; i++ ){

        floatType delta = eps * std::fabs( density_dot[ i ] ) + eps;

        std::array<floatType,nphases> xp = density_dot;

        std::array<floatType,nphases> xm = density_dot;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( density ), std::end( density ), std::begin( xp ), std::end( xp ),
                                                                         std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                         std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( density ), std::end( density ), std::begin( xm ), std::end( xm ),
                                                                         std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                         std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( j == i ){
                BOOST_TEST( dCdRhoDot[ i ] == grad );
            }
            else{
                BOOST_TEST( grad == 0 );
            }

        }

    }

    // Derivative w.r.t. density gradient
    for ( unsigned int i = 0; i < dim * nphases; i++ ){

        floatType delta = eps * std::fabs( density_gradient[ i ] ) + eps;

        std::array<floatType,nphases*dim> xp = density_gradient;

        std::array<floatType,nphases*dim> xm = density_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                         std::begin( xp ), std::end( xp ), std::begin( velocity ), std::end( velocity ),
                                                                         std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                         std::begin( xm ), std::end( xm ), std::begin( velocity ), std::end( velocity ),
                                                                         std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i / dim ) == j ){

                BOOST_TEST( dCdGradRho[ i ] == grad );

            }
            else{

                BOOST_TEST( grad == 0 );

            }

        }

    }

    // Derivative w.r.t. velocity
    for ( unsigned int i = 0; i < dim * nphases; i++ ){

        floatType delta = eps * std::fabs( velocity[ i ] ) + eps;

        std::array<floatType,nphases*dim> xp = velocity;

        std::array<floatType,nphases*dim> xm = velocity;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                         std::begin( density_gradient ), std::end( density_gradient ), std::begin( xp ), std::end( xp ),
                                                                         std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                         std::begin( density_gradient ), std::end( density_gradient ), std::begin( xm ), std::end( xm ),
                                                                         std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i / dim ) == j ){

                BOOST_TEST( dCdV[ i ] == grad );

            }
            else{

                BOOST_TEST( grad == 0 );

            }

        }

    }

    // Derivative w.r.t. velocity gradient
    for ( unsigned int i = 0; i < dim * dim * nphases; i++ ){

        floatType delta = eps * std::fabs( velocity_gradient[ i ] ) + eps;

        std::array<floatType,nphases*dim*dim> xp = velocity_gradient;

        std::array<floatType,nphases*dim*dim> xm = velocity_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                         std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                         std::begin( xp ), std::end( xp ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                         std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                         std::begin( xm ), std::end( xm ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i / ( dim * dim ) ) == j ){

                BOOST_TEST( dCdGradV[ i ] == grad );

            }
            else{

                BOOST_TEST( grad == 0 );

            }

        }

    }

}
