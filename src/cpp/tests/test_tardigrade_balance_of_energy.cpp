/**
  * \file test_tardigrade_balance_equations_balance_of_energy.cpp
  *
  * Tests for tardigrade_balance_equations_balance_of_energy
  */

#include<tardigrade_balance_of_energy.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_tardigrade_balance_equations_balance_of_energy
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

typedef tardigradeBalanceEquations::balanceOfEnergy::floatType floatType; //!< Define the float type to be the same as in the balance of linear momentum

typedef tardigradeBalanceEquations::balanceOfEnergy::floatVector floatVector; //!< Define the float vector type to be the same as in the balance of linear momentum

typedef tardigradeBalanceEquations::balanceOfEnergy::secondOrderTensor secondOrderTensor; //!< Define the second order tensor type to be the same as in the balance of linear momentum

BOOST_AUTO_TEST_CASE( test_computeBalanceOfEnergyNonDivergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    constexpr unsigned int dim = 3;

    floatType density = 0.69646919;

    floatType density_dot = 0.28613933;

    floatVector density_gradient = { 0.22685145, 0.55131477, 0.71946897 };

    floatType internal_energy = 0.42310646;

    floatType internal_energy_dot = 0.9807642;

    floatVector internal_energy_gradient = { 0.68482974, 0.4809319 , 0.39211752 };

    floatVector velocity = { 0.34317802, 0.72904971, 0.43857224 };

    secondOrderTensor velocity_gradient = { 0.0596779 , 0.39804426, 0.73799541, 0.18249173, 0.17545176,
                                            0.53155137, 0.53182759, 0.63440096, 0.84943179 };

    secondOrderTensor cauchy_stress = { 0.72445532, 0.61102351, 0.72244338, 0.32295891, 0.36178866,
                                        0.22826323, 0.29371405, 0.63097612, 0.09210494 };

    floatType volume_fraction = 0.43370117;

    floatType internal_heat_generation = 0.43086276;

    floatVector net_interphase_force = { 0.4936851 , 0.42583029, 0.31226122 };

    floatType answer = 0.87843988;

    floatType result;

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                      internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                      std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                      std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                      std::begin( net_interphase_force ), std::end( net_interphase_force ), result );

    BOOST_TEST( result == answer );

    result = 0;

    floatType dRdRho, dRdRhoDot, dRdE, dRdEDot, dRdPhi, dRdr;

    floatVector dRdGradRho, dRdGradE, dRdV, dRdpi;

    secondOrderTensor dRdGradV, dRdCauchy;

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                      internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                      std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                      std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                      std::begin( net_interphase_force ), std::end( net_interphase_force ), result,
                                                                                      dRdRho, dRdRhoDot, std::begin( dRdGradRho ), std::end( dRdGradRho ),
                                                                                      dRdE, dRdEDot, std::begin( dRdGradE ), std::end( dRdGradE ),
                                                                                      std::begin( dRdV ), std::end( dRdV ), std::begin( dRdGradV ), std::end( dRdGradV ),
                                                                                      std::begin( dRdCauchy ), std::end( dRdCauchy ), dRdPhi, dRdr,
                                                                                      std::begin( dRdpi ), std::end( dRdpi ) );

    BOOST_TEST( result == answer );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( density ) + eps;

        floatType xp = density;
        floatType xm = density;

        xp += delta;
        xm -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( xp, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( xm, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdRho == grad );

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( density_dot ) + eps;

        floatType xp = density_dot;
        floatType xm = density_dot;

        xp += delta;
        xm -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, xp, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, xm, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdRhoDot == grad );

    }

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::abs( density_gradient[ i ] ) + eps;

        floatVector xp = density_gradient;
        floatVector xm = density_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( xp ), std::end( xp ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( xm ), std::end( xm ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdGradRho[ i ] == grad );

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( internal_energy ) + eps;

        floatType xp = internal_energy;
        floatType xm = internal_energy;

        xp += delta;
        xm -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          xp, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          xm, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdE == grad );

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( internal_energy_dot ) + eps;

        floatType xp = internal_energy_dot;
        floatType xm = internal_energy_dot;

        xp += delta;
        xm -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, xp, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, xm, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdEDot == grad );

    }

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::abs( density_gradient[ i ] ) + eps;

        floatVector xp = density_gradient;
        floatVector xm = density_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( xp ), std::end( xp ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( xm ), std::end( xm ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdGradE[ i ] == grad );

    }

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::abs( velocity[ i ] ) + eps;

        floatVector xp = velocity;
        floatVector xm = velocity;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( xp ), std::end( xp ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( xm ), std::end( xm ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdV[ i ] == grad );

    }

    for ( unsigned int i = 0; i < dim * dim; i++ ){

        floatType delta = eps * std::abs( velocity_gradient[ i ] ) + eps;

        secondOrderTensor xp = velocity_gradient;
        secondOrderTensor xm = velocity_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( xp ), std::end( xp ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( xm ), std::end( xm ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdGradV[ i ] == grad );

    }

    for ( unsigned int i = 0; i < dim * dim; i++ ){

        floatType delta = eps * std::abs( cauchy_stress[ i ] ) + eps;

        secondOrderTensor xp = cauchy_stress;
        secondOrderTensor xm = cauchy_stress;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( xp ), std::end( xp ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( xm ), std::end( xm ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdCauchy[ i ] == grad );

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( volume_fraction ) + eps;

        floatType xp = volume_fraction;
        floatType xm = volume_fraction;

        xp += delta;
        xm -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), xp, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), xm, internal_heat_generation,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdPhi == grad );

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( volume_fraction ) + eps;

        floatType xp = volume_fraction;
        floatType xm = volume_fraction;

        xp += delta;
        xm -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, xp,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, xm,
                                                                                          std::begin( net_interphase_force ), std::end( net_interphase_force ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdr == grad );

    }

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::abs( net_interphase_force[ i ] ) + eps;

        floatVector xp = net_interphase_force;
        floatVector xm = net_interphase_force;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( xp ), std::end( xp ), vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
                                                                                          internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
                                                                                          std::begin( xm ), std::end( xm ), vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( dRdpi[ i ] == grad );

    }

}

BOOST_AUTO_TEST_CASE( test_multiphase_computeBalanceOfEnergyNonDivergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    constexpr unsigned int nphases = 5;

    constexpr unsigned int dim = 3;

    constexpr unsigned int sot_dim = dim * dim;

    std::array<floatType,nphases> density = { 0.69646919, 0.28613933, 0.22685145, 0.55131477, 0.71946897 };

    std::array<floatType,nphases> density_dot = { 0.42310646, 0.9807642 , 0.68482974, 0.4809319 , 0.39211752 };

    std::array<floatType,dim*nphases> density_gradient = { 0.34317802, 0.72904971, 0.43857224, 0.0596779 , 0.39804426,
                                                           0.73799541, 0.18249173, 0.17545176, 0.53155137, 0.53182759,
                                                           0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338 };

    std::array<floatType,nphases> internal_energy = { 0.32295891, 0.36178866, 0.22826323, 0.29371405, 0.63097612 };

    std::array<floatType,nphases> internal_energy_dot = { 0.09210494, 0.43370117, 0.43086276, 0.4936851 , 0.42583029 };

    std::array<floatType,dim*nphases> internal_energy_gradient = { 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668,
                                                                   0.62395295, 0.1156184 , 0.31728548, 0.41482621, 0.86630916,
                                                                   0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453 };

    std::array<floatType,dim*nphases> velocity = { 0.12062867, 0.8263408 , 0.60306013, 0.54506801, 0.34276383,
                                                   0.30412079, 0.41702221, 0.68130077, 0.87545684, 0.51042234,
                                                   0.66931378, 0.58593655, 0.6249035 , 0.67468905, 0.84234244 };

    std::array<floatType,dim*dim*nphases> velocity_gradient = { 0.08319499, 0.76368284, 0.24366637, 0.19422296, 0.57245696,
                                                                0.09571252, 0.88532683, 0.62724897, 0.72341636, 0.01612921,
                                                                0.59443188, 0.55678519, 0.15895964, 0.15307052, 0.69552953,
                                                                0.31876643, 0.6919703 , 0.55438325, 0.38895057, 0.92513249,
                                                                0.84167   , 0.35739757, 0.04359146, 0.30476807, 0.39818568,
                                                                0.70495883, 0.99535848, 0.35591487, 0.76254781, 0.59317692,
                                                                0.6917018 , 0.15112745, 0.39887629, 0.2408559 , 0.34345601,
                                                                0.51312815, 0.66662455, 0.10590849, 0.13089495, 0.32198061,
                                                                0.66156434, 0.84650623, 0.55325734, 0.85445249, 0.38483781 };

    std::array<floatType,dim*dim*nphases> cauchy_stress = { 0.3167879 , 0.35426468, 0.17108183, 0.82911263, 0.33867085,
                                                            0.55237008, 0.57855147, 0.52153306, 0.00268806, 0.98834542,
                                                            0.90534158, 0.20763586, 0.29248941, 0.52001015, 0.90191137,
                                                            0.98363088, 0.25754206, 0.56435904, 0.80696868, 0.39437005,
                                                            0.73107304, 0.16106901, 0.60069857, 0.86586446, 0.98352161,
                                                            0.07936579, 0.42834727, 0.20454286, 0.45063649, 0.54776357,
                                                            0.09332671, 0.29686078, 0.92758424, 0.56900373, 0.457412  ,
                                                            0.75352599, 0.74186215, 0.04857903, 0.7086974 , 0.83924335,
                                                            0.16593788, 0.78099794, 0.28653662, 0.30646975, 0.66526147 };

    std::array<floatType,nphases> volume_fraction = { 0.11139217, 0.66487245, 0.88785679, 0.69631127, 0.44032788 };

    std::array<floatType,nphases> internal_heat_generation = { 0.43821438, 0.7650961 , 0.565642  , 0.08490416, 0.58267109 };

    std::array<floatType,dim*nphases> net_interphase_force = { 0.8148437 , 0.33706638, 0.92757658, 0.750717  , 0.57406383,
                                                               0.75164399, 0.07914896, 0.85938908, 0.82150411, 0.90987166,
                                                               0.1286312 , 0.08178009, 0.13841557, 0.39937871, 0.42430686 };

    std::array<floatType,nphases> answer = { 0.68657614, -0.26333604, -1.83670758, -0.49695816,  0.19080456 };

    std::array<floatType,nphases> result;

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                      std::begin( density_dot ),              std::end( density_dot ),
                                                                                      std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                      std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                      std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                      std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                      std::begin( velocity ),                 std::end( velocity ),
                                                                                      std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                      std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                      std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                      std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                      std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                      std::begin( result ),                   std::end( result ) );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0. );

    std::array<floatType,nphases> dRdRho, dRdRhoDot, dRdE, dRdEDot, dRdPhi, dRdr;

    std::array<floatType,dim*nphases> dRdGradRho, dRdGradE, dRdV, dRdpi;

    std::array<floatType,sot_dim*nphases> dRdGradV, dRdCauchy;

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                      std::begin( density_dot ),              std::end( density_dot ),
                                                                                      std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                      std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                      std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                      std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                      std::begin( velocity ),                 std::end( velocity ),
                                                                                      std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                      std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                      std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                      std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                      std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                      std::begin( result ),                   std::end( result ),
                                                                                      std::begin( dRdRho ),                   std::end( dRdRho ),
                                                                                      std::begin( dRdRhoDot ),                std::end( dRdRhoDot ),
                                                                                      std::begin( dRdGradRho ),               std::end( dRdGradRho ),
                                                                                      std::begin( dRdE ),                     std::end( dRdE ),
                                                                                      std::begin( dRdEDot ),                  std::end( dRdEDot ),
                                                                                      std::begin( dRdGradE ),                 std::end( dRdGradE ),
                                                                                      std::begin( dRdV ),                     std::end( dRdV ),
                                                                                      std::begin( dRdGradV ),                 std::end( dRdGradV ),
                                                                                      std::begin( dRdCauchy ),                std::end( dRdCauchy ),
                                                                                      std::begin( dRdPhi ),                   std::end( dRdPhi ),
                                                                                      std::begin( dRdr ),                     std::end( dRdr ),
                                                                                      std::begin( dRdpi ),                    std::end( dRdpi ) );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < nphases; i++ ){

        floatType delta = eps * std::abs( density[ i ] ) + eps;

        std::array<floatType,nphases> xp = density;
        std::array<floatType,nphases> xm = density;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( i == j ){

                BOOST_TEST( dRdRho[ i ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < nphases; i++ ){

        floatType delta = eps * std::abs( density_dot[ i ] ) + eps;

        std::array<floatType,nphases> xp = density_dot;
        std::array<floatType,nphases> xm = density_dot;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( i == j ){

                BOOST_TEST( dRdRhoDot[ i ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < dim*nphases; i++ ){

        floatType delta = eps * std::abs( density_gradient[ i ] ) + eps;

        std::array<floatType,dim*nphases> xp = density_gradient;
        std::array<floatType,dim*nphases> xm = density_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i / dim ) == j ){

                unsigned int phase = i / dim;
                unsigned int row   = j;
                unsigned int col   = ( i - phase * dim );

                BOOST_TEST( dRdGradRho[ dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < nphases; i++ ){

        floatType delta = eps * std::abs( internal_energy[ i ] ) + eps;

        std::array<floatType,nphases> xp = internal_energy;
        std::array<floatType,nphases> xm = internal_energy;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( i == j ){

                BOOST_TEST( dRdE[ i ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < nphases; i++ ){

        floatType delta = eps * std::abs( internal_energy_dot[ i ] ) + eps;

        std::array<floatType,nphases> xp = internal_energy_dot;
        std::array<floatType,nphases> xm = internal_energy_dot;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( i == j ){

                BOOST_TEST( dRdEDot[ i ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < dim*nphases; i++ ){

        floatType delta = eps * std::abs( internal_energy_gradient[ i ] ) + eps;

        std::array<floatType,dim*nphases> xp = internal_energy_gradient;
        std::array<floatType,dim*nphases> xm = internal_energy_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i / dim ) == j ){

                unsigned int phase = i / dim;
                unsigned int row   = j;
                unsigned int col   = ( i - phase * dim );

                BOOST_TEST( dRdGradE[ dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < dim*nphases; i++ ){

        floatType delta = eps * std::abs( velocity[ i ] ) + eps;

        std::array<floatType,dim*nphases> xp = velocity;
        std::array<floatType,dim*nphases> xm = velocity;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i / dim ) == j ){

                unsigned int phase = i / dim;
                unsigned int row   = j;
                unsigned int col   = ( i - phase * dim );

                BOOST_TEST( dRdV[ dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < sot_dim*nphases; i++ ){

        floatType delta = eps * std::abs( velocity_gradient[ i ] ) + eps;

        std::array<floatType,sot_dim*nphases> xp = velocity_gradient;
        std::array<floatType,sot_dim*nphases> xm = velocity_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i / sot_dim ) == j ){

                unsigned int phase = i / sot_dim;
                unsigned int row   = j;
                unsigned int col   = ( i - phase * sot_dim );

                BOOST_TEST( dRdGradV[ sot_dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < sot_dim*nphases; i++ ){

        floatType delta = eps * std::abs( cauchy_stress[ i ] ) + eps;

        std::array<floatType,sot_dim*nphases> xp = cauchy_stress;
        std::array<floatType,sot_dim*nphases> xm = cauchy_stress;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i / sot_dim ) == j ){

                unsigned int phase = i / sot_dim;
                unsigned int row   = j;
                unsigned int col   = ( i - phase * sot_dim );

                BOOST_TEST( dRdCauchy[ sot_dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < nphases; i++ ){

        floatType delta = eps * std::abs( volume_fraction[ i ] ) + eps;

        std::array<floatType,nphases> xp = volume_fraction;
        std::array<floatType,nphases> xm = volume_fraction;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i ) == j ){

                unsigned int row   = j;

                BOOST_TEST( dRdPhi[ row ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < nphases; i++ ){

        floatType delta = eps * std::abs( internal_heat_generation[ i ] ) + eps;

        std::array<floatType,nphases> xp = internal_heat_generation;
        std::array<floatType,nphases> xm = internal_heat_generation;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( net_interphase_force ),     std::end( net_interphase_force ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i ) == j ){

                unsigned int row   = j;

                BOOST_TEST( dRdr[ row ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < dim*nphases; i++ ){

        floatType delta = eps * std::abs( velocity[ i ] ) + eps;

        std::array<floatType,dim*nphases> xp = velocity;
        std::array<floatType,dim*nphases> xm = velocity;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( xp ),                       std::end( xp ),
                                                                                          std::begin( vp ),                       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence( std::begin( density ),                  std::end( density ),
                                                                                          std::begin( density_dot ),              std::end( density_dot ),
                                                                                          std::begin( density_gradient ),         std::end( density_gradient ),
                                                                                          std::begin( internal_energy ),          std::end( internal_energy ),
                                                                                          std::begin( internal_energy_dot ),      std::end( internal_energy_dot ),
                                                                                          std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
                                                                                          std::begin( velocity ),                 std::end( velocity ),
                                                                                          std::begin( velocity_gradient ),        std::end( velocity_gradient ),
                                                                                          std::begin( cauchy_stress ),            std::end( cauchy_stress ),
                                                                                          std::begin( volume_fraction ),          std::end( volume_fraction ),
                                                                                          std::begin( internal_heat_generation ), std::end( internal_heat_generation ),
                                                                                          std::begin( xm ),                       std::end( xm ),
                                                                                          std::begin( vm ),                       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( i / dim ) == j ){

                unsigned int phase = i / dim;
                unsigned int row   = j;
                unsigned int col   = ( i - phase * dim );

                BOOST_TEST( dRdpi[ dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

}
