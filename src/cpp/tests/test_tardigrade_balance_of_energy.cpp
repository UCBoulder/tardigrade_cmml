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
