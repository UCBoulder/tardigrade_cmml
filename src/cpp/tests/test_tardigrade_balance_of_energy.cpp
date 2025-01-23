/**
  * \file test_tardigrade_balance_equations_balance_of_energy.cpp
  *
  * Tests for tardigrade_balance_equations_balance_of_energy
  */

#include<tardigrade_balance_of_energy.h>
#include<tardigrade_finite_element_utilities.h>
#define USE_EIGEN
#include<tardigrade_vector_tools.h>
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

typedef double floatType; //!< Define the float type

typedef std::array< floatType, 3 > floatVector; //!< Define the float vector type

typedef std::array< floatType, 9 > secondOrderTensor; //!< Define the second order tensor type

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

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
        density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
        internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
        std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
        std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
        std::begin( net_interphase_force ), std::end( net_interphase_force ), result
    );

    BOOST_TEST( result == answer );

    result = 0;

    floatType dRdRho, dRdRhoDot, dRdE, dRdEDot, dRdPhi, dRdr;

    floatVector dRdGradRho, dRdGradE, dRdV, dRdpi;

    secondOrderTensor dRdGradV, dRdCauchy;

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
        density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
        internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
        std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
        std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
        std::begin( net_interphase_force ), std::end( net_interphase_force ), result,
        dRdRho, dRdRhoDot, std::begin( dRdGradRho ), std::end( dRdGradRho ),
        dRdE, dRdEDot, std::begin( dRdGradE ), std::end( dRdGradE ),
        std::begin( dRdV ), std::end( dRdV ), std::begin( dRdGradV ), std::end( dRdGradV ),
        std::begin( dRdCauchy ), std::end( dRdCauchy ), dRdPhi, dRdr,
        std::begin( dRdpi ), std::end( dRdpi )
    );

    BOOST_TEST( result == answer );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::abs( density ) + eps;

        floatType xp = density;
        floatType xm = density;

        xp += delta;
        xm -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            xp, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            xm, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, xp, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, xm, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( xp ), std::end( xp ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( xm ), std::end( xm ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            xp, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            xm, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, xp, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, xm, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( xp ), std::end( xp ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( xm ), std::end( xm ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( xp ), std::end( xp ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( xm ), std::end( xm ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( xp ), std::end( xp ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( xm ), std::end( xm ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( xp ), std::end( xp ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( xm ), std::end( xm ), volume_fraction, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), xp, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), xm, internal_heat_generation,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
             density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
             internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
             std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
             std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, xp,
             std::begin( net_interphase_force ), std::end( net_interphase_force ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, xm,
            std::begin( net_interphase_force ), std::end( net_interphase_force ), vm
        );

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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( xp ), std::end( xp ), vp
        );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            internal_energy, internal_energy_dot, std::begin( internal_energy_gradient ), std::end( internal_energy_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ), volume_fraction, internal_heat_generation,
            std::begin( xm ), std::end( xm ), vm
        );

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

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( xp ),                       std::end( xp ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( xm ),                       std::end( xm ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyNonDivergence<dim>( std::begin( density ),                  std::end( density ),
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

BOOST_AUTO_TEST_CASE( test_computeBalanceOfEnergyDivergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    constexpr unsigned int dim = 3;

    floatVector grad_psi = { 0.69646919, 0.28613933, 0.22685145 };

    floatVector q = { 0.55131477, 0.71946897, 0.42310646 };

    floatType answer = -0.68582444;

    floatType result;

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( grad_psi ), std::end( grad_psi ),
                                                                                        std::begin( q ),        std::end( q ),
                                                                                        result );

    BOOST_TEST( result == answer );

    result = 0;

    floatVector dRdGradPsi, dRdq;

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( grad_psi ),   std::end( grad_psi ),
                                                                                        std::begin( q ),          std::end( q ),
                                                                                        result,
                                                                                        std::begin( dRdGradPsi ), std::end( dRdGradPsi ),
                                                                                        std::begin( dRdq ),       std::end( dRdq ) );

    BOOST_TEST( result == answer );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::fabs( grad_psi[ i ] ) + eps;

        floatVector xp = grad_psi;
        floatVector xm = grad_psi;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( xp ),       std::end( xp ),
                                                                                            std::begin( q ),        std::end( q ),
                                                                                            vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( xm ),       std::end( xm ),
                                                                                            std::begin( q ),        std::end( q ),
                                                                                            vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( grad == dRdGradPsi[ i ] );

    }

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::fabs( q[ i ] ) + eps;

        floatVector xp = q;
        floatVector xm = q;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatType vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( grad_psi ), std::end( grad_psi ),
                                                                                            std::begin( xp ),       std::end( xp ),
                                                                                            vp );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( grad_psi ), std::end( grad_psi ),
                                                                                            std::begin( xm ),       std::end( xm ),
                                                                                            vm );

        floatType grad = ( vp - vm ) / ( 2 * delta );

        BOOST_TEST( grad == dRdq[ i ] );

    }

}

BOOST_AUTO_TEST_CASE( test_multiphase_computeBalanceOfEnergyDivergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    constexpr unsigned int dim = 3;

    constexpr unsigned int nphases = 5;

    floatVector grad_psi = { 0.69646919, 0.28613933, 0.22685145 };

    std::array<floatType,dim*nphases> q = { 0.55131477, 0.71946897, 0.42310646, 0.9807642 , 0.68482974,
                                            0.4809319 , 0.39211752, 0.34317802, 0.72904971, 0.43857224,
                                            0.0596779 , 0.39804426, 0.73799541, 0.18249173, 0.17545176 };

    std::array<floatType,nphases> answer = { -0.68582444, -0.98812887, -0.53668048, -0.41282517, -0.60601061 };

    std::array<floatType,nphases> result;

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( grad_psi ), std::end( grad_psi ),
                                                                                        std::begin( q ),        std::end( q ),
                                                                                        std::begin( result ),   std::end( result ) );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0. );

    std::array<floatType,dim*nphases> dRdGradPsi, dRdq;

    tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( grad_psi ),   std::end( grad_psi ),
                                                                                        std::begin( q ),          std::end( q ),
                                                                                        std::begin( result ),     std::end( result ),
                                                                                        std::begin( dRdGradPsi ), std::end( dRdGradPsi ),
                                                                                        std::begin( dRdq ),       std::end( dRdq ) );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::fabs( grad_psi[ i ] ) + eps;

        floatVector xp = grad_psi;
        floatVector xm = grad_psi;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( xp ),       std::end( xp ),
                                                                                            std::begin( q ),        std::end( q ),
                                                                                            std::begin( vp ),       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( xm ),       std::end( xm ),
                                                                                            std::begin( q ),        std::end( q ),
                                                                                            std::begin( vm ),       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( grad == dRdGradPsi[ dim * j + i ] );

        }

    }

    for ( unsigned int i = 0; i < dim * nphases; i++ ){

        floatType delta = eps * std::fabs( q[ i ] ) + eps;

        std::array<floatType,dim*nphases> xp = q;
        std::array<floatType,dim*nphases> xm = q;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( grad_psi ), std::end( grad_psi ),
                                                                                            std::begin( xp ),       std::end( xp ),
                                                                                            std::begin( vp ),       std::end( vp ) );

        tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergyDivergence<dim>( std::begin( grad_psi ), std::end( grad_psi ),
                                                                                            std::begin( xm ),       std::end( xm ),
                                                                                            std::begin( vm ),       std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( i / dim ) == j ){

                unsigned int row = ( i - dim * j );

                BOOST_TEST( grad == dRdq[ dim * j + row ] );

            }
            else{

                BOOST_TEST( grad == 0 );

            }

        }

    }

}

template<
  typename dt_type, class v_t_in, class v_tp1_in,
  class vDot_t_in, typename alpha_type, class vDot_tp1_out,
  typename dVDotdV_type
>
void compute_current_rate_of_change(
    const dt_type &dt,
    const v_t_in &v_t_begin, const v_t_in &v_t_end,
    const v_tp1_in &v_tp1_begin, const v_tp1_in &v_tp1_end,
    const vDot_t_in &vDot_t_begin, const vDot_t_in &vDot_t_end,
    const alpha_type &alpha, vDot_tp1_out vDot_tp1_begin, vDot_tp1_out vDot_tp1_end,
    dVDotdV_type &dVDotdV
){

    dVDotdV = 1. / ( alpha * dt );

    for ( unsigned int i = 0; i < ( unsigned int )( v_t_end - v_t_begin ); ++i ){

        *( vDot_tp1_begin + i ) = ( ( *( v_tp1_begin + i ) ) - ( *( v_t_begin + i ) ) ) / ( alpha * dt ) - ( ( 1 - alpha ) / alpha ) * ( *( vDot_t_begin + i ) );

    }

}

template<
    int dim, int node_count, int nphases,
    class xi_in, typename dt_type, class density_t_in, class density_tp1_in, class e_t_in, class e_tp1_in,
    class u_t_in, class u_tp1_in, class umesh_t_in, class umesh_tp1_in, class density_dot_t_in, class e_dot_t_in, class u_dot_t_in,
    class X_in, class cauchy_stress_iter, class heat_flux_iter, class volume_fraction_iter, class internal_heat_generation_iter,
    class net_interphase_force_iter,
    typename alpha_type, typename beta_type, class value_out
>
void evaluate_at_nodes(
    const xi_in &xi_begin, const xi_in &xi_end, dt_type dt,
    const density_t_in &density_t_begin, const density_t_in &density_t_end,
    const density_tp1_in &density_tp1_begin, const density_tp1_in &density_tp1_end,
    const e_t_in &e_t_begin, const e_t_in &e_t_end,
    const e_tp1_in &e_tp1_begin, const e_tp1_in &e_tp1_end,
    const u_t_in &u_t_begin, const u_t_in &u_t_end,
    const u_tp1_in &u_tp1_begin, const u_tp1_in &u_tp1_end,
    const umesh_t_in &umesh_t_begin, const umesh_t_in &umesh_t_end,
    const umesh_tp1_in &umesh_tp1_begin, const umesh_tp1_in &umesh_tp1_end,
    const density_dot_t_in &density_dot_t_begin, const density_dot_t_in &density_dot_t_end,
    const e_dot_t_in &e_dot_t_begin, const e_dot_t_in &e_dot_t_end,
    const u_dot_t_in &u_dot_t_begin, const u_dot_t_in &u_dot_t_end,
    const X_in &X_begin, const X_in &X_end,
    const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
    const heat_flux_iter &heat_flux_begin, const heat_flux_iter &heat_flux_end,
    const volume_fraction_iter &volume_fraction_begin, const volume_fraction_iter &volume_fraction_end,
    const internal_heat_generation_iter &internal_heat_generation_begin, const internal_heat_generation_iter &internal_heat_generation_end,
    const net_interphase_force_iter &net_interphase_force_begin, const net_interphase_force_iter &net_interphase_force_end,
    const alpha_type &alpha, const beta_type &beta, value_out value_begin, value_out value_end
){

    // Update the mesh nodes
    std::array< typename std::iterator_traits<umesh_t_in  >::value_type, dim * node_count > x_t;
    std::array< typename std::iterator_traits<umesh_tp1_in>::value_type, dim * node_count > x_tp1;

    std::transform( X_begin, X_end,   umesh_t_begin,   std::begin( x_t ), std::plus<typename std::iterator_traits<umesh_t_in  >::value_type>( ) );
    std::transform( X_begin, X_end, umesh_tp1_begin, std::begin( x_tp1 ), std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>( ) );

    // Calculate the current rates of change
    std::array< typename std::iterator_traits<density_tp1_in>::value_type, node_count * nphases > density_dot_tp1;
    std::array< typename std::iterator_traits<e_tp1_in>::value_type, node_count * nphases > e_dot_tp1;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count * nphases > u_dot_tp1;

    floatType dRhoDotdRho, dEDotdE, dUDotdU;

    compute_current_rate_of_change(
        dt, density_t_begin, density_t_end, density_tp1_begin, density_tp1_end,
        density_dot_t_begin, density_dot_t_end, alpha,
        std::begin( density_dot_tp1 ), std::end( density_dot_tp1 ),
        dRhoDotdRho
    );

    compute_current_rate_of_change(
        dt, e_t_begin, e_t_end, e_tp1_begin, e_tp1_end,
        e_dot_t_begin, e_dot_t_end, alpha,
        std::begin( e_dot_tp1 ), std::end( e_dot_tp1 ),
        dEDotdE
    );

    compute_current_rate_of_change(
        dt, u_t_begin, u_t_end, u_tp1_begin, u_tp1_end,
        u_dot_t_begin, u_dot_t_end, alpha,
        std::begin( u_dot_tp1 ), std::end( u_dot_tp1 ),
        dUDotdU
    );

    // Instantiate the element
    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<
        floatType, typename std::array< floatType, 24 >::const_iterator,
        typename std::array< floatType, 3>::const_iterator,
        typename std::array< floatType, 8>::iterator,
        typename std::array< floatType, 24>::iterator
    > e(
        std::cbegin( x_tp1 ), std::cend( x_tp1 ), X_begin, X_end
    );

    std::array<
        typename std::iterator_traits<density_tp1_in>::value_type, nphases
    > density_tp1_p, density_dot_tp1_p;

    std::array<
        typename std::iterator_traits<e_tp1_in>::value_type, nphases
    > e_tp1_p, e_dot_tp1_p;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * nphases
    > u_dot_tp1_p;

    // Interpolate quantities to the local point

    e.InterpolateQuantity(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( density_tp1_p ), std::end( density_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::cbegin( density_dot_tp1 ), std::cend( density_dot_tp1 ),
        std::begin( density_dot_tp1_p ), std::end( density_dot_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, e_tp1_begin, e_tp1_end,
        std::begin( e_tp1_p ), std::end( e_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::cbegin( e_dot_tp1 ), std::cend( e_dot_tp1 ),
        std::begin( e_dot_tp1_p ), std::end( e_dot_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::cbegin( u_dot_tp1 ), std::cend( u_dot_tp1 ),
        std::begin( u_dot_tp1_p ), std::end( u_dot_tp1_p )
    );

    // Compute the gradients at the local point

    std::array<
        typename std::iterator_traits<density_tp1_in>::value_type, dim * nphases
    > grad_density_tp1;

    std::array<
        typename std::iterator_traits<e_tp1_in>::value_type, dim * nphases
    > grad_e_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * dim * nphases
    > grad_u_dot_tp1;

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( grad_density_tp1 ), std::end( grad_density_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, e_tp1_begin, e_tp1_end,
        std::begin( grad_e_tp1 ), std::end( grad_e_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, std::cbegin( u_dot_tp1 ), std::cend( u_dot_tp1 ),
        std::begin( grad_u_dot_tp1 ), std::end( grad_u_dot_tp1 )
    );

    // Get the Jacobian of transformation
    std::array< floatType, dim * dim > dxdxi;
    e.GetLocalQuantityGradient(
        xi_begin, xi_end, std::cbegin( x_tp1 ), std::cend( x_tp1 ),
        std::begin( dxdxi ), std::end( dxdxi )
    );

    floatType J = tardigradeVectorTools::determinant
    <
        typename std::array< floatType, dim * dim >::const_iterator,
        floatType, 3, 3
    >(
        std::cbegin( dxdxi ), std::cend( dxdxi ),
        3, 3
    );

    std::array< floatType, node_count> Ns;
    e.GetShapeFunctions( xi_begin, xi_end, std::begin( Ns ), std::end( Ns ) );

    std::array< floatType, node_count * dim > dNdx;
    e.GetGlobalShapeFunctionGradients( xi_begin, xi_end, std::cbegin( x_tp1 ), std::cend( x_tp1 ),
                                       std::begin( dNdx ), std::end( dNdx ) );

    if ( nphases == 1 ){

        for ( unsigned int i = 0; i < node_count; ++i ){

            tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<dim>(
                density_tp1_p[ 0 ], density_dot_tp1_p[ 0 ],
                std::cbegin( grad_density_tp1 ), std::cend( grad_density_tp1 ),
                e_tp1_p[ 0 ], e_dot_tp1_p[ 0 ],
                std::cbegin( grad_e_tp1 ), std::cend( grad_e_tp1 ),
                std::cbegin( u_dot_tp1_p ), std::cend( u_dot_tp1_p ),
                std::cbegin( grad_u_dot_tp1 ), std::cend( grad_u_dot_tp1 ),
                cauchy_stress_begin, cauchy_stress_end,
                *volume_fraction_begin,
                *internal_heat_generation_begin,
                net_interphase_force_begin, net_interphase_force_end,
                heat_flux_begin, heat_flux_end,
                Ns[ i ],
                std::begin( dNdx ) + dim * i, std::begin( dNdx ) + dim * ( i + 1 ),
                *( value_begin + i )
            );

            std::transform(
                value_begin + i, value_begin + ( i + 1 ), value_begin + i,
                std::bind(
                    std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                    std::placeholders::_1,
                    J
                )
            );

        }

    }
    else{

        for ( unsigned int i = 0; i < node_count; ++i ){

            tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<dim>(
                std::cbegin( density_tp1_p ),     std::cend( density_tp1_p ),
                std::cbegin( density_dot_tp1_p ), std::cend( density_dot_tp1_p ),
                std::cbegin( grad_density_tp1 ),  std::cend( grad_density_tp1 ),
                std::cbegin( e_tp1_p ),           std::cend( e_tp1_p ),
                std::cbegin( e_dot_tp1_p ),       std::cend( e_dot_tp1_p ),
                std::cbegin( grad_e_tp1 ),        std::cend( grad_e_tp1 ),
                std::cbegin( u_dot_tp1_p ),       std::cend( u_dot_tp1_p ),
                std::cbegin( grad_u_dot_tp1 ),    std::cend( grad_u_dot_tp1 ),
                cauchy_stress_begin,              cauchy_stress_end,
                volume_fraction_begin,            volume_fraction_end,
                internal_heat_generation_begin,   internal_heat_generation_end,
                net_interphase_force_begin,       net_interphase_force_end,
                heat_flux_begin,                  heat_flux_end,
                Ns[ i ],
                std::begin( dNdx ) + dim * i, std::begin( dNdx ) + 3 * ( i + 1 ),
                value_begin + nphases * i, value_begin + nphases * ( i + 1 )
            );

            std::transform(
                value_begin + nphases * i, value_begin + nphases * ( i + 1 ), value_begin + nphases * i,
                std::bind(
                    std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                    std::placeholders::_1,
                    J
                )
            );

        }

    }

}

template<
    int dim, int node_count, int nphases,
    class xi_in, typename dt_type, class density_t_in, class density_tp1_in, class e_t_in, class e_tp1_in,
    class u_t_in, class u_tp1_in, class umesh_t_in, class umesh_tp1_in, class density_dot_t_in, class e_dot_t_in, class u_dot_t_in,
    class X_in, class cauchy_stress_iter, class heat_flux_iter, class volume_fraction_iter, class internal_heat_generation_iter,
    class net_interphase_force_iter,
    typename alpha_type, typename beta_type, class value_out,
    class dRdRho_iter, class dRdE_iter, class dRdU_iter, class dRdCauchy_iter, class dRdq_iter, class dRdVolumeFraction_iter, class dRdr_iter,
    class dRdpi_iter, class dRdUMesh_iter
>
void evaluate_at_nodes(
    const xi_in &xi_begin, const xi_in &xi_end, dt_type dt,
    const density_t_in &density_t_begin, const density_t_in &density_t_end,
    const density_tp1_in &density_tp1_begin, const density_tp1_in &density_tp1_end,
    const e_t_in &e_t_begin, const e_t_in &e_t_end,
    const e_tp1_in &e_tp1_begin, const e_tp1_in &e_tp1_end,
    const u_t_in &u_t_begin, const u_t_in &u_t_end,
    const u_tp1_in &u_tp1_begin, const u_tp1_in &u_tp1_end,
    const umesh_t_in &umesh_t_begin, const umesh_t_in &umesh_t_end,
    const umesh_tp1_in &umesh_tp1_begin, const umesh_tp1_in &umesh_tp1_end,
    const density_dot_t_in &density_dot_t_begin, const density_dot_t_in &density_dot_t_end,
    const e_dot_t_in &e_dot_t_begin, const e_dot_t_in &e_dot_t_end,
    const u_dot_t_in &u_dot_t_begin, const u_dot_t_in &u_dot_t_end,
    const X_in &X_begin, const X_in &X_end,
    const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
    const heat_flux_iter &heat_flux_begin, const heat_flux_iter &heat_flux_end,
    const volume_fraction_iter &volume_fraction_begin, const volume_fraction_iter &volume_fraction_end,
    const internal_heat_generation_iter &internal_heat_generation_begin, const internal_heat_generation_iter &internal_heat_generation_end,
    const net_interphase_force_iter &net_interphase_force_begin, const net_interphase_force_iter &net_interphase_force_end,
    const alpha_type &alpha, const beta_type &beta, value_out value_begin, value_out value_end,
    dRdRho_iter dRdRho_begin,                       dRdRho_iter dRdRho_end,
    dRdE_iter dRdE_begin,                           dRdE_iter dRdE_end,
    dRdU_iter dRdU_begin,                           dRdU_iter dRdU_end,
    dRdCauchy_iter dRdCauchy_begin,                 dRdCauchy_iter dRdCauchy_end,
    dRdq_iter dRdq_begin,                           dRdq_iter dRdq_end,
    dRdVolumeFraction_iter dRdVolumeFraction_begin, dRdVolumeFraction_iter dRdVolumeFraction_end,
    dRdr_iter dRdr_begin,                           dRdr_iter dRdr_end,
    dRdpi_iter dRdpi_begin,                         dRdpi_iter dRdpi_end,
    dRdUMesh_iter dRdUMesh_begin,                   dRdUMesh_iter dRdUMesh_end
){

    // Update the mesh nodes
    std::array< typename std::iterator_traits<umesh_t_in  >::value_type, dim * node_count > x_t;
    std::array< typename std::iterator_traits<umesh_tp1_in>::value_type, dim * node_count > x_tp1;

    std::transform( X_begin, X_end,   umesh_t_begin,   std::begin( x_t ), std::plus<typename std::iterator_traits<umesh_t_in  >::value_type>( ) );
    std::transform( X_begin, X_end, umesh_tp1_begin, std::begin( x_tp1 ), std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>( ) );

    // Calculate the current rates of change
    std::array< typename std::iterator_traits<density_tp1_in>::value_type, node_count * nphases > density_dot_tp1;
    std::array< typename std::iterator_traits<e_tp1_in>::value_type, node_count * nphases > e_dot_tp1;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count * nphases > u_dot_tp1;

    floatType dRhoDotdRho, dEDotdE, dUDotdU;

    compute_current_rate_of_change(
        dt, density_t_begin, density_t_end, density_tp1_begin, density_tp1_end,
        density_dot_t_begin, density_dot_t_end, alpha,
        std::begin( density_dot_tp1 ), std::end( density_dot_tp1 ),
        dRhoDotdRho
    );

    compute_current_rate_of_change(
        dt, e_t_begin, e_t_end, e_tp1_begin, e_tp1_end,
        e_dot_t_begin, e_dot_t_end, alpha,
        std::begin( e_dot_tp1 ), std::end( e_dot_tp1 ),
        dEDotdE
    );

    compute_current_rate_of_change(
        dt, u_t_begin, u_t_end, u_tp1_begin, u_tp1_end,
        u_dot_t_begin, u_dot_t_end, alpha,
        std::begin( u_dot_tp1 ), std::end( u_dot_tp1 ),
        dUDotdU
    );

    // Instantiate the element
    tardigradeBalanceEquations::finiteElementUtilities::LinearHex<
        floatType, typename std::array< floatType, 24 >::const_iterator,
        typename std::array< floatType, 3>::const_iterator,
        typename std::array< floatType, 8>::iterator,
        typename std::array< floatType, 24>::iterator
    > e(
        std::cbegin( x_tp1 ), std::cend( x_tp1 ), X_begin, X_end
    );

    std::array<
        typename std::iterator_traits<density_tp1_in>::value_type, nphases
    > density_tp1_p, density_dot_tp1_p;

    std::array<
        typename std::iterator_traits<e_tp1_in>::value_type, nphases
    > e_tp1_p, e_dot_tp1_p;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * nphases
    > u_dot_tp1_p;

    // Interpolate quantities to the local point

    e.InterpolateQuantity(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( density_tp1_p ), std::end( density_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::cbegin( density_dot_tp1 ), std::cend( density_dot_tp1 ),
        std::begin( density_dot_tp1_p ), std::end( density_dot_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, e_tp1_begin, e_tp1_end,
        std::begin( e_tp1_p ), std::end( e_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::cbegin( e_dot_tp1 ), std::cend( e_dot_tp1 ),
        std::begin( e_dot_tp1_p ), std::end( e_dot_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::cbegin( u_dot_tp1 ), std::cend( u_dot_tp1 ),
        std::begin( u_dot_tp1_p ), std::end( u_dot_tp1_p )
    );

    // Compute the gradients at the local point

    std::array<
        typename std::iterator_traits<density_tp1_in>::value_type, dim * nphases
    > grad_density_tp1;

    std::array<
        typename std::iterator_traits<e_tp1_in>::value_type, dim * nphases
    > grad_e_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * dim * nphases
    > grad_u_dot_tp1;

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( grad_density_tp1 ), std::end( grad_density_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, e_tp1_begin, e_tp1_end,
        std::begin( grad_e_tp1 ), std::end( grad_e_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, std::cbegin( u_dot_tp1 ), std::cend( u_dot_tp1 ),
        std::begin( grad_u_dot_tp1 ), std::end( grad_u_dot_tp1 )
    );

    // Get the Jacobian of transformation
    std::array< floatType, dim * dim > dxdxi;
    e.GetLocalQuantityGradient(
        xi_begin, xi_end, std::cbegin( x_tp1 ), std::cend( x_tp1 ),
        std::begin( dxdxi ), std::end( dxdxi )
    );

    floatType J = tardigradeVectorTools::determinant
    <
        typename std::array< floatType, dim * dim >::const_iterator,
        floatType, 3, 3
    >(
        std::cbegin( dxdxi ), std::cend( dxdxi ),
        3, 3
    );

    std::array< floatType, node_count> Ns;
    e.GetShapeFunctions( xi_begin, xi_end, std::begin( Ns ), std::end( Ns ) );

    std::array< floatType, node_count * dim > dNdx;
    e.GetGlobalShapeFunctionGradients( xi_begin, xi_end, std::cbegin( x_tp1 ), std::cend( x_tp1 ),
                                       std::begin( dNdx ), std::end( dNdx ) );

    std::array< floatType, nphases >              value_p;
    std::array< floatType, nphases * 1   >        dRdRho_p;
    std::array< floatType, nphases * 1   >        dRdE_p;
    std::array< floatType, nphases * dim >        dRdU_p;
    std::array< floatType, nphases * dim * dim >  dRdCauchy_p;
    std::array< floatType, nphases * dim >        dRdq_p;
    std::array< floatType, nphases * 1   >        dRdVolumeFraction_p;
    std::array< floatType, nphases * 1   >        dRdr_p;
    std::array< floatType, nphases * dim >        dRdpi_p;
    std::array< floatType, nphases * dim >        dRdUMesh_p;

    std::fill( dRdRho_begin,            dRdRho_end,            0 );
    std::fill( dRdE_begin,              dRdE_end,              0 );
    std::fill( dRdU_begin,              dRdU_end,              0 );
    std::fill( dRdCauchy_begin,         dRdCauchy_end,         0 );
    std::fill( dRdq_begin,              dRdq_end,              0 );
    std::fill( dRdVolumeFraction_begin, dRdVolumeFraction_end, 0 );
    std::fill( dRdr_begin,              dRdr_end,              0 );
    std::fill( dRdpi_begin,             dRdpi_end,             0 );
    std::fill( dRdUMesh_begin,          dRdUMesh_end,          0 );

    if ( nphases == 1 ){

        for ( unsigned int i = 0; i < node_count; ++i ){

            tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<dim>(
                density_tp1_p[ 0 ], density_dot_tp1_p[ 0 ],
                std::cbegin( grad_density_tp1 ), std::cend( grad_density_tp1 ),
                e_tp1_p[ 0 ], e_dot_tp1_p[ 0 ],
                std::cbegin( grad_e_tp1 ), std::cend( grad_e_tp1 ),
                std::cbegin( u_dot_tp1_p ), std::cend( u_dot_tp1_p ),
                std::cbegin( grad_u_dot_tp1 ), std::cend( grad_u_dot_tp1 ),
                cauchy_stress_begin, cauchy_stress_end,
                *volume_fraction_begin,
                *internal_heat_generation_begin,
                net_interphase_force_begin, net_interphase_force_end,
                heat_flux_begin, heat_flux_end,
                Ns[ i ],
                std::begin( dNdx ) + dim * i, std::begin( dNdx ) + dim * ( i + 1 ),
                *( value_begin + i )
            );

            std::transform(
                value_begin + i, value_begin + ( i + 1 ), value_begin + i,
                std::bind(
                    std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                    std::placeholders::_1,
                    J
                )
            );

        }

        for ( unsigned int i = 0; i < node_count; ++i ){ // Loop over test functions

            for ( unsigned int j = 0; j < node_count; ++j ){ // Loop over interpolation functions

                tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<dim>(
                    density_tp1_p[ 0 ], density_dot_tp1_p[ 0 ],
                    std::cbegin( grad_density_tp1 ), std::cend( grad_density_tp1 ),
                    e_tp1_p[ 0 ], e_dot_tp1_p[ 0 ],
                    std::cbegin( grad_e_tp1 ), std::cend( grad_e_tp1 ),
                    std::cbegin( u_dot_tp1_p ), std::cend( u_dot_tp1_p ),
                    std::cbegin( grad_u_dot_tp1 ), std::cend( grad_u_dot_tp1 ),
                    cauchy_stress_begin, cauchy_stress_end,
                    *volume_fraction_begin,
                    *internal_heat_generation_begin,
                    net_interphase_force_begin, net_interphase_force_end,
                    heat_flux_begin, heat_flux_end,
                    Ns[ i ],
                    std::begin( dNdx ) + dim * i, std::begin( dNdx ) + dim * ( i + 1 ),
                    Ns[ j ],
                    std::begin( dNdx ) + dim * j, std::begin( dNdx ) + dim * ( j + 1 ),
                    dRhoDotdRho, dEDotdE, dUDotdU,
                    value_p[ 0 ],
                    dRdRho_p[ 0 ], dRdE_p[ 0 ],
                    std::begin( dRdU_p ),      std::end( dRdU_p ),
                    std::begin( dRdCauchy_p ), std::end( dRdCauchy_p ),
                    dRdVolumeFraction_p[ 0 ],
                    dRdr_p[ 0 ],
                    std::begin( dRdpi_p ),    std::end( dRdpi_p ),
                    std::begin( dRdq_p ),     std::end( dRdq_p ),
                    std::begin( dRdUMesh_p ), std::end( dRdUMesh_p )
                );

                value_p[ 0 ] *= J;

                BOOST_TEST( value_p[ 0 ] == ( *( value_begin + i ) ) );

                for ( unsigned int k = 0; k < 1; ++k ){

                    // node_count, 1, node_count, 1
                    // i, k, j, l
                    for ( unsigned int l = 0; l < 1; ++l ){

                        *( dRdRho_begin + 1 * node_count * 1 * i + node_count * 1 * k + 1 * j + l ) = dRdRho_p[ l ] * J;

                    }

                    for ( unsigned int l = 0; l < 1; ++l ){

                        *( dRdE_begin + 1 * node_count * 1 * i + node_count * 1 * k + 1 * j + l ) = dRdE_p[ l ] * J;

                    }

                    for ( unsigned int l = 0; l < dim; ++l ){

                        *( dRdU_begin + 1 * node_count * dim * i + node_count * dim * k + dim * j + l ) = dRdU_p[ l ] * J;

                        *( dRdUMesh_begin + 1 * node_count * dim * i + node_count * dim * k + dim * j + l ) = dRdUMesh_p[ l ] * J;

                    }

                }

            }

            // The quadrature point material model values aren't necessarily a quantity interpolated from the nodal quantities (i.e., no interpolation function)
            for ( unsigned int j = 0; j < 1; ++j ){

                for ( unsigned int k = 0; k < 1; ++k ){

                    for ( unsigned int l = 0; l < dim * dim; ++l ){

                        *( dRdCauchy_begin + 1 * 1 * dim * dim * i + 1 * dim * dim * k + dim * dim * j + l )
                            = dRdCauchy_p[ l ] * J;

                    }

                    for ( unsigned int l = 0; l < 1; ++l ){

                        *( dRdVolumeFraction_begin + 1 * 1 * 1 * i + 1 * 1 * k + 1 * j + l )
                            = dRdVolumeFraction_p[ l ] * J;

                        *( dRdr_begin + 1 * 1 * 1 * i + 1 * 1 * k + 1 * j + l )
                            = dRdr_p[ l ] * J;

                    }

                    for ( unsigned int l = 0; l < 3; ++l ){

                        *( dRdpi_begin + 1 * 1 * 3 * i + 1 * 3 * k + 3 * j + l )
                            = dRdpi_p[ l ] * J;

                        *( dRdq_begin + 1 * 1 * 3 * i + 1 * 3 * k + 3 * j + l )
                            = dRdq_p[ l ] * J;

                    }

                }

            }

        }

    }
    else{

        for ( unsigned int i = 0; i < node_count; ++i ){

            tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<dim>(
                std::cbegin( density_tp1_p ),     std::cend( density_tp1_p ),
                std::cbegin( density_dot_tp1_p ), std::cend( density_dot_tp1_p ),
                std::cbegin( grad_density_tp1 ),  std::cend( grad_density_tp1 ),
                std::cbegin( e_tp1_p ),           std::cend( e_tp1_p ),
                std::cbegin( e_dot_tp1_p ),       std::cend( e_dot_tp1_p ),
                std::cbegin( grad_e_tp1 ),        std::cend( grad_e_tp1 ),
                std::cbegin( u_dot_tp1_p ),       std::cend( u_dot_tp1_p ),
                std::cbegin( grad_u_dot_tp1 ),    std::cend( grad_u_dot_tp1 ),
                cauchy_stress_begin,              cauchy_stress_end,
                volume_fraction_begin,            volume_fraction_end,
                internal_heat_generation_begin,   internal_heat_generation_end,
                net_interphase_force_begin,       net_interphase_force_end,
                heat_flux_begin,                  heat_flux_end,
                Ns[ i ],
                std::begin( dNdx ) + dim * i, std::begin( dNdx ) + dim * ( i + 1 ),
                value_begin + nphases * i, value_begin + nphases * ( i + 1 )
            );

            std::transform(
                value_begin + nphases * i, value_begin + nphases * ( i + 1 ), value_begin + nphases * i,
                std::bind(
                    std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                    std::placeholders::_1,
                    J
                )
            );

        }

        for ( unsigned int i = 0; i < node_count; ++i ){ // Loop over test functions

            for ( unsigned int j = 0; j < node_count; ++j ){ // Loop over interpolation functions

                tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<dim>(
                    std::cbegin( density_tp1_p ),     std::cend( density_tp1_p ),
                    std::cbegin( density_dot_tp1_p ), std::cend( density_dot_tp1_p ),
                    std::cbegin( grad_density_tp1 ),  std::cend( grad_density_tp1 ),
                    std::cbegin( e_tp1_p ),           std::cend( e_tp1_p ),
                    std::cbegin( e_dot_tp1_p ),       std::cend( e_dot_tp1_p ),
                    std::cbegin( grad_e_tp1 ),        std::cend( grad_e_tp1 ),
                    std::cbegin( u_dot_tp1_p ),       std::cend( u_dot_tp1_p ),
                    std::cbegin( grad_u_dot_tp1 ),    std::cend( grad_u_dot_tp1 ),
                    cauchy_stress_begin,              cauchy_stress_end,
                    volume_fraction_begin,            volume_fraction_end,
                    internal_heat_generation_begin,   internal_heat_generation_end,
                    net_interphase_force_begin,       net_interphase_force_end,
                    heat_flux_begin,                  heat_flux_end,
                    Ns[ i ],
                    std::begin( dNdx ) + dim * i, std::begin( dNdx ) + dim * ( i + 1 ),
                    Ns[ j ],
                    std::begin( dNdx ) + dim * j, std::begin( dNdx ) + dim * ( j + 1 ),
                    dRhoDotdRho, dEDotdE, dUDotdU,
                    std::begin( value_p ),             std::end( value_p ),
                    std::begin( dRdRho_p ),            std::end( dRdRho_p ),
                    std::begin( dRdE_p ),              std::end( dRdE_p ),
                    std::begin( dRdU_p ),              std::end( dRdU_p ),
                    std::begin( dRdCauchy_p ),         std::end( dRdCauchy_p ),
                    std::begin( dRdVolumeFraction_p ), std::end( dRdVolumeFraction_p ),
                    std::begin( dRdr_p ),              std::end( dRdr_p ),
                    std::begin( dRdpi_p ),             std::end( dRdpi_p ),
                    std::begin( dRdq_p ),              std::end( dRdq_p ),
                    std::begin( dRdUMesh_p ),          std::end( dRdUMesh_p )
                );

                std::transform(
                    std::begin( value_p ), std::end( value_p ), std::begin( value_p ),
                    std::bind(
                        std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                        std::placeholders::_1,
                        J
                    )
                );

                for ( unsigned int k = 0; k < nphases; ++k ){

                    BOOST_TEST( value_p[ k ] == ( *( value_begin + nphases * i + k ) ) );

                    *( dRdRho_begin + nphases * node_count * nphases * 1 * i + node_count * nphases * 1 * k + nphases * 1 * j + k )
                        += dRdRho_p[ k ] * J;

                    *( dRdE_begin + nphases * node_count * nphases * 1 * i + node_count * nphases * 1 * k + nphases * 1 * j + k )
                        += dRdE_p[ k ] * J;

                    for ( unsigned int l = 0; l < dim; ++l ){

                        *( dRdU_begin + nphases * node_count * nphases * dim * i + node_count * nphases * dim * k + nphases * dim * j + dim * k + l )
                            += dRdU_p[ dim * k + l ] * J;

                        *( dRdUMesh_begin + nphases * node_count * dim * i + node_count * dim * k + dim * j + l )
                            += dRdUMesh_p[ dim * k + l ] * J;

                    }

                }

            }

            for ( unsigned int j = 0; j < nphases; ++j ){

                for ( unsigned int k = 0; k < dim * dim; ++k ){

                    *( dRdCauchy_begin + nphases * nphases * dim * dim * i + nphases * dim * dim * j + dim * dim * j + k )
                        += dRdCauchy_p[ dim * dim * j + k ] * J;

                }

                for ( unsigned int k = 0; k < 1; ++k ){

                    *( dRdVolumeFraction_begin + nphases * nphases * 1 * i + nphases * 1 * j + 1 * j + k )
                        += dRdVolumeFraction_p[ 1 * j + k ] * J;

                }

                for ( unsigned int k = 0; k < 1; ++k ){

                    *( dRdr_begin + nphases * nphases * 1 * i + nphases * 1 * j + 1 * j + k )
                        += dRdr_p[ 1 * j + k ] * J;

                }

                for ( unsigned int k = 0; k < dim; ++k ){

                    *( dRdq_begin + nphases * nphases * dim * i + nphases * dim * j + dim * j + k )
                        += dRdq_p[ dim * j + k ] * J;

                }

                for ( unsigned int k = 0; k < dim; ++k ){

                    *( dRdpi_begin + nphases * nphases * dim * i + nphases * dim * j + dim * j + k )
                        += dRdpi_p[ dim * j + k ] * J;

                }

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_computeBalanceOfEnergy_fea, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test computing the balance of mass in a finite element context
     */

    constexpr unsigned int nphases = 1;

    std::array< floatType, 8 > density_t = {
        0.61289453, 0.12062867, 0.8263408 , 0.60306013, 0.54506801,
        0.34276383, 0.30412079, 0.41702221
    };

    std::array< floatType, 8 > density_tp1 = {
        0.08319499, 0.76368284, 0.24366637, 0.19422296, 0.57245696,
        0.09571252, 0.88532683, 0.62724897
    };

    std::array< floatType, 8 > e_t = {
        0.72341636, 0.01612921, 0.59443188, 0.55678519, 0.15895964,
        0.15307052, 0.69552953, 0.31876643
    };

    std::array< floatType, 8 > e_tp1 = {
        0.39818568, 0.70495883, 0.99535848, 0.35591487, 0.76254781,
        0.59317692, 0.6917018 , 0.15112745
    };

    std::array< floatType, 24 > u_t = {
        -0.20224741, -0.5182882 , -0.31308797,  0.02625631,  0.3332491 ,
        -0.78818303, -0.7382101 , -0.35603879,  0.32312867,  0.69301245,
         0.10651469,  0.70890498, -0.23032438, -0.36642421, -0.29147065,
        -0.65783634,  0.65822527, -0.32265831,  0.10474015,  0.15710294,
         0.04306612, -0.99462387,  0.97669084,  0.81068315
    };

    std::array< floatType, 24 > u_tp1 = {
         0.50705198,  0.4837243 , -0.90284193,  0.41739479,  0.6784867 ,
        -0.66812423,  0.56199588, -0.42692677, -0.38706049,  0.33052293,
        -0.77721566,  0.3297449 ,  0.77571359,  0.39262254, -0.11934425,
        -0.12357123,  0.53019219,  0.131284  , -0.83019167,  0.16534218,
         0.62968741, -0.32586723,  0.85515316,  0.501434  
    };

    std::array< floatType, 24 > umesh_t = { 
        0.53182759, 0.63440096, 0.84943179, 0.51044152, 0.65634786,
        0.86791477, 0.48312667, 0.6486585 , 0.86600796, 0.50451273,
        0.6267116 , 0.84752498, 0.53695906, 0.68247738, 0.83864355,
        0.515573  , 0.70442428, 0.85712652, 0.48825814, 0.69673492,
        0.85521971, 0.50964421, 0.67478802, 0.83673674
    };

    std::array< floatType, 24 > umesh_tp1 = {
        0.72445532, 0.61102351, 0.72244338, 0.70877313, 0.5669913 ,
        0.69069256, 0.7316781 , 0.55679573, 0.65823773, 0.7473603 ,
        0.60082794, 0.68998856, 0.71831255, 0.63482305, 0.72559852,
        0.70263035, 0.59079084, 0.69384769, 0.72553532, 0.58059527,
        0.66139287, 0.74121752, 0.62462748, 0.6931437 
    };

    std::array< floatType, 8 > density_dot_t = {
        0.68130077, 0.87545684, 0.51042234, 0.66931378, 0.58593655,
        0.6249035 , 0.67468905, 0.84234244
    };

    std::array< floatType, 8 > e_dot_t = {
        0.6919703 , 0.55438325, 0.38895057, 0.92513249, 0.84167   ,
        0.35739757, 0.04359146, 0.30476807
    };

    std::array< floatType, 24 > u_dot_t = {
        -0.58472828, -0.41502117,  0.04002031,  0.80382275,  0.96726177,
        -0.48491587,  0.12871809,  0.61393737, -0.21125989,  0.46214607,
        -0.67786197,  0.20139714,  0.73172892,  0.96704322, -0.84126842,
        -0.14330545, -0.59091428, -0.09872702,  0.09552715, -0.81334658,
        -0.40627845,  0.85516848,  0.13800746, -0.085176  
    };

    std::array< floatType, 24 > X = {
        0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0., 0., 1., 1., 0.,
       1., 1., 1., 1., 0., 1., 1.
    };

    std::array< floatType, 9 > cauchy_stress = {
         0.14812765,  0.50328798, -0.84170208,  0.71877815,  0.64300823,
         0.81974332, -0.7427376 , -0.83643983, -0.72316885
    };

    std::array< floatType, 3 > heat_flux = {
        -0.20124258, -0.15138628,  0.12443676
    };

    std::array< floatType, 1 > volume_fractions = {
        0.12224355
    };

    std::array< floatType, 1 > internal_heat_generation = {
        -0.597201
    };

    std::array< floatType, 3 > net_interphase_force = {
         0.6232887 , -0.06402485,  0.61587642
    };

    std::array< floatType, 8 > answer = {
       0.0067579 , 0.0008028 , 0.00101024, 0.01134661, 0.0025969 ,
       0.01031854, 0.01300375, 0.04472838
    };

    std::array< floatType, 8 > result;

    std::array< floatType, 3 > local_point = {
        -0.98514724,  0.10318545,  0.8638643
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    floatType beta = 0.67;

    evaluate_at_nodes<3, 8, 1 >(
        std::cbegin( local_point ),              std::cend( local_point ), dt,
        std::cbegin( density_t ),                std::cend( density_t ),
        std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
        std::cbegin( e_t ),                      std::cend( e_t ),
        std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
        std::cbegin( u_t ),                      std::cend( u_t ),
        std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
        std::cbegin( umesh_t ),                  std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
        std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
        std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
        std::cbegin( X ),                        std::cend( X ),
        std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
        std::cbegin( heat_flux ),                std::cend( heat_flux ),
        std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
        std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
        std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
        alpha, beta,
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0 );

    std::array< floatType, 8 * 1 * 8 * 1 > dRdRho;
    std::array< floatType, 8 * 1 * 8 * 1 > dRdE;
    std::array< floatType, 8 * 1 * 8 * 3 > dRdU;
    std::array< floatType, 8 * 1 * 9 > dRdCauchy;
    std::array< floatType, 8 * 1 * 3 > dRdq;
    std::array< floatType, 8 * 1 * 1 > dRdVolumeFraction;
    std::array< floatType, 8 * 1 * 1 > dRdr;
    std::array< floatType, 8 * 1 * 3 > dRdpi;
    std::array< floatType, 8 * 1 * 8 * 3 > dRdUMesh;

    evaluate_at_nodes<3, 8, 1 >(
        std::cbegin( local_point ),              std::cend( local_point ), dt,
        std::cbegin( density_t ),                std::cend( density_t ),
        std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
        std::cbegin( e_t ),                      std::cend( e_t ),
        std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
        std::cbegin( u_t ),                      std::cend( u_t ),
        std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
        std::cbegin( umesh_t ),                  std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
        std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
        std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
        std::cbegin( X ),                        std::cend( X ),
        std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
        std::cbegin( heat_flux ),                std::cend( heat_flux ),
        std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
        std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
        std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
        alpha, beta,
        std::begin( result ),                    std::end( result ),
        std::begin( dRdRho ),                    std::end( dRdRho ),
        std::begin( dRdE ),                      std::end( dRdE ),
        std::begin( dRdU ),                      std::end( dRdU ),
        std::begin( dRdCauchy ),                 std::end( dRdCauchy ),
        std::begin( dRdq ),                      std::end( dRdq ),
        std::begin( dRdVolumeFraction ),         std::end( dRdVolumeFraction ),
        std::begin( dRdr ),                      std::end( dRdr ),
        std::begin( dRdpi ),                     std::end( dRdpi ),
        std::begin( dRdUMesh ),                  std::end( dRdUMesh )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    // Check the derivatives w.r.t. the density
    {

        constexpr unsigned int vardim = 1 * nphases * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( density_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = density_tp1;
            std::array< floatType, vardim > xm = density_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdRho[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the internal energy
    {

        constexpr unsigned int vardim = 1 * nphases * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( e_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = e_tp1;
            std::array< floatType, vardim > xm = e_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdE[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the displacement
    {

        constexpr unsigned int vardim = 3 * nphases * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( u_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = u_tp1;
            std::array< floatType, vardim > xm = u_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdU[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the Cauchy stress
    {

        constexpr unsigned int vardim = 9 * nphases;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( cauchy_stress[ i ] ) + eps;

            std::array< floatType, vardim > xp = cauchy_stress;
            std::array< floatType, vardim > xm = cauchy_stress;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdCauchy[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the volume fractions
    {

        constexpr unsigned int vardim = 1 * nphases;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( volume_fractions[ i ] ) + eps;

            std::array< floatType, vardim > xp = volume_fractions;
            std::array< floatType, vardim > xm = volume_fractions;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdVolumeFraction[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the internal heat generation
    {

        constexpr unsigned int vardim = 1 * nphases;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( internal_heat_generation[ i ] ) + eps;

            std::array< floatType, vardim > xp = internal_heat_generation;
            std::array< floatType, vardim > xm = internal_heat_generation;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdr[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the heat flux
    {

        constexpr unsigned int vardim = 3 * nphases;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( heat_flux[ i ] ) + eps;

            std::array< floatType, vardim > xp = heat_flux;
            std::array< floatType, vardim > xm = heat_flux;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdq[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the net interphase force
    {

        constexpr unsigned int vardim = 3 * nphases;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( net_interphase_force[ i ] ) + eps;

            std::array< floatType, vardim > xp = net_interphase_force;
            std::array< floatType, vardim > xm = net_interphase_force;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( xp ),                       std::cend( xp ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( xm ),                       std::cend( xm ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdpi[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the mesh displacement
    {

        constexpr unsigned int vardim = 3 * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( umesh_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = umesh_tp1;
            std::array< floatType, vardim > xm = umesh_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1 >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdUMesh[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_computeBalanceOfEnergy_multiphase_fea, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test computing the balance of mass in a finite element context
     */

    constexpr unsigned int nphases = 4;

    std::array< floatType, nphases * 8 > density_t = {
        0.20763586, 0.29248941, 0.52001015, 0.90191137, 0.98363088,
        0.25754206, 0.56435904, 0.80696868, 0.39437005, 0.73107304,
        0.16106901, 0.60069857, 0.86586446, 0.98352161, 0.07936579,
        0.42834727, 0.20454286, 0.45063649, 0.54776357, 0.09332671,
        0.29686078, 0.92758424, 0.56900373, 0.457412  , 0.75352599,
        0.74186215, 0.04857903, 0.7086974 , 0.83924335, 0.16593788,
        0.78099794, 0.28653662
    };

    std::array< floatType, nphases * 8 > density_tp1 = {
        0.80793821, 0.00742638, 0.55159273, 0.93193215, 0.58217546,
        0.20609573, 0.71775756, 0.37898585, 0.66838395, 0.02931972,
        0.63590036, 0.03219793, 0.74478066, 0.472913  , 0.12175436,
        0.54263593, 0.06677444, 0.65336487, 0.99608633, 0.76939734,
        0.57377411, 0.10263526, 0.69983407, 0.66116787, 0.04909713,
        0.7922993 , 0.51871659, 0.42586769, 0.78818717, 0.41156922,
        0.48102628, 0.18162884
    };

    std::array< floatType, nphases * 8 > e_t = {
        0.3213189 , 0.845533  , 0.18690375, 0.41729106, 0.98903451,
        0.23659981, 0.91683233, 0.91839747, 0.09129634, 0.46365272,
        0.50221634, 0.31366895, 0.04733954, 0.24168564, 0.09552964,
        0.23824991, 0.80779109, 0.89497829, 0.04322289, 0.30194684,
        0.9805822 , 0.53950482, 0.62630936, 0.00554541, 0.48490944,
        0.98832853, 0.37518553, 0.09703816, 0.46190876, 0.96300447,
        0.34183061, 0.79892273
    };

    std::array< floatType, nphases * 8 > e_tp1 = {
        0.01611863, 0.12695803, 0.77716246, 0.04589523, 0.71099869,
        0.97104614, 0.87168293, 0.71016165, 0.95850974, 0.42981334,
        0.87287891, 0.35595767, 0.92976365, 0.14877766, 0.94002901,
        0.8327162 , 0.84605484, 0.12392301, 0.5964869 , 0.01639248,
        0.72118437, 0.00773751, 0.08482228, 0.22549841, 0.87512453,
        0.36357632, 0.53995994, 0.56810321, 0.22546336, 0.57214677,
        0.6609518 , 0.29824539
    };

    std::array< floatType, nphases * 24 > u_t = {
        -0.16274628, -0.09382215,  0.86470132,  0.17498749,  0.89650474,
         0.11206951,  0.00112284, -0.99293558, -0.03822191,  0.85491   ,
        -0.60326862, -0.89581773, -0.18644221, -0.25520704,  0.71430612,
        -0.94677777,  0.84029846,  0.361806  ,  0.80845199,  0.21505814,
         0.62390662, -0.32891225, -0.30086754, -0.22025154,  0.50959416,
        -0.26141765, -0.51556039,  0.87533671,  0.81602217, -0.30240537,
         0.26927614, -0.45231558, -0.58776974, -0.32732094, -0.34580021,
         0.7645522 ,  0.64460763,  0.41924646,  0.91869045, -0.15491329,
        -0.50993392, -0.76520313, -0.39789328, -0.70947253, -0.81562781,
         0.20586439, -0.2716251 ,  0.12914069, -0.61732856,  0.35381172,
        -0.56898911, -0.44395281,  0.48352084,  0.11947579, -0.33032717,
         0.08597757,  0.38796941,  0.82426424,  0.16142643, -0.53462724,
         0.49339526,  0.55553804, -0.59919737,  0.64114844, -0.07013029,
         0.55953332, -0.52504356, -0.33483946,  0.90739424,  0.31563015,
         0.54575566,  0.37674869, -0.59139176, -0.0586225 ,  0.61792775,
         0.35007025, -0.98794423, -0.82518451, -0.30641056,  0.88873108,
        -0.01761904, -0.45964747, -0.27915256, -0.57869474, -0.15759989,
        -0.56392912,  0.69150501, -0.0874588 , -0.44039596,  0.8657833 ,
        -0.37129729,  0.81942932, -0.91316382,  0.41423012, -0.03222192,
        -0.11155788
    };

    std::array< floatType, nphases * 24 > u_tp1 = {
        -0.00698406, -0.14780869, -0.38870722,  0.83369757,  0.03524692,
         0.60805274,  0.71530357,  0.84476471, -0.39323853, -0.32037829,
         0.19014775, -0.11735173,  0.86568507, -0.2048719 , -0.0444439 ,
         0.23437218, -0.19052103,  0.98495687, -0.80229743, -0.55879336,
        -0.35468974, -0.70455431, -0.43156153,  0.55849059,  0.045784  ,
        -0.93209273,  0.96524517,  0.23201296, -0.88212104,  0.32233754,
        -0.24326126, -0.72865341,  0.12732919,  0.4541599 ,  0.34225321,
        -0.50497369,  0.04973244,  0.07532689,  0.43360673, -0.2802653 ,
         0.59546519,  0.2558437 , -0.92333679,  0.09295804,  0.72382419,
         0.13514833, -0.64834347,  0.02075274,  0.51389167, -0.77978961,
         0.63419816, -0.66503672,  0.06815298, -0.22851304, -0.50275246,
         0.29486503, -0.92521578,  0.52009161,  0.05388128,  0.75154242,
         0.04143664, -0.92993366, -0.71279806,  0.59120918, -0.0160479 ,
        -0.11624146, -0.36313044, -0.4309016 ,  0.93177262, -0.13406134,
         0.76800607,  0.29632625,  0.71685529,  0.70489909,  0.91262406,
         0.39588447,  0.61079387,  0.46625579,  0.21045367,  0.43470827,
         0.43150082, -0.91818441,  0.03222167,  0.58530272, -0.51407563,
        -0.06970403, -0.13002858, -0.19442566, -0.75632094,  0.05142308,
        -0.10750327,  0.32678551,  0.09882612, -0.94491414, -0.93616402,
         0.4027196 
    };

    std::array< floatType, 24 > umesh_t = {
        0.53182759, 0.63440096, 0.84943179, 0.51044152, 0.65634786,
        0.86791477, 0.48312667, 0.6486585 , 0.86600796, 0.50451273,
        0.6267116 , 0.84752498, 0.53695906, 0.68247738, 0.83864355,
        0.515573  , 0.70442428, 0.85712652, 0.48825814, 0.69673492,
        0.85521971, 0.50964421, 0.67478802, 0.83673674
    };

    std::array< floatType, 24 > umesh_tp1 = {
        0.72445532, 0.61102351, 0.72244338, 0.70877313, 0.5669913 ,
        0.69069256, 0.7316781 , 0.55679573, 0.65823773, 0.7473603 ,
        0.60082794, 0.68998856, 0.71831255, 0.63482305, 0.72559852,
        0.70263035, 0.59079084, 0.69384769, 0.72553532, 0.58059527,
        0.66139287, 0.74121752, 0.62462748, 0.6931437 
    };

    std::array< floatType, nphases * 8 > density_dot_t = {
        0.30646975, 0.66526147, 0.11139217, 0.66487245, 0.88785679,
        0.69631127, 0.44032788, 0.43821438, 0.7650961 , 0.565642  ,
        0.08490416, 0.58267109, 0.8148437 , 0.33706638, 0.92757658,
        0.750717  , 0.57406383, 0.75164399, 0.07914896, 0.85938908,
        0.82150411, 0.90987166, 0.1286312 , 0.08178009, 0.13841557,
        0.39937871, 0.42430686, 0.56221838, 0.12224355, 0.2013995 ,
        0.81164435, 0.46798757
    };

    std::array< floatType, nphases * 8 > e_dot_t = {
        0.79884633, 0.2082483 , 0.4433677 , 0.71560128, 0.41051979,
        0.19100696, 0.96749431, 0.65075037, 0.86545985, 0.02524236,
        0.26690581, 0.5020711 , 0.06744864, 0.99303326, 0.2364624 ,
        0.37429218, 0.21401191, 0.10544587, 0.23247979, 0.30061014,
        0.63444227, 0.28123478, 0.36227676, 0.00594284, 0.36571913,
        0.53388598, 0.16201584, 0.59743311, 0.29315247, 0.63205049,
        0.02619661, 0.88759346
    };

    std::array< floatType, nphases * 24 > u_dot_t = {
        -0.92735331, -0.91863362, -0.33449277,  0.89423908,  0.23531995,
        -0.26225032,  0.22395408, -0.58773693, -0.66986711, -0.27636547,
         0.7267067 ,  0.01880345, -0.40619697,  0.90050325,  0.63193218,
        -0.35405211,  0.94419649,  0.9747022 , -0.18267973,  0.31184621,
        -0.1886936 , -0.48530379, -0.83469465, -0.47277931, -0.45704029,
        -0.20272184, -0.63022794,  0.90763681, -0.79424023,  0.25041707,
        -0.11660522, -0.1529639 , -0.25601643,  0.73662942, -0.43904604,
        -0.95884769,  0.83619403,  0.72896056, -0.44619642,  0.0469751 ,
        -0.78182361, -0.81314586,  0.67493222, -0.17946856,  0.32343308,
         0.88640112, -0.50973882, -0.97368034, -0.95170319,  0.41877138,
         0.84910377, -0.06533945, -0.2497817 ,  0.08572085,  0.71783368,
         0.30430775, -0.53404021,  0.54916041, -0.73077301, -0.66888006,
         0.22536457, -0.52243319,  0.4095571 , -0.30096295, -0.44515208,
         0.99783681, -0.91876775,  0.29164504, -0.92260083,  0.52042052,
        -0.53982009, -0.82033627,  0.29689942,  0.46520243,  0.35619063,
        -0.89619811, -0.41138611, -0.09782331, -0.42579342,  0.62102691,
        -0.73776979,  0.22435872,  0.97642989,  0.80511308, -0.55568588,
        -0.99983622,  0.96119468,  0.76542597,  0.83894493, -0.1689929 ,
         0.48923092, -0.574337  , -0.21539186,  0.7030961 , -0.74477555,
         0.78773074
    };

    std::array< floatType, 24 > X = {
        0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0., 0., 1., 1., 0.,
       1., 1., 1., 1., 0., 1., 1.
    };

    std::array< floatType, nphases * 9 > cauchy_stress = {
         0.41516224,  0.91987827,  0.75340936, -0.06388067,  0.25181302,
        -0.08563654, -0.55410753, -0.24664601, -0.79223154,  0.33305424,
        -0.61593971, -0.04906443,  0.93487321, -0.93666214, -0.6965401 ,
        -0.40284163,  0.88361393,  0.81768359, -0.67599832,  0.96223555,
         0.50149505,  0.07995417,  0.86340577,  0.76121428, -0.21736701,
         0.31268639,  0.29477029, -0.34606363, -0.64121965, -0.06638025,
        -0.47343793, -0.28986975,  0.90828794, -0.07772426,  0.36978293,
        -0.32754021
    };

    std::array< floatType, nphases * 3 > heat_flux = {
         0.99172216,  0.31753521, -0.60798107, -0.803632  ,  0.88636114,
         0.88955566,  0.24265675, -0.966017  , -0.54893023,  0.60255357,
         0.75091966, -0.09202037
    };

    std::array< floatType, nphases * 1 > volume_fractions = {
        0.36552062, 0.27422501, 0.11697051, 0.11574454
    };

    std::array< floatType, nphases * 1 > internal_heat_generation = {
         0.9052054 ,  0.61725223, -0.67044128, -0.5858999 
    };

    std::array< floatType, nphases * 3 > net_interphase_force = {
         0.3111031 ,  0.52932843,  0.62062969, -0.67332462,  0.96825658,
        -0.54439587,  0.17883087,  0.17523151,  0.93472377,  0.31533489,
         0.16980853,  0.03754516
    };

    std::array< floatType, nphases * 8 > answer = {
        8.48561396e-02, -7.84057240e-02, -2.04859969e-02,  8.42953194e-02,
       -1.85890771e-01,  1.30399296e-01, -1.89514818e-01, -1.42221259e-03,
       -6.15974396e-02, -7.19744928e-02,  7.67919795e-02, -8.22857729e-02,
       -2.10349657e-03, -3.62151480e-02,  2.81241169e-02, -1.52729966e-02,
        1.46841956e-02, -2.30609745e-02,  1.46502219e-02,  1.40228979e-03,
        4.67431720e-02, -7.41206505e-02,  4.73426858e-02,  3.97909230e-03,
        5.46259217e-03, -8.97678887e-03,  5.82472207e-03,  2.99928234e-04,
        1.71658140e-03, -2.79231155e-03,  1.80211495e-03,  1.13383339e-04
    };

    std::array< floatType, nphases * 8 > result;

    std::array< floatType, 3 > local_point = {
        0.52931508, -0.78788948, -0.9958162
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    floatType beta = 0.67;

    evaluate_at_nodes<3, 8, nphases >(
        std::cbegin( local_point ),              std::cend( local_point ), dt,
        std::cbegin( density_t ),                std::cend( density_t ),
        std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
        std::cbegin( e_t ),                      std::cend( e_t ),
        std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
        std::cbegin( u_t ),                      std::cend( u_t ),
        std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
        std::cbegin( umesh_t ),                  std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
        std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
        std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
        std::cbegin( X ),                        std::cend( X ),
        std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
        std::cbegin( heat_flux ),                std::cend( heat_flux ),
        std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
        std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
        std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
        alpha, beta,
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0 );

    std::array< floatType, 8 * 1 * nphases * 8 * 1 * nphases > dRdRho;
    std::array< floatType, 8 * 1 * nphases * 8 * 1 * nphases > dRdE;
    std::array< floatType, 8 * 1 * nphases * 8 * 3 * nphases > dRdU;
    std::array< floatType, 8 * 1 * nphases * 9 * nphases > dRdCauchy;
    std::array< floatType, 8 * 1 * nphases * 3 * nphases > dRdq;
    std::array< floatType, 8 * 1 * nphases * 1 * nphases > dRdVolumeFraction;
    std::array< floatType, 8 * 1 * nphases * 1 * nphases > dRdr;
    std::array< floatType, 8 * 1 * nphases * 3 * nphases > dRdpi;
    std::array< floatType, 8 * 1 * nphases * 8 * 3 > dRdUMesh;

    evaluate_at_nodes<3, 8, nphases >(
        std::cbegin( local_point ),              std::cend( local_point ), dt,
        std::cbegin( density_t ),                std::cend( density_t ),
        std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
        std::cbegin( e_t ),                      std::cend( e_t ),
        std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
        std::cbegin( u_t ),                      std::cend( u_t ),
        std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
        std::cbegin( umesh_t ),                  std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
        std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
        std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
        std::cbegin( X ),                        std::cend( X ),
        std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
        std::cbegin( heat_flux ),                std::cend( heat_flux ),
        std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
        std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
        std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
        alpha, beta,
        std::begin( result ),                    std::end( result ),
        std::begin( dRdRho ),                    std::end( dRdRho ),
        std::begin( dRdE ),                      std::end( dRdE ),
        std::begin( dRdU ),                      std::end( dRdU ),
        std::begin( dRdCauchy ),                 std::end( dRdCauchy ),
        std::begin( dRdq ),                      std::end( dRdq ),
        std::begin( dRdVolumeFraction ),         std::end( dRdVolumeFraction ),
        std::begin( dRdr ),                      std::end( dRdr ),
        std::begin( dRdpi ),                     std::end( dRdpi ),
        std::begin( dRdUMesh ),                  std::end( dRdUMesh )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    floatType eps = 1e-4;

    // Check the derivatives w.r.t. the density
    {

        constexpr unsigned int vardim = 1 * nphases * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( density_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = density_tp1;
            std::array< floatType, vardim > xm = density_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdRho[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the internal energy
    {

        constexpr unsigned int vardim = 1 * nphases * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( e_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = e_tp1;
            std::array< floatType, vardim > xm = e_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdE[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the displacement
    {

        constexpr unsigned int vardim = 3 * nphases * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( u_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = u_tp1;
            std::array< floatType, vardim > xm = u_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdU[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the Cauchy stress
    {

        constexpr unsigned int vardim = 9 * nphases;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( cauchy_stress[ i ] ) + eps;

            std::array< floatType, vardim > xp = cauchy_stress;
            std::array< floatType, vardim > xm = cauchy_stress;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdCauchy[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the volume fractions
    {

        constexpr unsigned int vardim = 1 * nphases;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( volume_fractions[ i ] ) + eps;

            std::array< floatType, vardim > xp = volume_fractions;
            std::array< floatType, vardim > xm = volume_fractions;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdVolumeFraction[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the internal heat generation
    {

        constexpr unsigned int vardim = 1 * nphases;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( internal_heat_generation[ i ] ) + eps;

            std::array< floatType, vardim > xp = internal_heat_generation;
            std::array< floatType, vardim > xm = internal_heat_generation;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdr[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the heat flux
    {

        constexpr unsigned int vardim = 3 * nphases;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( heat_flux[ i ] ) + eps;

            std::array< floatType, vardim > xp = heat_flux;
            std::array< floatType, vardim > xm = heat_flux;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdq[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the net interphase force
    {

        constexpr unsigned int vardim = 3 * nphases;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( net_interphase_force[ i ] ) + eps;

            std::array< floatType, vardim > xp = net_interphase_force;
            std::array< floatType, vardim > xm = net_interphase_force;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( xp ),                       std::cend( xp ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),                std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( xm ),                       std::cend( xm ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdpi[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the mesh displacement
    {

        constexpr unsigned int vardim = 3 * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( umesh_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = umesh_tp1;
            std::array< floatType, vardim > xm = umesh_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( xp ),                       std::cend( xp ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases >(
                std::cbegin( local_point ),              std::cend( local_point ), dt,
                std::cbegin( density_t ),                std::cend( density_t ),
                std::cbegin( density_tp1 ),              std::cend( density_tp1 ),
                std::cbegin( e_t ),                      std::cend( e_t ),
                std::cbegin( e_tp1 ),                    std::cend( e_tp1 ),
                std::cbegin( u_t ),                      std::cend( u_t ),
                std::cbegin( u_tp1 ),                    std::cend( u_tp1 ),
                std::cbegin( umesh_t ),                  std::cend( umesh_t ),
                std::cbegin( xm ),                       std::cend( xm ),
                std::cbegin( density_dot_t ),            std::cend( density_dot_t ),
                std::cbegin( e_dot_t ),                  std::cend( e_dot_t ),
                std::cbegin( u_dot_t ),                  std::cend( u_dot_t ),
                std::cbegin( X ),                        std::cend( X ),
                std::cbegin( cauchy_stress ),            std::cend( cauchy_stress ),
                std::cbegin( heat_flux ),                std::cend( heat_flux ),
                std::cbegin( volume_fractions ),         std::cend( volume_fractions ),
                std::cbegin( internal_heat_generation ), std::cend( internal_heat_generation ),
                std::cbegin( net_interphase_force ),     std::cend( net_interphase_force ),
                alpha, beta,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dRdUMesh[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}
