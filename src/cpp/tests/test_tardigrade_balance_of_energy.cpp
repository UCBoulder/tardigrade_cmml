/**
  * \file test_tardigrade_balance_equations_balance_of_energy.cpp
  *
  * Tests for tardigrade_balance_equations_balance_of_energy
  */

#include<tardigrade_balance_of_energy.h>
#include<tardigrade_finite_element_utilities.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_hydraLinearTestMaterial.h>
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

std::vector< double > linear_test_params =
{
    0.37, 0.45, 0.22, 0.99, 0.32, 0.16, 0.12, 0.87, 0.25, 0.32, 0.30, 0.75, 0.88, 0.26, 0.14, 0.72, 0.47, 0.09, 0.42, 0.07, 0.98,
    0.91, 0.95, 0.86, 0.02, 0.13, 0.67, 0.14, 0.46, 0.20, 0.15, 0.11, 0.84, 0.84, 0.49, 0.61, 0.23, 0.92, 0.57, 0.64, 0.45, 0.67,
    0.65, 0.97, 0.60, 0.17, 0.72, 0.30, 0.23, 0.46, 0.54, 0.04, 0.59, 0.67, 0.91, 0.54, 0.10, 0.58, 0.86, 0.67, 0.85, 0.14, 0.50,
    0.55, 0.84, 0.33, 0.72, 0.60, 0.31, 0.58, 0.83, 0.48, 0.01, 0.22, 0.41, 0.68, 0.92, 0.36, 0.25, 0.10, 0.73, 0.71, 0.64, 0.81,
    0.45, 0.12, 0.33, 0.17, 0.94, 0.90, 0.63, 0.94, 0.95, 0.91, 0.30, 0.03, 0.52, 0.64, 0.47, 0.98, 0.70, 0.59, 0.78, 0.53, 0.87,
    0.47, 0.10, 0.20, 0.45, 0.44, 0.93, 0.43, 0.08, 0.48, 0.71, 0.23, 0.87, 0.01, 0.70, 0.25, 0.73, 0.61, 0.90, 0.09, 0.49, 0.61,
    0.86, 0.02, 0.79, 0.48, 0.90, 0.10, 0.23, 0.49, 0.69, 0.13, 0.61, 0.44, 0.35, 0.71, 0.57, 0.09, 0.39, 0.04, 0.38, 0.92, 0.63,
    0.71, 0.29, 0.68, 0.96, 0.74, 0.86, 0.15, 0.50, 0.20, 0.59, 0.63, 0.06, 0.98, 0.75, 0.65, 0.33, 0.37, 0.07, 0.12, 0.37, 0.77,
    0.73, 0.19, 0.28, 0.06, 0.72, 0.93, 0.43, 0.27, 0.23, 0.78, 0.65, 0.99, 0.95, 0.07, 0.18, 0.15, 0.50, 0.40, 0.79, 0.20, 0.10,
    0.63, 0.18, 0.81, 0.26, 0.39, 0.35, 0.95, 0.26, 0.42, 0.41, 0.85, 0.82, 0.85, 0.62, 0.25, 0.20, 0.12, 0.57, 0.39, 0.06, 0.71,
    0.46, 0.06, 0.44, 0.95, 0.18, 0.64, 0.43, 0.33, 0.15, 0.38, 0.50, 0.04, 0.32, 0.49, 0.57, 0.10, 0.65, 0.34, 0.18, 0.81, 0.07,
    0.93, 0.71, 0.48, 0.01, 0.39, 0.65, 0.86, 0.62, 0.40, 0.45, 0.87, 0.35, 0.07, 0.44, 1.00, 0.38, 0.19, 0.05, 0.17, 0.32, 0.57,
    0.67, 0.04, 0.94, 0.24, 0.15, 0.53, 0.68, 0.31, 0.67, 0.29, 0.90, 0.88, 0.07, 0.15, 0.70, 0.00, 0.82, 0.36, 0.74, 0.04, 0.21,
    0.07, 0.18, 0.38, 0.49, 0.92, 0.63, 0.71, 0.08, 0.29, 0.99, 0.56, 0.79, 0.03, 0.08, 0.11, 0.02, 0.92, 0.25, 1.00, 0.74, 0.25,
    0.99, 0.87, 0.16, 0.19, 0.68, 0.08, 0.58, 0.07, 0.22, 0.95, 0.50, 0.57, 0.27, 0.44, 0.13, 0.52, 0.05, 0.05, 0.68, 0.05, 0.52,
    0.27, 0.19, 0.57, 0.30, 0.46, 0.14, 0.13, 0.11, 0.07, 0.55, 0.76, 0.79, 0.96, 0.17, 0.14, 0.47, 0.86, 0.77, 0.56, 0.16, 0.39,
    0.72, 0.37, 0.04, 0.58, 0.72, 0.71, 0.16, 0.34, 0.41, 0.42, 0.11, 0.39, 0.51, 0.71, 0.57, 0.40, 0.93, 0.60, 0.97, 0.87, 0.67,
    0.90, 0.99, 0.22, 0.69, 0.92, 0.02, 0.80, 0.09, 0.92, 0.61, 0.16, 0.04, 0.83, 0.14, 0.94, 0.09, 0.84, 0.42, 0.87, 0.77, 0.85,
    0.47, 0.71, 0.32, 0.05, 0.23, 0.00, 0.36, 0.49, 0.02, 0.01, 0.69, 0.10, 0.07, 0.50, 0.72, 0.98, 0.30, 0.12, 0.07, 0.21, 0.10,
    0.00, 0.29, 0.16, 0.02, 0.35, 0.83, 0.30, 0.71, 0.73, 0.03, 0.96, 0.14, 0.72, 0.38, 0.87, 0.99, 0.19, 0.82, 0.51, 0.97, 0.66,
    0.59, 0.23, 0.72, 0.29, 0.85, 0.98, 0.01, 0.95, 0.36, 0.49, 0.92, 0.78, 0.24, 0.10, 0.84, 0.75, 0.86, 0.38, 0.49, 0.47, 0.31,
    0.34, 0.72, 0.75, 0.86, 0.84, 0.75, 0.33, 0.71, 0.97, 0.92, 0.97, 0.04, 0.61, 0.22, 0.25, 0.84, 0.75, 0.33, 0.52, 0.19, 0.00,
    0.90, 0.71, 0.53, 0.05, 0.42, 0.38, 0.17, 0.16, 0.17, 0.60, 0.12, 0.24, 0.95, 0.22, 0.76, 0.60, 0.41, 0.60, 0.13, 0.48, 0.97,
    0.43, 0.24, 0.93, 0.61, 0.73, 0.54, 0.01, 0.76, 0.98, 0.72, 0.81, 0.85, 1.00, 0.29, 0.45, 0.54, 0.15, 0.48, 0.28, 0.84, 0.35,
    0.08, 0.01, 0.86, 0.86, 0.93, 0.07, 0.39, 0.43, 0.86, 0.26, 0.30, 0.73, 0.18, 0.64, 0.74, 0.18, 0.22, 0.70, 0.79, 0.95, 0.31,
    0.55, 0.15, 0.44, 0.26, 0.35, 0.80, 0.80, 0.70, 0.12, 0.28, 0.94, 0.45, 0.62, 0.02, 0.52, 0.69, 0.06, 0.12, 0.37, 0.94, 0.25,
    0.02, 0.81, 0.85, 0.54, 0.89, 0.18, 0.27, 0.02, 0.58, 0.89, 0.61, 0.08, 0.69, 0.68, 0.17, 0.83, 0.91, 0.39, 0.31, 0.57, 0.97,
    0.72, 0.79, 0.77, 0.96, 0.82, 0.27, 0.23, 0.24, 0.21, 0.21, 0.93, 0.95, 0.09, 0.06, 0.86, 0.56, 0.21, 0.46, 0.07, 0.98, 0.14,
    0.81, 0.96, 0.71, 0.09, 0.59, 0.13, 0.87, 0.67, 0.94, 0.15, 0.45, 0.54, 0.11, 0.68, 0.57, 0.35, 0.41, 0.80, 0.36, 0.97, 0.92,
    0.08, 0.85, 0.67, 0.06, 0.96, 0.83, 0.83, 0.92, 0.92, 0.22, 0.94, 0.05, 0.88, 0.46, 0.83, 0.27, 0.75, 0.80, 0.08, 0.91, 0.82,
    0.60, 0.71, 0.35, 0.47, 0.15, 0.50, 0.78, 0.56, 0.13, 0.83, 0.74, 0.04, 0.80, 0.87, 0.51, 0.38, 0.08, 0.83, 0.00, 0.34, 0.33,
    0.91, 0.45, 0.89, 0.09, 0.85, 0.61, 0.54, 0.31, 0.49, 0.27, 0.46, 0.77, 0.63, 0.68, 0.89, 0.05, 0.20, 0.72, 0.41, 0.65, 0.90,
    0.74, 0.93, 0.15, 0.76, 0.70, 0.01, 0.04, 0.03, 0.53, 0.78, 0.41, 0.00, 0.17, 0.95, 0.24, 0.85, 0.37, 0.52, 0.49, 0.13, 0.37,
    1.00, 0.56, 0.47, 0.55, 0.26, 0.57, 0.73, 0.11, 0.28, 0.21, 0.74, 0.94, 0.42, 0.63, 0.03, 0.13, 0.38, 0.66, 0.64, 0.19, 0.13,
    0.94, 0.36, 0.46, 0.10, 0.91, 0.67, 0.51, 0.43, 0.05, 0.21, 0.86, 0.62, 0.19, 0.19, 0.35, 0.62, 0.54, 0.53, 0.87, 0.86, 0.85,
    0.70, 0.39, 0.27, 0.96, 0.05, 0.62, 0.85, 0.27, 0.98, 0.89, 0.20, 0.27, 0.90, 0.11, 0.19, 0.61, 0.15, 0.57, 0.58, 0.60, 0.25,
    0.82, 0.97, 0.28, 0.87, 0.29, 0.35, 0.79, 0.14, 0.06, 0.39, 0.13, 0.08, 0.59, 0.18, 0.30, 0.73, 0.09, 0.99, 0.31, 0.55, 0.79,
    0.24, 0.34, 0.21, 0.38, 0.73, 0.71, 0.18, 0.48, 0.69, 0.74, 0.64, 0.19, 0.87, 0.29, 0.37, 0.56, 0.86, 0.65, 0.70, 0.22, 0.35,
    0.14, 0.74, 0.95, 0.73, 0.12, 0.69, 0.43, 0.09, 1.00, 0.72, 0.62, 0.87, 0.33, 0.67, 0.47, 0.12, 0.54, 0.63, 0.47, 0.13, 0.08,
    0.66, 0.84, 0.92, 0.39, 0.94, 0.82, 0.42, 0.41, 0.33, 0.07, 0.01, 0.25, 0.25, 0.73, 0.16, 0.69, 0.73, 0.46, 0.12, 0.72, 0.91,
    0.37, 0.76, 0.74, 0.05, 0.06, 0.68, 0.83, 0.92, 0.14, 0.58, 0.28, 0.22, 0.74, 0.41, 0.44, 0.84, 0.24, 0.53, 0.83, 0.13, 0.60,
    0.84, 0.76, 0.88, 0.03, 0.69, 0.20, 0.56, 0.77, 0.42, 0.70, 0.84, 0.95, 0.55, 0.32, 0.39, 0.87, 0.24, 0.88, 0.06, 0.97, 0.12,
    0.88, 0.54, 0.56, 0.98, 0.95, 0.80, 0.65, 0.68, 0.13, 0.47, 0.43, 0.10, 0.01, 0.02, 0.74, 0.45, 0.49, 0.25, 0.54, 0.13, 0.16,
    0.82, 0.03, 0.92, 0.34, 0.93, 0.06, 0.31, 0.59, 0.94, 0.09, 0.81, 0.61, 0.64, 0.99, 0.32, 0.36, 0.54, 0.63, 0.31, 0.71, 0.08,
    0.19, 0.01, 0.75, 0.95, 0.77, 0.17, 0.28, 0.67, 0.31, 0.41, 0.26, 0.53, 0.87, 0.24, 0.53, 0.67, 0.01, 0.25, 0.57, 0.83, 0.23,
    0.90, 0.82, 0.39, 0.81, 0.85, 0.28, 0.12, 0.10, 0.88, 0.49, 0.67, 0.39, 0.78, 0.82, 0.89, 0.64, 0.89, 0.87, 0.27, 0.59, 0.26,
    0.64, 0.90, 0.48, 0.82, 0.02, 0.87, 0.94, 0.31, 0.71, 0.82, 0.12, 0.41, 0.21, 0.75, 0.41, 0.00, 0.26, 0.20, 0.07, 0.53, 0.34,
    0.22, 0.94, 0.32, 0.13, 0.52, 0.89, 0.78, 0.15, 0.30, 0.24, 0.93, 0.28, 0.17, 0.46, 0.61, 0.16, 0.02, 0.09, 0.65, 0.16, 0.26,
    0.68, 0.87, 0.95, 0.38, 0.23, 0.74, 0.71, 0.97, 0.10, 0.88, 0.18, 0.75, 1.00, 0.23, 0.74, 0.77, 0.28, 0.92, 0.07, 0.23, 0.29,
    0.09, 0.25, 0.79, 0.27, 0.51, 0.14, 0.59, 0.61, 0.22, 0.07, 0.70, 0.63, 0.48, 0.90, 0.98, 0.73, 0.96, 0.15, 0.16, 1.00, 0.38,
    0.56, 0.94, 0.60, 0.87, 0.71, 0.86, 0.47, 0.95, 0.13, 0.81, 0.24, 0.62, 0.78, 0.70, 0.44, 0.93, 0.16, 0.72, 0.25, 0.69, 0.19,
    0.05, 0.78, 0.82, 0.16, 0.81, 0.33, 0.57, 0.60, 0.91, 0.57, 0.87, 0.42, 0.65, 0.19, 0.71, 0.93, 0.43, 0.65, 0.54, 0.01, 0.91,
    0.81, 0.42, 0.09, 0.35, 0.53, 0.24, 0.83, 0.30, 0.20, 0.17, 0.73, 0.55, 0.70, 0.29, 0.40, 0.90, 0.33, 0.90, 0.97, 0.02, 0.49,
    0.62, 0.01, 0.35, 0.03, 0.62, 0.46, 0.26, 0.78, 0.54, 0.33, 0.94, 0.98, 0.18, 0.98, 0.59, 0.88, 0.29, 0.05, 0.01, 0.43, 0.00,
    0.45, 0.74, 0.67, 0.73, 0.59, 0.95, 0.24, 0.90, 0.98, 0.00, 0.43, 0.68, 0.61, 0.99, 0.66, 0.62, 0.93, 0.90, 0.83, 0.41, 0.70,
    0.96, 0.24, 0.49, 0.76, 0.60, 0.78, 0.26, 0.13, 0.47, 0.66, 0.60, 0.61, 0.09, 0.49, 0.20, 0.71, 0.56, 0.60, 0.60, 0.31, 0.08,
    0.41, 0.30, 0.57, 0.22, 0.22, 0.33, 0.17, 0.40, 0.82, 0.85, 0.11, 0.33, 0.61, 0.23, 0.92, 0.71, 0.22, 0.11, 0.07, 0.19, 0.40,
    0.75, 0.14, 0.04, 0.89, 0.34, 0.33, 0.53, 0.42, 0.30, 0.72, 0.79, 0.83, 0.77, 0.69, 0.19, 0.73, 0.27, 0.72, 0.10, 0.42, 0.39,
    0.04, 0.88, 0.68, 0.78, 0.67, 0.34, 0.51, 0.13, 0.16, 0.32, 0.64, 0.60, 0.78, 0.53, 0.68, 0.04, 0.93, 0.85, 0.94, 0.17, 0.15,
    0.54, 0.84, 0.74, 0.26, 0.48, 0.24, 0.56, 0.57, 0.08, 0.63, 0.10, 0.03, 0.86, 0.25, 0.47, 0.94, 0.63, 0.37, 0.85, 0.21, 0.49,
    0.57, 0.59, 0.02, 0.38, 0.51, 0.46, 0.82, 0.72, 0.93, 0.22, 0.16, 0.20, 0.66, 0.25, 0.66, 0.68, 0.08, 0.74, 0.31, 0.22, 0.26,
    0.90, 0.61, 0.65, 0.22, 0.36, 0.81, 0.34, 0.81, 0.20, 0.73, 0.12, 0.82, 0.67, 0.96, 0.03, 0.13, 0.03, 0.02, 0.85, 0.63, 0.64,
    0.38, 0.05, 0.43, 0.16, 0.93, 0.88, 0.30, 0.02, 0.02, 0.53, 0.82, 0.34, 0.66, 0.86, 0.64, 0.01, 0.55, 0.17, 0.80, 0.41, 0.44,
    0.30, 0.23, 0.16, 0.41, 0.54, 0.90, 0.40, 0.32, 0.82, 0.70, 0.11, 0.71, 0.95, 0.50, 0.49, 0.72, 0.20, 0.37, 0.32, 1.00, 0.16,
    0.17, 0.46, 0.27, 0.87, 0.44, 0.68, 0.82, 0.27, 0.73, 0.40, 0.25, 0.89, 0.57, 0.50, 0.18, 0.90, 0.64, 0.95, 0.04, 0.12, 0.88,
    0.45, 0.47, 0.22, 0.49, 0.44, 0.72, 0.64, 0.80, 0.78, 0.42, 0.30, 0.33, 0.99, 0.41, 0.76, 0.08, 0.43, 0.33, 0.16, 0.02, 0.25,
    0.79, 0.30, 0.87, 0.58, 0.80, 1.00, 0.64, 0.34, 0.31, 0.97, 0.71, 0.50, 0.21, 0.49, 0.76, 0.97, 0.67, 0.02, 0.47, 0.83, 0.51,
    0.34, 0.86, 0.58, 0.50, 0.49, 0.68, 0.46, 0.57, 0.29, 0.17, 0.82, 0.32, 0.13, 0.18, 0.41, 0.29, 0.70, 0.35, 0.59, 0.17, 0.14,
    0.91, 0.86, 0.90, 0.49, 0.07, 0.31, 0.88, 0.42, 0.33, 0.30, 0.55, 0.45, 0.37, 0.42, 0.27, 0.19, 0.20, 0.04, 0.95, 0.96, 0.10,
    0.07, 0.42, 0.76, 0.83, 0.02, 0.28, 0.54, 0.14, 0.97, 0.25, 0.53, 0.11, 0.34, 0.12, 0.91, 0.95, 0.54, 0.35, 0.32, 0.92, 0.62,
    0.41, 0.77, 0.54, 0.22, 0.24, 0.28, 0.02, 0.26, 0.61, 0.75, 0.11, 0.97, 0.49, 0.06, 0.94, 0.27, 0.31, 0.46, 0.32, 0.18, 0.13,
    0.51, 0.07, 0.71, 0.34, 0.99, 0.06, 0.05, 0.25, 0.35, 0.81, 0.83, 0.72, 0.83, 0.73, 0.87, 0.46, 0.28, 0.63, 0.36, 0.07, 0.86,
    0.83, 0.50, 0.47, 0.12, 0.13, 0.80, 0.85, 0.01, 0.78, 0.54, 0.44, 0.48, 0.71, 0.00, 0.06, 0.84, 1.00, 0.24, 0.40, 0.70, 0.66,
    0.78, 0.91, 0.26, 0.53, 0.14, 0.97, 0.76, 0.19, 0.74, 0.80, 0.14, 0.50, 0.82, 0.35, 0.05, 0.34, 0.78, 0.10, 0.31, 0.85, 0.16,
    0.39, 0.55, 0.18, 0.70, 0.78, 0.63, 0.86, 0.13, 0.00, 0.27, 0.84, 0.13, 0.76, 0.93, 0.56, 0.88, 0.25, 0.72, 0.78, 0.49, 0.60,
    0.09, 0.28, 0.08, 0.24, 0.89, 0.03, 0.50, 0.52, 0.29, 0.78, 0.59, 0.20, 0.47, 0.02, 0.02, 0.73, 0.17, 0.97, 0.36, 0.70, 0.46,
    0.41, 0.67, 0.63, 0.88, 0.08, 0.51, 0.61, 0.65, 0.63, 0.37, 0.03, 0.71, 0.04, 0.47, 0.72, 0.27, 0.03, 0.09, 0.14, 0.67, 0.65,
    0.98, 0.32, 0.13, 0.17, 0.03, 0.82, 0.08, 0.89, 0.05, 0.81, 0.20, 0.02, 0.70, 0.00, 0.30, 0.06, 0.12, 0.13, 0.72, 0.97, 0.24,
    0.40, 0.37, 0.26, 0.71, 0.27, 0.51, 0.13, 0.85, 0.69, 0.54, 0.61, 0.22, 0.84, 0.12, 0.88, 0.61, 0.94, 0.29, 0.94, 0.64, 0.01,
    0.10, 0.59, 0.10, 0.75, 0.67, 0.19, 0.64, 0.57, 0.25, 0.71, 0.25, 0.86, 0.40, 0.50, 0.66, 0.88, 0.46, 0.65, 0.60, 0.07, 0.63,
    0.05, 0.51, 0.11, 0.61, 0.72, 0.79, 0.59, 0.99, 0.36, 0.01, 0.82, 0.45, 0.54, 0.85, 0.26, 0.33, 0.10, 0.68, 0.21, 0.85, 0.61,
    0.91, 0.06, 0.99, 0.40, 0.55, 0.64, 0.49, 0.53, 0.42, 0.19, 0.97, 0.80, 0.60, 0.06, 0.89, 0.08, 0.78, 0.06, 0.76, 0.84, 0.45,
    0.55, 0.84, 0.97, 0.02, 0.13, 0.86, 0.68, 0.27, 0.45, 0.11, 0.03, 0.73, 0.72, 0.69, 0.88, 0.29, 0.42, 0.25, 0.63, 0.64, 0.41,
    0.94, 0.49, 0.92, 0.73, 0.94, 0.76, 0.22, 0.69, 0.24, 0.29, 0.22, 0.50, 0.43, 0.91, 0.59, 0.36, 0.62, 0.34, 0.84, 0.09, 0.22,
    0.09, 0.87, 0.27, 0.99, 0.15, 0.32, 0.89, 0.64, 0.49, 0.16, 0.79, 0.33, 0.75, 0.58, 0.62, 0.33, 0.08, 0.60, 0.74, 0.58, 0.60,
    0.70, 0.04, 0.22, 0.61, 0.14, 0.27, 0.05, 0.22, 0.01, 0.82, 0.59, 0.55, 0.01, 0.85, 0.00, 0.79, 0.22, 0.82, 0.13, 0.55, 0.03,
    0.55, 0.13, 0.81, 0.24, 0.12, 0.03, 0.57, 0.66, 0.89, 0.06, 0.39, 0.16, 0.50, 0.24, 0.41, 0.37, 0.69, 0.74, 0.42, 0.24, 0.79,
    0.85, 0.54, 0.25, 0.78, 0.36, 0.92, 0.27, 0.94, 0.43, 0.75, 0.21, 0.62, 0.58, 0.42, 0.74, 0.82, 0.73, 0.81, 0.75, 0.15, 0.77,
    0.19, 0.10, 0.90, 0.60, 0.18, 0.64, 0.33, 0.14, 0.80, 0.11, 0.96, 0.77, 0.71, 0.24, 0.76, 0.02, 0.55, 0.31, 0.83, 0.76, 0.49,
    0.47, 0.07, 0.23, 0.71, 0.07, 0.79, 0.98, 0.08, 0.46, 0.58, 0.81, 0.15, 0.36, 0.84, 0.52, 0.16, 0.87, 0.94, 0.74, 0.97, 0.93,
    0.12, 0.41, 0.60, 0.73, 0.54, 0.65, 0.12, 0.41, 0.01, 0.65, 0.16, 0.40, 0.50, 0.67, 0.76, 0.72, 0.57, 0.53, 0.86, 0.87, 0.06,
    0.27, 0.92, 0.59, 0.13, 0.07, 0.09, 0.77, 0.22, 0.24, 0.13, 0.60, 0.51, 0.28, 0.66, 0.78, 0.98, 0.48, 0.41, 0.29, 0.78, 0.97,
    0.60, 0.64, 0.69, 0.42, 0.22, 0.26, 0.62, 0.64, 0.46, 0.04, 0.88, 0.24, 0.58, 0.74, 0.18, 0.47, 0.00, 0.05, 0.25, 0.66, 0.89,
    0.42, 0.03, 0.29, 0.96, 0.11, 0.76, 0.75, 0.59, 0.93, 0.64, 0.97, 0.83, 0.79, 0.95, 0.62, 0.32, 0.10, 0.69, 0.52, 0.39, 0.95,
    0.91, 0.38, 0.95, 0.43, 0.21, 0.62, 0.09, 0.15, 0.81, 0.82, 0.24, 0.75, 0.97, 0.96, 0.44, 0.86, 0.56, 0.56, 0.24, 0.38, 0.21,
    0.74, 0.64, 0.80, 0.62, 0.28, 0.95, 0.65, 0.48, 0.78, 0.25, 0.81, 0.34, 0.07, 0.78, 0.19, 0.64, 0.27, 0.69, 0.13, 0.62, 0.19,
    0.50, 0.87, 0.48, 0.77, 0.15, 0.66, 0.15, 0.77, 0.71, 0.40, 0.47, 0.17, 0.64, 0.91, 0.56, 0.32, 0.41, 0.22, 0.44, 0.50, 0.94,
    0.43, 0.86, 0.17, 0.55, 0.59, 0.59, 0.50, 0.77, 0.17, 0.07, 0.14, 0.56, 0.25, 0.74, 0.06, 0.91, 0.95, 0.79, 0.88, 0.38, 0.16,
    0.20, 0.50, 0.66, 0.57, 0.10, 0.73, 0.79, 0.03, 0.46, 0.44, 0.07, 0.18, 0.81, 0.58, 0.22, 0.46, 0.85, 0.10, 0.96, 0.14, 0.29,
    0.18, 0.47, 0.31, 0.99, 0.02, 0.32, 0.59, 0.34, 0.46, 0.43, 0.09, 0.65, 0.83, 0.53, 0.15, 0.48, 0.01, 0.91, 0.24, 0.90, 0.58,
    0.98, 0.85, 0.01, 0.02, 0.74, 0.52, 0.76, 0.59, 0.23, 0.79, 0.03, 0.81, 0.67, 0.67, 0.04, 0.21, 0.49, 0.44, 0.03, 0.34, 0.97,
    0.31, 0.63, 0.68, 0.62, 0.29, 0.57, 0.70, 0.38, 0.28, 0.27, 0.92, 0.64, 0.75, 0.57, 0.63, 0.28, 0.89, 0.50, 0.81, 0.52, 0.12,
    0.11, 0.46, 0.43, 0.91, 0.13, 0.69, 0.92, 0.41, 0.88, 0.14, 0.11, 0.23, 0.52, 0.65, 0.77, 0.77, 0.11, 0.87, 0.11, 0.90, 0.96,
    0.58, 0.69, 0.37, 0.55, 0.83, 0.16, 0.50, 0.41, 0.87, 0.13, 0.89, 0.25, 0.61, 0.38, 0.93, 0.20, 0.50, 0.23, 0.80, 0.27, 0.95,
    0.92, 0.48, 0.01, 0.97, 0.69, 0.49, 0.60, 0.11, 0.16, 0.64, 0.05, 0.70, 0.92, 0.71, 0.56, 0.79, 0.75, 0.59, 0.20, 0.01, 0.46,
    0.24, 0.23, 0.51, 0.38, 0.50, 0.27, 0.61, 0.62, 0.21, 0.57, 0.45, 0.97, 0.72, 0.18, 0.88, 0.48, 0.85, 0.83, 0.02, 0.23, 0.91,
    0.77, 0.14, 0.21, 0.62, 0.28, 0.20, 0.98, 0.46, 0.52, 0.65, 0.08, 0.40, 0.91, 0.94, 0.50, 0.89, 0.50, 0.99, 0.79, 0.83, 0.19,
    0.00, 0.50, 0.98, 0.97, 0.61, 0.18, 0.01, 0.50, 0.30, 0.34, 0.83, 0.09, 0.09, 0.44, 0.41, 0.06, 0.55, 0.43, 0.99, 0.07, 0.12,
    0.16, 0.32, 0.66, 0.97, 0.91, 0.44, 0.07, 0.64, 0.42, 0.02, 0.82, 0.53, 0.44, 0.72, 0.86, 0.93, 0.33, 0.18, 0.33, 0.98, 0.10,
    0.95, 0.62, 0.01, 0.86, 0.68, 0.14, 0.34, 0.87, 0.73, 0.32, 0.79, 0.25, 0.55, 0.93, 0.44, 0.56, 0.78, 0.20, 0.45, 0.92, 0.65,
    0.70, 0.96, 0.59, 0.87, 0.99, 0.68, 0.68, 0.18, 0.46, 0.74, 0.34, 0.12, 0.02, 0.67, 0.03, 0.49, 0.88, 0.82, 0.72, 0.86, 0.75,
    0.02, 0.47, 0.01, 0.56, 0.24, 0.44, 0.20, 0.82, 0.18, 0.33, 0.02, 0.31, 0.03, 0.45, 0.03, 0.39, 0.61, 0.91, 0.55, 0.52, 0.05,
    0.37, 0.39, 0.32, 0.89, 0.40, 0.21, 0.79, 0.81, 0.35, 0.82, 0.71, 0.11, 0.37, 0.85, 0.37, 0.87, 0.45, 0.56, 0.39, 0.61, 0.52,
    0.84, 0.97, 0.29, 0.15, 0.33, 0.71, 0.21, 0.32, 0.93, 0.75, 0.95, 0.09, 0.49, 0.90, 0.44, 0.69, 0.17, 0.25, 0.42, 0.31, 0.64,
    0.83, 0.36, 0.76, 0.28, 0.40, 0.04, 0.14, 0.11, 0.58, 0.63, 0.06, 0.56, 0.75, 0.50, 0.46, 0.63, 0.38, 0.61, 0.92, 0.85, 0.74,
    0.87, 0.30, 0.75, 0.29, 0.52, 0.49, 0.05, 0.57, 0.50, 0.52, 0.03, 0.37, 0.42, 0.74, 0.73, 0.38, 0.40, 0.71, 0.42, 0.61, 0.27,
    0.05, 0.69, 0.67, 0.08, 0.02, 0.00, 0.13, 0.21, 0.86, 0.43, 0.65, 0.97, 0.44, 0.80, 0.39, 0.74, 0.26, 0.89, 0.05, 0.45, 0.99,
    0.14, 0.10, 0.33, 0.14, 0.68, 0.61, 0.53, 0.97, 0.32, 0.44, 0.46, 0.71, 0.23, 0.76, 0.19, 0.98, 0.95, 0.40, 0.88, 0.39, 0.54,
    0.75, 0.76, 0.89, 0.46, 0.59, 0.57, 0.51, 0.69, 0.26, 0.07, 0.52, 0.43, 0.53, 0.61, 0.36, 0.31, 0.68, 0.40, 0.83, 0.76, 0.90,
    0.10, 0.16, 0.61, 0.71, 0.21, 0.02, 0.86, 0.78, 0.92, 0.29, 0.51, 0.73, 0.63, 0.01, 0.51, 0.98, 0.23, 0.46, 0.81, 0.27, 0.41,
    0.69, 0.45, 0.95, 0.85, 0.29, 0.00, 0.54, 0.85, 0.91, 0.12, 0.89, 0.50, 0.98, 0.80, 0.17, 0.53, 0.88, 0.81, 0.97, 0.83, 0.99,
    0.52, 0.60, 0.07, 0.39, 0.32, 0.92, 0.33, 0.89, 0.24, 0.58, 0.53, 0.10, 0.32, 0.99, 0.59, 0.83, 0.72, 0.90, 0.34, 0.22, 0.07,
    0.84, 0.52, 0.72, 0.06, 0.85, 0.58, 0.12, 0.50, 0.60, 0.34, 0.71, 0.58, 0.84, 0.86, 0.06, 0.65, 0.72, 0.14, 0.88, 0.86, 0.14,
    0.88, 0.65, 0.79, 0.80, 0.16, 0.06, 0.64, 0.34, 0.26, 0.13, 0.86, 0.57, 0.84, 0.56, 0.05, 0.22, 0.82, 0.19, 0.91, 0.32, 0.42,
    0.66, 0.80, 0.77, 0.82, 0.51, 0.32, 0.21, 0.96, 0.57, 0.97, 0.49, 0.38, 0.61, 0.99, 0.86, 0.41, 0.04, 0.41, 0.76, 0.00, 0.83,
    0.26, 0.37, 0.92, 0.37, 0.48, 0.44, 0.36, 0.08, 0.42, 0.49, 0.22, 0.83, 0.49, 0.24, 0.10, 0.58, 0.48, 0.72, 0.99, 0.20, 0.73,
    0.43, 0.40, 0.29, 0.36, 0.12, 0.46, 0.59, 0.28, 0.45, 0.47, 0.66, 0.58, 0.94, 0.53, 0.61, 0.53, 0.35, 0.16, 0.02, 0.22, 0.56,
    0.03, 0.32, 0.04, 0.69, 0.51, 0.89, 0.63, 0.70, 0.88, 0.90, 0.37, 0.58, 0.40, 0.57, 0.94, 0.99, 0.47, 0.06, 0.96, 0.12, 0.02,
    0.70, 0.85, 0.62, 0.59, 0.33, 0.63, 0.82, 0.33, 0.46, 0.64, 0.78, 0.64, 0.16, 0.53, 0.44, 0.65, 0.65, 0.10, 0.14, 0.09, 0.10,
    0.52, 0.00, 0.78, 0.57, 0.72, 0.48, 0.91, 0.82, 0.50, 0.24, 0.13, 0.84, 0.68, 0.27, 0.73, 0.31, 0.25, 0.74, 0.74, 0.82, 0.85,
    0.55, 0.26, 0.61, 0.30, 0.79, 0.12, 0.33, 0.04, 0.45, 0.14, 0.97, 0.36, 0.79, 0.83, 0.09, 0.67, 0.34, 0.30, 0.75, 0.25, 0.31,
    0.76, 0.35, 0.39, 0.74, 0.14, 0.69, 0.32, 0.84, 0.62, 0.01, 0.13, 0.75, 0.22, 0.02, 0.74, 0.67, 0.22, 0.63, 0.63, 0.94, 0.10,
    0.17, 0.30, 0.77, 0.47, 0.37, 0.38, 0.51, 0.32, 0.68, 0.03, 0.60, 0.36, 0.48, 0.75, 0.04, 0.22, 0.51, 0.62, 0.19, 0.20, 0.62,
    0.72, 0.45, 0.45, 0.40, 0.25, 0.99, 0.63, 0.36, 0.26, 0.25, 0.32, 0.55, 0.58, 0.66, 0.02, 0.58, 0.58, 0.92, 0.23, 0.04, 0.95,
    0.56, 0.10, 0.33, 0.39, 0.77, 0.26, 0.95, 0.57, 0.10, 0.51, 0.58, 0.96, 0.91, 0.94, 0.65, 0.64, 0.68, 0.05, 0.26, 0.25, 0.78,
    0.17, 0.74, 0.42, 0.27, 0.90, 0.30, 0.77, 0.76, 0.38, 0.33, 0.00, 0.46, 0.18, 0.01, 0.62, 0.13, 0.31, 0.71, 0.75, 0.95, 0.91,
    0.52, 0.95, 0.34, 0.55, 0.54, 0.02, 0.04, 0.24, 0.93, 0.90, 0.64, 0.07, 0.95, 0.07, 0.51, 0.62, 0.71, 0.93, 0.34, 0.37, 0.54,
    0.25, 0.55, 0.35, 0.93, 0.31, 0.63, 0.39, 0.99, 0.89, 0.32, 0.01, 0.30, 0.49, 0.34, 0.66, 0.33, 0.33, 0.90, 0.29, 0.59, 0.52,
    0.26, 0.15, 0.35, 0.32, 0.36, 0.81, 0.53, 0.97, 0.30, 0.25, 0.04, 0.86, 0.48, 0.30, 0.03, 0.62, 0.30, 0.65, 0.11, 0.52, 0.14,
    0.65, 0.74, 0.80, 0.79, 0.88, 0.86, 0.81, 0.25, 0.88, 0.49, 0.14, 0.81, 0.81, 0.44, 0.71, 0.11, 0.12, 0.40, 0.24, 0.87, 0.86,
    0.14, 0.31, 0.45, 0.12, 0.12, 0.47, 0.61, 0.17, 0.15, 0.70, 0.92, 0.00, 0.77, 0.01, 0.94, 0.35, 0.77, 0.52, 0.37, 0.39, 0.69,
    0.88, 0.81, 0.34, 0.76, 0.17, 0.00, 0.69, 0.08, 0.40, 0.82, 0.38, 0.05, 0.02, 0.91, 0.82, 0.65, 0.01, 0.76, 0.43, 0.08, 0.76,
    0.69, 0.21, 0.32, 0.96, 0.12, 0.82, 0.65, 0.22, 0.57, 0.58, 0.05, 0.69, 0.05, 0.56, 0.04, 0.79, 0.97, 0.07, 0.88, 0.78, 0.41,
    0.17, 0.74, 0.73, 0.66, 0.31, 0.45, 0.01, 0.34, 0.49, 0.26, 0.56, 0.27, 0.65, 0.58, 0.86, 0.19, 0.76, 0.87, 0.04, 0.42, 0.35,
    0.02, 0.77, 0.88, 0.77, 0.72, 0.35, 0.24, 0.50, 0.10, 0.64, 0.61, 0.90, 0.91, 0.97, 0.07, 0.72, 0.31, 0.73, 0.56, 0.98, 0.47,
    0.12, 0.63, 0.44, 0.33, 0.60, 0.85, 0.24, 0.67, 0.74, 0.96, 0.22, 0.51, 0.70, 0.05, 0.87, 0.34, 0.68, 0.40, 0.90, 0.81, 0.47,
    0.95, 0.77, 0.43, 0.89, 0.53, 0.56, 0.11, 0.89, 0.39, 0.48, 0.99, 0.83, 0.12, 0.77, 0.24, 0.17, 0.82, 0.54, 0.47, 0.11, 0.09,
    0.40, 0.94, 0.22, 0.72, 0.54, 0.55, 0.95, 0.86, 0.85, 0.05, 0.99, 0.94, 0.11, 0.55, 0.03, 0.11, 0.87, 0.04, 0.26, 0.49, 0.26,
    0.93, 0.32, 0.71, 0.77, 0.05, 0.38, 0.20, 0.26, 0.41, 0.06, 0.61, 0.19, 0.08, 0.45, 0.34, 0.97, 0.80, 0.27, 0.59, 0.53, 0.18,
    0.76, 0.70, 0.13, 0.72, 0.63, 0.01, 0.05, 0.55, 0.09, 0.61, 0.17, 0.69, 0.37, 0.74, 0.19, 0.65, 0.97, 0.49, 0.28, 0.31, 0.19,
    0.57, 0.51, 0.55, 0.68, 0.02, 0.13, 0.21, 0.01, 0.19, 0.85, 0.65, 0.16, 0.04, 0.25, 0.71, 0.45, 0.54, 0.54, 0.24, 0.17, 0.16,
    0.96, 0.94, 0.75, 0.65, 0.22, 0.12, 0.93, 0.28, 0.19, 0.92, 0.16, 0.30, 0.18, 0.01, 0.18, 0.45, 0.53, 0.17, 0.72, 0.54, 0.60,
    0.32, 0.20, 0.80, 0.35, 0.52, 0.74, 0.88, 0.45, 0.97, 0.53, 0.79, 0.08, 0.09, 0.61, 0.39, 0.76, 0.57, 0.62, 0.17, 0.55, 0.54,
    0.74, 0.70, 0.62, 0.16, 0.97, 0.61, 0.16, 0.80, 0.14, 0.27, 0.56, 0.15, 0.09, 0.68, 0.04, 0.35, 0.25, 0.94, 0.58, 0.07, 0.41,
    0.37, 0.15, 0.73, 0.19, 0.85, 0.94, 0.09, 0.66, 0.63, 0.51, 0.64, 0.90, 0.76, 0.65, 0.43, 0.24, 0.20, 0.22, 0.86, 0.12, 0.12,
    0.02, 0.73, 0.01, 0.92, 0.97, 0.23, 0.04, 0.96, 0.21, 0.54, 0.65, 0.59, 0.85, 0.48, 0.52, 1.00, 0.56, 0.20, 0.17, 0.87, 0.36,
    0.95, 0.81, 0.94, 0.56, 0.04, 0.85, 0.58, 0.87, 0.28, 0.86, 0.20, 0.28, 0.14, 0.14, 0.10, 0.54, 0.68, 0.21, 0.91, 0.10, 0.11,
    0.88, 0.69, 0.18, 0.99, 0.03, 0.11, 0.28, 0.73, 0.71, 0.80, 0.55, 0.60, 0.35, 0.38, 0.10, 0.12, 0.70, 0.02, 0.45, 0.66, 0.48,
    0.07, 0.77, 0.48, 0.29, 0.17, 0.61, 0.92, 0.26, 0.87, 0.66, 0.02, 0.85, 0.40, 0.84, 0.20, 0.52, 0.21, 0.49, 0.92, 0.51, 0.99,
    0.27, 0.60, 0.97, 0.75, 0.69, 0.74, 0.33, 0.29, 0.14, 0.60, 0.92, 0.76, 0.99, 0.99, 0.43, 0.52, 0.82, 0.77, 0.77, 0.56, 0.59,
    0.84, 0.94, 0.78, 0.77, 0.16, 0.86, 0.95, 0.27, 0.19, 0.16, 0.02, 0.41, 0.18, 0.13, 0.73, 0.76, 0.70, 0.67, 0.52, 0.06, 0.98,
    0.99, 0.88, 0.76, 0.77, 0.60, 0.11, 0.82, 0.68, 0.08, 0.61, 0.19, 0.56, 0.25, 0.64, 0.78, 0.05, 0.07, 1.00, 0.68, 0.31, 0.24,
    0.42, 0.95, 0.69, 0.99, 0.55, 0.42, 0.21, 0.54, 0.64, 0.52, 0.13, 0.17, 0.79, 0.21, 0.11, 0.16, 0.90, 0.62, 0.77, 0.23, 0.63,
    0.19, 0.52, 0.07, 0.52, 0.23, 0.22, 0.50, 0.57, 0.70, 0.30, 0.29, 0.08, 0.32, 0.35, 0.49, 0.42, 0.70, 0.45, 0.84, 0.53, 0.59,
    0.70, 0.68, 0.93, 0.74, 0.96, 0.34, 0.97, 0.60, 0.31, 0.14, 0.52, 0.96, 0.84, 0.75, 0.05, 0.77, 0.69, 0.84, 0.10, 0.39, 0.34,
    0.77, 0.76, 0.50, 0.36, 0.10, 0.15, 0.01, 0.77, 0.92, 0.03, 0.80, 0.28, 0.51, 0.69, 0.34, 0.98, 0.99, 0.67, 0.98, 0.12, 0.01,
    0.85, 0.69, 0.31, 0.52, 0.40, 0.45, 0.57, 0.52, 0.61, 0.33, 0.56, 0.86, 0.36, 0.28, 0.49, 0.32, 0.90, 0.61, 0.74, 0.03, 0.10,
    0.40, 0.42, 0.07, 0.79, 0.56, 0.50, 0.68, 0.69, 0.33, 0.48, 0.58, 0.00, 0.14, 0.39, 0.04, 0.05, 0.51, 0.61, 0.24, 0.90, 0.81,
    0.86, 0.46, 0.82, 0.67, 0.09, 0.01, 0.94, 0.54, 0.73, 0.16, 0.70, 0.65, 0.90, 0.30, 0.36, 0.61, 0.88, 0.27, 0.10, 0.33, 0.64,
    1.00, 0.17, 0.81, 0.22, 0.51, 0.75, 0.11, 0.80, 0.08, 0.32, 0.03, 0.84, 0.60, 0.17, 0.12, 0.48, 0.01, 0.79, 0.97, 0.04, 0.59,
    0.61, 0.01, 0.77, 0.26, 0.65, 0.21, 0.09, 0.25, 0.88, 0.76, 0.49, 0.12, 0.21, 0.96, 0.14, 0.29, 0.55, 0.37, 0.27, 0.68, 0.16,
    0.47, 0.12, 0.29, 0.36, 0.64, 0.46, 0.28, 0.68, 0.47, 0.72, 0.39, 0.73, 0.57, 0.25, 0.54, 0.14, 0.01, 0.69, 0.18, 0.59, 0.47,
    0.44, 0.31, 0.52, 0.97, 0.77, 0.41, 0.87, 0.53, 0.99, 0.97, 0.03, 0.17, 0.84, 0.36, 0.25, 0.54, 0.67, 0.32, 0.86, 0.50, 0.23,
    0.67, 0.47, 0.15, 0.07, 0.60, 0.06, 0.97, 0.20, 0.99, 0.32, 0.93, 0.64, 0.10, 0.26, 0.39, 0.42, 0.56, 0.59, 0.31, 0.09, 0.69,
    0.46, 0.68, 0.57, 0.01, 0.80, 0.77, 0.40, 0.78, 0.02, 0.43, 0.42, 0.20, 0.22, 0.05, 0.66, 0.96, 0.93, 0.20, 0.54, 0.79, 0.33,
    0.18, 0.02, 0.74, 0.27, 0.19, 0.90, 0.83, 0.79, 0.06, 0.93, 0.45, 0.26, 0.73, 0.04, 0.15, 0.42, 0.97, 0.03, 0.06, 0.12, 0.30,
    0.26, 0.17, 0.70, 0.73, 0.39, 0.54, 0.58, 0.82, 0.42, 0.90, 0.81, 0.54, 0.29, 0.73, 0.70, 0.15, 0.90, 0.19, 0.06, 0.00, 0.77,
    0.00, 0.25, 0.97, 0.10, 0.59, 0.06, 0.46, 0.62, 0.53, 0.67, 0.15, 0.70, 0.03, 0.30, 0.93, 0.18, 0.15, 0.41, 0.43, 0.28, 0.19,
    0.23, 0.68, 0.11, 0.12, 0.16, 0.91, 0.85, 0.95, 0.94, 0.06, 0.81, 0.35, 0.81, 0.95, 0.71, 0.85, 0.56, 0.52, 0.95, 0.43, 0.10,
    0.22, 0.68, 0.41, 0.77, 0.53, 0.13, 0.66, 0.95, 0.51, 0.88, 0.97, 0.28, 0.53, 0.02, 0.52, 0.21, 0.41, 0.06, 0.78, 0.57, 0.64,
    0.39, 0.10, 0.37, 0.61, 0.56, 0.75, 0.35, 0.94, 0.68, 0.54, 0.06, 0.57, 0.82, 0.34, 0.81, 0.87, 0.76, 0.40, 0.08, 0.57, 0.41,
    0.37, 0.14, 0.67, 0.44, 0.36, 0.94, 0.39, 0.17, 0.06, 0.99, 0.38, 0.84, 0.08, 0.50, 0.40, 0.27, 0.74, 0.65, 0.07, 0.29, 0.27,
    0.91, 0.23, 0.51, 0.96, 0.99, 0.51
};

class hydraLinearTest : public tardigradeHydra::hydraBase{

    public:

        hydraLinearTest(
            const unsigned int &_nphases,               const unsigned int &_active_phase,
            const unsigned int &_num_phase_dof,         const unsigned int &_num_add_dof,
            const floatType &t,                         const floatType &dt,
            const std::vector< double > &additionalDOF, const std::vector< double > &previousAdditionalDOF
        ) : tardigradeHydra::hydraBase(
                t, dt,
                _getAdditionalDOFTemperature( _nphases, _active_phase, additionalDOF ),
                _getAdditionalDOFTemperature( _nphases, _active_phase, additionalDOF ),
                _getAdditionalDOFDeformationGradient( _nphases, _active_phase, _num_phase_dof, _num_add_dof, additionalDOF         ),
                _getAdditionalDOFDeformationGradient( _nphases, _active_phase, _num_phase_dof, _num_add_dof, previousAdditionalDOF ),
                additionalDOF, previousAdditionalDOF,
                std::vector< double >( 13, 0 ),
                linear_test_params,
                1, 13
            ){

            nphases = _nphases;

            active_phase = _active_phase;

            num_phase_dof = _num_phase_dof;

            num_add_dof = _num_add_dof;

        }

        //! The residual class
        tardigradeHydra::linearTestMaterial::residual residual; //!< The residual class

        virtual std::vector< double > getFullTangent( ){

            constexpr unsigned int dim = 3;

            computeTangents( );

            computedXdAdditionalDOF( );

            std::vector< double > dFdGradU( 81, 0 );

            for ( unsigned int i = 0; i < 3; ++i ){

                for ( unsigned int I = 0; I < 3; ++I ){

                    for ( unsigned int a = 0; a < 3; ++a ){

                        for ( unsigned int b = 0; b < 3; ++b ){

                            dFdGradU[ dim * dim * dim * i + dim * dim * I + dim * a + b ]
                                += ( *getDeformationGradient( ) )[ dim * i + a ] * ( *getDeformationGradient( ) )[ dim * b + I ];

                        }

                    }

                }

            }

            std::vector< double > full_jacobian( getFlatdXdAdditionalDOF( )->size( ), 0 );

            std::copy(
                std::begin( *getFlatdXdAdditionalDOF( ) ),
                std::end( *getFlatdXdAdditionalDOF( ) ),
                std::begin( full_jacobian )
            );

            // Incorporate the Jacobian of the temperature and deformation gradient to the full Jacobian
            unsigned int offset = nphases * num_phase_dof + num_add_dof + nphases * 3 + 9 * active_phase;

            for ( unsigned int I = 0; I < getNumUnknowns( ); ++I ){

                full_jacobian[ getAdditionalDOF( )->size( ) * I + nphases * ( 1 + 3 + 3 ) + active_phase ]
                    += ( *getFlatdXdT( ) )[ I ];

                for ( unsigned int ij = 0; ij < dim * dim; ++ij ){

                    for ( unsigned int ab = 0; ab < dim * dim; ++ab ){

                        full_jacobian[ getAdditionalDOF( )->size( ) * I + offset + ab ]
                            += ( *getFlatdXdF( ) )[ dim * dim * I + ij ]
                             * dFdGradU[ dim * dim * ij + ab ];

                    }

                }

            }

            return full_jacobian;

        }

    protected:

        unsigned int nphases;

        unsigned int active_phase;

        unsigned int num_phase_dof;

        unsigned int num_add_dof;

        virtual floatType _getAdditionalDOFTemperature(
            const unsigned int _nphases, const unsigned int _active_phase,
            const std::vector< floatType > &additional_dof
        ){

            unsigned int offset = _nphases * ( 1 + 3 + 3 ) + _active_phase;

            return additional_dof[ offset ];

        }

        virtual std::vector< floatType > _getAdditionalDOFDeformationGradient(
            const unsigned int _nphases, const unsigned int _active_phase,
            const unsigned int _num_phase_dof, const unsigned int _num_add_dof,
            const std::vector< floatType > &additional_dof
        ){

            unsigned int offset = _nphases * _num_phase_dof + _num_add_dof + _nphases * 3 + 9 * _active_phase;

            std::vector< floatType > gradW(
                std::begin( additional_dof ) + offset,
                std::begin( additional_dof ) + offset + 9
            );

            std::vector< floatType > F;

            tardigradeConstitutiveTools::computeDeformationGradient( gradW, F, true );

            return F;

        }

    private:

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses( ) override{
            /*!
             * Set the vector of residual classes (in this case, a single residual)
             */

            std::vector< tardigradeHydra::residualBase*  > residuals( 1 );

            TARDIGRADE_ERROR_TOOLS_CATCH( residual = tardigradeHydra::linearTestMaterial::residual( this, 22, linear_test_params ) );

            residuals[ 0 ] = &residual;

            setResidualClasses( residuals );

        }

};

BOOST_AUTO_TEST_CASE( test_linearHydraTest, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the linear hydra test model Jacobian
     */

    std::vector< floatType > additional_dof =
    {
        +0.392, -0.428, -0.546, +0.102, +0.438, -0.154,
        +0.962, +0.370, -0.038, -0.216, -0.314, +0.458,
        -0.122, -0.880, -0.204, +0.476, -0.636, -0.650,
        +0.064, +0.064, +0.268, +0.698, +0.448, +0.222,
        +0.444, -0.354, -0.276, -0.544, -0.412, +0.262,
        -0.816, -0.132, -0.138, -0.012, -0.148, -0.376,
        -0.148, +0.786, +0.888, +0.004, +0.248, -0.768,
        -0.366, -0.170, +0.732, -0.500, -0.034, +0.972,
        +0.038, +0.226, -0.758, +0.652, +0.206, +0.090,
        -0.314, -0.392, -0.166, +0.362, +0.750, +0.020,
        +0.338, +0.172, +0.250, +0.350, +0.684, -0.834,
        +0.528, -0.512, -0.612, +0.144, -0.808, +0.770,
        +0.254, +0.446, -0.968, +0.188, +0.114, -0.682,
        -0.694, +0.392, -0.362, +0.384, +0.108, -0.222,
        +0.850, +0.684, -0.286, -0.912, -0.390, -0.204,
        +0.410, +0.990, -0.288, +0.526, +0.186, +0.384,
        -0.698, -0.202, -0.518, -0.314, +0.026, +0.334,
        -0.788, -0.738, -0.356, +0.324, +0.694, +0.106,
        +0.708, -0.230, -0.366, -0.292, -0.658, +0.658,
        -0.322, +0.104, +0.158, +0.044, -0.994, +0.976,
        +0.810, -0.584, -0.416, +0.040, +0.804, +0.968,
        -0.484, +0.128, +0.614, -0.212, +0.462, -0.678,
        +0.202, +0.732, +0.968, -0.842, -0.144, -0.590,
        -0.098, +0.096, -0.814, -0.406, +0.856, +0.138,
        -0.086, +0.508, +0.484, -0.902, +0.418, +0.678,
        -0.668, +0.562, -0.426, -0.388, +0.330, -0.778,
        +0.330, +0.776, +0.392, -0.120, -0.124, +0.530,
        +0.132, -0.830
    };

    std::vector< floatType > previous_additional_dof( additional_dof.size( ), 0 );

    const unsigned int nphases = 4;

    const unsigned int active_phase = 3;

    const unsigned int num_phase_dof = 9;

    const unsigned int num_add_dof   = 5;

    std::vector< floatType > answer =
    {
        +3.137965460e+00, -9.070686731e-01, +5.603138559e+00, +1.108519558e+00, +2.602408209e+00, +9.284901523e-01,
        +2.837889005e+00, +1.476724760e+00, +3.349219055e+00, +8.644615772e-01, +2.088948758e+00, -1.135202759e+00,
        +4.019674120e+00, +1.078555126e+00, +2.870173714e+00, +1.248529766e+00, +1.673645827e+00, +1.566445004e+00,
        +7.420031308e-01, +2.596281643e+00, +4.319387829e-01, +2.005971987e+00
    };

    hydraLinearTest linearHydra(
        nphases, active_phase, num_phase_dof, num_add_dof,
        0, 0.1, additional_dof, previous_additional_dof
    );

    linearHydra.evaluate( );

    BOOST_TEST( answer == *linearHydra.getUnknownVector( ), CHECK_PER_ELEMENT );

    std::vector< double > dXdAdditionalDOF = linearHydra.getFullTangent( );

    // Check the Jacobian
    floatType eps = 1e-6;
    {

        constexpr unsigned int vardim = 164;
        constexpr unsigned int outdim = 22;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( additional_dof[ i ] ) + eps;

            std::vector< floatType > xp = additional_dof;
            std::vector< floatType > xm = additional_dof;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::vector< floatType > vp, vm;

            hydraLinearTest linearHydrap(
                nphases, active_phase, num_phase_dof, num_add_dof,
                0, 0.1, xp, previous_additional_dof
            );

            hydraLinearTest linearHydram(
                nphases, active_phase, num_phase_dof, num_add_dof,
                0, 0.1, xm, previous_additional_dof
            );

            linearHydrap.evaluate( );
            linearHydram.evaluate( );

            vp = *linearHydrap.getUnknownVector( );
            vm = *linearHydram.getUnknownVector( );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dXdAdditionalDOF[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

template<
    int dim, int node_count, int nphases, int num_additional_dof,
    class xi_in, typename dt_type,
    class density_t_in, class density_tp1_in,
    class u_t_in,       class u_tp1_in,
    class w_t_in,       class w_tp1_in,
    class theta_t_in,   class theta_tp1_in,
    class e_t_in,       class e_tp1_in,
    class z_t_in,       class z_tp1_in,
    class vf_t_in,      class vf_tp1_in,
    class umesh_t_in, class umesh_tp1_in, class density_dot_t_in, class v_t_in, class e_dot_t_in,
    class X_in, typename alpha_type, typename beta_type, class value_out,
    int material_response_size = 22
>
void evaluate_at_nodes(
    const xi_in &xi_begin, const xi_in &xi_end, dt_type dt,
    const density_t_in &density_t_begin,     const density_t_in &density_t_end,
    const density_tp1_in &density_tp1_begin, const density_tp1_in &density_tp1_end,
    const u_t_in       &u_t_begin,           const u_t_in       &u_t_end,
    const u_tp1_in     &u_tp1_begin,         const u_tp1_in     &u_tp1_end,
    const w_t_in       &w_t_begin,           const w_t_in       &w_t_end,
    const w_tp1_in     &w_tp1_begin,         const w_tp1_in     &w_tp1_end,
    const theta_t_in   &theta_t_begin,       const theta_t_in   &theta_t_end,
    const theta_tp1_in &theta_tp1_begin,     const theta_tp1_in &theta_tp1_end,
    const e_t_in       &e_t_begin,           const e_t_in       &e_t_end,
    const e_tp1_in     &e_tp1_begin,         const e_tp1_in     &e_tp1_end,
    const z_t_in       &z_t_begin,           const z_t_in       &z_t_end,
    const z_tp1_in     &z_tp1_begin,         const z_tp1_in     &z_tp1_end,
    const vf_t_in      &vf_t_begin,          const vf_t_in      &vf_t_end,
    const vf_tp1_in    &vf_tp1_begin,        const vf_tp1_in    &vf_tp1_end,
    const umesh_t_in &umesh_t_begin, const umesh_t_in &umesh_t_end,
    const umesh_tp1_in &umesh_tp1_begin, const umesh_tp1_in &umesh_tp1_end,
    const density_dot_t_in &density_dot_t_begin, const density_dot_t_in &density_dot_t_end,
    const v_t_in &v_t_begin, const v_t_in &v_t_end, const e_dot_t_in &e_dot_t_begin, const e_dot_t_in &e_dot_t_end,
    const X_in &X_begin, const X_in &X_end,
    const alpha_type &alpha, const beta_type &beta, value_out value_begin, value_out value_end,
    const int active_phase = -1
){

    // Update the mesh nodes
    std::array< typename std::iterator_traits<umesh_t_in  >::value_type, dim * node_count > x_t;
    std::array< typename std::iterator_traits<umesh_tp1_in>::value_type, dim * node_count > x_tp1;

    std::transform( X_begin, X_end,   umesh_t_begin,   std::begin( x_t ), std::plus<typename std::iterator_traits<umesh_t_in  >::value_type>( ) );
    std::transform( X_begin, X_end, umesh_tp1_begin, std::begin( x_tp1 ), std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>( ) );

    // Calculate the current rates of change
    std::array< typename std::iterator_traits<density_tp1_in>::value_type, node_count * nphases > density_dot_tp1, e_dot_tp1;

    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count * nphases > v_tp1;

    floatType dDensityDotdDensity;

    compute_current_rate_of_change(
        dt, density_t_begin, density_t_end, density_tp1_begin, density_tp1_end,
        density_dot_t_begin, density_dot_t_end, alpha,
        std::begin( density_dot_tp1 ), std::end( density_dot_tp1 ),
        dDensityDotdDensity
    );

    floatType dUDotdU;

    compute_current_rate_of_change(
        dt, u_t_begin, u_t_end, u_tp1_begin, u_tp1_end,
        v_t_begin, v_t_end, alpha,
        std::begin( v_tp1 ), std::end( v_tp1 ),
        dUDotdU
    );

    floatType dEDotdE;

    compute_current_rate_of_change(
        dt, e_t_begin, e_t_end, e_tp1_begin, e_tp1_end,
        e_dot_t_begin, e_dot_t_end, alpha,
        std::begin( e_dot_tp1 ), std::end( e_dot_tp1 ),
        dEDotdE
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
    > density_tp1_p, density_dot_tp1_p, theta_tp1_p, e_tp1_p, e_dot_tp1_p, vf_tp1_p;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * nphases
    > v_tp1_p, w_tp1_p;

    std::array<
         typename std::iterator_traits<z_tp1_in>::value_type, num_additional_dof
    > z_tp1_p;

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
        xi_begin, xi_end, std::cbegin( v_tp1 ), std::cend( v_tp1 ),
        std::begin( v_tp1_p ), std::end( v_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, w_tp1_begin, w_tp1_end,
        std::begin( w_tp1_p ), std::end( w_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, theta_tp1_begin, theta_tp1_end,
        std::begin( theta_tp1_p ), std::end( theta_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, e_tp1_begin, e_tp1_end,
        std::begin( e_tp1_p ), std::end( e_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::begin( e_dot_tp1 ), std::end( e_dot_tp1 ),
        std::begin( e_dot_tp1_p ), std::end( e_dot_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, z_tp1_begin, z_tp1_end,
        std::begin( z_tp1_p ), std::end( z_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, vf_tp1_begin, vf_tp1_end,
        std::begin( vf_tp1_p ), std::end( vf_tp1_p )
    );

   // Compute the gradients at the local point

    std::array<
        typename std::iterator_traits<density_tp1_in>::value_type, dim * nphases
    > grad_density_tp1, grad_theta_tp1, grad_e_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * dim * nphases
    > grad_velocity_tp1, grad_w_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type,dim * num_additional_dof
    > grad_z_tp1;

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( grad_density_tp1 ), std::end( grad_density_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, std::cbegin( v_tp1 ), std::cend( v_tp1 ),
        std::begin( grad_velocity_tp1 ), std::end( grad_velocity_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, w_tp1_begin, w_tp1_end,
        std::begin( grad_w_tp1 ), std::end( grad_w_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, theta_tp1_begin, theta_tp1_end,
        std::begin( grad_theta_tp1 ), std::end( grad_theta_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, e_tp1_begin, e_tp1_end,
        std::begin( grad_e_tp1 ), std::end( grad_e_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, z_tp1_begin, z_tp1_end,
        std::begin( grad_z_tp1 ), std::end( grad_z_tp1 )
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

   std::vector< floatType > dof_vector( nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 + 9 + 3 + 3 ) + num_additional_dof + 3 * num_additional_dof, 0 );

    std::copy(
        std::begin( density_tp1_p ),
        std::end(   density_tp1_p ),
        std::begin( dof_vector ) + nphases * 0
    );

    std::copy(
        std::begin( w_tp1_p ),
        std::end(   w_tp1_p ),
        std::begin( dof_vector ) + nphases * 1
    );

    std::copy(
        std::begin( v_tp1_p ),
        std::end(   v_tp1_p ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 )
    );

    std::copy(
        std::begin( theta_tp1_p ),
        std::end(   theta_tp1_p ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 )
    );

    std::copy(
        std::begin( e_tp1_p ),
        std::end(   e_tp1_p ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 )
    );

    std::copy(
        std::begin( z_tp1_p ),
        std::end(   z_tp1_p ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 )
    );

    std::copy(
        std::begin( grad_density_tp1 ),
        std::end(   grad_density_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 ) + num_additional_dof
    );

    std::copy(
        std::begin( grad_w_tp1 ),
        std::end(   grad_w_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 + 3 ) + num_additional_dof
    );

    std::copy(
        std::begin( grad_velocity_tp1 ),
        std::end(   grad_velocity_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 ) + num_additional_dof
    );

    std::copy(
        std::begin( grad_theta_tp1 ),
        std::end(   grad_theta_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 + 9 ) + num_additional_dof
    );

    std::copy(
        std::begin( grad_e_tp1 ),
        std::end(   grad_e_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 + 9 + 3 ) + num_additional_dof
    );

    std::copy(
        std::begin( grad_z_tp1 ),
        std::end(   grad_z_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 + 9 + 3 + 3 ) + num_additional_dof
    );

    std::vector< floatType > previous_dof_vector( dof_vector.size( ) );

    std::array< floatType, node_count> Ns;
    e.GetShapeFunctions( xi_begin, xi_end, std::begin( Ns ), std::end( Ns ) );

    std::array< floatType, node_count * dim > dNdx;
    e.GetGlobalShapeFunctionGradients( xi_begin, xi_end, std::cbegin( x_tp1 ), std::cend( x_tp1 ),
                                       std::begin( dNdx ), std::end( dNdx ) );

    unsigned int low_bound  = 0;
    unsigned int high_bound = nphases;

    if ( active_phase >= 0 ){

        low_bound  = active_phase;
        high_bound = active_phase + 1;

    }

    // Assemble the material response
    std::array< floatType, nphases * material_response_size > material_response;
    std::fill( std::begin( material_response ), std::end( material_response ), 0 );

    for ( unsigned int j = low_bound; j < high_bound; ++j ){

        hydraLinearTest linearTest(
            nphases, j, 9, num_additional_dof,
            0, 0.1, dof_vector, previous_dof_vector
        );

        linearTest.evaluate( );

        std::copy(
            std::begin( *linearTest.getUnknownVector( ) ),
            std::end(   *linearTest.getUnknownVector( ) ),
            std::begin( material_response ) + material_response_size * j
        );

    }

    std::fill( value_begin, value_end, 0 );

    for ( unsigned int i = 0; i < node_count; ++i ){

        if ( active_phase >= 0 ){

            unsigned int j = active_phase;

            tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<3,3,0,20,17,14,21>(
                density_tp1_p[ j ], density_dot_tp1_p[ j ],
                std::cbegin( grad_density_tp1 )  + 3 * j, std::cbegin( grad_density_tp1 )  + 3 * ( j + 1 ),
                e_tp1_p[ j ], e_dot_tp1_p[ j ],
                std::cbegin( grad_e_tp1 ) + 3 * j,        std::cbegin( grad_e_tp1 ) + 3 * ( j + 1 ),
                std::cbegin( v_tp1_p )           + 3 * j, std::cbegin( v_tp1_p )           + 3 * ( j + 1 ),
                std::cbegin( grad_velocity_tp1 ) + 9 * j, std::cbegin( grad_velocity_tp1 ) + 9 * ( j + 1 ),
                std::cbegin( material_response ) + material_response_size * j,
                std::cbegin( material_response ) + material_response_size * ( j + 1 ),
                vf_tp1_p[ j ],
                Ns[ i ], std::begin( dNdx ) + 3 * i, std::begin( dNdx ) + 3 * ( i + 1 ),
                *( value_begin + nphases * i + j )
            );

            std::transform(
                value_begin + nphases * i + j,
                value_begin + nphases * i + ( j + 1 ),
                value_begin + nphases * i + j,
                std::bind(
                    std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                    std::placeholders::_1,
                    J
                )
            );

        }
        else{

            tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<3,3,0,20,17,14,21>(
                std::cbegin( density_tp1_p ),     std::cend( density_tp1_p ),
                std::cbegin( density_dot_tp1_p ), std::cend( density_dot_tp1_p ),
                std::cbegin( grad_density_tp1 ),  std::cend( grad_density_tp1 ),
                std::cbegin( e_tp1_p ),           std::cend( e_tp1_p ),
                std::cbegin( e_dot_tp1_p ),       std::cend( e_dot_tp1_p ),
                std::cbegin( grad_e_tp1 ),        std::cend( grad_e_tp1 ),
                std::cbegin( v_tp1_p ),           std::cend( v_tp1_p ),
                std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
                std::cbegin( material_response ), std::cend( material_response ),
                std::cbegin( vf_tp1_p ),          std::cend( vf_tp1_p ),
                Ns[ i ], std::begin( dNdx ) + 3 * i, std::begin( dNdx ) + 3 * ( i + 1 ),
                value_begin + nphases * i,        value_begin + nphases * ( i + 1 )
            );

            std::transform(
                value_begin + nphases * i,
                value_begin + nphases * ( i + 1 ),
                value_begin + nphases * i,
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
    int dim, int node_count, int nphases, int num_additional_dof,
    class xi_in, typename dt_type,
    class density_t_in, class density_tp1_in,
    class u_t_in,       class u_tp1_in,
    class w_t_in,       class w_tp1_in,
    class theta_t_in,   class theta_tp1_in,
    class e_t_in,       class e_tp1_in,
    class z_t_in,       class z_tp1_in,
    class vf_t_in,      class vf_tp1_in,
    class umesh_t_in, class umesh_tp1_in, class density_dot_t_in, class v_t_in, class e_dot_t_in,
    class X_in, typename alpha_type, typename beta_type, class value_out,
    class dRdRho_iter,   class dRdU_iter, class dRdW_iter,
    class dRdTheta_iter, class dRdE_iter, class dRdZ_iter,
    class dRdVF_iter, class dRdUMesh_iter,
    int material_response_size = 22
>
void evaluate_at_nodes(
    const xi_in &xi_begin, const xi_in &xi_end, dt_type dt,
    const density_t_in &density_t_begin,     const density_t_in &density_t_end,
    const density_tp1_in &density_tp1_begin, const density_tp1_in &density_tp1_end,
    const u_t_in       &u_t_begin,           const u_t_in       &u_t_end,
    const u_tp1_in     &u_tp1_begin,         const u_tp1_in     &u_tp1_end,
    const w_t_in       &w_t_begin,           const w_t_in       &w_t_end,
    const w_tp1_in     &w_tp1_begin,         const w_tp1_in     &w_tp1_end,
    const theta_t_in   &theta_t_begin,       const theta_t_in   &theta_t_end,
    const theta_tp1_in &theta_tp1_begin,     const theta_tp1_in &theta_tp1_end,
    const e_t_in       &e_t_begin,           const e_t_in       &e_t_end,
    const e_tp1_in     &e_tp1_begin,         const e_tp1_in     &e_tp1_end,
    const z_t_in       &z_t_begin,           const z_t_in       &z_t_end,
    const z_tp1_in     &z_tp1_begin,         const z_tp1_in     &z_tp1_end,
    const vf_t_in      &vf_t_begin,          const vf_t_in      &vf_t_end,
    const vf_tp1_in    &vf_tp1_begin,        const vf_tp1_in    &vf_tp1_end,
    const umesh_t_in &umesh_t_begin, const umesh_t_in &umesh_t_end,
    const umesh_tp1_in &umesh_tp1_begin, const umesh_tp1_in &umesh_tp1_end,
    const density_dot_t_in &density_dot_t_begin, const density_dot_t_in &density_dot_t_end,
    const v_t_in &v_t_begin, const v_t_in &v_t_end, const e_dot_t_in &e_dot_t_begin, const e_dot_t_in &e_dot_t_end,
    const X_in &X_begin, const X_in &X_end,
    const alpha_type &alpha, const beta_type &beta, value_out value_begin, value_out value_end,
    dRdRho_iter   dRdRho_begin,   dRdRho_iter   dRdRho_end,
    dRdU_iter     dRdU_begin,     dRdU_iter     dRdU_end,
    dRdW_iter     dRdW_begin,     dRdW_iter     dRdW_end,
    dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end,
    dRdE_iter     dRdE_begin,     dRdE_iter     dRdE_end,
    dRdZ_iter     dRdZ_begin,     dRdZ_iter     dRdZ_end,
    dRdVF_iter    dRdVF_begin,    dRdVF_iter    dRdVF_end,
    dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end,
    const int active_phase = -1
){

    // Update the mesh nodes
    std::array< typename std::iterator_traits<umesh_t_in  >::value_type, dim * node_count > x_t;
    std::array< typename std::iterator_traits<umesh_tp1_in>::value_type, dim * node_count > x_tp1;

    std::transform( X_begin, X_end,   umesh_t_begin,   std::begin( x_t ), std::plus<typename std::iterator_traits<umesh_t_in  >::value_type>( ) );
    std::transform( X_begin, X_end, umesh_tp1_begin, std::begin( x_tp1 ), std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>( ) );

    // Calculate the current rates of change
    std::array< typename std::iterator_traits<density_tp1_in>::value_type, node_count * nphases > density_dot_tp1, e_dot_tp1;

    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count * nphases > v_tp1;

    floatType dDensityDotdDensity;

    compute_current_rate_of_change(
        dt, density_t_begin, density_t_end, density_tp1_begin, density_tp1_end,
        density_dot_t_begin, density_dot_t_end, alpha,
        std::begin( density_dot_tp1 ), std::end( density_dot_tp1 ),
        dDensityDotdDensity
    );

    floatType dUDotdU;

    compute_current_rate_of_change(
        dt, u_t_begin, u_t_end, u_tp1_begin, u_tp1_end,
        v_t_begin, v_t_end, alpha,
        std::begin( v_tp1 ), std::end( v_tp1 ),
        dUDotdU
    );

    floatType dEDotdE;

    compute_current_rate_of_change(
        dt, e_t_begin, e_t_end, e_tp1_begin, e_tp1_end,
        e_dot_t_begin, e_dot_t_end, alpha,
        std::begin( e_dot_tp1 ), std::end( e_dot_tp1 ),
        dEDotdE
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
    > density_tp1_p, density_dot_tp1_p, theta_tp1_p, e_tp1_p, e_dot_tp1_p, vf_tp1_p;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * nphases
    > v_tp1_p, w_tp1_p;

    std::array<
         typename std::iterator_traits<z_tp1_in>::value_type, num_additional_dof
    > z_tp1_p;

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
        xi_begin, xi_end, std::cbegin( v_tp1 ), std::cend( v_tp1 ),
        std::begin( v_tp1_p ), std::end( v_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, w_tp1_begin, w_tp1_end,
        std::begin( w_tp1_p ), std::end( w_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, theta_tp1_begin, theta_tp1_end,
        std::begin( theta_tp1_p ), std::end( theta_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, e_tp1_begin, e_tp1_end,
        std::begin( e_tp1_p ), std::end( e_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::begin( e_dot_tp1 ), std::end( e_dot_tp1 ),
        std::begin( e_dot_tp1_p ), std::end( e_dot_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, z_tp1_begin, z_tp1_end,
        std::begin( z_tp1_p ), std::end( z_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, vf_tp1_begin, vf_tp1_end,
        std::begin( vf_tp1_p ), std::end( vf_tp1_p )
    );

   // Compute the gradients at the local point

    std::array<
        typename std::iterator_traits<density_tp1_in>::value_type, dim * nphases
    > grad_density_tp1, grad_theta_tp1, grad_e_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * dim * nphases
    > grad_velocity_tp1, grad_w_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type,dim * num_additional_dof
    > grad_z_tp1;

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( grad_density_tp1 ), std::end( grad_density_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, std::cbegin( v_tp1 ), std::cend( v_tp1 ),
        std::begin( grad_velocity_tp1 ), std::end( grad_velocity_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, w_tp1_begin, w_tp1_end,
        std::begin( grad_w_tp1 ), std::end( grad_w_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, theta_tp1_begin, theta_tp1_end,
        std::begin( grad_theta_tp1 ), std::end( grad_theta_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, e_tp1_begin, e_tp1_end,
        std::begin( grad_e_tp1 ), std::end( grad_e_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, z_tp1_begin, z_tp1_end,
        std::begin( grad_z_tp1 ), std::end( grad_z_tp1 )
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

    constexpr unsigned int num_dof = nphases * ( 1 + 3 + 3 + 1 + 1 ) + num_additional_dof;

    constexpr unsigned int dof_vector_size = ( nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 + 9 + 3 + 3 ) + num_additional_dof + 3 * num_additional_dof );

    std::vector< floatType > dof_vector( nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 + 9 + 3 + 3 ) + num_additional_dof + 3 * num_additional_dof, 0 );

    std::copy(
        std::begin( density_tp1_p ),
        std::end(   density_tp1_p ),
        std::begin( dof_vector ) + nphases * 0
    );

    std::copy(
        std::begin( w_tp1_p ),
        std::end(   w_tp1_p ),
        std::begin( dof_vector ) + nphases * 1
    );

    std::copy(
        std::begin( v_tp1_p ),
        std::end(   v_tp1_p ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 )
    );

    std::copy(
        std::begin( theta_tp1_p ),
        std::end(   theta_tp1_p ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 )
    );

    std::copy(
        std::begin( e_tp1_p ),
        std::end(   e_tp1_p ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 )
    );

    std::copy(
        std::begin( z_tp1_p ),
        std::end(   z_tp1_p ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 )
    );

    std::copy(
        std::begin( grad_density_tp1 ),
        std::end(   grad_density_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 ) + num_additional_dof
    );

    std::copy(
        std::begin( grad_w_tp1 ),
        std::end(   grad_w_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 + 3 ) + num_additional_dof
    );

    std::copy(
        std::begin( grad_velocity_tp1 ),
        std::end(   grad_velocity_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 ) + num_additional_dof
    );

    std::copy(
        std::begin( grad_theta_tp1 ),
        std::end(   grad_theta_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 + 9 ) + num_additional_dof
    );

    std::copy(
        std::begin( grad_e_tp1 ),
        std::end(   grad_e_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 + 9 + 3 ) + num_additional_dof
    );

    std::copy(
        std::begin( grad_z_tp1 ),
        std::end(   grad_z_tp1 ),
        std::begin( dof_vector ) + nphases * ( 1 + 3 + 3 + 1 + 1 + 3 + 9 + 9 + 3 + 3 ) + num_additional_dof
    );

    std::vector< floatType > previous_dof_vector( dof_vector.size( ) );

    std::array< floatType, node_count> Ns;
    e.GetShapeFunctions( xi_begin, xi_end, std::begin( Ns ), std::end( Ns ) );

    std::array< floatType, node_count * dim > dNdx;
    e.GetGlobalShapeFunctionGradients( xi_begin, xi_end, std::cbegin( x_tp1 ), std::cend( x_tp1 ),
                                       std::begin( dNdx ), std::end( dNdx ) );

    unsigned int low_bound  = 0;
    unsigned int high_bound = nphases;

    if ( active_phase >= 0 ){

        low_bound  = active_phase;
        high_bound = active_phase + 1;

    }

    // Assemble the material response
    std::array< floatType, nphases * material_response_size > material_response;
    std::fill( std::begin( material_response ), std::end( material_response ), 0 );

    std::array< floatType, nphases * material_response_size * dof_vector_size > material_response_jacobian;
    std::fill( std::begin( material_response_jacobian ), std::end( material_response_jacobian ), 0 );

    for ( unsigned int j = low_bound; j < high_bound; ++j ){

        hydraLinearTest linearTest(
            nphases, j, 9, num_additional_dof,
            0, 0.1, dof_vector, previous_dof_vector
        );

        linearTest.evaluate( );

        std::copy(
            std::begin( *linearTest.getUnknownVector( ) ),
            std::end(   *linearTest.getUnknownVector( ) ),
            std::begin( material_response ) + material_response_size * j
        );

        std::vector full_jacobian = linearTest.getFullTangent( );

        std::copy(
            std::begin( full_jacobian ),
            std::end(   full_jacobian ),
            std::begin( material_response_jacobian ) + material_response_size * dof_vector_size * j
        );

    }

    std::fill( value_begin, value_end, 0 );

    std::array< floatType, node_count * nphases * 1 > value_n;

    std::array< floatType, nphases * 1 * nphases * 1 > dRdRho_n, dRdTheta_n, dRdE_n, dRdVolumeFraction_n;

    std::array< floatType, nphases * 1 * nphases * 3 > dRdU_n;

    std::array< floatType, nphases * 1 * nphases * 3 > dRdW_n;

    std::array< floatType, nphases * 1 * num_additional_dof > dRdZ_n;

    std::array< floatType, nphases * 1 * 3 > dRdUMesh_n;

    std::fill( dRdRho_begin,   dRdRho_end,   0 );

    std::fill( dRdU_begin,     dRdU_end,     0 );

    std::fill( dRdW_begin,     dRdW_end,     0 );

    std::fill( dRdTheta_begin, dRdTheta_end, 0 );

    std::fill( dRdE_begin,     dRdE_end,     0 );

    std::fill( dRdZ_begin,     dRdZ_end,     0 );

    std::fill( dRdVF_begin,    dRdVF_end,    0 );

    std::fill( dRdUMesh_begin, dRdUMesh_end, 0 );

    std::fill( std::begin( value_n ),             std::end( value_n ), 0 );

    std::fill( std::begin( dRdRho_n ),            std::end( dRdRho_n ), 0 );

    std::fill( std::begin( dRdU_n ),              std::end( dRdU_n ), 0 );

    std::fill( std::begin( dRdW_n ),              std::end( dRdW_n ), 0 );

    std::fill( std::begin( dRdTheta_n ),          std::end( dRdTheta_n ), 0 );

    std::fill( std::begin( dRdE_n ),              std::end( dRdE_n ), 0 );

    std::fill( std::begin( dRdZ_n ),              std::end( dRdZ_n ), 0 );

    std::fill( std::begin( dRdVolumeFraction_n ), std::end( dRdVolumeFraction_n ), 0 );

    std::fill( std::begin( dRdUMesh_n ),          std::end( dRdUMesh_n ), 0 );

    for ( unsigned int i = 0; i < node_count; ++i ){

        if ( active_phase >= 0 ){

            unsigned int j = active_phase;

            tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<3,3,0,20,17,14,21>(
                density_tp1_p[ j ], density_dot_tp1_p[ j ],
                std::cbegin( grad_density_tp1 )  + 3 * j, std::cbegin( grad_density_tp1 )  + 3 * ( j + 1 ),
                e_tp1_p[ j ], e_dot_tp1_p[ j ],
                std::cbegin( grad_e_tp1 ) + 3 * j,        std::cbegin( grad_e_tp1 ) + 3 * ( j + 1 ),
                std::cbegin( v_tp1_p )           + 3 * j, std::cbegin( v_tp1_p )           + 3 * ( j + 1 ),
                std::cbegin( grad_velocity_tp1 ) + 9 * j, std::cbegin( grad_velocity_tp1 ) + 9 * ( j + 1 ),
                std::cbegin( material_response ) + material_response_size * j,
                std::cbegin( material_response ) + material_response_size * ( j + 1 ),
                vf_tp1_p[ j ],
                Ns[ i ], std::begin( dNdx ) + 3 * i, std::begin( dNdx ) + 3 * ( i + 1 ),
                *( value_begin + nphases * i + j )
            );

            std::transform(
                value_begin + nphases * i + j,
                value_begin + nphases * i + ( j + 1 ),
                value_begin + nphases * i + j,
                std::bind(
                    std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                    std::placeholders::_1,
                    J
                )
            );

        }
        else{

            tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<3,3,0,20,17,14,21>(
                std::cbegin( density_tp1_p ),     std::cend( density_tp1_p ),
                std::cbegin( density_dot_tp1_p ), std::cend( density_dot_tp1_p ),
                std::cbegin( grad_density_tp1 ),  std::cend( grad_density_tp1 ),
                std::cbegin( e_tp1_p ),           std::cend( e_tp1_p ),
                std::cbegin( e_dot_tp1_p ),       std::cend( e_dot_tp1_p ),
                std::cbegin( grad_e_tp1 ),        std::cend( grad_e_tp1 ),
                std::cbegin( v_tp1_p ),           std::cend( v_tp1_p ),
                std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
                std::cbegin( material_response ), std::cend( material_response ),
                std::cbegin( vf_tp1_p ),          std::cend( vf_tp1_p ),
                Ns[ i ], std::begin( dNdx ) + 3 * i, std::begin( dNdx ) + 3 * ( i + 1 ),
                value_begin + nphases * i,        value_begin + nphases * ( i + 1 )
            );

            std::transform(
                value_begin + nphases * i,
                value_begin + nphases * ( i + 1 ),
                value_begin + nphases * i,
                std::bind(
                    std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                    std::placeholders::_1,
                    J
                )
            );

        }

        for ( unsigned int k = 0; k < node_count; ++k ){

            if ( active_phase >= 0 ){

                unsigned int j = active_phase;

                tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<3,3,0,20,17,14,21,num_dof>(
                    density_tp1_p[ j ], density_dot_tp1_p[ j ],
                    std::cbegin( grad_density_tp1 )  + 3 * j, std::cbegin( grad_density_tp1 )  + 3 * ( j + 1 ),
                    e_tp1_p[ j ], e_dot_tp1_p[ j ],
                    std::cbegin( grad_e_tp1 )        + 3 * j, std::cbegin( grad_e_tp1 )        + 3 * ( j + 1 ),
                    std::cbegin( v_tp1_p )           + 3 * j, std::cbegin( v_tp1_p )           + 3 * ( j + 1 ),
                    std::cbegin( grad_velocity_tp1 ) + 9 * j, std::cbegin( grad_velocity_tp1 ) + 9 * ( j + 1 ),
                    std::cbegin( material_response ) + material_response_size * j,
                    std::cbegin( material_response ) + material_response_size * ( j + 1 ),
                    std::cbegin( material_response_jacobian ) + material_response_size * dof_vector_size * j,
                    std::cbegin( material_response_jacobian ) + material_response_size * dof_vector_size * ( j + 1 ),
                    vf_tp1_p[ j ],
                    Ns[ i ], std::begin( dNdx ) + 3 * i, std::begin( dNdx ) + 3 * ( i + 1 ),
                    Ns[ k ], std::begin( dNdx ) + 3 * k, std::begin( dNdx ) + 3 * ( k + 1 ),
                    std::cbegin( dof_vector ) + num_dof,       std::cend( dof_vector ),
                    dDensityDotdDensity, dEDotdE, dUDotdU,
                    j,
                    *( std::begin( value_n ) + nphases * i + j ),
                    std::begin( dRdRho_n )            + nphases * 1 * j,        std::begin( dRdRho_n )            + nphases * 1 * ( j + 1 ),
                    std::begin( dRdU_n )              + nphases * 3 * j,        std::begin( dRdU_n )              + nphases * 3 * ( j + 1 ),
                    std::begin( dRdW_n )              + nphases * 3 * j,        std::begin( dRdW_n )              + nphases * 3 * ( j + 1 ),
                    std::begin( dRdTheta_n )          + nphases * 1 * j,        std::begin( dRdTheta_n )          + nphases * 1 * ( j + 1 ),
                    std::begin( dRdE_n )              + nphases * 1 * j,        std::begin( dRdE_n )              + nphases * 1 * ( j + 1 ),
                    std::begin( dRdZ_n )              + num_additional_dof * j, std::begin( dRdZ_n )              + num_additional_dof * ( j + 1 ),
                    std::begin( dRdVolumeFraction_n ) + nphases * 1 * j,        std::begin( dRdVolumeFraction_n ) + nphases * 1 * ( j + 1 ),
                    std::begin( dRdUMesh_n )          + dim * j,                std::begin( dRdUMesh_n )          + dim * ( j + 1 )
                );

                std::transform(
                    std::begin( value_n ) + nphases * i + j,
                    std::begin( value_n ) + nphases * i + ( j + 1 ),
                    std::begin( value_n ) + nphases * i + j,
                    std::bind(
                        std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                        std::placeholders::_1,
                        J
                    )
                );

            }
            else{

                tardigradeBalanceEquations::balanceOfEnergy::computeBalanceOfEnergy<3,3,0,20,17,14,21,num_dof>(
                    std::cbegin( density_tp1_p ),     std::cend( density_tp1_p ),
                    std::cbegin( density_dot_tp1_p ), std::cend( density_dot_tp1_p ),
                    std::cbegin( grad_density_tp1 ),  std::cend( grad_density_tp1 ),
                    std::cbegin( e_tp1_p ),           std::cend( e_tp1_p ),
                    std::cbegin( e_dot_tp1_p ),       std::cend( e_dot_tp1_p ),
                    std::cbegin( grad_e_tp1 ),        std::cend( grad_e_tp1 ),
                    std::cbegin( v_tp1_p ),           std::cend( v_tp1_p ),
                    std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
                    std::cbegin( material_response ), std::cend( material_response ),
                    std::cbegin( material_response_jacobian ),
                    std::cend(   material_response_jacobian ),
                    std::cbegin( vf_tp1_p ),          std::cend( vf_tp1_p ),
                    Ns[ i ], std::begin( dNdx ) + 3 * i, std::begin( dNdx ) + 3 * ( i + 1 ),
                    Ns[ k ], std::begin( dNdx ) + 3 * k, std::begin( dNdx ) + 3 * ( k + 1 ),
                    std::cbegin( dof_vector ) + num_dof,       std::cend( dof_vector ),
                    dDensityDotdDensity, dEDotdE, dUDotdU,
                    std::begin( value_n ) + nphases * i, std::begin( value_n ) + nphases * ( i + 1 ),
                    std::begin( dRdRho_n )             , std::end( dRdRho_n )           ,
                    std::begin( dRdU_n )               , std::end( dRdU_n )             ,
                    std::begin( dRdW_n )               , std::end( dRdW_n )             ,
                    std::begin( dRdTheta_n )           , std::end( dRdTheta_n )         ,
                    std::begin( dRdE_n )               , std::end( dRdE_n )             ,
                    std::begin( dRdZ_n )               , std::end( dRdZ_n )             ,
                    std::begin( dRdVolumeFraction_n )  , std::end( dRdVolumeFraction_n ),
                    std::begin( dRdUMesh_n )           , std::end( dRdUMesh_n )         
                );

                std::transform(
                    std::begin( value_n ) + nphases * i,
                    std::begin( value_n ) + nphases * ( i + 1 ),
                    std::begin( value_n ) + nphases * i,
                    std::bind(
                        std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                        std::placeholders::_1,
                        J
                    )
                );

            }

           for ( unsigned int j = 0; j < nphases; ++j ){

                BOOST_TEST( value_n[ nphases * i + j ] == *( value_begin + nphases * i + j ) );

//                // node, phase, node, x
//                //    i,     j,    k, l

                for ( unsigned int l = 0; l < nphases; ++l ){

                    *( dRdRho_begin + nphases * node_count * nphases * 1 * i + node_count * nphases * 1 * j + nphases * 1 * k + l )
                        += dRdRho_n[ nphases * j + l ] * J;

                }

                for ( unsigned int l = 0; l < nphases * dim; ++l ){

                    *( dRdU_begin + nphases * node_count * nphases * dim * i + node_count * nphases * dim * j + nphases * dim * k + l )
                        += dRdU_n[ nphases * dim * j + l ] * J;

                }

                for ( unsigned int l = 0; l < nphases * dim; ++l ){

                    *( dRdW_begin + nphases * node_count * nphases * dim * i + node_count * nphases * dim * j + nphases * dim * k + l )
                        += dRdW_n[ nphases * dim * j + l ] * J;

                }

                for ( unsigned int l = 0; l < nphases; ++l ){

                    *( dRdTheta_begin + nphases * node_count * nphases * 1 * i + node_count * nphases * 1 * j + nphases * 1 * k + l )
                        += dRdTheta_n[ nphases * 1 * j + l ] * J;

                }

                for ( unsigned int l = 0; l < nphases; ++l ){

                    *( dRdE_begin + nphases * node_count * nphases * 1 * i + node_count * nphases * 1 * j + nphases * 1 * k + l )
                        += dRdE_n[ nphases * 1 * j + l ] * J;

                }

                for ( unsigned int l = 0; l < num_additional_dof; ++l ){

                    *( dRdZ_begin + nphases * node_count * num_additional_dof * 1 * i + node_count * num_additional_dof * 1 * j + num_additional_dof * 1 * k + l )
                        += dRdZ_n[ num_additional_dof * 1 * j + l ] * J;

                }

                for ( unsigned int l = 0; l < nphases; ++l ){

                    *( dRdVF_begin + nphases * node_count * nphases * 1 * i + node_count * nphases * 1 * j + nphases * 1 * k + l )
                        += dRdVolumeFraction_n[ nphases * 1 * j + l ] * J;

                }

                for ( unsigned int l = 0; l < dim; ++l ){

                    *( dRdUMesh_begin + nphases * node_count * 3 * i + node_count * 3 * j + 3 * k + l )
                        += dRdUMesh_n[ 3 * j + l ] * J;

                }

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_computeBalanceOfLinearMomentum_hydra_fea, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test computing the balance of energy in a finite element context
     */

    constexpr unsigned int nphases = 4;

    constexpr unsigned int active_phase = 2;

    constexpr unsigned int num_additional_dof = 5;

    std::array< floatType, 8 * nphases > density_t = {
        0.3213189 , 0.845533  , 0.18690375, 0.41729106, 0.98903451,
        0.23659981, 0.91683233, 0.91839747, 0.09129634, 0.46365272,
        0.50221634, 0.31366895, 0.04733954, 0.24168564, 0.09552964,
        0.23824991, 0.80779109, 0.89497829, 0.04322289, 0.30194684,
        0.9805822 , 0.53950482, 0.62630936, 0.00554541, 0.48490944,
        0.98832853, 0.37518553, 0.09703816, 0.46190876, 0.96300447,
        0.34183061, 0.79892273
    };

    std::array< floatType, 8 * nphases > density_tp1 = {
        0.01611863, 0.12695803, 0.77716246, 0.04589523, 0.71099869,
        0.97104614, 0.87168293, 0.71016165, 0.95850974, 0.42981334,
        0.87287891, 0.35595767, 0.92976365, 0.14877766, 0.94002901,
        0.8327162 , 0.84605484, 0.12392301, 0.5964869 , 0.01639248,
        0.72118437, 0.00773751, 0.08482228, 0.22549841, 0.87512453,
        0.36357632, 0.53995994, 0.56810321, 0.22546336, 0.57214677,
        0.6609518 , 0.29824539
    };

    std::array< floatType, 24 * nphases > u_t = {
         0.17498749,  0.89650474,  0.11206951,  0.00112284, -0.99293558,
        -0.03822191,  0.85491   , -0.60326862, -0.89581773, -0.18644221,
        -0.25520704,  0.71430612, -0.94677777,  0.84029846,  0.361806  ,
         0.80845199,  0.21505814,  0.62390662, -0.32891225, -0.30086754,
        -0.22025154,  0.50959416, -0.26141765, -0.51556039,  0.87533671,
         0.81602217, -0.30240537,  0.26927614, -0.45231558, -0.58776974,
        -0.32732094, -0.34580021,  0.7645522 ,  0.64460763,  0.41924646,
         0.91869045, -0.15491329, -0.50993392, -0.76520313, -0.39789328,
        -0.70947253, -0.81562781,  0.20586439, -0.2716251 ,  0.12914069,
        -0.61732856,  0.35381172, -0.56898911, -0.44395281,  0.48352084,
         0.11947579, -0.33032717,  0.08597757,  0.38796941,  0.82426424,
         0.16142643, -0.53462724,  0.49339526,  0.55553804, -0.59919737,
         0.64114844, -0.07013029,  0.55953332, -0.52504356, -0.33483946,
         0.90739424,  0.31563015,  0.54575566,  0.37674869, -0.59139176,
        -0.0586225 ,  0.61792775,  0.35007025, -0.98794423, -0.82518451,
        -0.30641056,  0.88873108, -0.01761904, -0.45964747, -0.27915256,
        -0.57869474, -0.15759989, -0.56392912,  0.69150501, -0.0874588 ,
        -0.44039596,  0.8657833 , -0.37129729,  0.81942932, -0.91316382,
         0.41423012, -0.03222192, -0.11155788, -0.92735331, -0.91863362,
        -0.33449277
    };

    std::array< floatType, 24 * nphases > u_tp1 = {
         0.89423908,  0.23531995, -0.26225032,  0.22395408, -0.58773693,
        -0.66986711, -0.27636547,  0.7267067 ,  0.01880345, -0.40619697,
         0.90050325,  0.63193218, -0.35405211,  0.94419649,  0.9747022 ,
        -0.18267973,  0.31184621, -0.1886936 , -0.48530379, -0.83469465,
        -0.47277931, -0.45704029, -0.20272184, -0.63022794,  0.90763681,
        -0.79424023,  0.25041707, -0.11660522, -0.1529639 , -0.25601643,
         0.73662942, -0.43904604, -0.95884769,  0.83619403,  0.72896056,
        -0.44619642,  0.0469751 , -0.78182361, -0.81314586,  0.67493222,
        -0.17946856,  0.32343308,  0.88640112, -0.50973882, -0.97368034,
        -0.95170319,  0.41877138,  0.84910377, -0.06533945, -0.2497817 ,
         0.08572085,  0.71783368,  0.30430775, -0.53404021,  0.54916041,
        -0.73077301, -0.66888006,  0.22536457, -0.52243319,  0.4095571 ,
        -0.30096295, -0.44515208,  0.99783681, -0.91876775,  0.29164504,
        -0.92260083,  0.52042052, -0.53982009, -0.82033627,  0.29689942,
         0.46520243,  0.35619063, -0.89619811, -0.41138611, -0.09782331,
        -0.42579342,  0.62102691, -0.73776979,  0.22435872,  0.97642989,
         0.80511308, -0.55568588, -0.99983622,  0.96119468,  0.76542597,
         0.83894493, -0.1689929 ,  0.48923092, -0.574337  , -0.21539186,
         0.7030961 , -0.74477555,  0.78773074, -0.00698406, -0.14780869,
        -0.38870722
    };

    std::array< floatType, 24 * nphases > w_t = {
         0.08336976,  0.00352469,  0.06080527,  0.07153036,  0.08447647,
        -0.03932385, -0.03203783,  0.01901478, -0.01173517,  0.08656851,
        -0.02048719, -0.00444439,  0.02343722, -0.0190521 ,  0.09849569,
        -0.08022974, -0.05587934, -0.03546897, -0.07045543, -0.04315615,
         0.05584906,  0.0045784 , -0.09320927,  0.09652452,  0.0232013 ,
        -0.0882121 ,  0.03223375, -0.02432613, -0.07286534,  0.01273292,
         0.04541599,  0.03422532, -0.05049737,  0.00497324,  0.00753269,
         0.04336067, -0.02802653,  0.05954652,  0.02558437, -0.09233368,
         0.0092958 ,  0.07238242,  0.01351483, -0.06483435,  0.00207527,
         0.05138917, -0.07797896,  0.06341982, -0.06650367,  0.0068153 ,
        -0.0228513 , -0.05027525,  0.0294865 , -0.09252158,  0.05200916,
         0.00538813,  0.07515424,  0.00414366, -0.09299337, -0.07127981,
         0.05912092, -0.00160479, -0.01162415, -0.03631304, -0.04309016,
         0.09317726, -0.01340613,  0.07680061,  0.02963262,  0.07168553,
         0.07048991,  0.09126241,  0.03958845,  0.06107939,  0.04662558,
         0.02104537,  0.04347083,  0.04315008, -0.09181844,  0.00322217,
         0.05853027, -0.05140756, -0.0069704 , -0.01300286, -0.01944257,
        -0.07563209,  0.00514231, -0.01075033,  0.03267855,  0.00988261,
        -0.09449141, -0.0936164 ,  0.04027196,  0.04151622,  0.09198783,
         0.07534094
    };

    std::array< floatType, 24 * nphases > w_tp1 = {
        -0.00638807,  0.0251813 , -0.00856365, -0.05541075, -0.0246646 ,
        -0.07922315,  0.03330542, -0.06159397, -0.00490644,  0.09348732,
        -0.09366621, -0.06965401, -0.04028416,  0.08836139,  0.08176836,
        -0.06759983,  0.09622355,  0.05014951,  0.00799542,  0.08634058,
         0.07612143, -0.0217367 ,  0.03126864,  0.02947703, -0.03460636,
        -0.06412197, -0.00663802, -0.04734379, -0.02898697,  0.09082879,
        -0.00777243,  0.03697829, -0.03275402,  0.09917222,  0.03175352,
        -0.06079811, -0.0803632 ,  0.08863611,  0.08895557,  0.02426568,
        -0.0966017 , -0.05489302,  0.06025536,  0.07509197, -0.00920204,
        -0.02689588, -0.045155  , -0.0766059 , -0.07685109,  0.09052054,
         0.06172522, -0.06704413, -0.05858999,  0.03111031,  0.05293284,
         0.06206297, -0.06733246,  0.09682566, -0.05443959,  0.01788309,
         0.01752315,  0.09347238,  0.03153349,  0.01698085,  0.00375452,
         0.05293151, -0.07878895, -0.09958162,  0.09049777, -0.00026846,
        -0.03433292, -0.02638935,  0.06076866, -0.02352596,  0.05403383,
        -0.0119076 ,  0.06881549, -0.08475919, -0.00377433, -0.00663006,
        -0.0471344 ,  0.08872295,  0.08100569, -0.01128074, -0.08056808,
        -0.05864337, -0.04570163, -0.00315605, -0.03232458,  0.05482721,
        -0.00479468,  0.0740741 ,  0.09915635, -0.05603281,  0.02233428,
         0.06950046
    };

    std::array< floatType, 8 * nphases > theta_t = {
         0.89047327, -0.41982715,  0.45408549, -0.9699677 ,  0.75828488,
        -0.8721229 ,  0.4667908 ,  0.98922077,  0.00237956, -0.58133202,
         0.18928717,  0.24829996,  0.33614547, -0.65477652,  0.79742538,
         0.24198274, -0.91286259,  0.36808213, -0.60783191, -0.94531844,
         0.10190655,  0.62662728,  0.7198823 , -0.79295815,  0.32608557,
         0.42015045, -0.41096601,  0.942728  , -0.44262506, -0.86003562,
         0.03856072,  0.38862977
    };

    std::array< floatType, 8 * nphases > theta_tp1 = {
        -0.51068043, -0.3228356 ,  0.12725597,  0.77335633,  0.49465181,
        -0.58081608, -0.49644586,  0.04776138,  0.53791739,  0.23752356,
         0.00264853,  0.19425068,  0.51212006,  0.07415959,  0.79550549,
         0.89413498,  0.83070901,  0.50903668, -0.50735799, -0.22945711,
        -0.4400001 ,  0.31532048, -0.35155675,  0.50878321, -0.77298184,
         0.55072953,  0.17180395,  0.67077737, -0.1382487 ,  0.24992891,
         0.10882426,  0.95134253
    };

    std::array< floatType, 8 * nphases > e_t = {
         0.51094878,  0.0896265 , -0.65193582,  0.80822843, -0.58832443,
         0.30008651,  0.87294371, -0.55284074, -0.54815292,  0.70363782,
         0.65531005, -0.29659329, -0.46980742, -0.74522306,  0.97587221,
         0.6706862 ,  0.79878321,  0.02735865, -0.77123034, -0.89483932,
        -0.33883579,  0.84066087,  0.89516367,  0.68232773, -0.68264171,
        -0.16015366, -0.50751416, -0.58930045,  0.3696517 , -0.02777666,
        -0.35018071, -0.79957108
    };

    std::array< floatType, 8 * nphases > e_tp1 = {
         0.08952674, -0.30594969, -0.21780839, -0.37898252, -0.22560959,
         0.11171917, -0.97171239,  0.69529399,  0.84383972,  0.10105938,
        -0.46395777,  0.980478  , -0.23361194,  0.38731079,  0.37990509,
        -0.13138187, -0.60168369,  0.93315875, -0.87261828, -0.02970122,
        -0.55853858, -0.41205174,  0.65705464, -0.26546887, -0.83330346,
        -0.60738199,  0.72074682,  0.9540577 , -0.46403568,  0.35081798,
        -0.83760201,  0.44693118
    };

    std::array< floatType, 8 * num_additional_dof > z_t = {
         0.34787193,  0.94019709,  0.4022445 ,  0.64344147, -0.90992083,
         0.34539703,  0.30950529, -0.7965079 ,  0.68477499,  0.22834481,
        -0.80334382,  0.18893424, -0.0431683 , -0.53341286, -0.96048782,
        -0.26886544,  0.23970215, -0.34144173, -0.38549069,  0.50224248,
         0.5172493 ,  0.43753167, -0.79763609,  0.03233191,  0.11559732,
         0.48960905,  0.80635544, -0.26192227, -0.14267306,  0.46553498,
         0.32527284,  0.1157398 , -0.29972074, -0.60929531, -0.63238526,
        -0.83683342, -0.8375983 ,  0.69159645, -0.2326545 , -0.87852075
    };

    std::array< floatType, 8 * num_additional_dof > z_tp1 = {
         0.79285134, -0.55345905, -0.46375114, -0.61100432,  0.93500212,
        -0.77491983,  0.44432648,  0.86417749,  0.33600259,  0.71745322,
        -0.5151058 ,  0.34785596,  0.40174268, -0.08333497,  0.74109124,
         0.3887722 ,  0.78975558,  0.50640871,  0.04058084, -0.00262357,
        -0.09254475, -0.95670627,  0.0702828 , -0.15405353, -0.68493281,
        -0.7618606 , -0.10129624, -0.92017389,  0.97315978, -0.24375814,
        -0.23578161, -0.8977474 , -0.14665533, -0.96850913, -0.93981273,
        -0.32180153,  0.6419379 , -0.08235784, -0.97031884, -0.67355994
    };

    std::array< floatType, 8 * nphases > vf_t = {
        0.73992272, 0.73829374, 0.75452291, 0.35166938, 0.35227695,
        0.80207567, 0.39813789, 0.72719101, 0.58112301, 0.36434168,
        0.08000652, 0.11612538, 0.88955872, 0.45234051, 0.99400454,
        0.36389695, 0.2499543 , 0.35053932, 0.3430861 , 0.63735673,
        0.01273756, 0.76326864, 0.41641463, 0.43223919, 0.48111502,
        0.44921245, 0.4974709 , 0.34590431, 0.45334614, 0.40465134,
        0.51824272, 0.62326908
    };

    std::array< floatType, 8 * nphases> vf_tp1 = {
        0.24104061, 0.50843715, 0.59462192, 0.01694831, 0.52049374,
        0.23929325, 0.40453855, 0.82653016, 0.3262356 , 0.48321692,
        0.02474115, 0.30875088, 0.63972109, 0.31516177, 0.20579751,
        0.29065568, 0.95437831, 0.08680182, 0.46335776, 0.0583869 ,
        0.53865824, 0.14603574, 0.63408481, 0.26439746, 0.69091539,
        0.34714606, 0.00416849, 0.29489471, 0.08189444, 0.49503963,
        0.28889007, 0.63999187
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

    std::array< floatType, 8 * nphases > density_dot_t = {
        0.79884633, 0.2082483 , 0.4433677 , 0.71560128, 0.41051979,
        0.19100696, 0.96749431, 0.65075037, 0.86545985, 0.02524236,
        0.26690581, 0.5020711 , 0.06744864, 0.99303326, 0.2364624 ,
        0.37429218, 0.21401191, 0.10544587, 0.23247979, 0.30061014,
        0.63444227, 0.28123478, 0.36227676, 0.00594284, 0.36571913,
        0.53388598, 0.16201584, 0.59743311, 0.29315247, 0.63205049,
        0.02619661, 0.88759346
    };

    std::array< floatType, 24 * nphases > u_dot_t = {
        -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225,
        -0.81579012, -0.13259765, -0.13827447, -0.0126298 , -0.14833942,
        -0.37547755, -0.14729739,  0.78677833,  0.88832004,  0.00367335,
         0.2479059 , -0.76876321, -0.36542904, -0.17034758,  0.73261832,
        -0.49908927, -0.03393147,  0.97111957,  0.03897024,  0.22578905,
        -0.75874267,  0.6526816 ,  0.20612026,  0.09013601, -0.31447233,
        -0.39175842, -0.16595558,  0.36260153,  0.75091368,  0.02084467,
         0.33862757,  0.17187311,  0.249807  ,  0.3493781 ,  0.68468488,
        -0.83361002,  0.52736568, -0.51266725, -0.61155408,  0.14491391,
        -0.80857497,  0.77065365,  0.25449794,  0.44683272, -0.96774159,
         0.18886376,  0.11357038, -0.68208071, -0.69385897,  0.39105906,
        -0.36246715,  0.38394059,  0.1087665 , -0.22209885,  0.85026498,
         0.68333999, -0.28520487, -0.91281707, -0.39046385, -0.20362864,
         0.40991766,  0.99071696, -0.28817027,  0.52509563,  0.18635383,
         0.3834036 , -0.6977451 , -0.20224741, -0.5182882 , -0.31308797,
         0.02625631,  0.3332491 , -0.78818303, -0.7382101 , -0.35603879,
         0.32312867,  0.69301245,  0.10651469,  0.70890498, -0.23032438,
        -0.36642421, -0.29147065, -0.65783634,  0.65822527, -0.32265831,
         0.10474015,  0.15710294,  0.04306612, -0.99462387,  0.97669084,
         0.81068315
    };

    std::array< floatType, 8 * nphases > e_dot_t = {
        -0.16712678,  0.83631984, -0.37692766,  0.88293399,  0.00649485,
        -0.30221415,  0.29403919, -0.50050759, -0.54047281, -0.60730711,
         0.91979913, -0.01417254,  0.50322993, -0.05201624,  0.17508032,
         0.16827796,  0.95977262,  0.33686625, -0.52046106, -0.96960468,
        -0.56263572, -0.08896072, -0.21315933,  0.62465248,  0.57111351,
        -0.82180806,  0.90402145,  0.05491335,  0.19280793, -0.18988646,
         0.29900191,  0.74265261
    };

    std::array< floatType, 24 > X = {
        0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0., 0., 1., 1., 0.,
       1., 1., 1., 1., 0., 1., 1.
    };

    std::array< floatType, 3 > local_point = {
        -0.16274628, -0.09382215,  0.86470132
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    floatType beta = 0.67;

    std::array< floatType, 8 * nphases > answer = {
        0.3129818 ,  0.30842469,  0.2894097 ,  0.32023897,  0.21284998,
        0.20692824,  0.19621173,  0.21475773,  0.17067598,  0.16461001,
        0.15740845,  0.17082534,  0.25142996,  0.24604353,  0.23261885,
        0.25545469,  0.03614033, -0.15154134, -0.25696741, -0.00796672,
       -0.14649378, -0.31802617, -0.35293074, -0.22398137, -0.19936641,
       -0.35746103, -0.36327156, -0.28328484, -0.07838785, -0.25607707,
       -0.31132838, -0.14232
    };

    std::array< floatType, 8 * nphases > result;

    evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
        std::cbegin( local_point ),   std::cend( local_point ),
        dt,
        std::cbegin( density_t ),     std::cend( density_t ),
        std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
        std::cbegin( u_t ),           std::cend( u_t ),
        std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
        std::cbegin( w_t ),           std::cend( w_t ),
        std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
        std::cbegin( theta_t ),       std::cend( theta_t ),
        std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
        std::cbegin( e_t ),           std::cend( e_t ),
        std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
        std::cbegin( z_t ),           std::cend( z_t ),
        std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
        std::cbegin( vf_t ),          std::cend( vf_t ),
        std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
        std::cbegin( umesh_t ),       std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
        std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
        std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
        std::cbegin( X ),             std::cend( X ),
        alpha, beta,
        std::begin( result ),  std::end( result ),
        active_phase
    );

    for ( unsigned int i = 0; i < 8; ++i ){

            BOOST_TEST( result[ nphases * i + active_phase ] == answer[ nphases * i + active_phase ] );

    }

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 1 > dRdRho;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 3 > dRdU;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 3 > dRdW;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 1 > dRdTheta;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 1 > dRdE;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * num_additional_dof > dRdZ;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 1 > dRdVF;

    std::array< floatType, 8 * 1 * nphases * 8 * 3 > dRdUMesh;

    std::fill( std::begin( result ), std::end( result ), 0 );

    evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
        std::cbegin( local_point ),   std::cend( local_point ),
        dt,
        std::cbegin( density_t ),     std::cend( density_t ),
        std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
        std::cbegin( u_t ),           std::cend( u_t ),
        std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
        std::cbegin( w_t ),           std::cend( w_t ),
        std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
        std::cbegin( theta_t ),       std::cend( theta_t ),
        std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
        std::cbegin( e_t ),           std::cend( e_t ),
        std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
        std::cbegin( z_t ),           std::cend( z_t ),
        std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
        std::cbegin( vf_t ),          std::cend( vf_t ),
        std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
        std::cbegin( umesh_t ),       std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
        std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
        std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
        std::cbegin( X ),             std::cend( X ),
        alpha, beta,
        std::begin( result ),   std::end( result ),
        std::begin(   dRdRho ), std::end(   dRdRho ),
        std::begin(     dRdU ), std::end(     dRdU ),
        std::begin(     dRdW ), std::end(     dRdW ),
        std::begin( dRdTheta ), std::end( dRdTheta ),
        std::begin(     dRdE ), std::end(     dRdE ),
        std::begin(     dRdZ ), std::end(     dRdZ ),
        std::begin(    dRdVF ), std::end(    dRdVF ),
        std::begin( dRdUMesh ), std::end( dRdUMesh ),
        active_phase
    );

    for ( unsigned int i = 0; i < 8; ++i ){

            BOOST_TEST( result[ nphases * i + active_phase ] == answer[ nphases * i + active_phase ] );

    }

    const floatType eps = 1e-5;

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

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( xp ),            std::cend( xp ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vp ),             std::end( vp ),
                active_phase
            );

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( xm ),            std::cend( xm ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vm ),             std::end( vm ),
                active_phase
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                const unsigned int node = j / ( nphases );
                const unsigned int phase = ( j - nphases * node );

                if ( phase == active_phase ){

                    BOOST_TEST( dRdRho[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

                }

            }

        }

    }

    // Check the derivatives w.r.t. the velocity
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

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( xp ),            std::cend( xp ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vp ),             std::end( vp ),
                active_phase
            );

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( xm ),            std::cend( xm ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vm ),             std::end( vm ),
                active_phase
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                const unsigned int node = j / ( nphases );
                const unsigned int phase = ( j - nphases * node );

                if ( phase == active_phase ){

                    BOOST_TEST( dRdU[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

                }

            }

        }

    }

    // Check the derivatives w.r.t. the displacement
    {

        constexpr unsigned int vardim = 3 * nphases * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( w_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = w_tp1;
            std::array< floatType, vardim > xm = w_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( xp ),            std::cend( xp ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vp ),             std::end( vp ),
                active_phase
            );

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( xm ),            std::cend( xm ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vm ),             std::end( vm ),
                active_phase
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                const unsigned int node = j / ( nphases );
                const unsigned int phase = ( j - nphases * node );

                if ( phase == active_phase ){

                    BOOST_TEST( dRdW[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

                }

            }

        }

    }

    // Check the derivatives w.r.t. the temperature
    {

        constexpr unsigned int vardim = 1 * nphases * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( theta_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = theta_tp1;
            std::array< floatType, vardim > xm = theta_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( xp ),            std::cend( xp ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vp ),             std::end( vp ),
                active_phase
            );

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( xm ),            std::cend( xm ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vm ),             std::end( vm ),
                active_phase
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                const unsigned int node = j / ( nphases );
                const unsigned int phase = ( j - nphases * node );

                if ( phase == active_phase ){

                    BOOST_TEST( dRdTheta[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

                }

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

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( xp ),            std::cend( xp ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vp ),             std::end( vp ),
                active_phase
            );

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( xm ),            std::cend( xm ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vm ),             std::end( vm ),
                active_phase
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                const unsigned int node = j / ( nphases );
                const unsigned int phase = ( j - nphases * node );

                if ( phase == active_phase ){

                    BOOST_TEST( dRdE[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

                }

            }

        }

    }

    // Check the derivatives w.r.t. the additional dof
    {

        constexpr unsigned int vardim = num_additional_dof * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( z_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = z_tp1;
            std::array< floatType, vardim > xm = z_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( xp ),            std::cend( xp ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vp ),             std::end( vp ),
                active_phase
            );

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( xm ),            std::cend( xm ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vm ),             std::end( vm ),
                active_phase
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                const unsigned int node = j / ( nphases );
                const unsigned int phase = ( j - nphases * node );

                if ( phase == active_phase ){

                    BOOST_TEST( dRdZ[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

                }

            }

        }

    }

    // Check the derivatives w.r.t. the volume fraction
    {

        constexpr unsigned int vardim = 1 * nphases * 8;
        constexpr unsigned int outdim = 1 * nphases * 8;

        for ( unsigned int i = 0; i < vardim; ++i ){

            floatType delta = eps * std::fabs( vf_tp1[ i ] ) + eps;

            std::array< floatType, vardim > xp = vf_tp1;
            std::array< floatType, vardim > xm = vf_tp1;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( xp ),            std::cend( xp ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vp ),             std::end( vp ),
                active_phase
            );

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( xm ),            std::cend( xm ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vm ),             std::end( vm ),
                active_phase
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                const unsigned int node = j / ( nphases );
                const unsigned int phase = ( j - nphases * node );

                if ( phase == active_phase ){

                    BOOST_TEST( dRdVF[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

                }

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

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( xp ),            std::cend( xp ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vp ),             std::end( vp ),
                active_phase
            );

            evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
                std::cbegin( local_point ),   std::cend( local_point ),
                dt,
                std::cbegin( density_t ),     std::cend( density_t ),
                std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
                std::cbegin( u_t ),           std::cend( u_t ),
                std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
                std::cbegin( w_t ),           std::cend( w_t ),
                std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
                std::cbegin( theta_t ),       std::cend( theta_t ),
                std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
                std::cbegin( e_t ),           std::cend( e_t ),
                std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
                std::cbegin( z_t ),           std::cend( z_t ),
                std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
                std::cbegin( vf_t ),          std::cend( vf_t ),
                std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
                std::cbegin( umesh_t ),       std::cend( umesh_t ),
                std::cbegin( xm ),            std::cend( xm ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
                std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
                std::cbegin( X ),             std::cend( X ),
                alpha, beta,
                std::begin( vm ),             std::end( vm ),
                active_phase
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                const unsigned int node = j / ( nphases );
                const unsigned int phase = ( j - nphases * node );

                if ( phase == active_phase ){

                    BOOST_TEST( dRdUMesh[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

                }

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_computeBalanceOfLinearMomentum_hydra_fea_multiphase, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test computing the balance of energy in a finite element context
     */

    constexpr unsigned int nphases = 4;

    constexpr unsigned int num_additional_dof = 5;

    std::array< floatType, 8 * nphases > density_t = {
        0.3213189 , 0.845533  , 0.18690375, 0.41729106, 0.98903451,
        0.23659981, 0.91683233, 0.91839747, 0.09129634, 0.46365272,
        0.50221634, 0.31366895, 0.04733954, 0.24168564, 0.09552964,
        0.23824991, 0.80779109, 0.89497829, 0.04322289, 0.30194684,
        0.9805822 , 0.53950482, 0.62630936, 0.00554541, 0.48490944,
        0.98832853, 0.37518553, 0.09703816, 0.46190876, 0.96300447,
        0.34183061, 0.79892273
    };

    std::array< floatType, 8 * nphases > density_tp1 = {
        0.01611863, 0.12695803, 0.77716246, 0.04589523, 0.71099869,
        0.97104614, 0.87168293, 0.71016165, 0.95850974, 0.42981334,
        0.87287891, 0.35595767, 0.92976365, 0.14877766, 0.94002901,
        0.8327162 , 0.84605484, 0.12392301, 0.5964869 , 0.01639248,
        0.72118437, 0.00773751, 0.08482228, 0.22549841, 0.87512453,
        0.36357632, 0.53995994, 0.56810321, 0.22546336, 0.57214677,
        0.6609518 , 0.29824539
    };

    std::array< floatType, 24 * nphases > u_t = {
         0.17498749,  0.89650474,  0.11206951,  0.00112284, -0.99293558,
        -0.03822191,  0.85491   , -0.60326862, -0.89581773, -0.18644221,
        -0.25520704,  0.71430612, -0.94677777,  0.84029846,  0.361806  ,
         0.80845199,  0.21505814,  0.62390662, -0.32891225, -0.30086754,
        -0.22025154,  0.50959416, -0.26141765, -0.51556039,  0.87533671,
         0.81602217, -0.30240537,  0.26927614, -0.45231558, -0.58776974,
        -0.32732094, -0.34580021,  0.7645522 ,  0.64460763,  0.41924646,
         0.91869045, -0.15491329, -0.50993392, -0.76520313, -0.39789328,
        -0.70947253, -0.81562781,  0.20586439, -0.2716251 ,  0.12914069,
        -0.61732856,  0.35381172, -0.56898911, -0.44395281,  0.48352084,
         0.11947579, -0.33032717,  0.08597757,  0.38796941,  0.82426424,
         0.16142643, -0.53462724,  0.49339526,  0.55553804, -0.59919737,
         0.64114844, -0.07013029,  0.55953332, -0.52504356, -0.33483946,
         0.90739424,  0.31563015,  0.54575566,  0.37674869, -0.59139176,
        -0.0586225 ,  0.61792775,  0.35007025, -0.98794423, -0.82518451,
        -0.30641056,  0.88873108, -0.01761904, -0.45964747, -0.27915256,
        -0.57869474, -0.15759989, -0.56392912,  0.69150501, -0.0874588 ,
        -0.44039596,  0.8657833 , -0.37129729,  0.81942932, -0.91316382,
         0.41423012, -0.03222192, -0.11155788, -0.92735331, -0.91863362,
        -0.33449277
    };

    std::array< floatType, 24 * nphases > u_tp1 = {
         0.89423908,  0.23531995, -0.26225032,  0.22395408, -0.58773693,
        -0.66986711, -0.27636547,  0.7267067 ,  0.01880345, -0.40619697,
         0.90050325,  0.63193218, -0.35405211,  0.94419649,  0.9747022 ,
        -0.18267973,  0.31184621, -0.1886936 , -0.48530379, -0.83469465,
        -0.47277931, -0.45704029, -0.20272184, -0.63022794,  0.90763681,
        -0.79424023,  0.25041707, -0.11660522, -0.1529639 , -0.25601643,
         0.73662942, -0.43904604, -0.95884769,  0.83619403,  0.72896056,
        -0.44619642,  0.0469751 , -0.78182361, -0.81314586,  0.67493222,
        -0.17946856,  0.32343308,  0.88640112, -0.50973882, -0.97368034,
        -0.95170319,  0.41877138,  0.84910377, -0.06533945, -0.2497817 ,
         0.08572085,  0.71783368,  0.30430775, -0.53404021,  0.54916041,
        -0.73077301, -0.66888006,  0.22536457, -0.52243319,  0.4095571 ,
        -0.30096295, -0.44515208,  0.99783681, -0.91876775,  0.29164504,
        -0.92260083,  0.52042052, -0.53982009, -0.82033627,  0.29689942,
         0.46520243,  0.35619063, -0.89619811, -0.41138611, -0.09782331,
        -0.42579342,  0.62102691, -0.73776979,  0.22435872,  0.97642989,
         0.80511308, -0.55568588, -0.99983622,  0.96119468,  0.76542597,
         0.83894493, -0.1689929 ,  0.48923092, -0.574337  , -0.21539186,
         0.7030961 , -0.74477555,  0.78773074, -0.00698406, -0.14780869,
        -0.38870722
    };

    std::array< floatType, 24 * nphases > w_t = {
         0.08336976,  0.00352469,  0.06080527,  0.07153036,  0.08447647,
        -0.03932385, -0.03203783,  0.01901478, -0.01173517,  0.08656851,
        -0.02048719, -0.00444439,  0.02343722, -0.0190521 ,  0.09849569,
        -0.08022974, -0.05587934, -0.03546897, -0.07045543, -0.04315615,
         0.05584906,  0.0045784 , -0.09320927,  0.09652452,  0.0232013 ,
        -0.0882121 ,  0.03223375, -0.02432613, -0.07286534,  0.01273292,
         0.04541599,  0.03422532, -0.05049737,  0.00497324,  0.00753269,
         0.04336067, -0.02802653,  0.05954652,  0.02558437, -0.09233368,
         0.0092958 ,  0.07238242,  0.01351483, -0.06483435,  0.00207527,
         0.05138917, -0.07797896,  0.06341982, -0.06650367,  0.0068153 ,
        -0.0228513 , -0.05027525,  0.0294865 , -0.09252158,  0.05200916,
         0.00538813,  0.07515424,  0.00414366, -0.09299337, -0.07127981,
         0.05912092, -0.00160479, -0.01162415, -0.03631304, -0.04309016,
         0.09317726, -0.01340613,  0.07680061,  0.02963262,  0.07168553,
         0.07048991,  0.09126241,  0.03958845,  0.06107939,  0.04662558,
         0.02104537,  0.04347083,  0.04315008, -0.09181844,  0.00322217,
         0.05853027, -0.05140756, -0.0069704 , -0.01300286, -0.01944257,
        -0.07563209,  0.00514231, -0.01075033,  0.03267855,  0.00988261,
        -0.09449141, -0.0936164 ,  0.04027196,  0.04151622,  0.09198783,
         0.07534094
    };

    std::array< floatType, 24 * nphases > w_tp1 = {
        -0.00638807,  0.0251813 , -0.00856365, -0.05541075, -0.0246646 ,
        -0.07922315,  0.03330542, -0.06159397, -0.00490644,  0.09348732,
        -0.09366621, -0.06965401, -0.04028416,  0.08836139,  0.08176836,
        -0.06759983,  0.09622355,  0.05014951,  0.00799542,  0.08634058,
         0.07612143, -0.0217367 ,  0.03126864,  0.02947703, -0.03460636,
        -0.06412197, -0.00663802, -0.04734379, -0.02898697,  0.09082879,
        -0.00777243,  0.03697829, -0.03275402,  0.09917222,  0.03175352,
        -0.06079811, -0.0803632 ,  0.08863611,  0.08895557,  0.02426568,
        -0.0966017 , -0.05489302,  0.06025536,  0.07509197, -0.00920204,
        -0.02689588, -0.045155  , -0.0766059 , -0.07685109,  0.09052054,
         0.06172522, -0.06704413, -0.05858999,  0.03111031,  0.05293284,
         0.06206297, -0.06733246,  0.09682566, -0.05443959,  0.01788309,
         0.01752315,  0.09347238,  0.03153349,  0.01698085,  0.00375452,
         0.05293151, -0.07878895, -0.09958162,  0.09049777, -0.00026846,
        -0.03433292, -0.02638935,  0.06076866, -0.02352596,  0.05403383,
        -0.0119076 ,  0.06881549, -0.08475919, -0.00377433, -0.00663006,
        -0.0471344 ,  0.08872295,  0.08100569, -0.01128074, -0.08056808,
        -0.05864337, -0.04570163, -0.00315605, -0.03232458,  0.05482721,
        -0.00479468,  0.0740741 ,  0.09915635, -0.05603281,  0.02233428,
         0.06950046
    };

    std::array< floatType, 8 * nphases > theta_t = {
         0.89047327, -0.41982715,  0.45408549, -0.9699677 ,  0.75828488,
        -0.8721229 ,  0.4667908 ,  0.98922077,  0.00237956, -0.58133202,
         0.18928717,  0.24829996,  0.33614547, -0.65477652,  0.79742538,
         0.24198274, -0.91286259,  0.36808213, -0.60783191, -0.94531844,
         0.10190655,  0.62662728,  0.7198823 , -0.79295815,  0.32608557,
         0.42015045, -0.41096601,  0.942728  , -0.44262506, -0.86003562,
         0.03856072,  0.38862977
    };

    std::array< floatType, 8 * nphases > theta_tp1 = {
        -0.51068043, -0.3228356 ,  0.12725597,  0.77335633,  0.49465181,
        -0.58081608, -0.49644586,  0.04776138,  0.53791739,  0.23752356,
         0.00264853,  0.19425068,  0.51212006,  0.07415959,  0.79550549,
         0.89413498,  0.83070901,  0.50903668, -0.50735799, -0.22945711,
        -0.4400001 ,  0.31532048, -0.35155675,  0.50878321, -0.77298184,
         0.55072953,  0.17180395,  0.67077737, -0.1382487 ,  0.24992891,
         0.10882426,  0.95134253
    };

    std::array< floatType, 8 * nphases > e_t = {
         0.51094878,  0.0896265 , -0.65193582,  0.80822843, -0.58832443,
         0.30008651,  0.87294371, -0.55284074, -0.54815292,  0.70363782,
         0.65531005, -0.29659329, -0.46980742, -0.74522306,  0.97587221,
         0.6706862 ,  0.79878321,  0.02735865, -0.77123034, -0.89483932,
        -0.33883579,  0.84066087,  0.89516367,  0.68232773, -0.68264171,
        -0.16015366, -0.50751416, -0.58930045,  0.3696517 , -0.02777666,
        -0.35018071, -0.79957108
    };

    std::array< floatType, 8 * nphases > e_tp1 = {
         0.08952674, -0.30594969, -0.21780839, -0.37898252, -0.22560959,
         0.11171917, -0.97171239,  0.69529399,  0.84383972,  0.10105938,
        -0.46395777,  0.980478  , -0.23361194,  0.38731079,  0.37990509,
        -0.13138187, -0.60168369,  0.93315875, -0.87261828, -0.02970122,
        -0.55853858, -0.41205174,  0.65705464, -0.26546887, -0.83330346,
        -0.60738199,  0.72074682,  0.9540577 , -0.46403568,  0.35081798,
        -0.83760201,  0.44693118
    };

    std::array< floatType, 8 * num_additional_dof > z_t = {
         0.34787193,  0.94019709,  0.4022445 ,  0.64344147, -0.90992083,
         0.34539703,  0.30950529, -0.7965079 ,  0.68477499,  0.22834481,
        -0.80334382,  0.18893424, -0.0431683 , -0.53341286, -0.96048782,
        -0.26886544,  0.23970215, -0.34144173, -0.38549069,  0.50224248,
         0.5172493 ,  0.43753167, -0.79763609,  0.03233191,  0.11559732,
         0.48960905,  0.80635544, -0.26192227, -0.14267306,  0.46553498,
         0.32527284,  0.1157398 , -0.29972074, -0.60929531, -0.63238526,
        -0.83683342, -0.8375983 ,  0.69159645, -0.2326545 , -0.87852075
    };

    std::array< floatType, 8 * num_additional_dof > z_tp1 = {
         0.79285134, -0.55345905, -0.46375114, -0.61100432,  0.93500212,
        -0.77491983,  0.44432648,  0.86417749,  0.33600259,  0.71745322,
        -0.5151058 ,  0.34785596,  0.40174268, -0.08333497,  0.74109124,
         0.3887722 ,  0.78975558,  0.50640871,  0.04058084, -0.00262357,
        -0.09254475, -0.95670627,  0.0702828 , -0.15405353, -0.68493281,
        -0.7618606 , -0.10129624, -0.92017389,  0.97315978, -0.24375814,
        -0.23578161, -0.8977474 , -0.14665533, -0.96850913, -0.93981273,
        -0.32180153,  0.6419379 , -0.08235784, -0.97031884, -0.67355994
    };

    std::array< floatType, 8 * nphases > vf_t = {
        0.73992272, 0.73829374, 0.75452291, 0.35166938, 0.35227695,
        0.80207567, 0.39813789, 0.72719101, 0.58112301, 0.36434168,
        0.08000652, 0.11612538, 0.88955872, 0.45234051, 0.99400454,
        0.36389695, 0.2499543 , 0.35053932, 0.3430861 , 0.63735673,
        0.01273756, 0.76326864, 0.41641463, 0.43223919, 0.48111502,
        0.44921245, 0.4974709 , 0.34590431, 0.45334614, 0.40465134,
        0.51824272, 0.62326908
    };

    std::array< floatType, 8 * nphases> vf_tp1 = {
        0.24104061, 0.50843715, 0.59462192, 0.01694831, 0.52049374,
        0.23929325, 0.40453855, 0.82653016, 0.3262356 , 0.48321692,
        0.02474115, 0.30875088, 0.63972109, 0.31516177, 0.20579751,
        0.29065568, 0.95437831, 0.08680182, 0.46335776, 0.0583869 ,
        0.53865824, 0.14603574, 0.63408481, 0.26439746, 0.69091539,
        0.34714606, 0.00416849, 0.29489471, 0.08189444, 0.49503963,
        0.28889007, 0.63999187
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

    std::array< floatType, 8 * nphases > density_dot_t = {
        0.79884633, 0.2082483 , 0.4433677 , 0.71560128, 0.41051979,
        0.19100696, 0.96749431, 0.65075037, 0.86545985, 0.02524236,
        0.26690581, 0.5020711 , 0.06744864, 0.99303326, 0.2364624 ,
        0.37429218, 0.21401191, 0.10544587, 0.23247979, 0.30061014,
        0.63444227, 0.28123478, 0.36227676, 0.00594284, 0.36571913,
        0.53388598, 0.16201584, 0.59743311, 0.29315247, 0.63205049,
        0.02619661, 0.88759346
    };

    std::array< floatType, 24 * nphases > u_dot_t = {
        -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225,
        -0.81579012, -0.13259765, -0.13827447, -0.0126298 , -0.14833942,
        -0.37547755, -0.14729739,  0.78677833,  0.88832004,  0.00367335,
         0.2479059 , -0.76876321, -0.36542904, -0.17034758,  0.73261832,
        -0.49908927, -0.03393147,  0.97111957,  0.03897024,  0.22578905,
        -0.75874267,  0.6526816 ,  0.20612026,  0.09013601, -0.31447233,
        -0.39175842, -0.16595558,  0.36260153,  0.75091368,  0.02084467,
         0.33862757,  0.17187311,  0.249807  ,  0.3493781 ,  0.68468488,
        -0.83361002,  0.52736568, -0.51266725, -0.61155408,  0.14491391,
        -0.80857497,  0.77065365,  0.25449794,  0.44683272, -0.96774159,
         0.18886376,  0.11357038, -0.68208071, -0.69385897,  0.39105906,
        -0.36246715,  0.38394059,  0.1087665 , -0.22209885,  0.85026498,
         0.68333999, -0.28520487, -0.91281707, -0.39046385, -0.20362864,
         0.40991766,  0.99071696, -0.28817027,  0.52509563,  0.18635383,
         0.3834036 , -0.6977451 , -0.20224741, -0.5182882 , -0.31308797,
         0.02625631,  0.3332491 , -0.78818303, -0.7382101 , -0.35603879,
         0.32312867,  0.69301245,  0.10651469,  0.70890498, -0.23032438,
        -0.36642421, -0.29147065, -0.65783634,  0.65822527, -0.32265831,
         0.10474015,  0.15710294,  0.04306612, -0.99462387,  0.97669084,
         0.81068315
    };

    std::array< floatType, 8 * nphases > e_dot_t = {
        -0.16712678,  0.83631984, -0.37692766,  0.88293399,  0.00649485,
        -0.30221415,  0.29403919, -0.50050759, -0.54047281, -0.60730711,
         0.91979913, -0.01417254,  0.50322993, -0.05201624,  0.17508032,
         0.16827796,  0.95977262,  0.33686625, -0.52046106, -0.96960468,
        -0.56263572, -0.08896072, -0.21315933,  0.62465248,  0.57111351,
        -0.82180806,  0.90402145,  0.05491335,  0.19280793, -0.18988646,
         0.29900191,  0.74265261
    };

    std::array< floatType, 24 > X = {
        0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0., 0., 1., 1., 0.,
       1., 1., 1., 1., 0., 1., 1.
    };

    std::array< floatType, 3 > local_point = {
        -0.16274628, -0.09382215,  0.86470132
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    floatType beta = 0.67;

    std::array< floatType, 8 * nphases > answer = {
        0.3129818 ,  0.30842469,  0.2894097 ,  0.32023897,  0.21284998,
        0.20692824,  0.19621173,  0.21475773,  0.17067598,  0.16461001,
        0.15740845,  0.17082534,  0.25142996,  0.24604353,  0.23261885,
        0.25545469,  0.03614033, -0.15154134, -0.25696741, -0.00796672,
       -0.14649378, -0.31802617, -0.35293074, -0.22398137, -0.19936641,
       -0.35746103, -0.36327156, -0.28328484, -0.07838785, -0.25607707,
       -0.31132838, -0.14232
    };

    std::array< floatType, 8 * nphases > result;

    std::cout << "multiphase\n";

    evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
        std::cbegin( local_point ),   std::cend( local_point ),
        dt,
        std::cbegin( density_t ),     std::cend( density_t ),
        std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
        std::cbegin( u_t ),           std::cend( u_t ),
        std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
        std::cbegin( w_t ),           std::cend( w_t ),
        std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
        std::cbegin( theta_t ),       std::cend( theta_t ),
        std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
        std::cbegin( e_t ),           std::cend( e_t ),
        std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
        std::cbegin( z_t ),           std::cend( z_t ),
        std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
        std::cbegin( vf_t ),          std::cend( vf_t ),
        std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
        std::cbegin( umesh_t ),       std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
        std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
        std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
        std::cbegin( X ),             std::cend( X ),
        alpha, beta,
        std::begin( result ),  std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 1 > dRdRho;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 3 > dRdU;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 3 > dRdW;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 1 > dRdTheta;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 1 > dRdE;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * num_additional_dof > dRdZ;

    std::array< floatType, 8 * 1 * nphases * 8 * nphases * 1 > dRdVF;

    std::array< floatType, 8 * 1 * nphases * 8 * 3 > dRdUMesh;

    std::fill( std::begin( result ), std::end( result ), 0 );

    evaluate_at_nodes< 3, 8, nphases, num_additional_dof >(
        std::cbegin( local_point ),   std::cend( local_point ),
        dt,
        std::cbegin( density_t ),     std::cend( density_t ),
        std::cbegin( density_tp1 ),   std::cend( density_tp1 ),
        std::cbegin( u_t ),           std::cend( u_t ),
        std::cbegin( u_tp1 ),         std::cend( u_tp1 ),
        std::cbegin( w_t ),           std::cend( w_t ),
        std::cbegin( w_tp1 ),         std::cend( w_tp1 ),
        std::cbegin( theta_t ),       std::cend( theta_t ),
        std::cbegin( theta_tp1 ),     std::cend( theta_tp1 ),
        std::cbegin( e_t ),           std::cend( e_t ),
        std::cbegin( e_tp1 ),         std::cend( e_tp1 ),
        std::cbegin( z_t ),           std::cend( z_t ),
        std::cbegin( z_tp1 ),         std::cend( z_tp1 ),
        std::cbegin( vf_t ),          std::cend( vf_t ),
        std::cbegin( vf_tp1 ),        std::cend( vf_tp1 ),
        std::cbegin( umesh_t ),       std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),     std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
        std::cbegin( u_dot_t ),       std::cend( u_dot_t ),
        std::cbegin( e_dot_t ),       std::cend( e_dot_t ),
        std::cbegin( X ),             std::cend( X ),
        alpha, beta,
        std::begin( result ),   std::end( result ),
        std::begin(   dRdRho ), std::end(   dRdRho ),
        std::begin(     dRdU ), std::end(     dRdU ),
        std::begin(     dRdW ), std::end(     dRdW ),
        std::begin( dRdTheta ), std::end( dRdTheta ),
        std::begin(     dRdE ), std::end(     dRdE ),
        std::begin(     dRdZ ), std::end(     dRdZ ),
        std::begin(    dRdVF ), std::end(    dRdVF ),
        std::begin( dRdUMesh ), std::end( dRdUMesh )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

}
