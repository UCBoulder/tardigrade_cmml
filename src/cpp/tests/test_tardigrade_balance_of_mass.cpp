/**
  * \file test_tardigrade_balance_equations_balance_of_mass.cpp
  *
  * Tests for tardigrade_balance_equations_balance_of_mass
  */

#include<tardigrade_balance_of_mass.h>
#include<tardigrade_finite_element_utilities.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define USE_EIGEN
#include<tardigrade_vector_tools.h>

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

    result = 0.;

    tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, density_gradient, velocity, velocity_gradient, result,
                                                                     dCdRho,  dCdRhoDot,   dCdGradRho      , dCdV,     dCdGradV );

    BOOST_TEST( answer == result );

    floatType eps = 1e-6;

    // Derivative w.r.t. density
    for ( unsigned int i = 0; i < 1; ++i ){

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
    for ( unsigned int i = 0; i < 1; ++i ){

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
    for ( unsigned int i = 0; i < 3; ++i ){

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
    for ( unsigned int i = 0; i < 3; ++i ){

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
    for ( unsigned int i = 0; i < 9; ++i ){

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

template<
  typename dt_type, class v_t_in, class v_tp1_in,
  class vDot_t_in, typename alpha_type, class vDot_tp1_out
>
void compute_current_rate_of_change(
    const dt_type &dt,
    const v_t_in &v_t_begin, const v_t_in &v_t_end,
    const v_tp1_in &v_tp1_begin, const v_tp1_in &v_tp1_end,
    const vDot_t_in &vDot_t_begin, const vDot_t_in &vDot_t_end,
    alpha_type alpha, vDot_tp1_out vDot_tp1_begin, vDot_tp1_out vDot_tp1_end
){

    for ( unsigned int i = 0; i < ( unsigned int )( v_t_end - v_t_begin ); ++i ){

        *( vDot_tp1_begin + i ) = ( ( *( v_tp1_begin + i ) ) - ( *( v_t_begin + i ) ) ) / ( alpha * dt ) - ( ( 1 - alpha ) / alpha ) * ( *( vDot_t_begin + i ) );

    }

}

template<
    int dim, int node_count, int nphases,
    class xi_in, typename dt_type, class density_t_in, class density_tp1_in,
    class u_t_in, class u_tp1_in, class umesh_t_in, class umesh_tp1_in, class density_dot_t_in, class v_t_in,
    class X_in, typename alpha_type, class value_out
>
void evaluate_at_nodes(
    const xi_in &xi_begin, const xi_in &xi_end, dt_type dt, const density_t_in &density_t_begin,
    const density_t_in &density_t_end, const density_tp1_in &density_tp1_begin,
    const density_tp1_in &density_tp1_end,
    const u_t_in &u_t_begin, const u_t_in &u_t_end,
    const u_tp1_in &u_tp1_begin, const u_tp1_in &u_tp1_end,
    const umesh_t_in &umesh_t_begin, const umesh_t_in &umesh_t_end,
    const umesh_tp1_in &umesh_tp1_begin, const umesh_tp1_in &umesh_tp1_end,
    const density_dot_t_in &density_dot_t_begin, const density_dot_t_in &density_dot_t_end,
    const v_t_in &v_t_begin, const v_t_in &v_t_end, const X_in &X_begin, const X_in &X_end,
    const alpha_type &alpha, value_out value_begin, value_out value_end
){

    // Update the mesh nodes
    std::array< typename std::iterator_traits<umesh_t_in  >::value_type, dim * node_count > x_t;
    std::array< typename std::iterator_traits<umesh_tp1_in>::value_type, dim * node_count > x_tp1;

    std::transform( X_begin, X_end,   umesh_t_begin,   std::begin( x_t ), std::plus<typename std::iterator_traits<umesh_t_in  >::value_type>( ) );
    std::transform( X_begin, X_end, umesh_tp1_begin, std::begin( x_tp1 ), std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>( ) );

    // Calculate the current rates of change
    std::array< typename std::iterator_traits<density_tp1_in>::value_type, node_count * nphases > density_dot_tp1;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count * nphases > v_tp1;

    compute_current_rate_of_change(
        dt, density_t_begin, density_t_end, density_tp1_begin, density_tp1_end,
        density_dot_t_begin, density_dot_t_end, alpha,
        std::begin( density_dot_tp1 ), std::end( density_dot_tp1 )
    );

    compute_current_rate_of_change(
        dt, u_t_begin, u_t_end, u_tp1_begin, u_tp1_end,
        v_t_begin, v_t_end, alpha,
        std::begin( v_tp1 ), std::end( v_tp1 )
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
        typename std::iterator_traits<u_tp1_in>::value_type, dim * nphases
    > v_tp1_p;

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

    // Compute the gradients at the local point

    std::array<
        typename std::iterator_traits<density_tp1_in>::value_type, dim * nphases
    > grad_density_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * dim * nphases
    > grad_velocity_tp1;

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( grad_density_tp1 ), std::end( grad_density_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, std::cbegin( v_tp1 ), std::cend( v_tp1 ),
        std::begin( grad_velocity_tp1 ), std::end( grad_velocity_tp1 )
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

    if ( nphases == 1 ){

        for ( unsigned int i = 0; i < node_count; ++i ){

            tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>(
                density_tp1_p[ 0 ], density_dot_tp1_p[ 0 ],
                std::cbegin( grad_density_tp1 ), std::cend( grad_density_tp1 ),
                std::cbegin( v_tp1_p ), std::cend( v_tp1_p ),
                std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
                Ns[ i ], *( value_begin + i )
            );

            *( value_begin + i ) *= J;

        }

    }
    else{

        for ( unsigned int i = 0; i < node_count; ++i ){

            tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>(
                std::cbegin( density_tp1_p ),     std::cend( density_tp1_p ),
                std::cbegin( density_dot_tp1_p ), std::cend( density_dot_tp1_p ),
                std::cbegin( grad_density_tp1 ),  std::cend( grad_density_tp1 ),
                std::cbegin( v_tp1_p ),           std::cend( v_tp1_p ),
                std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
                Ns[ i ],
                value_begin + nphases * i,        value_begin + nphases * ( i + 1 )
            );

            std::transform(
                value_begin + nphases * i, value_begin + nphases * ( i + 1 ),
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
    int dim, int node_count, int nphases,
    class xi_in, typename dt_type, class density_t_in, class density_tp1_in,
    class u_t_in, class u_tp1_in, class umesh_t_in, class umesh_tp1_in,
    class density_dot_t_in, class v_t_in,
    class X_in, typename alpha_type, class value_out,
    class dCdRho_iter_out, class dCdU_iter_out, class dCdUMesh_iter_out
>
void evaluate_at_nodes(
    const xi_in &xi_begin, const xi_in &xi_end, dt_type dt, const density_t_in &density_t_begin,
    const density_t_in &density_t_end, const density_tp1_in &density_tp1_begin,
    const density_tp1_in &density_tp1_end,
    const u_t_in &u_t_begin, const u_t_in &u_t_end,
    const u_tp1_in &u_tp1_begin, const u_tp1_in &u_tp1_end,
    const umesh_t_in &umesh_t_begin, const umesh_t_in &umesh_t_end,
    const umesh_tp1_in &umesh_tp1_begin, const umesh_tp1_in &umesh_tp1_end,
    const density_dot_t_in &density_dot_t_begin, const density_dot_t_in &density_dot_t_end,
    const v_t_in &v_t_begin, const v_t_in &v_t_end, const X_in &X_begin, const X_in &X_end,
    const alpha_type &alpha, value_out value_begin, value_out value_end,
    dCdRho_iter_out dCdRho_begin,     dCdRho_iter_out dCdRho_end,
    dCdU_iter_out dCdU_begin,         dCdU_iter_out dCdU_end,
    dCdUMesh_iter_out dCdUMesh_begin, dCdUMesh_iter_out dCdUMesh_end
){

    // Update the mesh nodes
    std::array< typename std::iterator_traits<umesh_t_in  >::value_type, dim * node_count > x_t;
    std::array< typename std::iterator_traits<umesh_tp1_in>::value_type, dim * node_count > x_tp1;

    std::transform( X_begin, X_end,   umesh_t_begin,   std::begin( x_t ), std::plus<typename std::iterator_traits<umesh_t_in  >::value_type>( ) );
    std::transform( X_begin, X_end, umesh_tp1_begin, std::begin( x_tp1 ), std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>( ) );

    // Calculate the current rates of change
    std::array< typename std::iterator_traits<density_tp1_in>::value_type, node_count * nphases > density_dot_tp1;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count * nphases > v_tp1;

    compute_current_rate_of_change(
        dt, density_t_begin, density_t_end, density_tp1_begin, density_tp1_end,
        density_dot_t_begin, density_dot_t_end, alpha,
        std::begin( density_dot_tp1 ), std::end( density_dot_tp1 )
    );

    auto dDensityDotdDensity = 1. / ( alpha * dt );

    compute_current_rate_of_change(
        dt, u_t_begin, u_t_end, u_tp1_begin, u_tp1_end,
        v_t_begin, v_t_end, alpha,
        std::begin( v_tp1 ), std::end( v_tp1 )
    );

    auto dUDotdU = 1. / ( alpha * dt );

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
        typename std::iterator_traits<u_tp1_in>::value_type, dim * nphases
    > v_tp1_p;

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

    // Compute the gradients at the local point

    std::array<
        typename std::iterator_traits<density_tp1_in>::value_type, dim * nphases
    > grad_density_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * dim * nphases
    > grad_velocity_tp1;

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( grad_density_tp1 ), std::end( grad_density_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, std::cbegin( v_tp1 ), std::cend( v_tp1 ),
        std::begin( grad_velocity_tp1 ), std::end( grad_velocity_tp1 )
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

    std::array< floatType, dim * node_count> dNdxs;
    e.GetGlobalShapeFunctionGradients( xi_begin, xi_end, std::begin( x_tp1 ), std::end( x_tp1 ), std::begin( dNdxs ), std::end( dNdxs ) );

    std::fill( dCdRho_begin,   dCdRho_end,   0 );
    std::fill( dCdU_begin,     dCdU_end,     0 );
    std::fill( dCdUMesh_begin, dCdUMesh_end, 0 );

    if ( nphases == 1 ){

        typename std::iterator_traits<value_out>::value_type balance_of_mass;

        typename std::iterator_traits<dCdRho_iter_out>::value_type dCdRho_p;
        std::array< typename std::iterator_traits<dCdU_iter_out>::value_type, dim > dCdU_p;
        std::array< typename std::iterator_traits<dCdUMesh_iter_out>::value_type, dim > dCdUMesh_p;
    
        for ( unsigned int i = 0; i < node_count; ++i ){
    
            tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>(
                density_tp1_p[ 0 ], density_dot_tp1_p[ 0 ],
                std::cbegin( grad_density_tp1 ), std::cend( grad_density_tp1 ),
                std::cbegin( v_tp1_p ), std::cend( v_tp1_p ),
                std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
                Ns[ i ], *( value_begin + i )
            );
    
            *( value_begin + i ) *= J;
    
        }
    
        for ( unsigned int i = 0; i < node_count; ++i ){ //Loop over the test functions
    
            for ( unsigned int j = 0; j < node_count; ++j ){ //Loop over the interpolation functions
    
                tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>(
                    density_tp1_p[ 0 ], density_dot_tp1_p[ 0 ],
                    std::cbegin( grad_density_tp1 ), std::cend( grad_density_tp1 ),
                    std::cbegin( v_tp1_p ), std::cend( v_tp1_p ),
                    std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
                    Ns[ i ], Ns[ j ],
                    std::cbegin( dNdxs ) + dim * j, std::cbegin( dNdxs ) + dim * ( j + 1 ),
                    dDensityDotdDensity, dUDotdU,
                    balance_of_mass,
                    dCdRho_p,
                    std::begin( dCdU_p ), std::end( dCdU_p ),
                    std::begin( dCdUMesh_p ), std::end( dCdUMesh_p )
                );
    
                BOOST_TEST( balance_of_mass * J == ( *( value_begin + i ) ) );
    
                *( dCdRho_begin + node_count * i + j ) = dCdRho_p * J;
    
                std::transform(
                    std::begin( dCdU_p ), std::end( dCdU_p ), dCdU_begin + node_count * dim * i + dim * j,
                    std::bind(
                        std::multiplies< typename std::iterator_traits< dCdU_iter_out >::value_type >( ),
                        std::placeholders::_1, J
                    )
                );
    
                std::transform(
                    std::begin( dCdUMesh_p ), std::end( dCdUMesh_p ), dCdUMesh_begin + node_count * dim * i + dim * j,
                    std::bind(
                        std::multiplies< typename std::iterator_traits< dCdUMesh_iter_out >::value_type >( ),
                        std::placeholders::_1, J
                    )
                );
    
            }
    
        }

    }
    else{

        std::array< typename std::iterator_traits<value_out>::value_type, nphases > balance_of_mass;

        std::array< typename std::iterator_traits<dCdRho_iter_out>::value_type, nphases > dCdRho_p;
        std::array< typename std::iterator_traits<dCdU_iter_out>::value_type, dim * nphases > dCdU_p;
        std::array< typename std::iterator_traits<dCdUMesh_iter_out>::value_type, dim * nphases > dCdUMesh_p;

        std::array< dt_type, nphases > dDensityDotdDensity_array;
        std::array< dt_type, nphases > dUDotdU_array;
    
        for ( unsigned int n = 0; n < nphases; ++n ){

            dDensityDotdDensity_array[ n ] = dDensityDotdDensity;
            dUDotdU_array[ n ]             = dUDotdU;

        }

        for ( unsigned int i = 0; i < node_count; ++i ){

            tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>(
                std::cbegin( density_tp1_p ),     std::cend( density_tp1_p ),
                std::cbegin( density_dot_tp1_p ), std::cend( density_dot_tp1_p ),
                std::cbegin( grad_density_tp1 ),  std::cend( grad_density_tp1 ),
                std::cbegin( v_tp1_p ),           std::cend( v_tp1_p ),
                std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
                Ns[ i ],
                value_begin + nphases * i,        value_begin + nphases * ( i + 1 )
            );

            std::transform(
                value_begin + nphases * i, value_begin + nphases * ( i + 1 ),
                value_begin + nphases * i,
                std::bind(
                    std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                    std::placeholders::_1,
                    J
                )

            );
    
        }
    
        for ( unsigned int i = 0; i < node_count; ++i ){ //Loop over the test functions
    
            for ( unsigned int j = 0; j < node_count; ++j ){ //Loop over the interpolation functions
    
                tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>(
                    std::cbegin( density_tp1_p ),     std::cend( density_tp1_p ),
                    std::cbegin( density_dot_tp1_p ), std::cend( density_dot_tp1_p ),
                    std::cbegin( grad_density_tp1 ), std::cend( grad_density_tp1 ),
                    std::cbegin( v_tp1_p ), std::cend( v_tp1_p ),
                    std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
                    Ns[ i ], Ns[ j ],
                    std::cbegin( dNdxs ) + dim * j, std::cbegin( dNdxs ) + dim * ( j + 1 ),
                    std::cbegin( dDensityDotdDensity_array ), std::cend( dDensityDotdDensity_array ),
                    std::cbegin( dUDotdU_array ),             std::cend( dUDotdU_array ),
                    std::begin( balance_of_mass ),            std::end( balance_of_mass ),
                    std::begin( dCdRho_p ),                   std::end( dCdRho_p ),
                    std::begin( dCdU_p ),                     std::end( dCdU_p ),
                    std::begin( dCdUMesh_p ),                 std::end( dCdUMesh_p )
                );

                for ( unsigned int k = 0; k < nphases; ++k ){
    
                    BOOST_TEST( balance_of_mass[ k ] * J == ( *( value_begin + nphases * i + k ) ) );

                    *( dCdRho_begin + nphases * node_count * nphases * i + node_count * nphases * k + nphases * j + k )
                        = dCdRho_p[ k ] * J;

                    for ( unsigned int l = 0; l < dim; ++l ){

                        *( dCdU_begin + nphases * node_count * nphases * dim * i + node_count * nphases * dim * k + nphases * dim * j + dim * k + l )
                            = dCdU_p[ dim * k + l ] * J;

                        *( dCdUMesh_begin + nphases * node_count * dim * i + node_count * dim * k + dim * j + l )
                            = dCdUMesh_p[ dim * k + l ] * J;

                    }

                }

            }
    
        }

    }

}

BOOST_AUTO_TEST_CASE( test_computeBalanceOfMass_fea, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test computing the balance of mass in a finite element context
     */

    std::array< floatType, 8 > density_t = {
        0.61289453, 0.12062867, 0.8263408 , 0.60306013, 0.54506801,
        0.34276383, 0.30412079, 0.41702221
    };

    std::array< floatType, 8 > density_tp1 = {
        0.08319499, 0.76368284, 0.24366637, 0.19422296, 0.57245696,
        0.09571252, 0.88532683, 0.62724897
    };

    std::array< floatType, 24 > u_t = {
        0.11357038, -0.68208071, -0.69385897,  0.39105906, -0.36246715,
        0.38394059,  0.1087665 , -0.22209885,  0.85026498,  0.68333999,
       -0.28520487, -0.91281707, -0.39046385, -0.20362864,  0.40991766,
        0.99071696, -0.28817027,  0.52509563,  0.18635383,  0.3834036 ,
       -0.6977451 , -0.20224741, -0.5182882 , -0.31308797
    };

    std::array< floatType, 24 > u_tp1 = {
        0.02625631,  0.3332491 , -0.78818303, -0.7382101 , -0.35603879,
        0.32312867,  0.69301245,  0.10651469,  0.70890498, -0.23032438,
       -0.36642421, -0.29147065, -0.65783634,  0.65822527, -0.32265831,
        0.10474015,  0.15710294,  0.04306612, -0.99462387,  0.97669084,
        0.81068315, -0.58472828, -0.41502117,  0.04002031
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

    std::array< floatType, 24 > v_t = {
        -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225,
        -0.81579012, -0.13259765, -0.13827447, -0.0126298 , -0.14833942,
        -0.37547755, -0.14729739,  0.78677833,  0.88832004,  0.00367335,
         0.2479059 , -0.76876321, -0.36542904, -0.17034758,  0.73261832,
        -0.49908927, -0.03393147,  0.97111957,  0.03897024
    };

    std::array< floatType, 24 > X = {
        0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0., 0., 1., 1., 0.,
       1., 1., 1., 1., 0., 1., 1.
    };

    std::array< floatType, 8 > answer = {
        -0.01576343, -0.04122993, -0.00067591, -0.00025842, -0.0231041 ,
       -0.06042976, -0.00099066, -0.00037876
    };

    std::array< floatType, 8 > result;

    std::array< floatType, 3 > local_point = {
        0.44683272, -0.96774159,  0.18886376
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    evaluate_at_nodes<3, 8, 1>(
        std::cbegin( local_point ), std::cend( local_point ), dt,
        std::cbegin( density_t ), std::cend( density_t ),
        std::cbegin( density_tp1 ), std::cend( density_tp1 ),
        std::cbegin( u_t ), std::cend( u_t ),
        std::cbegin( u_tp1 ), std::cend( u_tp1 ),
        std::cbegin( umesh_t ), std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
        std::cbegin( v_t ), std::cend( v_t ),
        std::cbegin( X ), std::cend( X ),
        alpha,
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0 );

    std::array< floatType, 8 * 1 * 8 > dCdRho;

    std::array< floatType, 8 * 3 * 8 > dCdU, dCdUMesh;

    evaluate_at_nodes<3, 8, 1>(
        std::cbegin( local_point ), std::cend( local_point ), dt,
        std::cbegin( density_t ), std::cend( density_t ),
        std::cbegin( density_tp1 ), std::cend( density_tp1 ),
        std::cbegin( u_t ), std::cend( u_t ),
        std::cbegin( u_tp1 ), std::cend( u_tp1 ),
        std::cbegin( umesh_t ), std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
        std::cbegin( v_t ), std::cend( v_t ),
        std::cbegin( X ), std::cend( X ),
        alpha,
        std::begin( result ), std::end( result ),
        std::begin( dCdRho ), std::end( dCdRho ),
        std::begin( dCdU ), std::end( dCdU ),
        std::begin( dCdUMesh ), std::end( dCdUMesh )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    // Check the derivatives w.r.t. the density
    {

        constexpr unsigned int vardim = 1 * 8;
        constexpr unsigned int outdim = 8;

        for ( unsigned int i = 0; i < vardim; ++i ){
    
            floatType delta = eps * std::fabs( density_tp1[ i ] ) + eps;
    
            std::array< floatType, vardim > xp = density_tp1;
            std::array< floatType, vardim > xm = density_tp1;
    
            xp[ i ] += delta;
            xm[ i ] -= delta;
    
            std::array< floatType, outdim > vp, vm;
    
            evaluate_at_nodes<3, 8, 1>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( xp ), std::cend( xp ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( xm ), std::cend( xm ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dCdRho[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the deformation
    {

        constexpr unsigned int vardim = 3 * 8;
        constexpr unsigned int outdim = 8;

        for ( unsigned int i = 0; i < vardim; ++i ){
    
            floatType delta = eps * std::fabs( u_tp1[ i ] ) + eps;
    
            std::array< floatType, vardim > xp = u_tp1;
            std::array< floatType, vardim > xm = u_tp1;
    
            xp[ i ] += delta;
            xm[ i ] -= delta;
    
            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( xp ), std::cend( xp ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( xm ), std::cend( xm ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dCdU[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the mesh deformation
    {

        constexpr unsigned int vardim = 3 * 8;
        constexpr unsigned int outdim = 8;

        for ( unsigned int i = 0; i < vardim; ++i ){
    
            floatType delta = eps * std::fabs( umesh_tp1[ i ] ) + eps;
    
            std::array< floatType, vardim > xp = umesh_tp1;
            std::array< floatType, vardim > xm = umesh_tp1;
    
            xp[ i ] += delta;
            xm[ i ] -= delta;
    
            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, 1>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( xp ), std::cend( xp ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, 1>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( xm ), std::cend( xm ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dCdUMesh[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_computeBalanceOfMass_multiphase_fea, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test computing the multiphase balance of mass in a finite element context
     */

    constexpr int nphases = 4;

    std::array< floatType, 8 * nphases > density_t = {
        0.20763586, 0.29248941, 0.52001015, 0.90191137, 0.98363088,
        0.25754206, 0.56435904, 0.80696868, 0.39437005, 0.73107304,
        0.16106901, 0.60069857, 0.86586446, 0.98352161, 0.07936579,
        0.42834727, 0.20454286, 0.45063649, 0.54776357, 0.09332671,
        0.29686078, 0.92758424, 0.56900373, 0.457412  , 0.75352599,
        0.74186215, 0.04857903, 0.7086974 , 0.83924335, 0.16593788,
        0.78099794, 0.28653662
    };

    std::array< floatType, 8 * nphases > density_tp1 = {
        0.80793821, 0.00742638, 0.55159273, 0.93193215, 0.58217546,
        0.20609573, 0.71775756, 0.37898585, 0.66838395, 0.02931972,
        0.63590036, 0.03219793, 0.74478066, 0.472913  , 0.12175436,
        0.54263593, 0.06677444, 0.65336487, 0.99608633, 0.76939734,
        0.57377411, 0.10263526, 0.69983407, 0.66116787, 0.04909713,
        0.7922993 , 0.51871659, 0.42586769, 0.78818717, 0.41156922,
        0.48102628, 0.18162884
    };

    std::array< floatType, 24 * nphases > u_t = {
        -0.16541788,  0.97806901, -0.52680038,  0.83366467,  0.83679494,
        -0.81740732, -0.07269455,  0.00443267, -0.3726621 , -0.90532093,
        -0.51662873, -0.80894072, -0.52350019,  0.61558217,  0.78995658,
        -0.91355422, -0.39610633,  0.9611644 ,  0.07900965,  0.25261872,
        -0.98890918, -0.03018111,  0.97665707, -0.24962895, -0.80592368,
        -0.07618248,  0.92600893, -0.31633877,  0.59784547,  0.59769266,
        -0.58350341, -0.1132646 ,  0.43120255, -0.17896043, -0.61798609,
         0.93498861,  0.30150073,  0.7309197 , -0.94951528, -0.46618837,
         0.0041422 , -0.86510273,  0.98606652, -0.52707521, -0.25141564,
        -0.57197617, -0.78910827, -0.53504043, -0.39877973,  0.26888454,
        -0.43753044, -0.27544648, -0.98811431, -0.26856175,  0.06777196,
        -0.67596833,  0.19486622, -0.41369506,  0.26410099, -0.94760679,
         0.77518692, -0.96776274, -0.74608394,  0.55432492, -0.90820954,
         0.42199739,  0.94209228,  0.74336587,  0.4203233 ,  0.91701949,
        -0.14037332,  0.74575783, -0.28808466,  0.85952731, -0.70244469,
         0.88005803,  0.66543239,  0.69210968, -0.75215398,  0.1929738 ,
        -0.96721504,  0.44236873, -0.98452497, -0.83035545, -0.54900318,
         0.75024907, -0.27284736,  0.07991987,  0.13620643, -0.54907328,
         0.14429354,  0.32190359, -0.40350921, -0.16274628, -0.09382215,
         0.86470132
    };

    std::array< floatType, 24 * nphases > u_tp1 = {
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
        0.30646975, 0.66526147, 0.11139217, 0.66487245, 0.88785679,
        0.69631127, 0.44032788, 0.43821438, 0.7650961 , 0.565642  ,
        0.08490416, 0.58267109, 0.8148437 , 0.33706638, 0.92757658,
        0.750717  , 0.57406383, 0.75164399, 0.07914896, 0.85938908,
        0.82150411, 0.90987166, 0.1286312 , 0.08178009, 0.13841557,
        0.39937871, 0.42430686, 0.56221838, 0.12224355, 0.2013995 ,
        0.81164435, 0.46798757
    };

    std::array< floatType, 24 * nphases > v_t = {
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

    std::array< floatType, 24 > X = {
        0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0., 0., 1., 1., 0.,
       1., 1., 1., 1., 0., 1., 1.
    };

    std::array< floatType, 8 * nphases > answer = {
        0.01068099, -0.00525696,  0.00240617, -0.01230476,  0.00505687,
       -0.00248888,  0.00113919, -0.00582564,  0.02768069, -0.01362384,
        0.00623578, -0.03188882,  0.0584664 , -0.02877592,  0.01317105,
       -0.06735471,  0.0024552 , -0.0012084 ,  0.0005531 , -0.00282845,
        0.00116241, -0.00057211,  0.00026186, -0.00133912,  0.00636287,
       -0.00313167,  0.0014334 , -0.00733018,  0.01343948, -0.00661462,
        0.00302759, -0.01548261
    };

    std::array< floatType, 8 * nphases > result;

    std::array< floatType, 3 > local_point = {
        -0.3573622 ,  0.69106599, -0.6261925
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    evaluate_at_nodes<3, 8, nphases>(
        std::cbegin( local_point ), std::cend( local_point ), dt,
        std::cbegin( density_t ), std::cend( density_t ),
        std::cbegin( density_tp1 ), std::cend( density_tp1 ),
        std::cbegin( u_t ), std::cend( u_t ),
        std::cbegin( u_tp1 ), std::cend( u_tp1 ),
        std::cbegin( umesh_t ), std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
        std::cbegin( v_t ), std::cend( v_t ),
        std::cbegin( X ), std::cend( X ),
        alpha,
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0 );

    std::array< floatType, 8 * 1 * nphases * 8 * 1 * nphases > dCdRho;

    std::array< floatType, 8 * 1 * nphases * 8 * 3 * nphases > dCdU;

    std::array< floatType, 8 * 1 * nphases * 8 * 3 > dCdUMesh;

    evaluate_at_nodes<3, 8, nphases>(
        std::cbegin( local_point ), std::cend( local_point ), dt,
        std::cbegin( density_t ), std::cend( density_t ),
        std::cbegin( density_tp1 ), std::cend( density_tp1 ),
        std::cbegin( u_t ), std::cend( u_t ),
        std::cbegin( u_tp1 ), std::cend( u_tp1 ),
        std::cbegin( umesh_t ), std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
        std::cbegin( v_t ), std::cend( v_t ),
        std::cbegin( X ), std::cend( X ),
        alpha,
        std::begin( result ), std::end( result ),
        std::begin( dCdRho ), std::end( dCdRho ),
        std::begin( dCdU ), std::end( dCdU ),
        std::begin( dCdUMesh ), std::end( dCdUMesh )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    // Check the derivatives w.r.t. the density
    {

        constexpr unsigned int vardim = 1 * 8 * nphases;
        constexpr unsigned int outdim = 8 * nphases;

        for ( unsigned int i = 0; i < vardim; ++i ){
    
            floatType delta = eps * std::fabs( density_tp1[ i ] ) + eps;
    
            std::array< floatType, vardim > xp = density_tp1;
            std::array< floatType, vardim > xm = density_tp1;
    
            xp[ i ] += delta;
            xm[ i ] -= delta;
    
            std::array< floatType, outdim > vp, vm;
    
            evaluate_at_nodes<3, 8, nphases>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( xp ), std::cend( xp ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( xm ), std::cend( xm ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dCdRho[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the deformation
    {

        constexpr unsigned int vardim = 3 * 8 * nphases;
        constexpr unsigned int outdim = 8 * nphases;

        for ( unsigned int i = 0; i < vardim; ++i ){
    
            floatType delta = eps * std::fabs( u_tp1[ i ] ) + eps;
    
            std::array< floatType, vardim > xp = u_tp1;
            std::array< floatType, vardim > xm = u_tp1;
    
            xp[ i ] += delta;
            xm[ i ] -= delta;
    
            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( xp ), std::cend( xp ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( xm ), std::cend( xm ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dCdU[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // Check the derivatives w.r.t. the mesh deformation
    {

        constexpr unsigned int vardim = 3 * 8;
        constexpr unsigned int outdim = 8 * nphases;

        for ( unsigned int i = 0; i < vardim; ++i ){
    
            floatType delta = eps * std::fabs( umesh_tp1[ i ] ) + eps;
    
            std::array< floatType, vardim > xp = umesh_tp1;
            std::array< floatType, vardim > xm = umesh_tp1;
    
            xp[ i ] += delta;
            xm[ i ] -= delta;
    
            std::array< floatType, outdim > vp, vm;

            evaluate_at_nodes<3, 8, nphases>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( xp ), std::cend( xp ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8, nphases>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
                std::cbegin( umesh_t ), std::cend( umesh_t ),
                std::cbegin( xm ), std::cend( xm ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; ++j ){

                BOOST_TEST( dCdUMesh[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

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

    tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                          std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                          std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( result ), std::end( result ) );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    std::array<floatType,nphases> dCdRho, dCdRhoDot;

    std::array<floatType,nphases*dim> dCdGradRho, dCdV;

    std::array<floatType,nphases*dim*dim> dCdGradV;

    tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
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
    for ( unsigned int i = 0; i < nphases; ++i ){

        floatType delta = eps * std::fabs( density[ i ] ) + eps;

        std::array<floatType,nphases> xp = density;

        std::array<floatType,nphases> xm = density;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( xp ), std::end( xp ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( xm ), std::end( xm ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; ++j ){

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
    for ( unsigned int i = 0; i < nphases; ++i ){

        floatType delta = eps * std::fabs( density_dot[ i ] ) + eps;

        std::array<floatType,nphases> xp = density_dot;

        std::array<floatType,nphases> xm = density_dot;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( xp ), std::end( xp ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( xm ), std::end( xm ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; ++j ){

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
    for ( unsigned int i = 0; i < dim * nphases; ++i ){

        floatType delta = eps * std::fabs( density_gradient[ i ] ) + eps;

        std::array<floatType,nphases*dim> xp = density_gradient;

        std::array<floatType,nphases*dim> xm = density_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( xp ), std::end( xp ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( xm ), std::end( xm ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; ++j ){

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
    for ( unsigned int i = 0; i < dim * nphases; ++i ){

        floatType delta = eps * std::fabs( velocity[ i ] ) + eps;

        std::array<floatType,nphases*dim> xp = velocity;

        std::array<floatType,nphases*dim> xm = velocity;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( xp ), std::end( xp ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( xm ), std::end( xm ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; ++j ){

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
    for ( unsigned int i = 0; i < dim * dim * nphases; ++i ){

        floatType delta = eps * std::fabs( velocity_gradient[ i ] ) + eps;

        std::array<floatType,nphases*dim*dim> xp = velocity_gradient;

        std::array<floatType,nphases*dim*dim> xm = velocity_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( xp ), std::end( xp ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( xm ), std::end( xm ), std::begin( vm ), std::end( vm ) );

        for ( unsigned int j = 0; j < nphases; ++j ){

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
