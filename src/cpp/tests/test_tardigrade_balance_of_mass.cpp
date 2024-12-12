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
    int dim, int node_count,
    class xi_in, typename dt_type, class density_t_in, class density_tp1_in,
    class u_t_in, class u_tp1_in, class density_dot_t_in, class v_t_in,
    class X_in, typename alpha_type, class value_out
>
void evaluate_at_nodes(
    const xi_in &xi_begin, const xi_in &xi_end, dt_type dt, const density_t_in &density_t_begin,
    const density_t_in &density_t_end, const density_tp1_in &density_tp1_begin,
    const density_tp1_in &density_tp1_end, const u_t_in &u_t_begin, const u_t_in &u_t_end,
    const u_tp1_in &u_tp1_begin, const u_tp1_in &u_tp1_end,
    const density_dot_t_in &density_dot_t_begin, const density_dot_t_in &density_dot_t_end,
    const v_t_in &v_t_begin, const v_t_in &v_t_end, const X_in &X_begin, const X_in &X_end,
    const alpha_type &alpha, value_out value_begin, value_out value_end
){

    // Update the mesh nodes
    std::array< typename std::iterator_traits<u_t_in  >::value_type, dim * node_count > x_t;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count > x_tp1;

    std::transform( X_begin, X_end,   u_t_begin,   std::begin( x_t ), std::plus<typename std::iterator_traits<u_t_in  >::value_type>( ) );
    std::transform( X_begin, X_end, u_tp1_begin, std::begin( x_tp1 ), std::plus<typename std::iterator_traits<u_tp1_in>::value_type>( ) );

    // Calculate the current rates of change
    std::array< typename std::iterator_traits<density_tp1_in>::value_type, node_count > density_dot_tp1;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count > v_tp1;

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
        typename std::iterator_traits<density_tp1_in>::value_type, 1
    > density_tp1_p, density_dot_tp1_p;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim
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
        typename std::iterator_traits<density_tp1_in>::value_type, dim
    > grad_density_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * dim
    > grad_velocity_tp1;

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( grad_density_tp1 ), std::end( grad_density_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, std::cbegin( v_tp1 ), std::cend( v_tp1 ),
        std::begin( grad_velocity_tp1 ), std::end( grad_velocity_tp1 )
    );

    typename std::iterator_traits<value_out>::value_type balance_of_mass;

    tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>(
        density_tp1_p[ 0 ], density_dot_tp1_p[ 0 ],
        std::cbegin( grad_density_tp1 ), std::cend( grad_density_tp1 ),
        std::cbegin( v_tp1_p ), std::cend( v_tp1_p ),
        std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
        balance_of_mass
    );

    e.GetShapeFunctions( xi_begin, xi_end, value_begin, value_end );

    std::transform( value_begin, value_end, value_begin, std::bind( std::multiplies<typename std::iterator_traits<value_out>::value_type>( ), std::placeholders::_1, balance_of_mass ) );

}

template<
    int dim, int node_count,
    class xi_in, typename dt_type, class density_t_in, class density_tp1_in,
    class u_t_in, class u_tp1_in, class density_dot_t_in, class v_t_in,
    class X_in, typename alpha_type, class value_out,
    class dCdRho_iter_out, class dCdRhoDot_iter_out, class dCdGradRho_iter_out,
    class dCdV_iter_out, class dCdGradV_iter_out, class dCdU_iter_out
>
void evaluate_at_nodes(
    const xi_in &xi_begin, const xi_in &xi_end, dt_type dt, const density_t_in &density_t_begin,
    const density_t_in &density_t_end, const density_tp1_in &density_tp1_begin,
    const density_tp1_in &density_tp1_end, const u_t_in &u_t_begin, const u_t_in &u_t_end,
    const u_tp1_in &u_tp1_begin, const u_tp1_in &u_tp1_end,
    const density_dot_t_in &density_dot_t_begin, const density_dot_t_in &density_dot_t_end,
    const v_t_in &v_t_begin, const v_t_in &v_t_end, const X_in &X_begin, const X_in &X_end,
    const alpha_type &alpha, value_out value_begin, value_out value_end,
    dCdRho_iter_out dCdRho_begin, dCdRho_iter_out dCdRho_end,
    dCdRhoDot_iter_out dCdRhoDot_begin, dCdRhoDot_iter_out dCdRhoDot_end,
    dCdGradRho_iter_out dCdGradRho_begin, dCdGradRho_iter_out dCdGradRho_end,
    dCdV_iter_out dCdV_begin,             dCdV_iter_out dCdV_end,
    dCdGradV_iter_out dCdGradV_begin,     dCdGradV_iter_out dCdGradV_end,
    dCdU_iter_out dCdU_begin,             dCdU_iter_out dCdU_end
){

    // Update the mesh nodes
    std::array< typename std::iterator_traits<u_t_in  >::value_type, dim * node_count > x_t;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count > x_tp1;

    std::transform( X_begin, X_end,   u_t_begin,   std::begin( x_t ), std::plus<typename std::iterator_traits<u_t_in  >::value_type>( ) );
    std::transform( X_begin, X_end, u_tp1_begin, std::begin( x_tp1 ), std::plus<typename std::iterator_traits<u_tp1_in>::value_type>( ) );

    // Calculate the current rates of change
    std::array< typename std::iterator_traits<density_tp1_in>::value_type, node_count > density_dot_tp1;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count > v_tp1;

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
        typename std::iterator_traits<density_tp1_in>::value_type, 1
    > density_tp1_p, density_dot_tp1_p;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim
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
        typename std::iterator_traits<density_tp1_in>::value_type, dim
    > grad_density_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * dim
    > grad_velocity_tp1;

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( grad_density_tp1 ), std::end( grad_density_tp1 )
    );

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, std::cbegin( v_tp1 ), std::cend( v_tp1 ),
        std::begin( grad_velocity_tp1 ), std::end( grad_velocity_tp1 )
    );

    typename std::iterator_traits<value_out>::value_type balance_of_mass;

    typename std::iterator_traits<dCdRho_iter_out>::value_type dCdRho_p;
    typename std::iterator_traits<dCdRhoDot_iter_out>::value_type dCdRhoDot_p;
    std::array< typename std::iterator_traits<dCdGradRho_iter_out>::value_type, dim > dCdGradRho_p;
    std::array< typename std::iterator_traits<dCdV_iter_out>::value_type, dim > dCdV_p;
    std::array< typename std::iterator_traits<dCdGradV_iter_out>::value_type, dim * dim > dCdGradV_p;

    tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>(
        density_tp1_p[ 0 ], density_dot_tp1_p[ 0 ],
        std::cbegin( grad_density_tp1 ), std::cend( grad_density_tp1 ),
        std::cbegin( v_tp1_p ), std::cend( v_tp1_p ),
        std::cbegin( grad_velocity_tp1 ), std::cend( grad_velocity_tp1 ),
        balance_of_mass,
        dCdRho_p, dCdRhoDot_p,
        std::begin( dCdGradRho_p ), std::end( dCdGradRho_p ),
        std::begin( dCdV_p ), std::end( dCdV_p ),
        std::begin( dCdGradV_p ), std::end( dCdGradV_p )
    );

    std::array< floatType, node_count> Ns;
    e.GetShapeFunctions( xi_begin, xi_end, std::begin( Ns ), std::end( Ns ) );

    std::transform( std::begin( Ns ), std::end( Ns ), value_begin,     std::bind( std::multiplies<typename std::iterator_traits<value_out>::value_type>( ), std::placeholders::_1, balance_of_mass ) );

    for ( unsigned int i = 0; i < node_count; ++i ){ //Loop over the test functions

        for ( unsigned int j = 0; j < node_count; ++j ){ //Loop over the interpolation functions

            *( dCdRho_begin + dim * i + j ) = dCdRho_p * Ns[ i ] * Ns[ j ];

            *( dCdRhoDot_begin + dim * i + j ) = dCdRho_p * Ns[ i ] * Ns[ j ];

//            std::transform(
//                std::begin( dCdGradRho_p ), std::end( dCdGradRho_p ), dCdGradRho_begin + node_count * dim * i + dim * j,
//                std::bind(
//                    std::multiplies<
//                        typename std::iterator_traits<dCdGradRho_iter_out>::value_type
//                    >( ),
//                    std::placeholders::_1,
//                    Ns[ i ] * Ns[ j ]
//                )
//            );
//
//            std::transform(
//                std::begin( dCdV_p ), std::end( dCdV_p ), dCdV_begin + node_count * dim * i + dim * j,
//                std::bind(
//                    std::multiplies<
//                        typename std::iterator_traits<dCdV_iter_out>::value_type
//                    >( ),
//                    std::placeholders::_1,
//                    Ns[ i ] * Ns[ j ]
//                )
//            );
//
//            std::transform(
//                std::begin( dCdGradV_p ), std::end( dCdGradV_p ), dCdGradV_begin + node_count * dim * dim * i + dim * dim * j,
//                std::bind(
//                    std::multiplies<
//                        typename std::iterator_traits<dCdGradV_iter_out>::value_type
//                    >( ),
//                    std::placeholders::_1,
//                    Ns[ i ] * Ns[ j ]
//                )
//            );

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
        0.53182759, 0.63440096, 0.84943179, 0.51044152, 0.65634786,
        0.86791477, 0.48312667, 0.6486585 , 0.86600796, 0.50451273,
        0.6267116 , 0.84752498, 0.53695906, 0.68247738, 0.83864355,
        0.515573  , 0.70442428, 0.85712652, 0.48825814, 0.69673492,
        0.85521971, 0.50964421, 0.67478802, 0.83673674
    };

    std::array< floatType, 24 > u_tp1 = {
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
        -0.0830957 , -0.21734036, -0.003563  , -0.00136224, -0.12179146,
        -0.31855078, -0.0052222 , -0.0019966
    };

    std::array< floatType, 8 > result;

    std::array< floatType, 3 > local_point = {
        0.44683272, -0.96774159,  0.18886376
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    evaluate_at_nodes<3, 8>(
        std::cbegin( local_point ), std::cend( local_point ), dt,
        std::cbegin( density_t ), std::cend( density_t ),
        std::cbegin( density_tp1 ), std::cend( density_tp1 ),
        std::cbegin( u_t ), std::cend( u_t ),
        std::cbegin( u_tp1 ), std::cend( u_tp1 ),
        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
        std::cbegin( v_t ), std::cend( v_t ),
        std::cbegin( X ), std::cend( X ),
        alpha,
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0 );

    std::array< floatType, 1 * 1 * 8 * 8 > dCdRho, dCdRhoDot;

    std::array< floatType, 3 * 3 * 8 * 8 > dCdGradRho, dCdV, dCdU;

    std::array< floatType, 9 * 9 * 8 * 8 > dCdGradV;

    evaluate_at_nodes<3, 8>(
        std::cbegin( local_point ), std::cend( local_point ), dt,
        std::cbegin( density_t ), std::cend( density_t ),
        std::cbegin( density_tp1 ), std::cend( density_tp1 ),
        std::cbegin( u_t ), std::cend( u_t ),
        std::cbegin( u_tp1 ), std::cend( u_tp1 ),
        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
        std::cbegin( v_t ), std::cend( v_t ),
        std::cbegin( X ), std::cend( X ),
        alpha,
        std::begin( result ), std::end( result ),
        std::begin( dCdRho ), std::end( dCdRho ),
        std::begin( dCdRhoDot ), std::end( dCdRhoDot ),
        std::begin( dCdGradRho ), std::end( dCdGradRho ),
        std::begin( dCdV ), std::end( dCdV ),
        std::begin( dCdGradV ), std::end( dCdGradV ),
        std::begin( dCdU ), std::end( dCdU )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    // Check the derivatives w.r.t. the density
    {

        constexpr vardim = 8;
        constexpr outdim = 8;

        for ( unsigned int i = 0; i < vardim; i++ ){
    
            floatType delta = eps * std::fabs( rho_tp1[ i ] ) + eps;
    
            std::array< floatType, vardim > xp = density_tp1;
            std::array< floatType, vardim > xm = density_tp1;
    
            xp[ i ] += delta;
            xm[ i ] += delta;
    
            std::array< floatType, vardim > vp, vm;
    
            evaluate_at_nodes<3, 8>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( xp ), std::cend( xp ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vp ), std::end( vp )
            );

            evaluate_at_nodes<3, 8>(
                std::cbegin( local_point ), std::cend( local_point ), dt,
                std::cbegin( density_t ), std::cend( density_t ),
                std::cbegin( xm ), std::cend( xm ),
                std::cbegin( u_t ), std::cend( u_t ),
                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
                std::cbegin( v_t ), std::cend( v_t ),
                std::cbegin( X ), std::cend( X ),
                alpha,
                std::begin( vm ), std::end( vm )
            );

            for ( unsigned int j = 0; j < outdim; j++ ){

                BOOST_TEST( *( dCdRho_begin + vardim * j + i ) == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

//    // Check the derivatives w.r.t. the deformation
//    {
//
//        constexpr vardim = 3 * 8;
//        constexpr outdim = 3 * 8;
//
//        for ( unsigned int i = 0; i < vardim; i++ ){
//    
//            floatType delta = eps * std::fabs( u_tp1[ i ] ) + eps;
//    
//            std::array< floatType, vardim > xp = u_tp1;
//            std::array< floatType, vardim > xm = u_tp1;
//    
//            xp[ i ] += delta;
//            xm[ i ] += delta;
//    
//            std::array< floatType, vardim > vp, vm;
//    
//            evaluate_at_nodes<3, 8>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( xp ), std::cend( xp ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vp ), std::end( vp )
//            );
//
//            evaluate_at_nodes<3, 8>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( xm ), std::cend( xm ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vp ), std::end( vp )
//            );
//
//            for ( unsigned int j = 0; j < outdim; j++ ){
//
//                unsigned int node = outdim / 3;
//
//                unsigned int direction = j - 3 * node;
//
//                BOOST_TEST( *( dCdU_begin + vardim * j + i ) == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );
//
//            }
//
//        }
//
//    }

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
    for ( unsigned int i = 0; i < nphases; i++ ){

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

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( xp ), std::end( xp ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( xm ), std::end( xm ),
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

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( xp ), std::end( xp ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
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

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( xp ), std::end( xp ),
                                                                              std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
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

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
                                                                              std::begin( density_gradient ), std::end( density_gradient ), std::begin( velocity ), std::end( velocity ),
                                                                              std::begin( xp ), std::end( xp ), std::begin( vp ), std::end( vp ) );

        tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( std::begin( density ), std::end( density ), std::begin( density_dot ), std::end( density_dot ),
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
