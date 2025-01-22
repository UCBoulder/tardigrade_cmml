/**
  * \file test_tardigrade_balance_equations_balance_of_volume_fraction.cpp
  *
  * Tests for tardigrade_balance_equations_balance_of_volume_fraction
  */

#include<tardigrade_balance_of_volume_fraction.h>
#include<tardigrade_finite_element_utilities.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define USE_EIGEN
#include<tardigrade_vector_tools.h>

#define BOOST_TEST_MODULE test_tardigrade_balance_equations_balance_of_volume_fraction
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
    class xi_in, typename dt_type,
    class density_t_in, class density_tp1_in,
    class phi_t_in, class phi_tp1_in,
    class u_t_in, class u_tp1_in, class umesh_t_in, class umesh_tp1_in,
    class phi_dot_t_in, class u_dot_t_in,
    class X_in,
    class mass_change_rate_iter, class rest_density_iter, class trace_mass_change_velocity_gradient_iter,
    typename alpha_type, typename beta_type, class value_out
>
void evaluate_at_nodes(
    const xi_in &xi_begin, const xi_in &xi_end, dt_type dt,
    const density_t_in &density_t_begin, const density_t_in &density_t_end,
    const density_tp1_in &density_tp1_begin, const density_tp1_in &density_tp1_end,
    const phi_t_in &phi_t_begin, const phi_t_in &phi_t_end,
    const phi_tp1_in &phi_tp1_begin, const phi_tp1_in &phi_tp1_end,
    const u_t_in &u_t_begin, const u_t_in &u_t_end,
    const u_tp1_in &u_tp1_begin, const u_tp1_in &u_tp1_end,
    const umesh_t_in &umesh_t_begin, const umesh_t_in &umesh_t_end,
    const umesh_tp1_in &umesh_tp1_begin, const umesh_tp1_in &umesh_tp1_end,
    const phi_dot_t_in &phi_dot_t_begin, const phi_dot_t_in &phi_dot_t_end,
    const u_dot_t_in &u_dot_t_begin, const u_dot_t_in &u_dot_t_end,
    const X_in &X_begin, const X_in &X_end,
    const mass_change_rate_iter &mass_change_rate_begin, const mass_change_rate_iter &mass_change_rate_end,
    const rest_density_iter &rest_density_begin, const rest_density_iter &rest_density_end,
    const trace_mass_change_velocity_gradient_iter &trace_mass_change_velocity_gradient_begin,
    const trace_mass_change_velocity_gradient_iter &trace_mass_change_velocity_gradient_end,
    const alpha_type &alpha, const beta_type &beta, value_out value_begin, value_out value_end
){

    // Update the mesh nodes
    std::array< typename std::iterator_traits<umesh_t_in  >::value_type, dim * node_count > x_t;
    std::array< typename std::iterator_traits<umesh_tp1_in>::value_type, dim * node_count > x_tp1;

    std::transform( X_begin, X_end,   umesh_t_begin,   std::begin( x_t ), std::plus<typename std::iterator_traits<umesh_t_in  >::value_type>( ) );
    std::transform( X_begin, X_end, umesh_tp1_begin, std::begin( x_tp1 ), std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>( ) );

    // Calculate the current rates of change
    std::array< typename std::iterator_traits<phi_tp1_in>::value_type, node_count * nphases > phi_dot_tp1;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count * nphases > u_dot_tp1;

    floatType dPhiDotdPhi, dUDotdU;

    compute_current_rate_of_change(
        dt, phi_t_begin, phi_t_end, phi_tp1_begin, phi_tp1_end,
        phi_dot_t_begin, phi_dot_t_end, alpha,
        std::begin( phi_dot_tp1 ), std::end( phi_dot_tp1 ),
        dPhiDotdPhi
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
    > density_tp1_p;

    std::array<
        typename std::iterator_traits<phi_tp1_in>::value_type, nphases
    > phi_tp1_p, phi_dot_tp1_p;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * nphases
    > u_dot_tp1_p;

    // Interpolate quantities to the local point

    e.InterpolateQuantity(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( density_tp1_p ), std::end( density_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, phi_tp1_begin, phi_tp1_end,
        std::begin( phi_tp1_p ), std::end( phi_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::cbegin( phi_dot_tp1 ), std::cend( phi_dot_tp1 ),
        std::begin( phi_dot_tp1_p ), std::end( phi_dot_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::cbegin( u_dot_tp1 ), std::cend( u_dot_tp1 ),
        std::begin( u_dot_tp1_p ), std::end( u_dot_tp1_p )
    );

    // Compute the gradients at the local point

    std::array<
        typename std::iterator_traits<phi_tp1_in>::value_type, dim * nphases
    > grad_phi_tp1;

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, phi_tp1_begin, phi_tp1_end,
        std::begin( grad_phi_tp1 ), std::end( grad_phi_tp1 )
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

            tardigradeBalanceEquations::balanceOfVolumeFraction::computeBalanceOfVolumeFraction<dim>(
                density_tp1_p[ 0 ],
                std::cbegin( u_dot_tp1_p ), std::cend( u_dot_tp1_p ),
                phi_tp1_p[ 0 ], phi_dot_tp1_p[ 0 ],
                std::cbegin( grad_phi_tp1 ), std::cend( grad_phi_tp1 ),
                *( mass_change_rate_begin + 0 ),
                *( rest_density_begin + 0 ),
                *( trace_mass_change_velocity_gradient_begin + 0 ),
                Ns[ i ],
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

            tardigradeBalanceEquations::balanceOfVolumeFraction::computeBalanceOfVolumeFraction<dim>(
                std::cbegin( density_tp1_p ), std::cend( density_tp1_p ),
                std::cbegin( u_dot_tp1_p ),   std::cend( u_dot_tp1_p ),
                std::cbegin( phi_tp1_p ),     std::cend( phi_tp1_p ),
                std::cbegin( phi_dot_tp1_p ), std::cend( phi_dot_tp1_p ),
                std::cbegin( grad_phi_tp1 ),  std::cend( grad_phi_tp1 ),
                mass_change_rate_begin,       mass_change_rate_end,
                rest_density_begin,           rest_density_end,
                trace_mass_change_velocity_gradient_begin,
                trace_mass_change_velocity_gradient_end,
                Ns[ i ],
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

BOOST_AUTO_TEST_CASE( test_computeBalanceOfVolumeFraction_fea, * boost::unit_test::tolerance( 1e-5 ) ){

    constexpr unsigned int nphases = 1;

    std::array< floatType, 8 * nphases > density_t = {
        0.61289453, 0.12062867, 0.8263408 , 0.60306013, 0.54506801,
        0.34276383, 0.30412079, 0.41702221
    };

    std::array< floatType, 8 * nphases > density_tp1 = {
        0.68130077, 0.87545684, 0.51042234, 0.66931378, 0.58593655,
        0.6249035 , 0.67468905, 0.84234244
    };

    std::array< floatType, 8 * nphases > phi_t = {
        0.20454286, 0.45063649, 0.54776357, 0.09332671, 0.29686078,
        0.92758424, 0.56900373, 0.457412  
    };

    std::array< floatType, 8 * nphases > phi_tp1 = {
        0.30646975, 0.66526147, 0.11139217, 0.66487245, 0.88785679,
        0.69631127, 0.44032788, 0.43821438
    };

    std::array< floatType, 8 * 3 * nphases > u_t = {
        -0.83361002,  0.52736568, -0.51266725, -0.61155408,  0.14491391,
        -0.80857497,  0.77065365,  0.25449794,  0.44683272, -0.96774159,
         0.18886376,  0.11357038, -0.68208071, -0.69385897,  0.39105906,
        -0.36246715,  0.38394059,  0.1087665 , -0.22209885,  0.85026498,
         0.68333999, -0.28520487, -0.91281707, -0.39046385
    };

    std::array< floatType, 8 * 3 * nphases > u_tp1 = {
         0.65822527, -0.32265831,  0.10474015,  0.15710294,  0.04306612,
        -0.99462387,  0.97669084,  0.81068315, -0.58472828, -0.41502117,
         0.04002031,  0.80382275,  0.96726177, -0.48491587,  0.12871809,
         0.61393737, -0.21125989,  0.46214607, -0.67786197,  0.20139714,
         0.73172892,  0.96704322, -0.84126842, -0.14330545
    };

    std::array< floatType, 8 * 3 > umesh_t = {
        0.53182759, 0.63440096, 0.84943179, 0.51044152, 0.65634786,
        0.86791477, 0.48312667, 0.6486585 , 0.86600796, 0.50451273,
        0.6267116 , 0.84752498, 0.53695906, 0.68247738, 0.83864355,
        0.515573  , 0.70442428, 0.85712652, 0.48825814, 0.69673492,
        0.85521971, 0.50964421, 0.67478802, 0.83673674
    };

    std::array< floatType, 8 * 3 > umesh_tp1 = {
        0.72445532, 0.61102351, 0.72244338, 0.70877313, 0.5669913 ,
        0.69069256, 0.7316781 , 0.55679573, 0.65823773, 0.7473603 ,
        0.60082794, 0.68998856, 0.71831255, 0.63482305, 0.72559852,
        0.70263035, 0.59079084, 0.69384769, 0.72553532, 0.58059527,
        0.66139287, 0.74121752, 0.62462748, 0.6931437 
    };

    std::array< floatType, 8 * nphases > phi_dot_t = {
         0.50705198,  0.4837243 , -0.90284193,  0.41739479,  0.6784867 ,
        -0.66812423,  0.56199588, -0.42692677
    };

    std::array< floatType, 8 * 3 * nphases > u_dot_t = {
        -0.20362864,  0.40991766,  0.99071696, -0.28817027,  0.52509563,
         0.18635383,  0.3834036 , -0.6977451 , -0.20224741, -0.5182882 ,
        -0.31308797,  0.02625631,  0.3332491 , -0.78818303, -0.7382101 ,
        -0.35603879,  0.32312867,  0.69301245,  0.10651469,  0.70890498,
        -0.23032438, -0.36642421, -0.29147065, -0.65783634
    };

    std::array< floatType, 24 > X = {
        0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0., 0., 1., 1., 0.,
        1., 1., 1., 1., 0., 1., 1.
    };

    std::array< floatType, nphases > mass_change_rate = {
        0.53019219
    };

    std::array< floatType, nphases > rest_density = {
        0.565642
    };

    std::array< floatType, nphases > trace_mass_change_velocity_gradient = {
        0.08490416
    };

    std::array< floatType, 8 * nphases > answer = {
        -0.00277397, -0.00387299, -0.0170444 , -0.01220778, -0.00141041,
        -0.00196921, -0.00866617, -0.00620701
    };

    std::array< floatType, 8 > result;

    std::array< floatType, 3 > local_point = {
         0.16534218,  0.62968741, -0.32586723
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    floatType beta = 0.67;

    evaluate_at_nodes< 3, 8, nphases >(
        std::cbegin( local_point ), std::cend( local_point ),
        dt,
        std::cbegin( density_t ),        std::cend( density_t ),
        std::cbegin( density_tp1 ),      std::cend( density_tp1 ),
        std::cbegin( phi_t ),            std::cend( phi_t ),
        std::cbegin( phi_tp1 ),          std::cend( phi_tp1 ),
        std::cbegin( u_t ),              std::cend( u_t ),
        std::cbegin( u_tp1 ),            std::cend( u_tp1 ),
        std::cbegin( umesh_t ),          std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),        std::cend( umesh_tp1 ),
        std::cbegin( phi_dot_t ),        std::cend( phi_dot_t ),
        std::cbegin( u_dot_t ),          std::cend( u_dot_t ),
        std::cbegin( X ),                std::cend( X ),
        std::cbegin( mass_change_rate ), std::cend( mass_change_rate ),
        std::cbegin( rest_density ),     std::cend( rest_density ),
        std::cbegin( trace_mass_change_velocity_gradient ),
        std::cend(   trace_mass_change_velocity_gradient ),
        alpha, beta,
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0 );

}

BOOST_AUTO_TEST_CASE( test_computeBalanceOfVolumeFraction_multiphase_fea, * boost::unit_test::tolerance( 1e-5 ) ){

    constexpr unsigned int nphases = 4;

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
        0.30646975, 0.66526147, 0.11139217, 0.66487245, 0.88785679,
        0.69631127, 0.44032788, 0.43821438, 0.7650961 , 0.565642  ,
        0.08490416, 0.58267109, 0.8148437 , 0.33706638, 0.92757658,
        0.750717  , 0.57406383, 0.75164399, 0.07914896, 0.85938908,
        0.82150411, 0.90987166, 0.1286312 , 0.08178009, 0.13841557,
        0.39937871, 0.42430686, 0.56221838, 0.12224355, 0.2013995 ,
        0.81164435, 0.46798757
    };

    std::array< floatType, 8 * nphases > phi_t = {
        2.77423960e-01, 9.98918406e-01, 4.06161246e-02, 6.45822522e-01,
        3.86995850e-02, 7.60210258e-01, 2.30089957e-01, 8.98318671e-02,
        6.48449712e-01, 7.32601217e-01, 6.78095315e-01, 5.19009471e-02,
        2.94306946e-01, 4.51088346e-01, 2.87103290e-01, 8.10513456e-01,
        1.31115105e-01, 6.12179362e-01, 9.88214944e-01, 9.02556539e-01,
        2.22157062e-01, 8.18876137e-05, 9.80597342e-01, 8.82712985e-01,
        9.19472466e-01, 4.15503551e-01, 7.44615462e-01, 2.12831499e-01,
        3.92304071e-01, 8.51548051e-01, 1.27612224e-01, 8.93865368e-01
    };

    std::array< floatType, 8 * nphases > phi_tp1 = {
        0.56366459, 0.72707995, 0.6711266 , 0.24751315, 0.52486622,
        0.53766344, 0.71680336, 0.35986735, 0.7977326 , 0.62792185,
        0.03833161, 0.54647902, 0.8619121 , 0.56757416, 0.17582827,
        0.51037637, 0.75694584, 0.1101052 , 0.81709908, 0.16748164,
        0.53407649, 0.38574348, 0.24862377, 0.64743251, 0.03739211,
        0.76004581, 0.52694064, 0.87577121, 0.52071832, 0.03503317,
        0.14360097, 0.79560459
    };

    std::array< floatType, 8 * 3 * nphases > u_t = {
         0.61587642, -0.98514724,  0.10318545,  0.8638643 ,  0.16435092,
        -0.58780855,  0.43551512, -0.2420283 ,  0.33676789, -0.94136055,
         0.27180072, -0.93560413,  0.48956131, -0.054174  , -0.75649129,
         0.08527185, -0.86645111,  0.30672974,  0.99217265,  0.53879467,
         0.14754823, -0.79472948,  0.39966815,  0.32233573, -0.90180574,
         0.5845986 ,  0.03743318, -0.14826461,  0.57637435, -0.17686155,
        -0.03794745, -0.63674231, -0.3573622 ,  0.69106599, -0.6261925 ,
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
         0.77518692
    };

    std::array< floatType, 8 * 3 * nphases > u_tp1 = {
        -0.07013029,  0.55953332, -0.52504356, -0.33483946,  0.90739424,
         0.31563015,  0.54575566,  0.37674869, -0.59139176, -0.0586225 ,
         0.61792775,  0.35007025, -0.98794423, -0.82518451, -0.30641056,
         0.88873108, -0.01761904, -0.45964747, -0.27915256, -0.57869474,
        -0.15759989, -0.56392912,  0.69150501, -0.0874588 , -0.44039596,
         0.8657833 , -0.37129729,  0.81942932, -0.91316382,  0.41423012,
        -0.03222192, -0.11155788, -0.92735331, -0.91863362, -0.33449277,
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
        -0.30096295
    };

    std::array< floatType, 8 * 3 > umesh_t = {
        0.53182759, 0.63440096, 0.84943179, 0.51044152, 0.65634786,
        0.86791477, 0.48312667, 0.6486585 , 0.86600796, 0.50451273,
        0.6267116 , 0.84752498, 0.53695906, 0.68247738, 0.83864355,
        0.515573  , 0.70442428, 0.85712652, 0.48825814, 0.69673492,
        0.85521971, 0.50964421, 0.67478802, 0.83673674
    };

    std::array< floatType, 8 * 3 > umesh_tp1 = {
        0.72445532, 0.61102351, 0.72244338, 0.70877313, 0.5669913 ,
        0.69069256, 0.7316781 , 0.55679573, 0.65823773, 0.7473603 ,
        0.60082794, 0.68998856, 0.71831255, 0.63482305, 0.72559852,
        0.70263035, 0.59079084, 0.69384769, 0.72553532, 0.58059527,
        0.66139287, 0.74121752, 0.62462748, 0.6931437 
    };

    std::array< floatType, 8 * nphases > phi_dot_t = {
        -0.00698406, -0.14780869, -0.38870722,  0.83369757,  0.03524692,
         0.60805274,  0.71530357,  0.84476471, -0.39323853, -0.32037829,
         0.19014775, -0.11735173,  0.86568507, -0.2048719 , -0.0444439 ,
         0.23437218, -0.19052103,  0.98495687, -0.80229743, -0.55879336,
        -0.35468974, -0.70455431, -0.43156153,  0.55849059,  0.045784  ,
        -0.93209273,  0.96524517,  0.23201296, -0.88212104,  0.32233754,
        -0.24326126, -0.72865341
    };

    std::array< floatType, 8 * 3 * nphases > u_dot_t = {
        -0.96776274, -0.74608394,  0.55432492, -0.90820954,  0.42199739,
         0.94209228,  0.74336587,  0.4203233 ,  0.91701949, -0.14037332,
         0.74575783, -0.28808466,  0.85952731, -0.70244469,  0.88005803,
         0.66543239,  0.69210968, -0.75215398,  0.1929738 , -0.96721504,
         0.44236873, -0.98452497, -0.83035545, -0.54900318,  0.75024907,
        -0.27284736,  0.07991987,  0.13620643, -0.54907328,  0.14429354,
         0.32190359, -0.40350921, -0.16274628, -0.09382215,  0.86470132,
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
         0.64114844
    };

    std::array< floatType, 24 > X = {
        0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0., 0., 0., 1., 1., 0.,
        1., 1., 1., 1., 0., 1., 1.
    };

    std::array< floatType, nphases > mass_change_rate = {
        -0.0160479 , -0.11624146, -0.36313044, -0.4309016 
    };

    std::array< floatType, nphases > rest_density = {
        0.96588631, 0.43296933, 0.88400303, 0.64816312
    };

    std::array< floatType, nphases > trace_mass_change_velocity_gradient = {
        0.85842765, 0.85244954, 0.95631203, 0.69794224
    };

    std::array< floatType, 8 * nphases > answer = {
        0.000889  , 0.00297344, 0.00036272, 0.00416856, 0.00367928,
        0.01230606, 0.00150118, 0.01725227, 0.01010739, 0.03380613,
        0.00412391, 0.04739395, 0.00244219, 0.00816836, 0.00099643,
        0.01145151, 0.00136293, 0.00455858, 0.00055609, 0.00639082,
        0.0056407 , 0.01886642, 0.00230146, 0.02644946, 0.01549564,
        0.05182818, 0.00632236, 0.07265968, 0.00374412, 0.01252292,
        0.00152763, 0.01755631
    };

    std::array< floatType, 8 * nphases > result;

    std::array< floatType, 3 > local_point = {
        0.61079387, 0.46625579, 0.21045367
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    floatType beta = 0.67;

    evaluate_at_nodes< 3, 8, nphases >(
        std::cbegin( local_point ), std::cend( local_point ),
        dt,
        std::cbegin( density_t ),        std::cend( density_t ),
        std::cbegin( density_tp1 ),      std::cend( density_tp1 ),
        std::cbegin( phi_t ),            std::cend( phi_t ),
        std::cbegin( phi_tp1 ),          std::cend( phi_tp1 ),
        std::cbegin( u_t ),              std::cend( u_t ),
        std::cbegin( u_tp1 ),            std::cend( u_tp1 ),
        std::cbegin( umesh_t ),          std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),        std::cend( umesh_tp1 ),
        std::cbegin( phi_dot_t ),        std::cend( phi_dot_t ),
        std::cbegin( u_dot_t ),          std::cend( u_dot_t ),
        std::cbegin( X ),                std::cend( X ),
        std::cbegin( mass_change_rate ), std::cend( mass_change_rate ),
        std::cbegin( rest_density ),     std::cend( rest_density ),
        std::cbegin( trace_mass_change_velocity_gradient ),
        std::cend(   trace_mass_change_velocity_gradient ),
        alpha, beta,
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0 );

}
