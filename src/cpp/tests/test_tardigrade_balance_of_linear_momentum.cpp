/**
  * \file test_tardigrade_balance_equations_balance_of_linear_momentum.cpp
  *
  * Tests for tardigrade_balance_equations_balance_of_linear_momentum
  */

#include<tardigrade_balance_of_linear_momentum.h>
#include<tardigrade_finite_element_utilities.h>
#define USE_EIGEN
#include<tardigrade_vector_tools.h>
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

BOOST_AUTO_TEST_CASE( test_computeBalanceOfLinearMomentumNonDivergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    constexpr unsigned int dim = 3;

    constexpr unsigned int sot_dim = 9;

    floatType density = 0.69646919;

    floatType density_dot = 0.28613933;

    floatVector density_gradient = { 0.22685145, 0.55131477, 0.71946897 };

    floatVector velocity = { 0.42310646, 0.9807642,  0.68482974 };

    floatVector velocity_dot = { 0.4809319,  0.39211752, 0.34317802 };

    secondOrderTensor velocity_gradient = { 0.72904971, 0.43857224, 0.0596779,  0.39804426, 0.73799541, 0.18249173, 0.17545176, 0.53155137, 0.53182759 };

    floatVector body_force = { 0.63440096, 0.84943179, 0.72445532 };

    floatVector answer = { -1.62394626, -3.14362652, -2.32569958 };

    floatVector result;

    tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
        density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
        std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
        std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ), std::end( body_force ),
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0. );

    floatVector dRdRho, dRdRhoDot;

    secondOrderTensor dRdGradRho, dRdV, dRdVDot, dRdB;

    std::array< floatType, dim * dim * dim > dRdGradV;

    tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
        density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
        std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
        std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ), std::end( body_force ),
        std::begin( result ),     std::end( result ),
        std::begin( dRdRho ),     std::end( dRdRho ),
        std::begin( dRdRhoDot ),  std::end( dRdRhoDot ),
        std::begin( dRdGradRho ), std::end( dRdGradRho ),
        std::begin( dRdV ),       std::end( dRdV ),
        std::begin( dRdVDot ),    std::end( dRdVDot ),
        std::begin( dRdGradV ),   std::end( dRdGradV ),
        std::begin( dRdB ),       std::end( dRdB )
    );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( density ) + eps;

        floatType xp = density;
        floatType xm = density;

        xp += delta;
        xm -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            xp, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vp ), std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            xm, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vm ), std::end( vm )
        );

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

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, xp, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vp ), std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, xm, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vm ), std::end( vm )
        );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdRhoDot[ j ] == grad );

        }

    }

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::fabs( density_gradient[ i ] ) + eps;

        floatVector xp = density_gradient;
        floatVector xm = density_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, density_dot, std::begin( xp ), std::end( xp ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vp ), std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, density_dot, std::begin( xm ), std::end( xm ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vm ), std::end( vm )
        );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdGradRho[ dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::fabs( velocity[ i ] ) + eps;

        floatVector xp = velocity;
        floatVector xm = velocity;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( xp ), std::end( xp ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vp ), std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( xm ), std::end( xm ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vm ), std::end( vm )
        );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdV[ dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::fabs( velocity_dot[ i ] ) + eps;

        floatVector xp = velocity_dot;
        floatVector xm = velocity_dot;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( xp ), std::end( xp ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vp ), std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( xm ), std::end( xm ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vm ), std::end( vm )
        );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdVDot[ dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < sot_dim; i++ ){

        floatType delta = eps * std::fabs( velocity_gradient[ i ] ) + eps;

        secondOrderTensor xp = velocity_gradient;
        secondOrderTensor xm = velocity_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( xp ), std::end( xp ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vp ), std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( xm ), std::end( xm ), std::begin( body_force ),
            std::end( body_force ),
            std::begin( vm ), std::end( vm )
        );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdGradV[ dim * dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::fabs( body_force[ i ] ) + eps;

        floatVector xp = body_force;
        floatVector xm = body_force;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( xp ),
            std::end( xp ),
            std::begin( vp ), std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            density, density_dot, std::begin( density_gradient ), std::end( density_gradient ),
            std::begin( velocity ), std::end( velocity ), std::begin( velocity_dot ), std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ), std::begin( xm ),
            std::end( xm ),
            std::begin( vm ), std::end( vm )
        );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdB[ dim * j + i ] == grad );

        }

    }

}

BOOST_AUTO_TEST_CASE( test_computeBalanceOfLinearMomentumDivergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    constexpr unsigned int dim = 3;

    constexpr unsigned int sot_dim = 9;

    floatVector test_function_gradient = { 0.61102351, 0.72244338, 0.32295891 };

    secondOrderTensor cauchy_stress = { 0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,  0.42583029 };

    floatType volume_fraction = 0.31226122;

    floatVector answer = { -0.25482292, -0.11411737, -0.19682338 };

    floatVector result;

    tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
        std::begin( test_function_gradient ), std::end( test_function_gradient ),
        std::begin( cauchy_stress ),          std::end( cauchy_stress ),
        volume_fraction,
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0. );

    floatVector dRdPhi;

    secondOrderTensor dRdGradPsi;

    std::array< floatType, dim * dim * dim > dRdCauchy;

    tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
        std::begin( test_function_gradient ), std::end( test_function_gradient ),
        std::begin( cauchy_stress ),          std::end( cauchy_stress ),
        volume_fraction,
        std::begin( result ),     std::end( result ),
        std::begin( dRdGradPsi ), std::end( dRdGradPsi ),
        std::begin( dRdCauchy ),  std::end( dRdCauchy ),
        std::begin( dRdPhi ),     std::end( dRdPhi )
    );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::fabs( test_function_gradient[ i ] ) + eps;

        floatVector xp = test_function_gradient;
        floatVector xm = test_function_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( xp ), std::end( xp ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ),
            volume_fraction,
            std::begin( vp ), std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( xm ), std::end( xm ),
            std::begin( cauchy_stress ), std::end( cauchy_stress ),
            volume_fraction,
            std::begin( vm ), std::end( vm )
        );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdGradPsi[ dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < sot_dim; i++ ){

        floatType delta = eps * std::fabs( cauchy_stress[ i ] ) + eps;

        secondOrderTensor xp = cauchy_stress;
        secondOrderTensor xm = cauchy_stress;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( test_function_gradient ), std::end( test_function_gradient ),
            std::begin( xp ),                     std::end( xp ),
            volume_fraction,
            std::begin( vp ), std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( test_function_gradient ), std::end( test_function_gradient ),
            std::begin( xm ),          std::end( xm ),
            volume_fraction,
            std::begin( vm ), std::end( vm )
        );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdCauchy[ dim * dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( volume_fraction ) + eps;

        floatType xp = volume_fraction;
        floatType xm = volume_fraction;

        xp += delta;
        xm -= delta;

        floatVector vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( test_function_gradient ), std::end( test_function_gradient ),
            std::begin( cauchy_stress ),          std::end( cauchy_stress ),
            xp,
            std::begin( vp ), std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( test_function_gradient ), std::end( test_function_gradient ),
            std::begin( cauchy_stress ),          std::end( cauchy_stress ),
            xm,
            std::begin( vm ), std::end( vm )
        );

        for ( unsigned int j = 0; j < dim; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdPhi[ j ] == grad );

        }

    }

}

BOOST_AUTO_TEST_CASE( test_multiphase_computeBalanceOfLinearMomentumNonDivergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    constexpr unsigned int nphases = 5;

    constexpr unsigned int dim = 3;

    constexpr unsigned int sot_dim = dim * dim;

    std::array<floatType,nphases> density = { 0.69646919, 0.28613933, 0.22685145, 0.55131477, 0.71946897 };

    std::array<floatType,nphases> density_dot = { 0.42310646, 0.9807642,  0.68482974, 0.4809319 , 0.39211752 };

    std::array<floatType,dim*nphases> density_gradient = { 0.34317802, 0.72904971, 0.43857224, 0.0596779 , 0.39804426, 0.73799541,
                                                           0.18249173, 0.17545176, 0.53155137, 0.53182759, 0.63440096, 0.84943179,
                                                           0.72445532, 0.61102351, 0.72244338 };

    std::array<floatType,dim*nphases> velocity = { 0.32295891, 0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494,
                                                   0.43370117, 0.43086276, 0.4936851 , 0.42583029, 0.31226122, 0.42635131,
                                                   0.89338916, 0.94416002, 0.50183668 };

    std::array<floatType,dim*nphases> velocity_dot = { 0.62395295, 0.1156184 , 0.31728548, 0.41482621, 0.86630916, 0.25045537,
                                                       0.48303426, 0.98555979, 0.51948512, 0.61289453, 0.12062867, 0.8263408,
                                                       0.60306013, 0.54506801, 0.34276383 };

    std::array<floatType,dim*dim*nphases> velocity_gradient = { 0.30412079, 0.41702221, 0.68130077, 0.87545684, 0.51042234, 0.66931378,
                                                                0.58593655, 0.6249035 , 0.67468905, 0.84234244, 0.08319499, 0.76368284,
                                                                0.24366637, 0.19422296, 0.57245696, 0.09571252, 0.88532683, 0.62724897,
                                                                0.72341636, 0.01612921, 0.59443188, 0.55678519, 0.15895964, 0.15307052,
                                                                0.69552953, 0.31876643, 0.6919703 , 0.55438325, 0.38895057, 0.92513249,
                                                                0.84167   , 0.35739757, 0.04359146, 0.30476807, 0.39818568, 0.70495883,
                                                                0.99535848, 0.35591487, 0.76254781, 0.59317692, 0.6917018 , 0.15112745,
                                                                0.39887629, 0.2408559 , 0.34345601 };

    std::array<floatType,dim*nphases> body_force = { 0.51312815, 0.66662455, 0.10590849, 0.13089495, 0.32198061, 0.66156434,
                                                     0.84650623, 0.55325734, 0.85445249, 0.38483781, 0.3167879 , 0.35426468,
                                                     0.17108183, 0.82911263, 0.33867085 };

    std::array<floatType,dim*nphases> answer = { -0.98391364, -0.74824486, -0.98542636, -0.71396186, -1.35804429,
                                                 -0.2319744 , -0.68969665, -0.81421452, -0.82144783, -1.45965944,
                                                 -0.83539248, -1.48655174, -4.54064701, -3.94895894, -2.27310718 };

    std::array<floatType,dim*nphases> result;

    tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
        std::begin( density ),           std::end( density ),
        std::begin( density_dot ),       std::end( density_dot ),
        std::begin( density_gradient ),  std::end( density_gradient ),
        std::begin( velocity ),          std::end( velocity ),
        std::begin( velocity_dot ),      std::end( velocity_dot ),
        std::begin( velocity_gradient ), std::end( velocity_gradient ),
        std::begin( body_force ),        std::end( body_force ),
        std::begin( result ),            std::end( result )
    );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0. );

    std::array<floatType,dim*nphases> dRdRho, dRdRhoDot;

    std::array<floatType,dim*dim*nphases> dRdGradRho, dRdV, dRdVDot, dRdB;

    std::array< floatType, dim * dim * dim * nphases > dRdGradV;

    tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
        std::begin( density ),           std::end( density ),
        std::begin( density_dot ),       std::end( density_dot ),
        std::begin( density_gradient ),  std::end( density_gradient ),
        std::begin( velocity ),          std::end( velocity ),
        std::begin( velocity_dot ),      std::end( velocity_dot ),
        std::begin( velocity_gradient ), std::end( velocity_gradient ),
        std::begin( body_force ),        std::end( body_force ),
        std::begin( result ),            std::end( result ),
        std::begin( dRdRho ),            std::end( dRdRho ),
        std::begin( dRdRhoDot ),         std::end( dRdRhoDot ),
        std::begin( dRdGradRho ),        std::end( dRdGradRho ),
        std::begin( dRdV ),              std::end( dRdV ),
        std::begin( dRdVDot ),           std::end( dRdVDot ),
        std::begin( dRdGradV ),          std::end( dRdGradV ),
        std::begin( dRdB ),              std::end( dRdB )
    );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 1 * nphases; i++ ){

        floatType delta = eps * std::fabs( density[ i ] ) + eps;

        std::array<floatType,nphases> xp = density;
        std::array<floatType,nphases> xm = density;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,dim*nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( xp ),                std::end( xp ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vp ),                std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( xm ),                std::end( xm ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vm ),                std::end( vm )
        );

        for ( unsigned int j = 0; j < dim * nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( j / dim ) == i ){

                BOOST_TEST( dRdRho[ j ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < 1 * nphases; i++ ){

        floatType delta = eps * std::fabs( density_dot[ i ] ) + eps;

        std::array<floatType,nphases> xp = density_dot;
        std::array<floatType,nphases> xm = density_dot;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,dim*nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( xp ),                std::end( xp ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vp ),                std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( xm ),                std::end( xm ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vm ),                std::end( vm )
        );

        for ( unsigned int j = 0; j < dim * nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( j / dim ) == i ){

                BOOST_TEST( dRdRhoDot[ j ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < dim * nphases; i++ ){

        floatType delta = eps * std::fabs( density_gradient[ i ] ) + eps;

        std::array<floatType,dim*nphases> xp = density_gradient;
        std::array<floatType,dim*nphases> xm = density_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,dim*nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( xp ),                std::end( xp ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vp ),                std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( xm ),                std::end( xm ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vm ),                std::end( vm )
        );

        for ( unsigned int j = 0; j < dim * nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( j / dim ) == ( unsigned int )( i / dim ) ){

                unsigned int phase = ( unsigned int )( j / dim );
                unsigned int row   = ( j - phase * dim ) % nphases;
                unsigned int col   = ( i - phase * dim ) % dim;

                BOOST_TEST( dRdGradRho[ sot_dim * phase + dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < dim * nphases; i++ ){

        floatType delta = eps * std::fabs( velocity[ i ] ) + eps;

        std::array<floatType,dim*nphases> xp = velocity;
        std::array<floatType,dim*nphases> xm = velocity;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,dim*nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( xp ),                std::end( xp ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vp ),                std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( xm ),                std::end( xm ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vm ),                std::end( vm )
        );

        for ( unsigned int j = 0; j < dim * nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( j / dim ) == ( unsigned int )( i / dim ) ){

                unsigned int phase = ( unsigned int )( j / dim );
                unsigned int row   = ( j - phase * dim ) % nphases;
                unsigned int col   = ( i - phase * dim ) % dim;

                BOOST_TEST( dRdV[ sot_dim * phase + dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < dim * nphases; i++ ){

        floatType delta = eps * std::fabs( velocity_dot[ i ] ) + eps;

        std::array<floatType,dim*nphases> xp = velocity_dot;
        std::array<floatType,dim*nphases> xm = velocity_dot;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,dim*nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( xp ),                std::end( xp ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vp ),                std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( xm ),                std::end( xm ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vm ),                std::end( vm )
        );

        for ( unsigned int j = 0; j < dim * nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( j / dim ) == ( unsigned int )( i / dim ) ){

                unsigned int phase = ( unsigned int )( j / dim );
                unsigned int row   = ( j - phase * dim ) % nphases;
                unsigned int col   = ( i - phase * dim ) % dim;

                BOOST_TEST( dRdVDot[ sot_dim * phase + dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < sot_dim * nphases; i++ ){

        floatType delta = eps * std::fabs( velocity_gradient[ i ] ) + eps;

        std::array<floatType,dim*dim*nphases> xp = velocity_gradient;
        std::array<floatType,dim*dim*nphases> xm = velocity_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,dim*nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( xp ),                std::end( xp ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vp ),                std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( xm ),                std::end( xm ),
            std::begin( body_force ),        std::end( body_force ),
            std::begin( vm ),                std::end( vm )
        );

        for ( unsigned int j = 0; j < dim * nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( j / dim ) == ( unsigned int )( i / sot_dim ) ){

                unsigned int phase = ( unsigned int )( j / dim );
                unsigned int row   = ( j - phase * dim ) % nphases;
                unsigned int col   = ( i - phase * dim * dim ) % sot_dim;

                BOOST_TEST( dRdGradV[ sot_dim * dim * phase + sot_dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < dim * nphases; i++ ){

        floatType delta = eps * std::fabs( body_force[ i ] ) + eps;

        std::array<floatType,dim*nphases> xp = body_force;
        std::array<floatType,dim*nphases> xm = body_force;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,dim*nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( xp ),                std::end( xp ),
            std::begin( vp ),                std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumNonDivergence<dim>(
            std::begin( density ),           std::end( density ),
            std::begin( density_dot ),       std::end( density_dot ),
            std::begin( density_gradient ),  std::end( density_gradient ),
            std::begin( velocity ),          std::end( velocity ),
            std::begin( velocity_dot ),      std::end( velocity_dot ),
            std::begin( velocity_gradient ), std::end( velocity_gradient ),
            std::begin( xm ),                std::end( xm ),
            std::begin( vm ),                std::end( vm )
        );

        for ( unsigned int j = 0; j < dim * nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( j / dim ) == ( unsigned int )( i / dim ) ){

                unsigned int phase = ( unsigned int )( j / dim );
                unsigned int row   = ( j - phase * dim ) % nphases;
                unsigned int col   = ( i - phase * dim ) % dim;

                BOOST_TEST( dRdB[ sot_dim * phase + dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_multiphase_computeBalanceOfLinearMomentumDivergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    constexpr unsigned int nphases = 5;

    constexpr unsigned int dim = 3;

    constexpr unsigned int sot_dim = 9;

    std::array<floatType,dim> test_function_gradient = { 0.69646919, 0.28613933, 0.22685145 };

    std::array<floatType,sot_dim*nphases> cauchy_stress = { 0.4809319 , 0.39211752, 0.34317802, 0.72904971, 0.43857224,
                                                            0.0596779 , 0.39804426, 0.73799541, 0.18249173, 0.17545176,
                                                            0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                                            0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                                            0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                                            0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                                            0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                                            0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                                            0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013 };

    std::array<floatType,nphases> volume_fraction = { 0.55131477, 0.71946897, 0.42310646, 0.9807642 , 0.68482974 };

    std::array<floatType,dim*nphases> answer = { -0.34945691, -0.3120474 , -0.16400932, -0.31824658, -0.55913699,
                                                 -0.4683458 , -0.22435795, -0.12580069, -0.17993109, -0.50398514,
                                                 -0.50265385, -0.8776461 , -0.62506454, -0.34963036, -0.44417836 };

    std::array<floatType,dim*nphases> result;

    tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
        std::begin( test_function_gradient ), std::end( test_function_gradient ),
        std::begin( cauchy_stress ),          std::end( cauchy_stress ),
        std::begin( volume_fraction ),        std::end( volume_fraction ),
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    std::fill( std::begin( result ), std::end( result ), 0. );

    std::array<floatType,dim*nphases> dRdPhi;

    std::array<floatType,dim*dim*nphases> dRdGradPsi;

    std::array< floatType, dim * dim * dim * nphases > dRdCauchy;

    tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
        std::begin( test_function_gradient ), std::end( test_function_gradient ),
        std::begin( cauchy_stress ),          std::end( cauchy_stress ),
        std::begin( volume_fraction ),        std::end( volume_fraction ),
        std::begin( result ),     std::end( result ),
        std::begin( dRdGradPsi ), std::end( dRdGradPsi ),
        std::begin( dRdCauchy ),  std::end( dRdCauchy ),
        std::begin( dRdPhi ),     std::end( dRdPhi )
    );

    BOOST_TEST( answer == result, CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < dim; i++ ){

        floatType delta = eps * std::fabs( test_function_gradient[ i ] ) + eps;

        std::array<floatType,dim> xp = test_function_gradient;
        std::array<floatType,dim> xm = test_function_gradient;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,dim*nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( xp ),                     std::end( xp ),
            std::begin( cauchy_stress ),          std::end( cauchy_stress ),
            std::begin( volume_fraction ),        std::end( volume_fraction ),
            std::begin( vp ),                     std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( xm ),                     std::end( xm ),
            std::begin( cauchy_stress ),          std::end( cauchy_stress ),
            std::begin( volume_fraction ),        std::end( volume_fraction ),
            std::begin( vm ),                     std::end( vm )
        );

        for ( unsigned int j = 0; j < dim * nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            BOOST_TEST( dRdGradPsi[ dim * j + i ] == grad );

        }

    }

    for ( unsigned int i = 0; i < sot_dim * nphases; i++ ){

        floatType delta = eps * std::fabs( cauchy_stress[ i ] ) + eps;

        std::array<floatType,sot_dim*nphases> xp = cauchy_stress;
        std::array<floatType,sot_dim*nphases> xm = cauchy_stress;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,dim*nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( test_function_gradient ), std::end( test_function_gradient ),
            std::begin( xp ),                     std::end( xp ),
            std::begin( volume_fraction ),        std::end( volume_fraction ),
            std::begin( vp ),                     std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( test_function_gradient ), std::end( test_function_gradient ),
            std::begin( xm ),                     std::end( xm ),
            std::begin( volume_fraction ),        std::end( volume_fraction ),
            std::begin( vm ),                     std::end( vm )
        );

        for ( unsigned int j = 0; j < dim * nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( j / dim ) == ( unsigned int )( i / sot_dim ) ){

                unsigned int phase = ( unsigned int )( j / dim );
                unsigned int row   = ( j - phase * dim ) % nphases;
                unsigned int col   = ( i - phase * sot_dim ) % sot_dim;

                BOOST_TEST( dRdCauchy[ sot_dim * dim * phase + sot_dim * row + col ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

            }

        }

    }

    for ( unsigned int i = 0; i < nphases; i++ ){

        floatType delta = eps * std::fabs( volume_fraction[ i ] ) + eps;

        std::array<floatType,nphases> xp = volume_fraction;
        std::array<floatType,nphases> xm = volume_fraction;

        xp[ i ] += delta;
        xm[ i ] -= delta;

        std::array<floatType,dim*nphases> vp, vm;

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( test_function_gradient ), std::end( test_function_gradient ),
            std::begin( cauchy_stress ),          std::end( cauchy_stress ),
            std::begin( xp ),                     std::end( xp ),
            std::begin( vp ),                     std::end( vp )
        );

        tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentumDivergence<dim>(
            std::begin( test_function_gradient ), std::end( test_function_gradient ),
            std::begin( cauchy_stress ),          std::end( cauchy_stress ),
            std::begin( xm ),                     std::end( xm ),
            std::begin( vm ),                     std::end( vm )
        );

        for ( unsigned int j = 0; j < dim * nphases; j++ ){

            floatType grad = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

            if ( ( unsigned int )( j / dim ) == ( unsigned int )( i / 1 ) ){

                unsigned int phase = ( unsigned int )( j / dim );
                unsigned int row   = ( j - phase * dim ) % nphases;

                BOOST_TEST( dRdPhi[ dim * phase + row ] == grad );

            }
            else{

                BOOST_TEST( grad == 0. );

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
  typename dt_type, class v_t_in, class v_tp1_in,
  class vDot_t_in, class vDDot_t_in, typename alpha_type, typename beta_type, class vDDot_tp1_out,
  typename dVDDotdV_type
>
void compute_current_acceleration(
    const dt_type &dt,
    const v_t_in &v_t_begin, const v_t_in &v_t_end,
    const v_tp1_in &v_tp1_begin, const v_tp1_in &v_tp1_end,
    const vDot_t_in &vDot_t_begin, const vDot_t_in &vDot_t_end,
    const vDDot_t_in &vDDot_t_begin, const vDDot_t_in &vDDot_t_end,
    const alpha_type &alpha, const beta_type &beta,
    vDDot_tp1_out vDDot_tp1_begin, vDDot_tp1_out vDDot_tp1_end,
    dVDDotdV_type &dVDDotdV
){

    dVDDotdV = 1. / ( dt * dt * alpha * beta );

    for ( unsigned int i = 0; i < ( unsigned int )( v_t_end - v_t_begin ); ++i ){

        *( vDDot_tp1_begin + i ) = ( ( *( v_tp1_begin + i ) ) - ( *( v_t_begin + i ) ) - dt * ( *( vDot_t_begin + i ) ) - ( dt * dt ) * alpha * ( 1 - beta ) * ( *( vDDot_t_begin + i ) ) ) * dVDDotdV;

    }

}

template<
    int dim, int node_count, int nphases,
    class xi_in, typename dt_type, class density_t_in, class density_tp1_in,
    class u_t_in, class u_tp1_in, class umesh_t_in, class umesh_tp1_in, class density_dot_t_in, class u_dot_t_in, class u_ddot_t_in,
    class X_in, class cauchy_stress_iter, class body_force_iter, class volume_fraction_iter, typename alpha_type, typename beta_type, class value_out
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
    const u_dot_t_in &u_dot_t_begin, const u_dot_t_in &u_dot_t_end,
    const u_ddot_t_in &u_ddot_t_begin, const u_ddot_t_in &u_ddot_t_end,
    const X_in &X_begin, const X_in &X_end,
    const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
    const body_force_iter &body_force_begin, const body_force_iter &body_force_end,
    const volume_fraction_iter &volume_fraction_begin, const volume_fraction_iter &volume_fraction_end,
    const alpha_type &alpha, const beta_type &beta, value_out value_begin, value_out value_end
){

    // Update the mesh nodes
    std::array< typename std::iterator_traits<umesh_t_in  >::value_type, dim * node_count > x_t;
    std::array< typename std::iterator_traits<umesh_tp1_in>::value_type, dim * node_count > x_tp1;

    std::transform( X_begin, X_end,   umesh_t_begin,   std::begin( x_t ), std::plus<typename std::iterator_traits<umesh_t_in  >::value_type>( ) );
    std::transform( X_begin, X_end, umesh_tp1_begin, std::begin( x_tp1 ), std::plus<typename std::iterator_traits<umesh_tp1_in>::value_type>( ) );

    // Calculate the current rates of change
    std::array< typename std::iterator_traits<density_tp1_in>::value_type, node_count * nphases > density_dot_tp1;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count * nphases > u_dot_tp1;
    std::array< typename std::iterator_traits<u_tp1_in>::value_type, dim * node_count * nphases > u_ddot_tp1;

    floatType dRhoDotdRho, dUDotdU, dUDDotdU;

    compute_current_rate_of_change(
        dt, density_t_begin, density_t_end, density_tp1_begin, density_tp1_end,
        density_dot_t_begin, density_dot_t_end, alpha,
        std::begin( density_dot_tp1 ), std::end( density_dot_tp1 ),
        dRhoDotdRho
    );

    compute_current_rate_of_change(
        dt, u_t_begin, u_t_end, u_tp1_begin, u_tp1_end,
        u_dot_t_begin, u_dot_t_end, alpha,
        std::begin( u_dot_tp1 ), std::end( u_dot_tp1 ),
        dUDotdU
    );

    compute_current_acceleration(
        dt, u_t_begin, u_t_end, u_tp1_begin, u_tp1_end,
        u_dot_t_begin, u_dot_t_end,
        u_ddot_t_begin, u_ddot_t_end, alpha, beta,
        std::begin( u_ddot_tp1 ), std::end( u_ddot_tp1 ),
        dUDDotdU
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
    > u_dot_tp1_p;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * nphases
    > u_ddot_tp1_p;

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
        xi_begin, xi_end, std::cbegin( u_dot_tp1 ), std::cend( u_dot_tp1 ),
        std::begin( u_dot_tp1_p ), std::end( u_dot_tp1_p )
    );

    e.InterpolateQuantity(
        xi_begin, xi_end, std::cbegin( u_ddot_tp1 ), std::cend( u_ddot_tp1 ),
        std::begin( u_ddot_tp1_p ), std::end( u_ddot_tp1_p )
    );

    // Compute the gradients at the local point

    std::array<
        typename std::iterator_traits<density_tp1_in>::value_type, dim * nphases
    > grad_density_tp1;

    std::array<
        typename std::iterator_traits<u_tp1_in>::value_type, dim * dim * nphases
    > grad_u_dot_tp1;

    e.GetGlobalQuantityGradient(
        xi_begin, xi_end, density_tp1_begin, density_tp1_end,
        std::begin( grad_density_tp1 ), std::end( grad_density_tp1 )
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

            tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentum<dim>(
                density_tp1_p[ 0 ], density_dot_tp1_p[ 0 ],
                std::cbegin( grad_density_tp1 ), std::cend( grad_density_tp1 ),
                std::cbegin( u_dot_tp1_p ), std::cend( u_dot_tp1_p ),
                std::cbegin( u_ddot_tp1_p ), std::cend( u_ddot_tp1_p ),
                std::cbegin( grad_u_dot_tp1 ), std::cend( grad_u_dot_tp1 ),
                body_force_begin, body_force_end,
                cauchy_stress_begin, cauchy_stress_end,
                *volume_fraction_begin,
                Ns[ i ], std::begin( dNdx ) + 3 * i, std::end( dNdx ) + 3 * ( i + 1 ),
                value_begin + i * dim, value_begin + dim * ( i + 1 )
            );

            std::transform(
                value_begin + i * dim, value_begin + dim * ( i + 1 ), value_begin + i * dim,
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

            tardigradeBalanceEquations::balanceOfLinearMomentum::computeBalanceOfLinearMomentum<dim>(
                std::cbegin( density_tp1_p ),     std::cend( density_tp1_p ),
                std::cbegin( density_dot_tp1_p ), std::cend( density_dot_tp1_p ),
                std::cbegin( grad_density_tp1 ),  std::cend( grad_density_tp1 ),
                std::cbegin( u_dot_tp1_p ),       std::cend( u_dot_tp1_p ),
                std::cbegin( u_ddot_tp1_p ),      std::cend( u_ddot_tp1_p ),
                std::cbegin( grad_u_dot_tp1 ),    std::cend( grad_u_dot_tp1 ),
                body_force_begin,                 body_force_end,
                cauchy_stress_begin,              cauchy_stress_end,
                volume_fraction_begin,            volume_fraction_end,
                Ns[ i ], std::begin( dNdx ) + 3 * i, std::end( dNdx ) + 3 * ( i + 1 ),
                value_begin + dim * nphases * i, value_begin + dim * nphases * ( i + 1 )
            );

            std::transform(
                value_begin + dim * nphases * i, value_begin + dim * nphases * ( i + 1 ), value_begin + dim * nphases * i,
                std::bind(
                    std::multiplies< typename std::iterator_traits< value_out >::value_type >( ),
                    std::placeholders::_1,
                    J
                )
            );

        }

    }

}

BOOST_AUTO_TEST_CASE( test_computeBalanceOfLinearMomentum_fea, * boost::unit_test::tolerance( 1e-5 ) ){
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
         0.44683272, -0.96774159,  0.18886376,  0.11357038, -0.68208071,
        -0.69385897,  0.39105906, -0.36246715,  0.38394059,  0.1087665 ,
        -0.22209885,  0.85026498,  0.68333999, -0.28520487, -0.91281707,
        -0.39046385, -0.20362864,  0.40991766,  0.99071696, -0.28817027,
         0.52509563,  0.18635383,  0.3834036 , -0.6977451 
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

    std::array< floatType, 24 > u_dot_t = {
        -0.20224741, -0.5182882 , -0.31308797,  0.02625631,  0.3332491 ,
        -0.78818303, -0.7382101 , -0.35603879,  0.32312867,  0.69301245,
         0.10651469,  0.70890498, -0.23032438, -0.36642421, -0.29147065,
        -0.65783634,  0.65822527, -0.32265831,  0.10474015,  0.15710294,
         0.04306612, -0.99462387,  0.97669084,  0.81068315
    };

    std::array< floatType, 24 > u_ddot_t = {
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

    std::array< floatType, 3 > body_force = {
        -0.20124258, -0.15138628,  0.12443676
    };

    std::array< floatType, 1 > volume_fractions = {
        1.0
    };

    std::array< floatType, 8 * 3 > answer = {
       -0.04454743, -0.03780105, -0.0762556 , -0.008901  , -0.01554515,
        0.00779932, -0.0048873 , -0.0063567 , -0.00087836, -0.03020902,
       -0.02702708, -0.03966113,  0.132318  ,  0.19957351,  0.0004686 ,
        0.0068062 , -0.01650624,  0.07943597, -0.00967059, -0.01466125,
        0.00777241, -0.04839411, -0.02505299, -0.08791766
    };

    std::array< floatType, 8 * 3 > result;

    std::array< floatType, 3 > local_point = {
        -0.7555129, -0.597201 ,  0.6232887
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    floatType beta = 0.67;

    evaluate_at_nodes<3, 8, 1 >(
        std::cbegin( local_point ),      std::cend( local_point ), dt,
        std::cbegin( density_t ),        std::cend( density_t ),
        std::cbegin( density_tp1 ),      std::cend( density_tp1 ),
        std::cbegin( u_t ),              std::cend( u_t ),
        std::cbegin( u_tp1 ),            std::cend( u_tp1 ),
        std::cbegin( umesh_t ),          std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),        std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ),    std::cend( density_dot_t ),
        std::cbegin( u_dot_t ),          std::cend( u_dot_t ),
        std::cbegin( u_ddot_t ),         std::cend( u_ddot_t ),
        std::cbegin( X ),                std::cend( X ),
        std::cbegin( cauchy_stress ),    std::cend( cauchy_stress ),
        std::cbegin( body_force ),       std::cend( body_force ),
        std::cbegin( volume_fractions ), std::cend( volume_fractions ),
        alpha, beta,
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

//    std::fill( std::begin( result ), std::end( result ), 0 );
//
//    std::array< floatType, 8 * 1 * 8 > dCdRho;
//
//    std::array< floatType, 8 * 3 * 8 > dCdU, dCdUMesh;
//
//    evaluate_at_nodes<3, 8, 1>(
//        std::cbegin( local_point ), std::cend( local_point ), dt,
//        std::cbegin( density_t ), std::cend( density_t ),
//        std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//        std::cbegin( u_t ), std::cend( u_t ),
//        std::cbegin( u_tp1 ), std::cend( u_tp1 ),
//        std::cbegin( umesh_t ), std::cend( umesh_t ),
//        std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
//        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//        std::cbegin( v_t ), std::cend( v_t ),
//        std::cbegin( X ), std::cend( X ),
//        alpha,
//        std::begin( result ), std::end( result ),
//        std::begin( dCdRho ), std::end( dCdRho ),
//        std::begin( dCdU ), std::end( dCdU ),
//        std::begin( dCdUMesh ), std::end( dCdUMesh )
//    );
//
//    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );
//
//    floatType eps = 1e-6
//
//    // Check the derivatives w.r.t. the density
//    {
//
//        constexpr unsigned int vardim = 1 * 8;
//        constexpr unsigned int outdim = 8;
//
//        for ( unsigned int i = 0; i < vardim; ++i ){
//
//            floatType delta = eps * std::fabs( density_tp1[ i ] ) + eps;
//
//            std::array< floatType, vardim > xp = density_tp1;
//            std::array< floatType, vardim > xm = density_tp1;
//
//            xp[ i ] += delta;
//            xm[ i ] -= delta;
//
//            std::array< floatType, outdim > vp, vm;
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( xp ), std::cend( xp ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vp ), std::end( vp )
//            );
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( xm ), std::cend( xm ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vm ), std::end( vm )
//            );
//
//            for ( unsigned int j = 0; j < outdim; ++j ){
//
//                BOOST_TEST( dCdRho[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );
//
//            }
//
//        }
//
//    }
//
//    // Check the derivatives w.r.t. the deformation
//    {
//
//        constexpr unsigned int vardim = 3 * 8;
//        constexpr unsigned int outdim = 8;
//
//        for ( unsigned int i = 0; i < vardim; ++i ){
//
//            floatType delta = eps * std::fabs( u_tp1[ i ] ) + eps;
//
//            std::array< floatType, vardim > xp = u_tp1;
//            std::array< floatType, vardim > xm = u_tp1;
//
//            xp[ i ] += delta;
//            xm[ i ] -= delta;
//
//            std::array< floatType, outdim > vp, vm;
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( xp ), std::cend( xp ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vp ), std::end( vp )
//            );
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( xm ), std::cend( xm ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vm ), std::end( vm )
//            );
//
//            for ( unsigned int j = 0; j < outdim; ++j ){
//
//                BOOST_TEST( dCdU[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );
//
//            }
//
//        }
//
//    }
//
//    // Check the derivatives w.r.t. the mesh deformation
//    {
//
//        constexpr unsigned int vardim = 3 * 8;
//        constexpr unsigned int outdim = 8;
//
//        for ( unsigned int i = 0; i < vardim; ++i ){
//
//            floatType delta = eps * std::fabs( umesh_tp1[ i ] ) + eps;
//
//            std::array< floatType, vardim > xp = umesh_tp1;
//            std::array< floatType, vardim > xm = umesh_tp1;
//
//            xp[ i ] += delta;
//            xm[ i ] -= delta;
//
//            std::array< floatType, outdim > vp, vm;
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( xp ), std::cend( xp ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vp ), std::end( vp )
//            );
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( xm ), std::cend( xm ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vm ), std::end( vm )
//            );
//
//            for ( unsigned int j = 0; j < outdim; ++j ){
//
//                BOOST_TEST( dCdUMesh[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );
//
//            }
//
//        }
//
//    }

}

BOOST_AUTO_TEST_CASE( test_computeBalanceOfLinearMomentum_multiphase_fea, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test computing the balance of mass in a finite element context
     */

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
        0.80793821, 0.00742638, 0.55159273, 0.93193215, 0.58217546,
        0.20609573, 0.71775756, 0.37898585, 0.66838395, 0.02931972,
        0.63590036, 0.03219793, 0.74478066, 0.472913  , 0.12175436,
        0.54263593, 0.06677444, 0.65336487, 0.99608633, 0.76939734,
        0.57377411, 0.10263526, 0.69983407, 0.66116787, 0.04909713,
        0.7922993 , 0.51871659, 0.42586769, 0.78818717, 0.41156922,
        0.48102628, 0.18162884
    };

    std::array< floatType, 24 * nphases > u_t = {
        -0.3573622 ,  0.69106599, -0.6261925 , -0.16541788,  0.97806901,
        -0.52680038,  0.83366467,  0.83679494, -0.81740732, -0.07269455,
         0.00443267, -0.3726621 , -0.90532093, -0.51662873, -0.80894072,
        -0.52350019,  0.61558217,  0.78995658, -0.91355422, -0.39610633,
         0.9611644 ,  0.07900965,  0.25261872, -0.98890918, -0.03018111,
         0.97665707, -0.24962895, -0.80592368, -0.07618248,  0.92600893,
        -0.31633877,  0.59784547,  0.59769266, -0.58350341, -0.1132646 ,
         0.43120255, -0.17896043, -0.61798609,  0.93498861,  0.30150073,
         0.7309197 , -0.94951528, -0.46618837,  0.0041422 , -0.86510273,
         0.98606652, -0.52707521, -0.25141564, -0.57197617, -0.78910827,
        -0.53504043, -0.39877973,  0.26888454, -0.43753044, -0.27544648,
        -0.98811431, -0.26856175,  0.06777196, -0.67596833,  0.19486622,
        -0.41369506,  0.26410099, -0.94760679,  0.77518692, -0.96776274,
        -0.74608394,  0.55432492, -0.90820954,  0.42199739,  0.94209228,
         0.74336587,  0.4203233 ,  0.91701949, -0.14037332,  0.74575783,
        -0.28808466,  0.85952731, -0.70244469,  0.88005803,  0.66543239,
         0.69210968, -0.75215398,  0.1929738 , -0.96721504,  0.44236873,
        -0.98452497, -0.83035545, -0.54900318,  0.75024907, -0.27284736,
         0.07991987,  0.13620643, -0.54907328,  0.14429354,  0.32190359,
        -0.40350921
    };

    std::array< floatType, 24 * nphases > u_tp1 = {
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

    std::array< floatType, 8 * nphases > density_dot_t = {
        0.30646975, 0.66526147, 0.11139217, 0.66487245, 0.88785679,
        0.69631127, 0.44032788, 0.43821438, 0.7650961 , 0.565642  ,
        0.08490416, 0.58267109, 0.8148437 , 0.33706638, 0.92757658,
        0.750717  , 0.57406383, 0.75164399, 0.07914896, 0.85938908,
        0.82150411, 0.90987166, 0.1286312 , 0.08178009, 0.13841557,
        0.39937871, 0.42430686, 0.56221838, 0.12224355, 0.2013995 ,
        0.81164435, 0.46798757
    };

    std::array< floatType, 24 * nphases > u_dot_t = {
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

    std::array< floatType, 24 * nphases > u_ddot_t = {
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

    std::array< floatType, 9 * nphases > cauchy_stress = {
         0.41516224,  0.91987827,  0.75340936, -0.06388067,  0.25181302,
        -0.08563654, -0.55410753, -0.24664601, -0.79223154,  0.33305424,
        -0.61593971, -0.04906443,  0.93487321, -0.93666214, -0.6965401 ,
        -0.40284163,  0.88361393,  0.81768359, -0.67599832,  0.96223555,
         0.50149505,  0.07995417,  0.86340577,  0.76121428, -0.21736701,
         0.31268639,  0.29477029, -0.34606363, -0.64121965, -0.06638025,
        -0.47343793, -0.28986975,  0.90828794, -0.07772426,  0.36978293,
        -0.32754021
    };

    std::array< floatType, 3 * nphases > body_force = {
         0.99172216,  0.31753521, -0.60798107, -0.803632  ,  0.88636114,
         0.88955566,  0.24265675, -0.966017  , -0.54893023,  0.60255357,
         0.75091966, -0.09202037
    };

    std::array< floatType, 8 * 3 * nphases > answer = {
         2.18379402e-03,  2.23693322e-03,  8.64026322e-03,  2.95115578e-03,
        -4.60126608e-03, -5.50808623e-04, -1.15704288e-03,  1.76675546e-03,
         1.84501999e-03, -1.27344736e-03, -1.62460370e-03,  6.98276515e-04,
        -1.93341624e-02, -9.47497972e-02,  5.88068109e-02,  2.38255750e-02,
        -2.48680661e-02, -7.89457522e-03,  1.00283464e-02, -1.07706584e-02,
         1.32684268e-02, -9.31664038e-03, -2.00465229e-03,  1.83678079e-02,
        -7.55756741e-02, -4.56246815e-01,  2.54842047e-01, -3.42459433e-02,
         3.34295016e-02,  6.76410035e-02,  3.90534248e-02, -9.94236527e-02,
         9.47385772e-03, -1.04568757e-02,  1.09323634e-02,  2.33915722e-02,
         9.53176864e-03,  6.67086208e-03,  3.68248489e-02,  5.75681671e-03,
        -1.25506023e-02,  2.69786862e-03, -5.05413081e-03,  4.78269889e-03,
         5.47777579e-03, -3.94238611e-03, -5.89916116e-03,  2.52767464e-04,
         6.93212837e-04,  5.41785008e-04,  2.07705954e-03,  7.13548807e-04,
        -1.21016798e-03, -3.96846689e-04, -1.91477965e-04,  2.90835714e-04,
         3.11961353e-04, -2.35216863e-04, -3.72449821e-04,  1.83842332e-04,
         1.45888268e-03, -1.66738878e-02,  1.90872069e-02,  7.33981081e-03,
        -1.09837199e-02, -7.34939161e-03,  2.71795542e-03, -3.28509107e-03,
         1.57179969e-03, -1.51610027e-03, -1.43928359e-03,  4.54989231e-03,
         7.37140960e-03, -8.14805506e-02,  8.19055891e-02,  4.39570307e-03,
        -1.90846151e-02, -1.11280827e-02,  1.08293625e-02, -2.45172652e-02,
        -2.55026631e-03, -7.02598194e-04, -2.25356914e-03,  8.52815886e-03,
         2.98914230e-03,  1.74058169e-03,  8.83879597e-03,  1.69062907e-03,
        -3.75380530e-03, -6.85403624e-04, -8.41656488e-04,  6.99662789e-04,
         8.60814465e-04, -7.10098525e-04, -1.38327810e-03,  2.44575246e-04
    };

    std::array< floatType, nphases > volume_fractions = {
        0.36552062, 0.27422501, 0.11697051, 0.11574454
    };

    std::array< floatType, 8 * 3 * nphases > result;

    std::array< floatType, 3 > local_point = {
        0.9052054 ,  0.61725223, -0.67044128
    };

    floatType dt = 1.3929383711957233;

    floatType alpha = 0.56;

    floatType beta = 0.67;

    evaluate_at_nodes< 3, 8, nphases >(
        std::cbegin( local_point ),      std::cend( local_point ), dt,
        std::cbegin( density_t ),        std::cend( density_t ),
        std::cbegin( density_tp1 ),      std::cend( density_tp1 ),
        std::cbegin( u_t ),              std::cend( u_t ),
        std::cbegin( u_tp1 ),            std::cend( u_tp1 ),
        std::cbegin( umesh_t ),          std::cend( umesh_t ),
        std::cbegin( umesh_tp1 ),        std::cend( umesh_tp1 ),
        std::cbegin( density_dot_t ),    std::cend( density_dot_t ),
        std::cbegin( u_dot_t ),          std::cend( u_dot_t ),
        std::cbegin( u_ddot_t ),         std::cend( u_ddot_t ),
        std::cbegin( X ),                std::cend( X ),
        std::cbegin( cauchy_stress ),    std::cend( cauchy_stress ),
        std::cbegin( body_force ),       std::cend( body_force ),
        std::cbegin( volume_fractions ), std::cend( volume_fractions ),
        alpha, beta,
        std::begin( result ), std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

//    std::fill( std::begin( result ), std::end( result ), 0 );
//
//    std::array< floatType, 8 * 1 * 8 > dCdRho;
//
//    std::array< floatType, 8 * 3 * 8 > dCdU, dCdUMesh;
//
//    evaluate_at_nodes<3, 8, 1>(
//        std::cbegin( local_point ), std::cend( local_point ), dt,
//        std::cbegin( density_t ), std::cend( density_t ),
//        std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//        std::cbegin( u_t ), std::cend( u_t ),
//        std::cbegin( u_tp1 ), std::cend( u_tp1 ),
//        std::cbegin( umesh_t ), std::cend( umesh_t ),
//        std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
//        std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//        std::cbegin( v_t ), std::cend( v_t ),
//        std::cbegin( X ), std::cend( X ),
//        alpha,
//        std::begin( result ), std::end( result ),
//        std::begin( dCdRho ), std::end( dCdRho ),
//        std::begin( dCdU ), std::end( dCdU ),
//        std::begin( dCdUMesh ), std::end( dCdUMesh )
//    );
//
//    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );
//
//    floatType eps = 1e-6
//
//    // Check the derivatives w.r.t. the density
//    {
//
//        constexpr unsigned int vardim = 1 * 8;
//        constexpr unsigned int outdim = 8;
//
//        for ( unsigned int i = 0; i < vardim; ++i ){
//
//            floatType delta = eps * std::fabs( density_tp1[ i ] ) + eps;
//
//            std::array< floatType, vardim > xp = density_tp1;
//            std::array< floatType, vardim > xm = density_tp1;
//
//            xp[ i ] += delta;
//            xm[ i ] -= delta;
//
//            std::array< floatType, outdim > vp, vm;
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( xp ), std::cend( xp ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vp ), std::end( vp )
//            );
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( xm ), std::cend( xm ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vm ), std::end( vm )
//            );
//
//            for ( unsigned int j = 0; j < outdim; ++j ){
//
//                BOOST_TEST( dCdRho[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );
//
//            }
//
//        }
//
//    }
//
//    // Check the derivatives w.r.t. the deformation
//    {
//
//        constexpr unsigned int vardim = 3 * 8;
//        constexpr unsigned int outdim = 8;
//
//        for ( unsigned int i = 0; i < vardim; ++i ){
//
//            floatType delta = eps * std::fabs( u_tp1[ i ] ) + eps;
//
//            std::array< floatType, vardim > xp = u_tp1;
//            std::array< floatType, vardim > xm = u_tp1;
//
//            xp[ i ] += delta;
//            xm[ i ] -= delta;
//
//            std::array< floatType, outdim > vp, vm;
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( xp ), std::cend( xp ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vp ), std::end( vp )
//            );
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( xm ), std::cend( xm ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( umesh_tp1 ), std::cend( umesh_tp1 ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vm ), std::end( vm )
//            );
//
//            for ( unsigned int j = 0; j < outdim; ++j ){
//
//                BOOST_TEST( dCdU[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );
//
//            }
//
//        }
//
//    }
//
//    // Check the derivatives w.r.t. the mesh deformation
//    {
//
//        constexpr unsigned int vardim = 3 * 8;
//        constexpr unsigned int outdim = 8;
//
//        for ( unsigned int i = 0; i < vardim; ++i ){
//
//            floatType delta = eps * std::fabs( umesh_tp1[ i ] ) + eps;
//
//            std::array< floatType, vardim > xp = umesh_tp1;
//            std::array< floatType, vardim > xm = umesh_tp1;
//
//            xp[ i ] += delta;
//            xm[ i ] -= delta;
//
//            std::array< floatType, outdim > vp, vm;
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( xp ), std::cend( xp ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vp ), std::end( vp )
//            );
//
//            evaluate_at_nodes<3, 8, 1>(
//                std::cbegin( local_point ), std::cend( local_point ), dt,
//                std::cbegin( density_t ), std::cend( density_t ),
//                std::cbegin( density_tp1 ), std::cend( density_tp1 ),
//                std::cbegin( u_t ), std::cend( u_t ),
//                std::cbegin( u_tp1 ), std::cend( u_tp1 ),
//                std::cbegin( umesh_t ), std::cend( umesh_t ),
//                std::cbegin( xm ), std::cend( xm ),
//                std::cbegin( density_dot_t ), std::cend( density_dot_t ),
//                std::cbegin( v_t ), std::cend( v_t ),
//                std::cbegin( X ), std::cend( X ),
//                alpha,
//                std::begin( vm ), std::end( vm )
//            );
//
//            for ( unsigned int j = 0; j < outdim; ++j ){
//
//                BOOST_TEST( dCdUMesh[ vardim * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );
//
//            }
//
//        }
//
//    }

}
