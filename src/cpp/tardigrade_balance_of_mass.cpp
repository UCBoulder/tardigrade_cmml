/**
  ******************************************************************************
  * \file tardigrade_balance_of_mass.cpp
  ******************************************************************************
  * The source file for the equations associated with the balance of mass
  ******************************************************************************
  */

#include "tardigrade_balance_of_mass.h"
#include<numeric>

namespace tardigradeBalanceEquations{

    namespace balanceOfMass{

        void computeBalanceOfMass( const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
                                   const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &mass_change_rate ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient: The value of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity: The value of the velocity \f$ v_i \f$
             * \param &velocity_gradient: The value of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             */

            computeBalanceOfMass( density, density_dot,
                                  std::begin( density_gradient ),  std::end( density_gradient ),
                                  std::begin( velocity ),          std::end( velocity ),
                                  std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                  mass_change_rate );

        }

        void computeBalanceOfMass( const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
                                   const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &mass_change_rate, 
                                   floatType   &dCdRho, floatType         &dCdRhoDot, floatVector &dCdGradRho,
                                   floatVector &dCdV,   secondOrderTensor &dCdGradV ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate including the Jacobians
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient: The value of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity: The value of the velocity \f$ v_i \f$
             * \param &velocity_gradient: The value of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             * \param &dCdRho: The derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdRhoDot: The derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dCdGradRho: The derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dCdV: The derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dCdGradV: The derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             */

            computeBalanceOfMass( density,  density_dot,
                                  std::begin( density_gradient ),  std::end( density_gradient ),
                                  std::begin( velocity ),          std::end( velocity ),
                                  std::begin( velocity_gradient ), std::end( velocity_gradient ),
                                  mass_change_rate, dCdRho, dCdRhoDot,
                                  std::begin( dCdGradRho ), std::end( dCdGradRho ),
                                  std::begin( dCdV ),       std::end( dCdV ),
                                  std::begin( dCdGradV ),   std::end( dCdGradV ) );
        }

        template<typename T, class floatVector_iter, class secondOrderTensor_iter>
        void computeBalanceOfMass( const T   &density,  const T         &density_dot,
                                   const floatVector_iter &density_gradient_begin, const floatVector_iter &density_gradient_end,
                                   const floatVector_iter &velocity_begin, const floatVector_iter &velocity_end,
                                   const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                   T &mass_change_rate ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             */

            mass_change_rate = std::inner_product( density_gradient_begin, density_gradient_end,
                                                   velocity_begin, density_dot );

            for ( unsigned int i = 0; i < dim; i++ ){

                mass_change_rate += density * ( *( velocity_gradient_begin + dim * i + i ) );

            }

        }

        template<typename T, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out>
        void computeBalanceOfMass( const T   &density,  const T         &density_dot,
                                   const floatVector_iter &density_gradient_begin,        const floatVector_iter &density_gradient_end,
                                   const floatVector_iter &velocity_begin,                const floatVector_iter &velocity_end,
                                   const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                   T &mass_change_rate,
                                   T &dCdRho, T &dCdRhoDot,
                                   floatVector_iter_out dCdGradRho_begin,     floatVector_iter_out dCdGradRho_end,
                                   floatVector_iter_out dCdV_begin,           floatVector_iter_out dCdV_end,
                                   secondOrderTensor_iter_out dCdGradV_begin, secondOrderTensor_iter_out dCdGradV_end ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             * \param &dCdRho: The derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdRhoDot: The derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dCdGradRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dCdGradRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dCdV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dCdV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dCdGradV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &dCdGradV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             */

            computeBalanceOfMass( density, density_dot,
                                  density_gradient_begin, density_gradient_end,
                                  velocity_begin, velocity_end,
                                  velocity_gradient_begin, velocity_gradient_end, mass_change_rate );

            dCdRho = 0;

            dCdRhoDot = 1;

            std::copy( velocity_begin, velocity_end, dCdGradRho_begin );

            std::copy( density_gradient_begin, density_gradient_end, dCdV_begin );

            std::fill( dCdGradV_begin, dCdGradV_end, 0 );

            for ( unsigned int i = 0; i < dim; i++ ){

                dCdRho += ( *( velocity_gradient_begin + dim * i + i ) );

                ( *( dCdGradV_begin + dim * i + i ) ) = density;

            }

        }

        template<class scalarArray_iter, class floatVectorArray_iter, class secondOrderTensorArray_iter, class scalarArray_iter_out>
        void computeBalanceOfMass( const scalarArray_iter &density_begin,                      const scalarArray_iter &density_end,
                                   const scalarArray_iter &density_dot_begin,                  const scalarArray_iter &density_dot_end,
                                   const floatVectorArray_iter &density_gradient_begin,        const floatVectorArray_iter &density_gradient_end,
                                   const floatVectorArray_iter &velocity_begin,                const floatVectorArray_iter &velocity_end,
                                   const secondOrderTensorArray_iter &velocity_gradient_begin, const secondOrderTensorArray_iter &velocity_gradient_end,
                                   scalarArray_iter_out mass_change_rate_begin,                scalarArray_iter_out mass_change_rate_end ){
            /*!
             * Compute the balance of mass for a multi-phase continuum returning the values of the mass-change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density_begin: The starting iterator of the density \f$ \rho \f$
             * \param &density_end: The stopping iterator of the density \f$ \rho \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &mass_change_rate_begin: The starting iterator of the rate of change of the mass \f$ c \f$
             * \param &mass_change_rate_end: The stopping iterator of the rate of change of the mass \f$ c \f$
             */

            unsigned int phase;

            for ( auto rho = density_begin; rho != density_end; rho++ ){

                phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfMass( *( density_begin + phase ), *( density_dot_begin + phase ),
                                      density_gradient_begin + dim * phase,        density_gradient_begin + dim * ( phase + 1 ),
                                      velocity_begin + dim * phase,                velocity_begin + dim * ( phase + 1 ),
                                      velocity_gradient_begin + sot_dim * phase,   velocity_gradient_begin + sot_dim * ( phase + 1 ),
                                      *( mass_change_rate_begin + phase ) );

            }

        }

        template<class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class scalarArray_iter_out, class floatVector_iter_out, class secondOrderTensor_iter_out>
        void computeBalanceOfMass( const scalarArray_iter &density_begin,                 const scalarArray_iter &density_end,
                                   const scalarArray_iter &density_dot_begin,             const scalarArray_iter &density_dot_end,
                                   const floatVector_iter &density_gradient_begin,        const floatVector_iter &density_gradient_end,
                                   const floatVector_iter &velocity_begin,                const floatVector_iter &velocity_end,
                                   const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                   scalarArray_iter_out mass_change_rate_begin, scalarArray_iter_out mass_change_rate_end,
                                   scalarArray_iter_out dCdRho_begin,           scalarArray_iter_out dCdRho_end,
                                   scalarArray_iter_out dCdRhoDot_begin,        scalarArray_iter_out dCdRhoDot_end,
                                   floatVector_iter_out dCdGradRho_begin,       floatVector_iter_out dCdGradRho_end,
                                   floatVector_iter_out dCdV_begin,             floatVector_iter_out dCdV_end,
                                   secondOrderTensor_iter_out dCdGradV_begin,   secondOrderTensor_iter_out dCdGradV_end ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density_begin: The starting iterator of the value of the density \f$ \rho \f$
             * \param &density_end: The stopping iterator of the value of the density \f$ \rho \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &mass_change_rate_begin: The starting iterator of the rate of change of the mass \f$ c \f$
             * \param &mass_change_rate_end: The stopping iterator of the rate of change of the mass \f$ c \f$
             * \param &dCdRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdRhoDot_begin: The starting iterator of the derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dCdRhoDot_end: The stopping iterator of the derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dCdGradRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dCdGradRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dCdV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dCdV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dCdGradV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &dCdGradV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             */

            unsigned int phase;

            for ( auto rho = density_begin; rho != density_end; rho++ ){

                phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfMass( *( density_begin + phase ), *( density_dot_begin + phase ),
                                      density_gradient_begin + dim * phase,        density_gradient_begin + dim * ( phase + 1 ),
                                      velocity_begin + dim * phase,                velocity_begin + dim * ( phase + 1 ),
                                      velocity_gradient_begin + sot_dim * phase, velocity_gradient_begin + sot_dim * ( phase + 1 ),
                                      *( mass_change_rate_begin + phase ), *( dCdRho_begin + phase ), *( dCdRhoDot_begin + phase ),
                                      dCdGradRho_begin + dim * phase,   dCdGradRho_begin + dim * ( phase + 1 ),
                                      dCdV_begin + dim * phase,         dCdV_begin + dim * ( phase + 1 ),
                                      dCdGradV_begin + sot_dim * phase, dCdGradV_begin + sot_dim * ( phase + 1 ) );

            }

        }

    }

}
