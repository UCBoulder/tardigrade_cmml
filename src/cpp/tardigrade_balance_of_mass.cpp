/**
  ******************************************************************************
  * \file tardigrade_balance_of_mass.cpp
  ******************************************************************************
  * The source file for the equations associated with the balance of mass
  ******************************************************************************
  */

#include "tardigrade_balance_of_mass.h"

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
             * \param &density_dot: The value of the partial derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient: The value of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity: The value of the velocity \f$ v_i \f$
             * \param &velocity_gradient: The value of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             */

            mass_change_rate = density_dot;

            for ( unsigned int i = 0; i < dim; i++ ){

                mass_change_rate += density_gradient[ i ] * velocity[ i ] + density * velocity_gradient[ dim * i + i ];

            }

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
             * \param &density_dot: The value of the partial derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient: The value of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity: The value of the velocity \f$ v_i \f$
             * \param &velocity_gradient: The value of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             */

            computeBalanceOfMass( density, density_dot, density_gradient, velocity, velocity_gradient, mass_change_rate );

            dCdRho = 0;

            dCdRhoDot = 1;

            dCdGradRho = velocity;

            dCdV = density_gradient;

            std::fill( std::begin( dCdGradV ), std::end( dCdGradV ), 0 );

            for ( unsigned int i = 0; i < dim; i++ ){

                dCdRho += velocity_gradient[ dim * i + i ];

                dCdGradV[ dim * i + i ] = density;

            }

        }

    }

}

