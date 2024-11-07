/**
  ******************************************************************************
  * \file tardigrade_balance_of_mass.h
  ******************************************************************************
  * The header file for the equations associated with the balance of mass
  ******************************************************************************
  */

#ifndef TARDIGRADE_BALANCE_OF_MASS_H
#define TARDIGRADE_BALANCE_OF_MASS_H

#include<array>

#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace balanceOfMass{

        constexpr unsigned int dim = 3; //!< Set the dimension as 3D by default

        constexpr unsigned int sot_dim = dim * dim; //!< Set the dimensions of a standard second order tensor

        typedef double floatType; //!< Define the float type as a double

        typedef std::array< floatType, dim > floatVector; //!< Define a standard vector

        typedef std::array< floatType, sot_dim > secondOrderTensor; //!< Define a standard second-order tensor

        void computeBalanceOfMass( const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
                                   const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &mass_change_rate );

        void computeBalanceOfMass( const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
                                   const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &mass_change_rate,
                                   floatType   &dCdRho, floatType         &dCdRhoDot, floatVector &dCdGradRho,
                                   floatVector &dCdV,   secondOrderTensor &dCdGradV );

        template<class floatVector_iter, class secondOrderTensor_iter>
        void computeBalanceOfMass( const floatType   &density,  const floatType         &density_dot,
                                   const floatVector_iter &density_gradient_begin, const floatVector_iter &density_gradient_end,
                                   const floatVector_iter &velocity_begin, const floatVector_iter &velocity_end,
                                   const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                   floatType &mass_change_rate );

        template<class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out>
        void computeBalanceOfMass( const floatType   &density,  const floatType         &density_dot,
                                   const floatVector_iter &density_gradient_begin,        const floatVector_iter &density_gradient_end,
                                   const floatVector_iter &velocity_begin,                const floatVector_iter &velocity_end,
                                   const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                   floatType &mass_change_rate,
                                   floatType &dCdRho, floatType &dCdRhoDot,
                                   floatVector_iter_out dCdGradRho_begin,     floatVector_iter_out dCdGradRho_end,
                                   floatVector_iter_out dCdV_begin,           floatVector_iter_out dCdV_end,
                                   secondOrderTensor_iter_out dCdGradV_begin, secondOrderTensor_iter_out dCdGradV_end );

    }

}

#include "tardigrade_balance_of_mass.cpp"

#endif
