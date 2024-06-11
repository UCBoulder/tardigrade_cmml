/**
  ******************************************************************************
  * \file tardigrade_balance_of_mass.h
  ******************************************************************************
  * The header file for the equations associated with the balance of mass
  ******************************************************************************
  */

#ifndef TARDIGRADE_MASS_CHANGE_DEFORMATION_H
#define TARDIGRADE_MASS_CHANGE_DEFORMATION_H

#include<vector>
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

    }

}

#include "tardigrade_balance_of_mass.cpp"

#endif
