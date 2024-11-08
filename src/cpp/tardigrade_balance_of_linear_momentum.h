/**
  ******************************************************************************
  * \file tardigrade_balance_of_linear_momentum.h
  ******************************************************************************
  * The header file for the equations associated with the balance of linear
  * momentum
  ******************************************************************************
  */

#ifndef TARDIGRADE_BALANCE_OF_LINEAR_MOMENTUM_H
#define TARDIGRADE_BALANCE_OF_LINEAR_MOMENTUM_H

#include<array>

#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace balanceOfLinearMomentum{

        constexpr unsigned int dim = 3; //!< Set the dimension as 3D by default

        constexpr unsigned int sot_dim = dim * dim; //!< Set the dimensions of a standard second order tensor

        typedef double floatType; //!< Define the float type as a double

        typedef std::array< floatType, dim > floatVector; //!< Define a standard vector

        typedef std::array< floatType, sot_dim > secondOrderTensor; //!< Define a standard second-order tensor

        template<class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out>
        void computeBalanceOfLinearMomentumNonDivergence( const floatType &density, const floatType &density_dot,
                                                          const floatVector_iter &density_gradient_begin,        const floatVector_iter &density_gradient_end,
                                                          const floatVector_iter &velocity_begin,                const floatVector_iter &velocity_end,
                                                          const floatVector_iter &velocity_dot_begin,            const floatVector_iter &velocity_dot_end,
                                                          const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                                          const floatVector_iter &body_force_begin,              const floatVector_iter &body_force_end,
                                                          floatVector_iter_out result_begin,                     floatVector_iter_out result_end );

    }

}

#include "tardigrade_balance_of_linear_momentum.cpp"

#endif
