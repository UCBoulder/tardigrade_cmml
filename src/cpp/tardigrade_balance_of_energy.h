/**
  ******************************************************************************
  * \file tardigrade_balance_of_energy.h
  ******************************************************************************
  * The header file for the equations associated with the balance of energy
  ******************************************************************************
  */

#ifndef TARDIGRADE_BALANCE_OF_ENERGY_H
#define TARDIGRADE_BALANCE_OF_ENERGY_H

#include<array>

#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace balanceOfEnergy{

        constexpr unsigned int dim = 3; //!< Set the dimension as 3D by default

        constexpr unsigned int sot_dim = dim * dim; //!< Set the dimensions of a standard second order tensor

        typedef double floatType; //!< Define the float type as a double

        typedef std::array< floatType, dim > floatVector; //!< Define a standard vector

        typedef std::array< floatType, sot_dim > secondOrderTensor; //!< Define a standard second-order tensor

        template<class floatVector_iter, class secondOrderTensor_iter>
        void computeBalanceOfEnergyNonDivergence( const floatType &density, const floatType &density_dot,
                                                  const floatVector_iter &density_gradient_begin, const floatVector_iter &density_gradient_end,
                                                  const floatType &internal_energy, const floatType &internal_energy_dot,
                                                  const floatVector_iter &internal_energy_gradient_begin, const floatVector_iter &internal_energy_gradient_end,
                                                  const floatVector_iter &velocity_begin, const floatVector_iter &velocity_end,
                                                  const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                                  const secondOrderTensor_iter &cauchy_stress_begin, const secondOrderTensor_iter &cauchy_stress_end,
                                                  const floatType &volume_fraction,
                                                  const floatType &internal_heat_generation,
                                                  const floatVector_iter &net_interphase_force_begin, const floatVector_iter &net_interphase_force_end,
                                                  floatType &result );

    }

}

#include "tardigrade_balance_of_energy.cpp"

#endif
