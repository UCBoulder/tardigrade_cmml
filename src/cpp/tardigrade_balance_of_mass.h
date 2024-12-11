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

        constexpr unsigned int global_dim = 3; //!< Set the dimension as 3D by default

        constexpr unsigned int global_sot_dim = global_dim * global_dim; //!< Set the dimensions of a standard second order tensor

        typedef double floatType; //!< Define the float type as a double

        typedef std::array< floatType, global_dim > floatVector; //!< Define a standard vector

        typedef std::array< floatType, global_sot_dim > secondOrderTensor; //!< Define a standard second-order tensor

        void computeBalanceOfMass( const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
                                   const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &mass_change_rate );

        void computeBalanceOfMass( const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
                                   const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &mass_change_rate,
                                   floatType   &dCdRho, floatType         &dCdRhoDot, floatVector &dCdGradRho,
                                   floatVector &dCdV,   secondOrderTensor &dCdGradV );

        template<int dim, typename T, class densityGradient_iter, class velocity_iter, class velocityGradient_iter>
        void computeBalanceOfMass( const T   &density,  const T         &density_dot,
                                   const densityGradient_iter &density_gradient_begin, const densityGradient_iter &density_gradient_end,
                                   const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
                                   const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
                                   T &mass_change_rate );

        template<int dim, typename T, class densityGradient_iter, class velocity_iter, class velocityGradient_iter, class dCdGradRho_iter_out, class dCdV_iter_out, class dCdGradV_iter_out>
        void computeBalanceOfMass( const T   &density,  const T         &density_dot,
                                   const densityGradient_iter &density_gradient_begin,    const densityGradient_iter &density_gradient_end,
                                   const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
                                   const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
                                   T &mass_change_rate,
                                   T &dCdRho, T &dCdRhoDot,
                                   dCdGradRho_iter_out dCdGradRho_begin, dCdGradRho_iter_out dCdGradRho_end,
                                   dCdV_iter_out dCdV_begin,             dCdV_iter_out dCdV_end,
                                   dCdGradV_iter_out dCdGradV_begin,     dCdGradV_iter_out dCdGradV_end );

        template<int dim, class density_iter, class scalarArray_iter, class floatVectorArray_iter, class secondOrderTensorArray_iter, class scalarArray_iter_out>
        void computeBalanceOfMass( const density_iter &density_begin,                          const density_iter &density_end,
                                   const scalarArray_iter &density_dot_begin,                  const scalarArray_iter &density_dot_end,
                                   const floatVectorArray_iter &density_gradient_begin,        const floatVectorArray_iter &density_gradient_end,
                                   const floatVectorArray_iter &velocity_begin,                const floatVectorArray_iter &velocity_end,
                                   const secondOrderTensorArray_iter &velocity_gradient_begin, const secondOrderTensorArray_iter &velocity_gradient_end,
                                   scalarArray_iter_out mass_change_rate_start,                scalarArray_iter_out mass_change_rate_stop );

        template<int dim, class density_iter, class scalarArray_iter, class densityGradient_iter, class velocity_iter, class velocityGradient_iter, class scalarArray_iter_out, class dCdGradRho_iter_out, class dCdV_iter_out, class dCdGradV_iter_out>
        void computeBalanceOfMass( const density_iter &density_begin,                    const density_iter &density_end,
                                   const scalarArray_iter &density_dot_begin,            const scalarArray_iter &density_dot_end,
                                   const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
                                   const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
                                   const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
                                   scalarArray_iter_out mass_change_rate_begin, scalarArray_iter_out mass_change_rate_end,
                                   scalarArray_iter_out dCdRho_begin,    scalarArray_iter_out dCdRho_end,
                                   scalarArray_iter_out dCdRhoDot_begin, scalarArray_iter_out dCdRhoDot_end,
                                   dCdGradRho_iter_out dCdGradRho_begin, dCdGradRho_iter_out dCdGradRho_end,
                                   dCdV_iter_out dCdV_begin,             dCdV_iter_out dCdV_end,
                                   dCdGradV_iter_out dCdGradV_begin,     dCdGradV_iter_out dCdGradV_end );

    }

}

#include "tardigrade_balance_of_mass.cpp"

#endif
