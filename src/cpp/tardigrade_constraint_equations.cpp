/**
  ******************************************************************************
  * \file tardigrade_constraint_equations.cpp
  ******************************************************************************
  * The source file for the constraint equations associated with the balance
  * equations
  ******************************************************************************
  */

#include "tardigrade_constraint_equations.h"
#include<numeric>
#include<algorithm>
#include<functional>
#include<array>

namespace tardigradeBalanceEquations{

    namespace constraintEquations{

        template<
            int predicted_internal_energy_index,
            typename internal_energy_type,
            class material_response_iter,
            typename result_type
        >
        void computeInternalEnergyConstraint(
            const internal_energy_type &internal_energy,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            result_type &result
        ){
            /*!
             * Compute the value of the constraint on the internal energy.
             * 
             * Traditionally, the internal energy is related to the temperature via
             * 
             * \f$ e = c \theta \f$
             * 
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             * 
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             * 
             * \param &internal_energy: The internal energy degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &result: The resulting error between the internal energy DOF and the material response
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(
                predicted_internal_energy_index < ( unsigned int )( material_response_end - material_response_begin ),
                "The index for the predicted internal energy is out of range for the material response"
            )

            result = *( material_response_begin + predicted_internal_energy_index ) - internal_energy;

        }

        template<
            int predicted_internal_energy_index,
            class internal_energy_iter,
            class material_response_iter,
            class result_iter
        >
        void computeInternalEnergyConstraint(
            const internal_energy_iter &internal_energy_begin,
            const internal_energy_iter &internal_energy_end,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            result_iter result_begin,
            result_iter result_end
        ){
            /*!
             * Compute the value of the constraint on the internal energy for a multiphase problem.
             * 
             * Traditionally, the internal energy is related to the temperature via
             * 
             * \f$ e = c \theta \f$
             * 
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             * 
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             * 
             * \param &internal_energy_begin: The starting iterator of the internal energy degree of freedom
             * \param &internal_energy_end: The stopping iterator of the internal energy degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &result_begin: The starting iterator of the resulting error between the internal energy DOF and the material response
             * \param &result_end: The starting iterator of the resulting error between the internal energy DOF and the material response
             */

            const unsigned int nphases = ( unsigned int )( result_end - result_begin );
            const unsigned int material_response_size = ( unsigned int )( material_response_end - material_response_begin ) / nphases;

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( internal_energy_end - internal_energy_begin ),
                "The number of internal energy values must be equal to the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size == ( unsigned int )( material_response_end - material_response_begin ),
                "The material response vector must be a scalar multiple of the number of phases"
            )

            for ( auto v = std::pair< unsigned int, result_iter >( 0, result_begin ); v.second != result_end; ++v.first, ++v.second ){

                computeInternalEnergyConstraint<
                    predicted_internal_energy_index
                >
                (
                    *( internal_energy_begin + v.first ),
                    material_response_begin + material_response_size * v.first,
                    material_response_begin + material_response_size * ( v.first + 1 ),
                    *v.second
                );

            }

        }

    }

}
