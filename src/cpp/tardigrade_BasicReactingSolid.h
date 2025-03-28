/*!
=====================================================================
|                 tardigrade_BasicReactingSolid.h                   |
=====================================================================
| A basic reacting solid.                                           |
---------------------------------------------------------------------
| The output material vector is comprised of the following          |
| responses.                                                        |
|                                                                   |
| Stress response          : Linear, isotropic elasticity           |
| Predicted energy         : Linearly related to the temperature    |
|     via the specific heat.                                        |
| Mass generation rate     : Zero                                   |
| Body force               : Zero                                   |
| Interphasic force        : Zero                                   |
| Heat flux                : Linearly related to the temperature    |
| Internal heat generation : Zero                                   |
| Interphasic heat transfer: Zero                                   |
| trace mass change velocity gradient: Zero                         |
|                                                                   |
| The class can be used as the base class for more advanced         |
| behaviors                                                         |
=====================================================================
*/

#ifndef TARDIGRADE_GENERALREACTINGSOLID_H
#define TARDIGRADE_GENERALREACTINGSOLID_H

#include "tardigrade_error_tools.h"
#include "tardigrade_cmml.h"

namespace tardigradeCMML{

    /*!
     * A general reacting continuum solid
     */
    template<
        typename time_type,
        class current_dof_iter,
        class previous_dof_iter,
        class parameter_iter,
        class sdvs_iter,
        class result_iter,
        class jacobian_iter,
        class additional_iter
    >
    class BasicReactingSolid : public tardigradeCMML::CMMLMaterial<
        time_type, current_dof_iter, previous_dof_iter, parameter_iter,
        sdvs_iter, result_iter, jacobian_iter, additional_iter
    >{

        public:

            using tardigradeCMML::CMMLMaterialBase::setName;

            BasicReactingSolid( ){
                /*!
                 * A general reacting solid.
                 * 
                 * Stress response          : Linear, isotropic elasticity
                 * Predicted energy         : Linearly related to the temperature
                 *     via the specific heat.
                 * Mass generation rate     : Zero
                 * Body force               : Zero
                 * Interphasic force        : Zero
                 * Heat flux                : Linearly related to the temperature
                 * Internal heat generation : Zero
                 * Interphasic heat transfer: Zero
                 * trace mass change velocity gradient: Zero
                 *
                 * The class can be used as the base class for more advanced
                 * behaviors
                 * 
                 */

                setName( "BasicReactingSolid" );

            }

            virtual int evaluate_model(
                const time_type         &current_time,       const time_type &dt,
                const current_dof_iter  &current_dof_begin,  const current_dof_iter  &current_dof_end,
                const previous_dof_iter &previous_dof_begin, const previous_dof_iter &previous_dof_end,
                const parameter_iter    &parameters_begin,   const parameter_iter    &parameters_end,
                sdvs_iter               sdvs_begin,          sdvs_iter               sdvs_end,
                result_iter             result_begin,        result_iter             result_end,
                std::string             &output_message
            ){
            /*!
             * Evaluate the material model given the incoming data
             *
             * Returns:
             *     0: No errors. Solution converged.
             *     1: Convergence Error. Request timestep cutback.
             *     2: Fatal errors encountered. Terminate evaluation.
             *
             * \param &current_time: The current time of the simulation
             * \param &dt:           The change in time of the simulation
             * \param &current_dof_begin: The starting iterator of the current values of the degrees of freedom
             * \param &current_dof_end: The stopping iterator of the current values of the degrees of freedom
             * \param &previous_dof_begin: The starting iterator of the previous values of the degrees of freedom
             * \param &previous_dof_end: The stopping iterator of the previous values of the degrees of freedom
             * \param &parameters_begin: The starting iterator of the model parameters
             * \param &parameters_end: The stopping iterator of the model parameters
             * \param &sdvs_begin: The starting iterator of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param &sdvs_end: The stopping iterator of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param &result_begin: The starting iterator of the result vector.
             * \param &result_end: The stopping iterator of the result vector.
             * \param &output_message: An output string containing messages from the code.
             */

                TARDIGRADE_ERROR_TOOLS_CHECK( _expected_material_size == ( unsigned int )( result_end - result_begin ), "The output material size must be equal to " + std::to_string( _expected_material_size ) )

                output_message = "evaluate_model is not implemented";
                return 2;

            }

            virtual int evaluate_model(
                const time_type         &current_time,       const time_type &dt,
                const current_dof_iter  &current_dof_begin,  const current_dof_iter  &current_dof_end,
                const previous_dof_iter &previous_dof_begin, const previous_dof_iter &previous_dof_end,
                const parameter_iter    &parameters_begin,   const parameter_iter    &parameters_end,
                sdvs_iter               sdvs_begin,          sdvs_iter               sdvs_end,
                result_iter             result_begin,        result_iter             result_end,
                jacobian_iter           jacobian_begin,      jacobian_iter           jacobian_end,
                additional_iter         additional_begin,    additional_iter         additional_end,
                std::string             &output_message
            ){
            /*!
             * Evaluate the material model given the incoming data
             *
             * Returns:
             *     0: No errors. Solution converged.
             *     1: Convergence Error. Request timestep cutback.
             *     2: Fatal errors encountered. Terminate evaluation.
             *
             * \param &current_time: The current time of the simulation
             * \param &dt:           The change in time of the simulation
             * \param &current_dof_begin: The starting iterator of the current values of the degrees of freedom
             * \param &current_dof_end: The stopping iterator of the current values of the degrees of freedom
             * \param &previous_dof_begin: The starting iterator of the previous values of the degrees of freedom
             * \param &previous_dof_end: The stopping iterator of the previous values of the degrees of freedom
             * \param &parameters_begin: The starting iterator of the model parameters
             * \param &parameters_end: The stopping iterator of the model parameters
             * \param &sdvs_begin: The starting iterator of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param &sdvs_end: The stopping iterator of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param &result_begin: The starting iterator of the result vector.
             * \param &result_end: The stopping iterator of the result vector.
             * \param &jacobian_begin: The starting iterator of the jacobian of the result w.r.t. the dof.
             * \param &jacobian_end: The stopping iterator of the jacobian of the result w.r.t. the dof.
             * \param &additional_begin: The starting iterator of the additional information returned by the model
             * \param &additional_end: The stopping iterator of the additional information returned by the model
             * \param &output_message: An output string containing messages from the code.
             */

                output_message = "evaluate_model with additional information is not implemented";
                return 2;

            }

            unsigned int _expected_material_size = 23;

    };

    REGISTER_MATERIAL( tardigradeCMML::BasicReactingSolid );

}

#include "tardigrade_BasicReactingSolid.cpp"

#endif
