/*!
=====================================================================
|                 tardigrade_LinearTestMaterial.cpp                 |
=====================================================================
| A material which defines a response as being a linear product     |
| between the parameters and the degrees of freedom. To be used for |
| testing purposes.                                                 |
=====================================================================
*/

#include "tardigrade_LinearTestMaterial.h"

namespace tardigradeCMML {

    namespace LinearTestMaterial {

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
         * \param *current_dof_begin: A pointer to the starting element of the current values of the degrees of
         * freedom \param *previous_dof_begin: A pointer to the starting element of the previous values of the
         * degrees of freedom \param dof_size: The size of the dof arrays \param *parameters_begin: A pointer to the
         * starting element of the model parameters \param parameters_size: The size of the parameter vector \param
         * *sdvs_begin: A pointer to the starting element of the model state variables. These are initialized to the
         *     previous converged values and should be updated by the model
         * \param sdvs_size: The size of the SDVS vector
         * \param *result_begin: A pointer to the starting element of the result vector.
         * \param result_size: The size of the result
         * \param &output_message: An output string containing messages from the code.
         */
        int LinearTestMaterial::evaluate_model(const time_type &current_time, const time_type &dt,
                                               const current_dof_type  *current_dof_begin,
                                               const previous_dof_type *previous_dof_begin, const unsigned int dof_size,
                                               const parameter_type *parameters_begin,
                                               const unsigned int parameters_size, sdvs_type *sdvs_begin,
                                               const unsigned int sdvs_size, result_type *result_begin,
                                               const unsigned int result_size, std::string &output_message) {
            try {
                TARDIGRADE_ERROR_TOOLS_CHECK(parameters_size == (dof_size * result_size),
                                             "The parameter vector has a size of " + std::to_string(parameters_size) +
                                                 " but must have a size of " + std::to_string(dof_size * result_size))

                std::fill(result_begin, result_begin + result_size, result_type());
                for (unsigned int i = 0; i < result_size; ++i) {
                    for (unsigned int j = 0; j < dof_size; ++j) {
                        *(result_begin + i) += (*(parameters_begin + dof_size * i + j)) * (*(current_dof_begin + j));
                    }
                }

            } catch (tardigradeHydra::convergence_error &e) {
                return 1;

            } catch (std::exception &e) {
                tardigradeErrorTools::captureNestedExceptions(e, output_message);

                return 2;
            }

            return 0;
        }

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
         * \param *current_dof_begin: A pointer to the starting element of the current values of the degrees of
         * freedom \param *previous_dof_begin: A pointer to the starting element of the previous values of the
         * degrees of freedom \param dof_size: The size of the dof arrays \param *parameters_begin: A pointer to the
         * starting element of the model parameters \param parameters_size: The size of the parameter vector \param
         * *sdvs_begin: A pointer to the starting element of the model state variables. These are initialized to the
         *     previous converged values and should be updated by the model
         * \param sdvs_size: The size of the SDVS vector
         * \param *result_begin: A pointer to the starting element of the result vector.
         * \param result_size: The size of the result
         * \param *jacobian_begin: A pointer to the starting element of the jacobian of the result w.r.t. the dof.
         * (must have a size of result_size * dof_size) \param *additional_begin: A pointer to the starting element
         * of the additional information returned by the model \param additional_size: The size of the additional
         * information array \param &output_message: An output string containing messages from the code.
         */
        int LinearTestMaterial::evaluate_model(const time_type &current_time, const time_type &dt,
                                               const current_dof_type  *current_dof_begin,
                                               const previous_dof_type *previous_dof_begin, const unsigned int dof_size,
                                               const parameter_type *parameters_begin,
                                               const unsigned int parameters_size, sdvs_type *sdvs_begin,
                                               const unsigned int sdvs_size, result_type *result_begin,
                                               const unsigned int result_size, jacobian_type *jacobian_begin,
                                               additional_type *additional_begin, const unsigned int additional_size,
                                               std::string &output_message) {
            try {
                TARDIGRADE_ERROR_TOOLS_CHECK(parameters_size == (dof_size * result_size),
                                             "The parameter vector has a size of " + std::to_string(parameters_size) +
                                                 " but must have a size of " + std::to_string(dof_size * result_size))

                std::fill(result_begin, result_begin + result_size, result_type());
                for (unsigned int i = 0; i < result_size; ++i) {
                    for (unsigned int j = 0; j < dof_size; ++j) {
                        *(result_begin + i) += (*(parameters_begin + dof_size * i + j)) * (*(current_dof_begin + j));
                    }
                }

                std::copy(parameters_begin, parameters_begin + parameters_size, jacobian_begin);

            } catch (tardigradeHydra::convergence_error &e) {
                return 1;

            } catch (std::exception &e) {
                tardigradeErrorTools::captureNestedExceptions(e, output_message);

                return 2;
            }

            return 0;
        }

    }  // namespace LinearTestMaterial

}  // namespace tardigradeCMML
