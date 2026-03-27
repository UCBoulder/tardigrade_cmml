/*!
=====================================================================
|                      tardigrade_BasicSolid.h                      |
=====================================================================
| A basic solid.                                                    |
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

#ifndef TARDIGRADE_BASICSOLID_H
#define TARDIGRADE_BASICSOLID_H

#include <vector>

#include "tardigrade_error_tools.h"
#define TARDIGRADE_CMML_SKIP_HEADER_ONLY_MATERIALS_INCLUDE
#include "tardigrade_cmml.h"
#undef TARDIGRADE_CMML_SKIP_HEADER_ONLY_MATERIALS_INCLUDE
#include "tardigrade_hydraFourierHeatConduction.h"
#include "tardigrade_hydraLinearElasticity.h"
#include "tardigrade_hydraLinearInternalEnergy.h"

namespace tardigradeCMML {

    namespace BasicSolid {

        /*!
         * The basic solid material model hydra class
         */
        class basicSolidHydra : public tardigradeHydra::hydraBase {
           public:
            using tardigradeHydra::hydraBase::hydraBase;
            using tardigradeHydra::hydraBase::setResidualClasses;

            /*!
             * Set the parameters for the stress calculation
             *
             * \param &value_start: The starting value of the parameter iterator
             * \param &value_end: The stopping value of the parameter iterator
             */
            template <class parameter_iter>
            void setStressParameters(const parameter_iter &value_start, const parameter_iter &value_end) {
                _stress_parameters = std::vector<double>(value_start, value_end);
            }

            /*!
             * Set the parameters for the internal energy
             *
             * \param &value_start: The starting value of the parameter iterator
             * \param &value_end: The stopping value of the parameter iterator
             */
            template <class parameter_iter>
            void setInternalEnergyParameters(const parameter_iter &value_start, const parameter_iter &value_end) {
                _internal_energy_parameters = std::vector<double>(value_start, value_end);
            }

            /*!
             * Set the parameters for the heat conduction
             *
             * \param &value_start: The starting value of the parameter iterator
             * \param &value_end: The stopping value of the parameter iterator
             */
            template <class parameter_iter>
            void setHeatConductionParameters(const parameter_iter &value_start, const parameter_iter &value_end) {
                _heat_conduction_parameters = std::vector<double>(value_start, value_end);
            }

            /*!
             * Set the index of the temperature gradient
             *
             * \param &value: The value of the temperature gradient index
             */
            void setTemperatureGradientIndex(const unsigned int &value) {
                _temperature_gradient_index     = value;
                _temperature_gradient_index_set = true;
            }

            /*!
             * Get the temperature gradient index
             */
            const unsigned int getTemperatureGradientIndex() {
                TARDIGRADE_ERROR_TOOLS_CHECK(_temperature_gradient_index_set,
                                             "The temperature gradient index must be set before it is called");

                return _temperature_gradient_index;
            }

            /*!
             * Define all of the residual classes
             */
            virtual void setResidualClasses() override {
                stress = tardigradeHydra::linearElasticity::residual(this, getStressSize(), *getStressParameters());

                internal_energy = tardigradeHydra::linearInternalEnergy::residual(this, getInternalEnergySize(),
                                                                                  *getInternalEnergyParameters(), 1, 9);

                heat_conduction =
                    tardigradeHydra::fourierHeatConduction::residual(this, getHeatConductionSize(),
                                                                     *getHeatConductionParameters(),
                                                                     getTemperatureGradientIndex(), 1, 10);

                std::vector<tardigradeHydra::ResidualBase<> *> residuals(3);

                residuals[0] = &stress;

                residuals[1] = &internal_energy;

                residuals[2] = &heat_conduction;

                setResidualClasses(residuals);
            }

            /*! Get the stress parameters */
            const std::vector<double> *getStressParameters() { return &_stress_parameters; }

            /*! Get the internal energy parameters */
            const std::vector<double> *getInternalEnergyParameters() { return &_internal_energy_parameters; }

            /*! Get the heat conduction parameters */
            const std::vector<double> *getHeatConductionParameters() { return &_heat_conduction_parameters; }

            /*! Get the number of nonlinear unknowns in the stress response */
            const unsigned int getStressSize() { return _stress_size; }

            /*! Get the number of nonlinear unknowns in the internal energy response */
            const unsigned int getInternalEnergySize() { return _internal_energy_size; }

            /*! Get the number of nonlinear unknowns in the heat conduction response */
            const unsigned int getHeatConductionSize() { return _heat_conduction_size; }

           protected:
            tardigradeHydra::linearElasticity::residual stress;  //!< The residual that defines the stress
            tardigradeHydra::linearInternalEnergy::residual
                internal_energy;  //!< The residual that defines the internal energy
            tardigradeHydra::fourierHeatConduction::residual
                heat_conduction;  //!< The residual that defines the heat conduction

           private:
            //! The number of nonlinear unknowns in the stress response (must be at least 9)
            unsigned int _stress_size          = 9;
            //! The number of nonlinear unknowns in the internal energy response (must be at least 1)
            unsigned int _internal_energy_size = 1;
            //! The number of nonlinear unknowns in the heat conduction response (must be at least 3)
            unsigned int _heat_conduction_size = 3;
            //! The parameters associated with the calculation of the stress
            std::vector<double> _stress_parameters;
            //! The parameters associated with the calculation of the internal energy
            std::vector<double> _internal_energy_parameters;
            //! The parameters associated with the calculation of the heat conduction
            std::vector<double> _heat_conduction_parameters;
            //! The index of the temperature gradient
            unsigned int _temperature_gradient_index = 0;
            //! Flag for if the temperature gradient is set
            bool _temperature_gradient_index_set     = false;
        };

        /*!
         * A basic continuum solid
         */
        class BasicSolid : public tardigradeCMML::CMMLMaterial {
           public:
            /*!
             * A basic continuum solid.
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
            BasicSolid() {
                setName("BasicSolid");
                setEvaluateModelResultSize(23);
            }

            /*!
             * Set the index the displacement gradient is located in the dof vector
             *
             * \param value: The displacement gradient index
             */
            void setDisplacementGradientIndex(unsigned int value) { _displacement_gradient_index = value; }

            /*!
             * Set the index the temperature is located in the dof vector
             *
             * \param value: The temperature index
             */
            void setTemperatureIndex(unsigned int value) { _temperature_index = value; }

            virtual int evaluate_model(const time_type &current_time, const time_type &dt,
                                       const current_dof_type  *current_dof_begin,
                                       const previous_dof_type *previous_dof_begin, const unsigned int dof_size,
                                       const parameter_type *parameters_begin, const unsigned int parameters_size,
                                       sdvs_type *sdvs_begin, const unsigned int sdvs_size, result_type *result_begin,
                                       const unsigned int result_size, std::string &output_message) override;

            virtual int evaluate_model(const time_type &current_time, const time_type &dt,
                                       const current_dof_type  *current_dof_begin,
                                       const previous_dof_type *previous_dof_begin, const unsigned int dof_size,
                                       const parameter_type *parameters_begin, const unsigned int parameters_size,
                                       sdvs_type *sdvs_begin, const unsigned int sdvs_size, result_type *result_begin,
                                       const unsigned int result_size, jacobian_type *jacobian_begin,
                                       additional_type *additional_begin, const unsigned int additional_size,
                                       std::string &output_message) override;

            //! We assume that the material is being evaluated in 3D
            constexpr static unsigned int dim = 3;
            /*! Get if the dof gradients are with respect to the current configuration */
            const bool getIsCurrent() { return _isCurrent; }
            /*! Get the displacement gradient index */
            const unsigned int getDisplacementGradientIndex() { return _displacement_gradient_index; }
            /*! Get the temperature index */
            const unsigned int getTemperatureIndex() { return _temperature_index; }
            /*! Get the stress index in the material response vector */
            const unsigned int getStressIndex() { return _stress_index; }
            /*! Get the internal energy index in the material response vector */
            const unsigned int getInternalEnergyIndex() { return _internal_energy_index; }
            /*! Get the heat flux index in the material response vector */
            const unsigned int getHeatFluxIndex() { return _heat_flux_index; }

           protected:
            virtual void formDeformationGradients(const current_dof_type         *current_dof_begin,
                                                  const previous_dof_type        *previous_dof_begin,
                                                  std::vector<current_dof_type>  &deformationGradient,
                                                  std::vector<previous_dof_type> &previousDeformationGradient);

            virtual void formDeformationGradients(const current_dof_type         *current_dof_begin,
                                                  const previous_dof_type        *previous_dof_begin,
                                                  std::vector<current_dof_type>  &deformationGradient,
                                                  std::vector<previous_dof_type> &previousDeformationGradient,
                                                  std::vector<current_dof_type>  &dFdGradU);

            /*!
             * Extract the material parameters
             *
             * The parameter vector is assumed to be organized as
             *
             * displacement gradient index, temperature index, temperature gradient index, lambda, mu, specifc heat,
             * conductivity
             *
             * where lambda and mu are the linear elastic Lame parameters, specific heat is the linear
             * relationship between temperature and the internal energy, and the conductivity is the
             * conductivity parameter for the heat flux.
             *
             * \param *parameters_begin: A pointer to the initial element of the parameter vector
             * \param parameters_size: The size of the parameters vector
             */
            virtual void extract_parameters(const parameter_type *parameters_begin,
                                            const unsigned int    parameters_size) {
                TARDIGRADE_ERROR_TOOLS_CHECK(parameters_size == 7, "The parameters vector must have seven values")

                setDisplacementGradientIndex((unsigned int)(*(parameters_begin + 0) + 0.5));
                setTemperatureIndex((unsigned int)(*(parameters_begin + 1) + 0.5));
            }

           private:
            //! Whether the dof gradients are with respect to the current configuration or not
            bool _isCurrent = true;
            //! The index of the displacement gradient in the dof vector
            unsigned int _displacement_gradient_index;
            //! The index of the temperature in the dof vector
            unsigned int _temperature_index;

            //! The index of the stress in the material response vector
            unsigned int _stress_index          = 0;
            //! The index of the internal energy in the material response vector
            unsigned int _internal_energy_index = 9;
            //! The index of the heat flux in the material response vector
            unsigned int _heat_flux_index       = 17;
        };

        //! Register the material in the library
        REGISTER_MATERIAL(tardigradeCMML::BasicSolid::BasicSolid);

    }  // namespace BasicSolid

}  // namespace tardigradeCMML

#include "tardigrade_BasicSolid.cpp"

#endif
