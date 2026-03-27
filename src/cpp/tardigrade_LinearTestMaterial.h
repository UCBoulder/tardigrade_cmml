/*!
=====================================================================
|                  tardigrade_LinearTestMaterial.h                  |
=====================================================================
| A material which defines a response as being a linear product     |
| between the parameters and the degrees of freedom. To be used for |
| testing purposes.                                                 |
=====================================================================
*/

#ifndef TARDIGRADE_LINEARTESTMATERIAL
#define TARDIGRADE_LINEARTESTMATERIAL

#include "tardigrade_error_tools.h"
#include "tardigrade_constitutive_tools.h"
#include "tardigrade_hydra.h"
#include "tardigrade_hydraLinearTestMaterial.h"
#define TARDIGRADE_CMML_SKIP_HEADER_ONLY_MATERIALS_INCLUDE
#include "tardigrade_cmml.h"
#undef TARDIGRADE_CMML_SKIP_HEADER_ONLY_MATERIALS_INCLUDE

namespace tardigradeCMML {

    namespace LinearTestMaterial{

        /*!
         * A class which linearly relates the degrees of freedom to the response
         */
        class hydraLinearTest : public tardigradeHydra::hydraBase {
          public:
           /*!
            * A model that connects the degrees of freedom to the output response linearly
            *
            * \param &_nphases: The number of phases
            * \param &_active_phase: The currently active phase
            * \param &_num_phase_dof: The number of degrees of freedom in a given phase
            * \param &_num_add_dof: The number of additional degrees of freedom
            * \param &t: The current time
            * \param &dt: The change in time
            * \param &additionalDOF: The additional degrees of freedom
            * \param &previousAdditionalDOF: The previous additional degrees of freedom
            * \param &parameters: The parameter vector
            * \param &num_sdvs: The number of state variables
            * \param &num_nonlinear_sdvs: The number of nonlinear sdvs
            */
           hydraLinearTest(const unsigned int &_nphases, const unsigned int &_active_phase, const unsigned int &_num_phase_dof,
                           const unsigned int &_num_add_dof, const double &t, const double &dt,
                           const std::vector<double> &additionalDOF, const std::vector<double> &previousAdditionalDOF,
                           const std::vector<double> &parameters, const unsigned int num_sdvs, const unsigned int num_nonlinear_sdvs)
               : tardigradeHydra::hydraBase(dof_storage_class, model_configuration_class),
                 dof_storage_class(tardigradeHydra::DOFStorageBase(
                     t, dt, _getAdditionalDOFTemperature(_nphases, _active_phase, additionalDOF),
                     _getAdditionalDOFTemperature(_nphases, _active_phase, additionalDOF),
                     _getAdditionalDOFDeformationGradient(_nphases, _active_phase, _num_phase_dof, _num_add_dof,
                                                          additionalDOF),
                     _getAdditionalDOFDeformationGradient(_nphases, _active_phase, _num_phase_dof, _num_add_dof,
                                                          previousAdditionalDOF),
                     additionalDOF, previousAdditionalDOF)),
                 model_configuration_class(tardigradeHydra::ModelConfigurationBase(std::vector<double>(num_sdvs, 0),
                                                                                   parameters, 1, num_nonlinear_sdvs)) {
               nphases = _nphases;

               active_phase = _active_phase;

               num_phase_dof = _num_phase_dof;

               num_add_dof = _num_add_dof;
           }

           //! The degree of freedom storage class
           tardigradeHydra::DOFStorageBase dof_storage_class;

           //! The model configuration class
           tardigradeHydra::ModelConfigurationBase model_configuration_class;

           //! The residual class
           tardigradeHydra::linearTestMaterial::residual residual;  //!< The residual class

           virtual std::vector<double> getFullTangent();

          protected:
           //! The number of phases
           unsigned int nphases;

           //! The currently active phase
           unsigned int active_phase;

           //! The number of degrees of freedom in a phase
           unsigned int num_phase_dof;

           //! The number of additional degrees of freedom
           unsigned int num_add_dof;

           virtual double _getAdditionalDOFTemperature(const unsigned int _nphases, const unsigned int _active_phase,
                                                          const std::vector<double> &additional_dof) final;

           /*!
            * Get the deformation gradient from the additional DOF vector
            *
            * \param _nphases: The number of phases
            * \param _active_phase: The current active phase
            * \param &_num_phase_dof: The number of degrees of freedom in a given phase
            * \param &_num_add_dof: The number of additional degrees of freedom
            * \param additional_dof: The additional DOF vector
            */
           virtual std::vector<double> _getAdditionalDOFDeformationGradient(
               const unsigned int _nphases, const unsigned int _active_phase, const unsigned int _num_phase_dof,
               const unsigned int _num_add_dof, const std::vector<double> &additional_dof) final;

          private:
           using tardigradeHydra::hydraBase::setResidualClasses;

           /*!
            * Set the vector of residual classes (in this case, a single residual)
            */
           virtual void setResidualClasses() override {

               std::vector<tardigradeHydra::ResidualBase<> *> residuals(1);

               TARDIGRADE_ERROR_TOOLS_CATCH(residual =
                                                tardigradeHydra::linearTestMaterial::residual(this, 9 + model_configuration->_num_nonlinear_solve_state_variables, model_configuration->_parameters));

               residuals[0] = &residual;

               setResidualClasses(residuals);
           }
        };

    }

}

#include "tardigrade_LinearTestMaterial.cpp"

#endif
