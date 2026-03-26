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

        class hydraLinearTest : public tardigradeHydra::hydraBase {
          public:
           hydraLinearTest(const unsigned int &_nphases, const unsigned int &_active_phase, const unsigned int &_num_phase_dof,
                           const unsigned int &_num_add_dof, const double &t, const double &dt,
                           const std::vector<double> &additionalDOF, const std::vector<double> &previousAdditionalDOF,
                           const std::vector<double> &parameters)
               : tardigradeHydra::hydraBase(dof_storage_class, model_configuration_class),
                 dof_storage_class(tardigradeHydra::DOFStorageBase(
                     t, dt, _getAdditionalDOFTemperature(_nphases, _active_phase, additionalDOF),
                     _getAdditionalDOFTemperature(_nphases, _active_phase, additionalDOF),
                     _getAdditionalDOFDeformationGradient(_nphases, _active_phase, _num_phase_dof, _num_add_dof,
                                                          additionalDOF),
                     _getAdditionalDOFDeformationGradient(_nphases, _active_phase, _num_phase_dof, _num_add_dof,
                                                          previousAdditionalDOF),
                     additionalDOF, previousAdditionalDOF)),
                 model_configuration_class(tardigradeHydra::ModelConfigurationBase(std::vector<double>(13, 0),
                                                                                   parameters, 1, 13)) {
               nphases = _nphases;

               active_phase = _active_phase;

               num_phase_dof = _num_phase_dof;

               num_add_dof = _num_add_dof;
           }

           tardigradeHydra::DOFStorageBase dof_storage_class;

           tardigradeHydra::ModelConfigurationBase model_configuration_class;

           //! The residual class
           tardigradeHydra::linearTestMaterial::residual residual;  //!< The residual class

           virtual std::vector<double> getFullTangent() {
               constexpr unsigned int dim = 3;

               computeTangents();

               computedXdAdditionalDOF();

               std::vector<double> dFdGradU(81, 0);

               for (unsigned int i = 0; i < 3; ++i) {
                   for (unsigned int I = 0; I < 3; ++I) {
                       for (unsigned int a = 0; a < 3; ++a) {
                           for (unsigned int b = 0; b < 3; ++b) {
                               dFdGradU[dim * dim * dim * i + dim * dim * I + dim * a + b] +=
                                   (*getDeformationGradient())[dim * i + a] * (*getDeformationGradient())[dim * b + I];
                           }
                       }
                   }
               }

               std::vector<double> full_jacobian(getFlatdXdAdditionalDOF()->size(), 0);

               std::copy(std::begin(*getFlatdXdAdditionalDOF()), std::end(*getFlatdXdAdditionalDOF()),
                         std::begin(full_jacobian));

               // Incorporate the Jacobian of the temperature and deformation gradient to the full Jacobian
               unsigned int offset = nphases * num_phase_dof + num_add_dof + nphases * 3 + 9 * active_phase;

               for (unsigned int I = 0; I < getNumUnknowns(); ++I) {
                   full_jacobian[getAdditionalDOF()->size() * I + nphases * (1 + 3 + 3) + active_phase] += (*getFlatdXdT())[I];

                   for (unsigned int ij = 0; ij < dim * dim; ++ij) {
                       for (unsigned int ab = 0; ab < dim * dim; ++ab) {
                           full_jacobian[getAdditionalDOF()->size() * I + offset + ab] +=
                               (*getFlatdXdF())[dim * dim * I + ij] * dFdGradU[dim * dim * ij + ab];
                       }
                   }
               }

               return full_jacobian;
           }

          protected:
           unsigned int nphases;

           unsigned int active_phase;

           unsigned int num_phase_dof;

           unsigned int num_add_dof;

           virtual double _getAdditionalDOFTemperature(const unsigned int _nphases, const unsigned int _active_phase,
                                                          const std::vector<double> &additional_dof) final {
               unsigned int offset = _nphases * (1 + 3 + 3) + _active_phase;

               return additional_dof[offset];
           }

           virtual std::vector<double> _getAdditionalDOFDeformationGradient(
               const unsigned int _nphases, const unsigned int _active_phase, const unsigned int _num_phase_dof,
               const unsigned int _num_add_dof, const std::vector<double> &additional_dof) final {
               unsigned int offset = _nphases * _num_phase_dof + _num_add_dof + _nphases * 3 + 9 * _active_phase;

               std::vector<double> gradW(std::begin(additional_dof) + offset, std::begin(additional_dof) + offset + 9);

               std::vector<double> F;

               tardigradeConstitutiveTools::computeDeformationGradient(gradW, F, true);

               return F;
           }

          private:
           using tardigradeHydra::hydraBase::setResidualClasses;

           virtual void setResidualClasses() override {
               /*!
                * Set the vector of residual classes (in this case, a single residual)
                */

               std::vector<tardigradeHydra::ResidualBase<> *> residuals(1);

               TARDIGRADE_ERROR_TOOLS_CATCH(residual =
                                                tardigradeHydra::linearTestMaterial::residual(this, 22, model_configuration->_parameters));

               residuals[0] = &residual;

               setResidualClasses(residuals);
           }
        };

    }

}

#endif
