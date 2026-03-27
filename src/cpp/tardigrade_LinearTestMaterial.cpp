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

    namespace LinearTestMaterial{

        /*!
         * Get the full tangent matrix
         */
        std::vector<double> hydraLinearTest::getFullTangent() {
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

        /*!
         * Get the temperature from the additional DOF vector
         *
         * \param _nphases: The number of phases
         * \param _active_phase: The current active phase
         * \param additional_dof: The additional DOF vector
         */
        double hydraLinearTest::_getAdditionalDOFTemperature(const unsigned int _nphases, const unsigned int _active_phase,
                                                       const std::vector<double> &additional_dof) {
            unsigned int offset = _nphases * (1 + 3 + 3) + _active_phase;

            return additional_dof[offset];
        }

        /*!
         * Get the deformation gradient from the additional DOF vector
         *
         * \param _nphases: The number of phases
         * \param _active_phase: The current active phase
         * \param &_num_phase_dof: The number of degrees of freedom in a given phase
         * \param &_num_add_dof: The number of additional degrees of freedom
         * \param additional_dof: The additional DOF vector
         */
        std::vector<double> hydraLinearTest::_getAdditionalDOFDeformationGradient(
            const unsigned int _nphases, const unsigned int _active_phase, const unsigned int _num_phase_dof,
            const unsigned int _num_add_dof, const std::vector<double> &additional_dof) {
            unsigned int offset = _nphases * _num_phase_dof + _num_add_dof + _nphases * 3 + 9 * _active_phase;

            std::vector<double> gradW(std::begin(additional_dof) + offset, std::begin(additional_dof) + offset + 9);

            std::vector<double> F;

            tardigradeConstitutiveTools::computeDeformationGradient(gradW, F, true);

            return F;
        }
    }

}
