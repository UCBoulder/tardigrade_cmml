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
         * A linear test material
         */
        class LinearTestMaterial : public tardigradeCMML::CMMLMaterial {
           public:

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

           protected:

           private:

        };

        //! Register the material in the library
        REGISTER_MATERIAL(tardigradeCMML::LinearTestMaterial::LinearTestMaterial);

    }

}

#include "tardigrade_LinearTestMaterial.cpp"

#endif
