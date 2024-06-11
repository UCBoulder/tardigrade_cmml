/**
  ******************************************************************************
  * \file tardigrade_finite_element_utilities.cpp
  ******************************************************************************
  * The source file for utilities which can assist using the balance equations
  * in finite element codes. We here assume that unknown quantities are
  * computed at the evaluation point using interpolation functions (here called
  * interp) and are projected to the nodes using test functions (here called
  * test).
  ******************************************************************************
  */

namespace tardigradeBalanceEquations{

    namespace finiteElementUtilities{

        template<typename grad_iterator, typename output_iterator>
        void computeGradientSpatialJacobian( const grad_iterator &grad_a_start, const unsigned int grad_a_size,
                                             floatVector grad_test, const unsigned int index, output_iterator dgrad_adui_start ){
            /*!
             * Compute the derivative of the spatial gradient of a quantity a
             * in the current configuration w.r.t. the spatial degrees of freedom
             * 
             * \f$ \frac{D}{Du_a} \left( a_{i,j} \right) \f$
             * 
             * \param &grad_a_start: An iterator representing the start of the gradient of the quantity
             * \param &grad_a_size: The size of the gradient of a
             * \param &grad_test: The gradient of the test function
             * \param &index: The index of the spatial degree of freedom to compute the derivative for (0, 1 or 2)
             * \param &dgrad_adui_start: An iterator representing the start of the output
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( ( grad_a_size % dim ) == 0, "The incoming spatial gradient has a dimension of " + std::to_string( grad_a_size ) + " which is not a multiple of " + std::to_string( dim ) );

            const unsigned int a_size = grad_a_size / dim;

            std::fill( dgrad_adui_start, dgrad_adui_start + a_size, 0 );

            for ( unsigned int i = 0; i < a_size; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    *( dgrad_adui_start + dim * i + j ) -= *( grad_a_start + dim * i + index ) * grad_test[ j ];

                }

            }

        }

    }

}
