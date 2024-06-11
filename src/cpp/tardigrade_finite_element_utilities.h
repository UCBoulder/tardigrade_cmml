/**
  ******************************************************************************
  * \file tardigrade_finite_element_utilities.h
  ******************************************************************************
  * The header file for utilities which can assist using the balance equations
  * in finite element codes. We here assume that unknown quantities are
  * computed at the evaluation point using interpolation functions (here called
  * interp) and are projected to the nodes using test functions (here called
  * test).
  ******************************************************************************
  */

#ifndef TARDIGRADE_FINITE_ELEMENT_UTILITIES_H
#define TARDIGRADE_FINITE_ELEMENT_UTILITIES_H

#include<array>

#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace finiteElementUtilities{

        constexpr unsigned int dim = 3; //!< Set the dimension as 3D by default

        constexpr unsigned int sot_dim = dim * dim; //!< Set the dimensions of a standard second order tensor

        typedef double floatType; //!< Define the float type as a double

        typedef std::array< floatType, dim > floatVector; //!< Define a standard vector

        typedef std::array< floatType, sot_dim > secondOrderTensor; //!< Define a standard second-order tensor

        template<typename grad_iterator, typename output_iterator>
        void computeGradientSpatialJacobian( const grad_iterator &grad_a_start, const unsigned int grad_a_size, 
                                             floatVector grad_test, const unsigned int index, output_iterator dgrad_adui_start );

    }

}

#include "tardigrade_finite_element_utilities.cpp"

#endif
