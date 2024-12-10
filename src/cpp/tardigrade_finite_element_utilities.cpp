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

#include<Eigen/Dense>

namespace tardigradeBalanceEquations{

    namespace finiteElementUtilities{

        template<typename grad_iterator, typename output_iterator>
        void computeGradientSpatialJacobian( const grad_iterator &grad_a_start, const unsigned int grad_a_size,
                                             floatVector grad_interp, const unsigned int index, output_iterator dgrad_adui_start ){
            /*!
             * Compute the derivative of the spatial gradient of a quantity a
             * in the current configuration w.r.t. the spatial degrees of freedom
             * 
             * \f$ \frac{D}{Du_a} \left( a_{i,j} \right) \f$
             * 
             * \param &grad_a_start: An iterator representing the start of the gradient of the quantity
             * \param &grad_a_size: The size of the gradient of a
             * \param &grad_interp: The gradient of the interpolation function
             * \param &index: The index of the spatial degree of freedom to compute the derivative for (0, 1 or 2)
             * \param &dgrad_adui_start: An iterator representing the start of the output
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( ( grad_a_size % dim ) == 0, "The incoming spatial gradient has a dimension of " + std::to_string( grad_a_size ) + " which is not a multiple of " + std::to_string( dim ) );

            const unsigned int a_size = grad_a_size / dim;

            std::fill( dgrad_adui_start, dgrad_adui_start + a_size, 0 );

            for ( unsigned int i = 0; i < a_size; ++i ){

                for ( unsigned int j = 0; j < dim; ++j ){

                    *( dgrad_adui_start + dim * i + j ) -= *( grad_a_start + dim * i + index ) * grad_interp[ j ];

                }

            }

        }

        //! A base class for a simple finite element formulation useful for testing
        template<int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in, class shape_functions_out, class grad_shape_functions_out>
        FiniteElementBase<dim,local_dim,node_count,node_in,local_node_in,local_point_in,shape_functions_out,grad_shape_functions_out>::FiniteElementBase(
            const node_in &_x_begin,             const node_in &_x_end,
            const node_in &_X_begin,             const node_in &_X_end,
            const local_node_in &_local_node_xi_begin, const local_node_in &_local_node_xi_end )
          : x_begin( _x_begin ), x_end( _x_end ), X_begin( _X_begin ), X_end( _X_end ), local_node_xi_begin( _local_node_xi_begin ), local_node_xi_end( _local_node_xi_end ){
            /*!
             * Constructor for the base finite element class
             * 
             * \param &_x_begin: The starting iterator for the current node positions
             * \param &_x_end: The stopping iterator for the current node positions
             * \param &_X_begin: The starting iterator for the reference node positions
             * \param &_X_end: The stopping iterator for the reference node positions
             * \param &_local_node_xi_begin: The starting iterator for the local node positions
             * \param &_local_node_xi_end: The stopping iterator for the local node positions
             */

        }

        template<typename T, class node_in, class local_point_in, class shape_functions_out, class grad_shape_functions_out>
        LinearHex<T,node_in,local_point_in,shape_functions_out,grad_shape_functions_out>::LinearHex( const node_in &_x_begin, const node_in &_x_end, const node_in &_X_begin, const node_in &_X_end ) :
                              FiniteElementBase<3,3,8,node_in,typename std::array< T, 3 * 8 >::const_iterator,local_point_in,shape_functions_out,grad_shape_functions_out>( _x_begin, _x_end, _X_begin, _X_end, std::cbegin( local_nodes ), std::cend( local_nodes ) ){
            /*!
             * Constructor for the base finite element class
             * 
             * \param &_x_begin: The starting iterator for the current node positions
             * \param &_x_end: The stopping iterator for the current node positions
             * \param &_X_begin: The starting iterator for the reference node positions
             * \param &_X_end: The stopping iterator for the reference node positions
             */

        }

        template<typename T, class node_in, class local_point_in, class shape_functions_out, class grad_shape_functions_out>
        void LinearHex<T,node_in,local_point_in,shape_functions_out,grad_shape_functions_out>::GetShapeFunctions( const local_point_in &xi_begin, const local_point_in &xi_end, shape_functions_out N_begin, shape_functions_out N_end ){
            /*!
             * Get the shape functions of the linear hexahedral element
             *
             * \param &xi_begin: The starting iterator of the local coordinates (must have dimension 3)
             * \param &xi_end: The stopping iterator of the local coordinates (must have dimension 3)
             * \param N_begin: The starting iterator of the shape functions (must have dimension 8)
             * \param N_end: The stopping iterator of the shape functions (must have dimension 8)
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( ( size_type )( N_end - N_begin ) == 8, "The dimension of the shape-function iterator is " + std::to_string( ( size_type )( N_end - N_begin ) ) );

            for ( std::pair< size_type, shape_functions_out > i( 0, N_begin ); i.second != N_end; ++i.first, ++i.second ){

                *i.second = ( 1 + ( *( this->local_node_xi_begin + 3 * i.first + 0 ) ) * ( *( xi_begin + 0 ) ) ) * ( 1 + ( *( this->local_node_xi_begin + 3 * i.first + 1 ) ) * ( *( xi_begin + 1 ) ) ) * ( 1 + ( *( this->local_node_xi_begin + 3 * i.first + 2 ) ) * ( *( xi_begin + 2 ) ) ) / 8;

            }

        }

        template<typename T, class node_in, class local_point_in, class shape_functions_out, class grad_shape_functions_out>
        void LinearHex<T,node_in,local_point_in,shape_functions_out,grad_shape_functions_out>::GetLocalShapeFunctionGradients( const local_point_in &xi_begin, const local_point_in &xi_end, grad_shape_functions_out dNdxi_begin, grad_shape_functions_out dNdxi_end ){
            /*!
             * Get the shape functions of the linear hexahedral element
             *
             * \param &xi_begin: The starting iterator of the local coordinates (must have dimension 3)
             * \param &xi_end: The stopping iterator of the local coordinates (must have dimension 3)
             * \param dNdxi_begin: The starting iterator of the local gradient of the shape functions (must have dimension 24)
             * \param dNdxi_end: The stopping iterator of the local gradient of the shape functions (must have dimension 24)
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( ( size_type )( dNdxi_end - dNdxi_begin ) == 24, "The dimension of the shape-function iterator is " + std::to_string( ( size_type )( dNdxi_end - dNdxi_begin ) ) );

            for ( size_type i = 0; i < 8; ++i ){

                *( dNdxi_begin + 3 * i + 0 ) = ( *( this->local_node_xi_begin + 3 * i + 0 ) ) * ( 1 + ( *( this->local_node_xi_begin + 3 * i + 1 ) ) * ( *( xi_begin + 1 ) ) ) * ( 1 + ( *( this->local_node_xi_begin + 3 * i + 2 ) ) * ( *( xi_begin + 2 ) ) ) / 8;
                *( dNdxi_begin + 3 * i + 1 ) = ( 1 + ( *( this->local_node_xi_begin + 3 * i + 0 ) ) * ( *( xi_begin + 0 ) ) ) * ( *( this->local_node_xi_begin + 3 * i + 1 ) ) * ( 1 + ( *( this->local_node_xi_begin + 3 * i + 2 ) ) * ( *( xi_begin + 2 ) ) ) / 8;
                *( dNdxi_begin + 3 * i + 2 ) = ( 1 + ( *( this->local_node_xi_begin + 3 * i + 0 ) ) * ( *( xi_begin + 0 ) ) ) * ( 1 + ( *( this->local_node_xi_begin + 3 * i + 1 ) ) * ( *( xi_begin + 1 ) ) ) * ( *( this->local_node_xi_begin + 3 * i + 2 ) ) / 8;

            }

        }

        template<typename T, class node_in, class local_point_in, class shape_functions_out, class grad_shape_functions_out>
        void LinearHex<T,node_in,local_point_in,shape_functions_out,grad_shape_functions_out>::GetGlobalShapeFunctionGradients(
            const local_point_in &xi_begin, const local_point_in &xi_end,
            const node_in &node_positions_begin, const node_in &node_positions_end,
            grad_shape_functions_out value_begin, grad_shape_functions_out value_end ){
            /*!
             * Compute the global gradient of the quantity in the global coordinates
             *
             * \param &xi_begin: The starting iterator of the local point
             * \param &xi_end: The stopping iterator of the local point
             * \param &node_positions_begin: The starting iterator of the nodal positions (row major)
             * \param &node_positions_end: The stopping iterator of the nodal positions (row major)
             * \param &value_begin: The starting iterator of the shape function global gradient (row major)
             * \param &value_end: The stopping iterator of the shape function global gradient (row major)
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( ( size_type )( value_end - value_begin ) == 24, "The shape function global gradient must have a size of 24" );

            std::array< typename std::iterator_traits<node_in>::value_type, 9 > dxdxi;
            std::array< typename std::iterator_traits<node_in>::value_type, 9 > dxidx;

            TARDIGRADE_ERROR_TOOLS_CATCH( this->GetLocalQuantityGradient( xi_begin, xi_end, node_positions_begin, node_positions_end,
                                                                          std::begin( dxdxi ), std::end( dxdxi ) ) );

            Eigen::Map< const Eigen::Matrix< typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor > > _dxdxi( dxdxi.data( ) );

            Eigen::Map< Eigen::Matrix< typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor > > _dxidx( dxidx.data( ) );

            _dxidx = ( _dxdxi.inverse( ) ).eval( );

            std::fill( value_begin, value_end, 0 );

            for ( unsigned int node = 0; node < 8; ++node ){

                for ( unsigned int inner = 0; inner < 3; ++inner ){

                    for ( unsigned int outer = 0; outer < 3; ++outer ){

                        *( value_begin + 3 * node + outer ) += this->_local_gradshapefunctions[ 3 * node + inner ] * dxidx[ 3 * inner + outer ];

                    }

                }

            }

        }

        template<typename T, class node_in, class local_point_in, class shape_functions_out, class grad_shape_functions_out>
        void LinearHex<T,node_in,local_point_in,shape_functions_out,grad_shape_functions_out>::GetVolumeIntegralJacobianOfTransformation(
            const local_point_in &xi_begin, const local_point_in &xi_end,
            typename std::iterator_traits<node_in>::value_type &value, const bool configuration ){
            /*!
             * Compute the value of the Jacobian of transformation from the local coordinates to the configuration for volume integrals
             * 
             * \param &xi_begin: The starting iterator of the local coordinates
             * \param &xi_end: The stopping iterator of the local coordinates
             * \param value: The Jacobian of transformation going from the local coordinates to the indicated configuration
             * \param configuration: Compute the gradient w.r.t. the current configuration ( true ) or reference configuration ( false )
             */

            std::array< typename std::iterator_traits<node_in>::value_type, 9 > dxdxi;

            if ( configuration ){

                this->GetLocalQuantityGradient( xi_begin, xi_end, this->x_begin, this->x_end, std::begin( dxdxi ), std::end( dxdxi ) );

            }
            else{

                this->GetLocalQuantityGradient( xi_begin, xi_end, this->X_begin, this->X_end, std::begin( dxdxi ), std::end( dxdxi ) );

            }

            Eigen::Map< const Eigen::Matrix< typename std::iterator_traits<node_in>::value_type, 3, 3, Eigen::RowMajor > > _dxdxi( dxdxi.data( ) );

            value = _dxdxi.determinant( );

        }

    }

}
