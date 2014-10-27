/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_curve_piecewise_cst_airfoil_creator_hpp
#define eli_geom_curve_piecewise_cst_airfoil_creator_hpp

#include <vector>

#include "eli/code_eli.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise_polynomial_creator.hpp"
#include "eli/geom/curve/pseudo/cst_airfoil.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_cst_airfoil_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef eli::geom::curve::pseudo::cst_airfoil<data_type> cst_airfoil_type;
          typedef typename cst_airfoil_type::point_type cst_airfoil_point_type;
          typedef typename cst_airfoil_type::control_point_type cst_airfoil_control_point_type;
          typedef unsigned short dimension_type;

          piecewise_cst_airfoil_creator() : piecewise_creator_base<data_type, dim__, tolerance_type>(2, 0), cst(1) {}
          piecewise_cst_airfoil_creator(const piecewise_cst_airfoil_creator<data_type, dim__, tolerance_type> &pca)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(pca), cst(pca.cst) {}
          ~piecewise_cst_airfoil_creator() {}

          void set_airfoil(const cst_airfoil_type &ca)
          {
            cst = ca;
          }

          void get_airfoil(cst_airfoil_type &ca) const
          {
            ca = cst;
          }

          bool set_conditions(const cst_airfoil_type &ca)
          {
            set_airfoil(ca);

            return true;
          }


          /** This method creates an exact bezier representation of a CST airfoil upper and 
           *  lower surface. It uses the method from "Creating Exact Bezier Representations of CST
           *  Shapes" by Marshall, AIAA paper 2013-3077.
           */
          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;

            typename curve_type::monomial_coefficient_type bezu_mono_coef, bezl_mono_coef;
            typename cst_airfoil_type::monomial_coefficient_type cstu_mono_coef, cstl_mono_coef;

            // extract the monomial coefficients from the CST airfoil
            cst.get_upper_monomial_coefficients(cstu_mono_coef);
            cst.get_lower_monomial_coefficients(cstl_mono_coef);

            // convert the original coefficients to new monomial coefficients
            index_type i, nu(cstu_mono_coef.rows()-1), nl(cstl_mono_coef.rows()-1);
            data_type dte(cst.get_trailing_edge_thickness()/2);
            bezu_mono_coef.resize(2*nu+3+1, dim__);
            bezl_mono_coef.resize(2*nl+3+1, dim__);
            bezu_mono_coef.setZero();
            bezl_mono_coef.setZero();

            // upper
            i=0;
            bezu_mono_coef(2*i+1, 1)=cstu_mono_coef(i, 0);
            bezu_mono_coef(2*i+2, 0)=1;
            bezu_mono_coef(2*i+2, 1)=dte;
            for (i=1; i<=nu; ++i)
            {
              bezu_mono_coef(2*i+1, 1)=cstu_mono_coef(i, 0)-cstu_mono_coef(i-1, 0);
            }
            i=nu+1;
            bezu_mono_coef(2*i+1, 1)=-cstu_mono_coef(i-1);

            // lower
            i=0;
            bezl_mono_coef(2*i+1, 1)=cstl_mono_coef(i, 0);
            bezl_mono_coef(2*i+2, 0)=1;
            bezl_mono_coef(2*i+2, 1)=-dte;
            for (i=1; i<=nu; ++i)
            {
              bezl_mono_coef(2*i+1, 1)=cstl_mono_coef(i, 0)-cstl_mono_coef(i-1, 0);
            }
            i=nu+1;
            bezl_mono_coef(2*i+1, 1)=-cstl_mono_coef(i-1);

            // create the lower and upper curve
            piecewise_polynomial_creator<data_type, dim__, tolerance_type> poly_creator;
            pseudo::polynomial<data__, dim__> cu, cl;
            bool rtn_flag;

            for (dimension_type j=0; j<dim__; ++j)
            {
              cu.set_coefficients(bezu_mono_coef.col(j), j);
              cl.set_coefficients(bezl_mono_coef.col(j), j);
            }

            // lower
            rtn_flag=poly_creator.set_conditions(cl);
            if (!rtn_flag)
            {
              assert(false);
              return false;
            }
            poly_creator.set_t0(this->get_t0());
            poly_creator.set_segment_dt(this->get_segment_dt(0), 0);
            rtn_flag=poly_creator.create(pc);
            if (!rtn_flag)
            {
              assert(false);
              return false;
            }
            pc.reverse();

            // upper
            piecewise_curve_type pc_temp;
            rtn_flag=poly_creator.set_conditions(cu);
            if (!rtn_flag)
            {
              assert(false);
              return false;
            }
            poly_creator.set_t0(this->get_t0()+this->get_segment_dt(0));
            poly_creator.set_segment_dt(this->get_segment_dt(1), 0);
            rtn_flag=poly_creator.create(pc_temp);
            if (!rtn_flag)
            {
              assert(false);
              return false;
            }
            pc.push_back(pc_temp);

            return true;
          }

        private:
          cst_airfoil_type cst;
      };
    }
  }
}
#endif
