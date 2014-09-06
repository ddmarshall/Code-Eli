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

#ifndef eli_geom_curve_piecewise_polynomial_creator_hpp
#define eli_geom_curve_piecewise_polynomial_creator_hpp

#include <vector>

#include "eli/code_eli.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/pseudo/polynomial.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_polynomial_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef unsigned short dimension_type;
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef eli::geom::curve::pseudo::polynomial<data_type, dim__> polynomial_type;
          typedef typename polynomial_type::coefficient_type polynomial_coefficient_type;

          piecewise_polynomial_creator() : piecewise_creator_base<data_type, dim__, tolerance_type>(1, 0) {}
          piecewise_polynomial_creator(const piecewise_polynomial_creator<data_type, dim__, tolerance_type> &pp)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(pp) {}
          ~piecewise_polynomial_creator() {}

          void set_polynomial(const polynomial_type &p)
          {
            poly=p;
          }

          void get_polynomial(polynomial_type &p) const
          {
            p=poly;
          }

          bool set_conditions(const polynomial_type &p)
          {
            poly=p;

            return true;
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_point_collection_type;
            typedef typename piecewise_curve_type::error_code error_code;

            curve_type c;
            control_point_collection_type cp, coeff;
            polynomial_coefficient_type a[dim__];
            index_type deg[dim__], bez_deg;

            // make sure only have one segment
            if (this->get_number_segments()!=1)
            {
              assert(false);
              return false;
            }

            // extract all coefficients and degrees
            bez_deg=0;
            for (dimension_type j=0; j<dim__; ++j)
            {
              deg[j]=poly.degree(j);
              if (deg[j]>bez_deg)
              {
                bez_deg=deg[j];
              }
              poly.get_coefficients(a[j], j);
            }

            // build the coefficient matrix needed for conversion
            coeff.resize(bez_deg+1, dim__);
            coeff.setZero();
            for (dimension_type j=0; j<dim__; ++j)
            {
              for (index_type i=0; i<=deg[j]; ++i)
              {
                coeff(i, j)=a[j][i];
              }
            }

            // build the control points
            cp.resize(bez_deg+1, dim__);
            eli::geom::utility::monomial_to_bezier_control_points(cp, coeff);

            // create the curve
            c.resize(bez_deg);
            for (index_type i=0; i<=bez_deg; ++i)
            {
              c.set_control_point(cp.row(i), i);
            }

            // set the piecewise curve
            pc.clear();
            pc.set_t0(this->get_t0());
            error_code err = pc.push_back(c, this->get_segment_dt(0));

            if (err!=piecewise_curve_type::NO_ERRORS)
            {
              assert(false);
              return false;
            }

            return true;
          }

        private:
          polynomial_type poly;
      };
    }
  }
}
#endif
