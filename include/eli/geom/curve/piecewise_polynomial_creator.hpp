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

          void get_curve(polynomial_type &p) const
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
#if 0
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename curve_type::control_point_type control_point_type;
            typedef typename piecewise_curve_type::error_code error_code;

            // make sure only have one segment
            if (this->get_number_segments()!=1)
            {
              assert(false);
              return false;
            }

            curve_type c;
            control_point_type cp;
            eb_control_point_type eb_cp;
            index_type i, deg(eb.degree());

            // create the control points for the x-dimension
            polynomial_type x_eb(1);

            eb_cp << 0;
            x_eb.set_control_point(eb_cp, 0);
            eb_cp << 1;
            x_eb.set_control_point(eb_cp, 1);
            x_eb.degree_promote_to(deg);

            // make sure have same degree for x- and y-coordinates
            if (x_eb.degree()!=deg)
            {
              assert(false);
              return false;
            }

            // cycle through each control point and add to bezier curve
            c.resize(deg);
            for (i=0; i<=deg; ++i)
            {
              cp.setZero();
              cp(0) = x_eb.get_control_point(i).x();
              cp(1) = eb.get_control_point(i).x();
              c.set_control_point(cp, i);
            }
            if (rev)
            {
              c.reverse();
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
#endif
            return false;
          }

        private:
          polynomial_type poly;
      };
    }
  }
}
#endif
