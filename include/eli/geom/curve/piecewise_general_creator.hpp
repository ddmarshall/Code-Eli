/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_curve_piecewise_general_creator_hpp
#define eli_geom_curve_piecewise_general_creator_hpp

#include <iterator>

#include "eli/constants/math.hpp"

#include "eli/geom/point/distance.hpp"

#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {

      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_general_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;

        public:
          piecewise_general_creator() : piecewise_creator_base<data_type, dim__, tolerance_type>(4, 0)
          {
          }
          piecewise_general_creator(const index_type &ns)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(ns, 0)
          {
          }
          piecewise_general_creator(const piecewise_ellipse_creator_base<data_type, dim__, tolerance_type> &pcc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(pcc) {}
          virtual ~piecewise_general_creator() {};

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;
            typedef typename curve_type::control_point_type control_point_type;

            pc.clear();


            return false;
          }
      };
    }
  }
}
#endif
