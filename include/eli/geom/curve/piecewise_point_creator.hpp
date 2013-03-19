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

#ifndef eli_geom_curve_piecewise_point_creator_hpp
#define eli_geom_curve_piecewise_point_creator_hpp

#include <vector>

#include "Eigen/Eigen"

#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__>
      class piecewise_point_creator: public piecewise_creator_base<data__>
      {
        public:
          typedef data__  data_type;
          typedef int index_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;

          piecewise_point_creator() : piecewise_creator_base<data__>(4, 0) {}
          piecewise_point_creator(const index_type &ns) : piecewise_creator_base<data__>(ns, 0) {}
          piecewise_point_creator(const piecewise_point_creator<data__, dim__> &ppc)
            : piecewise_creator_base<data__>(ppc), point(ppc.point) {}
          ~piecewise_point_creator() {}

          void set_point(const point_type &p)
          {
            point=p;
          }
          point_type get_point() const
          {
            return point;
          }

          template<typename tol__>
          bool create(piecewise<bezier, data_type, dim__, tol__> &pc) const
          {
            typename piecewise<bezier, data_type, dim__, tol__>::curve_type c(1);

            // set the start parameter
            pc.set_t0(this->get_t0());

            // set the control points
            for (index_type i=0; i<this->get_number_segments(); ++i)
            {
              c.set_control_point(point, 0);
              c.set_control_point(point, 1);
              pc.push_back(c, this->get_segment_dt(i));
            }

            assert(pc.closed());
            return true;
          }

        private:
          point_type point;
      };
    }
  }
}
#endif

