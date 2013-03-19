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

#ifndef eli_geom_curve_piecewise_polygon_creator_hpp
#define eli_geom_curve_piecewise_polygon_creator_hpp

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
      class piecewise_polygon_creator : public piecewise_creator_base<data__>
      {
        public:
          typedef data__  data_type;
          typedef int index_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;

          piecewise_polygon_creator() : piecewise_creator_base<data__>(4, 0), corner(4) {}
          piecewise_polygon_creator(const index_type &ns) : piecewise_creator_base<data__>(ns, 0), corner(ns) {}
          piecewise_polygon_creator(const piecewise_polygon_creator<data__, dim__> &ppc)
            : piecewise_creator_base<data__>(ppc), corner(ppc.corner) {}
          ~piecewise_polygon_creator() {}

          void set_corner(const point_type &c, const index_type &i)
          {
            if ((i>=0) && (i<static_cast<index_type>(corner.size())))
              corner[i]=c;
            else
              assert(false);
          }
          point_type get_corner(const index_type &i) const
          {
            if ((i<0) || (i>=static_cast<index_type>(corner.size())))
            {
              return corner[0];
              assert(false);
            }
            return corner[i];
          }

          template<typename tol__>
          bool create(piecewise<bezier, data_type, dim__, tol__> &pc) const
          {
            typename piecewise<bezier, data_type, dim__, tol__>::curve_type c(1);
            index_type nsegs(this->get_number_segments());

            // do sanity check
            if (corner.size()!=static_cast<size_t>(nsegs))
            {
              assert(false);
              return false;
            }

            // set the start parameter
            pc.set_t0(this->get_t0());

            // set the first n-1 edges
            for (index_type i=0; i<(nsegs-1); ++i)
            {
              c.set_control_point(corner[i], 0);
              c.set_control_point(corner[i+1], 1);
              pc.push_back(c, this->get_segment_dt(i));
            }

            // set the last edge
            c.set_control_point(corner[corner.size()-1], 0);
            c.set_control_point(corner[0], 1);
            pc.push_back(c, this->get_segment_dt(nsegs-1));

            assert(pc.closed());
            return true;
          }

        private:
          void number_segments_changed() {corner.resize(this->get_number_segments());}

        private:
          std::vector<point_type, Eigen::aligned_allocator<point_type>> corner;
      };
    }
  }
}
#endif
