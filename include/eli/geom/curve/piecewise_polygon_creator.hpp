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

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__>
      class piecewise_polygon_creator
      {
        public:
          typedef data__  data_type;
          typedef int index_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;

          piecewise_polygon_creator() : corner(4), dt(4), t0(0)
          {
            set_t0(0);
            for (index_type i=0; i<static_cast<index_type>(corner.size()); ++i)
              dt[i]=1;
          }
          piecewise_polygon_creator(const index_type &ns) : corner(ns), dt(ns), t0(0)
          {
            set_t0(0);
            for (index_type i=0; i<ns; ++i)
              dt[i]=1;
          }
          piecewise_polygon_creator(const piecewise_polygon_creator<data__, dim__> &ppc)
            : corner(ppc.corner), dt(ppc.dt), t0(ppc.t0) {}
          ~piecewise_polygon_creator() {}

          index_type get_number_sides() const
          {
            assert(dt.size()==corner.size());

            return static_cast<index_type>(corner.size());
          }
          void set_number_sides(const index_type &ns)
          {
            corner.resize(ns);
            dt.resize(ns);
          }

          void set_t0(const data_type &tt0) {t0=tt0;}
          data_type get_t0() const {return t0;}

          void set_edge_dt(const data_type &dtt, const index_type &i)
          {
            if ((dtt>0) && (i>=0) && (i<static_cast<index_type>(dt.size())))
              dt[i]=dtt;
            else
              assert(false);
          }
          data_type get_edge_dt(const index_type &i) const
          {
            if ((i<0) || (i>=static_cast<index_type>(dt.size())))
            {
              assert(false);
              return static_cast<data_type>(-1);
            }

            return dt[i];
          }

          void set_corner(const point_type &c, const index_type &i)
          {
            if ((i>=0) && (i<static_cast<index_type>(corner.size())))
              corner[i]=c;
            else
              assert(false);
          }
          point_type get_corner(const index_type &i) const
          {
            if ((i<0) || (i>=static_cast<index_type>(dt.size())))
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

            // do sanity check
            if (corner.size()!=dt.size())
            {
              assert(false);
              return false;
            }

            // set the start parameter
            pc.set_t0(t0);

            // set the first n-1 edges
            for (index_type i=0; i<static_cast<index_type>(dt.size()-1); ++i)
            {
              c.set_control_point(corner[i], 0);
              c.set_control_point(corner[i+1], 1);
              pc.push_back(c, dt[i]);
            }

            // set the last edge
            c.set_control_point(corner[corner.size()-1], 0);
            c.set_control_point(corner[0], 1);
            pc.push_back(c, dt[dt.size()-1]);

            assert(pc.closed());
            return true;
          }

        private:
          std::vector<point_type, Eigen::aligned_allocator<point_type>> corner;
          std::vector<data_type> dt;
          data_type t0;
      };
    }
  }
}
#endif
