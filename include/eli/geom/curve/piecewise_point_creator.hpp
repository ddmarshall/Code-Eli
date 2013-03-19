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

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__>
      class piecewise_point_creator
      {
        public:
          typedef data__  data_type;
          typedef int index_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;

          piecewise_point_creator() : dt(4), t0(0)
          {
            set_t0(0);
            for (index_type i=0; i<static_cast<index_type>(dt.size()); ++i)
              dt[i]=1;
          }
          piecewise_point_creator(const index_type &ns) : dt(ns), t0(0)
          {
            set_t0(0);
            for (index_type i=0; i<ns; ++i)
              dt[i]=1;
          }
          piecewise_point_creator(const piecewise_point_creator<data__, dim__> &ppc)
            : point(ppc.point), dt(ppc.dt), t0(ppc.t0) {}
          ~piecewise_point_creator() {}

          index_type get_number_segments() const
          {
            return static_cast<index_type>(dt.size());
          }
          void set_number_segments(const index_type &ns)
          {
            dt.resize(ns);
            for (index_type i=0; i<ns; ++i)
              dt[i]=1;
          }

          void set_t0(const data_type &tt0) {t0=tt0;}
          data_type get_t0() const {return t0;}

          void set_segment_dt(const data_type &dtt, const index_type &i)
          {
            if ((dtt>0) && (i>=0) && (i<static_cast<index_type>(dt.size())))
              dt[i]=dtt;
            else
              assert(false);
          }
          data_type get_segment_dt(const index_type &i) const
          {
            if ((i<0) || (i>=static_cast<index_type>(dt.size())))
            {
              assert(false);
              return static_cast<data_type>(-1);
            }

            return dt[i];
          }

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
            pc.set_t0(t0);

            // set the control points
            for (index_type i=0; i<static_cast<index_type>(dt.size()); ++i)
            {
              c.set_control_point(point, 0);
              c.set_control_point(point, 1);
              pc.push_back(c, dt[i]);
            }

            assert(pc.closed());
            return true;
          }

        private:
          point_type point;
          std::vector<data_type> dt;
          data_type t0;
      };
    }
  }
}
#endif

