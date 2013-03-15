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

#ifndef eli_geom_curve_piecewise_circle_creator_hpp
#define eli_geom_curve_piecewise_circle_creator_hpp

#include <iterator>

#include "eli/geom/point/distance.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename point__, typename data__, typename tol__>
      bool create_circle_2(piecewise<bezier, data__, 2, tol__> &/*pc*/, const point__ &/*start*/, const point__ &/*origin*/)
      {
        // NOT IMPLEMENTED
        assert(false);
        return false;
      }

      template<typename point__, typename data__, typename tol__>
      bool create_circle_2(piecewise<bezier, data__, 2, tol__> &/*pc*/, const point__ &/*start*/, const point__ &/*middle*/, const point__ &/*last*/)
      {
        // NOT IMPLEMENTED
        assert(false);
        return false;
      }

      template<typename point__, typename data__, typename tol__>
      bool create_circular_arc_2(piecewise<bezier, data__, 2, tol__> &/*pc*/, const point__ &/*start*/, const point__ &/*origin*/, const point__ &/*normal*/, const data__ &/*angle*/)
      {
        // NOT IMPLEMENTED
        assert(false);
        return false;
      }

      template<typename point__, typename data__, typename tol__>
      bool create_circular_arc_2(piecewise<bezier, data__, 2, tol__> &/*pc*/, const point__ &/*start*/, const point__ &/*middle*/, const point__ &/*last*/)
      {
        // NOT IMPLEMENTED
        assert(false);
        return false;
      }

      template<typename point__, typename data__, typename tol__>
      bool create_circle_3(piecewise<bezier, data__, 3, tol__> &pc, const point__ &start, const point__ &origin, const point__ &normal)
      {
        typedef piecewise<bezier, data__, 3, tol__> piecewise_curve_type;
        typedef typename piecewise<bezier, data__, 3, tol__>::curve_type curve_type;
        typedef typename curve_type::control_point_type control_point_type;
        typedef typename curve_type::point_type point_type;
        typedef typename curve_type::data_type data_type;
        typedef typename curve_type::index_type index_type;
        typename piecewise_curve_type::error_code err;

        point_type x, y;
        data_type r, k;
        index_type i;

        // create local x and y vectors
        eli::geom::point::distance(r, start, origin);
        k=4*(eli::constants::math<data_type>::sqrt_two()-1)/3;

        pc.clear();

        tol__ tol;
        if (tol.approximately_equal(r, 0))
        {
          r=0;
          x << 1, 0, 0;
          y << 0, 1, 0;
        }
        else
        {
          x=start-origin;
          x.normalize();
          y=normal.cross(x);
          y.normalize();

#ifdef DEBUG
          point_type n(normal);
          n.normalize();
          assert(tol.approximately_equal(x.dot(n), 0));
#endif
        }

        curve_type c(3);
        control_point_type cp[4];

        // set 1st quadrant curve
        cp[0]=r*x+origin;
        cp[1]=r*(x+k*y)+origin;
        cp[2]=r*(k*x+y)+origin;
        cp[3]=r*y+origin;
        for (i=0; i<4; ++i)
        {
          c.set_control_point(cp[i], i);
        }
        err=pc.push_back(c, 0.25);
        if (err!=piecewise_curve_type::NO_ERROR)
          return false;

        // set 2nd quadrant curve
        cp[0]=r*y+origin;
        cp[1]=r*(y-k*x)+origin;
        cp[2]=r*(k*y-x)+origin;
        cp[3]=-r*x+origin;
        for (i=0; i<4; ++i)
        {
          c.set_control_point(cp[i], i);
        }
        err=pc.push_back(c, 0.25);
        if (err!=piecewise_curve_type::NO_ERROR)
          return false;

        // set 3rd quadrant curve
        cp[0]=-r*x+origin;
        cp[1]=-r*(x+k*y)+origin;
        cp[2]=-r*(k*x+y)+origin;
        cp[3]=-r*y+origin;
        for (i=0; i<4; ++i)
        {
          c.set_control_point(cp[i], i);
        }
        err=pc.push_back(c, 0.25);
        if (err!=piecewise_curve_type::NO_ERROR)
          return false;

        // set 4th quadrant curve
        cp[0]=-r*y+origin;
        cp[1]=-r*(y-k*x)+origin;
        cp[2]=-r*(k*y-x)+origin;
        cp[3]=r*x+origin;
        for (i=0; i<4; ++i)
        {
          c.set_control_point(cp[i], i);
        }
        err=pc.push_back(c, 0.25);
        if (err!=piecewise_curve_type::NO_ERROR)
          return false;

        return true;
      }

      template<typename point__, typename data__, typename tol__>
      bool create_circle_3(piecewise<bezier, data__, 2, tol__> &/*pc*/, const point__ &/*start*/, const point__ &/*middle*/, const point__ &/*last*/)
      {
        // NOT IMPLEMENTED
        assert(false);
        return false;
      }

      template<typename point__, typename data__, typename tol__>
      bool create_circular_arc_3(piecewise<bezier, data__, 3, tol__> &/*pc*/, const point__ &/*start*/, const point__ &/*origin*/, const point__ &/*normal*/, const data__ &/*angle*/)
      {
        // NOT IMPLEMENTED
        assert(false);
        return false;
      }

      template<typename point__, typename data__, typename tol__>
      bool create_circular_arc_3(piecewise<bezier, data__, 3, tol__> &/*pc*/, const point__ &/*start*/, const point__ &/*middle*/, const point__ &/*last*/)
      {
        // NOT IMPLEMENTED
        assert(false);
        return false;
      }
    }
  }
}
#endif
