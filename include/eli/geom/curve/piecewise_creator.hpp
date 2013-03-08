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

#ifndef eli_geom_curve_piecewise_creator_hpp
#define eli_geom_curve_piecewise_creator_hpp

#include <list>
#include <iterator>

#include "eli/util/tolerance.hpp"

#include "eli/mutil/fd/d1o2.hpp"

#include "eli/geom/point/distance.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"
namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename point_it__, typename dt_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_cubic_hermite_interpolating_polynomials(piecewise<bezier, data__, dim__, tol__> &pc, point_it__ itb, point_it__ ite, dt_it__ ditb)
      {
        typedef piecewise<bezier, data__, dim__, tol__> piecewise_curve_type;
        typedef typename piecewise<bezier, data__, dim__, tol__>::curve_type curve_type;
        typedef typename curve_type::data_type data_type;
        typedef typename curve_type::index_type index_type;

        index_type i, j, npts=static_cast<index_type>(std::distance(itb, ite));
        typename piecewise_curve_type::error_code err;

        if (npts<3)
          return false;

        // cycle through the points and create the segments
        point_it__ it(itb), itm1(itb), itp1(itb);
        dt_it__ dit(ditb), ditm1(ditb), ditp1(ditb);
        curve_type c(3);
        data_type dt, tmp[3];
        typename curve_type::control_point_type cp[4], m[2];
        eli::mutil::fd::d1o2<data_type> d1approx;

        // set the starting time
        pc.set_t0(*dit);

        // need to do first segment separately
        ++it; ++itp1; ++itp1; ++dit; ++ditp1; ++ditp1;
        dt=(*dit)-(*ditm1);
        if (dt<=0)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }
        d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::RIGHT);
        for (j=0; j<dim__; ++j)
        {
          tmp[0]=(*itm1)(j);
          tmp[1]=(*it)(j);
          tmp[2]=(*itp1)(j);
          d1approx.evaluate(m[0](j), tmp, ditm1);
        }
        d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::CENTER);
        for (j=0; j<dim__; ++j)
        {
          tmp[0]=(*itm1)(j);
          tmp[1]=(*it)(j);
          tmp[2]=(*itp1)(j);
          d1approx.evaluate(m[1](j), tmp, ditm1);
        }
        cp[0]=*itm1;
        cp[3]=*it;
        cp[1]=cp[0]+m[0]/3;
        cp[2]=cp[3]-m[1]/3;
        for (j=0; j<4; ++j)
        {
          c.set_control_point(cp[j], j);
        }
        pc.push_back(c, dt);
        m[0]=m[1];

        // do all interior segments
        for (i=1; i<npts-2; ++i, ++itm1, ++it, ++itp1, ++ditm1, ++dit, ++ditp1)
        {
          dt=(*ditp1)-(*dit);
          if (dt<=0)
          {
            pc.clear();
            pc.set_t0(0);
            return false;
          }
          for (j=0; j<dim__; ++j)
          {
            tmp[0]=(*itm1)(j);
            tmp[1]=(*it)(j);
            tmp[2]=(*itp1)(j);
            d1approx.evaluate(m[1](j), tmp, ditm1);
          }
          cp[0]=*it;
          cp[3]=*itp1;
          cp[1]=cp[0]+m[0]/3;
          cp[2]=cp[3]-m[1]/3;
          for (j=0; j<4; ++j)
          {
            c.set_control_point(cp[j], j);
          }
          pc.push_back(c, dt);
          m[0]=m[1];
        }

        // need to do last segment separately
        dt=(*ditp1)-(*dit);
        if (dt<=0)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }
        d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::LEFT);
        for (j=0; j<dim__; ++j)
        {
          tmp[0]=(*itm1)(j);
          tmp[1]=(*it)(j);
          tmp[2]=(*itp1)(j);
          d1approx.evaluate(m[1](j), tmp, ditm1);
        }
        cp[0]=*it;
        cp[3]=*itp1;
        cp[1]=cp[0]+m[0]/3;
        cp[2]=cp[3]-m[1]/3;
        for (j=0; j<4; ++j)
        {
          c.set_control_point(cp[j], j);
        }
        err=pc.push_back(c, dt);
        if (err!=piecewise_curve_type::NO_ERROR)
          return false;

        return true;
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

//       template<typename point_it__>
//       error_code create_piecewise_cubic(point_it__ itb, point_it__ ite)
//       {
//         size_type npts=std::distance(itb, ite);
//         std::vector<typename point_it__::T> dt(npts);
//         return create_piecewise_cubic(itb, ite, dt.begin());
//       }
    }
  }
}
#endif
