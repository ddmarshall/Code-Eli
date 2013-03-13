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

#ifndef eli_geom_curve_piecewise_spline_creator_hpp
#define eli_geom_curve_piecewise_spline_creator_hpp

#include <iterator>

#include "eli/mutil/fd/d1o2.hpp"

// #include "eli/geom/point/distance.hpp"

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

      template<typename point_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_cubic_spline(piecewise<bezier, data__, dim__, tol__> &/*pc*/, point_it__ /*itb*/, point_it__ /*ite*/)
      {
#if 0
        size_type npts=std::distance(itb, ite);
        std::vector<typename point_it__::T> dt(npts);
        return create_piecewise_cubic(itb, ite, dt.begin());
#else
        // NOT IMPLEMENTED
        assert(false);
        return false;
#endif
      }

      template<typename point_it__, typename dt_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_cardinal_spline(piecewise<bezier, data__, dim__, tol__> &/*pc*/, point_it__ /*itb*/, point_it__ /*ite*/, const data__ &c)
      {
        // check the parameter
        if ( (c<0) || (c>1) )
          return false;

        // NOT IMPLEMENTED
        assert(false);
        return false;
      }

      template<typename point_it__, typename dt_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_catmull_rom_spline(piecewise<bezier, data__, dim__, tol__> &pc, point_it__ itb, point_it__ ite)
      {
        return create_piecewise_cardinal_spline(pc, itb, ite, static_cast<data__>(0));
      }

      template<typename point_it__, typename dt_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_kochanek_bartels_spline(piecewise<bezier, data__, dim__, tol__> &/*pc*/, point_it__ /*itb*/, point_it__ /*ite*/, const data__ &tension, const data__ &bias, const data__ &continuity)
      {
        // check some parameters
        if ( (tension<0) || (tension>1) )
          return false;
        if ( (bias<0) || (bias>1) )
          return false;
        if ( (continuity<0) || (continuity>1) )
          return false;

        // NOT IMPLEMENTED
        assert(false);
        return false;
      }
    }
  }
}
#endif
