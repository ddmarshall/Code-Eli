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

#ifndef eli_geom_intersect_minimum_distance_curve_hpp
#define eli_geom_intersect_minimum_distance_curve_hpp

#include <cmath>
#include <vector>
#include <list>

#ifdef Success  // X11 #define collides with Eigen
#undef Success
#endif

#include "Eigen/Eigen"

#include "eli/mutil/nls/newton_raphson_constrained_method.hpp"

#include "eli/geom/point/distance.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      namespace internal
      {
        template <typename curve__>
        struct curve_g_functor
        {
          const curve__ *pc;
          typename curve__::point_type pt;

          typename curve__::data_type operator()(const typename curve__::data_type &t) const
          {
            typename curve__::data_type tt(t);

            assert((tt>=pc->get_t0()) && (tt<=pc->get_tmax()));

            return (pc->f(tt)-pt).dot(pc->fp(tt));
          }
        };

        template <typename curve__>
        struct curve_gp_functor
        {
          const curve__ *pc;
          typename curve__::point_type pt;

          typename curve__::data_type operator()(const typename curve__::data_type &t) const
          {
            typename curve__::data_type tt(t);

            assert((tt>=pc->get_t0()) && (tt<=pc->get_tmax()));

            typename curve__::point_type fp(pc->fp(tt));
            typename curve__::data_type rtn(fp.dot(fp)+pc->fpp(tt).dot(pc->f(tt)-pt));
            typename curve__::tolerance_type tol;

            if (tol.approximately_equal(rtn, 0))
            {
              curve_g_functor<curve__> g;

              g.pc=pc;
              g.pt=pt;
              if (t>=pc->get_tmax())
              {
                rtn=(g(pc->get_tmax())-g(static_cast<typename curve__::data_type>(pc->get_tmax()-.01)))/static_cast<typename curve__::data_type>(0.01);
              }
              else if (t<=pc->get_t0())
              {
                rtn=(g(pc->get_t0()+.01)-g(pc->get_t0()))/static_cast<typename curve__::data_type>(0.01);
              }
              else
              {
                rtn=(g(t+static_cast<typename curve__::data_type>(0.01))-g(t))/static_cast<typename curve__::data_type>(0.01);
              }
            }

            return rtn;
          }
        };
      }

      template<typename curve__>
      typename curve__::data_type minimum_distance(typename curve__::data_type &t, const curve__ &c, const typename curve__::point_type &pt, const typename curve__::data_type &t0)
      {
        eli::mutil::nls::newton_raphson_constrained_method<typename curve__::data_type> nrm;
        int stat;
        internal::curve_g_functor<curve__> g;
        internal::curve_gp_functor<curve__> gp;
        typename curve__::data_type dist0, dist;
        typename curve__::tolerance_type tol;

        // setup the functors
        g.pc=&c;
        g.pt=pt;
        gp.pc=&c;
        gp.pt=pt;

        // setup the solver
        nrm.set_absolute_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(10);
        if (c.open())
        {
          nrm.set_lower_condition(c.get_t0(), eli::mutil::nls::newton_raphson_constrained_method<typename curve__::data_type>::NRC_EXCLUSIVE);
          nrm.set_upper_condition(c.get_tmax(), eli::mutil::nls::newton_raphson_constrained_method<typename curve__::data_type>::NRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(c.get_t0(), c.get_tmax());
        }

        // set the initial guess
        nrm.set_initial_guess(t0);
        dist0=eli::geom::point::distance(c.f(t0), pt);

        // find the root
        stat = nrm.find_root(t, g, gp, 0);

        // if found root and it is within bounds and is closer than initial guess
        if (stat==eli::mutil::nls::newton_raphson_method<typename curve__::data_type>::converged)
        {
          assert((t>=c.get_t0()) && (t<=c.get_tmax()));

          dist = eli::geom::point::distance(c.f(t), pt);
          if  (dist<=dist0)
          {
            return dist;
          }
        }
        else
        {
//             std::cout << "# not converged!" << std::endl;
        }

        // couldn't find better answer so return initial guess
        t=t0;
        return dist0;
      }

      template<typename curve__>
      typename curve__::data_type minimum_distance(typename curve__::data_type &t, const curve__ &c, const typename curve__::point_type &pt)
      {
        typename curve__::tolerance_type tol;
        std::list<std::pair<typename curve__::data_type, typename curve__::data_type>> tinit;
        typename std::list<std::pair<typename curve__::data_type, typename curve__::data_type>>::iterator it;
        std::pair<typename curve__::data_type, typename curve__::data_type> cand_pair;

        // possible that end points are closest, so start by checking them
        typename curve__::data_type dist, tt, dd, tspan;


        typename curve__::index_type i, n;

        // Just a guess
        n=2*c.degree()+1;
        tspan = c.get_tmax()-c.get_t0();

        // Evenly spaced in parameter, don't repeat 0/1 if closed curve.
        typename curve__::data_type dt;
        if (c.open())
        {
          dt = tspan/(n-1);
        }
        else
        {
          dt = tspan/n;
        }

        // Find closest of evenly spaced points.
        tt = c.get_t0();
        dist = std::numeric_limits<typename curve__::data_type>::max();
        for (i = 0; i < n; i++)
        {
          dd=eli::geom::point::distance(c.f(tt), pt);

          if( dd < dist )
          {
            t=tt;
            dist=dd;
          }
          tt+=dt;
          if( tt >= c.get_tmax() )
          {
            tt=c.get_tmax();
          }
        }

        // Polish best point with Newton's method search.
        typename curve__::data_type t0(t);
        dist=minimum_distance(t, c, pt, t0);

        return dist;
      }
    }
  }
}
#endif
