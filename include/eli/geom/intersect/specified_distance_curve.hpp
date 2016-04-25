/*********************************************************************************
* Copyright (c) 2015 Rob McDonald <ramcdona@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - initial code and implementation
********************************************************************************/

#ifndef eli_geom_intersect_specified_distance_curve_hpp
#define eli_geom_intersect_specified_distance_curve_hpp

#include <cmath>
#include <vector>
#include <list>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/iterative_root_base_constrained.hpp"
#include "eli/mutil/nls/bisection_method.hpp"

#include "eli/geom/point/distance.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/intersect/minimum_distance_bounding_box.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__ >
      class piecewise;
    }

    namespace intersect
    {
      namespace internal
      {
        template <typename curve__>
        struct curve_spec_g_functor
        {
          const curve__ *pc;
          typename curve__::point_type pt;
          typename curve__::data_type r0;

          typename curve__::data_type operator()(const typename curve__::data_type &t) const
          {
            typename curve__::data_type tt(t);

            if ( !(tt>=pc->get_t0()) )
            {
              std::cout << "Specified distance curve g_functor, tt less than minimum.  tt: " << tt << " t0: " << pc->get_t0() << std::endl;
              tt=pc->get_t0();
            }
            if ( !(tt<=pc->get_tmax()) )
            {
              std::cout << "Specified distance curve g_functor, tt greater than maximum.  tt: " << tt << " tmax: " << pc->get_tmax() << std::endl;
              tt=pc->get_tmax();
            }

            assert((tt>=pc->get_t0()) && (tt<=pc->get_tmax()));

            typename curve__::point_type u = pc->f(tt)-pt;

            return u.dot(u) - r0*r0;
          }
        };

        template <typename curve__>
        struct curve_spec_gp_functor
        {
          const curve__ *pc;
          typename curve__::point_type pt;
          typename curve__::data_type r0;

          typename curve__::data_type operator()(const typename curve__::data_type &t) const
          {
            typename curve__::data_type tt(t);

            if ( !(tt>=pc->get_t0()) )
            {
              std::cout << "Specified distance curve gp_functor, tt less than minimum.  tt: " << tt << " t0: " << pc->get_t0() << std::endl;
              tt=pc->get_t0();
            }
            if ( !(tt<=pc->get_tmax()) )
            {
              std::cout << "Specified distance curve gp_functor, tt greater than maximum.  tt: " << tt << " tmax: " << pc->get_tmax() << std::endl;
              tt=pc->get_tmax();
            }

            assert((tt>=pc->get_t0()) && (tt<=pc->get_tmax()));

            typename curve__::point_type u = pc->f(tt)-pt;
            typename curve__::point_type du = pc->fp(tt);

            return 2.0 * u.dot(du);
          }
        };
      }

      template<typename curve__>
      typename curve__::data_type specified_distance(typename curve__::data_type &t, const curve__ &c, const typename curve__::point_type &pt, const typename curve__::data_type &r0, const typename curve__::data_type &t0)
      {
        eli::mutil::nls::newton_raphson_method<typename curve__::data_type> nrm;
        internal::curve_spec_g_functor<curve__> g;
        internal::curve_spec_gp_functor<curve__> gp;
        typename curve__::data_type dist0, dist;
        typename curve__::tolerance_type tol;

        // setup the functors
        g.pc=&c;
        g.pt=pt;
        g.r0=r0;
        gp.pc=&c;
        gp.pt=pt;
        gp.r0=r0;

        // setup the solver
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(10);
        if (c.open())
        {
          nrm.set_lower_condition(c.get_t0(), eli::mutil::nls::newton_raphson_method<typename curve__::data_type>::IRC_EXCLUSIVE);
          nrm.set_upper_condition(c.get_tmax(), eli::mutil::nls::newton_raphson_method<typename curve__::data_type>::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(c.get_t0(), c.get_tmax());
        }

        // set the initial guess
        nrm.set_initial_guess(t0);
        dist0=eli::geom::point::distance(c.f(t0), pt)-r0;

        // find the root
        nrm.find_root(t, g, gp, 0);

        // if root is within bounds and is closer than initial guess
        {
          assert((t>=c.get_t0()) && (t<=c.get_tmax()));

          dist = eli::geom::point::distance(c.f(t), pt)-r0;
          if ( abs(dist) <= abs(dist0) )
          {
            return dist;
          }
        }

        // couldn't find better answer so return initial guess
        t=t0;
        return dist0;
      }

      template<typename curve__>
      typename curve__::data_type specified_distance(typename curve__::data_type &t, const curve__ &c, const typename curve__::point_type &pt, const typename curve__::data_type &r0, const typename curve__::data_type &tmin, const typename curve__::data_type &tmax )
      {
        eli::mutil::nls::bisection_method<typename curve__::data_type> bm;
        internal::curve_spec_g_functor<curve__> g;
        typename curve__::data_type t0, dist0, dist;
        typename curve__::tolerance_type tol;

        // setup the functors
        g.pc=&c;
        g.pt=pt;
        g.r0=r0;

        // setup the solver
        bm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        bm.set_max_iteration(20);
        bm.set_bounds( tmin, tmax );

        // set the initial guess
        t0 = 0.5 * ( c.get_t0() + c.get_tmax() );
        dist0=eli::geom::point::distance(c.f(t0), pt)-r0;

        // find the root
        bm.find_root(t, g, 0);

        // if root is within bounds and is closer than initial guess
        {
          assert((t>=c.get_t0()) && (t<=c.get_tmax()));

          dist = eli::geom::point::distance(c.f(t), pt)-r0;
          if ( abs(dist) <= abs(dist0) )
          {
            return dist;
          }
        }

        // couldn't find better answer so return initial guess
        t=t0;
        return dist0;
      }

      template<typename curve__>
      typename curve__::data_type specified_distance(typename curve__::data_type &t, const curve__ &c, const typename curve__::point_type &pt, const typename curve__::data_type &r0)
      {
        return specified_distance( t, c, pt, r0, c.get_t0(), c.get_tmax() );
      }

    }
  }
}
#endif
