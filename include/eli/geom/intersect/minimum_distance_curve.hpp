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
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/newton_raphson_constrained_method.hpp"

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
        struct curve_g_functor
        {
          const curve__ *pc;
          typename curve__::point_type pt;

          typename curve__::data_type operator()(const typename curve__::data_type &t) const
          {
            typename curve__::data_type tt(t);

            if ( !(tt>=pc->get_t0()) )
            {
              std::cout << "Minimum distance curve g_functor, tt less than minimum.  tt: " << tt << " t0: " << pc->get_t0() << std::endl;
              tt=pc->get_t0();
            }
            if ( !(tt<=pc->get_tmax()) )
            {
              std::cout << "Minimum distance curve g_functor, tt greater than maximum.  tt: " << tt << " tmax: " << pc->get_tmax() << std::endl;
              tt=pc->get_tmax();
            }

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

            if ( !(tt>=pc->get_t0()) )
            {
              std::cout << "Minimum distance curve gp_functor, tt less than minimum.  tt: " << tt << " t0: " << pc->get_t0() << std::endl;
              tt=pc->get_t0();
            }
            if ( !(tt<=pc->get_tmax()) )
            {
              std::cout << "Minimum distance curve gp_functor, tt greater than maximum.  tt: " << tt << " tmax: " << pc->get_tmax() << std::endl;
              tt=pc->get_tmax();
            }

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
                rtn=(g(pc->get_t0()+static_cast<typename curve__::data_type>(0.01))-g(pc->get_t0()))/static_cast<typename curve__::data_type>(0.01);
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
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
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
        nrm.find_root(t, g, gp, 0);

        // if root is within bounds and is closer than initial guess
        {
          assert((t>=c.get_t0()) && (t<=c.get_tmax()));

          dist = eli::geom::point::distance(c.f(t), pt);
          if  (dist<=dist0)
          {
            return dist;
          }
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

      template< typename first__, typename second__>
      bool pairfirstcompare( const std::pair < first__, second__ > &a, const std::pair < first__, second__ > &b )
      {
          return ( a.first < b.first );
      }

      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__>
      typename curve::piecewise<curve__, data__, dim__, tol__>::data_type minimum_distance(typename curve::piecewise<curve__, data__, dim__, tol__>::data_type &t,
                                                                                    const curve::piecewise<curve__, data__, dim__, tol__> &pc,
                                                                                    const typename curve::piecewise<curve__, data__, dim__, tol__>::point_type &pt)
      {
        typedef curve::piecewise<curve__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::curve_type curve_type;
        typedef typename piecewise_type::data_type data_type;
        typedef typename piecewise_type::bounding_box_type bounding_box_type;

        typedef typename piecewise_type::segment_collection_type::const_iterator segit;

        typedef std::vector< std::pair<data_type,segit> > dvec;
        dvec minbbdist;

        // Find closest corner of bounding boxes, add them to vector
        // Simple linear search, would be more efficient with some sort of tree.
        for (segit seg=pc.segments.begin(); seg!=pc.segments.end(); ++seg)
        {
          bounding_box_type bb_local;
          seg->second.get_bounding_box(bb_local);

          data_type dbbmin;
          dbbmin = minimum_distance(bb_local, pt);

          minbbdist.push_back(std::make_pair(dbbmin, seg));
        }

        // Sort by nearest distance.
        std::sort( minbbdist.begin(), minbbdist.end(), pairfirstcompare<data_type, segit> );

        // Iterate over segments, starting with nearest bounding box
        data_type dmin(std::numeric_limits<data_type>::max());
        typename dvec::const_iterator it;
        for (it=minbbdist.begin(); it!=minbbdist.end(); ++it)
        {
          // If nearest bb distance is farther than current best, we're done.
          if(it->first < dmin )
          {
            segit seg = it->second;

            curve_type c(seg->second);

            data_type tlocal, d;
            d=minimum_distance(tlocal,c,pt);

            if(d < dmin)
            {
              data_type tstart(seg->first);
              data_type dt(pc.get_delta_t(seg));

              dmin = d;
              t=tstart+tlocal*dt;
            }
          }
          else
          {
            break;
          }

        }
        return dmin;
      }
    }
  }
}
#endif
