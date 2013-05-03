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

            assert((tt>=0) && (tt<=1));

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

            assert((tt>=0) && (tt<=1));

            typename curve__::point_type fp(pc->fp(tt));
            typename curve__::data_type rtn(fp.dot(fp)+pc->fpp(tt).dot(pc->f(tt)-pt));
            typename curve__::tolerance_type tol;

            if (tol.approximately_equal(rtn, 0))
            {
              curve_g_functor<curve__> g;

              g.pc=pc;
              g.pt=pt;
              if (t>=1)
              {
                rtn=(g(1)-g(0.99))/(0.01);
              }
              else if (t<=0)
              {
                rtn=(g(0.01)-g(0))/(0.01);
              }
              else
              {
                rtn=(g(t+0.01)-g(t))/(0.01);
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
        nrm.set_max_iteration(100);
        if (c.open())
        {
          nrm.set_lower_condition(0, eli::mutil::nls::newton_raphson_constrained_method<typename curve__::data_type>::NRC_EXCLUSIVE);
          nrm.set_upper_condition(1, eli::mutil::nls::newton_raphson_constrained_method<typename curve__::data_type>::NRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(0, 1);
        }

        // set the initial guess
        nrm.set_initial_guess(t0);
        dist0=eli::geom::point::distance(c.f(t0), pt);

        // find the root
        stat = nrm.find_root(t, g, gp, 0);

        // if found root and it is within bounds and is closer than initial guess
        if (stat==eli::mutil::nls::newton_raphson_method<typename curve__::data_type>::converged)
        {
          assert((t>=0) && (t<=1));

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
        typename curve__::data_type dist, tt, dd;

        // first check is start, middle and (if needed) end points
        {
          t=0;
          dist=eli::geom::point::distance(c.f(t), pt);
          tt=0.5;
          dd=eli::geom::point::distance(c.f(tt), pt);
          if (dd<dist)
          {
            t=tt;
            dist=dd;
          }
          if (c.open())
          {
            tt=1;
            dd=eli::geom::point::distance(c.f(tt), pt);
            if (dd<dist)
            {
              t=tt;
              dist=dd;
            }
          }
        }
        cand_pair.first=t;
        cand_pair.second=dist;
        tinit.push_back(cand_pair);


        // need to pick initial guesses
        typename curve__::index_type i, deg(c.degree()), ssize;
        std::vector<typename curve__::data_type> tsample(2*deg+1);
        typename curve__::point_type p0, p1;
        typename curve__::data_type temp, tlen;

        // determine the sample parameters from the control polygon points
        ssize=tsample.size();
        i=0;
        p1=c.get_control_point(i);
        tsample[i]=0;
        for (++i; i<=deg; ++i)
        {
          p0=p1;
          p1=c.get_control_point(i);
          temp=eli::geom::point::distance(p0, p1)/2;
          tsample[2*i-1]=tsample[2*i-2]+temp;
          tsample[2*i]=tsample[2*i-1]+temp;
        }
        tlen=tsample[tsample.size()-1];

        // add points that are minimums
        {
          // find candidate starting locations using distance between sampled points on curve and point
          for (i=0; i<ssize; ++i)
          {
            temp=eli::geom::point::distance(c.f(tsample[i]/tlen), pt);
//             std::cout << "point #=" << i << "\tdist_temp=" << temp << std::endl;
            if (temp<=1.01*dist)
            {
              cand_pair.first=tsample[i]/tlen;
              cand_pair.second=temp;
              tinit.push_back(cand_pair);
              if (temp<dist)
              {
                t=cand_pair.first;
                dist=cand_pair.second;
                it=tinit.begin();
                while (it!=tinit.end())
                {
                  // check to see if distance is beyond new threshold and remove if so
                  if (it->second>1.01*dist)
                  {
                    it=tinit.erase(it);
                  }
                  else
                  {
                    ++it;
                  }
                }
              }
//               std::cout << "% added point #=" << i << "\twith t=" << tsample[i]/tlen << std::endl;
            }
          }
        }

//         std::cout << "# t guesses=" << tinit.size() << std::endl;

        // make sure have some solutions to iterate
//         tinit.push_back(0);
//         if (c.open())
//         {
//           tinit.push_back(1);
//         }

        // cycle through all possible minima to find best
        for (it=tinit.begin(); it!=tinit.end(); ++it)
        {
          dd=minimum_distance(tt, c, pt, it->first);
//           std::cout << "% completed root starting at" << *it << std::endl;

          assert((tt>=0) && (tt<=1));

          dd=eli::geom::point::distance(c.f(tt), pt);

          // check to see if is closer than previous minimum
          if (dd<dist)
          {
            t=tt;
            dist=dd;
          }

//               std::cout << "# dd=" << dd << std::endl;
//               std::cout << "# j=" << j << "\tnj=" << tinit.size() << std::endl;
        }

//         std::cout << "# returning dist=" << dist << std::endl;
        return dist;
      }
    }
  }
}
#endif
