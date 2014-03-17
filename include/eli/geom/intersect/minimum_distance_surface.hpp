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

#ifndef eli_geom_intersect_minimum_distance_surface_hpp
#define eli_geom_intersect_minimum_distance_surface_hpp

#include <cmath>
#include <vector>
#include <list>

#ifdef Success  // X11 #define collides with Eigen
#undef Success
#endif

#include "Eigen/Eigen"

#include "eli/mutil/nls/newton_raphson_constrained_system_method.hpp"

#include "eli/geom/intersect/minimum_distance_curve.hpp"
#include "eli/geom/point/distance.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      namespace internal
      {
        template <typename surface__>
        struct surface_g_functor
        {
          const surface__ *ps;
          typename surface__::point_type pt;
          typedef typename Eigen::Matrix<typename surface__::data_type, 2, 1> vec;

          vec operator()(const vec &u) const
          {
            typename surface__::data_type uu(u[0]), vv(u[1]);
            vec rtn;

            typename surface__::data_type umin, umax, vmin, vmax;
            ps->get_parameter_min(umin,vmin);
            ps->get_parameter_max(umax,vmax);

            assert((uu>=umin) && (uu<=umax));
            assert((vv>=vmin) && (vv<=vmax));

            uu=std::min(std::max(uu, static_cast<typename surface__::data_type>(umin)), static_cast<typename surface__::data_type>(umax));
            vv=std::min(std::max(vv, static_cast<typename surface__::data_type>(vmin)), static_cast<typename surface__::data_type>(vmax));

            typename surface__::point_type tmp;

            tmp=ps->f(uu, vv)-pt;
            rtn(0)=tmp.dot(ps->f_u(uu, vv));
            rtn(1)=tmp.dot(ps->f_v(uu, vv));
            return rtn;
          }
        };

        template <typename surface__>
        struct surface_gp_functor
        {
          const surface__ *ps;
          typename surface__::point_type pt;
          typedef typename Eigen::Matrix<typename surface__::data_type, 2, 1> vec;
          typedef typename Eigen::Matrix<typename surface__::data_type, 2, 2> mat;

          mat operator()(const vec &u) const
          {
            typename surface__::data_type uu(u[0]), vv(u[1]);
            mat rtn;

            typename surface__::data_type umin, umax, vmin, vmax;
            ps->get_parameter_min(umin,vmin);
            ps->get_parameter_max(umax,vmax);

            assert((uu>=umin) && (uu<=umax));
            assert((vv>=vmin) && (vv<=vmax));

            uu=std::min(std::max(uu, static_cast<typename surface__::data_type>(umin)), static_cast<typename surface__::data_type>(umax));
            vv=std::min(std::max(vv, static_cast<typename surface__::data_type>(vmin)), static_cast<typename surface__::data_type>(vmax));

            typename surface__::point_type tmp, Su, Sv, Suu, Suv, Svv;

            tmp=ps->f(uu, vv)-pt;
            Su=ps->f_u(uu, vv);
            Sv=ps->f_v(uu, vv);
            Suu=ps->f_uu(uu, vv);
            Suv=ps->f_uv(uu, vv);
            Svv=ps->f_vv(uu, vv);

            rtn(0,0)=Su.dot(Su)+tmp.dot(Suu);
            rtn(0,1)=Su.dot(Sv)+tmp.dot(Suv);
            rtn(1,0)=rtn(0,1);
            rtn(1,1)=Sv.dot(Sv)+tmp.dot(Svv);

            // TODO: What to do if matrix becomes singular?

            return rtn;
          }
        };
      }

      template<typename surface__>
      typename surface__::data_type minimum_distance(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt,
                                                     const typename surface__::data_type &u0, const typename surface__::data_type &v0)
      {
        typedef eli::mutil::nls::newton_raphson_constrained_system_method<typename surface__::data_type, 2, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        int stat;
        internal::surface_g_functor<surface__> g;
        internal::surface_gp_functor<surface__> gp;
        typename surface__::data_type dist0, dist;
        typename surface__::tolerance_type tol;

        typename surface__::data_type umin, umax, vmin, vmax;
        s.get_parameter_min(umin,vmin);
        s.get_parameter_max(umax,vmax);

        // setup the functors
        g.ps=&s;
        g.pt=pt;
        gp.ps=&s;
        gp.pt=pt;

        // setup the solver
        nrm.set_absolute_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(20);
        nrm.set_norm_type(nonlinear_solver_type::max_norm);
        if (s.open_u())
        {
          nrm.set_lower_condition(0, umin, nonlinear_solver_type::NRC_EXCLUSIVE);
          nrm.set_upper_condition(0, umax, nonlinear_solver_type::NRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(0, umin, umax);
        }
        if (s.open_v())
        {
          nrm.set_lower_condition(1, vmin, nonlinear_solver_type::NRC_EXCLUSIVE);
          nrm.set_upper_condition(1, vmax, nonlinear_solver_type::NRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(1, vmin, vmax);
        }

        // set the initial guess
        typename nonlinear_solver_type::solution_matrix uinit, rhs, ans;

        uinit(0)=u0;
        uinit(1)=v0;
        nrm.set_initial_guess(uinit);
        rhs.setZero();
        dist0=eli::geom::point::distance(s.f(u0, v0), pt);

        // find the root
        stat = nrm.find_root(ans, g, gp, rhs);
        u=ans(0);
        v=ans(1);

        // if root is within bounds and is closer than initial guess
        {
          assert((u>=umin) && (u<=umax));
          assert((v>=vmin) && (v<=vmax));

          dist = eli::geom::point::distance(s.f(u, v), pt);
          if  (dist<=dist0)
          {
            return dist;
          }
        }
//         else
//         {
//             std::cout << "% not converged";
//             if (stat==nonlinear_solver_type::hit_constraint)
//               std::cout << " because hit constraint" << std::endl;
//             else if (stat==nonlinear_solver_type::max_iteration)
//               std::cout << " reached max iteration" << std::endl;
//             else
//               std::cout << " for out of range parameters (" << ans(0) << ", " << ans(1) << ")" << std::endl;
//         }

        // couldn't find better answer so return initial guess
        u=u0;
        v=v0;
        return dist0;
      }

      namespace internal
      {
        template<typename data__>
        struct surface_point_uv_pairs
        {
          data__ u, v, dist;
        };
      }

      template<typename surface__>
      typename surface__::data_type minimum_distance(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt)
      {
        typename surface__::tolerance_type tol;
        std::list<internal::surface_point_uv_pairs<typename surface__::data_type>> uvinit;
        typename std::list<internal::surface_point_uv_pairs<typename surface__::data_type>>::iterator it;
        internal::surface_point_uv_pairs<typename surface__::data_type> cand_match;

        // possible that end points are closest, so start by checking them
        typename surface__::data_type dist, uu, vv, dd;

        typename surface__::data_type umin, umax, vmin, vmax;
        s.get_parameter_min(umin,vmin);
        s.get_parameter_max(umax,vmax);

        // first check start and (if needed) end points of edges and middle of surface
        {
          bool openu=s.open_u(), openv=s.open_v();
          u=umin;
          v=vmin;
          dist=eli::geom::point::distance(s.f(u, v), pt);
          uu=0.5*(umin+umax);
          vv=0.5*(vmin+vmax);
          dd=eli::geom::point::distance(s.f(uu, vv), pt);
          if (dd<dist)
          {
            u=uu;
            v=vv;
            dist=dd;
          }
          if (openu)
          {
            uu=umax;
            vv=vmin;
            dd=eli::geom::point::distance(s.f(uu, vv), pt);
            if (dd<dist)
            {
              u=uu;
              v=vv;
              dist=dd;
            }
            if (openv)
            {
              uu=umax;
              vv=vmax;
              dd=eli::geom::point::distance(s.f(uu, vv), pt);
              if (dd<dist)
              {
                u=uu;
                v=vv;
                dist=dd;
              }
            }
          }
          if (s.open_v())
          {
            uu=umax;
            vv=vmin;
            dd=eli::geom::point::distance(s.f(uu, vv), pt);
            if (dd<dist)
            {
              u=uu;
              v=vv;
              dist=dd;
            }
          }

          // check center point
          uu=0.5*(umin+umax);
          vv=0.5*(vmin+vmax);
          dd=eli::geom::point::distance(s.f(uu, vv), pt);
          if (dd<dist)
          {
            u=uu;
            v=vv;
            dist=dd;
          }
        }

        // take best match so far and add it to cases
        cand_match.u=u;
        cand_match.v=v;
        cand_match.dist=dist;
        uvinit.push_back(cand_match);

        // next first check center and edges
        {
          typename surface__::boundary_curve_type bc;

          // check u-min edge
          u=umin;
          s.get_uconst_curve(bc, u);
          dist=eli::geom::intersect::minimum_distance(v, bc, pt);

          // if open in u-direction check u-max edge
          if (s.open_u())
          {
            uu=umax;
            s.get_uconst_curve(bc, uu);
            dd=eli::geom::intersect::minimum_distance(vv, bc, pt);
            if (dd<dist)
            {
              u=uu;
              v=vv;
              dist=dd;
            }
          }
          // check v-min edge
          vv=vmin;
          s.get_vconst_curve(bc, vv);
          dd=eli::geom::intersect::minimum_distance(uu, bc, pt);
          if (dd<dist)
          {
            u=uu;
            v=vv;
            dist=dd;
          }

          // if open in v-direction check v-max edge
          if (s.open_v())
          {
            vv=vmax;
            s.get_vconst_curve(bc, vv);
            dd=eli::geom::intersect::minimum_distance(uu, bc, pt);
            if (dd<dist)
            {
              u=uu;
              v=vv;
              dist=dd;
            }
          }
        }

        // take best match so far and add it to cases
        cand_match.u=u;
        cand_match.v=v;
        cand_match.dist=dist;
        uvinit.push_back(cand_match);

        // need to pick initial guesses
        typename surface__::index_type i, j, degu(s.degree_u()), degv(s.degree_v()), susize, svsize;
        std::vector<typename surface__::data_type> usample(1*degu+1), vsample(1*degv+1);
        typename surface__::point_type p0, p1;
        typename surface__::data_type temp;

        // determine the sample parameters from the control polygon points
// TODO: need to calculate distance between control points for better distribution of samples
#if 0
        ssize=tsample.size();
        i=0;
        p1=c.get_control_point(i);
        tsample[i]=0;
        for (++i; i<=degu; ++i)
        {
          p0=p1;
          p1=c.get_control_point(i);
          temp=eli::geom::point::distance(p0, p1)/2;
          tsample[2*i-1]=tsample[2*i-2]+temp;
          tsample[2*i]=tsample[2*i-1]+temp;
        }
        tlen=tsample[tsample.size()-1];
#else
          susize=usample.size();
          svsize=vsample.size();
          for (i=0; i<susize; ++i)
          {
            usample[i]=static_cast<typename surface__::data_type>(i)/(susize-1);
          }
          for (j=0; j<svsize; ++j)
          {
            vsample[j]=static_cast<typename surface__::data_type>(j)/(svsize-1);
          }
#endif

        // add points that are minimums
        {
          // find candidate starting locations using distance between sampled points on curve and point
          for (i=0; i<susize; ++i)
          {
            for (j=0; j<svsize; ++j)
            {
              temp=eli::geom::point::distance(s.f(usample[i], vsample[j]), pt);
//             std::cout << "% point #=" << i << ", " << j << "\tdist_temp=" << temp << std::endl;
              if (temp<=1.01*dist)
              {
                cand_match.u=usample[i];
                cand_match.v=vsample[j];
                cand_match.dist=temp;
                uvinit.push_back(cand_match);
                if (temp<dist)
                {
                  u=cand_match.u;
                  v=cand_match.v;
                  dist=cand_match.dist;
                  it=uvinit.begin();
                  while (it!=uvinit.end())
                  {
                    // check to see if distance is beyond new threshold and remove if so
                    if (it->dist>1.01*dist)
                    {
                      it=uvinit.erase(it);
                    }
                    else
                    {
                      ++it;
                    }
                  }
                }
//               std::cout << "% added point #=" << i << ", "<< j << "\tat (" << cand_match.u << ", " << cand_match.v << ")" << std::endl;
              }
            }
          }
        }

//         std::cout << "% t guesses=" << uvinit.size() << std::endl;

        // cycle through all possible minima to find best
        for (it=uvinit.begin(); it!=uvinit.end(); ++it)
        {
          dd=minimum_distance(uu, vv, s, pt, it->u, it->v);
//           std::cout << "% completed root starting at (" << it->u << ", " << it->v << ") and ending at (" << uu << ", " << vv << ") with distance=" << dd << std::endl;

          // if have valid solution
          if ((uu>=umin) && (uu<=umax) && (vv>=vmin) && (vv<=vmax))
          {
            dd=eli::geom::point::distance(s.f(uu, vv), pt);

//               std::cout << "# dd=" << dd << std::endl;
            // check to see if is closer than previous minimum
            if (dd<dist)
            {
              u=uu;
              v=vv;
              dist=dd;
            }
          }
        }

//         std::cout << "# returning dist=" << dist << std::endl;
        return dist;
      }
    }
  }
}
#endif
