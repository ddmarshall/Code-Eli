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
#include <algorithm>
#include <limits>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/iterative_system_root_base_constrained.hpp"
#include "eli/mutil/nls/newton_raphson_system_method.hpp"

#include "eli/geom/intersect/minimum_distance_curve.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/intersect/minimum_distance_bounding_box.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      namespace internal
      {

        template<typename surface__, size_t N__, size_t NSOL__=1>
        class tangent_plane_method : public mutil::nls::iterative_system_root_base_constrained<typename surface__::data_type, N__, NSOL__>
        {
          public:
            static const int hit_constraint = 101;
            surface__ s;
            typename surface__::point_type pt;
            typename surface__::data_type xtol;
            typename surface__::index_type maxit;

          public:
            tangent_plane_method()
            : mutil::nls::iterative_system_root_base_constrained<typename surface__::data_type, N__, NSOL__>()
            {
              x0.setConstant(static_cast<typename surface__::data_type>(0));
            }

            tangent_plane_method(const tangent_plane_method<surface__, N__, NSOL__> &tpm)
            : mutil::nls::iterative_system_root_base_constrained<typename surface__::data_type, N__, NSOL__>(tpm), x0(tpm.x0)
            {
            }

            ~tangent_plane_method()
            {
            }

            void set_initial_guess(const typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix &xg)
            {
              x0=xg;
            }

            const typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix & get_initial_guess() const
            {
              return x0;
            }

            int find_root(typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix &root) const
            {
              typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix dx, x(x0);
              typename surface__::data_type abs_x_norm, prev_dx;
              typename surface__::index_type count;

              typename surface__::point_type q;
              typename surface__::point_type Su, Sv, r, A, B, norm;

              typename surface__::tolerance_type tol;

              typename surface__::data_type umin, umax, vmin, vmax;
              s.get_parameter_min(umin,vmin);
              s.get_parameter_max(umax,vmax);

              bool divflag = false;

              bool all_zero = false;

              abs_x_norm = std::numeric_limits<typename surface__::data_type>::max();
              prev_dx = abs_x_norm;
              count = 0;
              while ( count < maxit && abs_x_norm > xtol && !all_zero)
              {
                q = s.f(x(0), x(1));
                Su = s.f_u(x(0), x(1));
                Sv = s.f_v(x(0), x(1));

                r = q - pt;

                A = Sv.cross( r );
                B = Su.cross( r );
                norm = Su.cross( Sv );

                typename surface__::data_type N = norm.dot( norm );

                if( fabs( N ) > tol.get_absolute_tolerance() )
                {
                    dx(0) = A.dot( norm ) / N;
                    dx(1) = -B.dot( norm ) / N;
                }
                else
                {
                    dx(0) = 0.0;
                    dx(1) = 0.0;
                }

                dx = this->calculate_delta_factor(x, dx);
                x+=dx;

                prev_dx = abs_x_norm;
                abs_x_norm = this->calculate_norm(dx);

                if ( divflag && abs_x_norm > prev_dx ) // Diverging twice in a row.
                {
                  // Solution diverging, return initial guess.
                  root = x0;
                  return this->no_root_found;
                }
                else if ( abs_x_norm > prev_dx ) // Diverging first time
                {
                  divflag = true;
                }
                else // Not diverging.
                {
                  divflag = false;
                }

                all_zero = true;
                for (size_t i=0; i<N__; ++i)
                {
                  // check if stuck and cannot move x anymore
                  if ( std::abs( dx(i) ) > std::numeric_limits<typename surface__::data_type>::epsilon() )
                  {
                    all_zero = false;
                    break;
                  }
                }

                ++count;
              }

              // Current solution good enough to keep.
              root = x;

              if ( all_zero )
              {
                return this->hit_constraint;
              }

              if ( count >= maxit )
              {
                return this->max_iteration;
              }

              return this->converged;
            }

          private:
            typename mutil::nls::iterative_system_root_base<typename surface__::data_type, N__, NSOL__>::solution_matrix x0;
        };

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

            if ( !(uu>=umin) )
            {
              std::cout << "Minimum distance surface g_functor, u less than minimum.  uu: " << uu << " umin: " << umin << std::endl;
              uu=umin;
            }
            if ( !(uu<=umax) )
            {
              std::cout << "Minimum distance surface g_functor, u greater than maximum.  uu: " << uu << " uamx: " << umax << std::endl;
              uu=umax;
            }

            if ( !(vv>=vmin) )
            {
              std::cout << "Minimum distance surface g_functor, v less than minimum.  vv: " << vv << " vmin: " << vmin << std::endl;
              vv=vmin;
            }
            if ( !(vv<=vmax) )
            {
              std::cout << "Minimum distance surface g_functor, v greater than maximum.  vv: " << vv << " vmax: " << vmax << std::endl;
              vv=vmax;
            }

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

            if ( !(uu>=umin) )
            {
              std::cout << "Minimum distance surface gp_functor, u less than minimum.  uu: " << uu << " umin: " << umin << std::endl;
              uu=umin;
            }
            if ( !(uu<=umax) )
            {
              std::cout << "Minimum distance surface gp_functor, u greater than maximum.  uu: " << uu << " uamx: " << umax << std::endl;
              uu=umax;
            }

            if ( !(vv>=vmin) )
            {
              std::cout << "Minimum distance surface gp_functor, v less than minimum.  vv: " << vv << " vmin: " << vmin << std::endl;
              vv=vmin;
            }
            if ( !(vv<=vmax) )
            {
              std::cout << "Minimum distance surface gp_functor, v greater than maximum.  vv: " << vv << " vmax: " << vmax << std::endl;
              vv=vmax;
            }

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
      typename surface__::data_type minimum_distance_tan(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt,
                                                     const typename surface__::data_type &u0, const typename surface__::data_type &v0, int & ret)
      {
        typedef internal::tangent_plane_method<surface__, 2, 1> nonlinear_solver_type;

        typename surface__::data_type umin, umax, vmin, vmax;
        s.get_parameter_min(umin,vmin);
        s.get_parameter_max(umax,vmax);

        nonlinear_solver_type tpsolve;

        typename surface__::tolerance_type tol;


        if (s.open_u())
        {
          tpsolve.set_lower_condition(0, umin, nonlinear_solver_type::IRC_EXCLUSIVE);
          tpsolve.set_upper_condition(0, umax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          tpsolve.set_periodic_condition(0, umin, umax);
        }

        if (s.open_v())
        {
          tpsolve.set_lower_condition(1, vmin, nonlinear_solver_type::IRC_EXCLUSIVE);
          tpsolve.set_upper_condition(1, vmax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          tpsolve.set_periodic_condition(1, vmin, vmax);
        }

        // setup the solver
        tpsolve.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        tpsolve.set_max_iteration(20);
        tpsolve.set_norm_type(nonlinear_solver_type::max_norm);

        typename nonlinear_solver_type::solution_matrix x, x0;

        typename surface__::point_type q;
        typename surface__::data_type dist, dist0;



        assert((u0>=umin) && (u0<=umax));
        assert((v0>=vmin) && (v0<=vmax));

        x0(0) = u0;
        x0(1) = v0;
        tpsolve.set_initial_guess( x0 );
        tpsolve.pt = pt;
        tpsolve.s = s;
        tpsolve.maxit = 20;
        tpsolve.xtol = tol.get_absolute_tolerance();

        q = s.f(x0(0), x0(1));
        dist0 = eli::geom::point::distance(q, pt);

        ret = tpsolve.find_root( x );

        q = s.f(x(0), x(1));
        dist = eli::geom::point::distance(q, pt);

        if ( dist > dist0 )
        {
          x = x0; // No progress made, restore initial guess.
          dist = dist0;
        }

        u = x(0);
        v = x(1);
        return dist;
      }

      template<typename surface__>
      typename surface__::data_type minimum_distance_nrm(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt,
                                                     const typename surface__::data_type &u0, const typename surface__::data_type &v0, int & ret)
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 2, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
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
        nrm.set_absolute_f_tolerance(tol.get_absolute_tolerance());
        nrm.set_max_iteration(20);
        nrm.set_norm_type(nonlinear_solver_type::max_norm);

        if (s.open_u())
        {
          nrm.set_lower_condition(0, umin, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(0, umax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(0, umin, umax);
        }

        if (s.open_v())
        {
          nrm.set_lower_condition(1, vmin, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(1, vmax, nonlinear_solver_type::IRC_EXCLUSIVE);
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
        ret = nrm.find_root(ans, g, gp, rhs);
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

      template<typename surface__>
      typename surface__::data_type minimum_distance(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt,
                                                     const typename surface__::data_type &u0, const typename surface__::data_type &v0)
      {
        internal::tangent_plane_method<surface__, 2, 1> tan_solver;
        typename surface__::data_type dist_tan, dist_nrmt, dist_nrm0;

        int rett = -1;
        dist_tan = minimum_distance_tan( u, v, s, pt, u0, v0, rett );

        if ( rett == tan_solver.converged || rett == tan_solver.hit_constraint )
        {
          return dist_tan;
        }

        typename surface__::data_type u0t, v0t;
        u0t = u;
        v0t = v;

        int retn = -1;
        dist_nrmt = minimum_distance_nrm( u, v, s, pt, u0t, v0t, retn );

        if ( retn == tan_solver.converged )
        {
          if ( dist_nrmt <= dist_tan )
          {
            return dist_nrmt;
          }
        }

        dist_nrm0 = minimum_distance_nrm( u, v, s, pt, u0, v0, retn );

        if ( dist_nrm0 < dist_tan )
        {
            return dist_nrm0;
        }

        u = u0t;
        v = v0t;
        return dist_tan;
      }

      template<typename surface__>
      typename surface__::data_type minimum_distance(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s, const typename surface__::point_type &pt)
      {
        typename surface__::tolerance_type tol;

        // possible that end points are closest, so start by checking them
        typename surface__::data_type dist, uu, vv, dd;

        typename surface__::data_type umin, umax, vmin, vmax, uspan, vspan;
        s.get_parameter_min(umin,vmin);
        s.get_parameter_max(umax,vmax);
        uspan=umax-umin;
        vspan=vmax-vmin;

        typename surface__::index_type i, j, nu, nv;
        typename surface__::data_type du, dv;

        nu=2*s.degree_u()+1;
        nv=2*s.degree_v()+1;

        // Evenly spaced in parameter, don't repeat 0/1 if closed curve.
        if (s.open_u())
        {
          du = uspan/(nu-1);
        }
        else
        {
          du = uspan/nu;
        }

        if (s.open_v())
        {
          dv = vspan/(nv-1);
        }
        else
        {
          dv = vspan/nv;
        }

        // Find closest of evenly spaced points.
        uu=umin;
        dist = std::numeric_limits<typename surface__::data_type>::max();
        for (i = 0; i < nu; i++)
        {
          vv=vmin;
          for (j = 0; j < nv; j++)
          {
            dd=eli::geom::point::distance(s.f(uu,vv), pt);

            if( dd < dist )
            {
              u=uu;
              v=vv;
              dist=dd;
            }
            vv+=dv;
            if(vv>=vmax)
            {
              vv=vmax;
            }
          }
          uu+=du;
          if(uu>=umax)
          {
            uu=umax;
          }
        }

        // Polish best point with Newton's method search.
        dd=minimum_distance(uu, vv, s, pt, u, v);

        if ((uu>=umin) && (uu<=umax) && (vv>=vmin) && (vv<=vmax))
        {
          if (dd<dist)
          {
            u=uu;
            v=vv;
            dist=dd;
          }
        }

        // next check edges
        // Since these are always edges, we could implement an edge curve extraction routine
        // that returned the control points directly instead of performing an arbitrary curve
        // extraction calculation.
        typename surface__::curve_type bc;
        if(u<=(umin+std::abs(umin)*2*std::numeric_limits<typename surface__::data_type>::epsilon()))
        {
          s.get_uconst_curve(bc, umin);
          dd=eli::geom::intersect::minimum_distance(vv, bc, pt, v);

          if (dd<dist)
          {
            u=umin;
            v=vv;
            dist=dd;
          }
        }

        if(u>=(umax-std::abs(umax)*2*std::numeric_limits<typename surface__::data_type>::epsilon()))
        {
          s.get_uconst_curve(bc, umax);
          dd=eli::geom::intersect::minimum_distance(vv, bc, pt, v);

          if (dd<dist)
          {
            u=umax;
            v=vv;
            dist=dd;
          }
        }

        if(v<=(vmin+std::abs(vmin)*2*std::numeric_limits<typename surface__::data_type>::epsilon()))
        {
          s.get_vconst_curve(bc, vmin);
          dd=eli::geom::intersect::minimum_distance(uu, bc, pt, u);

          if (dd<dist)
          {
            u=uu;
            v=vmin;
            dist=dd;
          }
        }

        if(v>=(vmax-std::abs(vmax)*2*std::numeric_limits<typename surface__::data_type>::epsilon()))
        {
          s.get_vconst_curve(bc, vmax);
          dd=eli::geom::intersect::minimum_distance(uu, bc, pt, u);

          if (dd<dist)
          {
            u=uu;
            v=vmax;
            dist=dd;
          }
        }

        return dist;

      }

// Defined for minimum_distance_curve.  Could be moved to util somewhere.
//      template< typename first__, typename second__>
//      bool pairfirstcompare( const std::pair < first__, second__ > &a, const std::pair < first__, second__ > &b )
//      {
//          return ( a.first < b.first );
//      }

      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__ >
      typename surface::piecewise<surface__, data__, dim__, tol__>::data_type minimum_distance(
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &u,
          typename surface::piecewise<surface__, data__, dim__, tol__>::data_type &v,
          const surface::piecewise<surface__, data__, dim__, tol__> &ps,
          const typename surface::piecewise<surface__, data__, dim__, tol__>::point_type &pt)
      {
        typedef surface::piecewise<surface__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::surface_type surface_type;
        typedef typename piecewise_type::index_type index_type;
        typedef typename piecewise_type::data_type data_type;
        typedef typename piecewise_type::bounding_box_type bounding_box_type;

        typedef typename piecewise_type::keymap_type keymap_type;
        typedef typename keymap_type::const_iterator keyit;

        typedef std::pair<keyit, keyit> itpair;
        typedef std::vector< std::pair<data_type, itpair > > dvec;
        dvec minbbdist;

        // Find closest corner of bounding boxes, add them to vector
        // Simple linear search, would be more efficient with some sort of tree.
        for(keyit uit = ps.ukey.key.begin(); uit != ps.ukey.key.end(); ++uit)
        {
          for(keyit vit = ps.vkey.key.begin(); vit != ps.vkey.key.end(); ++vit)
          {
            index_type uk = uit->second;
            index_type vk = vit->second;

            surface_type s = ps.patches[uk][vk];

            bounding_box_type bb_local;
            s.get_bounding_box(bb_local);

            data_type dbbmin;
            dbbmin = minimum_distance(bb_local, pt);

            minbbdist.push_back(std::make_pair(dbbmin, std::make_pair(uit, vit)));

          }
        }

        // Sort by nearest distance.
        std::sort( minbbdist.begin(), minbbdist.end(), pairfirstcompare<data_type, itpair > );


        // Iterate over segments, starting with nearest bounding box
        data_type dist(std::numeric_limits<data_type>::max());

        typename dvec::const_iterator it;
        for (it=minbbdist.begin(); it!=minbbdist.end(); ++it)
        {
          // If nearest bb distance is farther than current best, we're done.
          if(it->first < dist )
          {
            itpair itp = it->second;
            keyit uit = itp.first;
            keyit vit = itp.second;

            index_type uk = uit->second;
            index_type vk = vit->second;

            surface_type s = ps.patches[uk][vk];

            data_type uu, vv, d;
            d=minimum_distance(uu, vv, s, pt);

            if(d < dist)
            {
              data_type du(ps.ukey.get_delta_parm(uit));
              data_type dv(ps.vkey.get_delta_parm(vit));

              data_type ustart(uit->first);
              data_type vstart(vit->first);

              dist = d;
              u=ustart+uu*du;
              v=vstart+vv*dv;
            }
          }
          else
          {
            break;
          }

        }
        return dist;
      }

    }
  }
}
#endif
