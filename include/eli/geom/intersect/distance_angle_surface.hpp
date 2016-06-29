/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald
********************************************************************************/

#ifndef eli_geom_intersect_distance_angle_surface_hpp
#define eli_geom_intersect_distance_angle_surface_hpp

#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <limits>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/iterative_system_root_base_constrained.hpp"
#include "eli/mutil/nls/newton_raphson_system_method.hpp"

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
        struct surface_da_functor
        {
          const surface__ *ps;
          typename surface__::point_type pt;
          typename surface__::point_type dir;
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
              std::cout << "Distance angle surface g_functor, u less than minimum.  uu: " << uu << " umin: " << umin << std::endl;
              uu=umin;
            }
            if ( !(uu<=umax) )
            {
              std::cout << "Distance angle surface g_functor, u greater than maximum.  uu: " << uu << " uamx: " << umax << std::endl;
              uu=umax;
            }

            if ( !(vv>=vmin) )
            {
              std::cout << "Distance angle surface g_functor, v less than minimum.  vv: " << vv << " vmin: " << vmin << std::endl;
              vv=vmin;
            }
            if ( !(vv<=vmax) )
            {
              std::cout << "Distance angle surface g_functor, v greater than maximum.  vv: " << vv << " vmax: " << vmax << std::endl;
              vv=vmax;
            }

            assert((uu>=umin) && (uu<=umax));
            assert((vv>=vmin) && (vv<=vmax));

            uu=std::min(std::max(uu, static_cast<typename surface__::data_type>(umin)), static_cast<typename surface__::data_type>(umax));
            vv=std::min(std::max(vv, static_cast<typename surface__::data_type>(vmin)), static_cast<typename surface__::data_type>(vmax));

            typename surface__::point_type tmp;

            tmp=ps->f(uu, vv)-pt;
            rtn(0) = tmp.dot( tmp );
            rtn(1) = tmp.dot( dir ) / sqrt( rtn(0) );
            return rtn;
          }
        };

        template <typename surface__>
        struct surface_dap_functor
        {
          const surface__ *ps;
          typename surface__::point_type pt;
          typename surface__::point_type dir;
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
              std::cout << "Distance angle surface gp_functor, u less than minimum.  uu: " << uu << " umin: " << umin << std::endl;
              uu=umin;
            }
            if ( !(uu<=umax) )
            {
              std::cout << "Distance angle surface gp_functor, u greater than maximum.  uu: " << uu << " uamx: " << umax << std::endl;
              uu=umax;
            }

            if ( !(vv>=vmin) )
            {
              std::cout << "Distance angle surface gp_functor, v less than minimum.  vv: " << vv << " vmin: " << vmin << std::endl;
              vv=vmin;
            }
            if ( !(vv<=vmax) )
            {
              std::cout << "Distance angle surface gp_functor, v greater than maximum.  vv: " << vv << " vmax: " << vmax << std::endl;
              vv=vmax;
            }

            assert((uu>=umin) && (uu<=umax));
            assert((vv>=vmin) && (vv<=vmax));

            uu=std::min(std::max(uu, static_cast<typename surface__::data_type>(umin)), static_cast<typename surface__::data_type>(umax));
            vv=std::min(std::max(vv, static_cast<typename surface__::data_type>(vmin)), static_cast<typename surface__::data_type>(vmax));

            typename surface__::point_type tmp, Su, Sv;

            tmp=ps->f(uu, vv)-pt;
            Su=ps->f_u(uu, vv);
            Sv=ps->f_v(uu, vv);

            typename surface__::data_type k, krt;
            k = tmp.dot( tmp );
            krt = sqrt( k );

            typename surface__::data_type tmp2;
            tmp2 = tmp.dot(dir) / ( k * krt );

            rtn(0,0) = 2.0 * Su.dot(tmp);
            rtn(0,1) = 2.0 * Sv.dot(tmp);
            rtn(1,0) = Su.dot(dir) / krt - tmp2 * Su.dot(tmp);
            rtn(1,1) = Sv.dot(dir) / krt - tmp2 * Sv.dot(tmp);

            // TODO: What to do if matrix becomes singular?

            return rtn;
          }
        };
      }

      template<typename surface__>
      void distance_angle(typename surface__::data_type &u, typename surface__::data_type &v, const surface__ &s,
                                                   const typename surface__::point_type &pt, const typename surface__::point_type &dir,
												   const typename surface__::data_type &dsq, const typename surface__::data_type &dot,
                                                   const typename surface__::data_type &u0, const typename surface__::data_type &v0,
                                                   const typename surface__::data_type &vllim, const typename surface__::data_type &vulim, int & ret)
      {
        typedef eli::mutil::nls::newton_raphson_system_method<typename surface__::data_type, 2, 1> nonlinear_solver_type;
        nonlinear_solver_type nrm;
        internal::surface_da_functor<surface__> da;
        internal::surface_dap_functor<surface__> dap;
        typename surface__::data_type dist0, dist;
        typename surface__::tolerance_type tol;

        typename surface__::data_type umin, umax, vmin, vmax;
        s.get_parameter_min(umin,vmin);
        s.get_parameter_max(umax,vmax);

        // setup the functors
        da.ps=&s;
        da.pt=pt;
        da.dir=dir;

        dap.ps=&s;
        dap.pt=pt;
        dap.dir=dir;

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

        nrm.set_lower_condition(1, vllim, nonlinear_solver_type::IRC_EXCLUSIVE);
        nrm.set_upper_condition(1, vulim, nonlinear_solver_type::IRC_EXCLUSIVE);

        /*
        if (s.open_v())
        {
          nrm.set_lower_condition(1, vmin, nonlinear_solver_type::IRC_EXCLUSIVE);
          nrm.set_upper_condition(1, vmax, nonlinear_solver_type::IRC_EXCLUSIVE);
        }
        else
        {
          nrm.set_periodic_condition(1, vmin, vmax);
        }
        */

        // set the initial guess
        typename nonlinear_solver_type::solution_matrix uinit, rhs, ans, rhs0, rhs1;

        uinit(0)=u0;
        uinit(1)=v0;
        nrm.set_initial_guess(uinit);
        rhs(0)=dsq;
        rhs(1)=dot;

        rhs0 = da( uinit );

        // find the root
        ret = nrm.find_root(ans, da, dap, rhs);
        u=ans(0);
        v=ans(1);

        // if root is within bounds and is closer than initial guess
        {
          assert((u>=umin) && (u<=umax));
          assert((v>=vllim) && (v<=vulim));

          rhs1 = da( ans );

          if ( (rhs1-rhs).norm() <= (rhs0-rhs).norm() )
          {
              return;
          }
        }

        // couldn't find better answer so return initial guess
        u=u0;
        v=v0;
        return;
      }

    }
  }
}
#endif
