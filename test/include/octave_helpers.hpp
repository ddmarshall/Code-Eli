/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_test_octave_helpers_hpp
#define eli_geom_test_octave_helpers_hpp

#include <string>

#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/surface/bezier.hpp"
#include "eli/geom/surface/piecewise.hpp"

namespace eli
{
  namespace test
  {
    inline std::string random_string( size_t length )
    {
        auto randchar = []() -> char
        {
            const char charset[] =
            "0123456789"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            "abcdefghijklmnopqrstuvwxyz";
            const size_t max_index = (sizeof(charset) - 1);
            return charset[ rand() % max_index ];
        };
        std::string str(length,0);
        std::generate_n( str.begin(), length, randchar );
        return str;
    }

    inline void octave_start(int figno)
    {
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << "clf(" << figno << ", 'reset');" << std::endl;
      std::cout << "hold on;" << std::endl;
      std::cout << "grid on;" << std::endl;
    }

    inline void octave_finish(int figno)
    {
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << "hold off;" << std::endl;
      std::cout << "rotate3d on;" << std::endl;
      std::cout << "xlabel('x');" << std::endl;
      std::cout << "ylabel('y');" << std::endl;
      std::cout << "zlabel('z');" << std::endl;
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &pc,
                      const std::string &name="", bool show_control_points=true)
    {
      typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
      typedef typename piecewise_curve_type::curve_type curve_type;
      typedef typename piecewise_curve_type::data_type data_type;
      typedef typename piecewise_curve_type::index_type index_type;

      std::string nm, cpxbuf, cpybuf, cpzbuf, cxbuf, cybuf, czbuf;

      index_type i, pp, ns;
      data_type tmin(pc.get_parameter_min()), tmax(pc.get_parameter_max());

      ns=pc.number_segments();

      // build name
      if (name=="")
      {
        nm=random_string(5);
      }
      else
      {
        nm=name;
      }

      // set control points
      if (show_control_points)
      {
        cpxbuf=nm+"_curv_cp_x=[";
        cpybuf=nm+"_curv_cp_y=[";
        cpzbuf=nm+"_curv_cp_z=[";
        for (pp=0; pp<ns; ++pp)
        {
          index_type bez_deg;
          curve_type bez;
          pc.get(bez, pp);
          bez_deg=bez.degree();
          for (i=0; i<=bez_deg; ++i)
          {
            cpxbuf+=std::to_string(bez.get_control_point(i).x());
            cpybuf+=std::to_string(bez.get_control_point(i).y());
            cpzbuf+=std::to_string(bez.get_control_point(i).z());
            if ((pp<(ns-1)) || ((pp==(ns-1)) && (i<bez_deg)))
            {
              cpxbuf+=", ";
              cpybuf+=", ";
              cpzbuf+=", ";
            }
          }
        }
        cpxbuf+="];";
        cpybuf+="];";
        cpzbuf+="];";
      }

      // initialize the t parameters
      index_type nt(129);
      std::vector<data__> t(nt);
      for (i=0; i<nt; ++i)
      {
        t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(nt-1);
      }

      // set the curve points
      cxbuf=nm+"_curv_x=[";
      cybuf=nm+"_curv_y=[";
      czbuf=nm+"_curv_z=[";
      for (i=0; i<nt; ++i)
      {
        cxbuf+=std::to_string(pc.f(t[i]).x());
        cybuf+=std::to_string(pc.f(t[i]).y());
        czbuf+=std::to_string(pc.f(t[i]).z());
        if (i<nt)
        {
          cxbuf+=", ";
          cybuf+=", ";
          czbuf+=", ";
        }
      }
      cxbuf+="];";
      cybuf+="];";
      czbuf+="];";

      std::cout << "% curve: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << cxbuf << std::endl;
      std::cout << cybuf << std::endl;
      std::cout << czbuf << std::endl;
      if (show_control_points)
      {
        std::cout << cpxbuf << std::endl;
        std::cout << cpybuf << std::endl;
        std::cout << cpzbuf << std::endl;
      }
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot3(" << nm << "_curv_x, "
                            << nm << "_curv_y, "
                            << nm << "_curv_z, '-g');" << std::endl;
      if (show_control_points)
      {
        std::cout << "plot3(" << nm << "_curv_cp_x', "
                              << nm << "_curv_cp_y', "
                              << nm << "_curv_cp_z', '-o', 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0]);" << std::endl;
      }
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &pc,
                      const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &vec,
                      const std::string &name="")
    {
      typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
      typedef typename piecewise_curve_type::data_type data_type;
      typedef typename piecewise_curve_type::index_type index_type;
      typedef typename piecewise_curve_type::tolerance_type tolerance_type;

      std::string nm(random_string(5)), vecxbuf, vecybuf, veczbuf, cxbuf, cybuf, czbuf;

      index_type i, pp, ns;
      data_type tmin, tmax;

      tolerance_type tol;

      ns=pc.number_segments();
      pc.get_parameter_min(tmin);
      pc.get_parameter_max(tmax);

      // check parameterization of vec curve
      if (!tol.approximately_equal(vec.get_t0(), tmin))
      {
        return;
      }
      if (!tol.approximately_equal(vec.get_tmax(), tmax))
      {
        return;
      }

      // build name
      if (name=="")
      {
        nm=random_string(5);
      }
      else
      {
        nm=name;
      }

      // initialize the t parameters
      index_type nt(11);
      std::vector<data__> t(nt);
      for (i=0; i<nt; ++i)
      {
        t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(nt-1);
      }

      // set the surface points
      cxbuf=nm+"_loc_x=[";
      cybuf=nm+"_loc_y=[";
      czbuf=nm+"_loc_z=[";
      vecxbuf=nm+"_vec_x=[";
      vecybuf=nm+"_vec_y=[";
      veczbuf=nm+"_vec_z=[";
      for (i=0; i<nt; ++i)
      {
        cxbuf+=std::to_string(pc.f(t[i]).x());
        cybuf+=std::to_string(pc.f(t[i]).y());
        czbuf+=std::to_string(pc.f(t[i]).z());
        vecxbuf+=std::to_string(vec.f(t[i]).x());
        vecybuf+=std::to_string(vec.f(t[i]).y());
        veczbuf+=std::to_string(vec.f(t[i]).z());
        if (i<nt)
        {
          cxbuf+=", ";
          cybuf+=", ";
          czbuf+=", ";
          vecxbuf+=", ";
          vecybuf+=", ";
          veczbuf+=", ";
        }
      }
      cxbuf+="];";
      cybuf+="];";
      czbuf+="];";
      vecxbuf+="];";
      vecybuf+="];";
      veczbuf+="];";

      std::cout << "% curve: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << cxbuf << std::endl;
      std::cout << cybuf << std::endl;
      std::cout << czbuf << std::endl;
      std::cout << vecxbuf << std::endl;
      std::cout << vecybuf << std::endl;
      std::cout << veczbuf << std::endl;
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "quiver3(" << nm << "_loc_x, "
                              << nm << "_loc_y, "
                              << nm << "_loc_z, "
                              << nm << "_vec_x, "
                              << nm << "_vec_y, "
                              << nm << "_vec_z, "
                              << "'r');" << std::endl;
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> &ps,
                      const std::string &name="", bool show_control_points=true)
    {
      typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
      typedef typename piecewise_surface_type::surface_type surface_type;
      typedef typename piecewise_surface_type::data_type data_type;
      typedef typename piecewise_surface_type::index_type index_type;

      std::string nm, sxbuf, sybuf, szbuf;
      std::vector<std::string> cpxbuf, cpybuf, cpzbuf;

      index_type i, j, pp, qq, nup, nvp;
      data_type umin, vmin, umax, vmax;

      nup=ps.number_u_patches();
      nvp=ps.number_v_patches();
      ps.get_parameter_min(umin, vmin);
      ps.get_parameter_max(umax, vmax);


      // build name
      if (name=="")
      {
        nm=random_string(5);
      }
      else
      {
        nm=name;
      }

      // set control points
      if (show_control_points)
      {
        cpxbuf.resize(nup*nvp);
        cpybuf.resize(nup*nvp);
        cpzbuf.resize(nup*nvp);

        for (pp=0; pp<nup; ++pp)
        {
          for (qq=0; qq<nvp; ++qq)
          {
            surface_type bez;
            ps.get(bez, pp, qq);
            index_type ppqq;

            ppqq=qq*nup+pp;
            cpxbuf[ppqq]=nm+"_surf_cp"+std::to_string(pp)+std::to_string(qq)+"_x=[";
            cpybuf[ppqq]=nm+"_surf_cp"+std::to_string(pp)+std::to_string(qq)+"_y=[";
            cpzbuf[ppqq]=nm+"_surf_cp"+std::to_string(pp)+std::to_string(qq)+"_z=[";

            for (i=0; i<=bez.degree_u(); ++i)
            {
              for (j=0; j<=bez.degree_v(); ++j)
              {
                cpxbuf[ppqq]+=std::to_string(bez.get_control_point(i, j).x());
                cpybuf[ppqq]+=std::to_string(bez.get_control_point(i, j).y());
                cpzbuf[ppqq]+=std::to_string(bez.get_control_point(i, j).z());

                if ((i==bez.degree_u()) && (j==bez.degree_v()))
                {
                }
                else if ((j==bez.degree_v()) && (i<bez.degree_u()))
                {
                  cpxbuf[ppqq]+="; ";
                  cpybuf[ppqq]+="; ";
                  cpzbuf[ppqq]+="; ";
                }
                else
                {
                  cpxbuf[ppqq]+=", ";
                  cpybuf[ppqq]+=", ";
                  cpzbuf[ppqq]+=", ";
                }
              }
            }
            cpxbuf[ppqq]+="];";
            cpybuf[ppqq]+="];";
            cpzbuf[ppqq]+="];";
          }
        }
      }

      // initialize the u & v parameters
      index_type nu(10*nup+1), nv(10*nvp+1);
      std::vector<data__> u(nu), v(nv);
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        u[i]=umin+(umax-umin)*static_cast<data__>(i)/(u.size()-1);
      }
      for (j=0; j<static_cast<index_type>(v.size()); ++j)
      {
        v[j]=vmin+(vmax-vmin)*static_cast<data__>(j)/(v.size()-1);
      }

      // set the surface points
      sxbuf=nm+"_surf_x=[";
      sybuf=nm+"_surf_y=[";
      szbuf=nm+"_surf_z=[";
      for (i=0; i<nu; ++i)
      {
        for (j=0; j<nv; ++j)
        {
          sxbuf+=std::to_string(ps.f(u[i], v[j]).x());
          sybuf+=std::to_string(ps.f(u[i], v[j]).y());
          szbuf+=std::to_string(ps.f(u[i], v[j]).z());
          if (j==(nv-1))
          {
            if (i<(nu-1))
            {
              sxbuf+=";\n";
              sybuf+=";\n";
              szbuf+=";\n";
            }
          }
          else
          {
            sxbuf+=", ";
            sybuf+=", ";
            szbuf+=", ";
          }
        }
      }
      sxbuf+="];";
      sybuf+="];";
      szbuf+="];";

      std::cout << "% surface: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << sxbuf << std::endl;
      std::cout << sybuf << std::endl;
      std::cout << szbuf << std::endl;
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "mesh(" << nm << "_surf_x, "
                           << nm << "_surf_y, "
                           << nm << "_surf_z, "
                           << "'EdgeColor', [0 0 0]);" << std::endl;
      if (show_control_points)
      {
        for (pp=0; pp<nup; ++pp)
        {
          for (qq=0; qq<nvp; ++qq)
          {
            index_type ppqq;

            ppqq=qq*nup+pp;
            std::cout << cpxbuf[ppqq] << std::endl;
            std::cout << cpybuf[ppqq] << std::endl;
            std::cout << cpzbuf[ppqq] << std::endl;

            std::cout << "mesh(" << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_x', "
                                 << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_y', "
                                 << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_z', "
                                 << "'EdgeColor', [0.5 0.5 0.5], 'FaceColor', 'none');" << std::endl;
            std::cout << "plot3(" << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_x', "
                                  << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_y', "
                                  << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_z', "
                                  << "'o', 'MarkerEdgeColor', [0.5 0.5 0.5],'MarkerFaceColor', [0.5 0.5 0.5]);" << std::endl;
          }
        }
      }
    }
  }
}

#endif

