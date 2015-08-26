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
#include "eli/geom/curve/pseudo/explicit_bezier.hpp"
#include "eli/geom/curve/pseudo/polynomial.hpp"
#include "eli/geom/curve/pseudo/four_digit.hpp"
#include "eli/geom/curve/pseudo/cst_airfoil.hpp"
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

    inline void octave_finish(int figno, bool axis_equal)
    {
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << "hold off;" << std::endl;
      std::cout << "rotate3d on;" << std::endl;
      std::cout << "xlabel('x');" << std::endl;
      std::cout << "ylabel('y');" << std::endl;
      std::cout << "zlabel('z');" << std::endl;
      if (axis_equal)
      {
        std::cout << "axis equal;" << std::endl;
      }
    }

    template<typename data__, typename index__>
    void calculate_t(std::vector<data__> &t, std::vector<index__> &npseg, const std::vector<data__> &tseg, const data__ &dt)
    {
      typedef index__ index_type;
      typedef data__ data_type;

      index_type i, j, nt, nt_cumulative(0), nsegs(static_cast<index_type> (tseg.size())-1);
      nt=(tseg[nsegs]-tseg[0])/dt;

      // calulate the actual number of points to allocate taking into
      // account that small segments will get at least 2 points
      for (i=0; i<nsegs; ++i)
      {
        npseg[i]=std::max(static_cast<data_type>(2), std::ceil((tseg[i+1]-tseg[i])/dt));
        nt_cumulative+=npseg[i];
      }
      nt=nt_cumulative;

      // initialize the t parameters
      index_type nrun(0);
      data_type small_num(std::sqrt(std::numeric_limits<data_type>::epsilon()));

      t.resize(nt);
      for (i=0; i<nsegs; ++i)
      {
        if (i>0)
        {
          t[nrun]=tseg[i]*(1+small_num);
        }
        else
        {
          t[nrun]=tseg[i];
        }
        ++nrun;
        for (j=1; j<(npseg[i]-1); ++j, ++nrun)
        {
          t[nrun]=tseg[i]+(tseg[i+1]-tseg[i])*static_cast<data_type>(j)/(npseg[i]-1);
        }
        if (i<(nsegs-1))
        {
          t[nrun]=tseg[i+1]*(1-small_num);
        }
        else
        {
          t[nrun]=tseg[i+1];
        }
        ++nrun;
      }
    }

    template<typename point__>
    void octave_print(int figno, const point__ &pt, const std::string &name="")
    {
      std::string nm, ptxbuf, ptybuf, ptzbuf;
      bool d3(pt.cols()==3);

      // hack to ensure that we have either 2d or 3d points
      if ((!d3) && (pt.cols()!=2))
      {
        assert(false);
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

      // build point coordinates
      ptxbuf=nm+"_pt_x=[";
      ptybuf=nm+"_pt_y=[";
      if (d3)
      {
        ptzbuf=nm+"_pt_z=[";
      }

      ptxbuf+=std::to_string(pt.x());
      ptybuf+=std::to_string(pt.y());
      if (d3)
      {
        ptzbuf+=std::to_string(pt.z());
      }

      ptxbuf+="];";
      ptybuf+="];";
      if (d3)
      {
        ptzbuf+="];";
      }

      std::cout << "% point: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << ptxbuf << std::endl;
      std::cout << ptybuf << std::endl;
      if (d3)
      {
        std::cout << ptzbuf << std::endl;
      }
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      if (d3)
      {
        std::cout << "plot3(" << nm << "_pt_x, "
                              << nm << "_pt_y, "
                              << nm << "_pt_z, "
                              << "'d', 'Color', [0 0 0.75], 'MarkerFaceColor', [0 0 0.75]);" << std::endl;
      }
      else
      {
        std::cout << "plot(" << nm << "_pt_x, "
                             << nm << "_pt_y, "
                             << "'d', 'Color', [0 0 0.75], 'MarkerFaceColor', [0 0 0.75]);" << std::endl;
      }
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::bezier<data__, 2> &bc,
                      const std::string &name="", bool show_control_points=true)
    {
      typedef eli::geom::curve::bezier<data__, 2> curve_type;
      typedef typename curve_type::data_type data_type;
      typedef typename curve_type::point_type point_type;
      typedef typename curve_type::index_type index_type;

      std::string nm, cpxbuf, cpybuf, cxbuf, cybuf;

      index_type i;
      data_type tmin(0), tmax(1);

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
        {
          index_type bez_deg;
          bez_deg=bc.degree();
          for (i=0; i<=bez_deg; ++i)
          {
            point_type pt(bc.get_control_point(i));

            cpxbuf+=std::to_string(pt.x());
            cpybuf+=std::to_string(pt.y());
            if (i<bez_deg)
            {
              cpxbuf+=", ";
              cpybuf+=", ";
            }
          }
        }
        cpxbuf+="];";
        cpybuf+="];";
      }

      // initialize the t parameters
      index_type nt(201);
      std::vector<data__> t(nt);
      for (i=0; i<nt; ++i)
      {
        t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(nt-1);
      }

      // set the curve points
      cxbuf=nm+"_curv_x=[";
      cybuf=nm+"_curv_y=[";
      for (i=0; i<nt; ++i)
      {
        point_type pt(bc.f(t[i]));

        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        if (i<(nt-1))
        {
          cxbuf+=", ";
          cybuf+=", ";
        }
      }
      cxbuf+="];";
      cybuf+="];";

      std::cout << "% curve: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << cxbuf << std::endl;
      std::cout << cybuf << std::endl;
      if (show_control_points)
      {
        std::cout << cpxbuf << std::endl;
        std::cout << cpybuf << std::endl;
      }
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot(" << nm << "_curv_x, "
                           << nm << "_curv_y, "
                           << "'-g');" << std::endl;
      if (show_control_points)
      {
        std::cout << "plot(" << nm << "_curv_cp_x', "
                             << nm << "_curv_cp_y', "
                             << "'-o', 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0]);" << std::endl;
      }
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::bezier<data__, 3> &bc,
                      const std::string &name="", bool show_control_points=true)
    {
      typedef eli::geom::curve::bezier<data__, 3> curve_type;
      typedef typename curve_type::data_type data_type;
      typedef typename curve_type::point_type point_type;
      typedef typename curve_type::index_type index_type;

      std::string nm, cpxbuf, cpybuf, cpzbuf, cxbuf, cybuf, czbuf;

      index_type i;
      data_type tmin(0), tmax(1);

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
        {
          index_type bez_deg;
          bez_deg=bc.degree();
          for (i=0; i<=bez_deg; ++i)
          {
            point_type pt(bc.get_control_point(i));

            cpxbuf+=std::to_string(pt.x());
            cpybuf+=std::to_string(pt.y());
            cpzbuf+=std::to_string(pt.z());
            if (i<bez_deg)
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
      index_type nt(201);
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
        point_type pt(bc.f(t[i]));

        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        czbuf+=std::to_string(pt.z());
        if (i<(nt-1))
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
                            << nm << "_curv_z, "
                            << "'-g');" << std::endl;
      if (show_control_points)
      {
        std::cout << "plot3(" << nm << "_curv_cp_x', "
                              << nm << "_curv_cp_y', "
                              << nm << "_curv_cp_z', "
                              << "'-o', 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0]);" << std::endl;
      }
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::pseudo::explicit_bezier<data__> &ebc,
                      const std::string &name="")
    {
      typedef eli::geom::curve::pseudo::explicit_bezier<data__> explicit_bezier_curve_type;
      typedef typename explicit_bezier_curve_type::data_type data_type;
      typedef typename explicit_bezier_curve_type::point_type point_type;
      typedef typename explicit_bezier_curve_type::index_type index_type;

      std::string nm, cxbuf, cybuf;

      index_type i, pp, ns;
      data_type tmin(ebc.get_t0()), tmax(ebc.get_tmax());

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
      index_type nt(201);
      std::vector<data__> t(nt);
      for (i=0; i<nt; ++i)
      {
        t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(nt-1);
      }

      // set the curve points
      cxbuf=nm+"_curv_x=[";
      cybuf=nm+"_curv_y=[";
      for (i=0; i<nt; ++i)
      {
        point_type pt(ebc.f(t[i]));

        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        if (i<(nt-1))
        {
          cxbuf+=", ";
          cybuf+=", ";
        }
      }
      cxbuf+="];";
      cybuf+="];";

      std::cout << "% curve: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << cxbuf << std::endl;
      std::cout << cybuf << std::endl;
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot(" << nm << "_curv_x, "
                           << nm << "_curv_y, '-c');" << std::endl;
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::pseudo::polynomial<data__, 2> &pol,
                      const std::string &name="")
    {
      typedef eli::geom::curve::pseudo::polynomial<data__, 2> polynomial_curve_type;
      typedef typename polynomial_curve_type::data_type data_type;
      typedef typename polynomial_curve_type::point_type point_type;
      typedef typename polynomial_curve_type::index_type index_type;

      std::string nm, cxbuf, cybuf;

      index_type i, pp, ns;
      data_type tmin(0), tmax(1);

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
      index_type nt(201);
      std::vector<data__> t(nt);
      for (i=0; i<nt; ++i)
      {
        t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(nt-1);
      }

      // set the curve points
      cxbuf=nm+"_curv_x=[";
      cybuf=nm+"_curv_y=[";
      for (i=0; i<nt; ++i)
      {
        point_type pt(pol.f(t[i]));

        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        if (i<(nt-1))
        {
          cxbuf+=", ";
          cybuf+=", ";
        }
      }
      cxbuf+="];";
      cybuf+="];";

      std::cout << "% curve: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << cxbuf << std::endl;
      std::cout << cybuf << std::endl;
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot(" << nm << "_curv_x, "
                           << nm << "_curv_y, '-c');" << std::endl;
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::pseudo::polynomial<data__, 3> &pol,
                      const std::string &name="")
    {
      typedef eli::geom::curve::pseudo::polynomial<data__, 3> polynomial_curve_type;
      typedef typename polynomial_curve_type::data_type data_type;
      typedef typename polynomial_curve_type::point_type point_type;
      typedef typename polynomial_curve_type::index_type index_type;

      std::string nm, cxbuf, cybuf, czbuf;

      index_type i, pp, ns;
      data_type tmin(0), tmax(1);

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
      index_type nt(201);
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
        point_type pt(pol.f(t[i]));

        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        czbuf+=std::to_string(pt.z());
        if (i<(nt-1))
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
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot3(" << nm << "_curv_x, "
                            << nm << "_curv_y, "
                            << nm << "_curv_z, '-c');" << std::endl;
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::pseudo::four_digit<data__> &af,
                      const std::string &name="")
    {
      typedef eli::geom::curve::pseudo::four_digit<data__> four_digit_airfoil_curve_type;
      typedef typename four_digit_airfoil_curve_type::data_type data_type;
      typedef typename four_digit_airfoil_curve_type::point_type point_type;
      typedef typename four_digit_airfoil_curve_type::index_type index_type;

      std::string nm, cxbuf, cybuf;

      index_type i, pp, ns;
      data_type tmin(af.get_t0()), tmax(af.get_tmax());

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
      index_type nt(201);
      std::vector<data__> t(nt);
      for (i=0; i<nt; ++i)
      {
        t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(nt-1);
      }

      // set the curve points
      cxbuf=nm+"_curv_x=[";
      cybuf=nm+"_curv_y=[";
      for (i=0; i<nt; ++i)
      {
        point_type pt(af.f(t[i]));
        
        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        if (i<(nt-1))
        {
          cxbuf+=", ";
          cybuf+=", ";
        }
      }
      cxbuf+="];";
      cybuf+="];";

      std::cout << "% curve: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << cxbuf << std::endl;
      std::cout << cybuf << std::endl;
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot(" << nm << "_curv_x, "
                           << nm << "_curv_y, '-c');" << std::endl;
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::pseudo::cst_airfoil<data__> &cst,
                      const std::string &name="")
    {
      typedef eli::geom::curve::pseudo::cst_airfoil<data__> cst_airfoil_curve_type;
      typedef typename cst_airfoil_curve_type::data_type data_type;
      typedef typename cst_airfoil_curve_type::point_type point_type;
      typedef typename cst_airfoil_curve_type::index_type index_type;

      std::string nm, cxbuf, cybuf;

      index_type i, pp, ns;
      data_type tmin(cst.get_t0()), tmax(cst.get_tmax());

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
      index_type nt(201);
      std::vector<data__> t(nt);
      for (i=0; i<nt; ++i)
      {
        t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(nt-1);
      }

      // set the curve points
      cxbuf=nm+"_curv_x=[";
      cybuf=nm+"_curv_y=[";
      for (i=0; i<nt; ++i)
      {
        point_type pt(cst.f(t[i]));
        
        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        if (i<(nt-1))
        {
          cxbuf+=", ";
          cybuf+=", ";
        }
      }
      cxbuf+="];";
      cybuf+="];";

      std::cout << "% curve: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << cxbuf << std::endl;
      std::cout << cybuf << std::endl;
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot(" << nm << "_curv_x, "
                           << nm << "_curv_y, '-c');" << std::endl;
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 2> &pc,
                      const std::string &name="", bool show_control_points=true, bool show_segment_end_points=true)
    {
      typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 2> piecewise_curve_type;
      typedef typename piecewise_curve_type::curve_type curve_type;
      typedef typename piecewise_curve_type::data_type data_type;
      typedef typename piecewise_curve_type::point_type point_type;
      typedef typename piecewise_curve_type::index_type index_type;

      std::string nm, cpxbuf, cpybuf, cxbuf, cybuf, jointxbuf, jointybuf;

      index_type i, pp, ns;

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
      }
      if (show_segment_end_points)
      {
        jointxbuf=nm+"_curv_joint_x=[";
        jointybuf=nm+"_curv_joint_y=[";
      }
      for (pp=0; pp<ns; ++pp)
      {
        index_type bez_deg;
        curve_type bez;
        pc.get(bez, pp);
        bez_deg=bez.degree();
        for (i=0; i<=bez_deg; ++i)
        {
          point_type pt(bez.get_control_point(i));

          if (show_segment_end_points && ((i==0) || (i==bez_deg)))
          {
            jointxbuf+=std::to_string(pt.x());
            jointybuf+=std::to_string(pt.y());
            if ((pp<(ns-1)) || ((pp==(ns-1)) && (i==0)))
            {
              jointxbuf+=", ";
              jointybuf+=", ";
            }
          }
          if (show_control_points)
          {
            cpxbuf+=std::to_string(pt.x());
            cpybuf+=std::to_string(pt.y());
            if ((pp<(ns-1)) || ((pp==(ns-1)) && (i<bez_deg)))
            {
              cpxbuf+=", ";
              cpybuf+=", ";
            }
          }
        }
        if (show_control_points)
        {
          cpxbuf+="];";
          cpybuf+="];";
        }
        if (show_segment_end_points)
        {
          jointxbuf+="];";
          jointybuf+="];";
        }
      }

      // initialize the t parameters
      const index_type npts(129);
      index_type nsegs(pc.number_segments()), nt;
      data_type dt;
      std::vector<data_type> tseg(nsegs+1);
      std::vector<index_type> npseg(nsegs);
      std::vector<data_type> t;

      pc.get_parameters(tseg.begin());
      dt=(tseg[nsegs]-tseg[0])/npts;

      calculate_t(t, npseg, tseg, dt);
      nt=t.size();

      // set the curve points
      cxbuf=nm+"_curv_x=[";
      cybuf=nm+"_curv_y=[";
      for (i=0; i<nt; ++i)
      {
        point_type pt(pc.f(t[i]));

        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        if (i<(nt-1))
        {
          cxbuf+=", ";
          cybuf+=", ";
        }
      }
      cxbuf+="];";
      cybuf+="];";

      std::cout << "% curve: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << cxbuf << std::endl;
      std::cout << cybuf << std::endl;
      if (show_control_points)
      {
        std::cout << cpxbuf << std::endl;
        std::cout << cpybuf << std::endl;
      }
      if (show_segment_end_points)
      {
        std::cout << jointxbuf << std::endl;
        std::cout << jointybuf << std::endl;
      }
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot(" << nm << "_curv_x, "
                           << nm << "_curv_y, "
                           << "'-g');" << std::endl;
      if (show_control_points)
      {
        std::cout << "plot(" << nm << "_curv_cp_x', "
                             << nm << "_curv_cp_y', "
                             << "'-o', 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0]);" << std::endl;
      }
      if (show_segment_end_points)
      {
        std::cout << "plot(" << nm << "_curv_joint_x', "
                             << nm << "_curv_joint_y', "
                             << "'s', 'Color', [0.75 0 0], 'MarkerFaceColor', [0.75 0 0]);" << std::endl;
      }
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &pc,
                      const std::string &name="", bool show_control_points=true, bool show_segment_end_points=true)
    {
      typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
      typedef typename piecewise_curve_type::curve_type curve_type;
      typedef typename piecewise_curve_type::data_type data_type;
      typedef typename piecewise_curve_type::point_type point_type;
      typedef typename piecewise_curve_type::index_type index_type;

      std::string nm, cpxbuf, cpybuf, cpzbuf, cxbuf, cybuf, czbuf, jointxbuf, jointybuf, jointzbuf;

      index_type i, pp, ns;

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
      }
      if (show_segment_end_points)
      {
        jointxbuf=nm+"_curv_joint_x=[";
        jointybuf=nm+"_curv_joint_y=[";
        jointzbuf=nm+"_curv_joint_z=[";
      }
      for (pp=0; pp<ns; ++pp)
      {
        index_type bez_deg;
        curve_type bez;
        pc.get(bez, pp);
        bez_deg=bez.degree();
        for (i=0; i<=bez_deg; ++i)
        {
          point_type pt(bez.get_control_point(i));

          if (show_segment_end_points && ((i==0) || (i==bez_deg)))
          {
            jointxbuf+=std::to_string(pt.x());
            jointybuf+=std::to_string(pt.y());
            jointzbuf+=std::to_string(pt.z());
            if ((pp<(ns-1)) || ((pp==(ns-1)) && (i==0)))
            {
              jointxbuf+=", ";
              jointybuf+=", ";
              jointzbuf+=", ";
            }
          }
          if (show_control_points)
          {
            cpxbuf+=std::to_string(pt.x());
            cpybuf+=std::to_string(pt.y());
            cpzbuf+=std::to_string(pt.z());
            if ((pp<(ns-1)) || ((pp==(ns-1)) && (i<bez_deg)))
            {
              cpxbuf+=", ";
              cpybuf+=", ";
              cpzbuf+=", ";
            }
          }
        }
      }
      if (show_control_points)
      {
        cpxbuf+="];";
        cpybuf+="];";
        cpzbuf+="];";
      }
      if (show_segment_end_points)
      {
        jointxbuf+="];";
        jointybuf+="];";
        jointzbuf+="];";
      }

      // initialize the t parameters
      const index_type npts(129);
      index_type nsegs(pc.number_segments()), nt;
      data_type dt;
      std::vector<data_type> tseg(nsegs+1);
      std::vector<index_type> npseg(nsegs);
      std::vector<data_type> t;

      pc.get_parameters(tseg.begin());
      dt=(tseg[nsegs]-tseg[0])/npts;

      calculate_t(t, npseg, tseg, dt);
      nt=t.size();

      // set the curve points
      cxbuf=nm+"_curv_x=[";
      cybuf=nm+"_curv_y=[";
      czbuf=nm+"_curv_z=[";
      for (i=0; i<nt; ++i)
      {
        point_type pt(pc.f(t[i]));

        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        czbuf+=std::to_string(pt.z());
        if (i<(nt-1))
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
      if (show_segment_end_points)
      {
        std::cout << jointxbuf << std::endl;
        std::cout << jointybuf << std::endl;
        std::cout << jointzbuf << std::endl;
      }
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot3(" << nm << "_curv_x, "
                            << nm << "_curv_y, "
                            << nm << "_curv_z, "
                            << "'-g');" << std::endl;
      if (show_control_points)
      {
        std::cout << "plot3(" << nm << "_curv_cp_x', "
                              << nm << "_curv_cp_y', "
                              << nm << "_curv_cp_z', "
                              << "'-o', 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0]);" << std::endl;
      }
      if (show_segment_end_points)
      {
        std::cout << "plot3(" << nm << "_curv_joint_x', "
                              << nm << "_curv_joint_y', "
                              << nm << "_curv_joint_z', "
                              << "'s', 'Color', [0.75 0 0], 'MarkerFaceColor', [0.75 0 0]);" << std::endl;
      }
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 2> &pc,
                      const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 2> &vec,
                      const std::string &name="")
    {
      typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 2> piecewise_curve_type;
      typedef typename piecewise_curve_type::data_type data_type;
      typedef typename piecewise_curve_type::point_type point_type;
      typedef typename piecewise_curve_type::index_type index_type;
      typedef typename piecewise_curve_type::tolerance_type tolerance_type;

      std::string nm(random_string(5)), vecxbuf, vecybuf, cxbuf, cybuf;

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
      vecxbuf=nm+"_vec_x=[";
      vecybuf=nm+"_vec_y=[";
      for (i=0; i<nt; ++i)
      {
        point_type pt(pc.f(t[i])), v(vec.f(t[i]));

        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        vecxbuf+=std::to_string(v.x());
        vecybuf+=std::to_string(v.y());
        if (i<(nt-1))
        {
          cxbuf+=", ";
          cybuf+=", ";
          vecxbuf+=", ";
          vecybuf+=", ";
        }
      }
      cxbuf+="];";
      cybuf+="];";
      vecxbuf+="];";
      vecybuf+="];";

      std::cout << "% curve: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << cxbuf << std::endl;
      std::cout << cybuf << std::endl;
      std::cout << vecxbuf << std::endl;
      std::cout << vecybuf << std::endl;
      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "quiver(" << nm << "_loc_x, "
                             << nm << "_loc_y, "
                             << nm << "_vec_x, "
                             << nm << "_vec_y, "
                             << "'r');" << std::endl;
    }

    template<typename data__>
    void octave_print(int figno, const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &pc,
                      const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &vec,
                      const std::string &name="")
    {
      typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
      typedef typename piecewise_curve_type::data_type data_type;
      typedef typename piecewise_curve_type::point_type point_type;
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
        point_type pt(pc.f(t[i])), v(vec.f(t[i]));

        cxbuf+=std::to_string(pt.x());
        cybuf+=std::to_string(pt.y());
        czbuf+=std::to_string(pt.z());
        vecxbuf+=std::to_string(v.x());
        vecybuf+=std::to_string(v.y());
        veczbuf+=std::to_string(v.z());
        if (i<(nt-1))
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
      typedef typename piecewise_surface_type::point_type point_type;
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
                point_type pt(bez.get_control_point(i, j));

                cpxbuf[ppqq]+=std::to_string(pt.x());
                cpybuf[ppqq]+=std::to_string(pt.y());
                cpzbuf[ppqq]+=std::to_string(pt.z());

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
          point_type pt(ps.f(u[i], v[j]));

          sxbuf+=std::to_string(pt.x());
          sybuf+=std::to_string(pt.y());
          szbuf+=std::to_string(pt.z());
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

