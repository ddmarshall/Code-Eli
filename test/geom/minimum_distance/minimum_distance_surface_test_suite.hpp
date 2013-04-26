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

#ifndef minimum_distance_surface_test_suite_hpp
#define minimum_distance_surface_test_suite_hpp

#include <cassert>  // assert()
#include <cmath>    // cos(), sin()

#include <typeinfo> // typeid

#include "Eigen/Eigen"

#include "eli/code_eli.hpp"

#include "eli/geom/surface/bezier.hpp"
#include "eli/geom/intersect/minimum_distance_surface.hpp"

template<typename data__>
class minimum_distance_surface_test_suite : public Test::Suite
{
  private:
    typedef data__ data_type;

    typedef eli::geom::surface::bezier<data_type, 3> surface_type;
    typedef typename surface_type::point_type point_type;
    typedef typename surface_type::index_type index_type;
    typedef typename surface_type::tolerance_type tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(minimum_distance_surface_test_suite<float>::point_smooth_test);
      TEST_ADD(minimum_distance_surface_test_suite<float>::point_closed_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_surface_test_suite<double>::point_smooth_test);
      TEST_ADD(minimum_distance_surface_test_suite<double>::point_closed_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_surface_test_suite<long double>::point_smooth_test);
      TEST_ADD(minimum_distance_surface_test_suite<long double>::point_closed_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(minimum_distance_surface_test_suite<dd_real>::point_smooth_test);
      TEST_ADD(minimum_distance_surface_test_suite<dd_real>::point_closed_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(minimum_distance_surface_test_suite<qd_real>::point_smooth_test);
      TEST_ADD(minimum_distance_surface_test_suite<qd_real>::point_closed_test);
    }
#endif
  public:
    minimum_distance_surface_test_suite()
    {
      AddTests(data__());
    }
    ~minimum_distance_surface_test_suite()
    {
    }

  private:
    void octave_print(int figno, const std::vector<point_type> &pts, const surface_type &s) const
    {
      size_t i, j;

      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << "xpts=[" << pts[0].x();
      for (i=1; i<pts.size(); ++i)
        std::cout << ", " << pts[i].x();
      std::cout << "];" << std::endl;
      std::cout << "ypts=[" << pts[0].y();
      for (i=1; i<pts.size(); ++i)
        std::cout << ", " << pts[i].y();
      std::cout << "];" << std::endl;
      std::cout << "zpts=[" << pts[0].z();
      for (i=1; i<pts.size(); ++i)
        std::cout << ", " << pts[i].z();
      std::cout << "];" << std::endl;

      // initialize the u & v parameters
      std::vector<data__> u(51), v(51);
      for (i=0; i<u.size(); ++i)
      {
        u[i]=static_cast<data__>(i)/(u.size()-1);
      }
      for (j=0; j<v.size(); ++j)
      {
        v[j]=static_cast<data__>(j)/(v.size()-1);
      }

      // set the surface points
      std::cout << "surf_x=[";
      for (i=0; i<u.size(); ++i)
      {
        std::cout << s.f(u[i], v[0]).x();
        for (j=1; j<(v.size()-1); ++j)
        {
          std::cout << ", " << s.f(u[i], v[j]).x();
        }
        j=v.size()-1;
        std::cout << ", " << s.f(u[i], v[j]).x();
        if (i<(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_y=[";
      for (i=0; i<u.size(); ++i)
      {
        std::cout << s.f(u[i], v[0]).y();
        for (j=1; j<(v.size()-1); ++j)
        {
          std::cout << ", " << s.f(u[i], v[j]).y();
        }
        j=v.size()-1;
        std::cout << ", " << s.f(u[i], v[j]).y();
        if (i<(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_z=[";
      for (i=0; i<u.size(); ++i)
      {
        std::cout << s.f(u[i], v[0]).z();
        for (j=1; j<(v.size()-1); ++j)
        {
          std::cout << ", " << s.f(u[i], v[j]).z();
        }
        j=v.size()-1;
        std::cout << ", " << s.f(u[i], v[j]).z();
        if (i<(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "mesh(surf_x, surf_y, surf_z, zeros(size(surf_z)), 'EdgeColor', [0 0 0]);" << std::endl;
      std::cout << "hold on;" << std::endl;
      std::cout << "plot3(xpts, ypts, zpts, 'or', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');" << std::endl;
      std::cout << "hold off;" << std::endl;

      std::cout << "figure(" << figno+1 << ");" << std::endl;
      std::cout << "u=[" << u[0];
      for (i=1; i<u.size(); ++i)
        std::cout << ", " << u[i];
      std::cout << "];" << std::endl;
      std::cout << "v=[" << v[0];
      for (i=1; i<v.size(); ++i)
        std::cout << ", " << v[i];
      std::cout << "];" << std::endl;
      std::cout << "g_u=[";
      for (i=0; i<u.size(); ++i)
      {
        j=0;
        std::cout << (s.f(u[i], v[j])-pts[0]).dot(s.f_u(u[i], v[j]));
        for (j=1; j<(v.size()-1); ++j)
        {
          std::cout << ", " << (s.f(u[i], v[j])-pts[0]).dot(s.f_u(u[i], v[j]));
        }
        j=v.size()-1;
        std::cout << ", " << (s.f(u[i], v[j])-pts[0]).dot(s.f_u(u[i], v[j]));
        if (i<(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;
      std::cout << "g_v=[";
      for (i=0; i<u.size(); ++i)
      {
        j=0;
        std::cout << (s.f(u[i], v[j])-pts[0]).dot(s.f_v(u[i], v[j]));
        for (j=1; j<(v.size()-1); ++j)
        {
          std::cout << ", " << (s.f(u[i], v[j])-pts[0]).dot(s.f_v(u[i], v[j]));
        }
        j=v.size()-1;
        std::cout << ", " << (s.f(u[i], v[j])-pts[0]).dot(s.f_v(u[i], v[j]));
        if (i<(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "mesh(u, v, g_u, zeros(length(u), length(v)), 'EdgeColor', [0 0 0]);" << std::endl;
      std::cout << "hold on;" << std::endl;
      std::cout << "mesh(u, v, g_v, zeros(length(u), length(v)), 'EdgeColor', [0 0 0]);" << std::endl;
      std::cout << "surf(u, v, zeros(length(u), length(v)));" << std::endl;
      std::cout << "hold off;" << std::endl;
    }

    void point_smooth_test()
    {
      const index_type n(3), m(3);
      surface_type s(n, m);
      point_type cp[3+1][3+1], pt_out;
      point_type pt, norm, u_contra, v_contra;
      data_type  dist, u, v, dist_ref, u_ref, v_ref;
      index_type i, j;

      // create surface with specified control points
      cp[0][0] << -15, 0,  15;
      cp[1][0] <<  -5, 5,  15;
      cp[2][0] <<   5, 5,  15;
      cp[3][0] <<  15, 0,  15;
      cp[0][1] << -15, 5,   5;
      cp[1][1] <<  -5, 5,   5;
      cp[2][1] <<   5, 5,   5;
      cp[3][1] <<  15, 5,   5;
      cp[0][2] << -15, 5,  -5;
      cp[1][2] <<  -5, 5,  -5;
      cp[2][2] <<   5, 5,  -5;
      cp[3][2] <<  15, 5,  -5;
      cp[0][3] << -15, 0, -15;
      cp[1][3] <<  -5, 5, -15;
      cp[2][3] <<   5, 5, -15;
      cp[3][3] <<  15, 0, -15;

      // create surface with specified dimensions and set control points
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          s.set_control_point(cp[i][j], i, j);
        }
      }

      // test point on surface
      dist_ref=0;
      u_ref=0.25;
      v_ref=0.25;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point very near surface
      dist_ref=0.01;
      u_ref=0.25;
      v_ref=0.25;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=0.1;
      u_ref=0.25;
      v_ref=0.25;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near and on concave side of surface
      dist_ref=0.1;
      u_ref=0.25;
      v_ref=0.25;
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point far surface
      dist_ref=1.1;
      u_ref=0.25;
      v_ref=0.25;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point on surface
      dist_ref=0;
      u_ref=0.64;
      v_ref=0.32;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point very near surface
      dist_ref=0.01;
      u_ref=0.64;
      v_ref=0.32;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=0.1;
      u_ref=0.64;
      v_ref=0.32;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near and on concave side of surface
      dist_ref=0.1;
      u_ref=0.64;
      v_ref=0.32;
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=1.1;
      u_ref=0.64;
      v_ref=0.32;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=0.1;
      u_ref=0.01;
      v_ref=0.01;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=0.1;
      u_ref=0.99;
      v_ref=0.01;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=0.1;
      u_ref=0.01;
      v_ref=0.999;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=0.1;
      u_ref=0.999;
      v_ref=0.999;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point beyond edge of surface
      dist_ref=0.1;
      u_ref=1;
      v_ref=0.4;
      norm=s.normal(u_ref, v_ref);
      u_contra=-norm.cross(s.f_v(u_ref, v_ref));
      u_contra.normalize();
      pt=0.01*u_contra+s.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, s.f(u_ref, v_ref));
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point beyond edge of surface
      dist_ref=0.1;
      u_ref=0.6;
      v_ref=1;
      norm=s.normal(u_ref, v_ref);
      v_contra=norm.cross(s.f_u(u_ref, v_ref));
      v_contra.normalize();
      pt=0.01*v_contra+s.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, s.f(u_ref, v_ref));
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point beyond edge of surface
      dist_ref=0.5;
      u_ref=0;
      v_ref=0.2;
      norm=s.normal(u_ref, v_ref);
      u_contra=-norm.cross(s.f_v(u_ref, v_ref));
      u_contra.normalize();
      pt=0.1*u_contra+s.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, s.f(u_ref, v_ref));
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point beyond edge of surface
      dist_ref=0.5;
      u_ref=0.2;
      v_ref=0;
      norm=s.normal(u_ref, v_ref);
      v_contra=norm.cross(s.f_u(u_ref, v_ref));
      v_contra.normalize();
      pt=0.1*v_contra+s.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, s.f(u_ref, v_ref));
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
    }

    void point_closed_test()
    {
#if 0
      if (typeid(data_type)==typeid(double))
      {
        std::cout << "% s_f(" << u_ref << ", " << v_ref << ")=" << s.f(u_ref, v_ref)
                  << "\ts_f(" << u << ", " << v << ")=" << s.f(u, v) << "\tpt=" << pt << std::endl;
        std::cout << "% dist=" << dist << "\tdist_ref=" << dist_ref
                  << "\t(u, v)=(" << u << ", " << v << ")"
                  << "\t(u_ref, v_ref)=(" << u_ref << ", " << v_ref << ")" << std::endl;
        std::vector<point_type> vec(2);
        vec[0]=pt;
        vec[1]=s.f(u_ref, v_ref);
//         octave_print(1, vec, s);
      }
      point_type3 cntrl_in[13];

      // set control points
      cntrl_in[0]  <<   1, 0,   0;
      cntrl_in[1]  <<   1, 0.5,-0.5;
      cntrl_in[2]  << 0.5, 1,  -0.5;
      cntrl_in[3]  <<   0, 1,   0;
      cntrl_in[4]  <<-0.5, 1,   0;
      cntrl_in[5]  <<  -1, 0.5, 0.5;
      cntrl_in[6]  <<  -1, 0,   0.5;
      cntrl_in[7]  <<  -1,-0.5, 0;
      cntrl_in[8]  <<-0.5,-1,  -0.5;
      cntrl_in[9]  <<   0,-1,   0;
      cntrl_in[10] << 0.5,-1,   0.5;
      cntrl_in[11] <<   1,-0.5, 0.5;
      cntrl_in[12] <<   1,0,    0;

      curve_type3 c(12);
      point_type3 pt, norm, fp;
      data_type  dist, t, dist_ref, t_ref;

      // set control points
      for (index_type3 i=0; i<13; ++i)
      {
        c.set_control_point(cntrl_in[i], i);
      }

      // test point on curve
      dist_ref=0;
      t_ref=0.25;
      fp=c.fp(t_ref);
      norm << fp.z(), -fp.z(), -fp.x()+fp.y();
      norm.normalize();
      pt=c.f(t_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(t, c, pt);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near curve
      dist_ref=0.1;
      t_ref=0.25;
      fp=c.fp(t_ref);
      norm << fp.z(), -fp.z(), -fp.x()+fp.y();
      norm.normalize();
      pt=c.f(t_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(t, c, pt);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near center
      dist_ref=0.77;
      t_ref=0.25;
      fp=c.fp(t_ref);
      norm << fp.z(), -fp.z(), -fp.x()+fp.y();
      norm.normalize();
      pt=c.f(t_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(t, c, pt);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near and outside curve
      dist_ref=0.1;
      t_ref=0.7;
      fp=c.fp(t_ref);
      norm << fp.z(), -fp.z(), -fp.x()+fp.y();
      norm.normalize();
      pt=c.f(t_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(t, c, pt);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point far and outside curve
      dist_ref=3;
      t_ref=0.7;
      fp=c.fp(t_ref);
      norm << fp.z(), -fp.z(), -fp.x()+fp.y();
      norm.normalize();
      pt=c.f(t_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(t, c, pt);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near end of curve
      dist_ref=0.1;
      t_ref=0.999;
      fp=c.fp(t_ref);
      norm << fp.z(), -fp.z(), -fp.x()+fp.y();
      norm.normalize();
      pt=c.f(t_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(t, c, pt);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near end of and outside curve
      dist_ref=0.1;
      t_ref=0.999;
      fp=c.fp(t_ref);
      norm << -fp.z(), fp.z(), fp.x()-fp.y();
      norm.normalize();
      pt=c.f(t_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(t, c, pt);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near start of curve
      dist_ref=0.1;
      t_ref=0.001;
      fp=c.fp(t_ref);
      norm << fp.z(), -fp.z(), -fp.x()+fp.y();
      norm.normalize();
      pt=c.f(t_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(t, c, pt);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near start of and outside curve
      dist_ref=0.1;
      t_ref=0.001;
      fp=c.fp(t_ref);
      norm << -fp.z(), fp.z(), fp.x()-fp.y();
      norm.normalize();
      pt=c.f(t_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(t, c, pt);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near start of and far outside curve
      dist_ref=3;
      t_ref=0.001;
      fp=c.fp(t_ref);
      norm << -fp.z(), fp.z(), fp.x()-fp.y();
      norm.normalize();
      pt=c.f(t_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(t, c, pt);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
//       if (typeid(data_type)==typeid(double))
//       {
//         std::cout << "c_f(t_ref)=" << c.f(t_ref) << "\tpt=" << pt << std::endl;
//         std::cout << "dist=" << dist << "\tdist_ref=" << dist_ref << "\tt=" << t << "\tt_ref=" << t_ref << std::endl;
//         std::cout << "dist to end=" << (cntrl_in[0]-pt).norm() << std::endl;
//         std::vector<point_type3> vec(2);
//         vec[0]=pt;
//         vec[1]=c.f(t);
//         octave_print(1, vec, c);
//       }
#endif
    }
};

#endif

