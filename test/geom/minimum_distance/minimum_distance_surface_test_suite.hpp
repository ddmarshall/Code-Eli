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
#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/intersect/minimum_distance_surface.hpp"

template<typename data__>
class minimum_distance_surface_test_suite : public Test::Suite
{
  private:
    typedef data__ data_type;

    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data_type, 3> piecewise_surface_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
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
      TEST_ADD(minimum_distance_surface_test_suite<float>::point_piecewise_01_smooth_test);
      TEST_ADD(minimum_distance_surface_test_suite<float>::point_piecewise_trange_smooth_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_surface_test_suite<double>::point_smooth_test);
      TEST_ADD(minimum_distance_surface_test_suite<double>::point_closed_test);
      TEST_ADD(minimum_distance_surface_test_suite<double>::point_piecewise_01_smooth_test);
      TEST_ADD(minimum_distance_surface_test_suite<double>::point_piecewise_trange_smooth_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_surface_test_suite<long double>::point_smooth_test);
      TEST_ADD(minimum_distance_surface_test_suite<long double>::point_closed_test);
      TEST_ADD(minimum_distance_surface_test_suite<long double>::point_piecewise_01_smooth_test);
      TEST_ADD(minimum_distance_surface_test_suite<long double>::point_piecewise_trange_smooth_test);
    }

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
      std::cout << "axis equal" << std::endl;

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
      data_type dist, u, v, dist_ref, u_ref, v_ref;
      data_type u_off(0.2), v_off(0.2);
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
      TEST_ASSERT(s.open_u());
      TEST_ASSERT(s.open_v());

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

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point very near surface
      dist_ref=static_cast<data_type>(0.01);
      u_ref=0.25;
      v_ref=0.25;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=0.25;
      v_ref=0.25;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near and on concave side of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=0.25;
      v_ref=0.25;
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point far surface
      dist_ref=static_cast<data_type>(1.1);
      u_ref=0.25;
      v_ref=0.25;
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point on surface
      dist_ref=0;
      u_ref=static_cast<data_type>(0.64);
      v_ref=static_cast<data_type>(0.32);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point very near surface
      dist_ref=static_cast<data_type>(0.01);
      u_ref=static_cast<data_type>(0.64);
      v_ref=static_cast<data_type>(0.32);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.64);
      v_ref=static_cast<data_type>(0.32);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near and on concave side of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.64);
      v_ref=static_cast<data_type>(0.32);
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=static_cast<data_type>(1.1);
      u_ref=static_cast<data_type>(0.64);
      v_ref=static_cast<data_type>(0.32);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.01);
      v_ref=static_cast<data_type>(0.01);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.99);
      v_ref=static_cast<data_type>(0.01);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.01);
      v_ref=static_cast<data_type>(0.999);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.999);
      v_ref=static_cast<data_type>(0.999);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(1);
      v_ref=static_cast<data_type>(0.4);
      norm=s.normal(u_ref, v_ref);
      u_contra=-norm.cross(s.f_v(u_ref, v_ref));
      u_contra.normalize();
      pt=static_cast<data_type>(0.01)*u_contra+s.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, s.f(u_ref, v_ref));
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      data_type u_alt, v_alt, dist_alt;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref+v_off);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.6);
      v_ref=static_cast<data_type>(1);
      norm=s.normal(u_ref, v_ref);
      v_contra=norm.cross(s.f_u(u_ref, v_ref));
      v_contra.normalize();
      pt=static_cast<data_type>(0.01)*v_contra+s.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, s.f(u_ref, v_ref));
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.5);
      u_ref=static_cast<data_type>(0);
      v_ref=static_cast<data_type>(0.2);
      norm=s.normal(u_ref, v_ref);
      u_contra=-norm.cross(s.f_v(u_ref, v_ref));
      u_contra.normalize();
      pt=static_cast<data_type>(-0.1)*u_contra+s.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, s.f(u_ref, v_ref));
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.5);
      u_ref=static_cast<data_type>(0.2);
      v_ref=static_cast<data_type>(0);
      norm=s.normal(u_ref, v_ref);
      v_contra=norm.cross(s.f_u(u_ref, v_ref));
      v_contra.normalize();
      pt=static_cast<data_type>(-0.1)*v_contra+s.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, s.f(u_ref, v_ref));
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref+.02);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref+.02);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

    }

    void point_closed_test()
    {
      const index_type n(12), m(3);
      surface_type s(n, m);
      point_type cp[n+1], offset[m+1], pt_out;
      point_type pt, norm, u_contra, v_contra;
      data_type  z[m+1], dist, u, v, dist_ref, u_ref, v_ref;
      data_type u_off(0.1), v_off(0.1);
      index_type i, j;

      // set control points
      cp[0]  <<   1, 0,   0;
      cp[1]  <<   1, 0.5, 0;
      cp[2]  << 0.5, 1,   0;
      cp[3]  <<   0, 1,   0;
      cp[4]  <<-0.5, 1,   0;
      cp[5]  <<  -1, 0.5, 0;
      cp[6]  <<  -1, 0,   0;
      cp[7]  <<  -1,-0.5, 0;
      cp[8]  <<-0.5,-1,   0;
      cp[9]  <<   0,-1,   0;
      cp[10] << 0.5,-1,   0;
      cp[11] <<   1,-0.5, 0;
      cp[12] <<   1,0,    0;
      z[0]=-0.5;
      z[1]=0;
      z[2]=0.5;
      z[3]=1;
      offset[0] << 0, 0, 0;
      offset[1] << 0.5, 0, 0;
      offset[2] << -0.5, 0.5, 0;
      offset[3] << 0, 0, 0;

      // create surface with specified dimensions and set control points
      for (j=0; j<=m; ++j)
      {
        for (i=0; i<=n; ++i)
        {
          cp[i].z() = z[j];
          s.set_control_point(cp[i]+offset[j], i, j);
        }
      }
      TEST_ASSERT(s.closed_u());
      TEST_ASSERT(s.open_v());

      // test point on surface
      dist_ref=0;
      u_ref=static_cast<data_type>(0.35);
      v_ref=static_cast<data_type>(0.35);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point very near surface
      dist_ref=static_cast<data_type>(0.01);
      u_ref=static_cast<data_type>(0.35);
      v_ref=static_cast<data_type>(0.35);
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.35);
      v_ref=static_cast<data_type>(0.35);
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near center surface
      dist_ref=static_cast<data_type>(0.72);
      u_ref=static_cast<data_type>(0.35);
      v_ref=static_cast<data_type>(0.35);
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // +v_off cases seek local maximum.
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near and outside surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.35);
      v_ref=static_cast<data_type>(0.35);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point far and outside surface
      dist_ref=3;
      u_ref=static_cast<data_type>(0.25);
      v_ref=static_cast<data_type>(0.7);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near end of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.25);
      v_ref=static_cast<data_type>(0.999);
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near end and outside of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.25);
      v_ref=static_cast<data_type>(0.999);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near end of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.25);
      v_ref=static_cast<data_type>(0.001);
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near end and outside of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.25);
      v_ref=static_cast<data_type>(0.001);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near closed portion of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.001);
      v_ref=static_cast<data_type>(0.3);
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near closed portion and outside of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.001);
      v_ref=static_cast<data_type>(0.3);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 1-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near closed portion of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.999);
      v_ref=static_cast<data_type>(0.3);
      norm=-s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near closed portion and outside of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.999);
      v_ref=static_cast<data_type>(0.3);
      norm=s.normal(u_ref, v_ref);
      pt=s.f(u_ref, v_ref)+dist_ref*norm;
      dist=eli::geom::intersect::minimum_distance(u, v, s, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, 0+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, s, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

//       if (typeid(data_type)==typeid(double))
//       {
//         std::cout << "% s_f(" << u_ref << ", " << v_ref << ")=" << s.f(u_ref, v_ref)
//                   << "\ts_f(" << u << ", " << v << ")=" << s.f(u, v) << "\tpt=" << pt << std::endl;
//         std::cout << "% dist=" << dist << "\tdist_ref=" << dist_ref
//                   << "\t(u, v)=(" << u << ", " << v << ")"
//                   << "\t(u_ref, v_ref)=(" << u_ref << ", " << v_ref << ")" << std::endl;
//         std::vector<point_type> vec(3);
//         vec[0]=pt;
//         vec[1]=s.f(u_ref, v_ref);
//         vec[2]=s.f(u, v);
//         octave_print(1, vec, s);
//       }
    }

    void point_piecewise_01_smooth_test()
    {
      const index_type n(3), m(3);
      surface_type s(n, m);
      point_type cp[3+1][3+1], pt_out;
      point_type pt, norm, u_contra, v_contra;
      data_type dist, u, v, dist_ref, u_ref, v_ref;
      data_type u_off(0.2), v_off(0.2);
      index_type i, j;
      piecewise_surface_type pws;
      typename piecewise_surface_type::error_code err;

      pws.init_uv(1, 1);

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
      TEST_ASSERT(s.open_u());
      TEST_ASSERT(s.open_v());

      err=pws.set(s, 0, 0);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

      pws.split_u(0.2);
      pws.split_u(0.4);
      pws.split_u(0.6);
      pws.split_u(0.8);

      pws.split_v(0.2);
      pws.split_v(0.4);
      pws.split_v(0.6);
      pws.split_v(0.8);

      // test point on surface
      dist_ref=0;
      u_ref=0.25;
      v_ref=0.25;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point very near surface
      dist_ref=static_cast<data_type>(0.01);
      u_ref=0.25;
      v_ref=0.25;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=0.25;
      v_ref=0.25;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near and on concave side of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=0.25;
      v_ref=0.25;
      norm=-pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point far surface
      dist_ref=static_cast<data_type>(1.1);
      u_ref=0.25;
      v_ref=0.25;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point on surface
      dist_ref=0;
      u_ref=static_cast<data_type>(0.64);
      v_ref=static_cast<data_type>(0.32);
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point very near surface
      dist_ref=static_cast<data_type>(0.01);
      u_ref=static_cast<data_type>(0.64);
      v_ref=static_cast<data_type>(0.32);
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.64);
      v_ref=static_cast<data_type>(0.32);
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near and on concave side of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.64);
      v_ref=static_cast<data_type>(0.32);
      norm=-pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=static_cast<data_type>(1.1);
      u_ref=static_cast<data_type>(0.64);
      v_ref=static_cast<data_type>(0.32);
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.01);
      v_ref=static_cast<data_type>(0.01);
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.99);
      v_ref=static_cast<data_type>(0.01);
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.01);
      v_ref=static_cast<data_type>(0.999);
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.999);
      v_ref=static_cast<data_type>(0.999);
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(1);
      v_ref=static_cast<data_type>(0.4);
      norm=pws.normal(u_ref, v_ref);
      u_contra=-norm.cross(pws.f_v(u_ref, v_ref));
      u_contra.normalize();
      pt=static_cast<data_type>(0.01)*u_contra+pws.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, pws.f(u_ref, v_ref));
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      data_type u_alt, v_alt, dist_alt;
      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+v_off);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.6);
      v_ref=static_cast<data_type>(1);
      norm=pws.normal(u_ref, v_ref);
      v_contra=norm.cross(pws.f_u(u_ref, v_ref));
      v_contra.normalize();
      pt=static_cast<data_type>(0.01)*v_contra+pws.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, pws.f(u_ref, v_ref));
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.5);
      u_ref=static_cast<data_type>(0);
      v_ref=static_cast<data_type>(0.2);
      norm=pws.normal(u_ref, v_ref);
      u_contra=-norm.cross(pws.f_v(u_ref, v_ref));
      u_contra.normalize();
      pt=static_cast<data_type>(-0.1)*u_contra+pws.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, pws.f(u_ref, v_ref));
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.5);
      u_ref=static_cast<data_type>(0.2);
      v_ref=static_cast<data_type>(0);
      norm=pws.normal(u_ref, v_ref);
      v_contra=norm.cross(pws.f_u(u_ref, v_ref));
      v_contra.normalize();
      pt=static_cast<data_type>(-0.1)*v_contra+pws.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, pws.f(u_ref, v_ref));
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+.02);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+.02);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1, 1);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }
    }

    void point_piecewise_trange_smooth_test()
    {
      const index_type n(3), m(3);
      surface_type s(n, m);
      point_type cp[3+1][3+1], pt_out;
      point_type pt, norm, u_contra, v_contra;
      data_type dist, u, v, dist_ref, u_ref, v_ref;
      data_type u_off(0.2), v_off(0.2);
      index_type i, j;
      piecewise_surface_type pws;
      typename piecewise_surface_type::error_code err;

      pws.init_uv(1, 1);

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
      TEST_ASSERT(s.open_u());
      TEST_ASSERT(s.open_v());

      err=pws.set(s, 0, 0);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

      data_type u0(1.4);
      data_type v0(2.2);

      pws.set_u0(u0);
      pws.set_v0(v0);

      pws.split_u(u0+0.2);
      pws.split_u(u0+0.4);
      pws.split_u(u0+0.6);
      pws.split_u(u0+0.8);

      pws.split_v(v0+0.2);
      pws.split_v(v0+0.4);
      pws.split_v(v0+0.6);
      pws.split_v(v0+0.8);

      // test point on surface
      dist_ref=0;
      u_ref=0.25+u0;
      v_ref=0.25+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point very near surface
      dist_ref=static_cast<data_type>(0.01);
      u_ref=0.25+u0;
      v_ref=0.25+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=0.25+u0;
      v_ref=0.25+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near and on concave side of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=0.25+u0;
      v_ref=0.25+v0;
      norm=-pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point far surface
      dist_ref=static_cast<data_type>(1.1);
      u_ref=0.25+u0;
      v_ref=0.25+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point on surface
      dist_ref=0;
      u_ref=static_cast<data_type>(0.64)+u0;
      v_ref=static_cast<data_type>(0.32)+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point very near surface
      dist_ref=static_cast<data_type>(0.01);
      u_ref=static_cast<data_type>(0.64)+u0;
      v_ref=static_cast<data_type>(0.32)+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.64)+u0;
      v_ref=static_cast<data_type>(0.32)+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near and on concave side of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.64)+u0;
      v_ref=static_cast<data_type>(0.32)+v0;
      norm=-pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near surface
      dist_ref=static_cast<data_type>(1.1);
      u_ref=static_cast<data_type>(0.64)+u0;
      v_ref=static_cast<data_type>(0.32)+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.01)+u0;
      v_ref=static_cast<data_type>(0.01)+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.99)+u0;
      v_ref=static_cast<data_type>(0.01)+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.01)+u0;
      v_ref=static_cast<data_type>(0.999)+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point near corner of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.999)+u0;
      v_ref=static_cast<data_type>(0.999)+v0;
      norm=pws.normal(u_ref, v_ref);
      pt=pws.f(u_ref, v_ref)+dist_ref*norm;
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(1)+u0;
      v_ref=static_cast<data_type>(0.4)+v0;
      norm=pws.normal(u_ref, v_ref);
      u_contra=-norm.cross(pws.f_v(u_ref, v_ref));
      u_contra.normalize();
      pt=static_cast<data_type>(0.01)*u_contra+pws.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, pws.f(u_ref, v_ref));
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      data_type u_alt, v_alt, dist_alt;
      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+v_off);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.1);
      u_ref=static_cast<data_type>(0.6)+u0;
      v_ref=static_cast<data_type>(1)+v0;
      norm=pws.normal(u_ref, v_ref);
      v_contra=norm.cross(pws.f_u(u_ref, v_ref));
      v_contra.normalize();
      pt=static_cast<data_type>(0.01)*v_contra+pws.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, pws.f(u_ref, v_ref));
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.5);
      u_ref=static_cast<data_type>(0)+u0;
      v_ref=static_cast<data_type>(0.2)+v0;
      norm=pws.normal(u_ref, v_ref);
      u_contra=-norm.cross(pws.f_v(u_ref, v_ref));
      u_contra.normalize();
      pt=static_cast<data_type>(-0.1)*u_contra+pws.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, pws.f(u_ref, v_ref));
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref-v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }

      // test point beyond edge of surface
      dist_ref=static_cast<data_type>(0.5);
      u_ref=static_cast<data_type>(0.2)+u0;
      v_ref=static_cast<data_type>(0)+v0;
      norm=pws.normal(u_ref, v_ref);
      v_contra=norm.cross(pws.f_u(u_ref, v_ref));
      v_contra.normalize();
      pt=static_cast<data_type>(-0.1)*v_contra+pws.f(u_ref, v_ref)+dist_ref*norm;
      dist_ref=eli::geom::point::distance(pt, pws.f(u_ref, v_ref));
//      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt);
//      TEST_ASSERT(tol.approximately_equal(u, u_ref));
//      TEST_ASSERT(tol.approximately_equal(v, v_ref));
//      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // Surface solver converges to nearby apparent solution.
      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+.02);
      u_alt=u;
      v_alt=v;
      dist_alt=dist;

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref, v_ref+.02);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref+u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref+v_off);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, u_ref-u_off, v_ref);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      TEST_ASSERT(tol.approximately_equal(dist, dist_alt));

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 0+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 0+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }

      dist=eli::geom::intersect::minimum_distance(u, v, pws, pt, 1+u0, 1+v0);
      TEST_ASSERT(tol.approximately_equal(u, u_alt));
      TEST_ASSERT(tol.approximately_equal(v, v_alt));
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT(std::abs(dist-dist_alt) < 2e-3);
      }
      else
      {
        TEST_ASSERT(tol.approximately_equal(dist, dist_alt));
      }
    }


};

#endif

