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

#ifndef piecewise_surface_test_suite_hpp
#define piecewise_surface_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/util/tolerance.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/surface/bezier.hpp"
// #include "eli/geom/surface/area.hpp"
#include "eli/geom/surface/curvature.hpp"
#include "eli/geom/surface/piecewise.hpp"

template<typename data__>
class piecewise_surface_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;
    typedef typename piecewise_surface_type::tolerance_type tolerance_type;
    typedef typename piecewise_surface_type::piecewise_curve_type piecewise_curve_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_test_suite<float>::creation_test);
      TEST_ADD(piecewise_surface_test_suite<float>::bounding_box_test);
      TEST_ADD(piecewise_surface_test_suite<float>::reverse_test);
      TEST_ADD(piecewise_surface_test_suite<float>::swap_test);
      TEST_ADD(piecewise_surface_test_suite<float>::replace_test);
      TEST_ADD(piecewise_surface_test_suite<float>::transformation_test);
      TEST_ADD(piecewise_surface_test_suite<float>::evaluation_test);
      TEST_ADD(piecewise_surface_test_suite<float>::split_test);
      TEST_ADD(piecewise_surface_test_suite<float>::to_cubic_test);
      TEST_ADD(piecewise_surface_test_suite<float>::area_test);
      TEST_ADD(piecewise_surface_test_suite<float>::get_curve_test);
      TEST_ADD(piecewise_surface_test_suite<float>::continuity_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_test_suite<double>::creation_test);
      TEST_ADD(piecewise_surface_test_suite<double>::bounding_box_test);
      TEST_ADD(piecewise_surface_test_suite<double>::reverse_test);
      TEST_ADD(piecewise_surface_test_suite<double>::swap_test);
      TEST_ADD(piecewise_surface_test_suite<double>::replace_test);
      TEST_ADD(piecewise_surface_test_suite<double>::transformation_test);
      TEST_ADD(piecewise_surface_test_suite<double>::evaluation_test);
      TEST_ADD(piecewise_surface_test_suite<double>::split_test);
      TEST_ADD(piecewise_surface_test_suite<double>::to_cubic_test);
      TEST_ADD(piecewise_surface_test_suite<double>::area_test);
      TEST_ADD(piecewise_surface_test_suite<double>::get_curve_test);
      TEST_ADD(piecewise_surface_test_suite<double>::continuity_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_test_suite<long double>::creation_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::bounding_box_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::reverse_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::swap_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::replace_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::transformation_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::evaluation_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::split_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::to_cubic_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::area_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::get_curve_test);
      TEST_ADD(piecewise_surface_test_suite<long double>::continuity_test);
    }

  public:
    piecewise_surface_test_suite()
    {
      AddTests(data__());
    }
    ~piecewise_surface_test_suite()
    {
    }

  private:
    void octave_print(int figno, /*const point_type pts[][], */const piecewise_surface_type &ps) const
    {
      index_type i, j, pp, qq, nup, nvp;
      data_type umin, vmin, umax, vmax;

      nup=ps.number_u_patches();
      nvp=ps.number_v_patches();
      ps.get_parameter_min(umin, vmin);
      ps.get_parameter_max(umax, vmax);

      std::cout << "figure(" << figno << ");" << std::endl;
//      std::cout << "cp_x=[" << std::endl;
//      for (pp=0; pp<nup; ++pp)
//      {
//        for (qq=0; qq<nvp; ++qq)
//        {
//          surface_type bez;
//          ps.get(bez, pp, qq);
//          for (i=0; i<=bez.degree_u(); ++i)
//          {
//            std::cout << bez.get_control_point(i, 0).x();
//            for (j=1; j<bez.degree_v(); ++j)
//            {
//              std::cout << ", " << bez.get_control_point(i, j).x();
//            }
//            j=bez.degree_v();
//            std::cout << ", " << bez.get_control_point(i, j).x();
//            if (i<bez.degree_u())
//              std::cout << "; ";
//            else if (pp<nup-1)
//              std::cout << "; ";
//            else if (qq<nvp-1)
//              std::cout << "; ";
//          }
//          std::cout << std::endl;
//        }
//      }
//      std::cout << "];" << std::endl;
//
//      std::cout << "cp_y=[";
//      for (pp=0; pp<nup; ++pp)
//      {
//        for (qq=0; qq<nvp; ++qq)
//        {
//          surface_type bez;
//          ps.get(bez, pp, qq);
//          for (i=0; i<=bez.degree_u(); ++i)
//          {
//            std::cout << bez.get_control_point(i, 0).y();
//            for (j=1; j<bez.degree_v(); ++j)
//            {
//              std::cout << ", " << bez.get_control_point(i, j).y();
//            }
//            j=bez.degree_v();
//            std::cout << ", " << bez.get_control_point(i, j).y();
//            if (i<bez.degree_u())
//              std::cout << "; ";
//            else if (pp<nup-1)
//              std::cout << "; ";
//            else if (qq<nvp-1)
//              std::cout << "; ";
//          }
//          std::cout << std::endl;
//        }
//      }
//      std::cout << "];" << std::endl;
//
//      std::cout << "cp_z=[";
//      for (pp=0; pp<nup; ++pp)
//      {
//        for (qq=0; qq<nvp; ++qq)
//        {
//          surface_type bez;
//          ps.get(bez, pp, qq);
//          for (i=0; i<=bez.degree_u(); ++i)
//          {
//            std::cout << bez.get_control_point(i, 0).z();
//            for (j=1; j<bez.degree_v(); ++j)
//            {
//              std::cout << ", " << bez.get_control_point(i, j).z();
//            }
//            j=bez.degree_v();
//            std::cout << ", " << bez.get_control_point(i, j).z();
//            if (i<bez.degree_u())
//              std::cout << "; ";
//            else if (pp<nup-1)
//              std::cout << "; ";
//            else if (qq<nvp-1)
//              std::cout << "; ";
//          }
//          std::cout << std::endl;
//        }
//      }
//      std::cout << "];" << std::endl;

      // initialize the u & v parameters
      std::vector<data__> u(31), v(31);
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        u[i]=umin+(umax-umin)*static_cast<data__>(i)/(u.size()-1);
      }
      for (j=0; j<static_cast<index_type>(v.size()); ++j)
      {
        v[j]=vmin+(vmax-vmin)*static_cast<data__>(j)/(v.size()-1);
      }

      // set the surface points
      std::cout << "surf_x=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << ps.f(u[i], v[0]).x();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << ps.f(u[i], v[j]).x();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << ps.f(u[i], v[j]).x();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_y=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << ps.f(u[i], v[0]).y();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << ps.f(u[i], v[j]).y();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << ps.f(u[i], v[j]).y();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_z=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << ps.f(u[i], v[0]).z();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << ps.f(u[i], v[j]).z();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << ps.f(u[i], v[j]).z();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "mesh(surf_x, surf_y, surf_z, zeros(size(surf_z)), 'EdgeColor', [0 0 0]);" << std::endl;
//      std::cout << "hold on;" << std::endl;
//      std::cout << "plot3(cp_x, cp_y, cp_z, 'ok', 'MarkerFaceColor', [0 0 0]);" << std::endl;
//      std::cout << "hold off;" << std::endl;
    }

    void creation_test()
    {
      surface_type s, s1, s2, s3, s4, s5, s6, s_patches[6];
      piecewise_surface_type ps1, ps2;
      typename piecewise_surface_type::error_code err;

      // create 3x2 patches with unit spacing
      ps1.init_uv(3, 2);
      TEST_ASSERT(ps1.number_u_patches()==3);
      TEST_ASSERT(ps1.number_v_patches()==2);

      // create and set each surface
      index_type i, j, n(3), m(3);
      point_type pt[3+1][3+1], pt_out;

      // create surface with specified control points
      pt[0][0] << -15, 0,  15;
      pt[1][0] <<  -5, 5,  15;
      pt[2][0] <<   5, 5,  15;
      pt[3][0] <<  15, 0,  15;
      pt[0][1] << -15, 5,   5;
      pt[1][1] <<  -5, 5,   5;
      pt[2][1] <<   5, 5,   5;
      pt[3][1] <<  15, 5,   5;
      pt[0][2] << -15, 5,  -5;
      pt[1][2] <<  -5, 5,  -5;
      pt[2][2] <<   5, 5,  -5;
      pt[3][2] <<  15, 5,  -5;
      pt[0][3] << -15, 0, -15;
      pt[1][3] <<  -5, 5, -15;
      pt[2][3] <<   5, 5, -15;
      pt[3][3] <<  15, 0, -15;
      s.resize(n, m);
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          s.set_control_point(pt[i][j], i, j);
        }
      }

      s.split_v(s1, s2, 0.5);  // this splits surface into lower and upper
      s1.split_u(s3, s4, 0.5); // this splits lower into first segment and last two
      err=ps1.set(s3, 0, 0); s_patches[0]=s3;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      s2.split_u(s5, s6, 0.5); // this splits upper into first segment and last two
      err=ps1.set(s5, 0, 1); s_patches[3]=s5;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      s4.split_u(s1, s2, 0.5); // this splits lower end into final two pieces
      err=ps1.set(s1, 1, 0); s_patches[1]=s1;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      err=ps1.set(s2, 2, 0); s_patches[2]=s2;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      s6.split_u(s1, s2, 0.5); // this splits the upper end into final two pieces
      err=ps1.set(s1, 1, 1); s_patches[4]=s1;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      err=ps1.set(s2, 2, 1); s_patches[5]=s2;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

      // test getting parameter max
      data_type umax, vmax;
      ps1.get_parameter_max(umax, vmax);
      TEST_ASSERT(tol.approximately_equal(umax, 3));
      TEST_ASSERT(tol.approximately_equal(vmax, 2));

      // test copy ctr
      piecewise_surface_type ps1c(ps1);
      TEST_ASSERT(ps1==ps1c);

      // test changing u0 and v0
      ps1c.set_u0(-1.5);
      TEST_ASSERT(ps1c.get_u0()==-1.5);
      ps1c.set_v0(-1.5);
      TEST_ASSERT(ps1c.get_v0()==-1.5);

      // test assignment operator
      ps2=ps1;
      TEST_ASSERT(ps2==ps1);
    }

    void bounding_box_test()
    {
      surface_type s, s1, s2, s3, s4, s5, s6, s_patches[6];
      piecewise_surface_type ps1, ps2;
      typename piecewise_surface_type::error_code err;

      // create 3x2 patches with unit spacing
      ps1.init_uv(3, 2);
      TEST_ASSERT(ps1.number_u_patches()==3);
      TEST_ASSERT(ps1.number_v_patches()==2);

      // create and set each surface
      index_type i, j, n(3), m(3);
      point_type pt[3+1][3+1], pt_out;

      // create surface with specified control points
      pt[0][0] << -15, 0,  15;
      pt[1][0] <<  -5, 5,  15;
      pt[2][0] <<   5, 5,  15;
      pt[3][0] <<  15, 0,  15;
      pt[0][1] << -15, 5,   5;
      pt[1][1] <<  -5, 5,   5;
      pt[2][1] <<   5, 5,   5;
      pt[3][1] <<  15, 5,   5;
      pt[0][2] << -15, 5,  -5;
      pt[1][2] <<  -5, 5,  -5;
      pt[2][2] <<   5, 5,  -5;
      pt[3][2] <<  15, 5,  -5;
      pt[0][3] << -15, 0, -15;
      pt[1][3] <<  -5, 5, -15;
      pt[2][3] <<   5, 5, -15;
      pt[3][3] <<  15, 0, -15;
      s.resize(n, m);
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          s.set_control_point(pt[i][j], i, j);
        }
      }

      s.split_v(s1, s2, 0.5);  // this splits surface into lower and upper
      s1.split_u(s3, s4, 0.5); // this splits lower into first segment and last two
      err=ps1.set(s3, 0, 0); s_patches[0]=s3;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      s2.split_u(s5, s6, 0.5); // this splits upper into first segment and last two
      err=ps1.set(s5, 0, 1); s_patches[3]=s5;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      s4.split_u(s1, s2, 0.5); // this splits lower end into final two pieces
      err=ps1.set(s1, 1, 0); s_patches[1]=s1;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      err=ps1.set(s2, 2, 0); s_patches[2]=s2;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      s6.split_u(s1, s2, 0.5); // this splits the upper end into final two pieces
      err=ps1.set(s1, 1, 1); s_patches[4]=s1;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      err=ps1.set(s2, 2, 1); s_patches[5]=s2;
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

      // test the bounding box
      typename piecewise_surface_type::bounding_box_type bb;
      point_type pmin_ref, pmax_ref;

      ps1.get_bounding_box(bb);
      pmin_ref << -15, 0, -15;
      pmax_ref << 15, 4.6875, 15;
      TEST_ASSERT(bb.get_min()==pmin_ref);
      TEST_ASSERT(bb.get_max()==pmax_ref);
    }

    void reverse_test()
    {

      piecewise_surface_type ps1, ps2;
      data_type u, v, umax, vmax;

      // create 3x2 patches with unit spacing
      ps1.init_uv(3, 2);

      // create piecewise surface
      {
        typename piecewise_surface_type::error_code err;
        surface_type s, s1, s2, s3, s4, s5, s6;
        index_type i, j, n(3), m(3);
        point_type pt[3+1][3+1], pt_out;

        // create surface with specified control points
        pt[0][0] << -15, 0,  15;
        pt[1][0] <<  -5, 5,  15;
        pt[2][0] <<   5, 5,  15;
        pt[3][0] <<  15, 0,  15;
        pt[0][1] << -15, 5,   5;
        pt[1][1] <<  -5, 5,   5;
        pt[2][1] <<   5, 5,   5;
        pt[3][1] <<  15, 5,   5;
        pt[0][2] << -15, 5,  -5;
        pt[1][2] <<  -5, 5,  -5;
        pt[2][2] <<   5, 5,  -5;
        pt[3][2] <<  15, 5,  -5;
        pt[0][3] << -15, 0, -15;
        pt[1][3] <<  -5, 5, -15;
        pt[2][3] <<   5, 5, -15;
        pt[3][3] <<  15, 0, -15;
        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }

        // create and set each surface
        s.split_v(s1, s2, 0.5);  // this splits surface into lower and upper
        s1.split_u(s3, s4, 0.5); // this splits lower into first segment and last two
        err=ps1.set(s3, 0, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s2.split_u(s5, s6, 0.5); // this splits upper into first segment and last two
        err=ps1.set(s5, 0, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s4.split_u(s1, s2, 0.5); // this splits lower end into final two pieces
        err=ps1.set(s1, 1, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s6.split_u(s1, s2, 0.5); // this splits the upper end into final two pieces
        err=ps1.set(s1, 1, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      }
      ps1.get_parameter_max(umax, vmax);

      // reverse u-direction
      ps2=ps1;
      ps2.reverse_u();

      // test point
      u=1.25;
      v=0.75;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(umax-u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), -ps2.f_u(umax-u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), ps2.f_v(umax-u, v)));

      // test another point
      u=2.75;
      v=1.25;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(umax-u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), -ps2.f_u(umax-u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), ps2.f_v(umax-u, v)));

      // reverse v-direction
      ps2=ps1;
      ps2.reverse_v();

      // test point
      u=1.25;
      v=0.75;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(u, vmax-v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), ps2.f_u(u, vmax-v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), -ps2.f_v(u, vmax-v)));

      // test another point
      u=2.75;
      v=1.25;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(u, vmax-v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), ps2.f_u(u, vmax-v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), -ps2.f_v(u, vmax-v)));
    }

    void swap_test()
    {
      piecewise_surface_type ps1, ps2;
      data_type u, v;

      // create 3x2 patches with unit spacing
      ps1.init_uv(3, 2);
      TEST_ASSERT(ps1.number_u_patches()==3);
      TEST_ASSERT(ps1.number_v_patches()==2);

      // create piecewise surface
      {
        typename piecewise_surface_type::error_code err;
        surface_type s, s1, s2, s3, s4, s5, s6;
        index_type i, j, n(3), m(3);
        point_type pt[3+1][3+1], pt_out;

        // create surface with specified control points
        pt[0][0] << -15, 0,  15;
        pt[1][0] <<  -5, 5,  15;
        pt[2][0] <<   5, 5,  15;
        pt[3][0] <<  15, 0,  15;
        pt[0][1] << -15, 5,   5;
        pt[1][1] <<  -5, 5,   5;
        pt[2][1] <<   5, 5,   5;
        pt[3][1] <<  15, 5,   5;
        pt[0][2] << -15, 5,  -5;
        pt[1][2] <<  -5, 5,  -5;
        pt[2][2] <<   5, 5,  -5;
        pt[3][2] <<  15, 5,  -5;
        pt[0][3] << -15, 0, -15;
        pt[1][3] <<  -5, 5, -15;
        pt[2][3] <<   5, 5, -15;
        pt[3][3] <<  15, 0, -15;
        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }

        // create and set each surface
        s.split_v(s1, s2, 0.5);  // this splits surface into lower and upper
        s1.split_u(s3, s4, 0.5); // this splits lower into first segment and last two
        err=ps1.set(s3, 0, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s2.split_u(s5, s6, 0.5); // this splits upper into first segment and last two
        err=ps1.set(s5, 0, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s4.split_u(s1, s2, 0.5); // this splits lower end into final two pieces
        err=ps1.set(s1, 1, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s6.split_u(s1, s2, 0.5); // this splits the upper end into final two pieces
        err=ps1.set(s1, 1, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      }

      // create second piecewise surface and swap coordinate directions
      ps2=ps1;
      ps2.swap_uv();

      TEST_ASSERT(ps2.number_u_patches()==2);
      TEST_ASSERT(ps2.number_v_patches()==3);

      // test point
      u=1.25;
      v=0.75;
      TEST_ASSERT(ps1.f(u, v)==ps2.f(v, u));
      TEST_ASSERT(ps1.f_u(u, v)==ps2.f_v(v, u));
      TEST_ASSERT(ps1.f_v(u, v)==ps2.f_u(v, u));

      // test another point
      u=2.75;
      v=1.25;
      TEST_ASSERT(ps1.f(u, v)==ps2.f(v, u));
      TEST_ASSERT(ps1.f_u(u, v)==ps2.f_v(v, u));
      TEST_ASSERT(ps1.f_v(u, v)==ps2.f_u(v, u));
    }

    void replace_test()
    {
      typename piecewise_surface_type::error_code err;
      piecewise_surface_type ps1, ps2;
      data_type u, v;

      // create 3x2 patches with unit spacing
      ps1.init_uv(3, 2);

      // create piecewise surface
      {
        surface_type s, s1, s2, s3, s4, s5, s6;
        index_type i, j, n(3), m(3);
        point_type pt[3+1][3+1], pt_out;

        // create surface with specified control points
        pt[0][0] << -15, 0,  15;
        pt[1][0] <<  -5, 5,  15;
        pt[2][0] <<   5, 5,  15;
        pt[3][0] <<  15, 0,  15;
        pt[0][1] << -15, 5,   5;
        pt[1][1] <<  -5, 5,   5;
        pt[2][1] <<   5, 5,   5;
        pt[3][1] <<  15, 5,   5;
        pt[0][2] << -15, 5,  -5;
        pt[1][2] <<  -5, 5,  -5;
        pt[2][2] <<   5, 5,  -5;
        pt[3][2] <<  15, 5,  -5;
        pt[0][3] << -15, 0, -15;
        pt[1][3] <<  -5, 5, -15;
        pt[2][3] <<   5, 5, -15;
        pt[3][3] <<  15, 0, -15;
        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }

        // create and set each surface
        s.split_v(s1, s2, 0.5);  // this splits surface into lower and upper
        s1.split_u(s3, s4, 0.5); // this splits lower into first segment and last two
        err=ps1.set(s3, 0, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s2.split_u(s5, s6, 0.5); // this splits upper into first segment and last two
        err=ps1.set(s5, 0, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s4.split_u(s1, s2, 0.5); // this splits lower end into final two pieces
        err=ps1.set(s1, 1, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s6.split_u(s1, s2, 0.5); // this splits the upper end into final two pieces
        err=ps1.set(s1, 1, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      }

      surface_type s;

      // test adding valid patch
      ps2=ps1;
      ps2.get(s, 2, 1);
      s.promote_u();
      s.promote_v();
      err=ps2.replace(s, 2, 1);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

      u=2.25;
      v=1.75;
      TEST_ASSERT(ps1.f(u, v)==ps2.f(u, v));
      TEST_ASSERT(ps1.f_u(u, v)==ps2.f_u(u, v));
      TEST_ASSERT(ps1.f_v(u, v)==ps2.f_v(u, v));

      // test adding invalid patch
      point_type nudge;
      nudge << 0, static_cast<data_type>(0.01), 0;
      ps2=ps1;
      ps2.get(s, 2, 1);
      s.set_control_point(s.get_control_point(0, 0)+nudge, 0, 0);
      err=ps2.replace(s, 2, 1);
      TEST_ASSERT(err==piecewise_surface_type::PATCH_NOT_CONNECTED);
    }

    void transformation_test()
    {
      piecewise_surface_type ps1, ps2;
      data_type u, v;
      typename piecewise_surface_type::error_code err;

      // create 3x2 patches with unit spacing
      ps1.init_uv(3, 2);

      // create piecewise surface
      {
        surface_type s, s1, s2, s3, s4, s5, s6;
        index_type i, j, n(3), m(3);
        point_type pt[3+1][3+1], pt_out;

        // create surface with specified control points
        pt[0][0] << -15, 0,  15;
        pt[1][0] <<  -5, 5,  15;
        pt[2][0] <<   5, 5,  15;
        pt[3][0] <<  15, 0,  15;
        pt[0][1] << -15, 5,   5;
        pt[1][1] <<  -5, 5,   5;
        pt[2][1] <<   5, 5,   5;
        pt[3][1] <<  15, 5,   5;
        pt[0][2] << -15, 5,  -5;
        pt[1][2] <<  -5, 5,  -5;
        pt[2][2] <<   5, 5,  -5;
        pt[3][2] <<  15, 5,  -5;
        pt[0][3] << -15, 0, -15;
        pt[1][3] <<  -5, 5, -15;
        pt[2][3] <<   5, 5, -15;
        pt[3][3] <<  15, 0, -15;
        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }

        // create and set each surface
        s.split_v(s1, s2, 0.5);  // this splits surface into lower and upper
        s1.split_u(s3, s4, 0.5); // this splits lower into first segment and last two
        err=ps1.set(s3, 0, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s2.split_u(s5, s6, 0.5); // this splits upper into first segment and last two
        err=ps1.set(s5, 0, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s4.split_u(s1, s2, 0.5); // this splits lower end into final two pieces
        err=ps1.set(s1, 1, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s6.split_u(s1, s2, 0.5); // this splits the upper end into final two pieces
        err=ps1.set(s1, 1, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      }

      // test translation
      {
        point_type trans;

        // set up translation vector and apply
        u=1.5;
        v=0.75;
        ps2=ps1;
        trans << 2, 1, 3;
        ps2.translate(trans);
        TEST_ASSERT(tol.approximately_equal(ps2.f(u, v), ps1.f(u, v)+trans));
      }

      // test rotation about origin
      {
        typename piecewise_surface_type::rotation_matrix_type rmat;

        // set up rotation and apply
        data_type one(1);
        ps2=ps1;
        rmat << std::cos(one), 0,  -std::sin(one),
                0,             one, 0,
                std::sin(one), 0,   std::cos(one);
        ps2.rotate(rmat);
        TEST_ASSERT(tol.approximately_equal(ps2.f(u, v), ps1.f(u, v)*rmat.transpose()));
      }

      // test rotation about point
      {
        typename piecewise_surface_type::rotation_matrix_type rmat;
        point_type rorig;

        // set up rotation and apply
        data_type one(1);
        ps2=ps1;
        rorig << 2, 1, 3;
        rmat << std::cos(one), 0,  -std::sin(one),
                0,             one, 0,
                std::sin(one), 0,   std::cos(one);
        ps2.rotate(rmat, rorig);
        TEST_ASSERT(tol.approximately_equal(ps2.f(u, v), rorig+(ps1.f(u, v)-rorig)*rmat.transpose()));
      }
    }

    void evaluation_test()
    {
      surface_type s_patches[6];
      piecewise_surface_type ps1;
      data_type u, v, up, vp;
      typename piecewise_surface_type::error_code err;

      // create 3x2 patches with unit spacing
      ps1.init_uv(3, 2);

      // create piecewise surface
      {
        surface_type s, s1, s2, s3, s4, s5, s6;
        index_type i, j, n(3), m(3);
        point_type pt[3+1][3+1], pt_out;

        // create surface with specified control points
        pt[0][0] << -15, 0,  15;
        pt[1][0] <<  -5, 5,  15;
        pt[2][0] <<   5, 5,  15;
        pt[3][0] <<  15, 0,  15;
        pt[0][1] << -15, 5,   5;
        pt[1][1] <<  -5, 5,   5;
        pt[2][1] <<   5, 5,   5;
        pt[3][1] <<  15, 5,   5;
        pt[0][2] << -15, 5,  -5;
        pt[1][2] <<  -5, 5,  -5;
        pt[2][2] <<   5, 5,  -5;
        pt[3][2] <<  15, 5,  -5;
        pt[0][3] << -15, 0, -15;
        pt[1][3] <<  -5, 5, -15;
        pt[2][3] <<   5, 5, -15;
        pt[3][3] <<  15, 0, -15;
        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }

        // create and set each surface
        s.split_v(s1, s2, 0.5);  // this splits surface into lower and upper
        s1.split_u(s3, s4, 0.5); // this splits lower into first segment and last two
        err=ps1.set(s3, 0, 0); s_patches[0]=s3;
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s2.split_u(s5, s6, 0.5); // this splits upper into first segment and last two
        err=ps1.set(s5, 0, 1); s_patches[3]=s5;
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s4.split_u(s1, s2, 0.5); // this splits lower end into final two pieces
        err=ps1.set(s1, 1, 0); s_patches[1]=s1;
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 0); s_patches[2]=s2;
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s6.split_u(s1, s2, 0.5); // this splits the upper end into final two pieces
        err=ps1.set(s1, 1, 1); s_patches[4]=s1;
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 1); s_patches[5]=s2;
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      }

      // test evaluation on first patch
      u=0.25;
      v=0.75;
      up=u;
      vp=v;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), s_patches[0].f(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), s_patches[0].f_u(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), s_patches[0].f_v(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_uu(u, v), s_patches[0].f_uu(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_uv(u, v), s_patches[0].f_uv(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_vv(u, v), s_patches[0].f_vv(up, vp)));

      // test evaluation on a lower patch
      u=2.25;
      v=0.75;
      up=u-2;
      vp=v;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), s_patches[2].f(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), s_patches[2].f_u(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), s_patches[2].f_v(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_uu(u, v), s_patches[2].f_uu(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_uv(u, v), s_patches[2].f_uv(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_vv(u, v), s_patches[2].f_vv(up, vp)));

      // test evaluation on an upper patch
      u=2.25;
      v=1.75;
      up=u-2;
      vp=v-1;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), s_patches[5].f(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), s_patches[5].f_u(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), s_patches[5].f_v(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_uu(u, v), s_patches[5].f_uu(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_uv(u, v), s_patches[5].f_uv(up, vp)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_vv(u, v), s_patches[5].f_vv(up, vp)));
    }

    void split_test()
    {
      piecewise_surface_type ps1, ps2;
      data_type u, v;

      // create 3x2 patches with unit spacing
      ps1.init_uv(3, 2);

      // create piecewise surface
      {
        typename piecewise_surface_type::error_code err;
        surface_type s, s1, s2, s3, s4, s5, s6;
        index_type i, j, n(3), m(3);
        point_type pt[3+1][3+1], pt_out;

        // create surface with specified control points
        pt[0][0] << -15, 0,  15;
        pt[1][0] <<  -5, 5,  15;
        pt[2][0] <<   5, 5,  15;
        pt[3][0] <<  15, 0,  15;
        pt[0][1] << -15, 5,   5;
        pt[1][1] <<  -5, 5,   5;
        pt[2][1] <<   5, 5,   5;
        pt[3][1] <<  15, 5,   5;
        pt[0][2] << -15, 5,  -5;
        pt[1][2] <<  -5, 5,  -5;
        pt[2][2] <<   5, 5,  -5;
        pt[3][2] <<  15, 5,  -5;
        pt[0][3] << -15, 0, -15;
        pt[1][3] <<  -5, 5, -15;
        pt[2][3] <<   5, 5, -15;
        pt[3][3] <<  15, 0, -15;
        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }

        // create and set each surface
        s.split_v(s1, s2, 0.5);  // this splits surface into lower and upper
        s1.split_u(s3, s4, 0.5); // this splits lower into first segment and last two
        err=ps1.set(s3, 0, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s2.split_u(s5, s6, 0.5); // this splits upper into first segment and last two
        err=ps1.set(s5, 0, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s4.split_u(s1, s2, 0.5); // this splits lower end into final two pieces
        err=ps1.set(s1, 1, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        s6.split_u(s1, s2, 0.5); // this splits the upper end into final two pieces
        err=ps1.set(s1, 1, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
        err=ps1.set(s2, 2, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      }

      // split u-direction
      ps2=ps1;
      ps2.split_u(1.5);

      // test evaluation before split
      u=0.25;
      v=1.5;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), ps2.f_u(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), ps2.f_v(u, v)));

      // test evaluation before split
      u=1.25;
      v=1.75;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), ps2.f_u(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), ps2.f_v(u, v)));

      // test evaluation after split
      u=1.75;
      v=0.75;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), ps2.f_u(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), ps2.f_v(u, v)));

      // test evaluation after split
      u=2.25;
      v=0.25;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), ps2.f_u(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), ps2.f_v(u, v)));

      // split v-direction
      ps2=ps1;
      ps2.split_v(1.5);

      // test evaluation before split
      u=0.25;
      v=0.5;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), ps2.f_u(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), ps2.f_v(u, v)));

      // test evaluation before split
      u=1.25;
      v=0.75;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), ps2.f_u(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), ps2.f_v(u, v)));

      // test evaluation after split
      u=1.75;
      v=1.75;
      TEST_ASSERT(tol.approximately_equal(ps1.f(u, v), ps2.f(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_u(u, v), ps2.f_u(u, v)));
      TEST_ASSERT(tol.approximately_equal(ps1.f_v(u, v), ps2.f_v(u, v)));
    }

    void to_cubic_test()
    {
      piecewise_surface_type ps1, ps2;
      data_type u, v;

      // create 3x2 patches with unit spacing
      ps1.init_uv(1,1);

      // create piecewise surface
      {
        typename piecewise_surface_type::error_code err;
        surface_type s;
        index_type i, j, n(5), m(4);
        point_type pt[5+1][4+1], pt_out;

        // create surface with specified control points
        pt[0][0] <<  0, 0, 15;
        pt[1][0] <<  2, 6, 15;
        pt[2][0] <<  3, 0, 15;
        pt[3][0] <<  5, 4, 15;
        pt[4][0] <<  7, 1, 15;
        pt[5][0] <<  9, 1, 15;
        pt[0][1] <<  0, 0, 11;
        pt[1][1] <<  2, 6, 11;
        pt[2][1] <<  3, 0, 11;
        pt[3][1] <<  5, 4, 11;
        pt[4][1] <<  7, 1, 11;
        pt[5][1] <<  9, 1, 11;
        pt[0][2] <<  0, 0,  3;
        pt[1][2] <<  2, 6,  3;
        pt[2][2] <<  3, 0,  3;
        pt[3][2] <<  5, 4,  3;
        pt[4][2] <<  7, 1,  3;
        pt[5][2] <<  9, 1,  3;
        pt[0][3] <<  0, 0,  0;
        pt[1][3] <<  2, 6,  0;
        pt[2][3] <<  3, 0,  0;
        pt[3][3] <<  5, 4,  0;
        pt[4][3] <<  7, 1,  0;
        pt[5][3] <<  9, 1,  0;
        pt[0][4] <<  0, 0, -5;
        pt[1][4] <<  2, 6, -5;
        pt[2][4] <<  3, 0, -5;
        pt[3][4] <<  5, 4, -5;
        pt[4][4] <<  7, 1, -5;
        pt[5][4] <<  9, 1, -5;

        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }

        err=ps1.set(s, 0, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

      }

      data_type ttol = static_cast<data_type>(1e-3);

      // to_cubic u-direction
      ps2=ps1;

      ps2.to_cubic_u(ttol);

      index_type mind, maxd;

      ps2.degree_u(mind, maxd);

      TEST_ASSERT(mind==3);
      TEST_ASSERT(maxd==3);

      // printf("nu: %d nv: %d\n", ps2.number_u_patches(), ps2.number_v_patches());

      u=0.25;
      v=0.5;
      TEST_ASSERT((ps1.f(u, v)-ps2.f(u, v)).norm() < ttol);

      u=0.25;
      v=0.75;
      TEST_ASSERT((ps1.f(u, v)-ps2.f(u, v)).norm() < ttol);

      u=0.75;
      v=0.75;
      TEST_ASSERT((ps1.f(u, v)-ps2.f(u, v)).norm() < ttol);

      u=0.25;
      v=0.25;
      TEST_ASSERT((ps1.f(u, v)-ps2.f(u, v)).norm() < ttol);

      // to_cubic v-direction
      ps2=ps1;

      ps2.to_cubic_v(ttol);

      ps2.degree_v(mind, maxd);

      TEST_ASSERT(mind==3);
      TEST_ASSERT(maxd==3);

      // printf("nu: %d nv: %d\n", ps2.number_u_patches(), ps2.number_v_patches());

      u=0.25;
      v=0.5;
      TEST_ASSERT((ps1.f(u, v)-ps2.f(u, v)).norm() < ttol);

      u=0.25;
      v=0.75;
      TEST_ASSERT((ps1.f(u, v)-ps2.f(u, v)).norm() < ttol);

      u=0.75;
      v=0.75;
      TEST_ASSERT((ps1.f(u, v)-ps2.f(u, v)).norm() < ttol);

      // to_cubic u and v-direction
      ps2=ps1;

      ps2.to_cubic(ttol);

      ps2.degree_u(mind, maxd);

      TEST_ASSERT(mind==3);
      TEST_ASSERT(maxd==3);

      ps2.degree_v(mind, maxd);

      TEST_ASSERT(mind==3);
      TEST_ASSERT(maxd==3);

      // printf("nu: %d nv: %d\n", ps2.number_u_patches(), ps2.number_v_patches());

      u=0.25;
      v=0.5;
      TEST_ASSERT((ps1.f(u, v)-ps2.f(u, v)).norm() < ttol);

      u=0.25;
      v=0.75;
      TEST_ASSERT((ps1.f(u, v)-ps2.f(u, v)).norm() < ttol);

      u=0.75;
      v=0.75;
      TEST_ASSERT((ps1.f(u, v)-ps2.f(u, v)).norm() < ttol);
    }

    void area_test()
    {
#if 0
      data_type eps(std::numeric_limits<data__>::epsilon());
      piecewise_surface_type c1;
      curve_type bc[3];
      data_type dt[3], len, bc_len[3], ref_len, t0, t1, temp0, temp1;
      typename curve_type::control_point_type cntrl1_in[4], cntrl2_in[5], cntrl3_in[3];
      typename piecewise_surface_type::error_code err;
      data_type tol(std::sqrt(eps));

      // create piecewise curve
      cntrl1_in[0] << 2.0, 2.0, 0.0;
      cntrl1_in[1] << 1.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 0.0, 0.0;
      cntrl1_in[3] << 4.0, 1.0, 0.0;
      dt[0]=0.5;
      bc[0].resize(3);
      for (index_type i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      eli::geom::curve::length(bc_len[0], bc[0], tol);
      cntrl2_in[0] << 4.0, 1.0, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 1.0, 0.0;
      cntrl2_in[3] << 6.0, 0.0, 0.0;
      cntrl2_in[4] << 6.5,-0.5, 0.0;
      dt[1]=2.0;
      bc[1].resize(4);
      for (index_type i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      eli::geom::curve::length(bc_len[1], bc[1], tol);
      cntrl3_in[0] << 6.5,-0.5, 0.0;
      cntrl3_in[1] << 6.0,-1.0, 0.0;
      cntrl3_in[2] << 5.5,-2.0, 0.0;
      dt[2]=1.5;
      bc[2].resize(2);
      for (index_type i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      eli::geom::curve::length(bc_len[2], bc[2], tol);
      err=c1.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);
      TEST_ASSERT(c1.number_segments()==3);

      // create two segment curve calc length of each segment to compare
      eli::geom::curve::length(len, c1, tol);
      ref_len=bc_len[0]+bc_len[1]+bc_len[2];
      TEST_ASSERT(len==ref_len);

      // choose part of first segment to calc length and compare
      t0=0.125;
      t1=0.375;
      eli::geom::curve::length(len, c1, t0, t1, tol);
      eli::geom::curve::length(ref_len, bc[0], 0.25, 0.75, tol);
      TEST_ASSERT(len==ref_len);

      // choose part of second segment to calc length and compare
      t0=0.25;
      t1=1.5;
      eli::geom::curve::length(len, c1, t0, t1, tol);
      eli::geom::curve::length(temp0, bc[0], 0.5, 1, tol);
      eli::geom::curve::length(temp1, bc[1], 0, 0.5, tol);
      ref_len=temp0+temp1;
      TEST_ASSERT(std::abs(len-ref_len)<2*tol);

      // choose part of third segment to calc length and compare
      t0=0.25;
      t1=3.25;
      eli::geom::curve::length(len, c1, t0, t1, tol);
      eli::geom::curve::length(temp0, bc[0], 0.5, 1, tol);
      eli::geom::curve::length(temp1, bc[2], 0, 0.5, tol);
      ref_len=temp0+bc_len[1]+temp1;
      TEST_ASSERT(std::abs(len-ref_len)<166*tol);
#endif
    }

    void get_curve_test()
    {
      const index_type n(3), m(3);
      surface_type s(n, m);
      point_type cp[3+1][3+1], pt_out, pt_ref;
      point_type pt, norm, u_contra, v_contra;
      data_type u, v, d;
      index_type i, j;
      piecewise_surface_type pws;
      piecewise_curve_type pwc;
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
      pws.split_u(0.3);
      pws.split_u(0.6);
      pws.split_u(0.7);

      pws.split_v(0.2);
      pws.split_v(0.412);
      pws.split_v(0.6);
      pws.split_v(0.895);

      // Extract u curve
      u=0;
      pws.get_uconst_curve(pwc, u);

      v=0;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.2;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.4;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.6;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.8;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      // Extract u curve
      u=0.1;
      pws.get_uconst_curve(pwc, u);

      v=0;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.2;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.4;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.6;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.8;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      // Extract u curve
      u=0.2;
      pws.get_uconst_curve(pwc, u);

      v=0;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.2;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.4;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.6;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.8;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      // Extract u curve
      u=0.5;
      pws.get_uconst_curve(pwc, u);

      v=0;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.2;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.4;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.6;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.8;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      // Extract u curve
      u=1;
      pws.get_uconst_curve(pwc, u);

      v=0;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.2;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.4;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.6;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=0.8;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);

      v=1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(v);
      TEST_ASSERT(pt_out==pt_ref);



      // Extract v curve
      v=0;
      pws.get_vconst_curve(pwc, v);

      u=0;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.2;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.4;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.6;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.8;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      // Extract v curve
      v=0.1;
      pws.get_vconst_curve(pwc, v);

      u=0;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=0.2;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.4;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=0.6;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=0.8;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      // Extract v curve
      v=0.2;
      pws.get_vconst_curve(pwc, v);

      u=0;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.2;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.4;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=0.6;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=0.8;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      // Extract v curve
      v=0.4;
      pws.get_vconst_curve(pwc, v);

      u=0;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=0.2;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=0.4;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=0.6;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.8;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      d=(pt_ref - pt_out).norm();
      TEST_ASSERT(d<std::numeric_limits<data_type>::epsilon()*30);

      u=1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      // Extract v curve
      v=1;
      pws.get_vconst_curve(pwc, v);

      u=0;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.2;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.4;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.6;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=0.8;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);

      u=1;
      pt_ref=pws.f(u, v);
      pt_out=pwc.f(u);
      TEST_ASSERT(pt_out==pt_ref);
    }

    void continuity_test()
    {
      // create
      {
        surface_type s;
        piecewise_surface_type ps;
        std::vector<data_type> du(3), dv(2);
        typename piecewise_surface_type::error_code err;

        // set the du and dv
        du[0]=0.5;
        du[1]=1.5;
        du[2]=0.5;
        dv[0]=1.25;
        dv[1]=0.5;
        // create 3x2 patches with unit spacing
        ps.init_uv(du.begin(), du.end(), dv.begin(), dv.end(), -1, -2);
        TEST_ASSERT(ps.number_u_patches()==3);
        TEST_ASSERT(ps.number_v_patches()==2);

        // create and set each surface
        index_type i, j, n, m;
        point_type pt[3+1][3+1], pt_out;

        // create patch 0,0
        n=1;
        m=3;
        pt[0][0] << -20, 0,  15;
        pt[0][1] << -20, 5,   5;
        pt[0][2] << -20, 5,  -5;
        pt[0][3] << -20, 0, -15;
        pt[1][0] << -15, 0,  15;
        pt[1][1] << -15, 5,   5;
        pt[1][2] << -15, 5,  -5;
        pt[1][3] << -15, 0, -15;
        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }
        err=ps.set(s, 0, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

        // create patch 1, 0
        n=3;
        m=3;
        pt[0][0] << -15, 0,  15;
        pt[1][0] <<  -5, 5,  15;
        pt[2][0] <<   5, 5,  15;
        pt[3][0] <<  15, 0,  15;
        pt[0][1] << -15, 5,   5;
        pt[1][1] <<  -5, 5,   5;
        pt[2][1] <<   5, 5,   5;
        pt[3][1] <<  15, 5,   5;
        pt[0][2] << -15, 5,  -5;
        pt[1][2] <<  -5, 5,  -5;
        pt[2][2] <<   5, 5,  -5;
        pt[3][2] <<  15, 5,  -5;
        pt[0][3] << -15, 0, -15;
        pt[1][3] <<  -5, 5, -15;
        pt[2][3] <<   5, 5, -15;
        pt[3][3] <<  15, 0, -15;
        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }
        err=ps.set(s, 1, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

        // create patch 2,0
        n=1;
        m=3;
        pt[0][0] <<  15, 0,  15;
        pt[0][1] <<  15, 5,   5;
        pt[0][2] <<  15, 5,  -5;
        pt[0][3] <<  15, 0, -15;
        pt[1][0] <<  20, 0,  15;
        pt[1][1] <<  20, 5,   5;
        pt[1][2] <<  20, 5,  -5;
        pt[1][3] <<  20, 0, -15;
        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }
        err=ps.set(s, 2, 0);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

        // create patch 0,1
        n=1;
        m=2;
        pt[0][0] << -20,  0, -15;
        pt[0][1] << -20, -3,  -7;
        pt[0][2] << -20, -7, -10;
        pt[1][0] << -15,  0, -15;
        pt[1][1] << -15, -3,  -7;
        pt[1][2] << -15, -7, -10;

        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }
        err=ps.set(s, 0, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

        // create patch 1,1
        n=3;
        m=2;
        pt[0][0] << -15,  0, -15;
        pt[0][1] << -15, -3,  -7;
        pt[0][2] << -15, -7, -10;
        pt[1][0] <<  -5,  5, -15;
        pt[1][1] <<  -5, -1,  -7;
        pt[1][2] <<  -5, -4, -10;
        pt[2][0] <<   5,  5, -15;
        pt[2][1] <<   5, -1,  -7;
        pt[2][2] <<   5, -4, -10;
        pt[3][0] <<  15,  0, -15;
        pt[3][1] <<  15, -3,  -7;
        pt[3][2] <<  15, -7, -10;

        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }
        err=ps.set(s, 1, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

        // create patch 2,1
        n=1;
        m=2;
        pt[0][0] <<  15,  0, -15;
        pt[0][1] <<  15, -3,  -7;
        pt[0][2] <<  15, -7, -10;
        pt[1][0] <<  20,  0, -15;
        pt[1][1] <<  20, -3,  -7;
        pt[1][2] <<  20, -7, -10;

        s.resize(n, m);
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=m; ++j)
          {
            s.set_control_point(pt[i][j], i, j);
          }
        }
        err=ps.set(s, 2, 1);
        TEST_ASSERT(err==piecewise_surface_type::NO_ERRORS);

        std::vector<data_type> ujoints(2*ps.number_u_patches()+1), vjoints(2*ps.number_v_patches()+1);
        std::vector<data_type> disc_ujoints, disc_ujoints_ref(2), disc_vjoints, disc_vjoints_ref(1);
        tolerance_type tol;

        // store joints
        ujoints[0]=ps.get_u0();
        ujoints[1]=ujoints[0]+du[0]/2;
        ujoints[2]=ujoints[0]+du[0];
        ujoints[3]=ujoints[2]+du[1]/2;
        ujoints[4]=ujoints[2]+du[1];
        ujoints[5]=ujoints[4]+du[2]/2;
        ujoints[6]=ujoints[4]+du[2];
        vjoints[0]=ps.get_v0();
        vjoints[1]=vjoints[0]+dv[0]/2;
        vjoints[2]=vjoints[0]+dv[0];
        vjoints[3]=vjoints[2]+dv[1]/2;
        vjoints[4]=vjoints[2]+dv[1];

        // store the discontinuous joints
        disc_ujoints_ref[0]=ujoints[2];
        disc_ujoints_ref[1]=ujoints[4];
        disc_vjoints_ref[0]=vjoints[2];

        // split curve at 3 locations to create continuous joints
        ps.split_u(ujoints[1]);
        ps.split_u(ujoints[3]);
        ps.split_u(ujoints[5]);
        ps.split_v(vjoints[1]);
        ps.split_v(vjoints[3]);
        TEST_ASSERT(ps.number_u_patches()==6);
        TEST_ASSERT(ps.number_v_patches()==4);

//        if (typeid(data_type)==typeid(float))
//        {
//          octave_print(1, ps);
//        }

        ps.find_interior_C0_edges(disc_ujoints, disc_vjoints);
        TEST_ASSERT(disc_ujoints.size()==disc_ujoints_ref.size());
        TEST_ASSERT(disc_vjoints.size()==disc_vjoints_ref.size());
        TEST_ASSERT(tol.approximately_equal(disc_ujoints[0], disc_ujoints_ref[0]));
        TEST_ASSERT(tol.approximately_equal(disc_ujoints[1], disc_ujoints_ref[1]));
        TEST_ASSERT(tol.approximately_equal(disc_vjoints[0], disc_vjoints_ref[0]));
      }
    }
};

#endif

