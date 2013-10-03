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

#ifndef bezier_surface_test_suite_hpp
#define bezier_surface_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/surface/bezier.hpp"
#include "eli/geom/surface/curvature.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class bezier_surface_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::bezier<data__, 3> bezier_type;
    typedef typename bezier_type::control_point_type control_point_type;
    typedef typename bezier_type::point_type point_type;
    typedef typename bezier_type::data_type data_type;
    typedef typename bezier_type::index_type index_type;
    typedef typename bezier_type::tolerance_type tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(bezier_surface_test_suite<float>::assignment_test);
      TEST_ADD(bezier_surface_test_suite<float>::bounding_box_test);
      TEST_ADD(bezier_surface_test_suite<float>::reverse_test);
      TEST_ADD(bezier_surface_test_suite<float>::swap_test);
      TEST_ADD(bezier_surface_test_suite<float>::transformation_test);
      TEST_ADD(bezier_surface_test_suite<float>::evaluation_test);
      TEST_ADD(bezier_surface_test_suite<float>::derivative_1_test);
      TEST_ADD(bezier_surface_test_suite<float>::derivative_2_test);
      TEST_ADD(bezier_surface_test_suite<float>::derivative_3_test);
      TEST_ADD(bezier_surface_test_suite<float>::curvature_test);
      TEST_ADD(bezier_surface_test_suite<float>::promotion_test);
      TEST_ADD(bezier_surface_test_suite<float>::promotion_to_test);
      TEST_ADD(bezier_surface_test_suite<float>::demotion_test);
      TEST_ADD(bezier_surface_test_suite<float>::split_test);
      TEST_ADD(bezier_surface_test_suite<float>::normal_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(bezier_surface_test_suite<double>::assignment_test);
      TEST_ADD(bezier_surface_test_suite<double>::bounding_box_test);
      TEST_ADD(bezier_surface_test_suite<double>::reverse_test);
      TEST_ADD(bezier_surface_test_suite<double>::swap_test);
      TEST_ADD(bezier_surface_test_suite<double>::transformation_test);
      TEST_ADD(bezier_surface_test_suite<double>::evaluation_test);
      TEST_ADD(bezier_surface_test_suite<double>::derivative_1_test);
      TEST_ADD(bezier_surface_test_suite<double>::derivative_2_test);
      TEST_ADD(bezier_surface_test_suite<double>::derivative_3_test);
      TEST_ADD(bezier_surface_test_suite<double>::curvature_test);
      TEST_ADD(bezier_surface_test_suite<double>::promotion_test);
      TEST_ADD(bezier_surface_test_suite<double>::promotion_to_test);
      TEST_ADD(bezier_surface_test_suite<double>::demotion_test);
      TEST_ADD(bezier_surface_test_suite<double>::split_test);
      TEST_ADD(bezier_surface_test_suite<double>::normal_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(bezier_surface_test_suite<long double>::assignment_test);
      TEST_ADD(bezier_surface_test_suite<long double>::bounding_box_test);
      TEST_ADD(bezier_surface_test_suite<long double>::reverse_test);
      TEST_ADD(bezier_surface_test_suite<long double>::swap_test);
      TEST_ADD(bezier_surface_test_suite<long double>::transformation_test);
      TEST_ADD(bezier_surface_test_suite<long double>::evaluation_test);
      TEST_ADD(bezier_surface_test_suite<long double>::derivative_1_test);
      TEST_ADD(bezier_surface_test_suite<long double>::derivative_2_test);
      TEST_ADD(bezier_surface_test_suite<long double>::derivative_3_test);
      TEST_ADD(bezier_surface_test_suite<long double>::curvature_test);
      TEST_ADD(bezier_surface_test_suite<long double>::promotion_test);
      TEST_ADD(bezier_surface_test_suite<long double>::promotion_to_test);
      TEST_ADD(bezier_surface_test_suite<long double>::demotion_test);
      TEST_ADD(bezier_surface_test_suite<long double>::split_test);
      TEST_ADD(bezier_surface_test_suite<long double>::normal_test);
    }
#ifdef ELI_USING_QD
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(bezier_surface_test_suite<dd_real>::assignment_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::bounding_box_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::reverse_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::swap_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::transformation_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::evaluation_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::derivative_1_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::derivative_2_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::derivative_3_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::curvature_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::promotion_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::promotion_to_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::demotion_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::split_test);
      TEST_ADD(bezier_surface_test_suite<dd_real>::normal_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(bezier_surface_test_suite<qd_real>::assignment_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::bounding_box_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::reverse_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::swap_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::transformation_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::evaluation_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::derivative_1_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::derivative_2_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::derivative_3_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::curvature_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::promotion_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::promotion_to_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::demotion_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::split_test);
      TEST_ADD(bezier_surface_test_suite<qd_real>::normal_test);
    }
#endif

  public:
    bezier_surface_test_suite()
    {
      AddTests(data__());
    }
    ~bezier_surface_test_suite()
    {
    }

  private:
    void octave_print(int figno, /*const point_type pts[][], */const bezier_type &bez) const
    {
      index_type i, j;

      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << "cp_x=[";
      for (i=0; i<=bez.degree_u(); ++i)
      {
        std::cout << bez.get_control_point(i, 0).x();
        for (j=1; j<bez.degree_v(); ++j)
        {
          std::cout << ", " << bez.get_control_point(i, j).x();
        }
        j=bez.degree_v();
        std::cout << ", " << bez.get_control_point(i, j).x();
        if (i<bez.degree_u())
          std::cout << "; ";
      }
      std::cout << "];" << std::endl;

      std::cout << "cp_y=[";
      for (i=0; i<=bez.degree_u(); ++i)
      {
        std::cout << bez.get_control_point(i, 0).y();
        for (j=1; j<bez.degree_v(); ++j)
        {
          std::cout << ", " << bez.get_control_point(i, j).y();
        }
        j=bez.degree_v();
        std::cout << ", " << bez.get_control_point(i, j).y();
        if (i<bez.degree_u())
          std::cout << "; ";
      }
      std::cout << "];" << std::endl;

      std::cout << "cp_z=[";
      for (i=0; i<=bez.degree_u(); ++i)
      {
        std::cout << bez.get_control_point(i, 0).z();
        for (j=1; j<bez.degree_v(); ++j)
        {
          std::cout << ", " << bez.get_control_point(i, j).z();
        }
        j=bez.degree_v();
        std::cout << ", " << bez.get_control_point(i, j).z();
        if (i<bez.degree_u())
          std::cout << "; ";
      }
      std::cout << "];" << std::endl;

      // initialize the u & v parameters
      std::vector<data__> u(51), v(51);
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        u[i]=static_cast<data__>(i)/(u.size()-1);
      }
      for (j=0; j<static_cast<index_type>(v.size()); ++j)
      {
        v[j]=static_cast<data__>(j)/(v.size()-1);
      }

      // set the surface points
      std::cout << "surf_x=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << bez.f(u[i], v[0]).x();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << bez.f(u[i], v[j]).x();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << bez.f(u[i], v[j]).x();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_y=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << bez.f(u[i], v[0]).y();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << bez.f(u[i], v[j]).y();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << bez.f(u[i], v[j]).y();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_z=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << bez.f(u[i], v[0]).z();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << bez.f(u[i], v[j]).z();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << bez.f(u[i], v[j]).z();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "mesh(surf_x, surf_y, surf_z, zeros(size(surf_z)), 'EdgeColor', [0 0 0]);" << std::endl;
      std::cout << "hold on;" << std::endl;
      std::cout << "plot3(cp_x, cp_y, cp_z, '-ok', 'MarkerFaceColor', [0 0 0]);" << std::endl;
      std::cout << "plot3(cp_x', cp_y', cp_z', '-k');" << std::endl;
      std::cout << "hold off;" << std::endl;
    }

    void create_circle(std::vector<point_type> &/*pts*/)
    {
#if 0
      // NOTE: This will not create a closed circle, the last point will be
      //       the point just prior to the 2*pi point
      size_t n=pts.size();
      for (size_t i=0; i<n; ++i)
      {
        data__ theta(eli::constants::math<data__>::two_pi()*static_cast<data__>(i)/n);
        pts[i](0)=std::cos(theta);
        pts[i](1)=std::sin(theta);
        pts[i](2)=0;
      }
#endif
    }

    void assignment_test()
    {
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

      // create surface with specified dimensions and set control points
      bezier_type bez1(n, m);

      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          bez1.set_control_point(pt[i][j], i, j);
        }
      }
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          TEST_ASSERT(bez1.get_control_point(i, j) == pt[i][j]);
        }
      }

      // test order
      TEST_ASSERT(bez1.degree_u()==n);
      TEST_ASSERT(bez1.degree_v()==m);

      // test dimension
      TEST_ASSERT(bez1.dimension()==3);

      // test copy ctr
      bezier_type bez2(bez1);

      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          TEST_ASSERT(bez2.get_control_point(i, j) == bez1.get_control_point(i, j));
        }
      }

      // test equivalence operator
      TEST_ASSERT(bez2==bez1);

      // test assignment operator
      bezier_type bez3;

      bez3=bez1;
      TEST_ASSERT(bez3==bez1);

      // create surface then resize to needed dimensions
      bezier_type bez4;

      bez4.resize(n, m);
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          bez4.set_control_point(pt[i][j], i, j);
        }
      }
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          TEST_ASSERT(bez4.get_control_point(i, j) == pt[i][j]);
        }
      }
    }

    void bounding_box_test()
    {
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

      // create surface with specified dimensions and set control points
      bezier_type bez1(n, m);

      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          bez1.set_control_point(pt[i][j], i, j);
        }
      }
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          TEST_ASSERT(bez1.get_control_point(i, j) == pt[i][j]);
        }
      }

      // test bounding box
      typename bezier_type::bounding_box_type bb;
      point_type pmin_ref, pmax_ref;

      bez1.get_bounding_box(bb);
      pmin_ref << -15, 0, -15;
      pmax_ref <<  15, 5,  15;
      TEST_ASSERT(bb.get_min()==pmin_ref);
      TEST_ASSERT(bb.get_max()==pmax_ref);
    }

    void reverse_test()
    {
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

      // create surface with specified dimensions and set control points
      bezier_type bez1(n, m), bez2;

      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          bez1.set_control_point(pt[i][j], i, j);
        }
      }

      // reverse u-direction
      bez2=bez1;
      bez2.reverse_u();
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          TEST_ASSERT(bez2.get_control_point(i, j) == bez1.get_control_point(n-i, j));
        }
      }

      // reverse v-direction
      bez2=bez1;
      bez2.reverse_v();
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          TEST_ASSERT(bez2.get_control_point(i, j) == bez1.get_control_point(i, m-j));
        }
      }
    }

    void swap_test()
    {
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

      // create surface with specified dimensions and set control points
      bezier_type bez1(n, m), bez2;

      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          bez1.set_control_point(pt[i][j], i, j);
        }
      }

      // reverse swap u- & v-directions
      bez2=bez1;
      bez2.swap_uv();
      TEST_ASSERT(bez1.degree_u()==bez2.degree_v());
      TEST_ASSERT(bez1.degree_v()==bez2.degree_u());
      for (i=0; i<=n; ++i)
      {
        for (j=0; j<=m; ++j)
        {
          TEST_ASSERT(bez1.get_control_point(i, j) == bez2.get_control_point(j, i));
        }
      }
    }

    void transformation_test()
    {
      index_type n(3), m(3);
      point_type pt[3+1][3+1], pt_out, pt_ref;
      data_type u, v;

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

      // create surface with specified dimensions and set control points
      bezier_type bez(n, m), bez2;

      for (index_type i=0; i<=n; ++i)
      {
        for (index_type j=0; j<=m; ++j)
        {
          bez.set_control_point(pt[i][j], i, j);
        }
      }

      // test translation
      {
        point_type trans;

        // set up translation vector and apply
        bez2=bez;
        trans << 2, 1, 3;
        bez2.translate(trans);

        // test evaluation at corners
        u=0; v=0;
        pt_out=bez2.f(u, v);
        pt_ref=trans+pt[0][0];
        TEST_ASSERT(pt_out==pt_ref);
        u=1; v=0;
        pt_out=bez2.f(u, v);
        pt_ref=trans+pt[n][0];
        TEST_ASSERT(pt_out==pt_ref);
        u=0; v=1;
        pt_out=bez2.f(u, v);
        pt_ref=trans+pt[0][m];
        TEST_ASSERT(pt_out==pt_ref);
        u=1; v=1;
        pt_out=bez2.f(u, v);
        pt_ref=trans+pt[n][m];
        TEST_ASSERT(pt_out==pt_ref);

        // test evaluation at interior point u=v=1/2
        u=0.5; v=0.5;
        pt_ref << 0, 4.6875 , 0;
        pt_out=bez2.f(u, v);
        pt_ref=trans+pt_ref;
        TEST_ASSERT(pt_out==pt_ref);

        // test evaluation at interior point u=1/4, v=3/4
        u=0.25; v=0.75;
        pt_ref << -7.5, 4.04296875, -7.5;
        pt_out=bez2.f(u, v);
        pt_ref=trans+pt_ref;
        TEST_ASSERT(pt_out==pt_ref);
      }

      // test rotation about origin
      {
        typename bezier_type::rotation_matrix_type rmat;

        // set up rotation and apply
        bez2=bez;
        rmat << cos(1), 0, -sin(1),
                     0, 1,       0,
                sin(1), 0,  cos(1);
        bez2.rotate(rmat);

        // test evaluation at corners
        u=0; v=0;
        pt_out=bez2.f(u, v);
        pt_ref=pt[0][0]*rmat.transpose();
        TEST_ASSERT(pt_out==pt_ref);
        u=1; v=0;
        pt_out=bez2.f(u, v);
        pt_ref=pt[n][0]*rmat.transpose();
        TEST_ASSERT(pt_out==pt_ref);
        u=0; v=1;
        pt_out=bez2.f(u, v);
        pt_ref=pt[0][m]*rmat.transpose();
        TEST_ASSERT(pt_out==pt_ref);
        u=1; v=1;
        pt_out=bez2.f(u, v);
        pt_ref=pt[n][m]*rmat.transpose();
        TEST_ASSERT(pt_out==pt_ref);

        // test evaluation at interior point u=v=1/2
        u=0.5; v=0.5;
        pt_ref << 0, 4.6875 , 0;
        pt_out=bez2.f(u, v);
        pt_ref=pt_ref*rmat.transpose();
        TEST_ASSERT(pt_out==pt_ref);

        // test evaluation at interior point u=1/4, v=3/4
        u=0.25; v=0.75;
        pt_ref << -7.5, 4.04296875, -7.5;
        pt_out=bez2.f(u, v);
        pt_ref=pt_ref*rmat.transpose();
        TEST_ASSERT(tol.approximately_equal(pt_out,pt_ref));
      }

      // test rotation about point
      {
        point_type rorig;
        typename bezier_type::rotation_matrix_type rmat;

        // set up rotation and apply
        bez2=bez;
        rorig << 2, 1, 3;
        rmat << cos(1), 0, -sin(1),
                     0, 1,       0,
                sin(1), 0,  cos(1);
        bez2.rotate(rmat, rorig);

        // test evaluation at corners
        u=0; v=0;
        pt_out=bez2.f(u, v);
        pt_ref=rorig+(pt[0][0]-rorig)*rmat.transpose();
        TEST_ASSERT(pt_out==pt_ref);
        u=1; v=0;
        pt_out=bez2.f(u, v);
        pt_ref=rorig+(pt[n][0]-rorig)*rmat.transpose();
        TEST_ASSERT(pt_out==pt_ref);
        u=0; v=1;
        pt_out=bez2.f(u, v);
        pt_ref=rorig+(pt[0][m]-rorig)*rmat.transpose();
        TEST_ASSERT(pt_out==pt_ref);
        u=1; v=1;
        pt_out=bez2.f(u, v);
        pt_ref=rorig+(pt[n][m]-rorig)*rmat.transpose();
        TEST_ASSERT(pt_out==pt_ref);

        // test evaluation at interior point u=v=1/2
        u=0.5; v=0.5;
        pt_ref << 0, 4.6875 , 0;
        pt_out=bez2.f(u, v);
        pt_ref=rorig+(pt_ref-rorig)*rmat.transpose();
        TEST_ASSERT(tol.approximately_equal(pt_out,pt_ref));

        // test evaluation at interior point u=1/4, v=3/4
        u=0.25; v=0.75;
        pt_ref << -7.5, 4.04296875, -7.5;
        pt_out=bez2.f(u, v);
        pt_ref=rorig+(pt_ref-rorig)*rmat.transpose();
        TEST_ASSERT(tol.approximately_equal(pt_out,pt_ref));
      }
    }

    void evaluation_test()
    {
      index_type n(3), m(3);
      point_type pt[3+1][3+1], pt_out, pt_ref;
      data_type u, v;

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

      // create surface with specified dimensions and set control points
      bezier_type bez(n, m);

      for (index_type i=0; i<=n; ++i)
      {
        for (index_type j=0; j<=m; ++j)
        {
          bez.set_control_point(pt[i][j], i, j);
        }
      }

      // test evaluation at corners
      u=0; v=0;
      pt_out=bez.f(u, v);
      TEST_ASSERT(pt_out==pt[0][0]);
      u=1; v=0;
      pt_out=bez.f(u, v);
      TEST_ASSERT(pt_out==pt[n][0]);
      u=0; v=1;
      pt_out=bez.f(u, v);
      TEST_ASSERT(pt_out==pt[0][m]);
      u=1; v=1;
      pt_out=bez.f(u, v);
      TEST_ASSERT(pt_out==pt[n][m]);

      // test evaluation at interior point u=v=1/2
      u=0.5; v=0.5;
      pt_ref << 0, 4.6875 , 0;
      pt_out=bez.f(u, v);
      TEST_ASSERT(pt_out==pt_ref);

      // test evaluation at interior point u=1/4, v=3/4
      u=0.25; v=0.75;
      pt_ref << -7.5, 4.04296875, -7.5;
      pt_out=bez.f(u, v);
      TEST_ASSERT(pt_out==pt_ref);
    }

    void derivative_1_test()
    {
      index_type n(3), m(3);
      point_type pt[3+1][3+1], pt_out, pt_ref;
      data_type u, v;

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

      // create surface with specified dimensions and set control points
      bezier_type bez(n, m);

      for (index_type i=0; i<=n; ++i)
      {
        for (index_type j=0; j<=m; ++j)
        {
          bez.set_control_point(pt[i][j], i, j);
        }
      }

      // test evaluation at corners
      u=0; v=0;
      pt_out=bez.f_u(u, v);
      TEST_ASSERT(pt_out==n*(pt[1][0]-pt[0][0]));
      pt_out=bez.f_v(u, v);
      TEST_ASSERT(pt_out==m*(pt[0][1]-pt[0][0]));
      u=1; v=0;
      pt_out=bez.f_u(u, v);
      TEST_ASSERT(pt_out==n*(pt[3][0]-pt[2][0]));
      pt_out=bez.f_v(u, v);
      TEST_ASSERT(pt_out==m*(pt[3][1]-pt[3][0]));
      u=0; v=1;
      pt_out=bez.f_u(u, v);
      TEST_ASSERT(pt_out==n*(pt[1][3]-pt[0][3]));
      pt_out=bez.f_v(u, v);
      TEST_ASSERT(pt_out==m*(pt[0][3]-pt[0][2]));
      u=1; v=1;
      pt_out=bez.f_u(u, v);
      TEST_ASSERT(pt_out==n*(pt[3][3]-pt[2][3]));
      pt_out=bez.f_v(u, v);
      TEST_ASSERT(pt_out==m*(pt[3][3]-pt[3][2]));

      // test evaluation at interior point u=v=1/2
      u=0.5; v=0.5;
      pt_ref << 30, 0, 0;
      pt_out=bez.f_u(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, 0, -30;
      pt_out=bez.f_v(u, v);
      TEST_ASSERT(pt_out==pt_ref);

      // test evaluation at interior point u=1/4, v=3/4
      u=0.25; v=0.75;
      pt_ref << 30, 3.28125, 0;
      pt_out=bez.f_u(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, -3.28125, -30;
      pt_out=bez.f_v(u, v);
      TEST_ASSERT(pt_out==pt_ref);
    }

    void derivative_2_test()
    {
      index_type n(3), m(3);
      point_type pt[3+1][3+1], pt_out, pt_ref;
      data_type u, v;

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

      // create surface with specified dimensions and set control points
      bezier_type bez(n, m);

      for (index_type i=0; i<=n; ++i)
      {
        for (index_type j=0; j<=m; ++j)
        {
          bez.set_control_point(pt[i][j], i, j);
        }
      }

      // test evaluation at corners
      u=0; v=0;
      pt_out=bez.f_uu(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*(pt[2][0]-2*pt[1][0]+pt[0][0]));
      pt_out=bez.f_uv(u, v);
      TEST_ASSERT(pt_out==n*m*(pt[1][1]-pt[1][0]-pt[0][1]+pt[0][0]));
      pt_out=bez.f_vv(u, v);
      TEST_ASSERT(pt_out==m*(m-1)*(pt[0][2]-2*pt[0][1]+pt[0][0]));
      u=1; v=0;
      pt_out=bez.f_uu(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*(pt[3][0]-2*pt[2][0]+pt[1][0]));
      pt_out=bez.f_uv(u, v);
      TEST_ASSERT(pt_out==n*m*(pt[3][1]-pt[3][0]-pt[2][1]+pt[2][0]));
      pt_out=bez.f_vv(u, v);
      TEST_ASSERT(pt_out==m*(m-1)*(pt[3][2]-2*pt[3][1]+pt[3][0]));
      u=0; v=1;
      pt_out=bez.f_uu(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*(pt[2][3]-2*pt[1][3]+pt[0][3]));
      pt_out=bez.f_uv(u, v);
      TEST_ASSERT(pt_out==n*m*(pt[1][3]-pt[1][2]-pt[0][3]+pt[0][2]));
      pt_out=bez.f_vv(u, v);
      TEST_ASSERT(pt_out==m*(m-1)*(pt[0][3]-2*pt[0][2]+pt[0][1]));
      u=1; v=1;
      pt_out=bez.f_uu(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*(pt[3][3]-2*pt[2][3]+pt[1][3]));
      pt_out=bez.f_uv(u, v);
      TEST_ASSERT(pt_out==n*m*(pt[3][3]-pt[3][2]-pt[2][3]+pt[2][2]));
      pt_out=bez.f_vv(u, v);
      TEST_ASSERT(pt_out==m*(m-1)*(pt[3][3]-2*pt[3][2]+pt[3][1]));

      // test evaluation at interior point u=v=1/2
      u=0.5; v=0.5;
      pt_ref << 0, -7.5, 0;
      pt_out=bez.f_uu(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, 0, 0;
      pt_out=bez.f_uv(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, -7.5, 0;
      pt_out=bez.f_vv(u, v);
      TEST_ASSERT(pt_out==pt_ref);

      // test evaluation at interior point u=1/4, v=3/4
      u=0.25; v=0.75;
      pt_ref << 0, -13.125, 0;
      pt_out=bez.f_uu(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, 11.25, 0;
      pt_out=bez.f_uv(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, -13.125, 0;
      pt_out=bez.f_vv(u, v);
      TEST_ASSERT(pt_out==pt_ref);
    }

    void derivative_3_test()
    {
      index_type n(3), m(3);
      point_type pt[3+1][3+1], pt_out, pt_ref;
      data_type u, v;

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

      // create surface with specified dimensions and set control points
      bezier_type bez(n, m);

      for (index_type i=0; i<=n; ++i)
      {
        for (index_type j=0; j<=m; ++j)
        {
          bez.set_control_point(pt[i][j], i, j);
        }
      }

      // test evaluation at corners
      u=0; v=0;
      pt_out=bez.f_uuu(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*(n-2)*(pt[3][0]-3*pt[2][0]+3*pt[1][0]-pt[0][0]));
      pt_out=bez.f_uuv(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*m*((pt[2][1]-pt[2][0])-2*(pt[1][1]-pt[1][0])+(pt[0][1]-pt[0][0])));
      pt_out=bez.f_uvv(u, v);
      TEST_ASSERT(pt_out==n*m*(m-1)*((pt[1][2]-pt[0][2])-2*(pt[1][1]-pt[0][1])+(pt[1][0]-pt[0][0])));
      pt_out=bez.f_vvv(u, v);
      TEST_ASSERT(pt_out==m*(m-1)*(m-2)*(pt[0][3]-3*pt[0][2]+3*pt[0][1]-pt[0][0]));
      u=1; v=0;
      pt_out=bez.f_uuu(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*(n-2)*(pt[3][0]-3*pt[2][0]+3*pt[1][0]-pt[0][0]));
      pt_out=bez.f_uuv(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*m*((pt[3][1]-pt[3][0])-2*(pt[2][1]-pt[2][0])+(pt[1][1]-pt[1][0])));
      pt_out=bez.f_uvv(u, v);
      TEST_ASSERT(pt_out==n*m*(m-1)*((pt[3][2]-pt[2][2])-2*(pt[3][1]-pt[2][1])+(pt[3][0]-pt[2][0])));
      pt_out=bez.f_vvv(u, v);
      TEST_ASSERT(pt_out==m*(m-1)*(m-2)*(pt[3][3]-3*pt[3][2]+3*pt[3][1]-pt[3][0]));
      u=0; v=1;
      pt_out=bez.f_uuu(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*(n-2)*(pt[3][3]-3*pt[2][3]+3*pt[1][3]-pt[0][3]));
      pt_out=bez.f_uuv(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*m*((pt[2][3]-pt[2][2])-2*(pt[1][3]-pt[1][2])+(pt[0][3]-pt[0][2])));
      pt_out=bez.f_uvv(u, v);
      TEST_ASSERT(pt_out==n*m*(m-1)*((pt[1][3]-pt[0][3])-2*(pt[1][2]-pt[0][2])+(pt[1][1]-pt[0][1])));
      pt_out=bez.f_vvv(u, v);
      TEST_ASSERT(pt_out==m*(m-1)*(m-2)*(pt[0][3]-3*pt[0][2]+3*pt[0][1]-pt[0][0]));
      u=1; v=1;
      pt_out=bez.f_uuu(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*(n-2)*(pt[3][3]-3*pt[2][3]+3*pt[1][3]-pt[0][3]));
      pt_out=bez.f_uuv(u, v);
      TEST_ASSERT(pt_out==n*(n-1)*m*((pt[3][3]-pt[3][2])-2*(pt[2][3]-pt[2][2])+(pt[1][3]-pt[1][2])));
      pt_out=bez.f_uvv(u, v);
      TEST_ASSERT(pt_out==n*m*(m-1)*((pt[3][3]-pt[2][3])-2*(pt[3][2]-pt[2][2])+(pt[3][1]-pt[2][1])));
      pt_out=bez.f_vvv(u, v);
      TEST_ASSERT(pt_out==m*(m-1)*(m-2)*(pt[3][3]-3*pt[3][2]+3*pt[3][1]-pt[3][0]));

      // test evaluation at interior point u=v=1/2
      u=0.5; v=0.5;
      pt_ref << 0, 0, 0;
      pt_out=bez.f_uuu(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, 0, 0;
      pt_out=bez.f_uuv(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, 0, 0;
      pt_out=bez.f_uvv(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, 0, 0;
      pt_out=bez.f_vvv(u, v);
      TEST_ASSERT(pt_out==pt_ref);

      // test evaluation at interior point u=1/4, v=3/4
      u=0.25; v=0.75;
      pt_ref << 0, 0, 0;
      pt_out=bez.f_uuu(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, -45, 0;
      pt_out=bez.f_uuv(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, 45, 0;
      pt_out=bez.f_uvv(u, v);
      TEST_ASSERT(pt_out==pt_ref);
      pt_ref << 0, 0, 0;
      pt_out=bez.f_vvv(u, v);
      TEST_ASSERT(pt_out==pt_ref);
    }

    void curvature_test()
    {
      index_type n(3), m(3);
      point_type pt[3+1][3+1], pt_out, pt_ref;
      data_type u, v;

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

      // create surface with specified dimensions and set control points
      bezier_type bez(n, m);

      for (index_type i=0; i<=n; ++i)
      {
        for (index_type j=0; j<=m; ++j)
        {
          bez.set_control_point(pt[i][j], i, j);
        }
      }

      // test mean curvature
      {
        data_type curv_out, curv_ref;

        // test evaluation at corners
        u=0; v=0;
        curv_ref=-0.015876322406928;
        eli::geom::surface::mean_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-6);
        u=1; v=0;
        curv_ref=-0.015876322406928;
        eli::geom::surface::mean_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-6);
        u=0; v=1;
        curv_ref=-0.015876322406928;
        eli::geom::surface::mean_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-6);
        u=1; v=1;
        curv_ref=-0.015876322406928;
        eli::geom::surface::mean_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-6);

        // test at interior point u=v=1/2
        u=0.5; v=0.5;
        curv_ref=-0.008333333333333;
        eli::geom::surface::mean_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-6);

        // test at interior point u=1/4 & v=3/4
        u=0.25; v=0.75;
        curv_ref=-0.014099238414151;
        eli::geom::surface::mean_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-6);
      }

      // test Gaussian curvature
      {
        data_type curv_out, curv_ref;

        // test evaluation at corners
        u=0; v=0;
        curv_ref=-0.00061728395;
        eli::geom::surface::gaussian_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-10);
        u=1; v=0;
        curv_ref=-0.00061728395;
        eli::geom::surface::gaussian_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-10);
        u=0; v=1;
        curv_ref=-0.00061728395;
        eli::geom::surface::gaussian_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-10);
        u=1; v=1;
        curv_ref=-0.00061728395;
        eli::geom::surface::gaussian_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-10);

        // test at interior point u=v=1/2
        u=0.5; v=0.5;
        curv_ref=6.944444444444444e-5;
        eli::geom::surface::gaussian_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-10);

        // test at interior point u=1/4 & v=3/4
        u=0.25; v=0.75;
        curv_ref=5.381754978392456e-05;
        eli::geom::surface::gaussian_curvature(curv_out, bez, u, v);
        TEST_ASSERT(std::abs(curv_ref-curv_out)<1e-10);
      }

      // test principal curvatures without directions
      {
        data_type kmax_out, kmin_out, H_out, K_out;

        // test evaluation at corners
        u=0; v=0;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());
        u=1; v=0;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());
        u=0; v=1;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());
        u=1; v=1;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());

        // test at interior point u=v=1/2
        u=0.5; v=0.5;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());

        // test at interior point u=1/4 & v=3/4
        u=0.25; v=0.75;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());
      }

      // test principal curvatures with directions
      {
        data_type kmax_out, kmin_out, H_out, K_out;
        point_type kmax_dir_out, kmin_dir_out, n_out, kmax_dir_ref, kmin_dir_ref, n_ref;

        // test evaluation at corners
        u=0; v=0;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, kmax_dir_out, kmin_dir_out, n_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        kmax_dir_ref <<  0.70710678,          0,  0.70710678;
        kmin_dir_ref <<  0.57735027, 0.57735027, -0.57735027;
        n_ref        << -0.40824829, 0.81649658,  0.40824829;
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT((kmax_dir_out - kmax_dir_ref).norm()<1e-7);
        TEST_ASSERT((kmin_dir_out - kmin_dir_ref).norm()<1e-7);
        TEST_ASSERT((n_out - n_ref).norm()<1e-7);
        TEST_ASSERT(kmax_dir_out.dot(kmin_dir_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmax_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmin_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());
        u=1; v=0;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, kmax_dir_out, kmin_dir_out, n_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        kmax_dir_ref <<  0.70710678,          0, -0.70710678;
        kmin_dir_ref << -0.57735027, 0.57735027, -0.57735027;
        n_ref        <<  0.40824829, 0.81649658,  0.40824829;
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT((kmax_dir_out - kmax_dir_ref).norm()<1e-7);
        TEST_ASSERT((kmin_dir_out - kmin_dir_ref).norm()<1e-7);
        TEST_ASSERT((n_out - n_ref).norm()<1e-7);
        TEST_ASSERT(kmax_dir_out.dot(kmin_dir_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmax_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmin_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());
        u=0; v=1;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, kmax_dir_out, kmin_dir_out, n_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        kmax_dir_ref <<  0.70710678,          0, -0.70710678;
        kmin_dir_ref << -0.57735027,-0.57735027, -0.57735027;
        n_ref        << -0.40824829, 0.81649658, -0.40824829;
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT((kmax_dir_out - kmax_dir_ref).norm()<1e-7);
        TEST_ASSERT((kmin_dir_out - kmin_dir_ref).norm()<1e-7);
        TEST_ASSERT((n_out - n_ref).norm()<1e-7);
        TEST_ASSERT(kmax_dir_out.dot(kmin_dir_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmax_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmin_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());
        u=1; v=1;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, kmax_dir_out, kmin_dir_out, n_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        kmax_dir_ref <<  0.70710678,          0,  0.70710678;
        kmin_dir_ref <<  0.57735027,-0.57735027, -0.57735027;
        n_ref        <<  0.40824829, 0.81649658, -0.40824829;
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT((kmax_dir_out - kmax_dir_ref).norm()<1e-7);
        TEST_ASSERT((kmin_dir_out - kmin_dir_ref).norm()<1e-7);
        TEST_ASSERT((n_out - n_ref).norm()<1e-7);
        TEST_ASSERT(kmax_dir_out.dot(kmin_dir_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmax_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmin_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());

        // test at interior point u=v=1/2
        u=0.5; v=0.5;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, kmax_dir_out, kmin_dir_out, n_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        kmax_dir_ref <<  1, 0,  0;
        kmin_dir_ref <<  0, 0, -1;
        n_ref        <<  0, 1,  0;
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT((kmax_dir_out - kmax_dir_ref).norm()<1e-7);
        TEST_ASSERT((kmin_dir_out - kmin_dir_ref).norm()<1e-7);
        TEST_ASSERT((n_out - n_ref).norm()<1e-7);
        TEST_ASSERT(kmax_dir_out.dot(kmin_dir_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmax_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmin_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());

        // test at interior point u=1/4 & v=3/4
        u=0.25; v=0.75;
        eli::geom::surface::principal_curvature(kmax_out, kmin_out, kmax_dir_out, kmin_dir_out, n_out, bez, u, v);
        eli::geom::surface::mean_curvature(H_out, bez, u, v);
        eli::geom::surface::gaussian_curvature(K_out, bez, u, v);
        kmax_dir_ref <<  0.70710678,          0,  -0.70710678;
        kmin_dir_ref << -0.69879657, -0.15286175, -0.69879657;
        n_ref        << -0.10808958,  0.98824758, -0.10808958;
        TEST_ASSERT(std::abs((kmax_out+kmin_out)/2-H_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(std::abs((kmax_out*kmin_out)-K_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT((kmax_dir_out - kmax_dir_ref).norm()<1e-7);
        TEST_ASSERT((kmin_dir_out - kmin_dir_ref).norm()<1e-7);
        TEST_ASSERT((n_out - n_ref).norm()<1e-7);
        TEST_ASSERT(kmax_dir_out.dot(kmin_dir_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmax_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(kmin_dir_out.dot(n_out)<std::numeric_limits<data_type>::epsilon());
      }
    }

    void promotion_test()
    {
      index_type n(3), m(3);
      point_type pt[3+1][3+1];
      data_type u, v;

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

      // create surface with specified dimensions and set control points
      bezier_type bez(n, m);

      for (index_type i=0; i<=n; ++i)
      {
        for (index_type j=0; j<=m; ++j)
        {
          bez.set_control_point(pt[i][j], i, j);
        }
      }

      // test promote in u-direction
      {
        bezier_type bez2(bez);

        // promote
        bez2.promote_u();

        // test if degree increased
        TEST_ASSERT(bez2.degree_u()==bez.degree_u()+1);

        // test evaluation at u=v=1/2
        u=0.5; v=0.5;
        TEST_ASSERT(bez.f(u, v)==bez2.f(u, v));
        TEST_ASSERT(bez.f_u(u, v)==bez2.f_u(u, v));
        TEST_ASSERT(bez.f_v(u, v)==bez2.f_v(u, v));
        TEST_ASSERT(bez.f_uu(u, v)==bez2.f_uu(u, v));
        TEST_ASSERT(bez.f_uv(u, v)==bez2.f_uv(u, v));
        TEST_ASSERT(bez.f_vv(u, v)==bez2.f_vv(u, v));
        TEST_ASSERT(bez.f_uuu(u, v)==bez2.f_uuu(u, v));
        TEST_ASSERT(bez.f_uuv(u, v)==bez2.f_uuv(u, v));
        TEST_ASSERT(bez.f_uvv(u, v)==bez2.f_uvv(u, v));
        TEST_ASSERT(bez.f_vvv(u, v)==bez2.f_vvv(u, v));

        // test evaluation at interior point u=1/4, v=3/4
        u=0.25; v=0.75;
        TEST_ASSERT(bez.f(u, v)==bez2.f(u, v));
        TEST_ASSERT(bez.f_u(u, v)==bez2.f_u(u, v));
        TEST_ASSERT(bez.f_v(u, v)==bez2.f_v(u, v));
        TEST_ASSERT(bez.f_uu(u, v)==bez2.f_uu(u, v));
        TEST_ASSERT(bez.f_uv(u, v)==bez2.f_uv(u, v));
        TEST_ASSERT(bez.f_vv(u, v)==bez2.f_vv(u, v));
        TEST_ASSERT(bez.f_uuu(u, v)==bez2.f_uuu(u, v));
        TEST_ASSERT(bez.f_uuv(u, v)==bez2.f_uuv(u, v));
        TEST_ASSERT(bez.f_uvv(u, v)==bez2.f_uvv(u, v));
        TEST_ASSERT(bez.f_vvv(u, v)==bez2.f_vvv(u, v));
      }

      // test promote in v-direction
      {
        bezier_type bez2(bez);

        // promote
        bez2.promote_v();

        // test if degree increased
        TEST_ASSERT(bez2.degree_v()==bez.degree_v()+1);

        // test evaluation at u=v=1/2
        u=0.5; v=0.5;
        TEST_ASSERT(bez.f(u, v)==bez2.f(u, v));
        TEST_ASSERT(bez.f_u(u, v)==bez2.f_u(u, v));
        TEST_ASSERT(bez.f_v(u, v)==bez2.f_v(u, v));
        TEST_ASSERT(bez.f_uu(u, v)==bez2.f_uu(u, v));
        TEST_ASSERT(bez.f_uv(u, v)==bez2.f_uv(u, v));
        TEST_ASSERT(bez.f_vv(u, v)==bez2.f_vv(u, v));
        TEST_ASSERT(bez.f_uuu(u, v)==bez2.f_uuu(u, v));
        TEST_ASSERT(bez.f_uuv(u, v)==bez2.f_uuv(u, v));
        TEST_ASSERT(bez.f_uvv(u, v)==bez2.f_uvv(u, v));
        TEST_ASSERT(bez.f_vvv(u, v)==bez2.f_vvv(u, v));

        // test evaluation at interior point u=1/4, v=3/4
        u=0.25; v=0.75;
        TEST_ASSERT(bez.f(u, v)==bez2.f(u, v));
        TEST_ASSERT(bez.f_u(u, v)==bez2.f_u(u, v));
        TEST_ASSERT(bez.f_v(u, v)==bez2.f_v(u, v));
        TEST_ASSERT(bez.f_uu(u, v)==bez2.f_uu(u, v));
        TEST_ASSERT(bez.f_uv(u, v)==bez2.f_uv(u, v));
        TEST_ASSERT(bez.f_vv(u, v)==bez2.f_vv(u, v));
        TEST_ASSERT(bez.f_uuu(u, v)==bez2.f_uuu(u, v));
        TEST_ASSERT(bez.f_uuv(u, v)==bez2.f_uuv(u, v));
        TEST_ASSERT(bez.f_uvv(u, v)==bez2.f_uvv(u, v));
        TEST_ASSERT(bez.f_vvv(u, v)==bez2.f_vvv(u, v));
      }
    }

    void promotion_to_test()
    {
      index_type n(3), m(3);
      point_type pt[3+1][3+1];
      data_type u, v;

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

      // create surface with specified dimensions and set control points
      bezier_type bez(n, m);

      for (index_type i=0; i<=n; ++i)
      {
        for (index_type j=0; j<=m; ++j)
        {
          bez.set_control_point(pt[i][j], i, j);
        }
      }

      // test promote in u-direction
      {
        bezier_type bez2(bez);
        bezier_type bez3(bez);

        // promote
        bez2.promote_u();

        // promote to +1
        bez3.promote_u_to(bez.degree_u()+1);

        // test if degree increased
        TEST_ASSERT(bez3.degree_u()==bez.degree_u()+1);

        for (index_type i=0; i<=n+1; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            TEST_ASSERT(bez2.get_control_point(i, j)==bez3.get_control_point(i, j));
          }
        }

      }

      // test promote in v-direction
      {
        bezier_type bez2(bez);
        bezier_type bez3(bez);

        // promote
        bez2.promote_v();

        bez3.promote_v_to(bez.degree_v()+1);

        // test if degree increased
        TEST_ASSERT(bez3.degree_v()==bez.degree_v()+1);

        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m+1; ++j)
          {
            TEST_ASSERT(bez2.get_control_point(i, j)==bez3.get_control_point(i, j));
          }
        }
      }

      // test promote in u-direction
      {
        bezier_type bez2(bez);

        index_type inc = 3;

        // promote
        bez2.promote_u_to(bez.degree_u()+inc);

        // test if degree increased
        TEST_ASSERT(bez2.degree_u()==bez.degree_u()+inc);

        // test evaluation at u=v=1/2
        u=0.5; v=0.5;
        TEST_ASSERT(bez.f(u, v)==bez2.f(u, v));
        TEST_ASSERT(bez.f_u(u, v)==bez2.f_u(u, v));
        TEST_ASSERT(bez.f_v(u, v)==bez2.f_v(u, v));
        TEST_ASSERT(bez.f_uu(u, v)==bez2.f_uu(u, v));
        TEST_ASSERT(bez.f_uv(u, v)==bez2.f_uv(u, v));
        TEST_ASSERT(bez.f_vv(u, v)==bez2.f_vv(u, v));
        TEST_ASSERT(bez.f_uuu(u, v)==bez2.f_uuu(u, v));
        TEST_ASSERT(bez.f_uuv(u, v)==bez2.f_uuv(u, v));
        TEST_ASSERT(bez.f_uvv(u, v)==bez2.f_uvv(u, v));
        TEST_ASSERT(bez.f_vvv(u, v)==bez2.f_vvv(u, v));

        // test evaluation at interior point u=1/4, v=3/4
        u=0.25; v=0.75;
        TEST_ASSERT(bez.f(u, v)==bez2.f(u, v));
        TEST_ASSERT(bez.f_u(u, v)==bez2.f_u(u, v));
        TEST_ASSERT(bez.f_v(u, v)==bez2.f_v(u, v));
        TEST_ASSERT(bez.f_uu(u, v)==bez2.f_uu(u, v));
        TEST_ASSERT(bez.f_uv(u, v)==bez2.f_uv(u, v));
        TEST_ASSERT(bez.f_vv(u, v)==bez2.f_vv(u, v));
        TEST_ASSERT(bez.f_uuu(u, v)==bez2.f_uuu(u, v));
        TEST_ASSERT(bez.f_uuv(u, v)==bez2.f_uuv(u, v));
        TEST_ASSERT(bez.f_uvv(u, v)==bez2.f_uvv(u, v));
        TEST_ASSERT(bez.f_vvv(u, v)==bez2.f_vvv(u, v));
      }

      // test promote in v-direction
      {
        bezier_type bez2(bez);

        index_type inc = 3;

        // promote
        bez2.promote_v_to(bez.degree_v()+inc);

        // test if degree increased
        TEST_ASSERT(bez2.degree_v()==bez.degree_v()+inc);

        // test evaluation at u=v=1/2
        u=0.5; v=0.5;
        TEST_ASSERT(bez.f(u, v)==bez2.f(u, v));
        TEST_ASSERT(bez.f_u(u, v)==bez2.f_u(u, v));
        TEST_ASSERT(bez.f_v(u, v)==bez2.f_v(u, v));
        TEST_ASSERT(bez.f_uu(u, v)==bez2.f_uu(u, v));
        TEST_ASSERT(bez.f_uv(u, v)==bez2.f_uv(u, v));
        TEST_ASSERT(bez.f_vv(u, v)==bez2.f_vv(u, v));
        TEST_ASSERT(bez.f_uuu(u, v)==bez2.f_uuu(u, v));
        TEST_ASSERT(bez.f_uuv(u, v)==bez2.f_uuv(u, v));
        TEST_ASSERT(bez.f_uvv(u, v)==bez2.f_uvv(u, v));
        TEST_ASSERT(bez.f_vvv(u, v)==bez2.f_vvv(u, v));

        // test evaluation at interior point u=1/4, v=3/4
        u=0.25; v=0.75;
        TEST_ASSERT(bez.f(u, v)==bez2.f(u, v));
        TEST_ASSERT(bez.f_u(u, v)==bez2.f_u(u, v));
        TEST_ASSERT(bez.f_v(u, v)==bez2.f_v(u, v));
        TEST_ASSERT(bez.f_uu(u, v)==bez2.f_uu(u, v));
        TEST_ASSERT(bez.f_uv(u, v)==bez2.f_uv(u, v));
        TEST_ASSERT(bez.f_vv(u, v)==bez2.f_vv(u, v));
        TEST_ASSERT(bez.f_uuu(u, v)==bez2.f_uuu(u, v));
        TEST_ASSERT(bez.f_uuv(u, v)==bez2.f_uuv(u, v));
        TEST_ASSERT(bez.f_uvv(u, v)==bez2.f_uvv(u, v));
        TEST_ASSERT(bez.f_vvv(u, v)==bez2.f_vvv(u, v));
      }
    }

    void demotion_test()
    {
      // hand checked with "Degree Reduction of Bezier Curves" by Dave Morgan
      {
        index_type n(6), m(4);
        point_type pt[6+1][4+1], pt_ref[5+1][4+1];
        bool demoted;

        // create surface with specified control points
        pt[0][0] <<  0, 0, 15;
        pt[1][0] <<  2, 6, 15;
        pt[2][0] <<  3, 0, 15;
        pt[3][0] <<  5, 4, 15;
        pt[4][0] <<  7, 1, 15;
        pt[5][0] <<  5, 5, 15;
        pt[6][0] << 10, 6, 15;
        pt[0][1] <<  0, 0, 11;
        pt[1][1] <<  2, 6, 11;
        pt[2][1] <<  3, 0, 11;
        pt[3][1] <<  5, 4, 11;
        pt[4][1] <<  7, 1, 11;
        pt[5][1] <<  5, 5, 11;
        pt[6][1] << 10, 6, 11;
        pt[0][2] <<  0, 0,  3;
        pt[1][2] <<  2, 6,  3;
        pt[2][2] <<  3, 0,  3;
        pt[3][2] <<  5, 4,  3;
        pt[4][2] <<  7, 1,  3;
        pt[5][2] <<  5, 5,  3;
        pt[6][2] << 10, 6,  3;
        pt[0][3] <<  0, 0,  0;
        pt[1][3] <<  2, 6,  0;
        pt[2][3] <<  3, 0,  0;
        pt[3][3] <<  5, 4,  0;
        pt[4][3] <<  7, 1,  0;
        pt[5][3] <<  5, 5,  0;
        pt[6][3] << 10, 6,  0;
        pt[0][4] <<  0, 0, -5;
        pt[1][4] <<  2, 6, -5;
        pt[2][4] <<  3, 0, -5;
        pt[3][4] <<  5, 4, -5;
        pt[4][4] <<  7, 1, -5;
        pt[5][4] <<  5, 5, -5;
        pt[6][4] << 10, 6, -5;

        // create surface with specified dimensions and set control points
        bezier_type bez(n, m), bez2;

        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            bez.set_control_point(pt[i][j], i, j);
          }
        }
        bez2=bez;

        demoted=bez2.demote_u(eli::geom::general::NOT_CONNECTED);

        TEST_ASSERT(demoted);

        // test if the degree is correct
        TEST_ASSERT(bez.degree_u()==bez2.degree_u()+1);
        TEST_ASSERT(bez.degree_v()==bez2.degree_v());

        // set the reference control points
        pt_ref[0][0] <<  -0.00878906, 0.0610352,  15;
        pt_ref[1][0] <<   2.5177734,  6.3821289,  15;
        pt_ref[2][0] <<   2.8060547, -0.16982422, 15;
        pt_ref[3][0] <<   8.0060547,  2.5301758,  15;
        pt_ref[4][0] <<   4.1177734,  3.9821289,  15;
        pt_ref[5][0] <<   9.9912109,  6.0610352,  15;
        pt_ref[0][1] <<  -0.00878906, 0.0610352,  11;
        pt_ref[1][1] <<   2.5177734,  6.3821289,  11;
        pt_ref[2][1] <<   2.8060547, -0.16982422, 11;
        pt_ref[3][1] <<   8.0060547,  2.5301758,  11;
        pt_ref[4][1] <<   4.1177734,  3.9821289,  11;
        pt_ref[5][1] <<   9.9912109,  6.0610352,  11;
        pt_ref[0][2] <<  -0.00878906, 0.0610352,   3;
        pt_ref[1][2] <<   2.5177734,  6.3821289,   3;
        pt_ref[2][2] <<   2.8060547, -0.16982422,  3;
        pt_ref[3][2] <<   8.0060547,  2.5301758,   3;
        pt_ref[4][2] <<   4.1177734,  3.9821289,   3;
        pt_ref[5][2] <<   9.9912109,  6.0610352,   3;
        pt_ref[0][3] <<  -0.00878906, 0.0610352,   0;
        pt_ref[1][3] <<   2.5177734,  6.3821289,   0;
        pt_ref[2][3] <<   2.8060547, -0.16982422,  0;
        pt_ref[3][3] <<   8.0060547,  2.5301758,   0;
        pt_ref[4][3] <<   4.1177734,  3.9821289,   0;
        pt_ref[5][3] <<   9.9912109,  6.0610352,   0;
        pt_ref[0][4] <<  -0.00878906, 0.0610352,  -5;
        pt_ref[1][4] <<   2.5177734,  6.3821289,  -5;
        pt_ref[2][4] <<   2.8060547, -0.16982422, -5;
        pt_ref[3][4] <<   8.0060547,  2.5301758,  -5;
        pt_ref[4][4] <<   4.1177734,  3.9821289,  -5;
        pt_ref[5][4] <<   9.9912109,  6.0610352,  -5;

        // test if get the correct control points
        for (index_type j=0; j<=m; ++j)
        {
          for (index_type i=0; i<=n-1; ++i)
          {
            TEST_ASSERT((bez2.get_control_point(i,j)-pt_ref[i][j]).norm()<1e-6);
          }
        }
      }

      // test whether can promote and demote and get the same control points
      {
        index_type n(3), m(3);
        point_type pt[3+1][3+1];
        bool demoted;

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

        // create surface with specified dimensions and set control points
        bezier_type bez(n, m), bez2;

        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            bez.set_control_point(pt[i][j], i, j);
          }
        }
        bez2=bez;

        // promote -> demote u-direction
        bez2.promote_u();
        demoted=bez2.demote_u(eli::geom::general::NOT_CONNECTED);

        TEST_ASSERT(demoted);

        // test to see if degree has changed
        TEST_ASSERT(bez.degree_u()==bez2.degree_u());
        TEST_ASSERT(bez.degree_v()==bez2.degree_v());

        // test to see if get expected polygon
        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            TEST_ASSERT(bez2.get_control_point(i, j)==bez.get_control_point(i, j));
          }
        }

        // promote -> demote v-direction
        bez2=bez;
        bez2.promote_v();
        demoted=bez2.demote_v(eli::geom::general::NOT_CONNECTED);

        TEST_ASSERT(demoted);

        // test to see if degree has changed
        TEST_ASSERT(bez.degree_u()==bez2.degree_u());
        TEST_ASSERT(bez.degree_v()==bez2.degree_v());

        // test to see if get expected polygon
        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            TEST_ASSERT(bez2.get_control_point(i, j)==bez.get_control_point(i, j));
          }
        }
      }

      // test whether can demote and maintain continuity at end points
      {
        index_type n(6), m(4);
        point_type pt[6+1][4+1];

        // create surface with specified control points
        pt[0][0] <<  0, 0, 15;
        pt[1][0] <<  2, 6, 15;
        pt[2][0] <<  3, 0, 15;
        pt[3][0] <<  5, 4, 15;
        pt[4][0] <<  7, 1, 15;
        pt[5][0] <<  5, 5, 15;
        pt[6][0] << 10, 6, 15;
        pt[0][1] <<  0, 0, 11;
        pt[1][1] <<  2, 6, 11;
        pt[2][1] <<  3, 0, 11;
        pt[3][1] <<  5, 4, 11;
        pt[4][1] <<  7, 1, 11;
        pt[5][1] <<  5, 5, 11;
        pt[6][1] << 10, 6, 11;
        pt[0][2] <<  0, 0,  3;
        pt[1][2] <<  2, 6,  3;
        pt[2][2] <<  3, 0,  3;
        pt[3][2] <<  5, 4,  3;
        pt[4][2] <<  7, 1,  3;
        pt[5][2] <<  5, 5,  3;
        pt[6][2] << 10, 6,  3;
        pt[0][3] <<  0, 0,  0;
        pt[1][3] <<  2, 6,  0;
        pt[2][3] <<  3, 0,  0;
        pt[3][3] <<  5, 4,  0;
        pt[4][3] <<  7, 1,  0;
        pt[5][3] <<  5, 5,  0;
        pt[6][3] << 10, 6,  0;
        pt[0][4] <<  0, 0, -5;
        pt[1][4] <<  2, 6, -5;
        pt[2][4] <<  3, 0, -5;
        pt[3][4] <<  5, 4, -5;
        pt[4][4] <<  7, 1, -5;
        pt[5][4] <<  5, 5, -5;
        pt[6][4] << 10, 6, -5;

        // create surface with specified dimensions and set control points
        bezier_type bez(n, m);

        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            bez.set_control_point(pt[i][j], i, j);
          }
        }

        // No End Constraint
        {
          bezier_type bez2(bez);
          point_type eval_out[6], eval_ref[6];
          data_type d[6], u[6]={0, 1, 0, 1, 0.5, 0.25}, v[6]={0, 0, 1, 1, 0.5, 0.75};
          bool demoted;

          // demote u-direction
          demoted=bez2.demote_u(eli::geom::general::NOT_CONNECTED);
          TEST_ASSERT(demoted);
          TEST_ASSERT(bez.degree_u()==bez2.degree_u()+1);
          TEST_ASSERT(bez.degree_v()==bez2.degree_v());

          for (size_t k=0; k<6; ++k)
          {
            eval_out[0]=bez2.f(u[k], v[k]);
            eval_out[1]=bez2.f_u(u[k], v[k]);
            eval_out[2]=bez2.f_v(u[k], v[k]);
            eval_out[3]=bez2.f_uu(u[k], v[k]);
            eval_out[4]=bez2.f_uv(u[k], v[k]);
            eval_out[5]=bez2.f_vv(u[k], v[k]);
            eval_ref[0]=bez.f(u[k], v[k]);
            eval_ref[1]=bez.f_u(u[k], v[k]);
            eval_ref[2]=bez.f_v(u[k], v[k]);
            eval_ref[3]=bez.f_uu(u[k], v[k]);
            eval_ref[4]=bez.f_uv(u[k], v[k]);
            eval_ref[5]=bez.f_vv(u[k], v[k]);
            for (index_type i=0; i<6; ++i)
            {
              d[i]=eli::geom::point::distance(eval_out[i], eval_ref[i]);
            }

            // boundaries have different error than interior points
            if (k<4)
            {
              TEST_ASSERT_DELTA(d[0], 0.061665, 5e-6);
              TEST_ASSERT_DELTA(d[1], 4.43986, 5e-5);
              TEST_ASSERT_DELTA(d[2], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3], 103.597, 5e-3);
              TEST_ASSERT_DELTA(d[4], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 0, std::numeric_limits<data_type>::epsilon());
            }
            else if (k==4)
            {
              TEST_ASSERT_DELTA(d[0], 0.061665, 5e-6);
              TEST_ASSERT_DELTA(d[1], 0, 20*std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3], 8.8797199432, 5e-6);
              TEST_ASSERT_DELTA(d[4], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 0, std::numeric_limits<data_type>::epsilon());
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 0.061665, 5e-6);
              TEST_ASSERT_DELTA(d[1], 0, 20*std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3], 11.83962659, 5e-6);
              TEST_ASSERT_DELTA(d[4], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 0, std::numeric_limits<data_type>::epsilon());
            }
          }

          // demote v-direction
          bez2=bez;
          demoted=bez2.demote_v(eli::geom::general::NOT_CONNECTED);
          TEST_ASSERT(demoted);
          TEST_ASSERT(bez.degree_u()==bez2.degree_u());
          TEST_ASSERT(bez.degree_v()==bez2.degree_v()+1);

          for (size_t k=0; k<6; ++k)
          {
            eval_out[0]=bez2.f(u[k], v[k]);
            eval_out[1]=bez2.f_u(u[k], v[k]);
            eval_out[2]=bez2.f_v(u[k], v[k]);
            eval_out[3]=bez2.f_uu(u[k], v[k]);
            eval_out[4]=bez2.f_uv(u[k], v[k]);
            eval_out[5]=bez2.f_vv(u[k], v[k]);
            eval_ref[0]=bez.f(u[k], v[k]);
            eval_ref[1]=bez.f_u(u[k], v[k]);
            eval_ref[2]=bez.f_v(u[k], v[k]);
            eval_ref[3]=bez.f_uu(u[k], v[k]);
            eval_ref[4]=bez.f_uv(u[k], v[k]);
            eval_ref[5]=bez.f_vv(u[k], v[k]);
            for (index_type i=0; i<6; ++i)
            {
              d[i]=eli::geom::point::distance(eval_out[i], eval_ref[i]);
            }

            // boundaries have different error than interior points
            if (k<4)
            {
              TEST_ASSERT_DELTA(d[0],  0.125, 1e-6);
              TEST_ASSERT_DELTA(d[1],  0,     std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2],  4,     1e-5);
              TEST_ASSERT_DELTA(d[3],  0,     std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[4],  0,     std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 40,     1e-5);
            }
            else if (k==4)
            {
              TEST_ASSERT_DELTA(d[0],  0.125, 1e-6);
              TEST_ASSERT_DELTA(d[1],  0,     std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2],  0,     std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3],  0,     std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[4],  0,     std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5],  8,     1e-6);
            }
            else
            {
              TEST_ASSERT_DELTA(d[0],  0.0625, 1e-6);
              TEST_ASSERT_DELTA(d[1],  0,      std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2],  1,      1e-6);
              TEST_ASSERT_DELTA(d[3],  0,      std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[4],  0,      std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5],  4,      1e-6);
            }
          }
        }

        // Point End Constraint
        {
          bezier_type bez2(bez);
          point_type eval_out[6], eval_ref[6];
          data_type d[6], u[6]={0, 1, 0, 1, 0.5, 0.25}, v[6]={0, 0, 1, 1, 0.5, 0.75};
          bool demoted;

          // demote u-direction
          demoted=bez2.demote_u(eli::geom::general::C0);
          TEST_ASSERT(demoted);
          TEST_ASSERT(bez.degree_u()==bez2.degree_u()+1);
          TEST_ASSERT(bez.degree_v()==bez2.degree_v());

          for (size_t k=0; k<6; ++k)
          {
            eval_out[0]=bez2.f(u[k], v[k]);
            eval_out[1]=bez2.f_u(u[k], v[k]);
            eval_out[2]=bez2.f_v(u[k], v[k]);
            eval_out[3]=bez2.f_uu(u[k], v[k]);
            eval_out[4]=bez2.f_uv(u[k], v[k]);
            eval_out[5]=bez2.f_vv(u[k], v[k]);
            eval_ref[0]=bez.f(u[k], v[k]);
            eval_ref[1]=bez.f_u(u[k], v[k]);
            eval_ref[2]=bez.f_v(u[k], v[k]);
            eval_ref[3]=bez.f_uu(u[k], v[k]);
            eval_ref[4]=bez.f_uv(u[k], v[k]);
            eval_ref[5]=bez.f_vv(u[k], v[k]);
            for (index_type i=0; i<6; ++i)
            {
              d[i]=eli::geom::point::distance(eval_out[i], eval_ref[i]);
            }

            // boundaries have different error than interior points
            if (k<4)
            {
              TEST_ASSERT_DELTA(d[0], 0, 61*std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[1], 3.82695, 1e-5);
              TEST_ASSERT_DELTA(d[2], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3], 99.5007, 1e-4);
              TEST_ASSERT_DELTA(d[4], 0, 41*std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 0, std::numeric_limits<data_type>::epsilon());
            }
            else if (k==4)
            {
              TEST_ASSERT_DELTA(d[0], 0.0597961, 1e-6);
              TEST_ASSERT_DELTA(d[1], 0, 20*std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3], 9.08901, 1e-5);
              TEST_ASSERT_DELTA(d[4], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 0, std::numeric_limits<data_type>::epsilon());
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 0.0644677, 1e-6);
              TEST_ASSERT_DELTA(d[1], 0.0373726, 1e-6);
              TEST_ASSERT_DELTA(d[2], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3], 12.7067,  1e-4);
              TEST_ASSERT_DELTA(d[4], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 0, 9*std::numeric_limits<data_type>::epsilon());
            }
          }

          // demote v-direction
          bez2=bez;
          demoted=bez2.demote_v(eli::geom::general::C0);
          TEST_ASSERT(demoted);
          TEST_ASSERT(bez.degree_u()==bez2.degree_u());
          TEST_ASSERT(bez.degree_v()==bez2.degree_v()+1);

          for (size_t k=0; k<6; ++k)
          {
            eval_out[0]=bez2.f(u[k], v[k]);
            eval_out[1]=bez2.f_u(u[k], v[k]);
            eval_out[2]=bez2.f_v(u[k], v[k]);
            eval_out[3]=bez2.f_uu(u[k], v[k]);
            eval_out[4]=bez2.f_uv(u[k], v[k]);
            eval_out[5]=bez2.f_vv(u[k], v[k]);
            eval_ref[0]=bez.f(u[k], v[k]);
            eval_ref[1]=bez.f_u(u[k], v[k]);
            eval_ref[2]=bez.f_v(u[k], v[k]);
            eval_ref[3]=bez.f_uu(u[k], v[k]);
            eval_ref[4]=bez.f_uv(u[k], v[k]);
            eval_ref[5]=bez.f_vv(u[k], v[k]);
            for (index_type i=0; i<6; ++i)
            {
              d[i]=eli::geom::point::distance(eval_out[i], eval_ref[i]);
            }

            // boundaries have different error than interior points
            if (k<4)
            {
              TEST_ASSERT_DELTA(d[0],  0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[1],  0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2],  3.428571, 1e-5);
              TEST_ASSERT_DELTA(d[3],  0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[4],  0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 38.8571, 1e-4);
            }
            else if (k==4)
            {
              TEST_ASSERT_DELTA(d[0],  0.1428571, 1e-6);
              TEST_ASSERT_DELTA(d[1],  0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2],  0, 17*std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3],  0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[4],  0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 9.142857, 1e-5);
            }
            else
            {
              TEST_ASSERT_DELTA(d[0],  0.0803571, 1e-6);
              TEST_ASSERT_DELTA(d[1],  0,      std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2],  1.285714,      1e-6);
              TEST_ASSERT_DELTA(d[3],  0,      std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[4],  0,      std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5],  2.857143,      1e-5);
            }
          }
        }

        // Slope End Constraint
        {
          bezier_type bez2(bez);
          point_type eval_out[6], eval_ref[6];
          data_type d[6], u[6]={0, 1, 0, 1, 0.5, 0.25}, v[6]={0, 0, 1, 1, 0.5, 0.75};
          bool demoted;

          // demote u-direction
          demoted=bez2.demote_u(eli::geom::general::C1);
          TEST_ASSERT(demoted);
          TEST_ASSERT(bez.degree_u()==bez2.degree_u()+1);
          TEST_ASSERT(bez.degree_v()==bez2.degree_v());

          for (size_t k=0; k<6; ++k)
          {
            eval_out[0]=bez2.f(u[k], v[k]);
            eval_out[1]=bez2.f_u(u[k], v[k]);
            eval_out[2]=bez2.f_v(u[k], v[k]);
            eval_out[3]=bez2.f_uu(u[k], v[k]);
            eval_out[4]=bez2.f_uv(u[k], v[k]);
            eval_out[5]=bez2.f_vv(u[k], v[k]);
            eval_ref[0]=bez.f(u[k], v[k]);
            eval_ref[1]=bez.f_u(u[k], v[k]);
            eval_ref[2]=bez.f_v(u[k], v[k]);
            eval_ref[3]=bez.f_uu(u[k], v[k]);
            eval_ref[4]=bez.f_uv(u[k], v[k]);
            eval_ref[5]=bez.f_vv(u[k], v[k]);
            for (index_type i=0; i<6; ++i)
            {
              d[i]=eli::geom::point::distance(eval_out[i], eval_ref[i]);
            }

            // boundaries have different error than interior points
            if (k<4)
            {
              TEST_ASSERT_DELTA(d[0], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[1], 0, 20*std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3], 57.40425, 1e-4);
              TEST_ASSERT_DELTA(d[4], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 0, std::numeric_limits<data_type>::epsilon());
            }
            else if (k==4)
            {
              TEST_ASSERT_DELTA(d[0], 0.179388, 1e-6);
              TEST_ASSERT_DELTA(d[1], 0, 20*std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[2], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3], 18.65638, 1e-5);
              TEST_ASSERT_DELTA(d[4], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 0, std::numeric_limits<data_type>::epsilon());
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 0.1765853, 1e-6);
              TEST_ASSERT_DELTA(d[1], 1.278142,  1e-5);
              TEST_ASSERT_DELTA(d[2], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[3], 16.05525,  1e-5);
              TEST_ASSERT_DELTA(d[4], 0, std::numeric_limits<data_type>::epsilon());
              TEST_ASSERT_DELTA(d[5], 0, std::numeric_limits<data_type>::epsilon());
            }
          }

          // demote v-direction (cannot do, degree is to low)
          bez2=bez;
          demoted=bez2.demote_v(eli::geom::general::C1);
          TEST_ASSERT(!demoted);
        }

        // Second Derivative Constraint
        {
          bezier_type bez2(bez);
          bool demoted;

          // demote u-direction (cannot do, degree is to low)
          demoted=bez2.demote_u(eli::geom::general::C2);
          TEST_ASSERT(!demoted);

          // demote v-direction (cannot do, degree is to low)
          bez2=bez;
          demoted=bez2.demote_v(eli::geom::general::C2);
          TEST_ASSERT(!demoted);
        }
      }
    }

    void split_test()
    {
      // split with known control points
      {
        index_type n(3), m(3);
        point_type pt[3+1][3+1], pt_lo[3+1][3+1], pt_hi[3+1][3+1];

        // create surface with specified control points
        pt[0][0] << 0, 0, 0;
        pt[1][0] << 0, 2, 0;
        pt[2][0] << 8, 2, 0;
        pt[3][0] << 4, 0, 0;
        pt[0][1] << 0, 0, 1;
        pt[1][1] << 0, 2, 1;
        pt[2][1] << 8, 2, 1;
        pt[3][1] << 4, 0, 1;
        pt[0][2] << 0, 0, 3;
        pt[1][2] << 0, 2, 3;
        pt[2][2] << 8, 2, 3;
        pt[3][2] << 4, 0, 3;
        pt[0][3] << 0, 0, 4;
        pt[1][3] << 0, 2, 4;
        pt[2][3] << 8, 2, 4;
        pt[3][3] << 4, 0, 4;

        // create surface with specified dimensions and set control points
        bezier_type bez(n, m), bez_lo, bez_hi;

        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            bez.set_control_point(pt[i][j], i, j);
          }
        }

        // split the surface
        bez.split_u(bez_lo, bez_hi, 0.5);

        // set the reference control points
        pt_lo[0][0] << 0.0, 0.0, 0;
        pt_lo[1][0] << 0.0, 1.0, 0;
        pt_lo[2][0] << 2.0, 1.5, 0;
        pt_lo[3][0] << 3.5, 1.5, 0;
        pt_lo[0][1] << 0.0, 0.0, 1;
        pt_lo[1][1] << 0.0, 1.0, 1;
        pt_lo[2][1] << 2.0, 1.5, 1;
        pt_lo[3][1] << 3.5, 1.5, 1;
        pt_lo[0][2] << 0.0, 0.0, 3;
        pt_lo[1][2] << 0.0, 1.0, 3;
        pt_lo[2][2] << 2.0, 1.5, 3;
        pt_lo[3][2] << 3.5, 1.5, 3;
        pt_lo[0][3] << 0.0, 0.0, 4;
        pt_lo[1][3] << 0.0, 1.0, 4;
        pt_lo[2][3] << 2.0, 1.5, 4;
        pt_lo[3][3] << 3.5, 1.5, 4;

        pt_hi[0][0] << 3.5, 1.5, 0;
        pt_hi[1][0] << 5.0, 1.5, 0;
        pt_hi[2][0] << 6.0, 1.0, 0;
        pt_hi[3][0] << 4.0, 0.0, 0;
        pt_hi[0][1] << 3.5, 1.5, 1;
        pt_hi[1][1] << 5.0, 1.5, 1;
        pt_hi[2][1] << 6.0, 1.0, 1;
        pt_hi[3][1] << 4.0, 0.0, 1;
        pt_hi[0][2] << 3.5, 1.5, 3;
        pt_hi[1][2] << 5.0, 1.5, 3;
        pt_hi[2][2] << 6.0, 1.0, 3;
        pt_hi[3][2] << 4.0, 0.0, 3;
        pt_hi[0][3] << 3.5, 1.5, 4;
        pt_hi[1][3] << 5.0, 1.5, 4;
        pt_hi[2][3] << 6.0, 1.0, 4;
        pt_hi[3][3] << 4.0, 0.0, 4;

        // compare control points
        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            TEST_ASSERT(bez_lo.get_control_point(i,j)==pt_lo[i][j]);
            TEST_ASSERT(bez_hi.get_control_point(i,j)==pt_hi[i][j]);
          }
        }
      }

      // split in u-direction
      {
        index_type n(3), m(3);
        point_type pt[3+1][3+1];
        data_type u0(0.25), ul, uh, vl, vh;

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

        // create surface with specified dimensions and set control points
        bezier_type bez(n, m), bez_lo, bez_hi;
        point_type pt_out, pt_ref;

        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            bez.set_control_point(pt[i][j], i, j);
          }
        }

        // split surface
        bez.split_u(bez_lo, bez_hi, u0);

        // test lower side matches
        ul=0.125;
        vl=0.25;
        pt_ref=bez.f(ul, vl);
        pt_out=bez_lo.f(ul/u0, vl);
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_u(ul, vl);
        pt_out=bez_lo.f_u(ul/u0, vl)/u0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_v(ul, vl);
        pt_out=bez_lo.f_v(ul/u0, vl);
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_uu(ul, vl);
        pt_out=bez_lo.f_uu(ul/u0, vl)/u0/u0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_uv(ul, vl);
        pt_out=bez_lo.f_uv(ul/u0, vl)/u0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_vv(ul, vl);
        pt_out=bez_lo.f_vv(ul/u0, vl);
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_uuu(ul, vl);
        pt_out=bez_lo.f_uuu(ul/u0, vl)/u0/u0/u0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_uuv(ul, vl);
        pt_out=bez_lo.f_uuv(ul/u0, vl)/u0/u0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_uvv(ul, vl);
        pt_out=bez_lo.f_uvv(ul/u0, vl)/u0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_vvv(ul, vl);
        pt_out=bez_lo.f_vvv(ul/u0, vl);
        TEST_ASSERT(pt_ref==pt_out);

        // test upper side matches
        uh=0.75;
        vh=0.75;
        pt_ref=bez.f(uh, vh);
        pt_out=bez_hi.f((uh-u0)/(1-u0), vh);
        TEST_ASSERT((pt_ref - pt_out).norm()<5*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_u(uh, vh);
        pt_out=bez_hi.f_u((uh-u0)/(1-u0), vh)/(1-u0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_v(uh, vh);
        pt_out=bez_hi.f_v((uh-u0)/(1-u0), vh);
        TEST_ASSERT((pt_ref - pt_out).norm()<3*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_uu(uh, vh);
        pt_out=bez_hi.f_uu((uh-u0)/(1-u0), vh)/(1-u0)/(1-u0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_uv(uh, vh);
        pt_out=bez_hi.f_uv((uh-u0)/(1-u0), vh)/(1-u0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_vv(uh, vh);
        pt_out=bez_hi.f_vv((uh-u0)/(1-u0), vh);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_uuu(uh, vh);
        pt_out=bez_hi.f_uuu((uh-u0)/(1-u0), vh)/(1-u0)/(1-u0)/(1-u0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_uuv(uh, vh);
        pt_out=bez_hi.f_uuv((uh-u0)/(1-u0), vh)/(1-u0)/(1-u0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_uvv(uh, vh);
        pt_out=bez_hi.f_uvv((uh-u0)/(1-u0), vh)/(1-u0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_vvv(uh, vh);
        pt_out=bez_hi.f_vvv((uh-u0)/(1-u0), vh);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
      }

      // split in v-direction
      {
        index_type n(3), m(3);
        point_type pt[3+1][3+1];
        data_type v0(0.25), ul, uh, vl, vh;

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

        // create surface with specified dimensions and set control points
        bezier_type bez(n, m), bez_lo, bez_hi;
        point_type pt_out, pt_ref;

        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            bez.set_control_point(pt[i][j], i, j);
          }
        }

        // split surface
        bez.split_v(bez_lo, bez_hi, v0);

        // test lower side matches
        ul=0.25;
        vl=0.125;
        pt_ref=bez.f(ul, vl);
        pt_out=bez_lo.f(ul, vl/v0);
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_u(ul, vl);
        pt_out=bez_lo.f_u(ul, vl/v0);
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_v(ul, vl);
        pt_out=bez_lo.f_v(ul, vl/v0)/v0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_uu(ul, vl);
        pt_out=bez_lo.f_uu(ul, vl/v0);
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_uv(ul, vl);
        pt_out=bez_lo.f_uv(ul, vl/v0)/v0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_vv(ul, vl);
        pt_out=bez_lo.f_vv(ul, vl/v0)/v0/v0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_uuu(ul, vl);
        pt_out=bez_lo.f_uuu(ul, vl/v0);
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_uuv(ul, vl);
        pt_out=bez_lo.f_uuv(ul, vl/v0)/v0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_uvv(ul, vl);
        pt_out=bez_lo.f_uvv(ul, vl/v0)/v0/v0;
        TEST_ASSERT(pt_ref==pt_out);
        pt_ref=bez.f_vvv(ul, vl);
        pt_out=bez_lo.f_vvv(ul, vl/v0)/v0/v0/v0;
        TEST_ASSERT(pt_ref==pt_out);

        // test upper side matches
        uh=0.75;
        vh=0.75;
        pt_ref=bez.f(uh, vh);
        pt_out=bez_hi.f(uh, (vh-v0)/(1-v0));
        TEST_ASSERT((pt_ref - pt_out).norm()<5*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_u(uh, vh);
        pt_out=bez_hi.f_u(uh, (vh-v0)/(1-v0));
        TEST_ASSERT((pt_ref - pt_out).norm()<3*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_v(uh, vh);
        pt_out=bez_hi.f_v(uh, (vh-v0)/(1-v0))/(1-v0);
        TEST_ASSERT((pt_ref - pt_out).norm()<3*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_uu(uh, vh);
        pt_out=bez_hi.f_uu(uh, (vh-v0)/(1-v0));
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_uv(uh, vh);
        pt_out=bez_hi.f_uv(uh, (vh-v0)/(1-v0))/(1-v0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_vv(uh, vh);
        pt_out=bez_hi.f_vv(uh, (vh-v0)/(1-v0))/(1-v0)/(1-v0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_uuu(uh, vh);
        pt_out=bez_hi.f_uuu(uh, (vh-v0)/(1-v0));
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_uuv(uh, vh);
        pt_out=bez_hi.f_uuv(uh, (vh-v0)/(1-v0))/(1-v0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_uvv(uh, vh);
        pt_out=bez_hi.f_uvv(uh, (vh-v0)/(1-v0))/(1-v0)/(1-v0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
        pt_ref=bez.f_vvv(uh, vh);
        pt_out=bez_hi.f_vvv(uh, (vh-v0)/(1-v0))/(1-v0)/(1-v0)/(1-v0);
        TEST_ASSERT((pt_ref - pt_out).norm()<2*std::numeric_limits<data_type>::epsilon());
      }
    }

    void normal_test()
    {
      // simple surface test
      {
        index_type n(3), m(3);
        point_type pt[3+1][3+1], pt_out, pt_ref;
        data_type u, v;

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

        // create surface with specified dimensions and set control points
        bezier_type bez(n, m);

        for (index_type i=0; i<=n; ++i)
        {
          for (index_type j=0; j<=m; ++j)
          {
            bez.set_control_point(pt[i][j], i, j);
          }
        }

        // get the normal and the tangent vectors
        point_type normal, S_u, S_v, tmp;

        u=0.5;
        v=0.5;
        normal=bez.normal(u, v);
        S_u=bez.f_u(u, v);
        S_v=bez.f_v(u, v);

        // test the normal vector
        TEST_ASSERT(tol.approximately_equal(normal.norm(), 1));
        tmp=S_u.cross(S_v);
        tmp.normalize();
        TEST_ASSERT(tol.approximately_equal(tmp, normal));
      }

      // flat degenerate surface test
      {
        bezier_type bez(3, 3);
        point_type cp[4], x, y, origin;
        data_type k, xr, yr;
        index_type i;

        k=4*(eli::constants::math<data_type>::sqrt_two()-1)/3;
        x << 1, 0, 0;
        y << 0, 1, 0;
        origin << 0, 0, 0;

        // set the first section
        xr=0;
        yr=0;
        cp[0]=xr*x+origin;
        cp[1]=xr*x+yr*k*y+origin;
        cp[2]=xr*k*x+yr*y+origin;
        cp[3]=yr*y+origin;
        for (i=0; i<4; ++i)
        {
          bez.set_control_point(cp[i], i, 0);
        }

        // set the second section
        xr=0.5;
        yr=0.5;
        cp[0]=xr*x+origin;
        cp[1]=xr*x+yr*k*y+origin;
        cp[2]=xr*k*x+yr*y+origin;
        cp[3]=yr*y+origin;
        for (i=0; i<4; ++i)
        {
          bez.set_control_point(cp[i], i, 1);
        }

        // set the third section
        xr=1;
        yr=1;
        cp[0]=xr*x+origin;
        cp[1]=xr*x+yr*k*y+origin;
        cp[2]=xr*k*x+yr*y+origin;
        cp[3]=yr*y+origin;
        for (i=0; i<4; ++i)
        {
          bez.set_control_point(cp[i], i, 2);
        }

        // set the last section
        xr=1.5;
        yr=1.5;
        cp[0]=xr*x+origin;
        cp[1]=xr*x+yr*k*y+origin;
        cp[2]=xr*k*x+yr*y+origin;
        cp[3]=yr*y+origin;
        for (i=0; i<4; ++i)
        {
          bez.set_control_point(cp[i], i, 3);
        }

        // get the normal and the tangent vectors
        point_type normal, ref_normal, S_u, S_v;
        data_type u, v;

        u=0.5;
        v=0.5;
        normal=bez.normal(u, v);
        S_u=bez.f_u(u, v);
        S_v=bez.f_v(u, v);

        // test the normal vector
        TEST_ASSERT(tol.approximately_equal(normal.norm(), 1));
        ref_normal=S_u.cross(S_v);
        ref_normal.normalize();
        TEST_ASSERT(tol.approximately_equal(ref_normal, normal));

        u=0.5;
        v=0;
        normal=bez.normal(u, v);
        ref_normal=bez.normal(u, v+std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(tol.approximately_equal(ref_normal, normal));
      }

      // 3D degenerate surface test
      {
        bezier_type bez(3, 3);
        point_type cp[4], x, y, origin;
        data_type k, xr, yr;
        index_type i;

        k=4*(eli::constants::math<data_type>::sqrt_two()-1)/3;
        x << 1, 0, 0;
        y << 0, 1, 0;

        // set the first section
        xr=0;
        yr=0;
        origin << 0, 0, 0;
        cp[0]=xr*x+origin;
        cp[1]=xr*x+yr*k*y+origin;
        cp[2]=xr*k*x+yr*y+origin;
        cp[3]=yr*y+origin;
        for (i=0; i<4; ++i)
        {
          bez.set_control_point(cp[i], i, 0);
        }

        // set the second section
        xr=0.5;
        yr=0.5;
        origin << 0, 0, 1;
        cp[0]=xr*x+origin;
        cp[1]=xr*x+yr*k*y+origin;
        cp[2]=xr*k*x+yr*y+origin;
        cp[3]=yr*y+origin;
        for (i=0; i<4; ++i)
        {
          bez.set_control_point(cp[i], i, 1);
        }

        // set the third section
        xr=1;
        yr=1;
        origin << 0, 0, 2;
        cp[0]=xr*x+origin;
        cp[1]=xr*x+yr*k*y+origin;
        cp[2]=xr*k*x+yr*y+origin;
        cp[3]=yr*y+origin;
        for (i=0; i<4; ++i)
        {
          bez.set_control_point(cp[i], i, 2);
        }

        // set the last section
        xr=1.5;
        yr=1.5;
        origin << 0, 0, 3;
        cp[0]=xr*x+origin;
        cp[1]=xr*x+yr*k*y+origin;
        cp[2]=xr*k*x+yr*y+origin;
        cp[3]=yr*y+origin;
        for (i=0; i<4; ++i)
        {
          bez.set_control_point(cp[i], i, 3);
        }

        // get the normal and the tangent vectors
        point_type normal, ref_normal, S_u, S_v;
        data_type u, v;

        u=0.5;
        v=0.5;
        normal=bez.normal(u, v);
        S_u=bez.f_u(u, v);
        S_v=bez.f_v(u, v);

        // test the normal vector
        TEST_ASSERT(tol.approximately_equal(normal.norm(), 1));
        ref_normal=S_u.cross(S_v);
        ref_normal.normalize();
        TEST_ASSERT(tol.approximately_equal(ref_normal, normal));

        u=0.5;
        v=0;
        normal=bez.normal(u, v);
        ref_normal=bez.normal(u, v+std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(tol.approximately_equal(ref_normal, normal));
      }
    }
};

#endif

