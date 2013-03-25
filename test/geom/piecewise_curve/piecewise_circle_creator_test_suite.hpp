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

#ifndef piecewise_circle_creator_test_suite_hpp
#define piecewise_circle_creator_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_circle_creator.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_circle_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::piecewise_circle_creator<data__, 3, tolerance_type> circle_creator_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_circle_creator_test_suite<float>::create_circle_start_origin_test);
      TEST_ADD(piecewise_circle_creator_test_suite<float>::create_circle_start_origin_normal_test);
      TEST_ADD(piecewise_circle_creator_test_suite<float>::create_circle_3_point_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_circle_creator_test_suite<double>::create_circle_start_origin_test);
      TEST_ADD(piecewise_circle_creator_test_suite<double>::create_circle_start_origin_normal_test);
      TEST_ADD(piecewise_circle_creator_test_suite<double>::create_circle_3_point_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_circle_creator_test_suite<long double>::create_circle_start_origin_test);
      TEST_ADD(piecewise_circle_creator_test_suite<long double>::create_circle_start_origin_normal_test);
      TEST_ADD(piecewise_circle_creator_test_suite<long double>::create_circle_3_point_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_circle_creator_test_suite<dd_real>::create_circle_start_origin_test);
      TEST_ADD(piecewise_circle_creator_test_suite<dd_real>::create_circle_start_origin_normal_test);
      TEST_ADD(piecewise_circle_creator_test_suite<dd_real>::create_circle_3_point_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_circle_creator_test_suite<qd_real>::create_circle_start_origin_test);
      TEST_ADD(piecewise_circle_creator_test_suite<qd_real>::create_circle_start_origin_normal_test);
      TEST_ADD(piecewise_circle_creator_test_suite<qd_real>::create_circle_3_point_test);
    }
#endif
  public:
    piecewise_circle_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_circle_creator_test_suite()
    {
    }

  private:
    void octave_print(int figno, const piecewise_curve_type &pc) const
    {
      index_type i, pp, ns;
      data_type tmin, tmax;

      ns=pc.number_segments();
      pc.get_parameter_min(tmin);
      pc.get_parameter_max(tmax);

      std::cout << "figure(" << figno << ");" << std::endl;

      // get control points and print
      std::cout << "cp_x=[";
      for (pp=0; pp<ns; ++pp)
      {
        curve_type bez;
        pc.get(bez, pp);
        for (i=0; i<=bez.degree(); ++i)
        {
          std::cout << bez.get_control_point(i).x();
          if (i<bez.degree())
            std::cout << ", ";
          else if (pp<ns-1)
            std::cout << "; ";
        }
        std::cout << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "cp_y=[";
      for (pp=0; pp<ns; ++pp)
      {
        curve_type bez;
        pc.get(bez, pp);
        for (i=0; i<=bez.degree(); ++i)
        {
          std::cout << bez.get_control_point(i).y();
          if (i<bez.degree())
            std::cout << ", ";
          else if (pp<ns-1)
            std::cout << "; ";
        }
        std::cout << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "cp_z=[";
      for (pp=0; pp<ns; ++pp)
      {
        curve_type bez;
        pc.get(bez, pp);
        for (i=0; i<=bez.degree(); ++i)
        {
          std::cout << bez.get_control_point(i).z();
          if (i<bez.degree())
            std::cout << ", ";
          else if (pp<ns-1)
            std::cout << "; ";
        }
        std::cout << std::endl;
      }
      std::cout << "];" << std::endl;

      // initialize the t parameters
      std::vector<data__> t(129);
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(t.size()-1);
      }

      // set the surface points
      std::cout << "surf_x=[";
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        std::cout << pc.f(t[i]).x();
        if (i<static_cast<index_type>(t.size()-1))
          std::cout << ", ";
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_y=[";
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        std::cout << pc.f(t[i]).y();
        if (i<static_cast<index_type>(t.size()-1))
          std::cout << ", ";
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_z=[";
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        std::cout << pc.f(t[i]).z();
        if (i<static_cast<index_type>(t.size()-1))
          std::cout << ", ";
      }
      std::cout << "];" << std::endl;

      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot3(surf_x, surf_y, surf_z, '-k');" << std::endl;
      std::cout << "hold on;" << std::endl;
      std::cout << "plot3(cp_x', cp_y', cp_z', '-ok', 'MarkerFaceColor', [0 0 0]);" << std::endl;
      std::cout << "hold off;" << std::endl;
    }

    void create_circle_start_origin_test()
    {
      // create circle in x-z plane
      {
        circle_creator_type circle_creator;
        piecewise_curve_type pc;
        point_type start, origin, normal;

        // set the parameters for circle
        start  << 1, 0, 0;
        origin << 0, 0, 0;

        circle_creator.set(start, origin);

        // create the circle
        TEST_ASSERT(circle_creator.create(pc));
      }

      // create circle with zero radius
      {
        circle_creator_type circle_creator;
        piecewise_curve_type pc;
        point_type start, origin, normal;

        // set the parameters for circle
        start << 1, 0, 0;
        origin << 1, 0, 0;

        circle_creator.set(start, origin);

        // create the circle
        TEST_ASSERT(circle_creator.create(pc));

        TEST_ASSERT(pc.f(0.5)==pc.f(1.75));
      }
    }

    void create_circle_start_origin_normal_test()
    {
      // create circle in x-z plane
      {
        circle_creator_type circle_creator;
        piecewise_curve_type pc;
        point_type start, origin, normal;

        // set the parameters for circle
        start  << 1, 0, 0;
        origin << 0, 0, 0;
        normal << 0, 0, 1;

        circle_creator.set(start, origin, normal);

        // create the circle
        TEST_ASSERT(circle_creator.create(pc));
      }

      // create circle in 3d space
      {
        circle_creator_type circle_creator;
        piecewise_curve_type pc;
        point_type start, origin, normal;

        // set the parameters for circle
        start  << 2, 2, 1;
        origin << 1, 1, 1;
        normal << 1,-1, 2;

        circle_creator.set(start, origin, normal);

        // create the circle
        TEST_ASSERT(circle_creator.create(pc));
      }

      // create circle with zero radius
      {
        circle_creator_type circle_creator;
        piecewise_curve_type pc;
        point_type start, origin, normal;

        // set the parameters for circle
        start << 1, 0, 0;
        origin << 1, 0, 0;
        normal << 0, 0, 1;

        circle_creator.set(start, origin, normal);

        // create the circle
        TEST_ASSERT(circle_creator.create(pc));

        TEST_ASSERT(pc.f(0.5)==pc.f(1.75));
      }
    }

    void create_circle_3_point_test()
    {
      // test 2D
      {
        typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data_type, 2> piecewise_curve2_type;
        typedef typename piecewise_curve2_type::curve_type curve2_type;
        typedef typename piecewise_curve2_type::point_type point2_type;
        typedef eli::geom::curve::piecewise_circle_creator<data__, 2, tolerance_type> circle_creator2_type;

        point2_type origin(2, 4), start, middle, end, x, y, xref, yref;
        data_type radius(2), alpha1(1), alpha2(3), alpha3(5);
        piecewise_curve2_type pc;
        circle_creator2_type circle_creator;

        // set the three points
        start  << radius*std::cos(alpha1), radius*std::sin(alpha1);
        middle << radius*std::cos(alpha2), radius*std::sin(alpha2);
        end    << radius*std::cos(alpha3), radius*std::sin(alpha3);
        start+=origin;
        middle+=origin;
        end+=origin;

        circle_creator.set_3pt(start, middle, end);

        // test the circle creator parameters
        TEST_ASSERT(tol.approximately_equal(radius, circle_creator.get_radius()));
        TEST_ASSERT(tol.approximately_equal(origin, circle_creator.get_origin()));
        circle_creator.get_xy_directions(xref, yref);
        x=start-origin;
        x.normalize();
        y << -x.y(), x.x();
        TEST_ASSERT(tol.approximately_equal(x, xref));
        TEST_ASSERT(tol.approximately_equal(y, yref));

        // create the circle
        TEST_ASSERT(circle_creator.create(pc));
      }

      // test 3D
      {
        point_type origin(2, 4, 1), start, middle, end, x, y, xref, yref;
        data_type radius(2), alpha1(1), alpha2(3), alpha3(5);
        piecewise_curve_type pc;
        circle_creator_type circle_creator;
        Eigen::Matrix<data_type, 3, 3> rotx, roty;

        // set the three points
        start  << radius*std::cos(alpha1), radius*std::sin(alpha1), 0;
        middle << radius*std::cos(alpha2), radius*std::sin(alpha2), 0;
        end    << radius*std::cos(alpha3), radius*std::sin(alpha3), 0;
        rotx << 1,             0,              0,
                0, std::cos(0.5), -std::sin(0.5),
                0, std::sin(0.5),  std::cos(0.5);
        roty << std::cos(0.25), 0, std::sin(0.25),
                             0, 1,              0,
               -std::sin(0.25), 0, std::cos(0.25);
        start=start*rotx*roty+origin;
        middle=middle*rotx*roty+origin;
        end=end*rotx*roty+origin;

        circle_creator.set_3pt(start, middle, end);

        // test the circle creator parameters
        TEST_ASSERT(tol.approximately_equal(radius, circle_creator.get_radius()));
        TEST_ASSERT(tol.approximately_equal(origin, circle_creator.get_origin()));
        circle_creator.get_xy_directions(xref, yref);
        x=start-origin;
        x.normalize();
        y << -radius*std::sin(alpha1), radius*std::cos(alpha1), 0;
        y=y*rotx*roty;
        y.normalize();
        TEST_ASSERT(tol.approximately_equal(x, xref));
        TEST_ASSERT(tol.approximately_equal(y, yref));

        // create the circle
        TEST_ASSERT(circle_creator.create(pc));
      }
    }
};

#endif

