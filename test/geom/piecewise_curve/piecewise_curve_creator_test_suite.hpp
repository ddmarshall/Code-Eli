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

#ifndef piecewise_curve_creator_test_suite_hpp
#define piecewise_curve_creator_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_creator.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_curve_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_creator_test_suite<float>::create_piecewise_cubic_hermite_interpolating_polynomials_test);
      TEST_ADD(piecewise_curve_creator_test_suite<float>::create_circle_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_creator_test_suite<double>::create_piecewise_cubic_hermite_interpolating_polynomials_test);
      TEST_ADD(piecewise_curve_creator_test_suite<double>::create_circle_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_creator_test_suite<long double>::create_piecewise_cubic_hermite_interpolating_polynomials_test);
      TEST_ADD(piecewise_curve_creator_test_suite<long double>::create_circle_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_creator_test_suite<dd_real>::create_piecewise_cubic_hermite_interpolating_polynomials_test);
      TEST_ADD(piecewise_curve_creator_test_suite<dd_real>::create_circle_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_creator_test_suite<qd_real>::create_piecewise_cubic_hermite_interpolating_polynomials_test);
      TEST_ADD(piecewise_curve_creator_test_suite<qd_real>::create_circle_test);
    }
#endif
  public:
    piecewise_curve_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_curve_creator_test_suite()
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

    void create_piecewise_cubic_hermite_interpolating_polynomials_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc1;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_cubic_hermite_interpolating_polynomials(pc1, pts.begin(), pts.end(), t.begin()));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check the continuity at each point
        curve_type seg1, seg2;
        data_type t_out;
        pc1.get(seg1, t_out, 0);
        TEST_ASSERT(tol.approximately_equal(t_out, t[1]-t[0]));
        for (i=1; i<nseg; ++i)
        {
          pc1.get(seg2, t_out, i);
          TEST_ASSERT(tol.approximately_equal(t_out, t[i+1]-t[i]));
          TEST_ASSERT(tol.approximately_equal(seg1.fp(1), seg2.fp(0)));
          TEST_ASSERT(!tol.approximately_equal(seg1.fpp(1), seg2.fpp(0)));
          seg1=seg2;
        }
      }

      // create with default times
      {
      }
    }

    void create_circle_test()
    {
      // create circle in x-y plane
      {
        piecewise_curve_type pc1;
        point_type start, origin, normal;

        // set the parameters for circle
        start << 1, 0, 0;
        origin << 0, 0, 0;
        normal << 0, 0, 1;

        // create the circle
        TEST_ASSERT(create_circle_3(pc1, start, origin, normal));
      }

      // create circle in 3d space
      {
        piecewise_curve_type pc1;
        point_type start, origin, normal;

        // set the parameters for circle
        start << 2, 2, 1;
        origin << 1, 1, 1;
        normal << 1, -1, 2;

        // create the circle
        TEST_ASSERT(create_circle_3(pc1, start, origin, normal));
      }

      // create circle with zero radius
      {
        piecewise_curve_type pc1;
        point_type start, origin, normal;

        // set the parameters for circle
        start << 1, 0, 0;
        origin << 1, 0, 0;
        normal << 0, 0, 1;

        // create the circle
        TEST_ASSERT(create_circle_3(pc1, start, origin, normal));

        TEST_ASSERT(pc1.f(0.5)==pc1.f(0.75));
      }
    }

    void create_circular_arc_test()
    {
    }
};

#endif

