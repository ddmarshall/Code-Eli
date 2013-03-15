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

#ifndef piecewise_spline_creator_test_suite_hpp
#define piecewise_spline_creator_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_spline_creator.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_spline_creator_test_suite : public Test::Suite
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
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_closed_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_closed_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_closed_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_closed_piecewise_kochanek_bartels_spline_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_closed_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_closed_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_closed_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_closed_piecewise_kochanek_bartels_spline_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_closed_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_closed_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_closed_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_closed_piecewise_kochanek_bartels_spline_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_closed_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_closed_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_closed_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_closed_piecewise_kochanek_bartels_spline_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_closed_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_closed_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_closed_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_closed_piecewise_kochanek_bartels_spline_test);
    }
#endif
  public:
    piecewise_spline_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_spline_creator_test_suite()
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

    void create_piecewise_chip_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
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
        TEST_ASSERT(create_piecewise_chip(pc1, pts.begin(), pts.end(), t.begin(), eli::geom::general::NOT_CONNECTED));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }

      // create with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_chip(pc1, pts.begin(), pts.end(), eli::geom::general::NOT_CONNECTED));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }
    }

    void create_closed_piecewise_chip_test()
    {
      // create non-smooth with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type tend;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;
        tend=10;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_chip(pc1, pts.begin(), pts.end(), t.begin(), eli::geom::general::C0, tend));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(tend);
        pt_a=pc1.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(tend);
        pt_a=pc1.fp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(tend);
        pt_a=pc1.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type tend;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;
        tend=10;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_chip(pc1, pts.begin(), pts.end(), t.begin(), eli::geom::general::C1, tend));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(tend);
        pt_a=pc1.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(tend);
        pt_a=pc1.fp(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(tend);
        pt_a=pc1.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create non-smooth with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_chip(pc1, pts.begin(), pts.end(), eli::geom::general::C0));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(nseg);
        pt_a=pc1.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(nseg);
        pt_a=pc1.fp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(nseg);
        pt_a=pc1.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_chip(pc1, pts.begin(), pts.end(), eli::geom::general::C1));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(nseg);
        pt_a=pc1.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(nseg);
        pt_a=pc1.fp(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(nseg);
        pt_a=pc1.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }
    }

    void create_piecewise_cardinal_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
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
        TEST_ASSERT(create_piecewise_cardinal_spline(pc1, pts.begin(), pts.end(), t.begin(), static_cast<data_type>(0.75), eli::geom::general::NOT_CONNECTED));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }

      // create with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_cardinal_spline(pc1, pts.begin(), pts.end(), static_cast<data_type>(0.75), eli::geom::general::NOT_CONNECTED));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }
    }

    void create_closed_piecewise_cardinal_spline_test()
    {
      // create non-smooth with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type tend;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;
        tend=10;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_cardinal_spline(pc1, pts.begin(), pts.end(), t.begin(), static_cast<data_type>(0.75), eli::geom::general::C0, tend));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(tend);
        pt_a=pc1.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(tend);
        pt_a=pc1.fp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(tend);
        pt_a=pc1.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type tend;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;
        tend=10;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_cardinal_spline(pc1, pts.begin(), pts.end(), t.begin(), static_cast<data_type>(0.75), eli::geom::general::C1, tend));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));

//           if (typeid(data_type)==typeid(double))
//           {
//             std::cout << "joint " << i << std::endl;
//             std::cout << "    before=" << pc1.f(t[i]*(1-small)) << "\tafter=" << pc1.f(t[i]*(1+small)) << "\tdiff=" << pc1.f(t[i]*(1+small))-pc1.f(t[i]*(1-small)) << std::endl;
//             std::cout << "    before=" << pc1.fp(t[i]*(1-small)) << "\tafter=" << pc1.fp(t[i]*(1+small)) << "\tdiff=" << pc1.fp(t[i]*(1+small))-pc1.fp(t[i]*(1-small)) << std::endl;
//           }
        }

        // test the end conditions
        pt_b=pc1.f(tend);
        pt_a=pc1.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(tend);
        pt_a=pc1.fp(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(tend);
        pt_a=pc1.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));

//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << "end " << std::endl;
//           std::cout << "    before=" << pc1.f(t[0]) << "\tafter=" << pc1.f(tend) << "\tdiff=" << pc1.f(tend)-pc1.f(t[0]) << std::endl;
//           std::cout << "    before=" << pc1.fp(t[0]) << "\tafter=" << pc1.fp(tend) << "\tdiff=" << pc1.fp(tend)-pc1.fp(t[0]) << std::endl;
//         }
      }

      // create non-smooth with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_cardinal_spline(pc1, pts.begin(), pts.end(), static_cast<data_type>(0.75), eli::geom::general::C0));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(nseg);
        pt_a=pc1.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(nseg);
        pt_a=pc1.fp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(nseg);
        pt_a=pc1.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_cardinal_spline(pc1, pts.begin(), pts.end(), static_cast<data_type>(0.75), eli::geom::general::C1));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));

//           if (typeid(data_type)==typeid(double))
//           {
//             std::cout << "joint " << i << std::endl;
//             std::cout << "    before=" << pc1.f(t[i]*(1-small)) << "\tafter=" << pc1.f(t[i]*(1+small)) << "\tdiff=" << pc1.f(t[i]*(1+small))-pc1.f(t[i]*(1-small)) << std::endl;
//             std::cout << "    before=" << pc1.fp(t[i]*(1-small)) << "\tafter=" << pc1.fp(t[i]*(1+small)) << "\tdiff=" << pc1.fp(t[i]*(1+small))-pc1.fp(t[i]*(1-small)) << std::endl;
//           }
        }

        // test the end conditions
        pt_b=pc1.f(nseg);
        pt_a=pc1.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(nseg);
        pt_a=pc1.fp(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(nseg);
        pt_a=pc1.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));

//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << "end " << std::endl;
//           std::cout << "    before=" << pc1.f(t[0]) << "\tafter=" << pc1.f(tend) << "\tdiff=" << pc1.f(tend)-pc1.f(t[0]) << std::endl;
//           std::cout << "    before=" << pc1.fp(t[0]) << "\tafter=" << pc1.fp(tend) << "\tdiff=" << pc1.fp(tend)-pc1.fp(t[0]) << std::endl;
//         }
      }
    }

    void create_piecewise_catmull_rom_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
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
        TEST_ASSERT(create_piecewise_catmull_rom_spline(pc1, pts.begin(), pts.end(), t.begin(), eli::geom::general::NOT_CONNECTED));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }

      // create with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_catmull_rom_spline(pc1, pts.begin(), pts.end(), eli::geom::general::NOT_CONNECTED));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }
    }

    void create_closed_piecewise_catmull_rom_spline_test()
    {
      // create non-smooth with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type tend;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;
        tend=10;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_catmull_rom_spline(pc1, pts.begin(), pts.end(), t.begin(), eli::geom::general::C0, tend));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(tend);
        pt_a=pc1.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(tend);
        pt_a=pc1.fp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(tend);
        pt_a=pc1.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type tend;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;
        tend=10;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_catmull_rom_spline(pc1, pts.begin(), pts.end(), t.begin(), eli::geom::general::C1, tend));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(tend);
        pt_a=pc1.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(tend);
        pt_a=pc1.fp(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(tend);
        pt_a=pc1.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create non-smooth with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_catmull_rom_spline(pc1, pts.begin(), pts.end(), eli::geom::general::C0));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(nseg);
        pt_a=pc1.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(nseg);
        pt_a=pc1.fp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(nseg);
        pt_a=pc1.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_catmull_rom_spline(pc1, pts.begin(), pts.end(), eli::geom::general::C1));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(nseg);
        pt_a=pc1.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(nseg);
        pt_a=pc1.fp(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(nseg);
        pt_a=pc1.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }
    }

    void create_piecewise_kochanek_bartels_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
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
        TEST_ASSERT(create_piecewise_kochanek_bartels_spline(pc1, pts.begin(), pts.end(), t.begin(), static_cast<data_type>(0), static_cast<data_type>(0), static_cast<data_type>(0), eli::geom::general::NOT_CONNECTED));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }

      // create with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_kochanek_bartels_spline(pc1, pts.begin(), pts.end(), static_cast<data_type>(0.75), static_cast<data_type>(0.75), static_cast<data_type>(0), eli::geom::general::NOT_CONNECTED));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }
    }

    void create_closed_piecewise_kochanek_bartels_spline_test()
    {
      // create non-smooth with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type tend;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;
        tend=10;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_kochanek_bartels_spline(pc1, pts.begin(), pts.end(), t.begin(), static_cast<data_type>(0.75), static_cast<data_type>(0.75), static_cast<data_type>(0), eli::geom::general::C0, tend));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(tend);
        pt_a=pc1.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(tend);
        pt_a=pc1.fp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(tend);
        pt_a=pc1.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with specified times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type tend;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;
        tend=10;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_kochanek_bartels_spline(pc1, pts.begin(), pts.end(), t.begin(), static_cast<data_type>(0.75), static_cast<data_type>(0.75), static_cast<data_type>(0), eli::geom::general::C1, tend));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(t[i]*(1+small));
          pt_a=pc1.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(t[i]*(1+small));
          pt_a=pc1.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(t[i]*(1+small));
          pt_a=pc1.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));

//           if (typeid(data_type)==typeid(double))
//           {
//             std::cout << "joint " << i << std::endl;
//             std::cout << "    before=" << pc1.f(t[i]*(1-small)) << "\tafter=" << pc1.f(t[i]*(1+small)) << "\tdiff=" << pc1.f(t[i]*(1+small))-pc1.f(t[i]*(1-small)) << std::endl;
//             std::cout << "    before=" << pc1.fp(t[i]*(1-small)) << "\tafter=" << pc1.fp(t[i]*(1+small)) << "\tdiff=" << pc1.fp(t[i]*(1+small))-pc1.fp(t[i]*(1-small)) << std::endl;
//           }
        }

        // test the end conditions
        pt_b=pc1.f(tend);
        pt_a=pc1.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(tend);
        pt_a=pc1.fp(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(tend);
        pt_a=pc1.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));

//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << "end " << std::endl;
//           std::cout << "    before=" << pc1.f(t[0]) << "\tafter=" << pc1.f(tend) << "\tdiff=" << pc1.f(tend)-pc1.f(t[0]) << std::endl;
//           std::cout << "    before=" << pc1.fp(t[0]) << "\tafter=" << pc1.fp(tend) << "\tdiff=" << pc1.fp(tend)-pc1.fp(t[0]) << std::endl;
//         }
      }

      // create non-smooth with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_kochanek_bartels_spline(pc1, pts.begin(), pts.end(), static_cast<data_type>(0.75), static_cast<data_type>(0.75), static_cast<data_type>(0), eli::geom::general::C0));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc1.f(nseg);
        pt_a=pc1.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(nseg);
        pt_a=pc1.fp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(nseg);
        pt_a=pc1.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with default times
      {
        piecewise_curve_type pc1;
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_kochanek_bartels_spline(pc1, pts.begin(), pts.end(), static_cast<data_type>(0.75), static_cast<data_type>(0.75), static_cast<data_type>(0), eli::geom::general::C1));

        // check the number segments
        nseg=pc1.number_segments();
        TEST_ASSERT(nseg==(1+static_cast<index_type>(pts.size()-1)));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc1.get_t0()));

        // check if closed
        TEST_ASSERT(pc1.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc1.f(i*(1+small));
          pt_a=pc1.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fp(i*(1+small));
          pt_a=pc1.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc1.fpp(i*(1+small));
          pt_a=pc1.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));

//           if (typeid(data_type)==typeid(double))
//           {
//             std::cout << "joint " << i << std::endl;
//             std::cout << "    before=" << pc1.f(t[i]*(1-small)) << "\tafter=" << pc1.f(t[i]*(1+small)) << "\tdiff=" << pc1.f(t[i]*(1+small))-pc1.f(t[i]*(1-small)) << std::endl;
//             std::cout << "    before=" << pc1.fp(t[i]*(1-small)) << "\tafter=" << pc1.fp(t[i]*(1+small)) << "\tdiff=" << pc1.fp(t[i]*(1+small))-pc1.fp(t[i]*(1-small)) << std::endl;
//           }
        }

        // test the end conditions
        pt_b=pc1.f(nseg);
        pt_a=pc1.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fp(nseg);
        pt_a=pc1.fp(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc1.fpp(nseg);
        pt_a=pc1.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));

//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << "end " << std::endl;
//           std::cout << "    before=" << pc1.f(t[0]) << "\tafter=" << pc1.f(tend) << "\tdiff=" << pc1.f(tend)-pc1.f(t[0]) << std::endl;
//           std::cout << "    before=" << pc1.fp(t[0]) << "\tafter=" << pc1.fp(tend) << "\tdiff=" << pc1.fp(tend)-pc1.fp(t[0]) << std::endl;
//         }
      }
    }

#if 0
          if (typeid(data_type)==typeid(double))
          {
            std::cout << "joint " << i << std::endl;
            std::cout << "    before=" << pc1.f(t[i]*(1-small)) << "\tafter=" << pc1.f(t[i]*(1+small)) << "\tdiff=" << pc1.f(t[i]*(1+small))-pc1.f(t[i]*(1-small)) << std::endl;
            std::cout << "    before=" << pc1.fp(t[i]*(1-small)) << "\tafter=" << pc1.fp(t[i]*(1+small)) << "\tdiff=" << pc1.fp(t[i]*(1+small))-pc1.fp(t[i]*(1-small)) << std::endl;
          }
#endif
};

#endif

