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
    typedef eli::geom::curve::piecewise_cubic_spline_creator<data__, 3, tolerance_type> cubic_spline_creator_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_control_points_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_points_slopes_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_closed_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_closed_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_closed_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_closed_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_clamped_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_natural_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_periodic_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<float>::create_closed_cubic_spline_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_control_points_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_points_slopes_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_closed_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_closed_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_closed_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_closed_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_clamped_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_natural_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_periodic_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<double>::create_closed_cubic_spline_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_control_points_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_points_slopes_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_closed_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_closed_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_closed_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_closed_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_clamped_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_natural_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_periodic_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<long double>::create_closed_cubic_spline_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_control_points_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_points_slopes_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_closed_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_closed_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_closed_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_closed_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_clamped_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_natural_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_periodic_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<dd_real>::create_closed_cubic_spline_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_control_points_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_points_slopes_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_closed_piecewise_chip_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_closed_piecewise_cardinal_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_closed_piecewise_catmull_rom_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_closed_piecewise_kochanek_bartels_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_clamped_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_natural_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_periodic_cubic_spline_test);
      TEST_ADD(piecewise_spline_creator_test_suite<qd_real>::create_closed_cubic_spline_test);
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

    void create_control_points_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pt(16);
        point_type cp[4];
        std::vector<data_type> t(5);
        data_type dt;

        // set the points and times
        data_type k=4*(eli::constants::math<data_type>::sqrt_two()-1)/3;
        pt[0]  <<    2,  0,  0;
        pt[1]  <<    2,  k,  0;
        pt[2]  <<  2*k,  1,  0;
        pt[3]  <<    0,  1,  0;
        pt[4]  <<    0,  1,  0;
        pt[5]  << -2*k,  1,  0;
        pt[6]  <<   -2,  k,  0;
        pt[7]  <<   -2,  0,  0;
        pt[8]  <<   -2,  0,  0;
        pt[9]  <<   -2, -k,  0;
        pt[10] << -2*k, -1,  0;
        pt[11] <<    0, -1,  0;
        pt[12] <<    0, -1,  0;
        pt[13] <<  2*k, -1,  0;
        pt[14] <<    2, -k,  0;
        pt[15] <<    2,  0,  0;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;

        // set the control points
        spline_creator.set_segment_control_points(pt[0],  pt[1],  pt[2],  pt[3],  0);
        spline_creator.set_segment_control_points(pt[4],  pt[5],  pt[6],  pt[7],  1);
        spline_creator.set_segment_control_points(pt[8],  pt[9],  pt[10], pt[11], 2);
        spline_creator.set_segment_control_points(pt[12], pt[13], pt[14], pt[15], 3);

        // set the times
        spline_creator.set_t0(t[0]);
        spline_creator.set_segment_dt(t[1]-t[0], 0);
        spline_creator.set_segment_dt(t[2]-t[1], 1);
        spline_creator.set_segment_dt(t[3]-t[2], 2);
        spline_creator.set_segment_dt(t[4]-t[3], 3);

        // test control point settings
        spline_creator.get_segment_control_points(cp[0], cp[1], cp[2], cp[3], 0);
        TEST_ASSERT(cp[0]==pt[0]);
        TEST_ASSERT(cp[1]==pt[1]);
        TEST_ASSERT(cp[2]==pt[2]);
        TEST_ASSERT(cp[3]==pt[3]);
        spline_creator.get_segment_control_points(cp[0], cp[1], cp[2], cp[3], 1);
        TEST_ASSERT(cp[0]==pt[4]);
        TEST_ASSERT(cp[1]==pt[5]);
        TEST_ASSERT(cp[2]==pt[6]);
        TEST_ASSERT(cp[3]==pt[7]);
        spline_creator.get_segment_control_points(cp[0], cp[1], cp[2], cp[3], 2);
        TEST_ASSERT(cp[0]==pt[8]);
        TEST_ASSERT(cp[1]==pt[9]);
        TEST_ASSERT(cp[2]==pt[10]);
        TEST_ASSERT(cp[3]==pt[11]);
        spline_creator.get_segment_control_points(cp[0], cp[1], cp[2], cp[3], 3);
        TEST_ASSERT(cp[0]==pt[12]);
        TEST_ASSERT(cp[1]==pt[13]);
        TEST_ASSERT(cp[2]==pt[14]);
        TEST_ASSERT(cp[3]==pt[15]);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the spline
        TEST_ASSERT(spline_creator.create(pc));
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pt(16);
        point_type cp[4];
        data_type dt;

        // set the points and times
        data_type k=4*(eli::constants::math<data_type>::sqrt_two()-1)/3;
        pt[0]  <<    2,  0,  0;
        pt[1]  <<    2,  k,  0;
        pt[2]  <<  2*k,  1,  0;
        pt[3]  <<    0,  1,  0;
        pt[4]  <<    0,  1,  0;
        pt[5]  << -2*k,  1,  0;
        pt[6]  <<   -2,  k,  0;
        pt[7]  <<   -2,  0,  0;
        pt[8]  <<   -2,  0,  0;
        pt[9]  <<   -2, -k,  0;
        pt[10] << -2*k, -1,  0;
        pt[11] <<    0, -1,  0;
        pt[12] <<    0, -1,  0;
        pt[13] <<  2*k, -1,  0;
        pt[14] <<    2, -k,  0;
        pt[15] <<    2,  0,  0;

        // set the control points
        spline_creator.set_segment_control_points(pt[0],  pt[1],  pt[2],  pt[3],  0);
        spline_creator.set_segment_control_points(pt[4],  pt[5],  pt[6],  pt[7],  1);
        spline_creator.set_segment_control_points(pt[8],  pt[9],  pt[10], pt[11], 2);
        spline_creator.set_segment_control_points(pt[12], pt[13], pt[14], pt[15], 3);

        // test control point settings
        spline_creator.get_segment_control_points(cp[0], cp[1], cp[2], cp[3], 0);
        TEST_ASSERT(cp[0]==pt[0]);
        TEST_ASSERT(cp[1]==pt[1]);
        TEST_ASSERT(cp[2]==pt[2]);
        TEST_ASSERT(cp[3]==pt[3]);
        spline_creator.get_segment_control_points(cp[0], cp[1], cp[2], cp[3], 1);
        TEST_ASSERT(cp[0]==pt[4]);
        TEST_ASSERT(cp[1]==pt[5]);
        TEST_ASSERT(cp[2]==pt[6]);
        TEST_ASSERT(cp[3]==pt[7]);
        spline_creator.get_segment_control_points(cp[0], cp[1], cp[2], cp[3], 2);
        TEST_ASSERT(cp[0]==pt[8]);
        TEST_ASSERT(cp[1]==pt[9]);
        TEST_ASSERT(cp[2]==pt[10]);
        TEST_ASSERT(cp[3]==pt[11]);
        spline_creator.get_segment_control_points(cp[0], cp[1], cp[2], cp[3], 3);
        TEST_ASSERT(cp[0]==pt[12]);
        TEST_ASSERT(cp[1]==pt[13]);
        TEST_ASSERT(cp[2]==pt[14]);
        TEST_ASSERT(cp[3]==pt[15]);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the spline
        TEST_ASSERT(spline_creator.create(pc));
      }
    }

    void create_points_slopes_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pt(5), m0(4), m1(4);
        std::vector<data_type> t(5);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pt[0] << 2, 0, 0;
        pt[1] << 1, 1, 0;
        pt[2] << 2, 1, 1;
        pt[3] << 0, 2, 1;
        pt[4] << 0, 1, 0;
        m0[0] << 0, 1, 0;
        m1[0] << 1, 1, 1;
        m0[1] << 1, 1, 1;
        m1[1] << 0, 1, 0;
        m0[2] << 1, 0, 1;
        m1[2] << 1, 1, 1;
        m0[3] << 0, 1, 1;
        m1[3] << 1, 1, 0;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;

        // set the times
        spline_creator.set_t0(t[0]);
        spline_creator.set_segment_dt(t[1]-t[0], 0);
        spline_creator.set_segment_dt(t[2]-t[1], 1);
        spline_creator.set_segment_dt(t[3]-t[2], 2);
        spline_creator.set_segment_dt(t[4]-t[3], 3);

        // set the points and slopes
        spline_creator.set_segment_point_slope(pt[0], m0[0], pt[1], m1[0], 0);
        spline_creator.set_segment_point_slope(pt[1], m0[1], pt[2], m1[1], 1);
        spline_creator.set_segment_point_slope(pt[2], m0[2], pt[3], m1[2], 2);
        spline_creator.set_segment_point_slope(pt[3], m0[3], pt[4], m1[3], 3);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the spline
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pt.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // test the curve
        data_type small(std::numeric_limits<data_type>::epsilon());
        i=0;
        pt_a=pc.f(t[i]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt[i]));
        pt_b=pc.fp(t[i]);
        TEST_ASSERT(tol.approximately_equal(pt_b, m0[i]));
        for (i=1; i<nseg-1; ++i)
        {
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt[i]));
          TEST_ASSERT(tol.approximately_equal(pt_b, pt[i]));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, m1[i-1]));
          TEST_ASSERT(tol.approximately_equal(pt_b, m0[i]));
        }
        pt_a=pc.f(t[i+1]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt[i+1]));
        pt_b=pc.fp(t[i+1]);
        TEST_ASSERT(tol.approximately_equal(pt_b, m1[i]));
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pt(5), m0(4), m1(4);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pt[0] << 2, 0, 0;
        pt[1] << 1, 1, 0;
        pt[2] << 2, 1, 1;
        pt[3] << 0, 2, 1;
        pt[4] << 0, 1, 0;
        m0[0] << 0, 1, 0;
        m1[0] << 1, 1, 1;
        m0[1] << 1, 1, 1;
        m1[1] << 0, 1, 0;
        m0[2] << 1, 0, 1;
        m1[2] << 1, 1, 1;
        m0[3] << 0, 1, 1;
        m1[3] << 1, 1, 0;

        // set the points and slopes
        spline_creator.set_segment_point_slope(pt[0], m0[0], pt[1], m1[0], 0);
        spline_creator.set_segment_point_slope(pt[1], m0[1], pt[2], m1[1], 1);
        spline_creator.set_segment_point_slope(pt[2], m0[2], pt[3], m1[2], 2);
        spline_creator.set_segment_point_slope(pt[3], m0[3], pt[4], m1[3], 3);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the spline
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pt.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // test the curve
        data_type small(std::numeric_limits<data_type>::epsilon());
        i=0;
        pt_a=pc.f(i);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt[i]));
        pt_b=pc.fp(i);
        TEST_ASSERT(tol.approximately_equal(pt_b, m0[i]));
        for (i=1; i<nseg-1; ++i)
        {
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt[i]));
          TEST_ASSERT(tol.approximately_equal(pt_b, pt[i]));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, m1[i-1]));
          TEST_ASSERT(tol.approximately_equal(pt_b, m0[i]));
        }
        pt_a=pc.f(i+1);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt[i+1]));
        pt_b=pc.fp(i+1);
        TEST_ASSERT(tol.approximately_equal(pt_b, m1[i]));
      }
    }

    void create_piecewise_chip_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type dt;
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

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        spline_creator.set_chip(pts.begin(), eli::geom::general::NOT_CONNECTED);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        spline_creator.set_chip(pts.begin(), eli::geom::general::NOT_CONNECTED);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }
    }

    void create_closed_piecewise_chip_test()
    {
      // create non-smooth with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
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
        t[5]=10;

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }

        spline_creator.set_chip(pts.begin(), eli::geom::general::C0);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==t[5]-t[4]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(t[5]);
        pt_a=pc.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(t[5]);
        pt_a=pc.fp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(t[5]);
        pt_a=pc.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
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
        t[5]=10;

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }

        spline_creator.set_chip(pts.begin(), eli::geom::general::C1);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(t[5]);
        pt_a=pc.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(t[5]);
        pt_a=pc.fp(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(t[5]);
        pt_a=pc.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create non-smooth with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // set up the creator
        spline_creator.set_chip(pts.begin(), eli::geom::general::C0);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(nseg);
        pt_a=pc.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(nseg);
        pt_a=pc.fp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(nseg);
        pt_a=pc.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // set up the creator
        spline_creator.set_chip(pts.begin(), eli::geom::general::C1);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(nseg);
        pt_a=pc.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(nseg);
        pt_a=pc.fp(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(nseg);
        pt_a=pc.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }
    }

    void create_piecewise_cardinal_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type dt;
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

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        spline_creator.set_cardinal(pts.begin(), static_cast<data_type>(0.75), eli::geom::general::NOT_CONNECTED);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        spline_creator.set_cardinal(pts.begin(), static_cast<data_type>(0.75), eli::geom::general::NOT_CONNECTED);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }
    }

    void create_closed_piecewise_cardinal_spline_test()
    {
      // create non-smooth with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
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
        t[5]=10;

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }

        spline_creator.set_cardinal(pts.begin(), static_cast<data_type>(0.75), eli::geom::general::C0);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==t[5]-t[4]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(t[5]);
        pt_a=pc.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(t[5]);
        pt_a=pc.fp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(t[5]);
        pt_a=pc.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
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
        t[5]=10;

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }

        spline_creator.set_cardinal(pts.begin(), static_cast<data_type>(0.75), eli::geom::general::C1);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(t[5]);
        pt_a=pc.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(t[5]);
        pt_a=pc.fp(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(t[5]);
        pt_a=pc.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create non-smooth with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // set up the creator
        spline_creator.set_cardinal(pts.begin(), static_cast<data_type>(0.75), eli::geom::general::C0);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(nseg);
        pt_a=pc.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(nseg);
        pt_a=pc.fp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(nseg);
        pt_a=pc.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // set up the creator
        spline_creator.set_cardinal(pts.begin(), static_cast<data_type>(0.75), eli::geom::general::C1);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(nseg);
        pt_a=pc.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(nseg);
        pt_a=pc.fp(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(nseg);
        pt_a=pc.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }
    }

    void create_piecewise_catmull_rom_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type dt;
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

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        spline_creator.set_catmull_rom(pts.begin(), eli::geom::general::NOT_CONNECTED);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        spline_creator.set_catmull_rom(pts.begin(), eli::geom::general::NOT_CONNECTED);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }
    }

    void create_closed_piecewise_catmull_rom_spline_test()
    {
      // create non-smooth with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
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
        t[5]=10;

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }

        spline_creator.set_catmull_rom(pts.begin(), eli::geom::general::C0);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==t[5]-t[4]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(t[5]);
        pt_a=pc.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(t[5]);
        pt_a=pc.fp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(t[5]);
        pt_a=pc.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
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
        t[5]=10;

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }

        spline_creator.set_catmull_rom(pts.begin(), eli::geom::general::C1);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(t[5]);
        pt_a=pc.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(t[5]);
        pt_a=pc.fp(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(t[5]);
        pt_a=pc.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create non-smooth with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // set up the creator
        spline_creator.set_catmull_rom(pts.begin(), eli::geom::general::C0);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(nseg);
        pt_a=pc.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(nseg);
        pt_a=pc.fp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(nseg);
        pt_a=pc.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        // set up the creator
        spline_creator.set_catmull_rom(pts.begin(), eli::geom::general::C1);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(nseg);
        pt_a=pc.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(nseg);
        pt_a=pc.fp(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(nseg);
        pt_a=pc.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }
    }

    void create_piecewise_kochanek_bartels_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type dt, ten(0.75), bia(0.75), con(0);
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

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        spline_creator.set_kochanek_bartels(pts.begin(), ten, bia, con, eli::geom::general::NOT_CONNECTED);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt, ten(0.75), bia(0.75), con(0);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        spline_creator.set_kochanek_bartels(pts.begin(), ten, bia, con, eli::geom::general::NOT_CONNECTED);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }
      }
    }

    void create_closed_piecewise_kochanek_bartels_spline_test()
    {
      // create non-smooth with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt, ten(0.75), bia(0.75), con(0);
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
        t[5]=10;

        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        spline_creator.set_kochanek_bartels(pts.begin(), ten, bia, con, eli::geom::general::C0);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==t[5]-t[4]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==(static_cast<index_type>(pts.size())));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(t[5]);
        pt_a=pc.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(t[5]);
        pt_a=pc.fp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(t[5]);
        pt_a=pc.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt, ten(0.75), bia(0.75), con(0);
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
        t[5]=10;

        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        spline_creator.set_kochanek_bartels(pts.begin(), ten, bia, con, eli::geom::general::C1);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==t[5]-t[4]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==(static_cast<index_type>(pts.size())));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(t[5]);
        pt_a=pc.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(t[5]);
        pt_a=pc.fp(t[0]);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(t[5]);
        pt_a=pc.fpp(t[0]);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create non-smooth with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt, ten(0.75), bia(0.75), con(0);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        spline_creator.set_kochanek_bartels(pts.begin(), ten, bia, con, eli::geom::general::C0);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==(static_cast<index_type>(pts.size())));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(5);
        pt_a=pc.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(5);
        pt_a=pc.fp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(5);
        pt_a=pc.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }

      // create smooth with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt, ten(0.75), bia(0.75), con(0);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5, -1;
        pts[4] <<  1, 4,  1;

        spline_creator.set_kochanek_bartels(pts.begin(), ten, bia, con, eli::geom::general::C1);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==(static_cast<index_type>(pts.size())));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check if closed
        TEST_ASSERT(pc.closed());

        // check the continuity at each point
        for (i=1; i<(nseg-1); ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
        }

        // test the end conditions
        pt_b=pc.f(5);
        pt_a=pc.f(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fp(5);
        pt_a=pc.fp(0);
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        pt_b=pc.fpp(5);
        pt_a=pc.fpp(0);
        TEST_ASSERT(!tol.approximately_equal(pt_a, pt_b));
      }
    }

    void create_cubic_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type dt;
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

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        spline_creator.set_cubic_spline(pts.begin());

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        }

        // check the end points
        data_type small(std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(tol.approximately_equal(pts[0], pc.f(t[0])));
        pt_b=pc.fppp(t[1]*(1+small));
        pt_a=pc.fppp(t[1]*(1-small));
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        TEST_ASSERT(tol.approximately_equal(pts[4], pc.f(t[4])));
        pt_b=pc.fppp(t[3]*(1+small));
        pt_a=pc.fppp(t[3]*(1-small));
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        // set up the creator
        spline_creator.set_cubic_spline(pts.begin());

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        }

        // check the end points
        data_type small(std::numeric_limits<data_type>::epsilon());
        TEST_ASSERT(tol.approximately_equal(pts[0], pc.f(0)));
        pt_b=pc.fppp(1*(1+small));
        pt_a=pc.fppp(1*(1-small));
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        TEST_ASSERT(tol.approximately_equal(pts[4], pc.f(4)));
        pt_b=pc.fppp(3*(1+small));
        pt_a=pc.fppp(3*(1-small));
        TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
      }
    }

    void create_clamped_cubic_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type dt;
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

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        point_type m0, m1;
        m0 << -1, -1, 1;
        m1 << 1, 1, -1;
        spline_creator.set_clamped_cubic_spline(pts.begin(), m0, m1);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        }

        // check the end points
        TEST_ASSERT(tol.approximately_equal(pts[0], pc.f(t[0])));
        TEST_ASSERT(tol.approximately_equal(m0, pc.fp(t[0])));
        TEST_ASSERT(tol.approximately_equal(pts[4], pc.f(t[4])));
        TEST_ASSERT(tol.approximately_equal(m1, pc.fp(t[4])));
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        // set up the creator
        point_type m0, m1;
        m0 << -1, -1, 1;
        m1 << 1, 1, -1;
        spline_creator.set_clamped_cubic_spline(pts.begin(), m0, m1);

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        }

        // check the end points
        TEST_ASSERT(tol.approximately_equal(pts[0], pc.f(0)));
        TEST_ASSERT(tol.approximately_equal(m0, pc.fp(0)));
        TEST_ASSERT(tol.approximately_equal(pts[4], pc.f(4)));
        TEST_ASSERT(tol.approximately_equal(m1, pc.fp(4)));
      }
    }

    void create_natural_cubic_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type dt;
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

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        spline_creator.set_natural_cubic_spline(pts.begin());

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        }

        // check the end points
        pt_a.setZero();
        TEST_ASSERT(tol.approximately_equal(pts[0], pc.f(t[0])));
        TEST_ASSERT(tol.approximately_equal(pt_a, pc.fpp(t[0])));
        TEST_ASSERT(tol.approximately_equal(pts[4], pc.f(t[4])));
        TEST_ASSERT(tol.approximately_equal(pt_a, pc.fpp(t[4])));
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        // set up the creator
        spline_creator.set_natural_cubic_spline(pts.begin());

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        }

        // check the end points
        TEST_ASSERT(tol.approximately_equal(pts[0], pc.f(0)));
        pt_a.setZero();
        TEST_ASSERT(tol.approximately_equal(pt_a, pc.fpp(0)));
        TEST_ASSERT(tol.approximately_equal(pts[4], pc.f(4)));
        TEST_ASSERT(tol.approximately_equal(pt_a, pc.fpp(4)));
      }
    }

    void create_periodic_cubic_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        data_type dt;
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

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        spline_creator.set_periodic_cubic_spline(pts.begin());

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        }

        // check the end points
        TEST_ASSERT(tol.approximately_equal(pts[0], pc.f(t[0])));
        TEST_ASSERT(tol.approximately_equal(pc.fp(t[4]), pc.fp(t[0])));
        TEST_ASSERT(tol.approximately_equal(pts[4], pc.f(t[4])));
        TEST_ASSERT(tol.approximately_equal(pc.fpp(t[4]), pc.fpp(t[0])));
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(4);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        // set up the creator
        spline_creator.set_periodic_cubic_spline(pts.begin());

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        }

        // check the end points
        TEST_ASSERT(tol.approximately_equal(pts[0], pc.f(0)));
        TEST_ASSERT(tol.approximately_equal(pc.fp(0), pc.fp(4)));
        TEST_ASSERT(tol.approximately_equal(pts[4], pc.f(4)));
        TEST_ASSERT(tol.approximately_equal(pc.fpp(0), pc.fpp(4)));
      }
    }

    void create_closed_cubic_spline_test()
    {
      // create with specified times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(6);
        data_type dt;
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
        t[5]=10;

        // set up the creator
        spline_creator.set_t0(t[0]);
        for (i=0; i<spline_creator.get_number_segments(); ++i)
        {
          spline_creator.set_segment_dt(t[i+1]-t[i], i);
        }
        spline_creator.set_closed_cubic_spline(pts.begin());

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==t[0]);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==t[1]-t[0]);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==t[2]-t[1]);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==t[3]-t[2]);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==t[4]-t[3]);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==t[5]-t[4]);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(t[i]*(1+small));
          pt_a=pc.f(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(t[i]*(1+small));
          pt_a=pc.fp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(t[i]*(1+small));
          pt_a=pc.fpp(t[i]*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        }

        // check the end points
        TEST_ASSERT(tol.approximately_equal(pc.f(t[5]), pc.f(t[0])));
        TEST_ASSERT(tol.approximately_equal(pc.fp(t[5]), pc.fp(t[0])));
        TEST_ASSERT(tol.approximately_equal(pc.fpp(t[5]), pc.fpp(t[0])));
      }

      // create with default times
      {
        piecewise_curve_type pc;
        cubic_spline_creator_type spline_creator(5);
        point_type pt_b, pt_a;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        data_type dt;
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;

        // set up the creator
        spline_creator.set_closed_cubic_spline(pts.begin());

        // test time step settings
        TEST_ASSERT(spline_creator.get_t0()==0);
        dt=spline_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);
        dt=spline_creator.get_segment_dt(4);
        TEST_ASSERT(dt==1);

        // create the piecewise curve
        TEST_ASSERT(spline_creator.create(pc));

        // check the number segments
        nseg=pc.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(0, pc.get_t0()));

        // check the continuity at each point
        for (i=1; i<nseg; ++i)
        {
          data_type small(std::numeric_limits<data_type>::epsilon());
          pt_b=pc.f(i*(1+small));
          pt_a=pc.f(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fp(i*(1+small));
          pt_a=pc.fp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
          pt_b=pc.fpp(i*(1+small));
          pt_a=pc.fpp(i*(1-small));
          TEST_ASSERT(tol.approximately_equal(pt_a, pt_b));
        }

        // check the end points
        TEST_ASSERT(tol.approximately_equal(pc.f(5), pc.f(0)));
        TEST_ASSERT(tol.approximately_equal(pc.fp(5), pc.fp(0)));
        TEST_ASSERT(tol.approximately_equal(pc.fpp(5), pc.fpp(0)));
      }
    }
};

#endif

