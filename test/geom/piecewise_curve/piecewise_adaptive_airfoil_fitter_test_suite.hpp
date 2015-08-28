/*********************************************************************************
* Copyright (c) 2015 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef piecewise_adaptive_airfoil_fitter_test_suite_hpp
#define piecewise_adaptive_airfoil_fitter_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "octave_helpers.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_adaptive_airfoil_fitter.hpp"

template<typename data__>
class piecewise_adaptive_airfoil_fitter_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::control_point_type control_point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::piecewise_adaptive_airfoil_fitter<data_type, 3, tolerance_type> airfoil_fitter_type;
    typedef typename airfoil_fitter_type::subdivide_method subdivide_method;

    tolerance_type tol;

  private:
    void create_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &pts)
    {
      pts.resize(40);
      // lower surface
      pts[ 0] << 1.0000,-0.0020, 0;
      pts[ 1] << 0.9012,-0.0087, 0;
      pts[ 2] << 0.8078,-0.0247, 0;
      pts[ 3] << 0.7017,-0.0453, 0;
      pts[ 4] << 0.6178,-0.0572, 0;
      pts[ 5] << 0.4123,-0.0684, 0;
      pts[ 6] << 0.3545,-0.0678, 0;
      pts[ 7] << 0.2986,-0.0657, 0;
      pts[ 8] << 0.2453,-0.0621, 0;
      pts[ 9] << 0.1499,-0.0511, 0;
      pts[10] << 0.1029,-0.0441, 0;
      pts[11] << 0.0741,-0.0365, 0;
      pts[12] << 0.0451,-0.0284, 0;
      pts[13] << 0.0143,-0.0112, 0;
      pts[14] << 0.0076,-0.0116, 0;
      pts[15] << 0.0029,-0.0077, 0;
      pts[16] << 0.0000, 0.0000, 0;
      // upper surface
      pts[17] << 0.0013, 0.0030, 0;
      pts[18] << 0.0045, 0.0094, 0;
      pts[19] << 0.0096, 0.0138, 0;
      pts[20] << 0.0165, 0.0183, 0;
      pts[21] << 0.0252, 0.0228, 0;
      pts[22] << 0.0477, 0.0315, 0;
      pts[23] << 0.0767, 0.0397, 0;
      pts[24] << 0.1118, 0.0473, 0;
      pts[25] << 0.1524, 0.0539, 0;
      pts[26] << 0.1980, 0.0593, 0;
      pts[27] << 0.2979, 0.0636, 0;
      pts[28] << 0.3015, 0.0665, 0;
      pts[29] << 0.3578, 0.0680, 0;
      pts[30] << 0.4160, 0.0680, 0;
      pts[31] << 0.4455, 0.0675, 0;
      pts[32] << 0.5049, 0.0652, 0;
      pts[33] << 0.5930, 0.0585, 0;
      pts[34] << 0.6501, 0.0514, 0;
      pts[35] << 0.7050, 0.0416, 0;
      pts[36] << 0.7623, 0.0297, 0;
      pts[37] << 0.8168, 0.0221, 0;
      pts[38] << 0.9074, 0.0108, 0;
      pts[39] << 1.0000, 0.0050, 0;
    }

    void create_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &upt,
                               std::vector<point_type, Eigen::aligned_allocator<point_type>> &lpt, bool clean)
    {
      if (clean)
      {
        // upper surface
        upt.resize(23);
        upt[ 0] << 0.0000, 0.0000, 0;
        upt[ 1] << 0.0013, 0.0030, 0;
        upt[ 2] << 0.0045, 0.0094, 0;
        upt[ 3] << 0.0096, 0.0138, 0;
        upt[ 4] << 0.0165, 0.0183, 0;
        upt[ 5] << 0.0252, 0.0228, 0;
        upt[ 6] << 0.0477, 0.0315, 0;
        upt[ 7] << 0.0767, 0.0397, 0;
        upt[ 8] << 0.1118, 0.0473, 0;
        upt[ 9] << 0.1524, 0.0539, 0;
        upt[10] << 0.1980, 0.0593, 0;
        upt[11] << 0.3015, 0.0665, 0;
        upt[12] << 0.3578, 0.0680, 0;
        upt[13] << 0.4160, 0.0680, 0;
        upt[14] << 0.4455, 0.0675, 0;
        upt[15] << 0.5049, 0.0652, 0;
        upt[16] << 0.5930, 0.0585, 0;
        upt[17] << 0.6501, 0.0514, 0;
        upt[18] << 0.7050, 0.0416, 0;
        upt[19] << 0.7623, 0.0297, 0;
        upt[20] << 0.8168, 0.0221, 0;
        upt[21] << 0.9074, 0.0108, 0;
        upt[22] << 1.0000, 0.0050, 0;

        // lower surface
        lpt.resize(16);
        lpt[ 0] << 0.0000, 0.0000, 0;
        lpt[ 1] << 0.0029,-0.0077, 0;
        lpt[ 2] << 0.0076,-0.0116, 0;
        lpt[ 3] << 0.0451,-0.0284, 0;
        lpt[ 4] << 0.0741,-0.0365, 0;
        lpt[ 5] << 0.1029,-0.0441, 0;
        lpt[ 6] << 0.1499,-0.0511, 0;
        lpt[ 7] << 0.2453,-0.0621, 0;
        lpt[ 8] << 0.2986,-0.0657, 0;
        lpt[ 9] << 0.3545,-0.0678, 0;
        lpt[10] << 0.4123,-0.0684, 0;
        lpt[11] << 0.6178,-0.0572, 0;
        lpt[12] << 0.7017,-0.0453, 0;
        lpt[13] << 0.8078,-0.0247, 0;
        lpt[14] << 0.9012,-0.0087, 0;
        lpt[15] << 1.0000,-0.0020, 0;
      }
      else
      {
        // upper surface
        upt.resize(24);
        upt[ 0] << 0.0000, 0.0000, 0;
        upt[ 1] << 0.0013, 0.0030, 0;
        upt[ 2] << 0.0045, 0.0094, 0;
        upt[ 3] << 0.0096, 0.0138, 0;
        upt[ 4] << 0.0165, 0.0183, 0;
        upt[ 5] << 0.0252, 0.0228, 0;
        upt[ 6] << 0.0477, 0.0315, 0;
        upt[ 7] << 0.0767, 0.0397, 0;
        upt[ 8] << 0.1118, 0.0473, 0;
        upt[ 9] << 0.1524, 0.0539, 0;
        upt[10] << 0.1980, 0.0593, 0;
        upt[11] << 0.2979, 0.0636, 0;
        upt[12] << 0.3015, 0.0665, 0;
        upt[13] << 0.3578, 0.0680, 0;
        upt[14] << 0.4160, 0.0680, 0;
        upt[15] << 0.4455, 0.0675, 0;
        upt[16] << 0.5049, 0.0652, 0;
        upt[17] << 0.5930, 0.0585, 0;
        upt[18] << 0.6501, 0.0514, 0;
        upt[19] << 0.7050, 0.0416, 0;
        upt[20] << 0.7623, 0.0297, 0;
        upt[21] << 0.8168, 0.0221, 0;
        upt[22] << 0.9074, 0.0108, 0;
        upt[23] << 1.0000, 0.0050, 0;

        // lower surface
        lpt.resize(17);
        lpt[ 0] << 0.0000, 0.0000, 0;
        lpt[ 1] << 0.0029,-0.0077, 0;
        lpt[ 2] << 0.0076,-0.0116, 0;
        lpt[ 3] << 0.0143,-0.0112, 0;
        lpt[ 4] << 0.0451,-0.0284, 0;
        lpt[ 5] << 0.0741,-0.0365, 0;
        lpt[ 6] << 0.1029,-0.0441, 0;
        lpt[ 7] << 0.1499,-0.0511, 0;
        lpt[ 8] << 0.2453,-0.0621, 0;
        lpt[ 9] << 0.2986,-0.0657, 0;
        lpt[10] << 0.3545,-0.0678, 0;
        lpt[11] << 0.4123,-0.0684, 0;
        lpt[12] << 0.6178,-0.0572, 0;
        lpt[13] << 0.7017,-0.0453, 0;
        lpt[14] << 0.8078,-0.0247, 0;
        lpt[15] << 0.9012,-0.0087, 0;
        lpt[16] << 1.0000,-0.0020, 0;
      }
    }

    void create_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &upt,
                               std::vector<point_type, Eigen::aligned_allocator<point_type>> &lpt,
                               point_type & le_pt, point_type &te_pt)
    {
      // get the airfoil points
      create_airfoil_points(upt, lpt);

      // adjust the lower surface points to put trailing edge half-way
      // between upper and lower trailing edge points
      data_type y_offset(-0.003);

      for (index_type i=0; i<static_cast<index_type>(lpt.size()); ++i)
      {
        lpt[i].y()+=lpt[i].x()*y_offset;
      }

      // transform airfoil points
      data_type scale(2.5);
      data_type theta(30*eli::constants::math<data_type>::pi()/180);
      le_pt << 1, 2, 0;
      te_pt << le_pt.x()+std::cos(theta), le_pt.y()+std::sin(theta), 0;

      for (index_type i=0; i<static_cast<index_type>(upt.size()); ++i)
      {
        // (1) scale
        upt[i] = scale*upt[i];

        // (2) rotate
        data_type ptx(upt[i].x()), pty(upt[i].y()), ptz(upt[i].z()), ct(std::cos(theta)), st(std::sin(theta));
        upt[i] << (ptx*ct-pty*st), (ptx*st+pty*ct), ptz;

        // (3) translate
        upt[i] = upt[i]+le_pt;
      }
      for (index_type i=0; i<static_cast<index_type>(lpt.size()); ++i)
      {
        // (1) scale
        lpt[i] = scale*lpt[i];

        // (2) rotate
        data_type ptx(lpt[i].x()), pty(lpt[i].y()), ptz(lpt[i].z()), ct(std::cos(theta)), st(std::sin(theta));
        lpt[i] << (ptx*ct-pty*st), (ptx*st+pty*ct), ptz;

        // (3) translate
        lpt[i] = lpt[i]+le_pt;
      }
    }

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_adaptive_airfoil_fitter_test_suite<float>::fit_airfoil_simple_test);
      TEST_ADD(piecewise_adaptive_airfoil_fitter_test_suite<float>::fit_airfoil_transform_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_adaptive_airfoil_fitter_test_suite<double>::fit_airfoil_simple_test);
      TEST_ADD(piecewise_adaptive_airfoil_fitter_test_suite<double>::fit_airfoil_transform_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_adaptive_airfoil_fitter_test_suite<long double>::fit_airfoil_simple_test);
      TEST_ADD(piecewise_adaptive_airfoil_fitter_test_suite<long double>::fit_airfoil_transform_test);
    }

  public:
    piecewise_adaptive_airfoil_fitter_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_adaptive_airfoil_fitter_test_suite()
    {
    }

  private:

    void fit_airfoil_simple_test()
    {
      {
        airfoil_fitter_type pcaf;
        piecewise_curve_type af;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        point_type pt, pt_ref;
        data_type t0, t1, t2, angle, scale_factor;
        data_type fit_tol;
        bool rtn_flag;

        // get airfoil points
        create_airfoil_points(upt, lpt, true);

        // set the parameterization
        t0=0;
        t1=2;
        t2=4;

        // create curve
        subdivide_method sm(airfoil_fitter_type::BISECTION);
//        subdivide_method sm(airfoil_fitter_type::MAX_ERROR);

        fit_tol=1e-4;
        rtn_flag=pcaf.set_conditions(upt.begin(), static_cast<index_type>(upt.size()),
                                     lpt.begin(), static_cast<index_type>(lpt.size()),
                                     sm, fit_tol, false);
        TEST_ASSERT(rtn_flag);
        pcaf.set_t0(t0);
        pcaf.set_segment_dt(t1-t0, 0);
        pcaf.set_segment_dt(t2-t1, 1);
        rtn_flag=pcaf.create(af);
        TEST_ASSERT(rtn_flag);

#if 0
        // make sure got back correct angle and scale factor
        pt_ref.setZero();
        TEST_ASSERT(pt==pt_ref);
        TEST_ASSERT(tol.approximately_equal(angle, 0));
        TEST_ASSERT(tol.approximately_equal(scale_factor, 1));

        // cycle through CST control points to compare to the reference values
        cst_airfoil_control_point_type cp;
        index_type i;

        for (i=0; i<=cst.upper_degree(); ++i)
        {
          cp = cst.get_upper_control_point(i);
          std::string str("Error for Upper Control Point "+std::to_string(i));
          TEST_ASSERT_MSG(tol.approximately_equal(cp, cpu_ref[i]), str.data());
//          std::cout << "diff=" << (cp-cpu_ref[i]).norm() << std::endl;
        }

        for (i=0; i<=cst.lower_degree(); ++i)
        {
          cp = cst.get_lower_control_point(i);
          std::string str("Error for Lower Control Point "+std::to_string(i));
          TEST_ASSERT_MSG(tol.approximately_equal(cp, cpl_ref[i]), str.data());
//          std::cout << "diff=" << (cp-cpu_ref[i]).norm() << std::endl;
        }
#endif

        if (typeid(data_type)==typeid(float))
        {
          std::cout.flush();
          eli::test::octave_start(1);
          // print out the upper surface points
          for (size_t n=0; n<upt.size(); ++n)
          {
            std::string name("upt"); name+=std::to_string(n);
            eli::test::octave_print(1, upt[n], name);
          }
          // print out the lower surface points
          for (size_t n=0; n<lpt.size(); ++n)
          {
            std::string name("lpt"); name+=std::to_string(n);
            eli::test::octave_print(1, lpt[n], name);
          }
          eli::test::octave_print(1, af, "af");
          eli::test::octave_finish(1, true);
        }
      }

#if 0
      // fit to a known CST airfoil shape
      {
        airfoil_fitter_type pcaf;
        piecewise_curve_type pc;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        point_type pt, pt_ref;
        data_type t0, t1, t2;
        index_type degu, degl;
        bool rtn_flag;

        // get airfoil points
        create_cst_airfoil_points(upt, lpt);

        // set the parameterization
        t0=-1;
        t1=0;
        t2=1;

        // create curve
        degu=7;
        degl=5;
        rtn_flag=pcaf.set_conditions(upt.begin(), static_cast<index_type>(upt.size()), degu,
                                     lpt.begin(), static_cast<index_type>(lpt.size()), degl, false);
        TEST_ASSERT(rtn_flag);
        pcaf.set_t0(t0);
        pcaf.set_segment_dt(t1-t0, 0);
        pcaf.set_segment_dt(t2-t1, 1);
        rtn_flag=pcaf.create(pc);
        TEST_ASSERT(rtn_flag);

        // cycle through segments to get the control points to compare to the reference values
        control_point_type cp[32], cp_ref[32];
        curve_type crv;
        index_type i, ii;

        cp_ref[ 0] << static_cast<data_type>(1.00000000000000000), static_cast<data_type>(-0.00500000000000000), 0;
        cp_ref[ 1] << static_cast<data_type>(0.84615384615384600), static_cast<data_type>( 0.00346153846153849), 0;
        cp_ref[ 2] << static_cast<data_type>(0.70512820512820500), static_cast<data_type>(-0.00288461538461537), 0;
        cp_ref[ 3] << static_cast<data_type>(0.57692307692307700), static_cast<data_type>(-0.00288461538461520), 0;
        cp_ref[ 4] << static_cast<data_type>(0.46153846153846200), static_cast<data_type>( 0.01111888111888110), 0;
        cp_ref[ 5] << static_cast<data_type>(0.35897435897435900), static_cast<data_type>( 0.02306915306915300), 0;
        cp_ref[ 6] << static_cast<data_type>(0.26923076923076900), static_cast<data_type>( 0.01935314685314690), 0;
        cp_ref[ 7] << static_cast<data_type>(0.19230769230769200), static_cast<data_type>( 0.00146270396270403), 0;
        cp_ref[ 8] << static_cast<data_type>(0.12820512820512800), static_cast<data_type>(-0.01944444444444440), 0;
        cp_ref[ 9] << static_cast<data_type>(0.07692307692307690), static_cast<data_type>(-0.03283216783216790), 0;
        cp_ref[10] << static_cast<data_type>(0.03846153846153850), static_cast<data_type>(-0.03445804195804200), 0;
        cp_ref[11] << static_cast<data_type>(0.01282051282051280), static_cast<data_type>(-0.02621794871794880), 0;
        cp_ref[12] << static_cast<data_type>(0.00000000000000000), static_cast<data_type>(-0.01307692307692310), 0;
        cp_ref[13] << static_cast<data_type>(0.00000000000000000), static_cast<data_type>( 0.00000000000000000), 0;
        cp_ref[14] << static_cast<data_type>(0.00000000000000000), static_cast<data_type>( 0.00000000000000000), 0;
        cp_ref[15] << static_cast<data_type>(0.00000000000000000), static_cast<data_type>( 0.01000000000000000), 0;
        cp_ref[16] << static_cast<data_type>(0.00735294117647059), static_cast<data_type>( 0.02005147058823530), 0;
        cp_ref[17] << static_cast<data_type>(0.02205882352941180), static_cast<data_type>( 0.02980147058823530), 0;
        cp_ref[18] << static_cast<data_type>(0.04411764705882350), static_cast<data_type>( 0.03889705882352940), 0;
        cp_ref[19] << static_cast<data_type>(0.07352941176470590), static_cast<data_type>( 0.04703054298642530), 0;
        cp_ref[20] << static_cast<data_type>(0.11029411764705900), static_cast<data_type>( 0.05398472850678730), 0;
        cp_ref[21] << static_cast<data_type>(0.15441176470588200), static_cast<data_type>( 0.05961337926779100), 0;
        cp_ref[22] << static_cast<data_type>(0.20588235294117600), static_cast<data_type>( 0.06369210201563150), 0;
        cp_ref[23] << static_cast<data_type>(0.26470588235294100), static_cast<data_type>( 0.06573323735088460), 0;
        cp_ref[24] << static_cast<data_type>(0.33088235294117600), static_cast<data_type>( 0.06517508227067060), 0;
        cp_ref[25] << static_cast<data_type>(0.40441176470588200), static_cast<data_type>( 0.06229920814479620), 0;
        cp_ref[26] << static_cast<data_type>(0.48529411764705900), static_cast<data_type>( 0.05873642533936640), 0;
        cp_ref[27] << static_cast<data_type>(0.57352941176470600), static_cast<data_type>( 0.05454411764705930), 0;
        cp_ref[28] << static_cast<data_type>(0.66911764705882300), static_cast<data_type>( 0.04503676470588240), 0;
        cp_ref[29] << static_cast<data_type>(0.77205882352941200), static_cast<data_type>( 0.03525735294117730), 0;
        cp_ref[30] << static_cast<data_type>(0.88235294117647100), static_cast<data_type>( 0.02264705882353050), 0;
        cp_ref[31] << static_cast<data_type>(1.00000000000000000), static_cast<data_type>( 0.00700000000000000), 0;

        pc.get(crv, 0);
        ii=0;
        for (i=0; i<=crv.degree(); ++i, ++ii)
        {
          cp[ii]=crv.get_control_point(i);
        }
        pc.get(crv, 1);
        for (i=0; i<=crv.degree(); ++i, ++ii)
        {
          cp[ii]=crv.get_control_point(i);
        }

        for (i=0; i<32; ++i)
        {
          std::string str("Error for Control Point "+std::to_string(i));
          TEST_ASSERT_MSG(tol.approximately_equal(cp[i], cp_ref[i]), str.data());
//          std::cout << "cp_ref[" << std::setw(2) << i << "] << " << std::setprecision(15)
//                    << cp[i].x() << ", " << cp[i].y() << ", " << cp[i].z() << ";" << std::endl;
        }

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          // print out the upper surface points
//          for (size_t n=0; n<upt.size(); ++n)
//          {
//            std::string name("upt"); name+=std::to_string(n);
//            eli::test::octave_print(1, upt[n], name);
//          }
//          // print out the lower surface points
//          for (size_t n=0; n<lpt.size(); ++n)
//          {
//            std::string name("lpt"); name+=std::to_string(n);
//            eli::test::octave_print(1, lpt[n], name);
//          }
//          eli::test::octave_print(1, pc, "piecewise");
//          eli::test::octave_finish(1, true);
//        }
      }
#endif
    }

    void fit_airfoil_transform_test()
    {
#if 0
      // fit to simple two vector specification
      {
        airfoil_fitter_type pcaf;
        piecewise_curve_type pc;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        point_type pt, pt_ref;
        data_type t, t0, t1, t2;
        index_type degu, degl;
        bool rtn_flag;

        // get airfoil points
        create_airfoil_points(upt, lpt);

        // set the parameterization
        t0=-1;
        t1=0;
        t2=1;

        // create curve
        degu=8;
        degl=8;
        rtn_flag=pcaf.set_conditions(upt.begin(), static_cast<index_type>(upt.size()), degu,
                                     lpt.begin(), static_cast<index_type>(lpt.size()), degl, true);
        TEST_ASSERT(rtn_flag);
        pcaf.set_t0(t0);
        pcaf.set_segment_dt(t1-t0, 0);
        pcaf.set_segment_dt(t2-t1, 1);
        rtn_flag=pcaf.create(pc);
        TEST_ASSERT(rtn_flag);

        // test various points
        t  =-0.8;
        pt = pc.f(t);
        pt_ref << t*t, -0.05448192672, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);
        t  =-0.1;
        pt = pc.f(t);
        pt_ref << t*t, -0.0121366373, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);
        t  = 0.1;
        pt = pc.f(t);
        pt_ref << t*t, 0.01311318449, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);
        t  = 0.95;
        pt = pc.f(t);
        pt_ref << t*t, 0.01116413716, 0;
        TEST_ASSERT((pt-pt_ref).norm()<2e-4);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          // print out the upper surface points
//          for (size_t n=0; n<upt.size(); ++n)
//          {
//            std::string name("upt"); name+=std::to_string(n);
//            eli::test::octave_print(1, upt[n], name);
//          }
//          // print out the lower surface points
//          for (size_t n=0; n<lpt.size(); ++n)
//          {
//            std::string name("lpt"); name+=std::to_string(n);
//            eli::test::octave_print(1, lpt[n], name);
//          }
//          eli::test::octave_print(1, pc, "piecewise");
//          eli::test::octave_finish(1, true);
//        }
      }

      // scale, translate, and rotate points
      {
        airfoil_fitter_type pcaf;
        piecewise_curve_type pc;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        point_type pt, pt_ref, lept, tept;
        data_type t, t0, t1, t2;
        index_type degu, degl;
        bool rtn_flag;

        // get airfoil points
        create_airfoil_points(upt, lpt, lept, tept);

        // set the parameterization
        t0=-1;
        t1=0;
        t2=1;

        // create curve
        degu=8;
        degl=8;
        rtn_flag=pcaf.set_conditions(upt.begin(), static_cast<index_type>(upt.size()), degu,
                                     lpt.begin(), static_cast<index_type>(lpt.size()), degl, true);
        TEST_ASSERT(rtn_flag);
        pcaf.set_t0(t0);
        pcaf.set_segment_dt(t1-t0, 0);
        pcaf.set_segment_dt(t2-t1, 1);
        rtn_flag=pcaf.create(pc);
        TEST_ASSERT(rtn_flag);


        // test various points
        t  =-0.8;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 2.45614305445983, 2.67788624658704, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);

        t  =-0.1;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 1.03685893171540, 1.98615845755622, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);

        t  = 0.1;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 1.00525915448552, 2.04089087722623, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);

        t  = 0.95;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 2.93983922906606, 3.15259989674235, 0;
        TEST_ASSERT((pt-pt_ref).norm()<2e-4);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          // print out the upper surface points
//          for (size_t n=0; n<upt.size(); ++n)
//          {
//            std::string name("upt"); name+=std::to_string(n);
//            eli::test::octave_print(1, upt[n], name);
//          }
//          // print out the lower surface points
//          for (size_t n=0; n<lpt.size(); ++n)
//          {
//            std::string name("lpt"); name+=std::to_string(n);
//            eli::test::octave_print(1, lpt[n], name);
//          }
//          eli::test::octave_print(1, pc, "piecewise");
//          eli::test::octave_finish(1, true);
//        }
      }
#endif
    }
};

#endif
