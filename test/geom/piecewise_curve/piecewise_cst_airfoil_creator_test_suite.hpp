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

#ifndef piecewise_cst_airfoil_creator_test_suite_hpp
#define piecewise_cst_airfoil_creator_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"
#include "eli/mutil/fd/d1o2.hpp"
#include "eli/mutil/fd/d2o2.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_cst_airfoil_creator.hpp"
#include "eli/geom/curve/piecewise_cst_airfoil_fitter.hpp"
#include "eli/geom/curve/pseudo/cst_airfoil.hpp"

template<typename data__>
class piecewise_cst_airfoil_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::control_point_type control_point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::pseudo::cst_airfoil<data_type> cst_airfoil_type;
    typedef typename cst_airfoil_type::control_point_type cst_airfoil_control_point_type;

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

    void create_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &upts,
                               std::vector<point_type, Eigen::aligned_allocator<point_type>> &lpts)
    {
      // upper surface
      upts.resize(24);
      upts[ 0] << 0.0000, 0.0000, 0;
      upts[ 1] << 0.0013, 0.0030, 0;
      upts[ 2] << 0.0045, 0.0094, 0;
      upts[ 3] << 0.0096, 0.0138, 0;
      upts[ 4] << 0.0165, 0.0183, 0;
      upts[ 5] << 0.0252, 0.0228, 0;
      upts[ 6] << 0.0477, 0.0315, 0;
      upts[ 7] << 0.0767, 0.0397, 0;
      upts[ 8] << 0.1118, 0.0473, 0;
      upts[ 9] << 0.1524, 0.0539, 0;
      upts[10] << 0.1980, 0.0593, 0;
      upts[11] << 0.2979, 0.0636, 0;
      upts[12] << 0.3015, 0.0665, 0;
      upts[13] << 0.3578, 0.0680, 0;
      upts[14] << 0.4160, 0.0680, 0;
      upts[15] << 0.4455, 0.0675, 0;
      upts[16] << 0.5049, 0.0652, 0;
      upts[17] << 0.5930, 0.0585, 0;
      upts[18] << 0.6501, 0.0514, 0;
      upts[19] << 0.7050, 0.0416, 0;
      upts[20] << 0.7623, 0.0297, 0;
      upts[21] << 0.8168, 0.0221, 0;
      upts[22] << 0.9074, 0.0108, 0;
      upts[23] << 1.0000, 0.0050, 0;

      // lower surface
      lpts.resize(17);
      lpts[ 0] << 0.0000, 0.0000, 0;
      lpts[ 1] << 0.0029,-0.0077, 0;
      lpts[ 2] << 0.0076,-0.0116, 0;
      lpts[ 3] << 0.0143,-0.0112, 0;
      lpts[ 4] << 0.0451,-0.0284, 0;
      lpts[ 5] << 0.0741,-0.0365, 0;
      lpts[ 6] << 0.1029,-0.0441, 0;
      lpts[ 7] << 0.1499,-0.0511, 0;
      lpts[ 8] << 0.2453,-0.0621, 0;
      lpts[ 9] << 0.2986,-0.0657, 0;
      lpts[10] << 0.3545,-0.0678, 0;
      lpts[11] << 0.4123,-0.0684, 0;
      lpts[12] << 0.6178,-0.0572, 0;
      lpts[13] << 0.7017,-0.0453, 0;
      lpts[14] << 0.8078,-0.0247, 0;
      lpts[15] << 0.9012,-0.0087, 0;
      lpts[16] << 1.0000,-0.0020, 0;
    }

    void create_cst_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &upts,
                                   std::vector<point_type, Eigen::aligned_allocator<point_type>> &lpts)
    {
      // create an asymmetric CST airfoil
      typedef eli::geom::curve::pseudo::cst_airfoil<data__> cst_airfoil_type;
      typedef typename cst_airfoil_type::point_type cst_airfoil_point_type;
      typedef typename cst_airfoil_type::data_type cst_airfoil_data_type;
      typedef typename cst_airfoil_type::index_type cst_airfoil_index_type;
      typedef typename cst_airfoil_type::control_point_type cst_airfoil_control_point_type;

      cst_airfoil_type cst(7, 5);
      cst_airfoil_control_point_type cpu[8], cpl[6];
      cst_airfoil_data_type dteu(0.007), dtel(0.005);
      cst_airfoil_index_type i;

      // set the control points
      cpu[0] << static_cast<data_type>( 0.17);
      cpu[1] << static_cast<data_type>( 0.16);
      cpu[2] << static_cast<data_type>( 0.16);
      cpu[3] << static_cast<data_type>( 0.14);
      cpu[4] << static_cast<data_type>( 0.15);
      cpu[5] << static_cast<data_type>( 0.14);
      cpu[6] << static_cast<data_type>( 0.14);
      cpu[7] << static_cast<data_type>( 0.14);
      cpl[0] << static_cast<data_type>(-0.17);
      cpl[1] << static_cast<data_type>( 0.08);
      cpl[2] << static_cast<data_type>( 0.06);
      cpl[3] << static_cast<data_type>( 0.00);
      cpl[4] << static_cast<data_type>( 0.00);
      cpl[5] << static_cast<data_type>( 0.05);
      for (i=0; i<=cst.upper_degree(); ++i)
      {
        cst.set_upper_control_point(cpu[i], i);
      }
      for (i=0; i<=cst.lower_degree(); ++i)
      {
        cst.set_lower_control_point(cpl[i], i);
      }

      // set the trailing edge thickness of CST airfoil
      cst.set_trailing_edge_thickness(dteu, dtel);

      // sample the lower and upper surfaces
      upts.resize(24);
      for (i=0; i<upts.size(); ++i)
      {
        cst_airfoil_point_type pt(cst.f(static_cast<data_type>(i)/(upts.size()-1)));
        upts[i] << pt.x(), pt.y(), 0;
      }
      lpts.resize(20);
      for (i=0; i<lpts.size(); ++i)
      {
        cst_airfoil_point_type pt(cst.f(-static_cast<data_type>(i)/(lpts.size()-1)));
        lpts[i] << pt.x(), pt.y(), 0;
      }
    }

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<float>::create_airfoil_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<float>::fit_airfoil_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<double>::create_airfoil_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<double>::fit_airfoil_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<long double>::create_airfoil_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<long double>::fit_airfoil_test);
    }

  public:
    piecewise_cst_airfoil_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_cst_airfoil_creator_test_suite()
    {
    }

  private:

    void create_airfoil_test()
    {
      typedef eli::geom::curve::piecewise_cst_airfoil_creator<data__, 3, tolerance_type> airfoil_creator_type;

      airfoil_creator_type pcst;
      piecewise_curve_type pc;
      cst_airfoil_type cst(7);
      cst_airfoil_control_point_type cp[8];
      data_type dte(2*0.00126), ti, t0, t1, t2, t[6];
      point_type pt_out, pt_ref;
      typename cst_airfoil_type::point_type pt2_ref;
      index_type i;
      bool rtn_flag;

      // set the control points
      cp[0] << static_cast<data_type>(0.170987592880629);
      cp[1] << static_cast<data_type>(0.157286894410384);
      cp[2] << static_cast<data_type>(0.162311658384540);
      cp[3] << static_cast<data_type>(0.143623187913493);
      cp[4] << static_cast<data_type>(0.149218456400780);
      cp[5] << static_cast<data_type>(0.137218405082418);
      cp[6] << static_cast<data_type>(0.140720628655908);
      cp[7] << static_cast<data_type>(0.141104769355436);
      for (i=0; i<=cst.upper_degree(); ++i)
      {
        cst.set_upper_control_point(cp[i], i);
        cst.set_lower_control_point(-cp[i], i);
      }

      // set the trailing edge thickness of CST airfoil
      cst.set_trailing_edge_thickness(dte);

      // set the parameterization
      t0=-1;
      t1=0;
      t2=1;

      // set the parameters to evaluate the tests
      t[0] = t0+(t2-t0)*static_cast<data_type>(0);
      t[1] = t0+(t2-t0)*static_cast<data_type>(0.1);
      t[2] = t0+(t2-t0)*static_cast<data_type>(0.27);
      t[3] = t0+(t2-t0)*static_cast<data_type>(0.5);
      t[4] = t0+(t2-t0)*static_cast<data_type>(0.73);
      t[5] = t0+(t2-t0)*static_cast<data_type>(1);

      // create curve
      rtn_flag=pcst.set_conditions(cst);
      TEST_ASSERT(rtn_flag);
      pcst.set_t0(t0);
      pcst.set_segment_dt(t1-t0, 0);
      pcst.set_segment_dt(t2-t1, 1);
      rtn_flag=pcst.create(pc);
      TEST_ASSERT(rtn_flag);

      // evaluate the points (note need to transform parameterization to match points
      i=0;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=1;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=2;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=3;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=4;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=5;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));

//      if (typeid(data_type)==typeid(float))
//      {
//        std::cout.flush();
//        eli::test::octave_start(1);
//        eli::test::octave_print(1, cst, "cst");
//        eli::test::octave_print(1, pc, "piecewise");
//        eli::test::octave_finish(1, false);
//      }
    }

    void fit_airfoil_test()
    {
      typedef eli::geom::curve::piecewise_cst_airfoil_fitter<data__, 3, tolerance_type> airfoil_fitter_type;

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
        if (typeid(data_type)==typeid(float))
        {
          pt_ref << t*t, 0.01116413716, 0;
        }
        else
        {
          pt_ref << t*t, 0.01130447058, 0;
        }
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);

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
      }

      // fit to a known CST airfoil shape
      {
        airfoil_fitter_type pcaf;
        piecewise_curve_type pc;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        point_type pt, pt_ref;
        data_type t, t0, t1, t2;
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
    }
};

#endif
