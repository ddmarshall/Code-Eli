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
#include "eli/geom/curve/pseudo/cst_airfoil.hpp"

template<typename data__>
class piecewise_cst_airfoil_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::pseudo::cst_airfoil<data_type> cst_airfoil_type;
    typedef typename cst_airfoil_type::control_point_type cst_airfoil_control_point_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<float>::create_airfoil_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<double>::create_airfoil_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<long double>::create_airfoil_test);
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
//        eli::test::octave_finish(1);
//      }
    }
};

#endif
