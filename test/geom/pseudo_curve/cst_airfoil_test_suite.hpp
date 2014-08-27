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

#ifndef cst_airfoil_curve_test_suite_hpp
#define cst_airfoil_curve_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"
#include "eli/mutil/fd/d1o2.hpp"
#include "eli/mutil/fd/d2o2.hpp"

#include "eli/geom/point/distance.hpp"

#include "eli/geom/curve/length.hpp"
#include "eli/geom/curve/curvature.hpp"

#include "eli/geom/curve/pseudo/cst_airfoil.hpp"

#include "octave_helpers.hpp"

template<typename data__>
class cst_airfoil_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::cst_airfoil<data__> curve_type;
    typedef typename curve_type::point_type point_type;
    typedef typename curve_type::data_type data_type;
    typedef typename curve_type::index_type index_type;
    typedef typename curve_type::control_point_type control_point_type;
    typedef typename curve_type::tolerance_type tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(cst_airfoil_test_suite<float>::assignment_test);
      TEST_ADD(cst_airfoil_test_suite<float>::evaluation_test);
      TEST_ADD(cst_airfoil_test_suite<float>::derivative_test);
//      TEST_ADD(cst_airfoil_test_suite<float>::promotion_test);
//      TEST_ADD(cst_airfoil_test_suite<float>::demotion_test);
//      TEST_ADD(cst_airfoil_test_suite<float>::length_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(cst_airfoil_test_suite<double>::assignment_test);
      TEST_ADD(cst_airfoil_test_suite<double>::evaluation_test);
      TEST_ADD(cst_airfoil_test_suite<double>::derivative_test);
//      TEST_ADD(cst_airfoil_test_suite<double>::promotion_test);
//      TEST_ADD(cst_airfoil_test_suite<double>::demotion_test);
//      TEST_ADD(cst_airfoil_test_suite<double>::length_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(cst_airfoil_test_suite<long double>::assignment_test);
      TEST_ADD(cst_airfoil_test_suite<long double>::evaluation_test);
      TEST_ADD(cst_airfoil_test_suite<long double>::derivative_test);
//      TEST_ADD(cst_airfoil_test_suite<long double>::promotion_test);
//      TEST_ADD(cst_airfoil_test_suite<long double>::demotion_test);
//      TEST_ADD(cst_airfoil_test_suite<long double>::length_test);
    }

  public:
    cst_airfoil_test_suite()
    {
      AddTests(data__());
    }
    ~cst_airfoil_test_suite()
    {
    }

  private:
    void assignment_test()
    {
      curve_type ebc1, ebc2;
      control_point_type cntrl_in[4], cntrl_out;

      // test default constructor then set control points
      cntrl_in[0] << 2;
      cntrl_in[1] << 1;
      cntrl_in[2] << static_cast<data_type>(3.5);
      cntrl_in[3] << 4;

      ebc1.resize(3);
      for (index_type i=0; i<=ebc1.degree(); ++i)
      {
        ebc1.set_control_point(cntrl_in[i], i);
      }

      // test order
      TEST_ASSERT(ebc1.degree()==3);

      // test control points
      for (index_type i=0; i<=ebc1.degree(); ++i)
      {
        TEST_ASSERT(tol.approximately_equal(ebc1.get_control_point(i), cntrl_in[i]));
      }

      // test copy ctr
      curve_type ebcc2(ebc1);
      TEST_ASSERT(ebcc2==ebc1);

      // test assignment operator
      ebc2=ebc1;
      TEST_ASSERT(ebc2==ebc1);
    }

    void evaluation_test()
    {
      typedef eli::geom::curve::pseudo::explicit_bezier<data_type> explicit_bezier_curve_type;

      explicit_bezier_curve_type ebc(7);
      curve_type cst(7);
      control_point_type cp[8];
      data_type dte(0.01), t[6];
      point_type pt_out, pt_ref;
      index_type i;

      // set the control points
      cp[0] << static_cast<data_type>(0.170987592880629);
      cp[1] << static_cast<data_type>(0.157286894410384);
      cp[2] << static_cast<data_type>(0.162311658384540);
      cp[3] << static_cast<data_type>(0.143623187913493);
      cp[4] << static_cast<data_type>(0.149218456400780);
      cp[5] << static_cast<data_type>(0.137218405082418);
      cp[6] << static_cast<data_type>(0.140720628655908);
      cp[7] << static_cast<data_type>(0.141104769355436);
      for (i=0; i<=cst.degree(); ++i)
      {
        ebc.set_control_point(cp[i], i);
        cst.set_control_point(cp[i], i);
      }

      // set the trailing edge thickness of CST airfoil
      cst.set_trailing_edge_thickness(dte);

      // set the parameters to evaluate the tests
      t[0] = static_cast<data_type>(0);
      t[1] = static_cast<data_type>(0.1);
      t[2] = static_cast<data_type>(0.27);
      t[3] = static_cast<data_type>(0.5);
      t[4] = static_cast<data_type>(0.73);
      t[5] = static_cast<data_type>(1);

      // evaluate the point
      point_type dte_vec;
      i=0;
      pt_out = cst.f(t[i]);
      dte_vec << 0, t[i]*dte;
      pt_ref = std::sqrt(t[i])*(1-t[i])*ebc.f(t[i])+dte_vec;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=1;
      pt_out = cst.f(t[i]);
      dte_vec << 0, t[i]*dte;
      pt_ref = std::sqrt(t[i])*(1-t[i])*ebc.f(t[i])+dte_vec;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=2;
      pt_out = cst.f(t[i]);
      dte_vec << 0, t[i]*dte;
      pt_ref = std::sqrt(t[i])*(1-t[i])*ebc.f(t[i])+dte_vec;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=3;
      pt_out = cst.f(t[i]);
      dte_vec << 0, t[i]*dte;
      pt_ref = std::sqrt(t[i])*(1-t[i])*ebc.f(t[i])+dte_vec;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=4;
      pt_out = cst.f(t[i]);
      dte_vec << 0, t[i]*dte;
      pt_ref = std::sqrt(t[i])*(1-t[i])*ebc.f(t[i])+dte_vec;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=5;
      pt_out = cst.f(t[i]);
      dte_vec << 0, t[i]*dte;
      pt_ref = std::sqrt(t[i])*(1-t[i])*ebc.f(t[i])+dte_vec;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
    }

    void derivative_test()
    {
      typedef eli::geom::curve::pseudo::explicit_bezier<data_type> explicit_bezier_curve_type;

      explicit_bezier_curve_type ebc(7);
      curve_type cst(7);
      control_point_type cp[8];
      data_type dte(0.01), t[6];
      point_type pt, xp, xp_ref, xpp, xpp_ref;
      index_type i;

      // set the control points
      cp[0] << static_cast<data_type>(0.170987592880629);
      cp[1] << static_cast<data_type>(0.157286894410384);
      cp[2] << static_cast<data_type>(0.162311658384540);
      cp[3] << static_cast<data_type>(0.143623187913493);
      cp[4] << static_cast<data_type>(0.149218456400780);
      cp[5] << static_cast<data_type>(0.137218405082418);
      cp[6] << static_cast<data_type>(0.140720628655908);
      cp[7] << static_cast<data_type>(0.141104769355436);
      for (i=0; i<=cst.degree(); ++i)
      {
        ebc.set_control_point(cp[i], i);
        cst.set_control_point(cp[i], i);
      }

      // set the trailing edge thickness of CST airfoil
      cst.set_trailing_edge_thickness(dte);

      // set the parameters to evaluate the tests
      t[0] = static_cast<data_type>(0.1);
      t[1] = static_cast<data_type>(0.24);
      t[2] = static_cast<data_type>(0.42);
      t[3] = static_cast<data_type>(0.5);
      t[4] = static_cast<data_type>(0.73);
      t[5] = static_cast<data_type>(0.95);

      data_type x[3], y[3], dt(static_cast<data_type>(1e2)*std::sqrt(std::numeric_limits<data_type>::epsilon()));
      eli::mutil::fd::d1o2<data_type> d1_calc;
      eli::mutil::fd::d2o2<data_type> d2_calc;

      for (i=0; i<6; ++i)
      {
        pt=cst.f(t[i]-dt); x[0]=pt(0); y[0]=pt(1);
        pt=cst.f(t[i]);    x[1]=pt(0); y[1]=pt(1);
        pt=cst.f(t[i]+dt); x[2]=pt(0); y[2]=pt(1);
        d1_calc.evaluate(xp_ref(0), x, dt);
        d1_calc.evaluate(xp_ref(1), y, dt);
        d2_calc.evaluate(xpp_ref(0), x, dt);
        d2_calc.evaluate(xpp_ref(1), y, dt);
        xp=cst.fp(t[i]);
        xpp=cst.fpp(t[i]);
        if (typeid(data_type)==typeid(float))
        {
          TEST_ASSERT((xp-xp_ref).norm()<6e-3);
          if (i==0)
          {
            TEST_ASSERT((xpp-xpp_ref).norm()<7e-2);
          }
          else if (i==1)
          {
            TEST_ASSERT((xpp-xpp_ref).norm()<4e-3);
          }
          else
          {
            TEST_ASSERT((xpp-xpp_ref).norm()<3e-3);
          }
        }
        else
        {
          TEST_ASSERT(tol.approximately_equal(xp, xp_ref));
          TEST_ASSERT((xpp-xpp_ref).norm()<3e-5);
        }
      }
    }

    void promotion_test()
    {
      typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
      data_type eps(std::numeric_limits<data__>::epsilon());

      control_point_type cntrl_in(5, 1);
      typename bezier_curve_type::control_point_type bez_cntrl(5,2);
      curve_type ebc;
      bezier_curve_type bc;
      typename curve_type::point_type eval_out, eval_ref;
      typename curve_type::data_type t, curv_out, curv_ref;

      // set control points and create curves
      cntrl_in << 2.0,
                  1.5,
                  0.0,
                  1.0,
                  0.5;
      bez_cntrl.col(0) << 0, 0.25, 0.5, 0.75, 1;
      bez_cntrl.col(1)=cntrl_in;
      ebc.set_control_points(cntrl_in);
      bc.set_control_points(bez_cntrl);

      ebc.degree_promote();
      bc.degree_promote();

      // test to see if degree has increased
      TEST_ASSERT(ebc.degree()==bc.degree());

      // test evaluation at end points
      t=0;
      eval_out=ebc.f(t);
      eval_ref=bc.f(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.f(t);
      eval_ref=bc.f(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test evaluation at interior point
      t=0.45;
      eval_out=ebc.f(t);
      eval_ref=bc.f(t);
      if (typeid(data_type)==typeid(double))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<2.1*eps);
      }
      else
      {
        TEST_ASSERT(eval_out==eval_ref);
      }

      // test 1st derivative at end points
      t=0;
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      if (typeid(data_type)==typeid(double))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<24.1*eps);
      }
      else
      {
        TEST_ASSERT(eval_out==eval_ref);
      }

      // test 1st derivative at interior point
      t=0.45;
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      if ((typeid(data_type)==typeid(float)) || (typeid(data_type)==typeid(double)))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<4.1*eps);
      }
      else
      {
        TEST_ASSERT(eval_out==eval_ref);
      }

      // test 2nd derivative at end points
      t=0;
      eval_out=ebc.fpp(t);
      eval_ref=bc.fpp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.fpp(t);
      eval_ref=bc.fpp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 2nd derivative at interior point
      t=0.45;
      eval_out=ebc.fpp(t);
      eval_ref=bc.fpp(t);
      if (typeid(data_type)==typeid(double))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<9.1*eps);
      }
      else
      {
        TEST_ASSERT(eval_out==eval_ref);
      }

      // test 3rd derivative at end points
      t=0;
      eval_out=ebc.fppp(t);
      eval_ref=bc.fppp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.fppp(t);
      eval_ref=bc.fppp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 3rd derivative at interior point
      t=0.45;
      eval_out=ebc.fppp(t);
      eval_ref=bc.fppp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test curvature at end points
      t=0;
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      TEST_ASSERT(curv_out==curv_ref);
      t=1;
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      if (typeid(data_type)==typeid(double))
      {
        TEST_ASSERT(std::abs(curv_out-curv_ref)<47.1*eps);
      }
      else
      {
        TEST_ASSERT(curv_out==curv_ref);
      }

      // test curvature at interior point
      t=0.45;
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      if ((typeid(data_type)==typeid(float)) || (typeid(data_type)==typeid(double)))
      {
        TEST_ASSERT(std::abs(curv_out-curv_ref)<6.1*eps);
      }
      else
      {
        TEST_ASSERT(curv_out==curv_ref);
      }
    }

    void demotion_test()
    {
      // no constraint
      {
        typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;

        control_point_type cntrl_in(9, 1), cntrl_out;
        typename bezier_curve_type::control_point_type bez_cntrl(9,2), bez_cntrl_out;
        curve_type ebc;
        bezier_curve_type bc;

        // set control points and create curves
        cntrl_in << 2.0,
                    1.5,
                    1.0,
                    0.5,
                    0.0,
                   -0.5,
                   -1.0,
                    1.0,
                    0.5;
        bez_cntrl.col(0) << 0, 0.125, 0.25, 0.325, 0.5, 0.625, 0.75, 0.875, 1;
        bez_cntrl.col(1)=cntrl_in;
        ebc.set_control_points(cntrl_in);
        bc.set_control_points(bez_cntrl);

        ebc.degree_demote(eli::geom::general::NOT_CONNECTED);
        bc.degree_demote(eli::geom::general::NOT_CONNECTED);
        bc.get_control_points(bez_cntrl_out);
        ebc.get_control_points(cntrl_out);
        TEST_ASSERT(bez_cntrl_out.col(1)==cntrl_out);
      }

      // C0 constraint
      {
        typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;

        control_point_type cntrl_in(9, 1), cntrl_out;
        typename bezier_curve_type::control_point_type bez_cntrl(9,2), bez_cntrl_out;
        curve_type ebc;
        bezier_curve_type bc;

        // set control points and create curves
        cntrl_in << 2.0,
                    1.5,
                    1.0,
                    0.5,
                    0.0,
                   -0.5,
                   -1.0,
                    1.0,
                    0.5;
        bez_cntrl.col(0) << 0, 0.125, 0.25, 0.325, 0.5, 0.625, 0.75, 0.875, 1;
        bez_cntrl.col(1)=cntrl_in;
        ebc.set_control_points(cntrl_in);
        bc.set_control_points(bez_cntrl);

        ebc.degree_demote(eli::geom::general::C0);
        bc.degree_demote(eli::geom::general::C0);
        bc.get_control_points(bez_cntrl_out);
        ebc.get_control_points(cntrl_out);
        TEST_ASSERT(bez_cntrl_out.col(1)==cntrl_out);
      }

      // C1 constraint
      {
        typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;

        control_point_type cntrl_in(9, 1), cntrl_out;
        typename bezier_curve_type::control_point_type bez_cntrl(9,2), bez_cntrl_out;
        curve_type ebc;
        bezier_curve_type bc;

        // set control points and create curves
        cntrl_in << 2.0,
                    1.5,
                    1.0,
                    0.5,
                    0.0,
                   -0.5,
                   -1.0,
                    1.0,
                    0.5;
        bez_cntrl.col(0) << 0, 0.125, 0.25, 0.325, 0.5, 0.625, 0.75, 0.875, 1;
        bez_cntrl.col(1)=cntrl_in;
        ebc.set_control_points(cntrl_in);
        bc.set_control_points(bez_cntrl);

        ebc.degree_demote(eli::geom::general::C1);
        bc.degree_demote(eli::geom::general::C1);
        bc.get_control_points(bez_cntrl_out);
        ebc.get_control_points(cntrl_out);
        TEST_ASSERT(bez_cntrl_out.col(1)==cntrl_out);
      }

      // C2 constraint
      {
        typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;

        control_point_type cntrl_in(9, 1), cntrl_out;
        typename bezier_curve_type::control_point_type bez_cntrl(9,2), bez_cntrl_out;
        curve_type ebc;
        bezier_curve_type bc;

        // set control points and create curves
        cntrl_in << 2.0,
                    1.5,
                    1.0,
                    0.5,
                    0.0,
                   -0.5,
                   -1.0,
                    1.0,
                    0.5;
        bez_cntrl.col(0) << 0, 0.125, 0.25, 0.325, 0.5, 0.625, 0.75, 0.875, 1;
        bez_cntrl.col(1)=cntrl_in;
        ebc.set_control_points(cntrl_in);
        bc.set_control_points(bez_cntrl);

        ebc.degree_demote(eli::geom::general::C2);
        bc.degree_demote(eli::geom::general::C2);
        bc.get_control_points(bez_cntrl_out);
        ebc.get_control_points(cntrl_out);
        TEST_ASSERT(bez_cntrl_out.col(1)==cntrl_out);
      }
    }

    void length_test()
    {
      typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
      data_type eps(std::numeric_limits<data__>::epsilon());

      control_point_type cntrl_in(5, 1);
      typename bezier_curve_type::control_point_type bez_cntrl(5,2);
      curve_type ebc;
      bezier_curve_type bc;
      typename curve_type::data_type length_cal, length_ref;

      // set control points and create curves
      cntrl_in << 2.0,
                  1.5,
                  0.0,
                  1.0,
                  0.5;
      bez_cntrl.col(0) << 0, 0.25, 0.5, 0.75, 1;
      bez_cntrl.col(1)=cntrl_in;
      ebc.set_control_points(cntrl_in);
      bc.set_control_points(bez_cntrl);

      // calculate the length of curve
      data_type tol(std::sqrt(eps));
      length(length_cal, ebc, tol);
      length(length_ref, bc, tol);
      TEST_ASSERT(length_cal==length_ref);

      // test computing some segment length
      typename curve_type::data_type t0(0.2), t1(0.7);

      length(length_cal, ebc, t0, t1, tol);
      length(length_ref, bc, t0, t1, tol);
      TEST_ASSERT(length_cal==length_ref);
    }
};
#endif

