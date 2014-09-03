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

#ifndef polynomial_curve_test_suite_hpp
#define polynomial_curve_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/curve/length.hpp"
#include "eli/geom/curve/curvature.hpp"
#include "eli/geom/curve/pseudo/polynomial.hpp"

#include "octave_helpers.hpp"

template<typename data__>
class polynomial_curve_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::pseudo::polynomial<data__, 3> curve_type;
    typedef typename curve_type::point_type point_type;
    typedef typename curve_type::data_type data_type;
    typedef typename curve_type::index_type index_type;
    typedef typename curve_type::coefficient_type coefficient_type;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(polynomial_curve_test_suite<float>::assignment_test);
      TEST_ADD(polynomial_curve_test_suite<float>::evaluation_test);
      TEST_ADD(polynomial_curve_test_suite<float>::derivative_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(polynomial_curve_test_suite<double>::assignment_test);
      TEST_ADD(polynomial_curve_test_suite<double>::evaluation_test);
      TEST_ADD(polynomial_curve_test_suite<double>::derivative_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(polynomial_curve_test_suite<long double>::assignment_test);
      TEST_ADD(polynomial_curve_test_suite<long double>::evaluation_test);
      TEST_ADD(polynomial_curve_test_suite<long double>::derivative_test);
    }

  public:
    polynomial_curve_test_suite()
    {
      AddTests(data__());
    }
    ~polynomial_curve_test_suite()
    {
    }

  private:
    void octave_print(int figno, const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts, const curve_type &bez) const
    {
      size_t i;

      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << "xpts=[" << pts[0].x();
      for (i=1; i<pts.size(); ++i)
        std::cout << ", " << pts[i].x();
      std::cout << "];" << std::endl;
      std::cout << "ypts=[" << pts[0].y();
      for (i=1; i<pts.size(); ++i)
        std::cout << ", " << pts[i].y();
      std::cout << "];" << std::endl;

      std::vector<data_type> t(101);
      for (i=0; i<t.size(); ++i)
        t[i]=static_cast<data_type>(i)/(t.size()-1);

      std::cout << "xint=[" << bez.f(t[0])(0);
      for (i=1; i<t.size(); ++i)
        std::cout << ", " << bez.f(t[i])(0);
      std::cout << "];" << std::endl;
      std::cout << "yint=[" << bez.f(t[0])(1);
      for (i=1; i<t.size(); ++i)
        std::cout << ", " << bez.f(t[i])(1);
      std::cout << "];" << std::endl;

      std::cout << "plot(xpts, ypts, 'bo', xint, yint, 'k-');" << std::endl;
    }

    void assignment_test()
    {
      curve_type c1, c2;
      coefficient_type coef1(6), coef2(1), coef_out;

      // test setting coefficients
      coef1 << 2, 0, 1, 4, 1, 8;
      coef2 << 0;
      c1.set_coefficients(coef1, 0);
      c1.set_coefficients(coef2, 1);
      c1.set_coefficients(coef2, 2);
      c1.get_coefficients(coef_out, 0);
      TEST_ASSERT(coef1==coef_out);
      coef_out.setZero();
      c1.get_coefficients(coef_out, 1);
      TEST_ASSERT(coef2==coef_out);
      coef_out.setZero();
      c1.get_coefficients(coef_out, 2);
      TEST_ASSERT(coef2==coef_out);
      coef_out.setZero();

      // test assignment operator
      c2=c1;
      c2.get_coefficients(coef_out, 0);
      TEST_ASSERT(coef1==coef_out);
      coef_out.setZero();
      c2.get_coefficients(coef_out, 1);
      TEST_ASSERT(coef2==coef_out);
      coef_out.setZero();
      c2.get_coefficients(coef_out, 2);
      TEST_ASSERT(coef2==coef_out);
      coef_out.setZero();

      // test equivalence operator
      TEST_ASSERT(c1==c2);
      coef2 << 1;
      c2.set_coefficients(coef2, 2);
      TEST_ASSERT(c1!=c2);

      // test copy ctr
      curve_type c3(c1);
      TEST_ASSERT(c3==c1);
    }

    void evaluation_test()
    {
      curve_type c;
      coefficient_type coef1(5), coef2(2);
      point_type eval_out, eval_ref;
      data_type t;

      // set coefficients
      coef1 << 2, 4, 3, 1, 2;
      coef2 << 1, 1;
      c.set_coefficients(coef1, 0);
      c.set_coefficients(coef2, 1);

      // test evaluation at points
      t=0;
      eval_out=c.f(t);
      eval_ref << coef1(0), coef2(0), 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=c.f(t);
      eval_ref << coef1.sum(), coef2.sum(), 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=10;
      eval_out=c.f(t);
      eval_ref << 21342, 11, 0;
      TEST_ASSERT(eval_out==eval_ref);

//      if (typeid(data_type)==typeid(float))
//      {
//        std::cout.flush();
//        eli::test::octave_start(1);
//        eli::test::octave_print(1, c, "poly");
//        eli::test::octave_finish(1);
//      }
    }

    void derivative_test()
    {
#if 0
      typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
      data_type eps(std::numeric_limits<data__>::epsilon());
      control_point_type cntrl_in[5];
      typename bezier_curve_type::control_point_type bez_cntrl[5];
      curve_type ebc;
      bezier_curve_type bc;
      typename curve_type::point_type eval_out, eval_ref;
      typename curve_type::data_type t;

      // set control points and create curves
      cntrl_in[0] << 2.0;
      cntrl_in[1] << 1.5;
      cntrl_in[2] << 0.0;
      cntrl_in[3] << 1.0;
      cntrl_in[4] << 0.5;
      bez_cntrl[0] << 0,    2;
      bez_cntrl[1] << 0.25, 1.5;
      bez_cntrl[2] << 0.5,  0;
      bez_cntrl[3] << 0.75, 1;
      bez_cntrl[4] << 1,    0.5;

      ebc.resize(4);
      bc.resize(4);
      for (index_type i=0; i<5; ++i)
      {
        ebc.set_control_point(cntrl_in[i], i);
        bc.set_control_point(bez_cntrl[i], i);
      }

      // test 1st derivative at end points
      t=0;
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 1st derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      if (typeid(data__)==typeid(float))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
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
      t=static_cast<data__>(0.45);
      eval_out=ebc.fpp(t);
      eval_ref=bc.fpp(t);
      if (typeid(data__)==typeid(float))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<17*eps);
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
      t=static_cast<data__>(0.45);
      eval_out=ebc.fppp(t);
      eval_ref=bc.fppp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test curvature at end points
      data_type curv_out, curv_ref;
      t=0;
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      TEST_ASSERT(curv_out==curv_ref);
      t=1;
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      TEST_ASSERT(curv_out==curv_ref);

      // test curvature at interior point
      t=static_cast<data__>(0.45);
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      if (typeid(data__)==typeid(float))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
      }
      else
      {
        TEST_ASSERT(curv_out==curv_ref);
      }
#endif
    }
};
#endif

