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
//        eli::test::octave_finish(1, false);
//      }
    }

    void derivative_test()
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

      // test evaluation of f' at points
      t=0;
      eval_out=c.fp(t);
      eval_ref << coef1(1), coef2(1), 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=c.fp(t);
      eval_ref << (4+6+3+8), 1, 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=10;
      eval_out=c.fp(t);
      eval_ref << (4+60+300+8000), 1, 0;
      TEST_ASSERT(eval_out==eval_ref);

      // test evaluation of f'' at points
      t=0;
      eval_out=c.fpp(t);
      eval_ref << 2*coef1(2), 0, 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=c.fpp(t);
      eval_ref << (6+6+24), 0, 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=10;
      eval_out=c.fpp(t);
      eval_ref << (6+60+2400), 0, 0;
      TEST_ASSERT(eval_out==eval_ref);

      // test evaluation of f''' at points
      t=0;
      eval_out=c.fppp(t);
      eval_ref << 6*coef1(3), 0, 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=c.fppp(t);
      eval_ref << (6+48), 0, 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=10;
      eval_out=c.fppp(t);
      eval_ref << (6+480), 0, 0;
      TEST_ASSERT(eval_out==eval_ref);
    }
};
#endif

