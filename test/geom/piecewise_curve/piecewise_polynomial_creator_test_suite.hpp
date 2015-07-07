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

#ifndef piecewise_polynomial_creator_test_suite_hpp
#define piecewise_polynomial_creator_test_suite_hpp

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
#include "eli/geom/curve/piecewise_polynomial_creator.hpp"

#include "octave_helpers.hpp"

template<typename data__>
class piecewise_polynomial_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::pseudo::polynomial<data__, 3> polynomial_type;
    typedef typename polynomial_type::coefficient_type coefficient_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_polynomial_creator_test_suite<float>::create_curve_test);
      TEST_ADD(piecewise_polynomial_creator_test_suite<float>::full_cycle_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_polynomial_creator_test_suite<double>::create_curve_test);
      TEST_ADD(piecewise_polynomial_creator_test_suite<double>::full_cycle_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_polynomial_creator_test_suite<long double>::create_curve_test);
      TEST_ADD(piecewise_polynomial_creator_test_suite<long double>::full_cycle_test);
    }

  public:
    piecewise_polynomial_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_polynomial_creator_test_suite()
    {
    }

  private:

    void create_curve_test()
    {
      eli::geom::curve::piecewise_polynomial_creator<data_type, 3, tolerance_type> ppc;
      piecewise_curve_type pc;
      polynomial_type c;
      coefficient_type coef1(5), coef2(2);
      point_type eval_out, eval_ref;
      data_type t, t0, t1;
      bool rtn_flag;

      // set coefficients
      coef1 << 2, 4, 3, 1, 2;
      coef2 << 1, 1;
      c.set_coefficients(coef1, 0);
      c.set_coefficients(coef2, 1);

      // set the parameterization
      t0=-1;
      t1=1;

      // create curve
      rtn_flag=ppc.set_conditions(c);
      TEST_ASSERT(rtn_flag);
      ppc.set_t0(t0);
      ppc.set_segment_dt(t1-t0, 0);
      rtn_flag=ppc.create(pc);
      TEST_ASSERT(rtn_flag);

      // test evaluation at ends
      t=t0;
      eval_out=pc.f(t);
      eval_ref=c.f((t-t0)/(t1-t0));
      TEST_ASSERT(tol.approximately_equal(eval_out, eval_ref));
      t=t1;
      eval_out=pc.f(t);
      eval_ref=c.f((t-t0)/(t1-t0));
      TEST_ASSERT(tol.approximately_equal(eval_out, eval_ref));

      // test evaluation at interior point
      t=(t1-t0)*static_cast<data_type>(0.2)+t0;
      eval_out=pc.f(t);
      eval_ref=c.f((t-t0)/(t1-t0));
      TEST_ASSERT(tol.approximately_equal(eval_out, eval_ref));

//      if (typeid(data_type)==typeid(float))
//      {
//        std::cout.flush();
//        eli::test::octave_start(1);
//        eli::test::octave_print(1, c, "poly");
//        eli::test::octave_print(1, pc, "piecewise");
//        eli::test::octave_finish(1, false);
//      }
    }

    void full_cycle_test()
    {
      eli::geom::curve::piecewise_polynomial_creator<data_type, 3, tolerance_type> ppc;
      piecewise_curve_type pc;
      polynomial_type c;
      coefficient_type coef1(5), coef2(2);
      point_type eval_out, eval_ref;
      data_type t0, t1;
      bool rtn_flag;

      // set coefficients
      coef1 << 2, 4, 3, 1, 2;
      coef2 << 1, 1;
      c.set_coefficients(coef1, 0);
      c.set_coefficients(coef2, 1);

      // set the parameterization
      t0=-1;
      t1=1;

      // create curve
      rtn_flag=ppc.set_conditions(c);
      TEST_ASSERT(rtn_flag);
      ppc.set_t0(t0);
      ppc.set_segment_dt(t1-t0, 0);
      rtn_flag=ppc.create(pc);
      TEST_ASSERT(rtn_flag);

      // extract the curve and obtain the monomial coefficients
      curve_type crv;
      typename curve_type::monomial_coefficient_type coef_out;

      pc.get(crv, 0);
      crv.get_monomial_coefficients(coef_out);

      // test the coefficients
      index_type i;

      i=0;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 0), coef1(i)));
      i=1;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 0), coef1(i)));
      i=2;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 0), coef1(i)));
      i=3;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 0), coef1(i)));
      i=4;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 0), coef1(i)));

      i=0;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 1), coef2(i)));
      i=1;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 1), coef2(i)));
      i=2;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 1), 0));
      i=3;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 1), 0));
      i=4;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 1), 0));

      i=0;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 2), 0));
      i=1;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 2), 0));
      i=2;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 2), 0));
      i=3;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 2), 0));
      i=4;
      TEST_ASSERT(tol.approximately_equal(coef_out(i, 2), 0));
    }
};

#endif
