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

#ifndef poly_basic_test_suite_hpp
#define poly_basic_test_suite_hpp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/code_eli.hpp"

#include "eli/poly/polynomial.hpp"

template<typename data__>
class poly_basic_test_suite : public Test::Suite
{
  protected:
    void AddTests(const int &)
    {
      TEST_ADD(poly_basic_test_suite<int>::assignment_test);
      TEST_ADD(poly_basic_test_suite<int>::evaluation_test);
    }

    void AddTests(const float &)
    {
      TEST_ADD(poly_basic_test_suite<float>::assignment_test);
      TEST_ADD(poly_basic_test_suite<float>::evaluation_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(poly_basic_test_suite<double>::assignment_test);
      TEST_ADD(poly_basic_test_suite<double>::evaluation_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(poly_basic_test_suite<long double>::assignment_test);
      TEST_ADD(poly_basic_test_suite<long double>::evaluation_test);
    }

#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      TEST_ADD(poly_basic_test_suite<dd_real>::assignment_test);
      TEST_ADD(poly_basic_test_suite<dd_real>::evaluation_test);
    }

    void AddTests(const qd_real &)
    {
      TEST_ADD(poly_basic_test_suite<qd_real>::assignment_test);
      TEST_ADD(poly_basic_test_suite<qd_real>::evaluation_test);
    }
#endif

  public:
    poly_basic_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~poly_basic_test_suite()
    {
    }

  private:

    void assignment_test()
    {
      eli::poly::polynomial<data__> p1, p2;
      typename eli::poly::polynomial<data__>::coefficient_type coef_in(6), coef_out;

      // test default constructor then set coefficients
      coef_in << 2, 0, 1, 4, 1, 8;
      p1.set_coefficients(coef_in);
      p1.get_coefficients(coef_out);
      TEST_ASSERT(coef_in==coef_out);
      coef_out.setZero();

      // test constructor with vector of coefficients
      eli::poly::polynomial<data__> pc1(coef_in);
      pc1.get_coefficients(coef_out);
      TEST_ASSERT(coef_in==coef_out);
      coef_out.setZero();

      // test copy ctr
      eli::poly::polynomial<data__> p1c(p1);
      p1c.get_coefficients(coef_out);
      TEST_ASSERT(coef_in==coef_out);
      coef_out.setZero();

      // test assignment operator
      p2=p1;
      p2.get_coefficients(coef_out);
      TEST_ASSERT(coef_in==coef_out);
      coef_out.setZero();

      // test order
      TEST_ASSERT(p2.degree()==coef_in.rows()-1);

      // test data assignment operator
      data__ zero(0);
      p2=zero;
      p2.get_coefficients(coef_out);
      TEST_ASSERT(p2.degree()==0);
      TEST_ASSERT(coef_out(0)==zero);
    }

    void evaluation_test()
    {
      typename eli::poly::polynomial<data__>::coefficient_type coef_in(5);

      // set coefficients
      coef_in << 2, 4, 3, 1, 2;

      eli::poly::polynomial<data__> p1(coef_in), p2;
      typename eli::poly::polynomial<data__>::data_type eval_out, eval_ref;
      typename eli::poly::polynomial<data__>::data_type t;

      // test evaluation at points
      t=static_cast<data__>(0);
      eval_out=p1.f(t);
      eval_ref=coef_in(0);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=p1.f(t);
      eval_ref = coef_in.sum();
      TEST_ASSERT(eval_out==eval_ref);
      t=10;
      eval_out=p1.f(t);
      eval_ref = 21342;

      TEST_ASSERT(eval_out==eval_ref);
    }
};

#endif
