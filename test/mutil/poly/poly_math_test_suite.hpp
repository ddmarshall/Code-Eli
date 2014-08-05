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

#ifndef poly_math_test_suite_hpp
#define poly_math_test_suite_hpp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/mutil/poly/polynomial.hpp"

template<typename data__>
class poly_math_test_suite : public Test::Suite
{
  protected:
    void AddTests(const int &)
    {
      TEST_ADD(poly_math_test_suite<int>::add_test);
      TEST_ADD(poly_math_test_suite<int>::subtract_test);
      TEST_ADD(poly_math_test_suite<int>::multiply_test);
      TEST_ADD(poly_math_test_suite<int>::divide_test);
      TEST_ADD(poly_math_test_suite<int>::creation_test);
      TEST_ADD(poly_math_test_suite<int>::derivative_test);
    }

    void AddTests(const float &)
    {
      TEST_ADD(poly_math_test_suite<float>::add_test);
      TEST_ADD(poly_math_test_suite<float>::subtract_test);
      TEST_ADD(poly_math_test_suite<float>::multiply_test);
      TEST_ADD(poly_math_test_suite<float>::divide_test);
      TEST_ADD(poly_math_test_suite<float>::creation_test);
      TEST_ADD(poly_math_test_suite<float>::derivative_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(poly_math_test_suite<double>::add_test);
      TEST_ADD(poly_math_test_suite<double>::subtract_test);
      TEST_ADD(poly_math_test_suite<double>::multiply_test);
      TEST_ADD(poly_math_test_suite<double>::divide_test);
      TEST_ADD(poly_math_test_suite<double>::creation_test);
      TEST_ADD(poly_math_test_suite<double>::derivative_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(poly_math_test_suite<long double>::add_test);
      TEST_ADD(poly_math_test_suite<long double>::subtract_test);
      TEST_ADD(poly_math_test_suite<long double>::multiply_test);
      TEST_ADD(poly_math_test_suite<long double>::divide_test);
      TEST_ADD(poly_math_test_suite<long double>::creation_test);
      TEST_ADD(poly_math_test_suite<long double>::derivative_test);
    }

  public:
    poly_math_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~poly_math_test_suite()
    {
    }

  private:

    void add_test()
    {
      typename eli::mutil::poly::polynomial<data__>::coefficient_type coef1(4), coef2(3);

      // set coefficients
      coef1 << 0, 0, 8, 4;
      coef2 << 1, 2, 3;

      // create two polynomials
      eli::mutil::poly::polynomial<data__> p1(coef1), p2(coef2), pr;
      typename eli::mutil::poly::polynomial<data__>::data_type eval_out, eval_ref;
      typename eli::mutil::poly::polynomial<data__>::data_type x[5] = {0, 1, 6, 9, 15};

      // create polynomial that is addition of two other polynomials
      pr.add(p1, p2);

      // test resulting polynomial at several points
      for (size_t i=0; i<5; ++i)
      {
        // test point
        eval_out=pr.f(x[i]);
        eval_ref=p1.f(x[i]) + p2.f(x[i]);
        TEST_ASSERT(eval_out==eval_ref);
      }
    }

    void subtract_test()
    {
      typename eli::mutil::poly::polynomial<data__>::coefficient_type coef1(4), coef2(3);

      // set coefficients
      coef1 << 0, 0, 8, 4;
      coef2 << 1, 2, 3;

      // create polynomial that is addition of two other polynomials
      eli::mutil::poly::polynomial<data__> p1(coef1), p2(coef2), pr;
      typename eli::mutil::poly::polynomial<data__>::data_type eval_out, eval_ref;
      typename eli::mutil::poly::polynomial<data__>::data_type x[5] = {0, 1, 6, 9, 15};

      // create polynomial that is subtraction of two other polynomials
      pr.subtract(p1, p2);

      // test resulting polynomial at several points
      for (size_t i=0; i<5; ++i)
      {
        // test point
        eval_out=pr.f(x[i]);
        eval_ref=p1.f(x[i]) - p2.f(x[i]);
        TEST_ASSERT(eval_out==eval_ref);
      }
    }

    void multiply_test()
    {
      typename eli::mutil::poly::polynomial<data__>::coefficient_type coef1(4), coef2(3);

      // set coefficients
      coef1 << 0, 0, 8, 4;
      coef2 << 1, 2, 3;

      // create polynomial that is addition of two other polynomials
      eli::mutil::poly::polynomial<data__> p1(coef1), p2(coef2), pr;
      typename eli::mutil::poly::polynomial<data__>::data_type eval_out, eval_ref;
      typename eli::mutil::poly::polynomial<data__>::data_type x[5] = {0, 1, 6, 9, 15};

      // create polynomial that is multiplication of two other polynomials
      pr.multiply(p1, p2);

      // test resulting polynomial at several points
      for (size_t i=0; i<5; ++i)
      {
        // test point
        eval_out=pr.f(x[i]);
        eval_ref=p1.f(x[i]) * p2.f(x[i]);
        TEST_ASSERT(eval_out==eval_ref);
      }
    }

    void divide_test()
    {
      // test against known results: simple test
      {
        typename eli::mutil::poly::polynomial<data__>::coefficient_type coef1(4), coef2(3), coef_q(2), coef_r(2), coef_out;

        // set coefficients
        coef1  << -3,  4, -2, 3;
        coef2  <<  3,  3,  1;
        coef_q <<-11,  3;
        coef_r << 30, 28;

        // create polynomials
        eli::mutil::poly::polynomial<data__> p1(coef1), p2(coef2), pr, prem;

        // create polynomial that is division of two other polynomials
        pr.divide(prem, p1, p2);

        // test answers
        pr.get_coefficients(coef_out);
        if (coef_out.rows()==coef_q.rows())
        {
          TEST_ASSERT(coef_out==coef_q);
        }
        else
        {
          TEST_ASSERT_MSG(false, "Quotient not correct degree");
        }
        prem.get_coefficients(coef_out);
        if (coef_out.rows()==coef_r.rows())
        {
          TEST_ASSERT(coef_out==coef_r);
        }
        else
        {
          TEST_ASSERT_MSG(false, "Remainder not correct degree");
        }
      }

      // test against known results: this tests zeros in the numerator
      {
        typename eli::mutil::poly::polynomial<data__>::coefficient_type coef1(4), coef2(2), coef_q(3), coef_r(1), coef_out;

        // set coefficients
        coef1  <<-1, 0, 0, 1;
        coef2  << 2, 1;
        coef_q << 4,-2, 1;
        coef_r <<-9;

        // create polynomials
        eli::mutil::poly::polynomial<data__> p1(coef1), p2(coef2), pr, prem;

        // create polynomial that is division of two other polynomials
        pr.divide(prem, p1, p2);

        // test answers
        pr.get_coefficients(coef_out);
        if (coef_out.rows()==coef_q.rows())
        {
          TEST_ASSERT(coef_out==coef_q);
        }
        else
        {
          TEST_ASSERT_MSG(false, "Quotient not correct degree");
        }
        prem.get_coefficients(coef_out);
        if (coef_out.rows()==coef_r.rows())
        {
          TEST_ASSERT(coef_out==coef_r);
        }
        else
        {
          TEST_ASSERT_MSG(false, "Remainder not correct degree");
        }
      }

      // test against know results: this tests multiple steps in the division
      {
        typename eli::mutil::poly::polynomial<data__>::coefficient_type coef1(6), coef2(2), coef_q(5), coef_r(1), coef_out;

        // set coefficients
        coef1  <<-1, 0, 0, 0, 0,-5;
        coef2  << 5, 5;
        coef_q <<-1, 1,-1, 1,-1;
        coef_r << 4;

        // create polynomials
        eli::mutil::poly::polynomial<data__> p1(coef1), p2(coef2), pr, prem;

        // create polynomial that is division of two other polynomials
        pr.divide(prem, p1, p2);

        // test answers
        pr.get_coefficients(coef_out);
        if (coef_out.rows()==coef_q.rows())
        {
          TEST_ASSERT(coef_out==coef_q);
        }
        else
        {
          TEST_ASSERT_MSG(false, "Quotient not correct degree");
        }
        prem.get_coefficients(coef_out);
        if (coef_out.rows()==coef_r.rows())
        {
          TEST_ASSERT(coef_out==coef_r);
        }
        else
        {
          TEST_ASSERT_MSG(false, "Remainder not correct degree");
        }
      }

      // test against know results: this tests with no remainder
      {
        typename eli::mutil::poly::polynomial<data__>::coefficient_type coef1(4), coef2(3), coef_q(2), coef_r(2), coef_out;

        // set coefficients
        coef1  <<-15, 3,-5, 1;
        coef2  <<  3, 0, 1;
        coef_q << -5, 1;
        coef_r <<  0, 0;

        // create polynomials
        eli::mutil::poly::polynomial<data__> p1(coef1), p2(coef2), pr, prem;

        // create polynomial that is division of two other polynomials
        pr.divide(prem, p1, p2);

        // test answers
        pr.get_coefficients(coef_out);
        if (coef_out.rows()==coef_q.rows())
        {
          TEST_ASSERT(coef_out==coef_q);
        }
        else
        {
          TEST_ASSERT_MSG(false, "Quotient not correct degree");
        }
        prem.get_coefficients(coef_out);
        if (coef_out.rows()==coef_r.rows())
        {
          TEST_ASSERT(coef_out==coef_r);
        }
        else
        {
          TEST_ASSERT_MSG(false, "Remainder not correct degree");
        }
      }

      // test against know results: this tests leading zeros in divisor
      {
        typename eli::mutil::poly::polynomial<data__>::coefficient_type coef1(6), coef2(3), coef_q(5), coef_r(2), coef_out;

        // set coefficients
        coef1  <<-1, 0, 0, 0, 0,-5;
        coef2  << 5, 5, 0;
        coef_q <<-1, 1,-1, 1,-1;
        coef_r << 4, 0;

        // create polynomials
        eli::mutil::poly::polynomial<data__> p1(coef1), p2(coef2), pr, prem;

        // create polynomial that is division of two other polynomials
        pr.divide(prem, p1, p2);

        // test answers
        pr.get_coefficients(coef_out);
        if (coef_out.rows()==coef_q.rows())
        {
          TEST_ASSERT(coef_out==coef_q);
        }
        else
        {
          TEST_ASSERT_MSG(false, "Quotient not correct degree");
        }
        prem.get_coefficients(coef_out);
        if (coef_out.rows()==coef_r.rows())
        {
          TEST_ASSERT(coef_out==coef_r);
        }
        else
        {
          TEST_ASSERT_MSG(false, "Remainder not correct degree");
        }
      }

      {
        typename eli::mutil::poly::polynomial<data__>::coefficient_type coef1(4), coef2(3);

        // set control points
        coef1 << 1, 1, 1, 1;
        coef2 << 1, 1, 0;

        // create two polynomials
        eli::mutil::poly::polynomial<data__> p1(coef1), p2(coef2), pr, prem;
        typename eli::mutil::poly::polynomial<data__>::data_type eval_out, eval_ref;
        typename eli::mutil::poly::polynomial<data__>::data_type x[5] = {0, 1, 6, 9, 15};

        // create polynomial that is division of two other polynomials
        pr.divide(prem, p1, p2);

        // test resulting polynomials at several points
        for (size_t i=0; i<5; ++i)
        {
          // test point
          eval_out=pr.f(x[i])+prem.f(x[i])/p2.f(x[i]);
          eval_ref=p1.f(x[i])/p2.f(x[i]);
          TEST_ASSERT(eval_out==eval_ref);
        }
      }
    }

    void creation_test()
    {
      // create polynomial with 1 root
      {
        data__ roots(2);
        eli::mutil::poly::polynomial<data__> p(roots);

        TEST_ASSERT(p.f(roots)==0);
      }

      // create polynomial with 2 roots
      {
        data__ roots[2] = {2, -3};
        eli::mutil::poly::polynomial<data__> p(roots[0], roots[1]);

        TEST_ASSERT(p.f(roots[0])==0);
        TEST_ASSERT(p.f(roots[1])==0);
      }

      // create polynomial with 3 roots
      {
        data__ roots[3] = {2, -3, 5};
        eli::mutil::poly::polynomial<data__> p(roots[0], roots[1], roots[2]);

        TEST_ASSERT(p.f(roots[0])==0);
        TEST_ASSERT(p.f(roots[1])==0);
        TEST_ASSERT(p.f(roots[2])==0);
      }

      // create polynomial with 4 roots
      {
        data__ roots[4] = {2, -3, 5, 7};
        eli::mutil::poly::polynomial<data__> p(roots[0], roots[1], roots[2], roots[3]);

        typename eli::mutil::poly::polynomial<data__>::coefficient_type c;
        p.get_coefficients(c);

        TEST_ASSERT(p.f(roots[0])==0);
        TEST_ASSERT(p.f(roots[1])==0);
        TEST_ASSERT(p.f(roots[2])==0);
        TEST_ASSERT(p.f(roots[3])==0);
      }

      // create polynomial with 4 roots
      {
        data__ roots[4] = {2, -3, 5, 7};
        eli::mutil::poly::polynomial<data__> p(roots, roots+4);

        TEST_ASSERT(p.f(roots[0])==0);
        TEST_ASSERT(p.f(roots[1])==0);
        TEST_ASSERT(p.f(roots[2])==0);
        TEST_ASSERT(p.f(roots[3])==0);
      }

      // create polynomial with many roots
      {
        int i, nr;

        // NOTE: different types fail the below tests at different number of roots
        if (typeid(data__)==typeid(float))
          nr=5;
        else if (typeid(data__)==typeid(double))
          nr=7;
        else if (typeid(data__)==typeid(long double))
          nr=9;
        else if (typeid(data__)==typeid(int)) // for computational time reasons this is limited
          nr=9;                               //  but it worked for 30 and took 12 minutes :(
        else
          nr=15;

        std::vector<data__> roots(nr);
        roots[0]=2;
        for (i=1; i<static_cast<int>(roots.size()); ++i)
          roots[i]=3*i-roots[i-1];

        eli::mutil::poly::polynomial<data__> p(roots.begin(), roots.end());
        for (i=0; i<static_cast<int>(roots.size()); ++i)
        {
          TEST_ASSERT(p.f(roots[i])==0);
        }
      }
    }

    void derivative_test()
    {
      typename eli::mutil::poly::polynomial<data__>::coefficient_type coef_in(4);

      // set coefficients
      coef_in << 3, 5, 1, 8;

      eli::mutil::poly::polynomial<data__> p1(coef_in);
      typename eli::mutil::poly::polynomial<data__>::data_type eval_out, eval_ref;
      typename eli::mutil::poly::polynomial<data__>::data_type t;

      {
        // test 1st derivative function
        eli::mutil::poly::polynomial<data__> *p2=p1.f();

        if (p2==0)
        {
          TEST_ASSERT_MSG(false, "Null derivative function");
        }
        else
        {
          typename eli::mutil::poly::polynomial<data__>::coefficient_type coef_ref, coef_out;

          coef_ref = coef_in;

          // compare coefficients
          p2->get_coefficients(coef_out);
          if (coef_out.rows()==coef_ref.rows())
          {
            TEST_ASSERT(coef_out==coef_ref);
          }
          else
          {
            TEST_ASSERT_MSG(false, "Degree not correct degree");
          }
        }
        delete p2;
      }

      {
        // test 1st derivative at points
        t=static_cast<data__>(0);
        eval_out=p1.fp(t);
        eval_ref=1*coef_in(1);
        TEST_ASSERT(eval_out==eval_ref);
        t=static_cast<data__>(1);
        eval_out=p1.fp(t);
        eval_ref=static_cast<data__>(31);
        TEST_ASSERT(eval_out==eval_ref);
        t=static_cast<data__>(11);
        eval_out=p1.fp(t);
        eval_ref=static_cast<data__>(2931);
        TEST_ASSERT(eval_out==eval_ref);

        // test 1st derivative function
        eli::mutil::poly::polynomial<data__> *p2=p1.fp();

        if (p2==0)
        {
          TEST_ASSERT_MSG(false, "Null derivative function");
        }
        else
        {
          typename eli::mutil::poly::polynomial<data__>::coefficient_type coef_ref(3), coef_out;

          coef_ref << 1*5, 2*1, 3*8;

          // compare coefficients
          p2->get_coefficients(coef_out);
          if (coef_out.rows()==coef_ref.rows())
          {
            TEST_ASSERT(coef_out==coef_ref);
          }
          else
          {
            TEST_ASSERT_MSG(false, "Degree not correct degree");
          }
        }
        delete p2;
      }

      {
        // test 2nd derivative at points
        t=static_cast<data__>(0);
        eval_out=p1.fpp(t);
        eval_ref=2*1*coef_in(2);
        TEST_ASSERT(eval_out==eval_ref);
        t=static_cast<data__>(1);
        eval_out=p1.fpp(t);
        eval_ref=static_cast<data__>(50);
        TEST_ASSERT(eval_out==eval_ref);
        t=static_cast<data__>(11);
        eval_out=p1.fpp(t);
        eval_ref=static_cast<data__>(530);
        TEST_ASSERT(eval_out==eval_ref);

        // test 2nd derivative function
        eli::mutil::poly::polynomial<data__> *p2=p1.fpp();

        if (p2==0)
        {
          TEST_ASSERT_MSG(false, "Null derivative function");
        }
        else
        {
          typename eli::mutil::poly::polynomial<data__>::coefficient_type coef_ref(2), coef_out;

          coef_ref << (2*1)*1, (3*2)*8;

          // compare coefficients
          p2->get_coefficients(coef_out);
          if (coef_out.rows()==coef_ref.rows())
          {
            TEST_ASSERT(coef_out==coef_ref);
          }
          else
          {
            TEST_ASSERT_MSG(false, "Degree not correct degree");
          }
        }
        delete p2;
      }

      {
        // test 3rd derivative at points
        t=static_cast<data__>(0);
        eval_out=p1.fppp(t);
        eval_ref=3*2*1*coef_in(3);
        TEST_ASSERT(eval_out==eval_ref);
        t=static_cast<data__>(1);
        eval_out=p1.fppp(t);
        eval_ref=static_cast<data__>(48);
        TEST_ASSERT(eval_out==eval_ref);
        t=static_cast<data__>(11);
        eval_out=p1.fppp(t);
        eval_ref=static_cast<data__>(48);
        TEST_ASSERT(eval_out==eval_ref);

        // test 3rd derivative function
        eli::mutil::poly::polynomial<data__> *p2=p1.fppp();

        if (p2==0)
        {
          TEST_ASSERT_MSG(false, "Null derivative function");
        }
        else
        {
          typename eli::mutil::poly::polynomial<data__>::coefficient_type coef_ref(1), coef_out;

          coef_ref << (3*2*1)*8;

          // compare coefficients
          p2->get_coefficients(coef_out);
          if (coef_out.rows()==coef_ref.rows())
          {
            TEST_ASSERT(coef_out==coef_ref);
          }
          else
          {
            TEST_ASSERT_MSG(false, "Degree not correct degree");
          }
        }
        delete p2;
      }
    }

};

#endif
