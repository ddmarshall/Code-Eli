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

#ifndef tolerance_test_suite_hpp
#define tolerance_test_suite_hpp

#include "Eigen/Eigen"

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include <limits>     // numeric_limits

template<typename data__>
class tolerance_test_suite : public Test::Suite
{
  protected:
    void AddTests(const float &)
    {
      TEST_ADD(tolerance_test_suite<float>::exactly_equal_test);
      TEST_ADD(tolerance_test_suite<float>::approximately_equal_test);
      TEST_ADD(tolerance_test_suite<float>::approximately_less_than_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(tolerance_test_suite<double>::exactly_equal_test);
      TEST_ADD(tolerance_test_suite<double>::approximately_equal_test);
      TEST_ADD(tolerance_test_suite<double>::approximately_less_than_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(tolerance_test_suite<long double>::exactly_equal_test);
      TEST_ADD(tolerance_test_suite<long double>::approximately_equal_test);
      TEST_ADD(tolerance_test_suite<long double>::approximately_less_than_test);
    }
#ifdef ELI_USING_QD
    void AddTests(const dd_real &)
    {
      TEST_ADD(tolerance_test_suite<dd_real>::exactly_equal_test);
      TEST_ADD(tolerance_test_suite<dd_real>::approximately_equal_test);
      TEST_ADD(tolerance_test_suite<dd_real>::approximately_less_than_test);
    }

    void AddTests(const qd_real &)
    {
      TEST_ADD(tolerance_test_suite<qd_real>::exactly_equal_test);
      TEST_ADD(tolerance_test_suite<qd_real>::approximately_equal_test);
      TEST_ADD(tolerance_test_suite<qd_real>::approximately_less_than_test);
    }
#endif

  public:
    tolerance_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~tolerance_test_suite()
    {
    }

  private:
    typedef data__ data_type;
    typedef std::numeric_limits<data_type> limits_type;

  private:
    void exactly_equal_test()
    {
      eli::util::tolerance<data_type> tol(1e-6, 1e-3);
      data__ rel_tol(tol.get_relative_tolerance());

      // test same type comparisons
      {
        data_type d1, d2;

        d1=1;
        d2=1;
        TEST_ASSERT(tol.exactly_equal(d1, d2));
        TEST_ASSERT(tol.exactly_equal(d2, d1));
        d1=1;
        d2=1+rel_tol;
        TEST_ASSERT(!tol.exactly_equal(d1, d2));
        TEST_ASSERT(!tol.exactly_equal(d2, d1));
        d1=2;
        d2=1;
        TEST_ASSERT(!tol.exactly_equal(d1, d2));
        TEST_ASSERT(!tol.exactly_equal(d2, d1));
      }

      // test mixed type comparisons
      {
        data_type d2;
        int i1;

        i1=1;
        d2=1;
        TEST_ASSERT(tol.exactly_equal(i1, d2));
        TEST_ASSERT(tol.exactly_equal(d2, i1));
        i1=1;
        d2=1+rel_tol;
        TEST_ASSERT(!tol.exactly_equal(i1, d2));
        TEST_ASSERT(!tol.exactly_equal(d2, i1));
        i1=2;
        d2=1;
        TEST_ASSERT(!tol.exactly_equal(i1, d2));
        TEST_ASSERT(!tol.exactly_equal(d2, i1));
      }
    }

    void approximately_equal_test()
    {
      eli::util::tolerance<data_type> tol(1e-6, 1e-3);
      data__ rel_tol(tol.get_relative_tolerance()), abs_tol(tol.get_absolute_tolerance());

      // test same type comparisons
      {
        data_type d1, d2;

        // less than absolute error
        d1=1;
        d2=1+static_cast<data_type>(0.09)*abs_tol;
        TEST_ASSERT(tol.approximately_equal(d1, d2));
        TEST_ASSERT(tol.approximately_equal(d2, d1));
        d1=1;
        d2=1+static_cast<data_type>(0.99)*abs_tol;
        TEST_ASSERT(tol.approximately_equal(d1, d2));
        TEST_ASSERT(tol.approximately_equal(d2, d1));
        d1=0;
        d2=static_cast<data_type>(0.99)*abs_tol;
        TEST_ASSERT(tol.approximately_equal(d1, d2));
        TEST_ASSERT(tol.approximately_equal(d2, d1));
        d1=0;
        d2=static_cast<data_type>(1.01)*abs_tol;
        TEST_ASSERT(!tol.approximately_equal(d1, d2));
        TEST_ASSERT(!tol.approximately_equal(d2, d1));

        // less than relative error
        d1=1;
        d2=1+static_cast<data_type>(1.1)*abs_tol;
        TEST_ASSERT(tol.approximately_equal(d1, d2));
        TEST_ASSERT(tol.approximately_equal(d2, d1));
        d1=1;
        d2=1+static_cast<data_type>(0.99)*rel_tol;
        TEST_ASSERT(tol.approximately_equal(d1, d2));
        TEST_ASSERT(tol.approximately_equal(d2, d1));
        d1=1;
        d2=1+static_cast<data_type>(1.01)*rel_tol;
        TEST_ASSERT(!tol.approximately_equal(d1, d2));
        TEST_ASSERT(!tol.approximately_equal(d2, d1));
      }

      // test mixed type comparisons
      {
        data_type d2;
        int i1;

        // less than absolute error
        i1=1;
        d2=1+static_cast<data_type>(0.09)*abs_tol;
        TEST_ASSERT(tol.approximately_equal(i1, d2));
        TEST_ASSERT(tol.approximately_equal(d2, i1));
        i1=1;
        d2=1+static_cast<data_type>(0.99)*abs_tol;
        TEST_ASSERT(tol.approximately_equal(i1, d2));
        TEST_ASSERT(tol.approximately_equal(d2, i1));
        i1=0;
        d2=static_cast<data_type>(0.99)*abs_tol;
        TEST_ASSERT(tol.approximately_equal(i1, d2));
        TEST_ASSERT(tol.approximately_equal(d2, i1));
        i1=0;
        d2=static_cast<data_type>(1.01)*abs_tol;
        TEST_ASSERT(!tol.approximately_equal(i1, d2));
        TEST_ASSERT(!tol.approximately_equal(d2, i1));

        // less than relative error
        i1=1;
        d2=1+static_cast<data_type>(1.1)*abs_tol;
        TEST_ASSERT(tol.approximately_equal(i1, d2));
        TEST_ASSERT(tol.approximately_equal(d2, i1));
        i1=1;
        d2=1+static_cast<data_type>(0.99)*rel_tol;
        TEST_ASSERT(tol.approximately_equal(i1, d2));
        TEST_ASSERT(tol.approximately_equal(d2, i1));
        i1=1;
        d2=1+static_cast<data_type>(1.01)*rel_tol;
        TEST_ASSERT(!tol.approximately_equal(i1, d2));
        TEST_ASSERT(!tol.approximately_equal(d2, i1));
      }
    }

    void approximately_less_than_test()
    {
      eli::util::tolerance<data_type> tol(1e-6, 1e-3);
      data__ rel_tol(tol.get_relative_tolerance());

      // test same type comparisons
      {
        data_type d1, d2;

        d1=1;
        d2=1+static_cast<data_type>(0.99)*rel_tol;
        TEST_ASSERT(!tol.approximately_less_than(d1, d2));
        TEST_ASSERT(!tol.approximately_less_than(d2, d1));
        d1=1;
        d2=2;
        TEST_ASSERT(tol.approximately_less_than(d1, d2));
        TEST_ASSERT(!tol.approximately_less_than(d2, d1));
        d1=1;
        d2=1+static_cast<data_type>(1.01)*rel_tol;
        TEST_ASSERT(tol.approximately_less_than(d1, d2));
        TEST_ASSERT(!tol.approximately_less_than(d2, d1));
      }

      // test mixed type comparisons
      {
        data_type d2;
        int i1;

        i1=1;
        d2=1+static_cast<data_type>(0.99)*rel_tol;
        TEST_ASSERT(!tol.approximately_less_than(i1, d2));
        TEST_ASSERT(!tol.approximately_less_than(d2, i1));
        i1=1;
        d2=2;
        TEST_ASSERT(tol.approximately_less_than(i1, d2));
        TEST_ASSERT(!tol.approximately_less_than(d2, i1));
        i1=1;
        d2=1+static_cast<data_type>(1.01)*rel_tol;
        TEST_ASSERT(tol.approximately_less_than(i1, d2));
        TEST_ASSERT(!tol.approximately_less_than(d2, i1));
      }
    }
};

#endif
