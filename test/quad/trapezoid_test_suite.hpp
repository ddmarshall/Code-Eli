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

#ifndef trapezoid_test_suite_hpp
#define trapezoid_test_suite_hpp

#include <cmath>    // std::exp
#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <vector>   // std::vector

#include "eli/code_eli.hpp"

#include "eli/quad/trapezoid.hpp"
#include "eli/constants/math.hpp"

template<typename data__>
class trapezoid_test_suite : public Test::Suite
{
  protected:
    void AddTests(const float &)
    {
      TEST_ADD(trapezoid_test_suite<float>::uniform_points_test);
      TEST_ADD(trapezoid_test_suite<float>::nonuniform_points_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(trapezoid_test_suite<double>::uniform_points_test);
      TEST_ADD(trapezoid_test_suite<double>::nonuniform_points_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(trapezoid_test_suite<long double>::uniform_points_test);
      TEST_ADD(trapezoid_test_suite<long double>::nonuniform_points_test);
    }

#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      TEST_ADD(trapezoid_test_suite<dd_real>::uniform_points_test);
      TEST_ADD(trapezoid_test_suite<dd_real>::nonuniform_points_test);
    }

    void AddTests(const qd_real &)
    {
      TEST_ADD(trapezoid_test_suite<qd_real>::uniform_points_test);
      TEST_ADD(trapezoid_test_suite<qd_real>::nonuniform_points_test);
    }
#endif

  public:
    trapezoid_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~trapezoid_test_suite()
    {
    }

  private:
    void uniform_points_test()
    {
      data__ x1(1.0), x2(3.0), F_exact(std::exp(x2)-std::exp(x1)), F_quad;
      size_t n, npts(1001);
      std::vector<data__> x(npts), f(npts);

      // set the x points
      for (n=0; n<npts; ++n)
        x[n]=(n)*(x2-x1)/(npts-1)+x1;

      // set the function points
      for (n=0; n<npts; ++n)
        f[n]=std::exp(x[n]);

      // find the integral approximation
      eli::quad::trapezoid<data__> quad;

      F_quad=quad(x[1]-x[0], f.begin(), f.end());

      if (typeid(data__)==typeid(float))
      {
        TEST_ASSERT_DELTA(1, F_quad/F_exact, 2e-5);
      }
      else
      {
        TEST_ASSERT_DELTA(1, F_quad/F_exact, 4e-7);
      }
    }

    void nonuniform_points_test()
    {
      data__ x1(1.0), x2(3.0), F_exact(std::exp(x2)-std::exp(x1)), F_quad;
      size_t n, npts(1001);
      std::vector<data__> x(npts), f(npts);

      // set the x points
      for (n=0; n<npts; ++n)
        x[n]=sin(eli::constants::math<data__>::pi_by_two()*n/(npts-1))*(x2-x1)+x1;

      // set the function points
      for (n=0; n<npts; ++n)
        f[n]=std::exp(x[n]);

      // find the integral approximation
      eli::quad::trapezoid<data__> quad;

      F_quad=quad(x.begin(), f.begin(), f.end());

      TEST_ASSERT_DELTA(1, F_quad/F_exact, 5e-7);
    }
};

#endif
