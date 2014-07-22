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

#ifndef simpson_test_suite_hpp
#define simpson_test_suite_hpp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/code_eli.hpp"

#include "eli/mutil/quad/simpson.hpp"
#include "eli/constants/math.hpp"


template <typename data__>
struct exp_functor
{
  data__ operator()(const data__ &x) {return std::exp(x);}
};

// default is for bad values
template <typename data__>
void get_reference_adaptive_params(typename eli::mutil::quad::simpson<data__>::adaptive_params &ap)
{
  ap.recursion_depth=0;
  ap.function_count=0;
  ap.coarse_value=-1;
  ap.fine_value=-1;
  ap.approximate_error=-1;
}

template <>
void get_reference_adaptive_params<float>(typename eli::mutil::quad::simpson<float>::adaptive_params &ap)
{
  ap.recursion_depth=2;
  ap.function_count=9;
  ap.coarse_value=17.3731117248535f;
  ap.fine_value=17.3676300048828f;
#ifdef _MSC_VER
# ifdef _WIN64
  ap.approximate_error=5.48267344129272e-05f;
# else
  ap.approximate_error=5.4817199e-05f;
# endif
#else
  ap.approximate_error=5.48267344129272e-05f;
#endif
}

template <>
void get_reference_adaptive_params<double>(typename eli::mutil::quad::simpson<double>::adaptive_params &ap)
{
  ap.recursion_depth=5;
  ap.function_count=65;
  ap.coarse_value=17.367256566284723363;
  ap.fine_value=17.367255186732954542;
  ap.approximate_error=1.3795517685433190629e-08;
}

template <>
void get_reference_adaptive_params<long double>(typename eli::mutil::quad::simpson<long double>::adaptive_params &ap)
{
#ifdef _MSC_VER
  ap.recursion_depth=5;
  ap.function_count=65;
  ap.coarse_value=17.367256566284723363L;
  ap.fine_value=17.367255186732954542L;
  ap.approximate_error=1.3795517685433190629e-08L;
#else
  ap.recursion_depth=6;
  ap.function_count=129;
  ap.coarse_value=17.367255186732954673478612761527983820997178554534912109375L;
  ap.fine_value=17.36725510047939464303157208746597461868077516555786133L;
# ifdef __INTEL_COMPILER
  ap.approximate_error=8.6253560028725869727e-10L;
# else
  ap.approximate_error=8.62535600286716596198481699594144637825189614321885756e-10L;
# endif
#endif
}

template<typename data__>
class simpson_test_suite : public Test::Suite
{
  protected:
    void AddTests(const float &)
    {
      TEST_ADD(simpson_test_suite<float>::uniform_points_test);
      TEST_ADD(simpson_test_suite<float>::nonuniform_points_test);
      TEST_ADD(simpson_test_suite<float>::adaptive_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(simpson_test_suite<double>::uniform_points_test);
      TEST_ADD(simpson_test_suite<double>::nonuniform_points_test);
      TEST_ADD(simpson_test_suite<double>::adaptive_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(simpson_test_suite<long double>::uniform_points_test);
      TEST_ADD(simpson_test_suite<long double>::nonuniform_points_test);
      TEST_ADD(simpson_test_suite<long double>::adaptive_test);
    }

  public:
    simpson_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~simpson_test_suite()
    {
    }

  private:
    void uniform_points_test()
    {
      data__ x1(1.0), x2(3.0), F_exact(std::exp(x2)-std::exp(x1)), F_quad;

      {
        size_t n, npts(1001);
        std::vector<data__> x(npts), f(npts);

        // set the x points
        for (n=0; n<npts; ++n)
          x[n]=(n)*(x2-x1)/(npts-1)+x1;

        // set the function points
        for (n=0; n<npts; ++n)
          f[n]=std::exp(x[n]);

        // find the integral approximation
        eli::mutil::quad::simpson<data__> quad;

        F_quad=quad(x[1]-x[0], f.begin(), f.end());

        if (typeid(data__)==typeid(float))
        {
          TEST_ASSERT_DELTA(1, F_quad/F_exact, 2e-5);
        }
        else
        {
          TEST_ASSERT_DELTA(1, F_quad/F_exact, 1e-13);
        }
      }

// bug in release build for clang current version (3.3) and older
#ifdef NDEBUG
# ifdef __clang__
#   if ( (__clang_major__ < 3) || ((__clang_major__==3) && (__clang_minor__<=3)) )
      TEST_ASSERT_MSG(false, "Clang release build bug cannot build this test.");
      return;
#   endif
# endif
#endif
      {
        size_t n, npts(1000);
        std::vector<data__> x(npts), f(npts);

        // set the x points
        for (n=0; n<npts; ++n)
          x[n]=(n)*(x2-x1)/(npts-1)+x1;

        // set the function points
        for (n=0; n<npts; ++n)
          f[n]=std::exp(x[n]);

        // find the integral approximation
        eli::mutil::quad::simpson<data__> quad;

        F_quad=quad(x[1]-x[0], f.begin(), f.end());

        if (typeid(data__)==typeid(float))
        {
          TEST_ASSERT_DELTA(1, F_quad/F_exact, 1e-5);
        }
        else
        {
          TEST_ASSERT_DELTA(1, F_quad/F_exact, 1e-12);
        }
      }
    }

    void nonuniform_points_test()
    {
      data__ x1(1.0), x2(3.0), F_exact(std::exp(x2)-std::exp(x1)), F_quad;

      {
        size_t n, npts(1001);
        std::vector<data__> x(npts), f(npts);

        // set the x points
        for (n=0; n<npts; ++n)
          x[n]=sin(eli::constants::math<data__>::pi_by_two()*n/(npts-1))*(x2-x1)+x1;

        // set the function points
        for (n=0; n<npts; ++n)
          f[n]=std::exp(x[n]);

        // find the integral approximation
        eli::mutil::quad::simpson<data__> quad;

        F_quad=quad(x.begin(), f.begin(), f.end());

        if (typeid(data__)==typeid(float))
        {
          TEST_ASSERT_DELTA(1, F_quad/F_exact, 3e-7);
        }
        else
        {
          TEST_ASSERT_DELTA(1, F_quad/F_exact, 5e-13);
        }
      }

      {
        size_t n, npts(1000);
        std::vector<data__> x(npts), f(npts);

        // set the x points
        for (n=0; n<npts; ++n)
          x[n]=sin(eli::constants::math<data__>::pi_by_two()*n/(npts-1))*(x2-x1)+x1;

        // set the function points
        for (n=0; n<npts; ++n)
          f[n]=std::exp(x[n]);

        // find the integral approximation
        eli::mutil::quad::simpson<data__> quad;

        F_quad=quad(x.begin(), f.begin(), f.end());

        if (typeid(data__)==typeid(float))
        {
          TEST_ASSERT_DELTA(1, F_quad/F_exact, 9e-7);
        }
        else
        {
          TEST_ASSERT_DELTA(1, F_quad/F_exact, 5e-13);
        }
      }
    }

    void adaptive_test()
    {
      data__ x0(1.0), x1(3.0), F_exact(std::exp(x1)-std::exp(x0)), F_quad, tol(std::numeric_limits<data__>::epsilon());

      // find the integral approximation
      eli::mutil::quad::simpson<data__> quad;
      typename eli::mutil::quad::simpson<data__>::adaptive_params ap, ap2, ap_ref;

      // integrate with depth as termination case
      ap.tolerance=tol;
      ap.max_depth=2;
      F_quad=quad(exp_functor<data__>(), x0, x1, ap);
      TEST_ASSERT(ap.recursion_depth==ap.max_depth);
      TEST_ASSERT(std::abs(1 - F_quad/F_exact)>ap.tolerance);

      // integrate with tolerance as termination case
      ap.tolerance=tol;
      ap.max_depth=23;
      F_quad=quad(exp_functor<data__>(), x0, x1, ap);
      TEST_ASSERT(ap.recursion_depth<ap.max_depth);
      TEST_ASSERT_DELTA(1, F_quad/F_exact, ap.tolerance);

      // integrate to a known state
      ap=ap2;
      F_quad=quad(exp_functor<data__>(), x0, x1, ap);
      get_reference_adaptive_params<data__>(ap_ref);

      TEST_ASSERT(ap.max_depth==30);
      TEST_ASSERT(ap.tol_factor==1.25);
      TEST_ASSERT(ap.error_factor==100.0);
      TEST_ASSERT(ap.tolerance==ap_ref.tolerance);
      TEST_ASSERT_DELTA(1, F_quad/F_exact, ap.tolerance);
      TEST_ASSERT(ap.recursion_depth==ap_ref.recursion_depth);
      TEST_ASSERT(ap.function_count==ap_ref.function_count);
      TEST_ASSERT_DELTA(1, ap.coarse_value/ap_ref.coarse_value, std::numeric_limits<data__>::epsilon());
      TEST_ASSERT_DELTA(1, ap.fine_value/ap_ref.fine_value, std::numeric_limits<data__>::epsilon());
      TEST_ASSERT_DELTA(1, ap.approximate_error/ap_ref.approximate_error, std::numeric_limits<data__>::epsilon());
    }
};

#endif
