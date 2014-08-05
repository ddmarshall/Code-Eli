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

#ifndef mutil_fd_testsuite_hpp
#define mutil_fd_testsuite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <vector>   // std::vector

#include "eli/mutil/fd.hpp"

// NOTE: Derive and implememt d1o4, d1o5, d1o6, d2o3, d2o4, d2o5, d2o6
// NOTE: How to extend this to multiple dimensional derivatives?
// NOTE: Can implement the MUSCL, FV discretization, (structured?) etc. as same type

namespace sveli
{
  namespace post
  {
    template<typename itn__, typename ite__, typename data__>
    int calculate_order(itn__ itn, ite__ iteb, ite__ itee, const data__ &refrat);
//
//    template<typename itn__, typename ite__, typename itr__>
//    int calculate_order(itn__ itn, ite__ iteb, ite__ itee, itr__ itr);
  }
}

#include <cmath> // std::log

template<typename itn__, typename ite__, typename data__>
int calculate_order(itn__ itn, ite__ iteb, ite__ itee, const data__ &refrat)
{
  ite__ ite(iteb), itep(iteb);

  for (++ite; ite!=itee; ++ite, ++itep, ++itn)
    (*itn)=std::log(std::abs((*itep)/(*ite)))/std::log(refrat);

  return 0;
}

template<typename data__>
class fd_test_suite : public Test::Suite
{
  protected:
    void AddTests(const float &)
    {
      TEST_ADD(fd_test_suite<float>::d1o1_test);
      TEST_ADD(fd_test_suite<float>::d1o2_test);
      TEST_ADD(fd_test_suite<float>::d1o3_test);
      TEST_ADD(fd_test_suite<float>::d2o1_test);
      TEST_ADD(fd_test_suite<float>::d2o2_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(fd_test_suite<double>::d1o1_test);
      TEST_ADD(fd_test_suite<double>::d1o2_test);
      TEST_ADD(fd_test_suite<double>::d1o3_test);
      TEST_ADD(fd_test_suite<double>::d2o1_test);
      TEST_ADD(fd_test_suite<double>::d2o2_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(fd_test_suite<long double>::d1o1_test);
      TEST_ADD(fd_test_suite<long double>::d1o2_test);
      TEST_ADD(fd_test_suite<long double>::d1o3_test);
      TEST_ADD(fd_test_suite<long double>::d2o1_test);
      TEST_ADD(fd_test_suite<long double>::d2o2_test);
    }

  public:
    fd_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~fd_test_suite()
    {
    }

  private:
    void d1o1_test()
    {
      typedef eli::mutil::fd::d1o1<data__> d1o1_type;

      std::vector<data__> dx(1); dx[0]=static_cast<data__>(1.0);
      d1_driver("*** 1st Derivative, 1st Order Left Uniform Spacing ***", dx, d1o1_type(d1o1_type::LEFT));
      d1_driver("*** 1st Derivative, 1st Order Right Uniform Spacing ***", dx, d1o1_type(d1o1_type::RIGHT));

      std::vector<data__> dx_factor(2); dx_factor[0]=static_cast<data__>(1.0); dx_factor[1]=static_cast<data__>(1.0);
      d1_driver("*** 1st Derivative, 1st Order Left Uniform Spacing ***", dx_factor, d1o1_type(d1o1_type::LEFT));
      d1_driver("*** 1st Derivative, 1st Order Right Uniform Spacing ***", dx_factor, d1o1_type(d1o1_type::RIGHT));
      dx_factor[0]=static_cast<data__>(0.95); dx_factor[1]=static_cast<data__>(1.00);
      d1_driver("*** 1st Derivative, 1st Order Left Mildly Non-Uniform Spacing ***", dx_factor, d1o1_type(d1o1_type::LEFT));
      dx_factor[0]=static_cast<data__>(1.00); dx_factor[1]=static_cast<data__>(0.95);
      d1_driver("*** 1st Derivative, 1st Order Right Mildly Non-Uniform Spacing ***", dx_factor, d1o1_type(d1o1_type::RIGHT));
      dx_factor[0]=static_cast<data__>(0.1); dx_factor[1]=static_cast<data__>(1.00);
      d1_driver("*** 1st Derivative, 1st Order Left Extremely Non-Uniform Spacing ***", dx_factor, d1o1_type(d1o1_type::LEFT));
      dx_factor[0]=static_cast<data__>(1.0); dx_factor[1]=static_cast<data__>(0.1);
      d1_driver("*** 1st Derivative, 1st Order Right Extremely Non-Uniform Spacing ***", dx_factor, d1o1_type(d1o1_type::RIGHT));
    }

    void d1o2_test()
    {
      typedef eli::mutil::fd::d1o2<data__> d1o2_type;

      std::vector<data__> dx(1); dx[0]=static_cast<data__>(1.0);
      d1_driver("*** 1st Derivative, 2nd Order Center Uniform Spacing ***", dx, d1o2_type(d1o2_type::CENTER));
      d1_driver("*** 1st Derivative, 2nd Order Left Uniform Spacing ***", dx, d1o2_type(d1o2_type::LEFT));
      d1_driver("*** 1st Derivative, 2nd Order Right Uniform Spacing ***", dx, d1o2_type(d1o2_type::RIGHT));

      std::vector<data__> dx_factor(3); dx_factor[0]=static_cast<data__>(1.0); dx_factor[1]=static_cast<data__>(1.0); dx_factor[2]=static_cast<data__>(1.0);
      d1_driver("*** 1st Derivative, 2nd Order Center Uniform Spacing ***", dx_factor, d1o2_type(d1o2_type::CENTER));
      d1_driver("*** 1st Derivative, 2nd Order Left Uniform Spacing ***", dx_factor, d1o2_type(d1o2_type::LEFT));
      d1_driver("*** 1st Derivative, 2nd Order Right Uniform Spacing ***", dx_factor, d1o2_type(d1o2_type::RIGHT));
      dx_factor[0]=static_cast<data__>(0.95); dx_factor[1]=static_cast<data__>(1.0); dx_factor[2]=static_cast<data__>(1.05);
      d1_driver("*** 1st Derivative, 2nd Order Center Mildly Non-Uniform Spacing ***", dx_factor, d1o2_type(d1o2_type::CENTER));
      dx_factor[0]=static_cast<data__>(1.025); dx_factor[1]=static_cast<data__>(0.95); dx_factor[2]=static_cast<data__>(1.00);
      d1_driver("*** 1st Derivative, 2nd Order Left Mildly Non-Uniform Spacing ***", dx_factor, d1o2_type(d1o2_type::LEFT));
      dx_factor[0]=static_cast<data__>(1.00); dx_factor[1]=static_cast<data__>(0.95); dx_factor[2]=static_cast<data__>(1.025);
      d1_driver("*** 1st Derivative, 2nd Order Right Mildly Non-Uniform Spacing ***", dx_factor, d1o2_type(d1o2_type::RIGHT));
      dx_factor[0]=static_cast<data__>(0.1); dx_factor[1]=static_cast<data__>(1.0); dx_factor[2]=static_cast<data__>(1.8);
      d1_driver("*** 1st Derivative, 2nd Order Center Extremely Non-Uniform Spacing ***", dx_factor, d1o2_type(d1o2_type::CENTER));
      dx_factor[0]=static_cast<data__>(0.95); dx_factor[1]=static_cast<data__>(0.1); dx_factor[2]=static_cast<data__>(1.00);
      d1_driver("*** 1st Derivative, 2nd Order Left Extremely Non-Uniform Spacing ***", dx_factor, d1o2_type(d1o2_type::LEFT));
      dx_factor[0]=static_cast<data__>(1.0); dx_factor[1]=static_cast<data__>(0.1); dx_factor[2]=static_cast<data__>(0.95);
      d1_driver("*** 1st Derivative, 2nd Order Right Extremely Non-Uniform Spacing ***", dx_factor, d1o2_type(d1o2_type::RIGHT));
    }

    void d1o3_test()
    {
      typedef eli::mutil::fd::d1o3<data__> d1o3_type;

      std::vector<data__> dx(1); dx[0]=static_cast<data__>(1.0);
      d1_driver("*** 1st Derivative, 3rd Order Left Uniform Spacing ***", dx, d1o3_type(d1o3_type::LEFT));
      d1_driver("*** 1st Derivative, 3rd Order Left-Biased Uniform Spacing ***", dx, d1o3_type(d1o3_type::LEFT_BIASED));
      d1_driver("*** 1st Derivative, 3rd Order Right-Biased Uniform Spacing ***", dx, d1o3_type(d1o3_type::RIGHT_BIASED));
      d1_driver("*** 1st Derivative, 3rd Order Right Uniform Spacing ***", dx, d1o3_type(d1o3_type::RIGHT));

      std::vector<data__> dx_factor(4); dx_factor[0]=static_cast<data__>(1.0); dx_factor[1]=static_cast<data__>(1.0); dx_factor[2]=static_cast<data__>(1.0); dx_factor[3]=static_cast<data__>(1.0);
      d1_driver("*** 1st Derivative, 3rd Order Left Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::LEFT));
      d1_driver("*** 1st Derivative, 3rd Order Left-Biased Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::LEFT_BIASED));
      d1_driver("*** 1st Derivative, 3rd Order Right-Biased Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::RIGHT_BIASED));
      d1_driver("*** 1st Derivative, 3rd Order Right Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::RIGHT));
      dx_factor[0]=static_cast<data__>(2.95/3); dx_factor[1]=static_cast<data__>(1.025); dx_factor[2]=static_cast<data__>(0.95); dx_factor[3]=static_cast<data__>(1.00);
      d1_driver("*** 1st Derivative, 3rd Order Left Mildly Non-Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::RIGHT_BIASED));
      dx_factor[0]=static_cast<data__>(2.95/3); dx_factor[1]=static_cast<data__>(1.025); dx_factor[2]=static_cast<data__>(1.00); dx_factor[3]=static_cast<data__>(0.95);
      d1_driver("*** 1st Derivative, 3rd Order Left-Biased Mildly Non-Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::LEFT_BIASED));
      dx_factor[0]=static_cast<data__>(0.95); dx_factor[1]=static_cast<data__>(1.00); dx_factor[2]=static_cast<data__>(1.025); dx_factor[3]=static_cast<data__>(2.95/3);
      d1_driver("*** 1st Derivative, 3rd Order Right-Biased Mildly Non-Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::RIGHT_BIASED));
      dx_factor[0]=static_cast<data__>(1.00); dx_factor[1]=static_cast<data__>(0.95); dx_factor[2]=static_cast<data__>(1.025); dx_factor[3]=static_cast<data__>(2.95/3);
      d1_driver("*** 1st Derivative, 3rd Order Right Mildly Non-Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::RIGHT));
      dx_factor[0]=static_cast<data__>(1.2); dx_factor[1]=static_cast<data__>(0.95); dx_factor[2]=static_cast<data__>(0.1); dx_factor[3]=static_cast<data__>(1.00);
      d1_driver("*** 1st Derivative, 3rd Order Left Extremely Non-Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::LEFT));
      dx_factor[0]=static_cast<data__>(1.2); dx_factor[1]=static_cast<data__>(0.95); dx_factor[2]=static_cast<data__>(1.0); dx_factor[3]=static_cast<data__>(0.1);
      d1_driver("*** 1st Derivative, 3rd Order Left-Biased Extremely Non-Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::LEFT_BIASED));
      dx_factor[0]=static_cast<data__>(0.1); dx_factor[1]=static_cast<data__>(1.0); dx_factor[2]=static_cast<data__>(0.95); dx_factor[3]=static_cast<data__>(1.2);
      d1_driver("*** 1st Derivative, 3rd Order Right-Biased Extremely Non-Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::RIGHT_BIASED));
      dx_factor[0]=static_cast<data__>(1.0); dx_factor[1]=static_cast<data__>(0.1); dx_factor[2]=static_cast<data__>(0.95); dx_factor[3]=static_cast<data__>(1.2);
      d1_driver("*** 1st Derivative, 3rd Order Right Extremely Non-Uniform Spacing ***", dx_factor, d1o3_type(d1o3_type::RIGHT));
    }

    void d2o1_test()
    {
      typedef eli::mutil::fd::d2o1<data__> d2o1_type;

      std::vector<data__> dx(1); dx[0]=static_cast<data__>(1.0);
      d1_driver("*** 2nd Derivative, 1st Order Center Uniform Spacing ***", dx, d2o1_type(d2o1_type::CENTER));
      d1_driver("*** 2nd Derivative, 1st Order Left Uniform Spacing ***", dx, d2o1_type(d2o1_type::LEFT));
      d1_driver("*** 2nd Derivative, 1st Order Right Uniform Spacing ***", dx, d2o1_type(d2o1_type::RIGHT));

      std::vector<data__> dx_factor(3); dx_factor[0]=static_cast<data__>(1.0); dx_factor[1]=static_cast<data__>(1.0); dx_factor[2]=static_cast<data__>(1.0);
      d1_driver("*** 2nd Derivative, 1st Order Center Uniform Spacing ***", dx_factor, d2o1_type(d2o1_type::CENTER));
      d1_driver("*** 2nd Derivative, 1st Order Left Uniform Spacing ***", dx_factor, d2o1_type(d2o1_type::LEFT));
      d1_driver("*** 2nd Derivative, 1st Order Right Uniform Spacing ***", dx_factor, d2o1_type(d2o1_type::RIGHT));
      dx_factor[0]=static_cast<data__>(0.95); dx_factor[1]=static_cast<data__>(1.0); dx_factor[2]=static_cast<data__>(1.05);
      d1_driver("*** 2nd Derivative, 1st Order Center Mildly Non-Uniform Spacing ***", dx_factor, d2o1_type(d2o1_type::CENTER));
      dx_factor[0]=static_cast<data__>(1.025); dx_factor[1]=static_cast<data__>(0.95); dx_factor[2]=static_cast<data__>(1.00);
      d1_driver("*** 2nd Derivative, 1st Order Left Mildly Non-Uniform Spacing ***", dx_factor, d2o1_type(d2o1_type::LEFT));
      dx_factor[0]=static_cast<data__>(1.00); dx_factor[1]=static_cast<data__>(0.95); dx_factor[2]=static_cast<data__>(1.025);
      d1_driver("*** 2nd Derivative, 1st Order Right Mildly Non-Uniform Spacing ***", dx_factor, d2o1_type(d2o1_type::RIGHT));
      dx_factor[0]=static_cast<data__>(0.1); dx_factor[1]=static_cast<data__>(1.0); dx_factor[2]=static_cast<data__>(1.8);
      d1_driver("*** 2nd Derivative, 1st Order Center Extremely Non-Uniform Spacing ***", dx_factor, d2o1_type(d2o1_type::CENTER));
      dx_factor[0]=static_cast<data__>(0.95); dx_factor[1]=static_cast<data__>(0.1); dx_factor[2]=static_cast<data__>(1.00);
      d1_driver("*** 2nd Derivative, 1st Order Left Extremely Non-Uniform Spacing ***", dx_factor, d2o1_type(d2o1_type::LEFT));
      dx_factor[0]=static_cast<data__>(1.0); dx_factor[1]=static_cast<data__>(0.1); dx_factor[2]=static_cast<data__>(0.95);
      d1_driver("*** 2nd Derivative, 1st Order Right Extremely Non-Uniform Spacing ***", dx_factor, d2o1_type(d2o1_type::RIGHT));
    }

    void d2o2_test()
    {
      typedef eli::mutil::fd::d2o2<data__> d2o2_type;

      std::vector<data__> dx(1); dx[0]=static_cast<data__>(1.0);
      d1_driver("*** 2nd Derivative, 2nd Order Left Uniform Spacing ***", dx, d2o2_type(d2o2_type::LEFT));
      d1_driver("*** 2nd Derivative, 2nd Order Left-Biased Uniform Spacing ***", dx, d2o2_type(d2o2_type::LEFT_BIASED));
      d1_driver("*** 2nd Derivative, 2nd Order Right-Biased Uniform Spacing ***", dx, d2o2_type(d2o2_type::RIGHT_BIASED));
      d1_driver("*** 2nd Derivative, 2nd Order Right Uniform Spacing ***", dx, d2o2_type(d2o2_type::RIGHT));

      std::vector<data__> dx_factor(4); dx_factor[0]=static_cast<data__>(1.0); dx_factor[1]=static_cast<data__>(1.0); dx_factor[2]=static_cast<data__>(1.0); dx_factor[3]=static_cast<data__>(1.0);
      d1_driver("*** 2nd Derivative, 2nd Order Left Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::LEFT));
      d1_driver("*** 2nd Derivative, 2nd Order Left-Biased Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::LEFT_BIASED));
      d1_driver("*** 2nd Derivative, 2nd Order Right-Biased Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::RIGHT_BIASED));
      d1_driver("*** 2nd Derivative, 2nd Order Right Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::RIGHT));
      dx_factor[0]=static_cast<data__>(2.95/3); dx_factor[1]=static_cast<data__>(1.025); dx_factor[2]=static_cast<data__>(0.95); dx_factor[3]=static_cast<data__>(1.00);
      d1_driver("*** 2nd Derivative, 2nd Order Left Mildly Non-Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::LEFT));
      dx_factor[0]=static_cast<data__>(2.95/3); dx_factor[1]=static_cast<data__>(1.025); dx_factor[2]=static_cast<data__>(1.00); dx_factor[3]=static_cast<data__>(0.95);
      d1_driver("*** 2nd Derivative, 2nd Order Left-Biased Mildly Non-Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::LEFT_BIASED));
      dx_factor[0]=static_cast<data__>(0.95); dx_factor[1]=static_cast<data__>(1.00); dx_factor[2]=static_cast<data__>(1.025); dx_factor[3]=static_cast<data__>(2.95/3);
      d1_driver("*** 2nd Derivative, 2nd Order Right-Biased Mildly Non-Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::RIGHT_BIASED));
      dx_factor[0]=static_cast<data__>(1.00); dx_factor[1]=static_cast<data__>(0.95); dx_factor[2]=static_cast<data__>(1.025); dx_factor[3]=static_cast<data__>(2.95/3);
      d1_driver("*** 2nd Derivative, 2nd Order Right Mildly Non-Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::RIGHT));
      dx_factor[0]=static_cast<data__>(1.2); dx_factor[1]=static_cast<data__>(0.95); dx_factor[2]=static_cast<data__>(0.1); dx_factor[3]=static_cast<data__>(1.00);
      d1_driver("*** 2nd Derivative, 2nd Order Left Extremely Non-Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::LEFT));
      dx_factor[0]=static_cast<data__>(1.2); dx_factor[1]=static_cast<data__>(0.95); dx_factor[2]=static_cast<data__>(1.0); dx_factor[3]=static_cast<data__>(0.1);
      d1_driver("*** 2nd Derivative, 2nd Order Left-Biased Extremely Non-Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::LEFT_BIASED));
      dx_factor[0]=static_cast<data__>(0.1); dx_factor[1]=static_cast<data__>(1.0); dx_factor[2]=static_cast<data__>(0.95); dx_factor[3]=static_cast<data__>(1.2);
      d1_driver("*** 2nd Derivative, 2nd Order Right-Biased Extremely Non-Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::RIGHT_BIASED));
      dx_factor[0]=static_cast<data__>(1.0); dx_factor[1]=static_cast<data__>(0.1); dx_factor[2]=static_cast<data__>(0.95); dx_factor[3]=static_cast<data__>(1.2);
      d1_driver("*** 2nd Derivative, 2nd Order Right Extremely Non-Uniform Spacing ***", dx_factor, d2o2_type(d2o2_type::RIGHT));
    }

    template<typename d1__>
    void d1_driver(const std::string &str, const std::vector<data__> &dx_factor, const d1__ &dop)
    {
      const size_t nnodes(dop.number_nodes());
      const size_t nref(5);
      const data__ x0(0.5), rrat(3.0), dx0(0.5);
      data__ d(0);
      size_t c;
      std::ptrdiff_t ind0;
      std::vector<data__> phi(nref);
      std::vector<data__> phi_err(nref);
      std::vector<data__> phi_te(nref);
      std::vector<data__> dx(nref);
      std::vector<data__> ord(nref);
      std::vector<data__> x(nnodes);
      std::vector<data__> dx_term(nnodes);
      std::vector<std::ptrdiff_t> ind(nnodes);
      bool uniform;

      // capture the mesh refinement
      if (dx_factor.size() == 1)
      {
        uniform=true;
        for (size_t i=0; i<nnodes; ++i)
          dx_term[i]=static_cast<data__>(1.0);
      }
      else
      {
        dx_term=dx_factor;
        uniform=true;
        for (size_t i=1; i<dx_factor.size(); ++i)
        {
          if (dx_factor[i]!=dx_factor[i-1])
            uniform=false;
        }
      }

      // get the index offsets
      ind0=dop.index(ind.begin());

      for (c=0; c<nref; ++c)
      {
        dx[c]=dx0*std::pow(rrat,-static_cast<int>(c));
        for (size_t i=0; i<nnodes; ++i)
        {
          x[i]=x0+ind[i]*dx[c]*dx_term[i];
          phi[i]=std::exp(x[i]);
        }

        TEST_ASSERT(x0==x[ind0]);

        if (dx_factor.size() == 1)
        {
          int tertn;
          dop.evaluate(d, phi.begin(), dx[c]);
          tertn=dop.truncation_error(phi_te[c], std::exp(x[ind0]), dx[c]);
          if (tertn>0)
          {
            // need to evaluate the order+tertn derivative
            phi_te[c]*=std::exp(x[ind0]);
          }
        }
        else
        {
          dop.evaluate(d, phi.begin(), x.begin());
          dop.truncation_error(phi_te[c], std::exp(x[ind0]), x.begin());
        }
        phi_err[c]=std::exp(x[ind0])-d;
      }

      // calculate the order of accuracy
      calculate_order(ord.begin()+1, phi_err.begin(), phi_err.end(), rrat);
      if (typeid(data__)==typeid(float))
      {
        // for float, if truncation error is too small, then will hit numerical precision limits
        if (std::abs(phi_te[(nref-1)/2])<5e-3)
        {
          TEST_ASSERT(true);
        }
        else
        {
          TEST_ASSERT_DELTA_MSG(ord[(nref-1)/2], dop.order(uniform), static_cast<data__>(5.0e-1), str.c_str());
        }
      }
      else
      {
        TEST_ASSERT_DELTA_MSG(ord[nref-1], dop.order(uniform), static_cast<data__>(5.0e-2), str.c_str());
      }
    }
};

#endif

