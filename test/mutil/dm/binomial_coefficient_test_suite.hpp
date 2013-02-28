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

#ifndef dm_binomial_coefficient_test_suite_hpp
#define dm_binomial_coefficient_test_suite_hpp

#include <cassert>  // assert()

#include <typeinfo>   // typeid
#include <string>     // std::string
#include <sstream>    // std::stringstream
#include <iomanip>    // std::setw
#include <vector>     // std::vector
#include <functional> // std::less

#include "eli/code_eli.hpp"

#include "eli/mutil/dm/factorial.hpp"
#include "eli/mutil/dm/binomial_coefficient.hpp"

#include <string>
#include <vector>

class binomial_coefficient_test_suite : public Test::Suite
{
  protected:
    void AddTests()
    {
      TEST_ADD(binomial_coefficient_test_suite::nchoosek_test);
      TEST_ADD(binomial_coefficient_test_suite::binomial_coefficient_test);
    }

  public:
    binomial_coefficient_test_suite()
    {
      // add the tests
      AddTests();
    }
    ~binomial_coefficient_test_suite()
    {
    }

  private:
    void nchoosek_test()
    {
      int iv;
      double dv;

      eli::mutil::dm::n_choose_k(iv, 10, 7);
      TEST_ASSERT(iv==120);

      eli::mutil::dm::n_choose_k(dv, 10, 7);
#if defined(__INTEL_COMPILER) && defined(NDEBUG)
      TEST_ASSERT(std::abs(1-dv/120)<std::numeric_limits<double>::epsilon());
#else
      TEST_ASSERT(dv==120);
#endif

      eli::mutil::dm::n_choose_k(iv, 20, 11);
      TEST_ASSERT(iv==167960);

      eli::mutil::dm::n_choose_k(dv, 20, 11);
      TEST_ASSERT(dv==167960);
    }

    void binomial_coefficient_test()
    {
      double dv;

      eli::mutil::dm::binomial_coefficient(dv, 10, 7);
#if defined(__INTEL_COMPILER) && defined(NDEBUG)
      TEST_ASSERT(std::abs(1-dv/120)<std::numeric_limits<double>::epsilon());
#else
      TEST_ASSERT(dv==120);
#endif

      eli::mutil::dm::binomial_coefficient(dv, 20, 11);
#if defined(__INTEL_COMPILER) && defined(NDEBUG)
      TEST_ASSERT(std::abs(1-dv/167960)<2*std::numeric_limits<double>::epsilon());
#else
      TEST_ASSERT(dv==167960);
#endif

      eli::mutil::dm::binomial_coefficient(dv, 10.5, 7);
      TEST_ASSERT(std::abs(dv-202.98)<1e-2);
    }
};

#endif

