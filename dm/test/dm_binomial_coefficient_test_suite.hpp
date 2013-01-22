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

#include "code_eli.hpp"
#include "dm/factorial.hpp"
#include "dm/binomial_coefficient.hpp"

#include <string>
#include <vector>

class dm_binomial_coefficient_test_suite : public Test::Suite
{
  protected:
    void AddTests()
    {
      TEST_ADD(dm_binomial_coefficient_test_suite::nchoosek_test);
      TEST_ADD(dm_binomial_coefficient_test_suite::binomial_coefficient_test);
    }

  public:
    dm_binomial_coefficient_test_suite()
    {
      // add the tests
      AddTests();
    }
    ~dm_binomial_coefficient_test_suite()
    {
    }

  private:
    void nchoosek_test()
    {
      int iv;
      double dv;

      eli::dm::n_choose_k(iv, 10, 7);
      TEST_ASSERT(iv==120);

      eli::dm::n_choose_k(dv, 10, 7);
      TEST_ASSERT(dv==120.0);

      eli::dm::n_choose_k(iv, 20, 11);
      TEST_ASSERT(iv==167960);

      eli::dm::n_choose_k(dv, 20, 11);
      TEST_ASSERT(dv==167960.0);
    }

    void binomial_coefficient_test()
    {
      double dv;

      eli::dm::binomial_coefficient(dv, 10, 7);
      TEST_ASSERT(dv==120.0);

      eli::dm::binomial_coefficient(dv, 20, 11);
      TEST_ASSERT(dv==167960.0);

      eli::dm::binomial_coefficient(dv, 10.5, 7);
      TEST_ASSERT(std::abs(dv-202.98)<1e-2);
    }
};

#endif

