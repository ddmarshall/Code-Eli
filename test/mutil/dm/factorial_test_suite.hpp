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

#ifndef dm_factorial_test_suite_hpp
#define dm_factorial_test_suite_hpp

#include <typeinfo>   // typeid
#include <string>     // std::string
#include <sstream>    // std::stringstream
#include <iomanip>    // std::setw
#include <vector>     // std::vector
#include <functional> // std::less

#include "eli/mutil/dm/factorial.hpp"

#include <string>
#include <vector>

class factorial_test_suite : public Test::Suite
{
  protected:
    void AddTests()
    {
      TEST_ADD(factorial_test_suite::factorial_test);
    }

  public:
    factorial_test_suite()
    {
      // add the tests
      AddTests();
    }
    ~factorial_test_suite()
    {
    }

  private:
    void factorial_test()
    {
      int iv;
      double dv;

      eli::mutil::dm::factorial(iv, 6);
      TEST_ASSERT(iv==720);

      eli::mutil::dm::factorial(dv, 6);
      TEST_ASSERT(dv==720.0);

      eli::mutil::dm::factorial(iv, 12);
      TEST_ASSERT(iv==479001600);

      eli::mutil::dm::factorial(dv, 12);
      TEST_ASSERT(dv==479001600.0);
    }
};

#endif

