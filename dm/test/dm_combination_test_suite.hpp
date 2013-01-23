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

#ifndef dm_combination_test_suite_hpp
#define dm_combination_test_suite_hpp

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo>   // typeid
#include <string>     // std::string
#include <sstream>    // std::stringstream
#include <iomanip>    // std::setw
#include <vector>     // std::vector
#include <functional> // std::less

#include "eli/code_eli.hpp"
#include "eli/dm/combination.hpp"

#include <string>
#include <vector>

class dm_combination_test_suite : public Test::Suite
{
  protected:
    void AddTests()
    {
      TEST_ADD(dm_combination_test_suite::string_test);
      TEST_ADD(dm_combination_test_suite::vector_test);
    }

  public:
    dm_combination_test_suite()
    {
      // add the tests
      AddTests();
    }
    ~dm_combination_test_suite()
    {
    }

  private:
    void string_test()
    {
      std::string s = "1234";
      std::size_t comb_size = 3;

      TEST_ASSERT(eli::dm::next_combination(s.begin(),s.begin() + comb_size,s.end()));
      TEST_ASSERT(s=="1243");
      TEST_ASSERT(eli::dm::next_combination(s.begin(),s.begin() + comb_size,s.end()));
      TEST_ASSERT(s=="1342");
      TEST_ASSERT(eli::dm::next_combination(s.begin(),s.begin() + comb_size,s.end()));
      TEST_ASSERT(s=="2341");
      TEST_ASSERT(!eli::dm::next_combination(s.begin(),s.begin() + comb_size,s.end()));
      TEST_ASSERT(s=="1234");
    }

    void vector_test()
    {
      size_t i;
      std::vector<int> s(4);
      std::size_t comb_size = 3;

      // set s values
      for (i=0; i<s.size(); ++i)
        s[i]=i+1;

      TEST_ASSERT(eli::dm::next_combination(s.begin(),s.begin() + comb_size,s.end(), std::less<int>()));
      TEST_ASSERT((s[0]==1) && (s[1]==2) && (s[2]==4) && (s[3]==3));
      TEST_ASSERT(eli::dm::next_combination(s.begin(),s.begin() + comb_size,s.end(), std::less<int>()));
      TEST_ASSERT((s[0]==1) && (s[1]==3) && (s[2]==4) && (s[3]==2));
      TEST_ASSERT(eli::dm::next_combination(s.begin(),s.begin() + comb_size,s.end(), std::less<int>()));
      TEST_ASSERT((s[0]==2) && (s[1]==3) && (s[2]==4) && (s[3]==1));
      TEST_ASSERT(!eli::dm::next_combination(s.begin(),s.begin() + comb_size,s.end(), std::less<int>()));
      TEST_ASSERT((s[0]==1) && (s[1]==2) && (s[2]==3) && (s[3]==4));
    }
};

#endif

