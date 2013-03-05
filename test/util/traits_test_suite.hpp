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

#ifndef traits_test_suite_hpp
#define traits_test_suite_hpp

#include "Eigen/Eigen"

#include "eli/code_eli.hpp"

#include "eli/util/traits.hpp"

#include <limits>     // numeric_limits

template<typename data__>
class traits_test_suite : public Test::Suite
{
  protected:
    void AddTests(const float &)
    {
      TEST_ADD(traits_test_suite<float>::cast_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(traits_test_suite<double>::cast_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(traits_test_suite<long double>::cast_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      TEST_ADD(traits_test_suite<dd_real>::cast_test);
    }

    void AddTests(const qd_real &)
    {
      TEST_ADD(traits_test_suite<qd_real>::cast_test);
    }
#endif

  public:
    traits_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~traits_test_suite()
    {
    }

  private:
    typedef data__ data_type;

  private:
    void cast_test()
    {
    }
};

#endif
