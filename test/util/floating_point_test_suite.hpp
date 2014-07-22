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

#ifndef floating_point_test_suite_hpp
#define floating_point_test_suite_hpp

#include <string>
#include <sstream>

#include "Eigen/Eigen"

#include "eli/code_eli.hpp"

#include "eli/util/floating_point.hpp"

template<typename data__>
class floating_point_test_suite : public Test::Suite
{
  protected:
    void AddTests(const float &)
    {
      TEST_ADD(floating_point_test_suite<float>::printing_test);
      TEST_ADD(floating_point_test_suite<float>::increment_ulp_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(floating_point_test_suite<double>::printing_test);
      TEST_ADD(floating_point_test_suite<double>::increment_ulp_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(floating_point_test_suite<long double>::printing_test);
      TEST_ADD(floating_point_test_suite<long double>::increment_ulp_test);
    }

  public:
    floating_point_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~floating_point_test_suite()
    {
    }

  private:
    typedef data__ data_type;

  private:
// TODO: Need to have methods that create data for each test for each data type

    void printing_test()
    {
      printing_test(data_type());
    }

    void printing_test(const float &)
    {
      float f(2);
      const float cf(f);
      eli::util::float_type *pft=eli::util::set_floating_point_type(&f);
      const eli::util::float_type * pcft=eli::util::set_floating_point_type(&cf);
      std::stringstream ostr;
      std::string ref;

      ref="0x0 0x0 0x80";
      ostr << (*pft);
      TEST_ASSERT(ostr.str()==ref);
      ostr.str("");
      ostr << (*pcft);
      TEST_ASSERT(ostr.str()==ref);
      ostr.str("");
    }

    void printing_test(const double &)
    {
      double f(2);
      const double cf(f);
      eli::util::double_type *pft=eli::util::set_floating_point_type(&f);
      const eli::util::double_type * pcft=eli::util::set_floating_point_type(&cf);
      std::stringstream ostr;
      std::string ref;

      ref="0x0 0x0 0x400";
      ostr << (*pft);
      TEST_ASSERT(ostr.str()==ref);
      ostr.str("");
      ostr << (*pcft);
      TEST_ASSERT(ostr.str()==ref);
      ostr.str("");
    }

    void printing_test(const long double &)
    {
      long double f(2);
      const long double cf(f);
      eli::util::long_double_type *pft=eli::util::set_floating_point_type(&f);
      const eli::util::long_double_type * pcft=eli::util::set_floating_point_type(&cf);
      std::stringstream ostr;
      std::string ref;

#if defined(_WIN32)
      ref="0x0 0x0 0x400";
#else
      ref="0x0 0x1 0x0 0x4000";
#endif
      ostr << (*pft);
      TEST_ASSERT(ostr.str()==ref);
      ostr.str("");
      ostr << (*pcft);
      TEST_ASSERT(ostr.str()==ref);
      ostr.str("");
    }

    void increment_ulp_test()
    {
      // increment from zero by one ulp
      {
        data_type f, f_inc, f_inc_ref;

        get_zero_and_one_ulp(f, f_inc_ref);
        f_inc=eli::util::increment_ulp(f, 1);
        TEST_ASSERT(f_inc==f_inc_ref);
      }

      // increment from one by one ulp
      {
        data_type f, f_inc, f_inc_ref;

        get_one_and_many_ulp(f, f_inc_ref, 1);
        f_inc=eli::util::increment_ulp(f, 1);
        TEST_ASSERT(f_inc==f_inc_ref);
      }

      // increment from one by multiple ulps
      {
        data_type f, f_inc, f_inc_ref;

        get_one_and_many_ulp(f, f_inc_ref, 1000);
        f_inc=eli::util::increment_ulp(f, 1000);
        TEST_ASSERT(f_inc==f_inc_ref);
      }


      // increment from largest mantissa by one ulp
      {
        data_type f, f_inc, f_inc_ref;

        get_max_mantissa_and_many_ulp(f, f_inc_ref, 1);
        f_inc=eli::util::increment_ulp(f, 1);
        TEST_ASSERT(f_inc==f_inc_ref);
      }

      // increment from largest mantissa by multiple ulps
      {
        data_type f, f_inc, f_inc_ref;

        get_max_mantissa_and_many_ulp(f, f_inc_ref, 1000);
        f_inc=eli::util::increment_ulp(f, 1000);
        TEST_ASSERT(f_inc==f_inc_ref);
      }
    }

    void get_zero_and_one_ulp(float &f, float &f_ulp)
    {
      int32_t *pi;
      f=0;
      f_ulp=f;

      pi=static_cast<int32_t *>(static_cast<void *>(&f_ulp));
      ++(*pi);
    }

    void get_zero_and_one_ulp(double &f, double &f_ulp)
    {
      int64_t *pi;
      f=0;
      f_ulp=f;

      pi=static_cast<int64_t *>(static_cast<void *>(&f_ulp));
      ++(*pi);
    }

    void get_zero_and_one_ulp(long double &f, long double &f_ulp)
    {
      int64_t *pi;
      f=0;
      f_ulp=f;

      pi=static_cast<int64_t *>(static_cast<void *>(&f_ulp));
      ++(*pi);
    }

    void get_one_and_many_ulp(data_type &f, data_type &f_ulp, int ulps)
    {
      f=1;
      f_ulp=f+ulps*std::numeric_limits<data_type>::epsilon();
    }

    void get_max_mantissa_and_many_ulp(data_type &f, data_type &f_ulp, int ulps)
    {
      f_ulp=2;
      f=2-ulps*std::numeric_limits<data_type>::epsilon();
    }
};

#endif
