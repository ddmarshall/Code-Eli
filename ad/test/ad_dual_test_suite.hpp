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

#ifndef ad_dual_test_suite_hpp
#define ad_dual_test_suite_hpp

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/ad/dual.hpp"

template<typename data__>
class ad_dual_test_suite : public Test::Suite
{
  protected:
    void AddTests(const float &)
    {
      TEST_ADD(ad_dual_test_suite<float>::dual_assignment_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_addition_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_subtraction_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_multiplication_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_division_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_exact_comparison_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_inexact_comparison_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_trig_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_trig_nostd_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_hyp_trig_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_hyp_trig_nostd_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_exp_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_exp_nostd_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_power_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_power_nostd_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_misc_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_misc_nostd_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_erf_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_erf_nostd_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_gamma_test);
//       TEST_ADD(ad_dual_test_suite<float>::dual_gamma_nostd_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(ad_dual_test_suite<double>::dual_assignment_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_addition_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_subtraction_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_multiplication_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_division_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_exact_comparison_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_inexact_comparison_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_trig_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_trig_nostd_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_hyp_trig_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_hyp_trig_nostd_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_exp_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_exp_nostd_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_power_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_power_nostd_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_misc_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_misc_nostd_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_erf_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_erf_nostd_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_gamma_test);
//       TEST_ADD(ad_dual_test_suite<double>::dual_gamma_nostd_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(ad_dual_test_suite<long double>::dual_assignment_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_addition_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_subtraction_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_multiplication_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_division_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_exact_comparison_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_inexact_comparison_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_trig_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_trig_nostd_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_hyp_trig_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_hyp_trig_nostd_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_exp_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_exp_nostd_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_power_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_power_nostd_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_misc_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_misc_nostd_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_erf_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_erf_nostd_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_gamma_test);
//       TEST_ADD(ad_dual_test_suite<long double>::dual_gamma_nostd_test);
    }

#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      TEST_ADD(ad_dual_test_suite<dd_real>::dual_assignment_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_addition_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_subtraction_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_multiplication_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_division_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_exact_comparison_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_inexact_comparison_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_trig_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_trig_nostd_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_hyp_trig_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_hyp_trig_nostd_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_exp_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_exp_nostd_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_power_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_power_nostd_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_misc_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_misc_nostd_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_erf_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_erf_nostd_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_gamma_test);
//       TEST_ADD(ad_dual_test_suite<dd_real>::dual_gamma_nostd_test);
    }

    void AddTests(const qd_real &)
    {
      TEST_ADD(ad_dual_test_suite<qd_real>::dual_assignment_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_addition_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_subtraction_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_multiplication_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_division_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_exact_comparison_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_inexact_comparison_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_trig_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_trig_nostd_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_hyp_trig_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_hyp_trig_nostd_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_exp_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_exp_nostd_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_power_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_power_nostd_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_misc_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_misc_nostd_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_erf_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_erf_nostd_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_gamma_test);
//       TEST_ADD(ad_dual_test_suite<qd_real>::dual_gamma_nostd_test);
    }
#endif

  public:
    ad_dual_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~ad_dual_test_suite()
    {
    }

  private:
    void dual_assignment_test()
    {
      eli::ad::dual<data__, false> d1, d2(1,2), d4(-2,4), d5;
      data__ v1(2);

      // test copy ctr
      eli::ad::dual<data__, false> d3(d2);
      TEST_ASSERT(d3.exact(d2));

      // test constructor with expression
      eli::ad::dual<data__, false> de(d2+d4);
      d5.set_real(d2.real()+d4.real());
      d5.set_nonreal(d2.nonreal()+d4.nonreal());

      // test assignment operator
      d1=d2;
      TEST_ASSERT(d1.exact(d2));

      // test A=v
      d5.set_real(v1);
      d5.set_nonreal(0.0);
      d4=v1;
      TEST_ASSERT(d4.exact(d5));

      // test A=int
      d5.set_real(2);
      d5.set_nonreal(0.0);
      d4=2;
      TEST_ASSERT(d4.exact(d5));
    }

    void dual_exact_comparison_test()
    {
      eli::ad::dual<data__, false> d1(1,3), d2(1,2), d3(1,2), d4(-2,4), d5(2,4), d6(2,0), d7(1,-3), d8(3,-1);
      data__ v1(2);

      // compare dual & dual
      TEST_ASSERT(d2==d3); // both equal
      TEST_ASSERT(d1!=d3); // real equal
      TEST_ASSERT(d4!=d5); // nonreal equal
      TEST_ASSERT(d1!=d4); // neither equal

      // compare dual & v
      TEST_ASSERT(d6==v1); // equal real, zero nonreal
      TEST_ASSERT(d5!=v1); // equal real, nonzero nonreal
      TEST_ASSERT(d1!=v1); // neither equal
      TEST_ASSERT(v1==d6); // equal real, zero nonreal
      TEST_ASSERT(v1!=d5); // equal real, nonzero nonreal
      TEST_ASSERT(v1!=d1); // neither equal

      // compare dual & int
      TEST_ASSERT(d6==2); // equal real, zero nonreal
      TEST_ASSERT(d5!=2); // equal real, nonzero nonreal
      TEST_ASSERT(d1!=2); // neither equal
      TEST_ASSERT(2==d6); // equal real, zero nonreal
      TEST_ASSERT(2!=d5); // equal real, nonzero nonreal
      TEST_ASSERT(2!=d1); // neither equal

      // compare dual & expression
      TEST_ASSERT(d5==(d2+d3)); // both equal
      TEST_ASSERT(d5!=(d1+d3)); // real equal
      TEST_ASSERT(d4!=(d5+d6)); // nonreal equal
      TEST_ASSERT(d5!=(d1+d4)); // neither equal
      TEST_ASSERT((d2+d3)==d5); // both equal
      TEST_ASSERT((d1+d3)!=d5); // real equal
      TEST_ASSERT((d5+d6)!=d4); // nonreal equal
      TEST_ASSERT((d1+d4)!=d5); // neither equal

      // compare expression & v
      TEST_ASSERT((d7+d1)==v1); // equal real, zero nonreal
      TEST_ASSERT((d1+d2)!=v1); // equal real, nonzero nonreal
      TEST_ASSERT((d1+d4)!=v1); // neither equal
      TEST_ASSERT(v1==(d7+d1)); // equal real, zero nonreal
      TEST_ASSERT(v1!=(d1+d2)); // equal real, nonzero nonreal
      TEST_ASSERT(v1!=(d1+d4)); // neither equal

      // compare expression & int
      TEST_ASSERT((d7+d1)==2); // equal real, zero nonreal
      TEST_ASSERT((d1+d2)!=2); // equal real, nonzero nonreal
      TEST_ASSERT((d1+d4)!=2); // neither equal
      TEST_ASSERT(2==(d7+d1)); // equal real, zero nonreal
      TEST_ASSERT(2!=(d1+d2)); // equal real, nonzero nonreal
      TEST_ASSERT(2!=(d1+d4)); // neither equal

      // compare expression & expression
      TEST_ASSERT((d2+d3)==(d2+d3)); // both equal
      TEST_ASSERT((d2+d3)!=(d1+d3)); // real equal
      TEST_ASSERT((d4+d6)!=(d5+d6)); // nonreal equal
      TEST_ASSERT((d1+d5)!=(d1+d4)); // neither equal

      // TODO: DO NOT KNOW HOW TO DO GREATER/LESS THAN WITH EXACT COMPARISONS
    }

    void dual_inexact_comparison_test()
    {
      eli::ad::dual<data__, true> d1(1,3), d2(1,2), d3(1,2), d4(-2,4), d5(2,4), d6(2,0), d7(1,-3), d8(3,-1);
      data__ v1(2);

      // compare dual & dual
      TEST_ASSERT(d2==d3); // both equal
      TEST_ASSERT(d1==d3); // real equal
      TEST_ASSERT(d4!=d5); // nonreal equal
      TEST_ASSERT(d1!=d4); // neither equal

      // compare dual & v
      TEST_ASSERT(d6==v1); // equal real, zero nonreal
      TEST_ASSERT(d5==v1); // equal real, nonzero nonreal
      TEST_ASSERT(d1!=v1); // neither equal
      TEST_ASSERT(v1==d6); // equal real, zero nonreal
      TEST_ASSERT(v1==d5); // equal real, nonzero nonreal
      TEST_ASSERT(v1!=d1); // neither equal

      // compare dual & int
      TEST_ASSERT(d6==2); // equal real, zero nonreal
      TEST_ASSERT(d5==2); // equal real, nonzero nonreal
      TEST_ASSERT(d1!=2); // neither equal
      TEST_ASSERT(2==d6); // equal real, zero nonreal
      TEST_ASSERT(2==d5); // equal real, nonzero nonreal
      TEST_ASSERT(2!=d1); // neither equal

      // compare dual & expression
      TEST_ASSERT(d5==(d2+d3)); // both equal
      TEST_ASSERT(d5==(d1+d3)); // real equal
      TEST_ASSERT(d4!=(d5+d6)); // nonreal equal
      TEST_ASSERT(d5!=(d1+d4)); // neither equal
      TEST_ASSERT((d2+d3)==d5); // both equal
      TEST_ASSERT((d1+d3)==d5); // real equal
      TEST_ASSERT((d5+d6)!=d4); // nonreal equal
      TEST_ASSERT((d1+d4)!=d5); // neither equal

      // compare expression & v
      TEST_ASSERT((d7+d1)==v1); // equal real, zero nonreal
      TEST_ASSERT((d1+d2)==v1); // equal real, nonzero nonreal
      TEST_ASSERT((d1+d4)!=v1); // neither equal
      TEST_ASSERT(v1==(d7+d1)); // equal real, zero nonreal
      TEST_ASSERT(v1==(d1+d2)); // equal real, nonzero nonreal
      TEST_ASSERT(v1!=(d1+d4)); // neither equal

      // compare expression & INT
      TEST_ASSERT((d7+d1)==2); // equal real, zero nonreal
      TEST_ASSERT((d1+d2)==2); // equal real, nonzero nonreal
      TEST_ASSERT((d1+d4)!=2); // neither equal
      TEST_ASSERT(2==(d7+d1)); // equal real, zero nonreal
      TEST_ASSERT(2==(d1+d2)); // equal real, nonzero nonreal
      TEST_ASSERT(2!=(d1+d4)); // neither equal

      // compare expression & expression
      TEST_ASSERT((d2+d3)==(d2+d3)); // both equal
      TEST_ASSERT((d2+d3)==(d1+d3)); // real equal
      TEST_ASSERT((d4+d6)!=(d5+d6)); // nonreal equal
      TEST_ASSERT((d1+d5)!=(d1+d4)); // neither equal

      // compare dual & dual
      TEST_ASSERT(d2<=d3); // both equal
      TEST_ASSERT(d1<=d3); // real equal
      TEST_ASSERT(d4<=d5); // nonreal equal
      TEST_ASSERT(d4<=d1); // neither equal

      // compare dual & v
      TEST_ASSERT(d6<=v1); // equal real, zero nonreal
      TEST_ASSERT(d5<=v1); // equal real, nonzero nonreal
      TEST_ASSERT(d1<=v1); // neither equal
      TEST_ASSERT(v1<=d6); // equal real, zero nonreal
      TEST_ASSERT(v1<=d5); // equal real, nonzero nonreal
      TEST_ASSERT(v1<=d8); // neither equal

      // compare dual & int
      TEST_ASSERT(d6<=2); // equal real, zero nonreal
      TEST_ASSERT(d5<=2); // equal real, nonzero nonreal
      TEST_ASSERT(d1<=2); // neither equal
      TEST_ASSERT(2<=d6); // equal real, zero nonreal
      TEST_ASSERT(2<=d5); // equal real, nonzero nonreal
      TEST_ASSERT(2<=d8); // neither equal

      // compare dual & expression
      TEST_ASSERT(d5<=(d2+d3)); // both equal
      TEST_ASSERT(d5<=(d1+d3)); // real equal
      TEST_ASSERT(d4<=(d5+d6)); // nonreal equal
      TEST_ASSERT(d4<=(d1+d2)); // neither equal
      TEST_ASSERT((d2+d3)<=d5); // both equal
      TEST_ASSERT((d1+d3)<=d5); // real equal
      TEST_ASSERT((d2-d3)<=d6); // nonreal equal
      TEST_ASSERT((d3+d4)<=d5); // neither equal

      // compare expression & v
      TEST_ASSERT((d7+d1)<=v1); // equal real, zero nonreal
      TEST_ASSERT((d1+d2)<=v1); // equal real, nonzero nonreal
      TEST_ASSERT((d1+d4)<=v1); // neither equal
      TEST_ASSERT(v1<=(d7+d1)); // equal real, zero nonreal
      TEST_ASSERT(v1<=(d1+d2)); // equal real, nonzero nonreal
      TEST_ASSERT(v1<=(d5+d6)); // neither equal

      // compare expression & INT
      TEST_ASSERT((d7+d1)<=2); // equal real, zero nonreal
      TEST_ASSERT((d1+d2)<=2); // equal real, nonzero nonreal
      TEST_ASSERT((d1+d4)<=2); // neither equal
      TEST_ASSERT(2<=(d7+d1)); // equal real, zero nonreal
      TEST_ASSERT(2<=(d1+d2)); // equal real, nonzero nonreal
      TEST_ASSERT(2<=(d5+d6)); // neither equal

      // compare expression & expression
      TEST_ASSERT((d2+d3)<=(d2+d3)); // both equal
      TEST_ASSERT((d2+d3)<=(d1+d3)); // real equal
      TEST_ASSERT((d4+d6)<=(d5+d6)); // nonreal equal
      TEST_ASSERT((d1+d4)<=(d1+d5)); // neither equal

      // compare dual & dual
      TEST_ASSERT(d3>=d2); // both equal
      TEST_ASSERT(d3>=d1); // real equal
      TEST_ASSERT(d5>=d4); // nonreal equal
      TEST_ASSERT(d1>=d4); // neither equal

      // compare dual & v
      TEST_ASSERT(d6>=v1); // equal real, zero nonreal
      TEST_ASSERT(d5>=v1); // equal real, nonzero nonreal
      TEST_ASSERT(d8>=v1); // neither equal
      TEST_ASSERT(v1>=d6); // equal real, zero nonreal
      TEST_ASSERT(v1>=d5); // equal real, nonzero nonreal
      TEST_ASSERT(v1>=d1); // neither equal

      // compare dual & int
      TEST_ASSERT(d6>=2); // equal real, zero nonreal
      TEST_ASSERT(d5>=2); // equal real, nonzero nonreal
      TEST_ASSERT(d8>=2); // neither equal
      TEST_ASSERT(2>=d6); // equal real, zero nonreal
      TEST_ASSERT(2>=d5); // equal real, nonzero nonreal
      TEST_ASSERT(2>=d1); // neither equal

      // compare dual & expression
      TEST_ASSERT(d5>=(d2+d3)); // both equal
      TEST_ASSERT(d5>=(d1+d3)); // real equal
      TEST_ASSERT(d6>=(d2-d3)); // nonreal equal
      TEST_ASSERT(d5>=(d3+d4)); // neither equal
      TEST_ASSERT((d2+d3)>=d5); // both equal
      TEST_ASSERT((d1+d3)>=d5); // real equal
      TEST_ASSERT((d5+d6)>=d4); // nonreal equal
      TEST_ASSERT((d1+d2)>=d4); // neither equal

      // compare expression & v
      TEST_ASSERT((d7+d1)>=v1); // equal real, zero nonreal
      TEST_ASSERT((d1+d2)>=v1); // equal real, nonzero nonreal
      TEST_ASSERT((d5+d6)>=v1); // neither equal
      TEST_ASSERT(v1>=(d7+d1)); // equal real, zero nonreal
      TEST_ASSERT(v1>=(d1+d2)); // equal real, nonzero nonreal
      TEST_ASSERT(v1>=(d1+d4)); // neither equal

      // compare expression & INT
      TEST_ASSERT((d7+d1)>=2); // equal real, zero nonreal
      TEST_ASSERT((d1+d2)>=2); // equal real, nonzero nonreal
      TEST_ASSERT((d5+d6)>=2); // neither equal
      TEST_ASSERT(2>=(d7+d1)); // equal real, zero nonreal
      TEST_ASSERT(2>=(d1+d2)); // equal real, nonzero nonreal
      TEST_ASSERT(2>=(d1+d4)); // neither equal

      // compare expression & expression
      TEST_ASSERT((d2+d3)>=(d2+d3)); // both equal
      TEST_ASSERT((d2+d3)>=(d1+d3)); // real equal
      TEST_ASSERT((d5+d6)>=(d4+d6)); // nonreal equal
      TEST_ASSERT((d1+d5)>=(d1+d4)); // neither equal

      // compare dual & dual
      TEST_ASSERT(d4<d5); // nonreal equal
      TEST_ASSERT(d4<d1); // neither equal

      // compare dual & v
      TEST_ASSERT(d1<v1); // neither equal
      TEST_ASSERT(v1<d8); // neither equal

      // compare dual & int
      TEST_ASSERT(d1<2); // neither equal
      TEST_ASSERT(2<d8); // neither equal

      // compare dual & expression
      TEST_ASSERT(d4<(d5+d6)); // nonreal equal
      TEST_ASSERT(d4<(d1+d2)); // neither equal
      TEST_ASSERT((d2-d3)<d6); // nonreal equal
      TEST_ASSERT((d3+d4)<d5); // neither equal

      // compare expression & v
      TEST_ASSERT((d1+d4)<v1); // neither equal
      TEST_ASSERT(v1<(d5+d6)); // neither equal

      // compare expression & INT
      TEST_ASSERT((d1+d4)<2); // neither equal
      TEST_ASSERT(2<(d5+d6)); // neither equal

      // compare expression & expression
      TEST_ASSERT((d4+d6)<(d5+d6)); // nonreal equal
      TEST_ASSERT((d1+d4)<(d1+d5)); // neither equal

      // compare dual & dual
      TEST_ASSERT(d5>d4); // nonreal equal
      TEST_ASSERT(d1>d4); // neither equal

      // compare dual & v
      TEST_ASSERT(d8>v1); // neither equal
      TEST_ASSERT(v1>d1); // neither equal

      // compare dual & int
      TEST_ASSERT(d8>2); // neither equal
      TEST_ASSERT(2>d1); // neither equal

      // compare dual & expression
      TEST_ASSERT(d6>(d2-d3)); // nonreal equal
      TEST_ASSERT(d5>(d3+d4)); // neither equal
      TEST_ASSERT((d5+d6)>d4); // nonreal equal
      TEST_ASSERT((d1+d2)>d4); // neither equal

      // compare expression & v
      TEST_ASSERT((d5+d6)>v1); // neither equal
      TEST_ASSERT(v1>(d1+d4)); // neither equal

      // compare expression & INT
      TEST_ASSERT((d5+d6)>2); // neither equal
      TEST_ASSERT(2>(d1+d4)); // neither equal

      // compare expression & expression
      TEST_ASSERT((d5+d6)>(d4+d6)); // nonreal equal
      TEST_ASSERT((d1+d5)>(d1+d4)); // neither equal
    }

    void dual_addition_test()
    {
      eli::ad::dual<data__, false> d1(3,-4), d2(1,2), d3(d2), d4, d5;
      data__ v1(2);

      // test C=A+B
      d5.set_real(d2.real()+d3.real());
      d5.set_nonreal(d2.nonreal()+d3.nonreal());
      d4=d2+d3;
      TEST_ASSERT(d4.exact(d5));

      // test B=A+B
      d4=d3;
      d5.set_real(d2.real()+d4.real());
      d5.set_nonreal(d2.nonreal()+d4.nonreal());
      d4=d2+d4;
      TEST_ASSERT(d4.exact(d5));

      // test B+=A
      d4=d3;
      d5.set_real(d2.real()+d4.real());
      d5.set_nonreal(d2.nonreal()+d4.nonreal());
      d4+=d2;
      TEST_ASSERT(d4.exact(d5));

      // test B+=v
      d4=d3;
      d5.set_real(d2.real()+v1);
      d5.set_nonreal(d2.nonreal()+0.0);
      d4+=v1;
      TEST_ASSERT(d4.exact(d5));

       // test B+=int
      d4=d3;
      d5.set_real(d2.real()+2);
      d5.set_nonreal(d2.nonreal()+0.0);
      d4+=2;
      TEST_ASSERT(d4.exact(d5));

      // test C+=A+B
      d4=d3;
      d5.set_real(d4.real()+(d1.real()+d2.real()));
      d5.set_nonreal(d4.nonreal()+(d1.nonreal()+d2.nonreal()));
      d4+=d1+d2;
      TEST_ASSERT(d4.exact(d5));

      // test D=(A+B)+C
      d5.set_real((d1.real()+d2.real())+d3.real());
      d5.set_nonreal((d1.nonreal()+d2.nonreal())+d3.nonreal());
      d4=(d1+d2)+d3;
      TEST_ASSERT(d4.exact(d5));

      // test D=A+(B+C)
      d5.set_real(d1.real()+(d2.real()+d3.real()));
      d5.set_nonreal(d1.nonreal()+(d2.nonreal()+d3.nonreal()));
      d4=d1+(d2+d3);
      TEST_ASSERT(d4.exact(d5));

      // test D=A+B+C
      d5.set_real(d1.real()+d2.real()+d3.real());
      d5.set_nonreal(d1.nonreal()+d2.nonreal()+d3.nonreal());
      d4=d1+d2+d3;
      TEST_ASSERT(d4.exact(d5));

      // test D=(A+B)+(C+C)
      d5.set_real((d1.real()+d2.real())+(d3.real()+d3.real()));
      d5.set_nonreal((d1.nonreal()+d2.nonreal())+(d3.nonreal()+d3.nonreal()));
      d4=(d1+d2)+(d3+d3);
      TEST_ASSERT(d4.exact(d5));

      // test C=A+v
      d5.set_real(d2.real()+v1);
      d5.set_nonreal(d2.nonreal()+0.0);
      d4=d2+v1;
      TEST_ASSERT(d4.exact(d5));

      // test C=v+A
      d5.set_real(v1+d2.real());
      d5.set_nonreal(0.0+d2.nonreal());
      d4=v1+d2;
      TEST_ASSERT(d4.exact(d5));

      // test C=A+B+v
      d5.set_real(d2.real()+d3.real()+v1);
      d5.set_nonreal(d2.nonreal()+d3.nonreal()+0.0);
      d4=d2+d3+v1;
      TEST_ASSERT(d4.exact(d5));

      // test C=v+A+B
      d5.set_real(v1+d2.real()+d3.real());
      d5.set_nonreal(0.0+d2.nonreal()+d3.nonreal());
      d4=v1+d2+d3;
      TEST_ASSERT(d4.exact(d5));

      // test C=A+int
      d5.set_real(d2.real()+2);
      d5.set_nonreal(d2.nonreal()+0.0);
      d4=d2+2;
      TEST_ASSERT(d4.exact(d5));

      // test C=int+A
      d5.set_real(2+d2.real());
      d5.set_nonreal(d2.nonreal()+0.0);
      d4=2+d2;
      TEST_ASSERT(d4.exact(d5));
    }

    void dual_subtraction_test()
    {
      eli::ad::dual<data__, false> d1(3,-4), d2(1,2), d3(d2), d4, d5;
      data__ v1(2);

      // test C=A-B
      d5.set_real(d2.real()-d3.real());
      d5.set_nonreal(d2.nonreal()-d3.nonreal());
      d4=d2-d3;
      TEST_ASSERT(d4.exact(d5));

      // test B=A-B
      d4=d3;
      d5.set_real(d2.real()-d4.real());
      d5.set_nonreal(d2.nonreal()-d4.nonreal());
      d4=d2-d4;
      TEST_ASSERT(d4.exact(d5));

      // test B-=A
      d4=d3;
      d5.set_real(d2.real()-d4.real());
      d5.set_nonreal(d2.nonreal()-d4.nonreal());
      d4-=d2;
      TEST_ASSERT(d4.exact(d5));

      // test B-=v
      d4=d3;
      d5.set_real(d2.real()-v1);
      d5.set_nonreal(d2.nonreal()-0.0);
      d4-=v1;
      TEST_ASSERT(d4.exact(d5));

       // test B-=int
      d4=d3;
      d5.set_real(d2.real()-2);
      d5.set_nonreal(d2.nonreal()-0.0);
      d4-=2;
      TEST_ASSERT(d4.exact(d5));

      // test C-=A-B
      d4=d3;
      d5.set_real(d4.real()-(d1.real()-d2.real()));
      d5.set_nonreal(d4.nonreal()-(d1.nonreal()-d2.nonreal()));
      d4-=d1-d2;
      TEST_ASSERT(d4.exact(d5));

      // test D=(A-B)-C
      d5.set_real((d1.real()-d2.real())-d3.real());
      d5.set_nonreal((d1.nonreal()-d2.nonreal())-d3.nonreal());
      d4=(d1-d2)-d3;
      TEST_ASSERT(d4.exact(d5));

      // test D=A-(B-C)
      d5.set_real(d1.real()-(d2.real()-d3.real()));
      d5.set_nonreal(d1.nonreal()-(d2.nonreal()-d3.nonreal()));
      d4=d1-(d2-d3);
      TEST_ASSERT(d4.exact(d5));

      // test D=A-B-C
      d5.set_real(d1.real()-d2.real()-d3.real());
      d5.set_nonreal(d1.nonreal()-d2.nonreal()-d3.nonreal());
      d4=d1-d2-d3;
      TEST_ASSERT(d4.exact(d5));

      // test D=(A-B)-(C+C)
      d5.set_real((d1.real()-d2.real())-(d3.real()+d3.real()));
      d5.set_nonreal((d1.nonreal()-d2.nonreal())-(d3.nonreal()+d3.nonreal()));
      d4=(d1-d2)-(d3+d3);
      TEST_ASSERT(d4.exact(d5));

      // test C=A-v
      d5.set_real(d2.real()-v1);
      d5.set_nonreal(d2.nonreal()-0.0);
      d4=d2-v1;
      TEST_ASSERT(d4.exact(d5));

      // test C=v-A
      d5.set_real(v1-d2.real());
      d5.set_nonreal(0.0-d2.nonreal());
      d4=v1-d2;
      TEST_ASSERT(d4.exact(d5));

      // test C=A-B-v
      d5.set_real(d2.real()-d3.real()-v1);
      d5.set_nonreal(d2.nonreal()-d3.nonreal()-0.0);
      d4=d2-d3-v1;
      TEST_ASSERT(d4.exact(d5));

      // test C=v-A-B
      d5.set_real(v1-d2.real()-d3.real());
      d5.set_nonreal(0.0-d2.nonreal()-d3.nonreal());
      d4=v1-d2-d3;
      TEST_ASSERT(d4.exact(d5));

      // test C=A-int
      d5.set_real(d2.real()-2);
      d5.set_nonreal(d2.nonreal()-0.0);
      d4=d2-2;
      TEST_ASSERT(d4.exact(d5));

      // test C=int-A
      d5.set_real(2-d2.real());
      d5.set_nonreal(0.0-d2.nonreal());
      d4=2-d2;
      TEST_ASSERT(d4.exact(d5));
    }

    void dual_multiplication_test()
    {
      eli::ad::dual<data__, false> d1(3,-4), d2(1,2), d3(d2), d4, d5, d6;
      data__ v1(2);

      // test C=A*B
      d5.set_real(d2.real()*d3.real());
      d5.set_nonreal(d2.real()*d3.nonreal()+d2.nonreal()*d3.real());
      d4=d2*d3;
      TEST_ASSERT(d4.exact(d5));

      // test B=A*B
      d4=d3;
      d5.set_real(d2.real()*d4.real());
      d5.set_nonreal(d2.real()*d4.nonreal()+d2.nonreal()*d4.real());
      d4=d2*d4;
      TEST_ASSERT(d4.exact(d5));

      // test B*=A
      d4=d3;
      d5.set_real(d2.real()*d4.real());
      d5.set_nonreal(d2.real()*d4.nonreal()+d2.nonreal()*d4.real());
      d4*=d2;
      TEST_ASSERT(d4.exact(d5));

      // test B*=v
      d4=d3;
      d5.set_real(d2.real()*v1);
      d5.set_nonreal(d2.nonreal()*v1);
      d4*=v1;
      TEST_ASSERT(d4.exact(d5));

       // test B*=int
      d4=d3;
      d5.set_real(d2.real()*2);
      d5.set_nonreal(d2.nonreal()*2);
      d4*=2;
      TEST_ASSERT(d4.exact(d5));

      // test C*=A-B
      d4=d3;
      d6=d1-d2;
      d5.set_real(d4.real()*d6.real());
      d5.set_nonreal(d4.real()*d6.nonreal()+d4.nonreal()*d6.real());
      d4*=d1-d2;
      TEST_ASSERT(d4.exact(d5));

      // test D=(A-B)*C
      d6=d1-d2;
      d5.set_real(d6.real()*d3.real());
      d5.set_nonreal(d6.real()*d3.nonreal()+d6.nonreal()*d3.real());
      d4=(d1-d2)*d3;
      TEST_ASSERT(d4.exact(d5));

      // test D=A*(B-C)
      d6=d2-d3;
      d5.set_real(d1.real()*d6.real());
      d5.set_nonreal(d1.real()*d6.nonreal()+d1.nonreal()*d6.real());
      d4=d1*(d2-d3);
      TEST_ASSERT(d4.exact(d5));

      // test D=A*B*C
      d6.set_real(d1.real()*d2.real());
      d6.set_nonreal(d1.real()*d2.nonreal()+d1.nonreal()*d2.real());
      d5.set_real(d6.real()*d3.real());
      d5.set_nonreal(d6.real()*d3.nonreal()+d6.nonreal()*d3.real());
      d4=d1*d2*d3;
      TEST_ASSERT(d4.exact(d5));

      // test D=(A*B)-(C*C)
      d5.set_real(d1.real()*d2.real());
      d5.set_nonreal(d1.real()*d2.nonreal()+d1.nonreal()*d2.real());
      d6.set_real(d3.real()*d3.real());
      d6.set_nonreal(d3.real()*d3.nonreal()+d3.nonreal()*d3.real());
      d5-=d6;
      d4=(d1*d2)-(d3*d3);
      TEST_ASSERT(d4.exact(d5));

      // test C=A*v
      d5.set_real(d2.real()*v1);
      d5.set_nonreal(d2.nonreal()*v1);
      d4=d2*v1;
      TEST_ASSERT(d4.exact(d5));

      // test C=v*A
      d5.set_real(v1*d2.real());
      d5.set_nonreal(v1*d2.nonreal());
      d4=v1*d2;
      TEST_ASSERT(d4.exact(d5));

      // test C=A*B*v
      d5.set_real(d2.real()*d3.real()*v1);
      d5.set_nonreal((d2.real()*d3.nonreal()+d2.nonreal()*d3.real())*v1);
      d4=d2*d3*v1;
      TEST_ASSERT(d4.exact(d5));

      // test C=v*A*B
      d5.set_real(v1*d2.real()*d3.real());
      d5.set_nonreal(v1*(d2.real()*d3.nonreal()+d2.nonreal()*d3.real()));
      d4=v1*d2*d3;
      TEST_ASSERT(d4.exact(d5));

      // test C=A*int
      d5.set_real(d2.real()*2);
      d5.set_nonreal(d2.nonreal()*2);
      d4=d2*2;
      TEST_ASSERT(d4.exact(d5));

      // test C=int*A
      d5.set_real(2*d2.real());
      d5.set_nonreal(2*d2.nonreal());
      d4=2*d2;
      TEST_ASSERT(d4.exact(d5));
    }

    void dual_division_test()
    {
      eli::ad::dual<data__, false> d1(3,-4), d2(1,2), d3(d2), d4, d5, d6, d7;
      data__ v1(2);

      // test C=A/B
      d5.set_real(d2.real()/d3.real());
      d5.set_nonreal((d2.nonreal()*d3.real()-d2.real()*d3.nonreal())/d3.real()/d3.real());
      d4=d2/d3;
      TEST_ASSERT(d4.exact(d5));

      // test B=A/B
      d4=d3;
      d5.set_real(d2.real()/d4.real());
      d5.set_nonreal((d2.nonreal()*d4.real()-d2.real()*d4.nonreal())/d4.real()/d4.real());
      d4=d2/d4;
      TEST_ASSERT(d4.exact(d5));

      // test B/=A
      d4=d3;
      d5.set_real(d2.real()/d4.real());
      d5.set_nonreal((d2.nonreal()/d4.real()-d2.real()*d4.nonreal())/d4.real()/d4.real());
      d4/=d2;
      TEST_ASSERT(d4.exact(d5));

      // test B/=v
      d4=d3;
      d5.set_real(d2.real()/v1);
      d5.set_nonreal(d2.nonreal()/v1);
      d4/=v1;
      TEST_ASSERT(d4.exact(d5));

       // test B/=int
      d4=d3;
      d5.set_real(d2.real()/2);
      d5.set_nonreal(d2.nonreal()/2);
      d4/=2;
      TEST_ASSERT(d4.exact(d5));

      // test C/=A-B
      d4=d3;
      d6=d1-d2;
      d5.set_real(d4.real()/d6.real());
      d5.set_nonreal((d4.nonreal()*d6.real()-d4.real()*d6.nonreal())/d6.real()/d6.real());
      d4/=d1-d2;
      TEST_ASSERT(d4.exact(d5));

      // test D=(A-B)/C
      d6=d1-d2;
      d5.set_real(d6.real()/d3.real());
      d5.set_nonreal((d6.nonreal()*d3.real()-d6.real()*d3.nonreal())/d3.real()/d3.real());
      d4=(d1-d2)/d3;
      TEST_ASSERT(d4.exact(d5));

      // test D=A/(B-C)
      d6=d1-d2;
      d5.set_real(d1.real()/d6.real());
      d5.set_nonreal((d1.nonreal()*d6.real()-d1.real()*d6.nonreal())/d6.real()/d6.real());
      d4=d1/(d1-d2);
      TEST_ASSERT(d4.exact(d5));

      // test D=A*B*C
      d6.set_real(d1.real()/d2.real());
      d6.set_nonreal((d1.nonreal()*d2.real()-d1.real()*d2.nonreal())/d2.real()/d2.real());
      d5.set_real(d6.real()/d3.real());
      d5.set_nonreal((d6.nonreal()*d3.real()-d6.real()*d3.nonreal())/d3.real()/d3.real());
      d4=d1/d2/d3;
      TEST_ASSERT(d4.exact(d5));

      // test E=(A/B)-(C/D)
      d5.set_real(d1.real()/d2.real());
      d5.set_nonreal((d1.nonreal()*d2.real()-d1.real()*d2.nonreal())/d2.real()/d2.real());
      TEST_ASSERT(d5.exact(d1/d2));
      d6.set_real(d3.real()/d4.real());
      d6.set_nonreal((d3.nonreal()*d4.real()-d3.real()*d4.nonreal())/d4.real()/d4.real());
      TEST_ASSERT(d6.exact(d3/d4));
      d5-=d6;
      d7=(d1/d2)-(d3/d4);
      TEST_ASSERT(d7.exact(d5));

      // test C=A/v
      d5.set_real(d2.real()/v1);
      d5.set_nonreal(d2.nonreal()/v1);
      d4=d2/v1;
      TEST_ASSERT(d4.exact(d5));

      // test C=v/A
      d5.set_real(v1/d2.real());
      d5.set_nonreal(-v1*d2.nonreal()/d2.real()/d2.real());
      d4=v1/d2;
      TEST_ASSERT(d4.exact(d5));

      // test C=A/B/v
      d5.set_real(d2.real()/d3.real()/v1);
      d5.set_nonreal((d2.nonreal()*d3.real()-d2.real()*d3.nonreal())/d3.real()/d3.real()/v1);
      d4=d2/d3/v1;
      TEST_ASSERT(d4.exact(d5));

      // test C=v/A/B
      d6.set_real(v1/d2.real());
      d6.set_nonreal(-v1*d2.nonreal()/d2.real()/d2.real());
      d5.set_real(d6.real()/d3.real());
      d5.set_nonreal((d6.nonreal()*d3.real()-d6.real()*d3.nonreal())/d3.real()/d3.real());
      d4=v1/d2/d3;
      TEST_ASSERT(d4.exact(d5));

      // test C=A/int
      d5.set_real(d2.real()/2);
      d5.set_nonreal(d2.nonreal()/2);
      d4=d2/2;
      TEST_ASSERT(d4.exact(d5));

      // test C=int/A
      d5.set_real(2/d2.real());
      d5.set_nonreal(-2*d2.nonreal()/d2.real()/d2.real());
      d4=2/d2;
      TEST_ASSERT(d4.exact(d5));
    }

    void dual_trig_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ v1(4);
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // sin(dual) test
      d3=d1;
      dref.set_real(std::sin(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::cos(d3.real()));
      d0=std::sin(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sin(expression) test
      d3=d1+d2;
      dref.set_real(std::sin(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::cos(d3.real()));
      d0=std::sin(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sin(sin(expression)) test
      d3=d1+d2;
      d4.set_real(std::sin(d3.real()));
      d4.set_nonreal(d3.nonreal()*std::cos(d3.real()));
      dref.set_real(std::sin(d4.real()));
      dref.set_nonreal(d4.nonreal()*std::cos(d4.real()));
      d0=std::sin(std::sin(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // cos(dual) test
      d3=d1;
      dref.set_real(std::cos(d3.real()));
      dref.set_nonreal(-d3.nonreal()*std::sin(d3.real()));
      d0=std::cos(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cos(expression) test
      d3=d1+d2;
      dref.set_real(std::cos(d3.real()));
      dref.set_nonreal(-d3.nonreal()*std::sin(d3.real()));
      d0=std::cos(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cos(cos(expression)) test
      d3=d1+d2;
      d4.set_real(std::cos(d3.real()));
      d4.set_nonreal(-d3.nonreal()*std::sin(d3.real()));
      dref.set_real(std::cos(d4.real()));
      dref.set_nonreal(-d4.nonreal()*std::sin(d4.real()));
      d0=std::cos(std::cos(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // tan(dual) test
      d3=d1;
      dref.set_real(std::tan(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::cos(d3.real())/std::cos(d3.real()));
      d0=std::tan(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // tan(expression) test
      d3=d1+d2;
      dref.set_real(std::tan(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::cos(d3.real())/std::cos(d3.real()));
      d0=std::tan(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // tan(tan(expression)) test
      d3=d1+d2;
      d4.set_real(std::tan(d3.real()));
      d4.set_nonreal(d3.nonreal()/std::cos(d3.real())/std::cos(d3.real()));
      dref.set_real(std::tan(d4.real()));
      dref.set_nonreal(d4.nonreal()/std::cos(d4.real())/std::cos(d4.real()));
      d0=std::tan(std::tan(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // asin(dual) test
      d3=d1/v1;
      dref.set_real(std::asin(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::sqrt(1-d3.real()*d3.real()));
      d0=std::asin(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // asin(expression) test
      d3=(d1+d2)/v1;
      dref.set_real(std::asin(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::sqrt(1-d3.real()*d3.real()));
      d0=std::asin((d1+d2)/v1);
      TEST_ASSERT(d0.nearly(dref, tol));

      // asin(sin(expression)) test
      d3=(d1+d2)/v1;
      dref=d3;
      d0=std::asin(std::sin((d1+d2)/v1));
      TEST_ASSERT(d0.nearly(dref, tol));

      // acos(dual) test
      d3=d1/v1;
      dref.set_real(std::acos(d3.real()));
      dref.set_nonreal(-d3.nonreal()/std::sqrt(1-d3.real()*d3.real()));
      d0=std::acos(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // acos(expression) test
      d3=(d1+d2)/v1;
      dref.set_real(std::acos(d3.real()));
      dref.set_nonreal(-d3.nonreal()/std::sqrt(1-d3.real()*d3.real()));
      d0=std::acos((d1+d2)/v1);
      TEST_ASSERT(d0.nearly(dref, tol));

      // acos(cos(expression)) test
      d3=(d1+d2)/v1;
      dref=d3;
      d0=std::acos(std::cos((d1+d2)/v1));
      TEST_ASSERT(d0.nearly(dref, tol));

      // atan(dual) test
      d3=d1/v1;
      dref.set_real(std::atan(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1+d3.real()*d3.real()));
      d0=std::atan(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // atan(expression) test
      d3=(d1+d2)/v1;
      dref.set_real(std::atan(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1+d3.real()*d3.real()));
      d0=std::atan((d1+d2)/v1);
      TEST_ASSERT(d0.nearly(dref, tol));

      // atan(tan(expression)) test
      d3=(d1+d2)/v1;
      dref=d3;
      d0=std::atan(std::tan((d1+d2)/v1));
      TEST_ASSERT(d0.nearly(dref, tol));

      // TODO: ADD atan2 test when bin_fun implemented
    }

    void dual_trig_nostd_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ v1(4);
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // sin(dual) test
      d3=d1;
      dref.set_real(sin(d3.real()));
      dref.set_nonreal(d3.nonreal()*cos(d3.real()));
      d0=sin(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sin(expression) test
      d3=d1+d2;
      dref.set_real(sin(d3.real()));
      dref.set_nonreal(d3.nonreal()*cos(d3.real()));
      d0=sin(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sin(sin(expression)) test
      d3=d1+d2;
      d4.set_real(sin(d3.real()));
      d4.set_nonreal(d3.nonreal()*cos(d3.real()));
      dref.set_real(sin(d4.real()));
      dref.set_nonreal(d4.nonreal()*cos(d4.real()));
      d0=sin(sin(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // cos(dual) test
      d3=d1;
      dref.set_real(cos(d3.real()));
      dref.set_nonreal(-d3.nonreal()*sin(d3.real()));
      d0=cos(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cos(expression) test
      d3=d1+d2;
      dref.set_real(cos(d3.real()));
      dref.set_nonreal(-d3.nonreal()*sin(d3.real()));
      d0=cos(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cos(cos(expression)) test
      d3=d1+d2;
      d4.set_real(cos(d3.real()));
      d4.set_nonreal(-d3.nonreal()*sin(d3.real()));
      dref.set_real(cos(d4.real()));
      dref.set_nonreal(-d4.nonreal()*sin(d4.real()));
      d0=cos(cos(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // tan(dual) test
      d3=d1;
      dref.set_real(tan(d3.real()));
      dref.set_nonreal(d3.nonreal()/cos(d3.real())/cos(d3.real()));
      d0=tan(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // tan(expression) test
      d3=d1+d2;
      dref.set_real(tan(d3.real()));
      dref.set_nonreal(d3.nonreal()/cos(d3.real())/cos(d3.real()));
      d0=tan(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // tan(tan(expression)) test
      d3=d1+d2;
      d4.set_real(tan(d3.real()));
      d4.set_nonreal(d3.nonreal()/cos(d3.real())/cos(d3.real()));
      dref.set_real(tan(d4.real()));
      dref.set_nonreal(d4.nonreal()/cos(d4.real())/cos(d4.real()));
      d0=tan(tan(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // asin(dual) test
      d3=d1/v1;
      dref.set_real(asin(d3.real()));
      dref.set_nonreal(d3.nonreal()/sqrt(1-d3.real()*d3.real()));
      d0=asin(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // asin(expression) test
      d3=(d1+d2)/v1;
      dref.set_real(asin(d3.real()));
      dref.set_nonreal(d3.nonreal()/sqrt(1-d3.real()*d3.real()));
      d0=asin((d1+d2)/v1);
      TEST_ASSERT(d0.nearly(dref, tol));

      // asin(sin(expression)) test
      d3=(d1+d2)/v1;
      dref=d3;
      d0=asin(sin((d1+d2)/v1));
      TEST_ASSERT(d0.nearly(dref, tol));

      // acos(dual) test
      d3=d1/v1;
      dref.set_real(acos(d3.real()));
      dref.set_nonreal(-d3.nonreal()/sqrt(1-d3.real()*d3.real()));
      d0=acos(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // acos(expression) test
      d3=(d1+d2)/v1;
      dref.set_real(acos(d3.real()));
      dref.set_nonreal(-d3.nonreal()/sqrt(1-d3.real()*d3.real()));
      d0=acos((d1+d2)/v1);
      TEST_ASSERT(d0.nearly(dref, tol));

      // acos(cos(expression)) test
      d3=(d1+d2)/v1;
      dref=d3;
      d0=acos(cos((d1+d2)/v1));
      TEST_ASSERT(d0.nearly(dref, tol));

      // atan(dual) test
      d3=d1/v1;
      dref.set_real(atan(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1+d3.real()*d3.real()));
      d0=atan(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // atan(expression) test
      d3=(d1+d2)/v1;
      dref.set_real(atan(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1+d3.real()*d3.real()));
      d0=atan((d1+d2)/v1);
      TEST_ASSERT(d0.nearly(dref, tol));

      // tan(tan(expression)) test
      d3=(d1+d2)/v1;
      dref=d3;
      d0=atan(std::tan((d1+d2)/v1));
      TEST_ASSERT(d0.nearly(dref, tol));

      // TODO: ADD atan2 test when bin_fun implemented
    }

    void dual_hyp_trig_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // sinh(dual) test
      d3=d1;
      dref.set_real(std::sinh(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::cosh(d3.real()));
      d0=std::sinh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sinh(expression) test
      d3=d1+d2;
      dref.set_real(std::sinh(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::cosh(d3.real()));
      d0=std::sinh(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sinh(sinh(expression)) test
      d3=d1+d2;
      d4.set_real(std::sinh(d3.real()));
      d4.set_nonreal(d3.nonreal()*std::cosh(d3.real()));
      dref.set_real(std::sinh(d4.real()));
      dref.set_nonreal(d4.nonreal()*std::cosh(d4.real()));
      d0=std::sinh(std::sinh(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // cosh(dual) test
      d3=d1;
      dref.set_real(std::cosh(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::sinh(d3.real()));
      d0=std::cosh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cosh(expression) test
      d3=d1+d2;
      dref.set_real(std::cosh(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::sinh(d3.real()));
      d0=std::cosh(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cosh(cosh(expression)) test
      d3=d1+d2;
      d4.set_real(std::cosh(d3.real()));
      d4.set_nonreal(d3.nonreal()*std::sinh(d3.real()));
      dref.set_real(std::cosh(d4.real()));
      dref.set_nonreal(d4.nonreal()*std::sinh(d4.real()));
      d0=std::cosh(std::cosh(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // tanh(dual) test
      d3=d1;
      dref.set_real(std::tanh(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::cosh(d3.real())/std::cosh(d3.real()));
      d0=std::tanh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // tanh(expression) test
      d3=d1+d2;
      dref.set_real(std::tanh(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::cosh(d3.real())/std::cosh(d3.real()));
      d0=std::tanh(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // tanh(tanh(expression)) test
      d3=d1+d2;
      d4.set_real(std::tanh(d3.real()));
      d4.set_nonreal(d3.nonreal()/std::cosh(d3.real())/std::cosh(d3.real()));
      dref.set_real(std::tanh(d4.real()));
      dref.set_nonreal(d4.nonreal()/std::cosh(d4.real())/std::cosh(d4.real()));
      d0=std::tanh(std::tanh(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // asinh(dual) test
      d3=d1;
      dref.set_real(std::asinh(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::sqrt(1+d3.real()*d3.real()));
      d0=std::asinh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // asinh(expression) test
      d3=(d1+d2);
      dref.set_real(std::asinh(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::sqrt(1+d3.real()*d3.real()));
      d0=std::asinh((d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // asinh(sinh(expression)) test
      d3=(d1+d2);
      dref=d3;
      d0=std::asinh(std::sinh(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // acosh(dual) test
      d3=d1;
      dref.set_real(std::acosh(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::sqrt(d3.real()*d3.real()-1));
      d0=std::acosh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // acosh(expression) test
      d3=(d1+d2);
      dref.set_real(std::acosh(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::sqrt(d3.real()*d3.real()-1));
      d0=std::acosh(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // acosh(cosh(expression)) test
      d3=(d1+d2);
      dref=d3;
      d0=std::acos(std::cosh(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // atanh(dual) test
      d3=d1/4;
      dref.set_real(std::atanh(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1-d3.real()*d3.real()));
      d0=std::atanh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // atanh(expression) test
      d3=(d1+d2)/4;
      dref.set_real(std::atanh(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1-d3.real()*d3.real()));
      d0=std::atanh((d1+d2)/4);
      TEST_ASSERT(d0.nearly(dref, tol));

      // atanh(tanh(expression)) test
      d3=(d1+d2)/4;
      dref=d3;
      d0=std::atanh(std::tanh((d1+d2)/4));
      TEST_ASSERT(d0.nearly(dref, tol));
    }

    void dual_hyp_trig_nostd_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // sinh(dual) test
      d3=d1;
      dref.set_real(sinh(d3.real()));
      dref.set_nonreal(d3.nonreal()*cosh(d3.real()));
      d0=sinh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sinh(expression) test
      d3=d1+d2;
      dref.set_real(sinh(d3.real()));
      dref.set_nonreal(d3.nonreal()*cosh(d3.real()));
      d0=sinh(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sinh(sinh(expression)) test
      d3=d1+d2;
      d4.set_real(sinh(d3.real()));
      d4.set_nonreal(d3.nonreal()*cosh(d3.real()));
      dref.set_real(sinh(d4.real()));
      dref.set_nonreal(d4.nonreal()*cosh(d4.real()));
      d0=sinh(sinh(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // cosh(dual) test
      d3=d1;
      dref.set_real(cosh(d3.real()));
      dref.set_nonreal(d3.nonreal()*sinh(d3.real()));
      d0=cosh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cosh(expression) test
      d3=d1+d2;
      dref.set_real(cosh(d3.real()));
      dref.set_nonreal(d3.nonreal()*sinh(d3.real()));
      d0=cosh(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cosh(cosh(expression)) test
      d3=d1+d2;
      d4.set_real(cosh(d3.real()));
      d4.set_nonreal(d3.nonreal()*sinh(d3.real()));
      dref.set_real(cosh(d4.real()));
      dref.set_nonreal(d4.nonreal()*sinh(d4.real()));
      d0=cosh(cosh(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // tanh(dual) test
      d3=d1;
      dref.set_real(tanh(d3.real()));
      dref.set_nonreal(d3.nonreal()/cosh(d3.real())/cosh(d3.real()));
      d0=tanh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // tanh(expression) test
      d3=d1+d2;
      dref.set_real(tanh(d3.real()));
      dref.set_nonreal(d3.nonreal()/cosh(d3.real())/cosh(d3.real()));
      d0=tanh(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // tanh(tanh(expression)) test
      d3=d1+d2;
      d4.set_real(tanh(d3.real()));
      d4.set_nonreal(d3.nonreal()/cosh(d3.real())/cosh(d3.real()));
      dref.set_real(tanh(d4.real()));
      dref.set_nonreal(d4.nonreal()/cosh(d4.real())/cosh(d4.real()));
      d0=tanh(tanh(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // asinh(dual) test
      d3=d1;
      dref.set_real(asinh(d3.real()));
      dref.set_nonreal(d3.nonreal()/sqrt(1+d3.real()*d3.real()));
      d0=asinh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // asinh(expression) test
      d3=(d1+d2);
      dref.set_real(asinh(d3.real()));
      dref.set_nonreal(d3.nonreal()/sqrt(1+d3.real()*d3.real()));
      d0=asinh((d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // asinh(sinh(expression)) test
      d3=(d1+d2);
      dref=d3;
      d0=asinh(sinh(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // acosh(dual) test
      d3=d1;
      dref.set_real(acosh(d3.real()));
      dref.set_nonreal(d3.nonreal()/sqrt(d3.real()*d3.real()-1));
      d0=acosh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // acosh(expression) test
      d3=(d1+d2);
      dref.set_real(acosh(d3.real()));
      dref.set_nonreal(d3.nonreal()/sqrt(d3.real()*d3.real()-1));
      d0=acosh(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // acosh(cosh(expression)) test
      d3=(d1+d2);
      dref=d3;
      d0=acos(cosh(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // atanh(dual) test
      d3=d1/4;
      dref.set_real(atanh(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1-d3.real()*d3.real()));
      d0=atanh(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // atanh(expression) test
      d3=(d1+d2)/4;
      dref.set_real(atanh(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1-d3.real()*d3.real()));
      d0=atanh((d1+d2)/4);
      TEST_ASSERT(d0.nearly(dref, tol));

      // atanh(tanh(expression)) test
      d3=(d1+d2)/4;
      dref=d3;
      d0=atanh(tanh((d1+d2)/4));
      TEST_ASSERT(d0.nearly(dref, tol));
    }

    void dual_exp_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // exp(dual) test
      d3=d1;
      dref.set_real(std::exp(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::exp(d3.real()));
      d0=std::exp(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // exp(expression) test
      d3=d1+d2;
      dref.set_real(std::exp(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::exp(d3.real()));
      d0=std::exp(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // exp(exp(expression)) test
      d3=d1+d2;
      d4.set_real(std::exp(d3.real()));
      d4.set_nonreal(d3.nonreal()*std::exp(d3.real()));
      dref.set_real(std::exp(d4.real()));
      dref.set_nonreal(d4.nonreal()*std::exp(d4.real()));
      d0=std::exp(std::exp(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // expm1(dual) test
      d3=d1;
      dref.set_real(std::expm1(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::exp(d3.real()));
      d0=std::expm1(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // expm1(expression) test
      d3=d1+d2;
      dref.set_real(std::expm1(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::exp(d3.real()));
      d0=std::expm1(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // expm1(expm1(expression)) test
      d3=d1+d2;
      d4.set_real(std::expm1(d3.real()));
      d4.set_nonreal(d3.nonreal()*std::exp(d3.real()));
      dref.set_real(std::expm1(d4.real()));
      dref.set_nonreal(d4.nonreal()*std::exp(d4.real()));
      d0=std::expm1(std::expm1(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // exp2(dual) test
      d3=d1;
      dref.set_real(std::exp2(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::log(2)*std::exp2(d3.real()));
      d0=std::exp2(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // exp2(expression) test
      d3=d1+d2;
      dref.set_real(std::exp2(d3.real()));
      dref.set_nonreal(d3.nonreal()*std::log(2)*std::exp2(d3.real()));
      d0=std::exp2(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // exp2(exp2(expression)) test
      d3=d1+d2;
      d4.set_real(std::exp2(d3.real()));
      d4.set_nonreal(d3.nonreal()*std::log(2)*std::exp2(d3.real()));
      dref.set_real(std::exp2(d4.real()));
      dref.set_nonreal(d4.nonreal()*std::log(2)*std::exp2(d4.real()));
      d0=std::exp2(std::exp2(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // log(dual) test
      d3=d1;
      dref.set_real(std::log(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real());
      d0=std::log(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log(expression) test
      d3=d1+d2;
      dref.set_real(std::log(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real());
      d0=std::log(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log(log(expression)) test
      d3=d1+d2;
      d4.set_real(std::log(d3.real()));
      d4.set_nonreal(d3.nonreal()/d3.real());
      dref.set_real(std::log(d4.real()));
      dref.set_nonreal(d4.nonreal()/d4.real());
      d0=std::log(std::log(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // log10(dual) test
      d3=d1;
      dref.set_real(std::log10(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real()/std::log(10));
      d0=std::log10(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log10(expression) test
      d3=d1+d2;
      dref.set_real(std::log10(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real()/std::log(10));
      d0=std::log10(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log10(log10(expression)) test
      d3=d1+d2;
      d4.set_real(std::log10(d3.real()));
      d4.set_nonreal(d3.nonreal()/d3.real()/std::log(10));
      dref.set_real(std::log10(d4.real()));
      dref.set_nonreal(d4.nonreal()/d4.real()/std::log(10));
      d0=std::log10(std::log10(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // log2(dual) test
      d3=d1;
      dref.set_real(std::log2(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real()/std::log(2));
      d0=std::log2(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log2(expression) test
      d3=d1+d2;
      dref.set_real(std::log2(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real()/std::log(2));
      d0=std::log2(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log2(log2(expression)) test
      d3=d1+d2;
      d4.set_real(std::log2(d3.real()));
      d4.set_nonreal(d3.nonreal()/d3.real()/std::log(2));
      dref.set_real(std::log2(d4.real()));
      dref.set_nonreal(d4.nonreal()/d4.real()/std::log(2));
      d0=std::log2(std::log2(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // log1p(dual) test
      d3=d1;
      dref.set_real(std::log1p(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1+d3.real()));
      d0=std::log1p(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log1p(expression) test
      d3=d1+d2;
      dref.set_real(std::log1p(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1+d3.real()));
      d0=std::log1p(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log1p(log1p(expression)) test
      d3=d1+d2;
      d4.set_real(std::log1p(d3.real()));
      d4.set_nonreal(d3.nonreal()/(1+d3.real()));
      dref.set_real(std::log1p(d4.real()));
      dref.set_nonreal(d4.nonreal()/(1+d4.real()));
      d0=std::log1p(std::log1p(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));
    }

    void dual_exp_nostd_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // exp(dual) test
      d3=d1;
      dref.set_real(exp(d3.real()));
      dref.set_nonreal(d3.nonreal()*exp(d3.real()));
      d0=exp(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // exp(expression) test
      d3=d1+d2;
      dref.set_real(exp(d3.real()));
      dref.set_nonreal(d3.nonreal()*exp(d3.real()));
      d0=exp(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // exp(exp(expression)) test
      d3=d1+d2;
      d4.set_real(exp(d3.real()));
      d4.set_nonreal(d3.nonreal()*exp(d3.real()));
      dref.set_real(exp(d4.real()));
      dref.set_nonreal(d4.nonreal()*exp(d4.real()));
      d0=exp(exp(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // expm1(dual) test
      d3=d1;
      dref.set_real(expm1(d3.real()));
      dref.set_nonreal(d3.nonreal()*exp(d3.real()));
      d0=expm1(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // expm1(expression) test
      d3=d1+d2;
      dref.set_real(expm1(d3.real()));
      dref.set_nonreal(d3.nonreal()*exp(d3.real()));
      d0=expm1(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // expm1(expm1(expression)) test
      d3=d1+d2;
      d4.set_real(expm1(d3.real()));
      d4.set_nonreal(d3.nonreal()*exp(d3.real()));
      dref.set_real(expm1(d4.real()));
      dref.set_nonreal(d4.nonreal()*exp(d4.real()));
      d0=expm1(expm1(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // exp2(dual) test
      d3=d1;
      dref.set_real(exp2(d3.real()));
      dref.set_nonreal(d3.nonreal()*log(2)*exp2(d3.real()));
      d0=exp2(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // exp2(expression) test
      d3=d1+d2;
      dref.set_real(exp2(d3.real()));
      dref.set_nonreal(d3.nonreal()*log(2)*exp2(d3.real()));
      d0=exp2(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // exp2(exp2(expression)) test
      d3=d1+d2;
      d4.set_real(exp2(d3.real()));
      d4.set_nonreal(d3.nonreal()*log(2)*exp2(d3.real()));
      dref.set_real(exp2(d4.real()));
      dref.set_nonreal(d4.nonreal()*log(2)*exp2(d4.real()));
      d0=exp2(exp2(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // log(dual) test
      d3=d1;
      dref.set_real(log(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real());
      d0=log(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log(expression) test
      d3=d1+d2;
      dref.set_real(log(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real());
      d0=log(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log(log(expression)) test
      d3=d1+d2;
      d4.set_real(log(d3.real()));
      d4.set_nonreal(d3.nonreal()/d3.real());
      dref.set_real(log(d4.real()));
      dref.set_nonreal(d4.nonreal()/d4.real());
      d0=log(log(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // log10(dual) test
      d3=d1;
      dref.set_real(log10(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real()/log(10));
      d0=log10(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log10(expression) test
      d3=d1+d2;
      dref.set_real(log10(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real()/log(10));
      d0=log10(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log10(log10(expression)) test
      d3=d1+d2;
      d4.set_real(log10(d3.real()));
      d4.set_nonreal(d3.nonreal()/d3.real()/log(10));
      dref.set_real(log10(d4.real()));
      dref.set_nonreal(d4.nonreal()/d4.real()/log(10));
      d0=log10(log10(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // log2(dual) test
      d3=d1;
      dref.set_real(log2(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real()/log(2));
      d0=log2(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log2(expression) test
      d3=d1+d2;
      dref.set_real(log2(d3.real()));
      dref.set_nonreal(d3.nonreal()/d3.real()/log(2));
      d0=log2(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log2(log2(expression)) test
      d3=d1+d2;
      d4.set_real(log2(d3.real()));
      d4.set_nonreal(d3.nonreal()/d3.real()/log(2));
      dref.set_real(log2(d4.real()));
      dref.set_nonreal(d4.nonreal()/d4.real()/log(2));
      d0=log2(log2(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // log1p(dual) test
      d3=d1;
      dref.set_real(log1p(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1+d3.real()));
      d0=log1p(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log1p(expression) test
      d3=d1+d2;
      dref.set_real(log1p(d3.real()));
      dref.set_nonreal(d3.nonreal()/(1+d3.real()));
      d0=log1p(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // log1p(log1p(expression)) test
      d3=d1+d2;
      d4.set_real(log1p(d3.real()));
      d4.set_nonreal(d3.nonreal()/(1+d3.real()));
      dref.set_real(log1p(d4.real()));
      dref.set_nonreal(d4.nonreal()/(1+d4.real()));
      d0=log1p(log1p(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));
    }

    void dual_power_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // sqrt(dual) test
      d3=d1;
      dref.set_real(std::sqrt(d3.real()));
      dref.set_nonreal(0.5*d3.nonreal()/std::sqrt(d3.real()));
      d0=std::sqrt(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sqrt(expression) test
      d3=d1+d2;
      dref.set_real(std::sqrt(d3.real()));
      dref.set_nonreal(0.5*d3.nonreal()/std::sqrt(d3.real()));
      d0=std::sqrt(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sqrt(sqrt(expression)) test
      d3=d1+d2;
      d4.set_real(std::sqrt(d3.real()));
      d4.set_nonreal(0.5*d3.nonreal()/std::sqrt(d3.real()));
      dref.set_real(std::sqrt(d4.real()));
      dref.set_nonreal(0.5*d4.nonreal()/std::sqrt(d4.real()));
      d0=std::sqrt(std::sqrt(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // cbrt(dual) test
      d3=d1;
      dref.set_real(std::cbrt(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::cbrt(d3.real()*d3.real())/3);
      d0=std::cbrt(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cbrt(expression) test
      d3=d1+d2;
      dref.set_real(std::cbrt(d3.real()));
      dref.set_nonreal(d3.nonreal()/std::cbrt(d3.real()*d3.real())/3);
      d0=std::cbrt(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cbrt(cbrt(expression)) test
      d3=d1+d2;
      d4.set_real(std::cbrt(d3.real()));
      d4.set_nonreal(d3.nonreal()/std::cbrt(d3.real()*d3.real())/3);
      dref.set_real(std::cbrt(d4.real()));
      dref.set_nonreal(d4.nonreal()/std::cbrt(d4.real()*d4.real())/3);
      d0=std::cbrt(std::cbrt(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // TODO: ADD hypot test when bin_fun implemented

      // TODO: ADD pow test when bin_fun implemented
    }

    void dual_power_nostd_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // sqrt(dual) test
      d3=d1;
      dref.set_real(sqrt(d3.real()));
      dref.set_nonreal(0.5*d3.nonreal()/sqrt(d3.real()));
      d0=sqrt(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sqrt(expression) test
      d3=d1+d2;
      dref.set_real(sqrt(d3.real()));
      dref.set_nonreal(0.5*d3.nonreal()/sqrt(d3.real()));
      d0=sqrt(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // sqrt(sqrt(expression)) test
      d3=d1+d2;
      d4.set_real(sqrt(d3.real()));
      d4.set_nonreal(0.5*d3.nonreal()/sqrt(d3.real()));
      dref.set_real(sqrt(d4.real()));
      dref.set_nonreal(0.5*d4.nonreal()/sqrt(d4.real()));
      d0=sqrt(sqrt(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // cbrt(dual) test
      d3=d1;
      dref.set_real(cbrt(d3.real()));
      dref.set_nonreal(d3.nonreal()/cbrt(d3.real()*d3.real())/3);
      d0=cbrt(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cbrt(expression) test
      d3=d1+d2;
      dref.set_real(cbrt(d3.real()));
      dref.set_nonreal(d3.nonreal()/cbrt(d3.real()*d3.real())/3);
      d0=cbrt(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // cbrt(cbrt(expression)) test
      d3=d1+d2;
      d4.set_real(cbrt(d3.real()));
      d4.set_nonreal(d3.nonreal()/cbrt(d3.real()*d3.real())/3);
      dref.set_real(cbrt(d4.real()));
      dref.set_nonreal(d4.nonreal()/cbrt(d4.real()*d4.real())/3);
      d0=cbrt(cbrt(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // TODO: ADD hypot test when bin_fun implemented

      // TODO: ADD pow test when bin_fun implemented
    }

    void dual_misc_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // abs(dual) test
      d3=d1;
      dref.set_real(std::abs(d3.real()));
      dref.set_nonreal(d3.nonreal()*((d3.real()<0)?(-1):(1)));
      d0=std::abs(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // abs(expression) test
      d3=d1+d2;
      dref.set_real(std::abs(d3.real()));
      dref.set_nonreal(d3.nonreal()*((d3.real()<0)?(-1):(1)));
      d0=std::abs(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // abs(abs(expression)) test
      d3=d1+d2;
      d4.set_real(std::abs(d3.real()));
      d4.set_nonreal(d3.nonreal()*((d3.real()<0)?(-1):(1)));
      dref.set_real(std::abs(d4.real()));
      dref.set_nonreal(d4.nonreal()*((d4.real()<0)?(-1):(1)));
      d0=std::abs(std::abs(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // ceil(dual) test
      d3=d1;
      dref.set_real(std::ceil(d3.real()));
      dref.set_nonreal(d3.nonreal()*(0));
      d0=std::ceil(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // ceil(expression) test
      d3=d1+d2;
      dref.set_real(std::ceil(d3.real()));
      dref.set_nonreal(d3.nonreal()*(0));
      d0=std::ceil(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // ceil(ceil(expression)) test
      d3=d1+d2;
      d4.set_real(std::ceil(d3.real()));
      d4.set_nonreal(d3.nonreal()*(0));
      dref.set_real(std::ceil(d4.real()));
      dref.set_nonreal(d4.nonreal()*(0));
      d0=std::ceil(std::ceil(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // floor(dual) test
      d3=d1;
      dref.set_real(std::floor(d3.real()));
      dref.set_nonreal(d3.nonreal()*(0));
      d0=std::floor(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // floor(expression) test
      d3=d1+d2;
      dref.set_real(std::floor(d3.real()));
      dref.set_nonreal(d3.nonreal()*(0));
      d0=std::floor(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // floor(floor(expression)) test
      d3=d1+d2;
      d4.set_real(std::floor(d3.real()));
      d4.set_nonreal(d3.nonreal()*(0));
      dref.set_real(std::floor(d4.real()));
      dref.set_nonreal(d4.nonreal()*(0));
      d0=std::floor(std::floor(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // TODO: ADD fmod test when bin_fun implemented

      // TODO: ADD fmin test when bin_fun implemented

      // TODO: ADD fmax test when bin_fun implemented

      // TODO: ADD min test when bin_fun implemented

      // TODO: ADD max test when bin_fun implemented

      // TODO: ADD modf test when bin_fun implemented

      // TODO: ADD ldexp test when bin_fun implemented
    }

    void dual_misc_nostd_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // abs(dual) test
      d3=d1;
      dref.set_real(abs(d3.real()));
      dref.set_nonreal(d3.nonreal()*((d3.real()<0)?(-1):(1)));
      d0=abs(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // abs(expression) test
      d3=d1+d2;
      dref.set_real(abs(d3.real()));
      dref.set_nonreal(d3.nonreal()*((d3.real()<0)?(-1):(1)));
      d0=abs(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // abs(abs(expression)) test
      d3=d1+d2;
      d4.set_real(abs(d3.real()));
      d4.set_nonreal(d3.nonreal()*((d3.real()<0)?(-1):(1)));
      dref.set_real(abs(d4.real()));
      dref.set_nonreal(d4.nonreal()*((d4.real()<0)?(-1):(1)));
      d0=abs(abs(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // ceil(dual) test
      d3=d1;
      dref.set_real(ceil(d3.real()));
      dref.set_nonreal(d3.nonreal()*(0));
      d0=ceil(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // ceil(expression) test
      d3=d1+d2;
      dref.set_real(ceil(d3.real()));
      dref.set_nonreal(d3.nonreal()*(0));
      d0=ceil(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // ceil(ceil(expression)) test
      d3=d1+d2;
      d4.set_real(ceil(d3.real()));
      d4.set_nonreal(d3.nonreal()*(0));
      dref.set_real(ceil(d4.real()));
      dref.set_nonreal(d4.nonreal()*(0));
      d0=ceil(ceil(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // floor(dual) test
      d3=d1;
      dref.set_real(floor(d3.real()));
      dref.set_nonreal(d3.nonreal()*(0));
      d0=floor(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // floor(expression) test
      d3=d1+d2;
      dref.set_real(floor(d3.real()));
      dref.set_nonreal(d3.nonreal()*(0));
      d0=floor(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // floor(floor(expression)) test
      d3=d1+d2;
      d4.set_real(floor(d3.real()));
      d4.set_nonreal(d3.nonreal()*(0));
      dref.set_real(floor(d4.real()));
      dref.set_nonreal(d4.nonreal()*(0));
      d0=floor(floor(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // TODO: ADD fmod test when bin_fun implemented

      // TODO: ADD fmin test when bin_fun implemented

      // TODO: ADD fmax test when bin_fun implemented

      // TODO: ADD min test when bin_fun implemented

      // TODO: ADD max test when bin_fun implemented

      // TODO: ADD modf test when bin_fun implemented

      // TODO: ADD ldexp test when bin_fun implemented
    }

    void dual_erf_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // erf(dual) test
      d3=d1;
      dref.set_real(std::erf(d3.real()));
      dref.set_nonreal(d3.nonreal()*(2*std::exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=std::erf(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // erf(expression) test
      d3=d1+d2;
      dref.set_real(std::erf(d3.real()));
      dref.set_nonreal(d3.nonreal()*(2*std::exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=std::erf(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // erf(erf(expression)) test
      d3=d1+d2;
      d4.set_real(std::erf(d3.real()));
      d4.set_nonreal(d3.nonreal()*(2*std::exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      dref.set_real(std::erf(d4.real()));
      dref.set_nonreal(d4.nonreal()*(2*std::exp(-d4.real()*d4.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=std::erf(std::erf(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // erfc(dual) test
      d3=d1;
      dref.set_real(std::erfc(d3.real()));
      dref.set_nonreal(d3.nonreal()*(-2*std::exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=std::erfc(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // erfc(expression) test
      d3=d1+d2;
      dref.set_real(std::erfc(d3.real()));
      dref.set_nonreal(d3.nonreal()*(-2*std::exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=std::erfc(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // erfc(erfc(expression)) test
      d3=d1+d2;
      d4.set_real(std::erfc(d3.real()));
      d4.set_nonreal(d3.nonreal()*(-2*std::exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      dref.set_real(std::erfc(d4.real()));
      dref.set_nonreal(d4.nonreal()*(-2*std::exp(-d4.real()*d4.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=std::erfc(std::erfc(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));
    }

    void dual_erf_nostd_test()
    {
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // erf(dual) test
      d3=d1;
      dref.set_real(erf(d3.real()));
      dref.set_nonreal(d3.nonreal()*(2*exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=erf(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // erf(expression) test
      d3=d1+d2;
      dref.set_real(erf(d3.real()));
      dref.set_nonreal(d3.nonreal()*(2*exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=erf(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // erf(erf(expression)) test
      d3=d1+d2;
      d4.set_real(erf(d3.real()));
      d4.set_nonreal(d3.nonreal()*(2*exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      dref.set_real(erf(d4.real()));
      dref.set_nonreal(d4.nonreal()*(2*exp(-d4.real()*d4.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=erf(erf(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // erfc(dual) test
      d3=d1;
      dref.set_real(erfc(d3.real()));
      dref.set_nonreal(d3.nonreal()*(-2*exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=erfc(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // erfc(expression) test
      d3=d1+d2;
      dref.set_real(erfc(d3.real()));
      dref.set_nonreal(d3.nonreal()*(-2*exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=erfc(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // erfc(erfc(expression)) test
      d3=d1+d2;
      d4.set_real(erfc(d3.real()));
      d4.set_nonreal(d3.nonreal()*(-2*exp(-d3.real()*d3.real())/eli::constants::math<data__>::sqrt_pi()));
      dref.set_real(erfc(d4.real()));
      dref.set_nonreal(d4.nonreal()*(-2*exp(-d4.real()*d4.real())/eli::constants::math<data__>::sqrt_pi()));
      d0=erfc(erfc(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));
    }

      // error & gamma functions: lgamma, tgamma

    void dual_gamma_test()
    {
      // TODO: ADD lgamma test when know how to implement derivative

      // TODO: ADD tgamma test when know how to implement derivative
#if 0
      eli::ad::dual<data__, false> d1(2,-2), d2(2,3), d3, d4, d0, dref;
      data__ tol(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // lgamma(dual) test
      d3=d1;
      dref.set_real(std::lgamma(d3.real()));
      dref.set_nonreal(d3.nonreal()*(xxx));
      d0=std::lgamma(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // lgamma(expression) test
      d3=d1+d2;
      dref.set_real(std::lgamma(d3.real()));
      dref.set_nonreal(d3.nonreal()*(xxx));
      d0=std::lgamma(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // lgamma(lgamma(expression)) test
      d3=d1+d2;
      d4.set_real(std::lgamma(d3.real()));
      d4.set_nonreal(d3.nonreal()*(xxx));
      dref.set_real(std::lgamma(d4.real()));
      dref.set_nonreal(d4.nonreal()*(xxx));
      d0=std::erf(std::lgamma(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));

      // erfc(dual) test
      d3=d1;
      dref.set_real(std::erfc(d3.real()));
      dref.set_nonreal(d3.nonreal()*(xxx));
      d0=std::erfc(d3);
      TEST_ASSERT(d0.nearly(dref, tol));

      // erfc(expression) test
      d3=d1+d2;
      dref.set_real(std::erfc(d3.real()));
      dref.set_nonreal(d3.nonreal()*(xxx));
      d0=std::erfc(d1+d2);
      TEST_ASSERT(d0.nearly(dref, tol));

      // erfc(erfc(expression)) test
      d3=d1+d2;
      d4.set_real(std::erfc(d3.real()));
      d4.set_nonreal(d3.nonreal()*(xxx));
      dref.set_real(std::erfc(d4.real()));
      dref.set_nonreal(d4.nonreal()*(xxx));
      d0=std::erfc(std::erfc(d1+d2));
      TEST_ASSERT(d0.nearly(dref, tol));
#endif
    }

    void dual_gamma_nostd_test()
    {
      // TODO: ADD lgamma test when know how to implement derivative

      // TODO: ADD tgamma test when know how to implement derivative
    }

    // TODO: Implement text write
    // TODO: Implement binary read
    // TODO: Implement binary write
};

#endif

