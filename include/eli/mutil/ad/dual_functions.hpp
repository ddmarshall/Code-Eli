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

#ifndef eli_mutil_ad_dual_functions_hpp
#define eli_mutil_ad_dual_functions_hpp

#include "eli/util/traits.hpp"

#include "eli/mutil/ad/dual_number.hpp"
#include "eli/mutil/ad/dual_operators.hpp"

namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      namespace dual_number
      {
        template <typename data__, typename fun_ob__>
        class unary_fun
        {
          public:
            typedef typename data__::data_type data_type;

          private:
            typename eli::util::traits<data__>::const_expr_ref val;

          public:
            unary_fun(const data__ &v) : val(v) {}

            data_type real() const
            {
              return fun_ob__::f(val.real());
            }

            data_type nonreal() const
            {
              return val.nonreal()*fun_ob__::fp(val.real());
            }
        };
      }
    }
  }
}

#define ELI_AD_DUAL_UNARY_OP_HELPER(x)                                                                                    \
template <typename data__, bool comp_real_only>                                                                             \
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::unary_fun<eli::mutil::ad::dual<data__, comp_real_only>,                     \
                                                                               eli::mutil::ad::x ## _fun<data__> >, comp_real_only>              \
  x (const eli::mutil::ad::dual<data__, comp_real_only> &v)                                                                        \
{                                                                                                                           \
  typedef eli::mutil::ad::dual_number::unary_fun<eli::mutil::ad::dual<data__, comp_real_only>,eli::mutil::ad::x ## _fun<data__> > un_fun_t;      \
                                                                                                                            \
  return eli::mutil::ad::dual_number::expression<un_fun_t, comp_real_only>(un_fun_t(v));                                           \
}                                                                                                                           \
template <typename data__, bool comp_real_only>                                                                             \
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::unary_fun<eli::mutil::ad::dual_number::expression<data__, comp_real_only>,  \
                                                                               eli::mutil::ad::x ## _fun<typename eli::mutil::ad::dual_number::expression<data__, comp_real_only>::data_type > >, \
                                                                               comp_real_only>                                                                          \
  x (const eli::mutil::ad::dual_number::expression<data__, comp_real_only> &v)                                                     \
{                                                                                                                           \
  typedef eli::mutil::ad::dual_number::unary_fun<eli::mutil::ad::dual_number::expression<data__, comp_real_only>,                         \
                                                 eli::mutil::ad::x ## _fun<typename eli::mutil::ad::dual_number::expression<data__, comp_real_only>::data_type > > un_fun_t; \
                                                                                                                            \
  return eli::mutil::ad::dual_number::expression<un_fun_t, comp_real_only>(un_fun_t(v));                                           \
}

#endif
