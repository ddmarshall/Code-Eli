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

#ifndef eli_mutil_ad_dual_operators_hpp
#define eli_mutil_ad_dual_operators_hpp

#include "eli/code_eli.hpp"

#include "eli/util/traits.hpp"

#include "eli/mutil/ad/dual_number.hpp"

//
// start add operation
//
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      namespace dual_number
      {
        template <typename left__, typename right__>
        class add
        {
          public:
            typedef typename left__::data_type left_data_type;
            typedef typename right__::data_type right_data_type;
            typedef typename eli::util::promote_traits<left_data_type, right_data_type>::promote_t data_type;

          private:
            typename eli::util::traits<left__>::const_expr_ref l_val;
            typename eli::util::traits<right__>::const_expr_ref r_val;

          public:
            add(const left__ &l, const right__ &r) : l_val(l), r_val(r) {}

            data_type real() const
            {
              return l_val.real()+r_val.real();
            }

            data_type nonreal() const
            {
              return l_val.nonreal()+r_val.nonreal();
            }
        };
      }
    }
  }
}

// dual + dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::add<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                         eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator+(const eli::mutil::ad::dual<left__, comp_real_only> &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::add<eli::mutil::ad::dual<left__, comp_real_only>,
                                           eli::mutil::ad::dual<right__, comp_real_only> > add_t;

  return eli::mutil::ad::dual_number::expression<add_t, comp_real_only>(add_t(l, r));
}

// dual + constant
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::add<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                         eli::mutil::ad::dual_number::constant<right__> >, comp_real_only>
  operator+(const eli::mutil::ad::dual<left__, comp_real_only> &l, const right__ &r)
{
  typedef eli::mutil::ad::dual_number::add<eli::mutil::ad::dual<left__, comp_real_only>,
                                           eli::mutil::ad::dual_number::constant<right__> > add_t;

  return eli::mutil::ad::dual_number::expression<add_t, comp_real_only>(add_t(l, eli::mutil::ad::dual_number::constant<right__>(r)));
}

// constant + dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::add<eli::mutil::ad::dual_number::constant<left__>,
                                                                         eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator+(const left__ &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::add<eli::mutil::ad::dual_number::constant<left__>,
                                           eli::mutil::ad::dual<right__, comp_real_only> > add_t;

  return eli::mutil::ad::dual_number::expression<add_t, comp_real_only>(add_t(l, r));
}

// expression + dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::add<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                         eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator+(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::add<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                           eli::mutil::ad::dual<right__, comp_real_only> > add_t;

  return eli::mutil::ad::dual_number::expression<add_t, comp_real_only>(add_t(l, r));
}

// dual + expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::add<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                         eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator+(const eli::mutil::ad::dual<left__, comp_real_only> &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::add<eli::mutil::ad::dual<left__, comp_real_only>,
                                           eli::mutil::ad::dual_number::expression<right__, comp_real_only> > add_t;

  return eli::mutil::ad::dual_number::expression<add_t, comp_real_only>(add_t(l, r));
}

// expression + constant
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::add<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                         eli::mutil::ad::dual_number::constant<right__> >, comp_real_only>
  operator+(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const right__ &r)
{
  typedef eli::mutil::ad::dual_number::add<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                           eli::mutil::ad::dual_number::constant<right__> > add_t;

  return eli::mutil::ad::dual_number::expression<add_t, comp_real_only>(add_t(l, r));
}

// constant + expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::add<eli::mutil::ad::dual_number::constant<left__>,
                                                                         eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator+(const left__ &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::add<eli::mutil::ad::dual_number::constant<left__>,
                                           eli::mutil::ad::dual_number::expression<right__, comp_real_only> > add_t;

  return eli::mutil::ad::dual_number::expression<add_t, comp_real_only>(add_t(l, r));
}

// expression + expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::add<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                         eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator+(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::add<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                           eli::mutil::ad::dual_number::expression<right__, comp_real_only> > add_t;

  return eli::mutil::ad::dual_number::expression<add_t, comp_real_only>(add_t(l, r));
}
//
// end add operation
//

//
// start subtract operation
//
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      namespace dual_number
      {
        template <typename left__, typename right__>
        class subtract
        {
          public:
            typedef typename left__::data_type left_data_type;
            typedef typename right__::data_type right_data_type;
            typedef typename eli::util::promote_traits<left_data_type, right_data_type>::promote_t data_type;

          private:
            typename eli::util::traits<left__>::const_expr_ref l_val;
            typename eli::util::traits<right__>::const_expr_ref r_val;

          public:
            subtract(const left__ &l, const right__ &r) : l_val(l), r_val(r) {}

            data_type real() const
            {
              return l_val.real()-r_val.real();
            }

            data_type nonreal() const
            {
              return l_val.nonreal()-r_val.nonreal();
            }
        };
      }
    }
  }
}

// dual - dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator-(const eli::mutil::ad::dual<left__, comp_real_only> &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual<left__, comp_real_only>,
                                                eli::mutil::ad::dual<right__, comp_real_only> > sub_t;

  return eli::mutil::ad::dual_number::expression<sub_t, comp_real_only>(sub_t(l, r));
}

// dual - constant
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual_number::constant<right__> >, comp_real_only>
  operator-(const eli::mutil::ad::dual<left__, comp_real_only> &l, const right__ &r)
{
  typedef eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual<left__, comp_real_only>,
                                                eli::mutil::ad::dual_number::constant<right__> > sub_t;

  return eli::mutil::ad::dual_number::expression<sub_t, comp_real_only>(sub_t(l, eli::mutil::ad::dual_number::constant<right__>(r)));
}

// constant - dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual_number::constant<left__>,
                                                                              eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator-(const left__ &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual_number::constant<left__>,
                                                eli::mutil::ad::dual<right__, comp_real_only> > sub_t;

  return eli::mutil::ad::dual_number::expression<sub_t, comp_real_only>(sub_t(l, r));
}

// expression - dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator-(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                eli::mutil::ad::dual<right__, comp_real_only> > sub_t;

  return eli::mutil::ad::dual_number::expression<sub_t, comp_real_only>(sub_t(l, r));
}

// dual - expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator-(const eli::mutil::ad::dual<left__, comp_real_only> &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual<left__, comp_real_only>,
                                                eli::mutil::ad::dual_number::expression<right__, comp_real_only> > sub_t;

  return eli::mutil::ad::dual_number::expression<sub_t, comp_real_only>(sub_t(l, r));
}

// expression - constant
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual_number::constant<right__> >, comp_real_only>
  operator-(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const right__ &r)
{
  typedef eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                eli::mutil::ad::dual_number::constant<right__> > sub_t;

  return eli::mutil::ad::dual_number::expression<sub_t, comp_real_only>(sub_t(l, r));
}

// constant - expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual_number::constant<left__>,
                                                                              eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator-(const left__ &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual_number::constant<left__>,
                                                eli::mutil::ad::dual_number::expression<right__, comp_real_only> > sub_t;

  return eli::mutil::ad::dual_number::expression<sub_t, comp_real_only>(sub_t(l, r));
}

// expression - expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator-(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::subtract<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                eli::mutil::ad::dual_number::expression<right__, comp_real_only> > sub_t;

  return eli::mutil::ad::dual_number::expression<sub_t, comp_real_only>(sub_t(l, r));
}
//
// end subtract operation
//

//
// start multiply operation
//
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      namespace dual_number
      {
        template <typename left__, typename right__>
        class multiply
        {
          public:
            typedef typename left__::data_type left_data_type;
            typedef typename right__::data_type right_data_type;
            typedef typename eli::util::promote_traits<left_data_type, right_data_type>::promote_t data_type;

          private:
            typename eli::util::traits<left__>::const_expr_ref l_val;
            typename eli::util::traits<right__>::const_expr_ref r_val;

          public:
            multiply(const left__ &l, const right__ &r) : l_val(l), r_val(r) {}

            data_type real() const
            {
              return l_val.real()*r_val.real();
            }

            data_type nonreal() const
            {
              return l_val.real()*r_val.nonreal()+l_val.nonreal()*r_val.real();
            }
        };
      }
    }
  }
}

// dual * dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator*(const eli::mutil::ad::dual<left__, comp_real_only> &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual<left__, comp_real_only>,
                                                eli::mutil::ad::dual<right__, comp_real_only> > mult_t;

  return eli::mutil::ad::dual_number::expression<mult_t, comp_real_only>(mult_t(l, r));
}

// dual * constant
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual_number::constant<right__> >, comp_real_only>
  operator*(const eli::mutil::ad::dual<left__, comp_real_only> &l, const right__ &r)
{
  typedef eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual<left__, comp_real_only>,
                                                eli::mutil::ad::dual_number::constant<right__> > mult_t;

  return eli::mutil::ad::dual_number::expression<mult_t, comp_real_only>(mult_t(l, eli::mutil::ad::dual_number::constant<right__>(r)));
}

// constant * dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual_number::constant<left__>,
                                                                              eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator*(const left__ &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual_number::constant<left__>,
                                                eli::mutil::ad::dual<right__, comp_real_only> > mult_t;

  return eli::mutil::ad::dual_number::expression<mult_t, comp_real_only>(mult_t(l, r));
}

// expression * dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator*(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                eli::mutil::ad::dual<right__, comp_real_only> > mult_t;

  return eli::mutil::ad::dual_number::expression<mult_t, comp_real_only>(mult_t(l, r));
}

// dual * expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator*(const eli::mutil::ad::dual<left__, comp_real_only> &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual<left__, comp_real_only>,
                                                eli::mutil::ad::dual_number::expression<right__, comp_real_only> > mult_t;

  return eli::mutil::ad::dual_number::expression<mult_t, comp_real_only>(mult_t(l, r));
}

// expression * constant
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual_number::constant<right__> >, comp_real_only>
  operator*(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const right__ &r)
{
  typedef eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                eli::mutil::ad::dual_number::constant<right__> > mult_t;

  return eli::mutil::ad::dual_number::expression<mult_t, comp_real_only>(mult_t(l, r));
}

// constant * expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual_number::constant<left__>,
                                                                              eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator*(const left__ &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual_number::constant<left__>,
                                                eli::mutil::ad::dual_number::expression<right__, comp_real_only> > mult_t;

  return eli::mutil::ad::dual_number::expression<mult_t, comp_real_only>(mult_t(l, r));
}

// expression * expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                              eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator*(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::multiply<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                eli::mutil::ad::dual_number::expression<right__, comp_real_only> > mult_t;

  return eli::mutil::ad::dual_number::expression<mult_t, comp_real_only>(mult_t(l, r));
}
//
// end multiply operation
//

//
// start divide operation
//
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      namespace dual_number
      {
        template <typename left__, typename right__>
        class divide
        {
          public:
            typedef typename left__::data_type left_data_type;
            typedef typename right__::data_type right_data_type;
            typedef typename eli::util::promote_traits<left_data_type, right_data_type>::promote_t data_type;

          private:
            typename eli::util::traits<left__>::const_expr_ref l_val;
            typename eli::util::traits<right__>::const_expr_ref r_val;

          public:
            divide(const left__ &l, const right__ &r) : l_val(l), r_val(r) {}

            data_type real() const
            {
              return l_val.real()/r_val.real();
            }

            data_type nonreal() const
            {
              right_data_type c(r_val.real());
              return (l_val.nonreal()*c-l_val.real()*r_val.nonreal())/c/c;
            }
        };
      }
    }
  }
}

// dual / dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                            eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator/(const eli::mutil::ad::dual<left__, comp_real_only> &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual<left__, comp_real_only>,
                                              eli::mutil::ad::dual<right__, comp_real_only> > div_t;

  return eli::mutil::ad::dual_number::expression<div_t, comp_real_only>(div_t(l, r));
}

// dual / constant
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                            eli::mutil::ad::dual_number::constant<right__> >, comp_real_only>
  operator/(const eli::mutil::ad::dual<left__, comp_real_only> &l, const right__ &r)
{
  typedef eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual<left__, comp_real_only>,
                                              eli::mutil::ad::dual_number::constant<right__> > div_t;

  return eli::mutil::ad::dual_number::expression<div_t, comp_real_only>(div_t(l, eli::mutil::ad::dual_number::constant<right__>(r)));
}

// constant / dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual_number::constant<left__>,
                                                                            eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator/(const left__ &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual_number::constant<left__>,
                                              eli::mutil::ad::dual<right__, comp_real_only> > div_t;

  return eli::mutil::ad::dual_number::expression<div_t, comp_real_only>(div_t(l, r));
}

// expression / dual
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                            eli::mutil::ad::dual<right__, comp_real_only> >, comp_real_only>
  operator/(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const eli::mutil::ad::dual<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                              eli::mutil::ad::dual<right__, comp_real_only> > div_t;

  return eli::mutil::ad::dual_number::expression<div_t, comp_real_only>(div_t(l, r));
}

// dual / expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual<left__, comp_real_only>,
                                                                            eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator/(const eli::mutil::ad::dual<left__, comp_real_only> &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual<left__, comp_real_only>,
                                              eli::mutil::ad::dual_number::expression<right__, comp_real_only> > div_t;

  return eli::mutil::ad::dual_number::expression<div_t, comp_real_only>(div_t(l, r));
}

// expression / constant
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                            eli::mutil::ad::dual_number::constant<right__> >, comp_real_only>
  operator/(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const right__ &r)
{
  typedef eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                              eli::mutil::ad::dual_number::constant<right__> > div_t;

  return eli::mutil::ad::dual_number::expression<div_t, comp_real_only>(div_t(l, r));
}

// constant / expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual_number::constant<left__>,
                                                                            eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator/(const left__ &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual_number::constant<left__>,
                                              eli::mutil::ad::dual_number::expression<right__, comp_real_only> > div_t;

  return eli::mutil::ad::dual_number::expression<div_t, comp_real_only>(div_t(l, r));
}

// expression / expression
template <typename left__, typename right__, bool comp_real_only>
eli::mutil::ad::dual_number::expression<eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                                                            eli::mutil::ad::dual_number::expression<right__, comp_real_only> >, comp_real_only>
  operator/(const eli::mutil::ad::dual_number::expression<left__, comp_real_only> &l, const eli::mutil::ad::dual_number::expression<right__, comp_real_only> &r)
{
  typedef eli::mutil::ad::dual_number::divide<eli::mutil::ad::dual_number::expression<left__, comp_real_only>,
                                              eli::mutil::ad::dual_number::expression<right__, comp_real_only> > div_t;

  return eli::mutil::ad::dual_number::expression<div_t, comp_real_only>(div_t(l, r));
}
//
// end divide operation
//
#endif
