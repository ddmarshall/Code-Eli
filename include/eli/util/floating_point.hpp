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

#ifndef eli_util_floating_point_hpp
#define eli_util_floating_point_hpp

#include <iostream>

#include "traits.hpp"

namespace eli
{
  namespace util
  {
    union float_type
    {
      typedef uint32_t integer_type;
      integer_type i;
      struct
      {
        uint64_t mantissa : 23;
        uint32_t exponent : 8;
        uint32_t sign : 1;
      } parts;
    };

    struct double_type
    {
      typedef uint64_t integer_type;
      integer_type i;
      struct
      {
        uint64_t mantissa : 52;
        uint32_t exponent : 11;
        uint32_t sign : 1;
      } parts;
    };

    union long_double_type
    {
#if defined(_WIN32)
      typedef uint64_t integer_type;
      integer_type i;
      struct
      {
        uint64_t mantissa : 52;
        uint32_t exponent : 11;
        uint32_t sign : 1;
      } parts;
#else
      typedef uint64_t integer_type;
      integer_type i[2];
      struct
      {
        uint64_t mantissa : 63;
        uint32_t integer_bit : 1;
        uint32_t exponent : 15;
        uint32_t sign : 1;
      } parts;
#endif
    };

    const float_type * set_floating_point_type(const float * pcf)
    {
      return static_cast<const float_type *>(static_cast<const void *>(pcf));
    }

    float_type * set_floating_point_type(float *pf)
    {
      return static_cast<float_type *>(static_cast<void *>(pf));
    }

    const double_type * set_floating_point_type(const double * pcf)
    {
      return static_cast<const double_type *>(static_cast<const void *>(pcf));
    }

    double_type * set_floating_point_type(double *pf)
    {
      return static_cast<double_type *>(static_cast<void *>(pf));
    }

    const long_double_type * set_floating_point_type(const long double * pcf)
    {
      return static_cast<const long_double_type *>(static_cast<const void *>(pcf));
    }

    long_double_type * set_floating_point_type(long double *pf)
    {
      return static_cast<long_double_type *>(static_cast<void *>(pf));
    }

    std::ostream & operator<<(std::ostream &ostr, const float_type &ft)
    {
      float_type::integer_type sign_shift(31);
      float_type::integer_type mantissa_mask(1); mantissa_mask=(mantissa_mask<<23) - 1;
      float_type::integer_type exponent_shift(23);
      float_type::integer_type exponent_mask(0xff);

      ostr << std::hex << "0x"  << (ft.i >> sign_shift)
                       << " 0x" << (ft.i & mantissa_mask)
                       << " 0x" << ((ft.i >> exponent_shift) & exponent_mask)
           << std::dec;

      return ostr;
    }

    std::ostream & operator<<(std::ostream &ostr, const double_type &ft)
    {
      double_type::integer_type sign_shift(63);
      double_type::integer_type mantissa_mask(1); mantissa_mask=(mantissa_mask<<52) - 1;
      double_type::integer_type exponent_shift(52);
      double_type::integer_type exponent_mask(0x7ffb);

      ostr << std::hex << "0x"  << (ft.i >> sign_shift)
                       << " 0x" << (ft.i & mantissa_mask)
                       << " 0x" << ((ft.i >> exponent_shift) & exponent_mask)
           << std::dec;

      return ostr;
    }

    std::ostream & operator<<(std::ostream &ostr, const long_double_type &ft)
    {
#if defined(_WIN32)
      double_type::integer_type sign_shift(63);
      double_type::integer_type mantissa_mask(1); mantissa_mask=(mantissa_mask<<52) - 1;
      double_type::integer_type exponent_shift(52);
      double_type::integer_type exponent_mask(0x7ffb);

      ostr << std::hex << "0x"  << (ft.i >> sign_shift)
                       << " 0x" << (ft.i & mantissa_mask)
                       << " 0x" << ((ft.i >> exponent_shift) & exponent_mask)
           << std::dec;

      return ostr;
#else
      double_type::integer_type mantissa_mask(1); mantissa_mask=0x7FFFFFFFFFFFFFFF;
      double_type::integer_type exponent_shift(0);
      double_type::integer_type exponent_mask(0x7fff);

      ostr << std::hex << "0x"  << ft.parts.sign
                       << " 0x" << ft.parts.integer_bit
                       << " 0x" << (ft.i[0] & mantissa_mask)
                       << " 0x" << ((ft.i[1] >> exponent_shift) & exponent_mask)
           << std::dec;

      return ostr;
#endif
    }

    template<typename data__>
    data__ increment_ulp(const data__ &/*d*/, const int &/*n_ulp*/)
    {
      static_assert(always_false<data__>::value, "Function not specialized for given type");
    }

    template<>
    float increment_ulp<float>(const float &f, const int &n_ulp)
    {
#if defined(__GNUC__) && defined(NDEBUG) && !defined(__clang__) && (__GNUC__==4) && (__GNUC_MINOR__==5)
      volatile
#endif
      int32_t *pi;
      float fr(f);

      pi=static_cast<int32_t *>(static_cast<void *>(&fr));
      (*pi)+=n_ulp;
      return fr;
    }

    template<>
    double increment_ulp<double>(const double &f, const int &n_ulp)
    {
#if defined(__GNUC__) && defined(NDEBUG) && !defined(__clang__) && (__GNUC__==4) && (__GNUC_MINOR__==5)
      volatile
#endif
      int64_t *pi;
      double fr(f);

      pi=static_cast<int64_t *>(static_cast<void *>(&fr));
      (*pi)+=n_ulp;
      return fr;
    }

    template<>
    long double increment_ulp<long double>(const long double &f, const int &n_ulp)
    {
#if defined(_WIN32)
      int64_t *pi;
      double fr(f);

      pi=static_cast<int64_t *>(static_cast<void *>(&fr));
      (*pi)+=n_ulp;
      return fr;
#else
      int32_t orig_integer_bit;
#if defined(__GNUC__) && defined(NDEBUG) && !defined(__clang__) && (__GNUC__==4) && (__GNUC_MINOR__==5)
      volatile
#endif
      int64_t *pi;
#if defined(__GNUC__) && defined(NDEBUG) && !defined(__clang__) && (__GNUC__==4) && (__GNUC_MINOR__==5)
      volatile
#endif
      long_double_type *pft;
      long double fr(f);

      pi=static_cast<int64_t *>(static_cast<void *>(&fr));
      pft=static_cast<long_double_type *>(static_cast<void *>(&fr));

      orig_integer_bit=pft->parts.integer_bit;
      pi[0]+=n_ulp;

      // check if wrapped mantissa and correct
      if ((pft->parts.integer_bit!=orig_integer_bit) && (orig_integer_bit==1))
      {
        pft->parts.integer_bit=1;
        if (n_ulp>0)
        {
          ++(pft->parts.exponent);
        }
        else
        {
          --(pft->parts.exponent);
        }
      }
      return fr;
#endif
    }
  }
}
#endif
