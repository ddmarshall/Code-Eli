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
    struct float_type
    {
      uint64_t mantissa : 23;
      uint32_t exponent : 8;
      uint32_t sign : 1;
    };

    struct double_type
    {
      uint64_t mantissa : 52;
      uint32_t exponent : 11;
      uint32_t sign : 1;
    };

    struct long_double_type
    {
      uint64_t mantissa : 63;
      uint32_t integer_bit : 1;
      uint32_t exponent : 15;
      uint32_t sign : 1;
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
      ostr << std::hex << "0x" << ft.sign << " 0x" << ft.mantissa << " 0x" << ft.exponent << std::dec;

      return ostr;
    }

    std::ostream & operator<<(std::ostream &ostr, const double_type &ft)
    {
      ostr << std::hex << "0x" << ft.sign << " 0x" << ft.mantissa << " 0x" << ft.exponent << std::dec;

      return ostr;
    }

    std::ostream & operator<<(std::ostream &ostr, const long_double_type &ft)
    {
      ostr << std::hex << "0x" << ft.sign << " 0x" << ft.integer_bit << " 0x" << ft.mantissa << " 0x" << ft.exponent << std::dec;

      return ostr;
    }

    template<typename data__>
    data__ increment_ulp(const data__ &/*d*/, const int &/*n_ulp*/)
    {
      static_assert(always_false<data__>::value, "Function not specialized for given type");
    }

    template<>
    float increment_ulp<float>(const float &f, const int &n_ulp)
    {
      int32_t *pi;
      float fr(f);

      pi=static_cast<int32_t *>(static_cast<void *>(&fr));
      (*pi)+=n_ulp;
      return fr;
    }

    template<>
    double increment_ulp<double>(const double &f, const int &n_ulp)
    {
      int64_t *pi;
      double fr(f);

      pi=static_cast<int64_t *>(static_cast<void *>(&fr));
      (*pi)+=n_ulp;
      return fr;
    }

    template<>
    long double increment_ulp<long double>(const long double &f, const int &n_ulp)
    {
      int32_t orig_integer_bit;
      int64_t *pi;
      long_double_type *pft;
      long double fr(f);

      pi=static_cast<int64_t *>(static_cast<void *>(&fr));
      pft=static_cast<long_double_type *>(static_cast<void *>(&fr));

      orig_integer_bit=pft->integer_bit;
      pi[0]+=n_ulp;

      // check if wrapped mantissa and correct
      if ((pft->integer_bit!=orig_integer_bit) && (orig_integer_bit==1))
      {
        pft->integer_bit=1;
        if (n_ulp>0)
        {
          ++(pft->exponent);
        }
        else
        {
          --(pft->exponent);
        }
      }
      return fr;
    }
  }
}
#endif
