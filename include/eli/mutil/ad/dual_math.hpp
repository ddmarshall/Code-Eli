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

#ifndef eli_mutil_ad_dual_math_hpp
#define eli_mutil_ad_dual_math_hpp

#include "eli/constants/math.hpp"
#include "eli/mutil/ad/dual_functions.hpp"

#include <cmath>

// sin and std::sin
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct sin_fun
      {
        static T__ f(const T__ &t)  {return std::sin(t);}
        static T__ fp(const T__ &t) {return std::cos(t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n%4)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            case(2):
            {
              return -std::sin(t);
              break;
            }
            default:
            case(3):
            {
              return -std::cos(t);
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(sin) }
using std::sin;

// cos and std::cos
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct cos_fun
      {
        static T__ f(const T__ &t)  {return std::cos(t);}
        static T__ fp(const T__ &t) {return -std::sin(t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n%4)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            case(2):
            {
              return -std::cos(t);
              break;
            }
            default:
            case(3):
            {
              return std::sin(t);
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(cos) }
using std::cos;

// tan and std::tan
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct tan_fun
      {
        static T__ f(const T__ &t)  {return std::tan(t);}
        static T__ fp(const T__ &t) {T__ v(std::cos(t)); return 1/v/v;}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              // TODO: IMPLEMENT THIS ALGORITHM
              // d(tan)=1+tan^2
              // d^n(tan)=n*tan^(n-1)+n*tan^(n+1)
              assert(false);
              return 0;
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(tan) }
using std::tan;

// asin and std::asin
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct asin_fun
      {
        static T__ f(const T__ &t)  {return std::asin(t);}
        static T__ fp(const T__ &t) {return 1/std::sqrt(1-t*t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            case(2):
            {
              // TODO: IMPLEMENT THIS
              assert(false);
              return 0;
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(asin) }
using std::asin;

// acos and std::acos
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct acos_fun
      {
        static T__ f(const T__ &t)  {return std::acos(t);}
        static T__ fp(const T__ &t) {return -1/std::sqrt(1-t*t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            case(2):
            {
              // TODO: IMPLEMENT THIS
              assert(false);
              return 0;
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(acos) }
using std::acos;

// atan and std::atan
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct atan_fun
      {
        static T__ f(const T__ &t)  {return std::atan(t);}
        static T__ fp(const T__ &t) {return 1/(1+t*t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            case(2):
            {
              // TODO: IMPLEMENT THIS
              assert(false);
              return 0;
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(atan) }
using std::atan;

// atan2 and std::atan2
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct atan2_fun
      {
        static T__ f(const T__ &y, const T__ &x)  {return std::atan2(y, x);}
        static T__ fx(const T__ &y, const T__ &x) {return -y/(x*x+y*y);}
        static T__ fy(const T__ &y, const T__ &x) {return x/(x*x+y*y);}
        static T__ fxy(const T__ &y, const T__ &x, const size_t &nx, const size_t &ny)
        {
          switch(nx)
          {
            case(0):
            {
              switch(ny)
              {
                case(0):
                {
                  return f(y,x);
                  break;
                }
                case(1):
                {
                  return fy(y,x);
                  break;
                }
                default:
                {
                  // TODO: NEED TO IMPLEMENT THIS
                  assert(false);
                  return 0;
                  break;
                }
              }
            }
            case(1):
            {
              switch(ny)
              {
                case(0):
                {
                  return fx(y,x);
                  break;
                }
                default:
                {
                  // TODO: NEED TO IMPLEMENT THIS
                  assert(false);
                  return 0;
                  break;
                }
              }
            }
            default:
            case(2):
            {
              // TODO: IMPLEMENT THIS
              assert(false);
              return 0;
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(atan2) }
using std::atan2;


// sinh and std::sinh
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct sinh_fun
      {
        static T__ f(const T__ &t)  {return std::sinh(t);}
        static T__ fp(const T__ &t) {return std::cosh(t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          if ((n%2)==0)
            return f(t);
          else
            return fp(t);
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(sinh) }
using std::sinh;

// cosh and std::cosh
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct cosh_fun
      {
        static T__ f(const T__ &t)  {return std::cosh(t);}
        static T__ fp(const T__ &t) {return std::sinh(t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          if ((n%2)==0)
            return f(t);
          else
            return fp(t);
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(cosh) }
using std::cosh;

// tanh and std::tanh
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct tanh_fun
      {
        static T__ f(const T__ &t)  {return std::tanh(t);}
        static T__ fp(const T__ &t) {T__ v(std::cosh(t)); return 1/v/v;}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              // TODO: IMPLEMENT THIS ALGORITHM
              // d(tanh)=1-tanh^2
              // d(tanh^n)=n*tanh^(n-1)-n*tanh^(n+1)
              assert(false);
              return 0;
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(tanh) }
using std::tanh;

// asinh and std::asinh
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct asinh_fun
      {
        static T__ f(const T__ &t)  {return std::asinh(t);}
        static T__ fp(const T__ &t) {return 1/std::sqrt(1+t*t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            case(2):
            {
              // TODO: IMPLEMENT THIS
              assert(false);
              return 0;
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(asinh) }
using std::asinh;

// acos and std::acosh
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct acosh_fun
      {
        static T__ f(const T__ &t)  {return std::acosh(t);}
        static T__ fp(const T__ &t) {return 1/std::sqrt(t*t-1);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            case(2):
            {
              // TODO: IMPLEMENT THIS
              assert(false);
              return 0;
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(acosh) }
using std::acosh;

// atanh and std::atanh
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct atanh_fun
      {
        static T__ f(const T__ &t)  {return std::atanh(t);}
        static T__ fp(const T__ &t) {return 1/(1-t*t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            case(2):
            {
              // TODO: IMPLEMENT THIS
              assert(false);
              return 0;
              break;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(atanh) }
using std::atanh;

// exp and std::exp
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct exp_fun
      {
        static T__ f(const T__ &t)  {return std::exp(t);}
        static T__ fp(const T__ &t) {return std::exp(t);}
        static T__ fp(const T__ &t, const size_t &/*n*/) {return std::exp(t);}
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(exp) }
using std::exp;

// expm1 and std::expm1
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct expm1_fun
      {
        static T__ f(const T__ &t)  {return std::expm1(t);}
        static T__ fp(const T__ &t) {return std::exp(t);}
        static T__ fp(const T__ &t, const size_t &/*n*/) {return std::exp(t);}
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(expm1) }
using std::expm1;

// exp2 and std::exp2
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct exp2_fun
      {
        static T__ f(const T__ &t)  {return std::exp2(t);}
        static T__ fp(const T__ &t) {return std::log(static_cast<T__>(2))*std::exp2(t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch (n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              return std::pow(std::log(static_cast<T__>(2)),n-1)*std::exp2(t);
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(exp2) }
using std::exp2;

// log and std::log
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct log_fun
      {
        static T__ f(const T__ &t)  {return std::log(t);}
        static T__ fp(const T__ &t) {return 1/t;}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              T__ val(fp(t));
              for (size_t i=2; i<=n; ++i)
                val/=-t/i;
              return val;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(log) }
using std::log;

// log10 and std::log10
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct log10_fun
      {
        static T__ f(const T__ &t)  {return std::log10(t);}
        static T__ fp(const T__ &t) {return 1/t/std::log(static_cast<T__>(10));}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              T__ val(fp(t)/std::log(static_cast<T__>(10)));
              for (size_t i=2; i<=n; ++i)
                val/=-t/i;
              return val;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(log10) }
using std::log10;

// log2 and std::log2
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct log2_fun
      {
        static T__ f(const T__ &t)  {return std::log2(t);}
        static T__ fp(const T__ &t) {return 1/t/std::log(static_cast<T__>(2));}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              T__ val(fp(t)/std::log(static_cast<T__>(2)));
              for (size_t i=2; i<=n; ++i)
                val/=-t/i;
              return val;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(log2) }
using std::log2;

// log1p and std::log1p
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct log1p_fun
      {
        static T__ f(const T__ &t)  {return std::log1p(t);}
        static T__ fp(const T__ &t) {return 1/(1+t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              T__ val(fp(t));
              for (size_t i=2; i<=n; ++i)
                val/=-(1+t)/i;
              return val;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(log1p) }
using std::log1p;

// sqrt and std::sqrt
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct sqrt_fun
      {
        static T__ f(const T__ &t)  {return std::sqrt(t);}
        static T__ fp(const T__ &t) {return static_cast<T__>(0.5)/std::sqrt(t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              T__ val(fp(t));
              for (size_t i=1; i<n; ++i)
                val*=0.5*(1-2*i)/t;

              return val;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(sqrt) }
using std::sqrt;

// cbrt and std::cbrt
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct cbrt_fun
      {
        static T__ f(const T__ &t)  {return std::cbrt(t);}
        static T__ fp(const T__ &t) {return (static_cast<T__>(1)/static_cast<T__>(3))/std::cbrt(t*t);}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              T__ val(fp(t));
              for (size_t i=1; i<n; ++i)
                val*=(1-3*i)/t/3;

              return val;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(cbrt) }
using std::cbrt;

// abs and std::abs
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct abs_fun
      {
        static T__ f(const T__ &t)  {return std::abs(t);}
        static T__ fp(const T__ &t) {return (t<0)?(static_cast<T__>(-1)):(static_cast<T__>(1));}
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              return 0;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(abs) }
using std::abs;

// ceil and std::ceil
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct ceil_fun
      {
        static T__ f(const T__ &t)  {return std::ceil(t);}
        static T__ fp(const T__ &/*t*/) {return 0;}
        static T__ fp(const T__ &t, const size_t &n)
        {
          if (n==0)
            return f(t);
          else
            return 0;
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(ceil) }
using std::ceil;

// floor and std::floor
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct floor_fun
      {
        static T__ f(const T__ &t)  {return std::floor(t);}
        static T__ fp(const T__ &/*t*/) {return 0;}
        static T__ fp(const T__ &t, const size_t &n)
        {
          if (n==0)
            return f(t);
          else
            return 0;
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(floor) }
using std::floor;

// erf and std::erf
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct erf_fun
      {
        static T__ f(const T__ &t)  {return std::erf(t);}
        static T__ fp(const T__ &t)
        {
          return static_cast<T__>(2)*std::exp(-t*t)/eli::constants::math<T__>::sqrt_pi();
        }
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              // TODO: IMPLEMENT THIS
              assert(false);
              return 0;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(erf) }
using std::erf;

// erfc and std::erfc
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename T__>
      struct erfc_fun
      {
        static T__ f(const T__ &t)  {return std::erfc(t);}
        static T__ fp(const T__ &t)
        {
          return -2*std::exp(-t*t)/eli::constants::math<T__>::sqrt_pi();
        }
        static T__ fp(const T__ &t, const size_t &n)
        {
          switch(n)
          {
            case(0):
            {
              return f(t);
              break;
            }
            case(1):
            {
              return fp(t);
              break;
            }
            default:
            {
              // TODO: IMPLEMENT THIS
              assert(false);
              return 0;
            }
          }
        }
      };
    }
  }
}
namespace std { ELI_AD_DUAL_UNARY_OP_HELPER(erfc) }
using std::erfc;

#endif
