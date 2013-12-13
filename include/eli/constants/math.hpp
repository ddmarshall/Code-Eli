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

#ifndef eli_constants_math_hpp
#define eli_constants_math_hpp

#include <cmath>

#include "eli/code_eli.hpp"

namespace eli
{
  namespace constants
  {
    template <typename T__>
    class math
    {
    };

    template<>
    class math <float>
    {
      public:
        static float exp()     {return 2.7182818f;}
        static float ln_two()  {return 0.6931472f;}
        static float log_exp() {return 0.43429447f;}

        static float pi()             {return 3.1415927f;}
        static float two_pi()         {return 6.2831853f;}
        static float pi_by_two()      {return 1.57079633f;}
        static float pi_by_four()     {return 0.78539816f;}
        static float pi_squared()     {return 9.8696046f;}
        static float pi_cubed()       {return 31.006279f;}
        static float sqrt_pi()        {return 1.7724539f;}
        static float cbrt_pi()        {return 1.4645919f;}
        static float one_by_pi()      {return 0.31830988f;}
        static float two_by_pi()      {return 0.63661977f;}
        static float one_by_sqrt_pi() {return 0.56418958f;}
        static float two_by_sqrt_pi() {return 1.12837916f;}

        static float sqrt_two()        {return 1.41421356f;}
        static float sqrt_two_by_two() {return 0.70710678f;}
    };

    template <>
    class math <double>
    {
      public:
        static const double exp1;

      public:
        static double exp()     {return 2.718281828459045;}
        static double ln_two()  {return 0.6931471805599453;}
        static double log_exp() {return 0.4342944819032518;}

        static double pi()             {return 3.141592653589793;}
        static double two_pi()         {return 6.283185307179586;}
        static double pi_by_two()      {return 1.5707963267948966;}
        static double pi_by_four()     {return 0.7853981633974483;}
        static double pi_squared()     {return 9.8696044010893586;}
        static double pi_cubed()       {return 31.006276680299817;}
        static double sqrt_pi()        {return 1.7724538509055159;}
#if defined(__clang__)
        static double cbrt_pi()        {return 1.4645918875615233;}
#elif defined(__INTEL_COMPILER)
        static double cbrt_pi()        {return 1.4645918875615231;}
#elif defined(__GNUC__)
# if (__GNUC__==4) && (__GNUC_MINOR__>=7) && (defined NDEBUG)
        static double cbrt_pi()        {return 1.4645918875615231;}
# else
        static double cbrt_pi()        {return 1.4645918875615233;}
# endif
#else
        static double cbrt_pi()        {return 1.4645918875615231;}
#endif
        static double one_by_pi()      {return 0.31830988618379067;}
        static double two_by_pi()      {return 0.63661977236758134;}
        static double one_by_sqrt_pi() {return 0.56418958354775629;}
        static double two_by_sqrt_pi() {return 1.1283791670955126;}

        static double sqrt_two()        {return 1.41421356237309504;}
        static double sqrt_two_by_two() {return 0.70710678118654752;}
    };

    template<>
    class math <long double>
    {
      public:
        static long double exp()     {return 2.7182818284590452354L;}
        static long double ln_two()  {return 0.69314718055994530942L;}
        static long double log_exp() {return 0.43429448190325182766L;}

        static long double pi()             {return 3.1415926535897932385L;}
        static long double two_pi()         {return 6.2831853071795864770L;}
        static long double pi_by_two()      {return 1.57079632679489661923L;}
        static long double pi_by_four()     {return 0.78539816339744830962L;}
        static long double pi_squared()     {return 9.869604401089358619L;}
#ifdef _MSC_VER
        static long double pi_cubed()       {return 31.006276680299817L;}
        static long double sqrt_pi()        {return 1.7724538509055159L;}
        static long double cbrt_pi()        {return 1.4645918875615231L;}
#else
        static long double pi_cubed()       {return 31.0062766802998201763L;}
        static long double sqrt_pi()        {return 1.7724538509055160273L;}
        static long double cbrt_pi()        {return 1.464591887561523263L;}
#endif
        static long double one_by_pi()      {return 0.31830988618379067154L;}
        static long double two_by_pi()      {return 0.6366197723675813431L;}
        static long double one_by_sqrt_pi() {return 0.5641895835477562869L;}
        static long double two_by_sqrt_pi() {return 1.1283791670955125738L;}

        static long double sqrt_two()        {return 1.4142135623730950488L;}
        static long double sqrt_two_by_two() {return 0.7071067811865475244L;}
    };
  }
}

#endif
