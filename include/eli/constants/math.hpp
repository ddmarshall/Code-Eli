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

#ifdef ELI_QD_FOUND
    template<>
    class math <dd_real>
    {
      public:
        static dd_real exp()     {return dd_real::_e;}
        static dd_real ln_two()  {return dd_real(0.6931471805599453, 2.3190468138462999e-17);}
        static dd_real log_exp() {return dd_real(0.4342944819032518, 1.0983196502167658e-17);}

        static dd_real pi()             {return dd_real::_pi;}
        static dd_real two_pi()         {return dd_real::_2pi;}
        static dd_real pi_by_two()      {return dd_real::_pi2;}
        static dd_real pi_by_four()     {return dd_real::_pi4;}
        static dd_real pi_squared()     {return dd_real(9.8696044010893586,  6.265295508739711e-16);}
        static dd_real pi_cubed()       {return dd_real(31.006276680299820,  4.1641946985288283e-16);}
        static dd_real sqrt_pi()        {return dd_real(1.7724538509055161, -7.666586499825800e-17);}
        static dd_real cbrt_pi()        {return dd_real(1.4645918875615234, -1.0150473108658175e-16);}
        static dd_real one_by_pi()      {return dd_real(0.31830988618379067,-1.9678676675182486e-17);}
        static dd_real two_by_pi()      {return dd_real(0.63661977236758134,-3.9357353350364972e-17);}
        static dd_real one_by_sqrt_pi() {return dd_real(0.56418958354775629, 7.6677298065829422e-18);}
        static dd_real two_by_sqrt_pi() {return dd_real(1.1283791670955126,  1.5335459613165884e-17);}

        static dd_real sqrt_two()        {return dd_real(1.41421356237309504, -9.667293313452916e-17);}
        static dd_real sqrt_two_by_two() {return dd_real(0.70710678118654752, -4.833646656726458e-17);}
    };

    template<>
    class math <qd_real>
    {
      public:
        static qd_real exp()     {return qd_real::_e;}
        static qd_real ln_two()  {return qd_real::_log2;}
        static qd_real log_exp() {return qd_real(0.4342944819032518,     1.0983196502167651e-17,
                                                 3.7171812331109590e-34, 7.7344843465042927e-51);}

        static qd_real pi()             {return qd_real::_pi;}
        static qd_real two_pi()         {return qd_real::_2pi;}
        static qd_real pi_by_two()      {return qd_real::_pi2;}
        static qd_real pi_by_four()     {return qd_real::_pi4;}
        static qd_real pi_squared()     {return qd_real(9.8696044010893586,     6.265295508739711e-16,
                                                        3.730017701459809e-32, -1.3336037774250846e-48);}
        static qd_real pi_cubed()       {return qd_real(31.006276680299820,     4.1641946985288347e-16,
                                                        2.1094349902320112e-32, 9.749332643144236e-49);}
        static qd_real sqrt_pi()        {return qd_real(1.7724538509055161,    -7.666586499825799e-17,
                                                       -1.3058334907945429e-33,-2.6110142087827117e-50);}
        static qd_real cbrt_pi()        {return qd_real(1.4645918875615234,    -1.0150473108658177e-16,
                                                       -1.3864695806293404e-33, 1.8314339704237657e-50);}
        static qd_real one_by_pi()      {return qd_real(0.31830988618379067,   -1.9678676675182486e-17,
                                                       -1.0721436282893004e-33, 8.0535639265941101e-50);}
        static qd_real two_by_pi()      {return qd_real(0.63661977236758134,   -3.9357353350364972e-17,
                                                       -2.1442872565786008e-33, 1.6107127853188220e-49);}
        static qd_real one_by_sqrt_pi() {return qd_real(0.56418958354775629,    7.6677298065829406e-18,
                                                       -2.3828422983468431e-34,-1.0038973308276328e-50);}
        static qd_real two_by_sqrt_pi() {return qd_real(1.1283791670955126,     1.5335459613165881e-17,
                                                       -4.7656845966936863e-34,-2.0077946616552656e-50);}

        static qd_real sqrt_two()        {return qd_real(1.41421356237309504,   -9.667293313452913e-17,
                                                         4.1386753086994136e-33, 4.935546991468351e-50);}
        static qd_real sqrt_two_by_two() {return qd_real(0.70710678118654752,   -4.833646656726457e-17,
                                                         2.0693376543497068e-33, 2.4677734957341755e-50);}
    };
#endif
  }
}

#endif
