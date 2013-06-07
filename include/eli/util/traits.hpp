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

#ifndef eli_util_traits_hpp
#define eli_util_traits_hpp

//
// traits class
//
namespace eli
{
  namespace util
  {
    // always false class
    template<typename T__>
    struct always_false
    {
      enum {value=false};
    };

    // always true class
    template<typename T__>
    struct always_true
    {
      enum {value=true};
    };

    // generic traits class
    template <typename data__>
    class traits
    {
      public:
        // how to refer to a constant reference
        typedef const data__ & const_expr_ref;

        // how to refer to a reference
        typedef data__ & expr_ref;
    };

    struct promoted_type_undefined {};
    template <typename left__, typename right__>
    struct promote_traits
    {
      typedef promoted_type_undefined promote_t;
    };

#define ELI_UTIL_PROMOTE_TRAITS_INT_HELPER(x_type)    \
    template <>                                       \
    struct promote_traits<x_type, x_type>             \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<x_type, char>               \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<char, x_type>               \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<x_type, unsigned char>      \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<unsigned char, x_type>      \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<x_type, short int>          \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<short int, x_type>          \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<x_type, unsigned short int> \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<unsigned short int, x_type> \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<x_type, int>                \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<int, x_type>                \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<x_type, unsigned int>       \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<unsigned int, x_type>       \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<x_type, long int>           \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<long int, x_type>           \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<x_type, unsigned long int>  \
    {                                                 \
      typedef x_type promote_t;                       \
    };                                                \
    template <>                                       \
    struct promote_traits<unsigned long int, x_type>  \
    {                                                 \
      typedef x_type promote_t;                       \
    };

#define ELI_UTIL_PROMOTE_TRAITS_HELPER(x_type)      \
    ELI_UTIL_PROMOTE_TRAITS_INT_HELPER(x_type)      \
    template <>                                     \
    struct promote_traits<x_type, float>            \
    {                                               \
      typedef x_type promote_t;                     \
    };                                              \
    template <>                                     \
    struct promote_traits<float, x_type>            \
    {                                               \
      typedef x_type promote_t;                     \
    };                                              \
    template <>                                     \
    struct promote_traits<x_type, double>           \
    {                                               \
      typedef x_type promote_t;                     \
    };                                              \
    template <>                                     \
    struct promote_traits<double, x_type>           \
    {                                               \
      typedef x_type promote_t;                     \
    };                                              \
    template <>                                     \
    struct promote_traits<x_type, long double>      \
    {                                               \
      typedef x_type promote_t;                     \
    };                                              \
    template <>                                     \
    struct promote_traits<long double, x_type>      \
    {                                               \
      typedef x_type promote_t;                     \
    };

    //
    // traits promotion definitions for float
    //
    ELI_UTIL_PROMOTE_TRAITS_INT_HELPER(float)

    //
    // traits promotion definitions for double
    //
    ELI_UTIL_PROMOTE_TRAITS_INT_HELPER(double)
    // double-float promotions
    template <>
    struct promote_traits<double, float>
    {
      typedef double promote_t;
    };
    template <>
    struct promote_traits<float, double>
    {
      typedef double promote_t;
    };

    //
    // traits promotion definitions for long double
    //
    ELI_UTIL_PROMOTE_TRAITS_INT_HELPER(long double)
    // long double-float promotions
    template <>
    struct promote_traits<long double, float>
    {
      typedef long double promote_t;
    };
    template <>
    struct promote_traits<float, long double>
    {
      typedef long double promote_t;
    };
    // long double-double promotions
    template <>
    struct promote_traits<long double, double>
    {
      typedef long double promote_t;
    };
    template <>
    struct promote_traits<double, long double>
    {
      typedef long double promote_t;
    };

#ifdef ELI_USING_QD
    //
    // traits promotion definitions for QD types
    //
    ELI_UTIL_PROMOTE_TRAITS_HELPER(dd_real)
    ELI_UTIL_PROMOTE_TRAITS_HELPER(qd_real)
#endif
  }
}
#endif
