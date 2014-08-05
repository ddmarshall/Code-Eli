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

#ifndef eli_util_tolerance_hpp
#define eli_util_tolerance_hpp

#include <cmath>  // std::cmath
#include <limits> // std::numeric_limits

#include "eli/code_eli.hpp"

namespace eli
{
  namespace util
  {
    template<typename data__>
    class tolerance
    {
      public:
        tolerance() : abs_tol(10000*std::numeric_limits<data__>::epsilon()),
                      rel_tol(std::sqrt(std::numeric_limits<data__>::epsilon()))
        {
        }
        tolerance(const data__ &abs, const data__ &rel) : abs_tol(abs), rel_tol(rel) {}
        tolerance(const tolerance<data__> &tol) : abs_tol(tol.abs_tol), rel_tol(tol.rel_tol) {}
        ~tolerance() {}

        tolerance & operator=(const tolerance<data__> &tol)
        {
          if (this == &tol)
            return *this;

          abs_tol=tol.abs_tol;
          rel_tol=tol.rel_tol;

          return *this;
        }

        bool operator==(const tolerance<data__> &tol) const
        {
          if (this == &tol)
            return true;

          if (abs_tol!=tol.abs_tol)
            return false;

          if (rel_tol!=tol.rel_tol)
            return false;

          return true;
        }

        bool operator!=(const tolerance<data__> &tol) const
        {
          return !operator==(tol);
        }

        data__ get_relative_tolerance() const {return rel_tol;}
        data__ get_absolute_tolerance() const {return abs_tol;}

        template<typename Derived1, typename Derived2>
        bool exactly_equal(const Eigen::MatrixBase<Derived1> &m1, const Eigen::MatrixBase<Derived2> &m2) const
        {
          if (m1.rows()!=m2.rows())
            return false;
          if (m1.cols()!=m2.cols())
            return false;

          return (m1==m2);
        }

        bool exactly_equal(const data__ &t1, const data__ &t2) const
        {
          return (t1==t2);
        }

        template<typename type2__>
        bool exactly_equal(const data__ &t1, const type2__ &t2) const
        {
          return (t1==static_cast<data__>(t2));
        }

        template<typename type1__>
        bool exactly_equal(const type1__ &t1, const data__ &t2) const
        {
          return exactly_equal(t2, t1);
        }

        template<typename Derived1, typename Derived2>
        bool approximately_equal(const Eigen::MatrixBase<Derived1> &m1, const Eigen::MatrixBase<Derived2> &m2) const
        {
          typename Derived1::Index i, j;

          if (m1.rows()!=m2.rows())
            return false;
          if (m1.cols()!=m2.cols())
            return false;

          for (i=0; i<m1.rows(); ++i)
          {
            for (j=0; j<m1.cols(); ++j)
            {
              if (!approximately_equal(m1(i,j), m2(i,j)))
                return false;
            }
          }

          return true;
        }

        bool approximately_equal(const data__ &t1, const data__ &t2) const
        {
          data__ diff(std::abs(t1-t2));
          if (diff<=abs_tol)
          {
            return true;
          }

          data__ denom(std::max(std::abs(t1), std::abs(t2)));
          if (diff/denom<=rel_tol)
          {
            return true;
          }

          return false;
        }

        template<typename type2__>
        bool approximately_equal(const data__ &t1, const type2__ &ti2) const
        {
          data__ t2=static_cast<data__>(ti2);

          return approximately_equal(t1, t2);
        }

        template<typename type1__>
        bool approximately_equal(const type1__ &ti1, const data__ &t2) const
        {
          return approximately_equal(t2, ti1);
        }

        bool approximately_less_than(const data__ &t1, const data__ &t2) const
        {
          if (approximately_equal(t1, t2))
            return false;

          return (t1<t2);
        }

        template<typename type2__>
        bool approximately_less_than(const data__ &t1, const type2__ &ti2) const
        {
          if (approximately_equal(t1, ti2))
            return false;

          data__ t2=static_cast<data__>(ti2);
          return (t1<t2);
        }

        template<typename type1__>
        bool approximately_less_than(const type1__ &ti1, const data__ &t2) const
        {
          if (approximately_equal(ti1, t2))
            return false;

          data__ t1=static_cast<data__>(ti1);
          return (t1<t2);
        }

      private:
        data__ abs_tol;
        data__ rel_tol;
    };
  }
}
#endif
