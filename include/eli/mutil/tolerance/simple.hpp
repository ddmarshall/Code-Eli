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

#ifndef eli_geom_tolerance_simple_hpp
#define eli_geom_tolerance_simple_hpp

#include <Eigen/Eigen>  // Eigen::MatrixBase
#include <cmath>  // std::cmath
#include <limits> // std::numeric_limits

namespace eli
{
  namespace mutil
  {
    namespace tolerance
    {
      template<typename data__>
      class simple
      {
        private:
          data__ rel_tol;
          data__ abs_tol;
        public:
          simple() : rel_tol(std::sqrt(std::numeric_limits<data__>::epsilon())),
                     abs_tol(std::numeric_limits<data__>::epsilon()) {}
          simple(const data__ &rel, const data__ &abs) : rel_tol(rel), abs_tol(abs) {}
          ~simple() {}

          template<typename Derived>
          bool operator()(const Eigen::MatrixBase<Derived> &t1, const Eigen::MatrixBase<Derived> &t2)
          {
            data__ diff((t1-t2).norm());
            if (diff<abs_tol)
              return true;

            data__ denom(std::max(t1.norm(), t2.norm()));
            if (diff/denom<rel_tol)
              return true;

            return false;
          }

          template<typename type2__>
          bool operator()(const data__ &t1, const type2__ &t2)
          {
            data__ diff(std::abs(t1-t2));
            if (diff<abs_tol)
              return true;

            data__ denom(std::max(std::abs(t1), std::abs(t2)));
            if (diff/denom<rel_tol)
              return true;

            return false;
          }

          template<typename type1__>
          bool operator()(const type1__ &t1, const data__ &t2)
          {
            data__ diff(std::abs(t1-t2));
            if (diff<abs_tol)
              return true;

            data__ denom(std::max(std::abs(t1), std::abs(t2)));
            if (diff/denom<rel_tol)
              return true;

            return false;
          }

          bool operator()(const data__ &t1, const data__ &t2)
          {
            data__ diff(std::abs(t1-t2));
            if (diff<abs_tol)
              return true;

            data__ denom(std::max(std::abs(t1), std::abs(t2)));
            if (diff/denom<rel_tol)
              return true;

            return false;
          }
      };
    }
  }
}
#endif
