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

#ifndef eli_geom_intersect_minimum_distance_line_hpp
#define eli_geom_intersect_minimum_distance_line_hpp

#include <cmath>

#ifdef Success  // X11 #define collides with Eigen
#undef Success
#endif

#include "Eigen/Eigen"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      template<typename Derived1__, typename Derived2__, typename Derived3__>
      typename Derived1__::Scalar minimum_distance(typename Derived1__::Scalar &t, const Eigen::MatrixBase<Derived1__> &a0,
                                                   const Eigen::MatrixBase<Derived2__> &a1, const Eigen::MatrixBase<Derived3__> &pt)
      {
        Eigen::Matrix<typename Derived1__::Scalar, 1, Eigen::Dynamic> pma0;
        typename Derived1__::Scalar a1a1;

        a1a1=a1.dot(a1);
        pma0=pt-a0;

        if (a1a1==0)
        {
          t=0;
          return eli::geom::point::distance(a0, pt);
        }

        t=pma0.dot(a1)/a1a1;
        return std::sqrt(std::max(static_cast<typename Derived1__::Scalar>(0), pma0.dot(pma0)-a1a1*t*t));
      }
    }
  }
}
#endif
