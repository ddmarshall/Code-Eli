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

#ifndef eli_geom_point_distance
#define eli_geom_point_distance

#include <cmath>

#include "eli/code_eli.hpp"

namespace eli
{
  namespace geom
  {
    namespace point
    {
      template<typename Derived1__, typename Derived2__>
      typename Derived1__::Scalar distance2(const Eigen::MatrixBase<Derived1__> &p1, const Eigen::MatrixBase<Derived2__> &p2)
      {
        return (p1-p2).squaredNorm();
      }

      template<typename Derived1__, typename Derived2__>
      typename Derived1__::Scalar distance(const Eigen::MatrixBase<Derived1__> &p1, const Eigen::MatrixBase<Derived2__> &p2)
      {
        return (p1-p2).norm();
      }
    }
  }
}
#endif
