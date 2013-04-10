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

#ifndef eli_geom_intersect_minimum_distance_point_hpp
#define eli_geom_intersect_minimum_distance_point_hpp

#include "Eigen/Eigen"

#include "eli/geom/point/distance.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      template<typename Derived1__, typename Derived2__>
      typename Derived1__::Scalar minimum_distance(const Eigen::MatrixBase<Derived1__> &pt1, const Eigen::MatrixBase<Derived2__> &pt2)
      {
        return eli::geom::point::distance(pt1, pt2);
      }
    }
  }
}
#endif
