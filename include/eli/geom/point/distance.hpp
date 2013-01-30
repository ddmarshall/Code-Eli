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

#include <Eigen/Eigen>

namespace eli
{
  namespace geom
  {
    namespace point
    {
      template<typename data__, typename point_type1__, typename point_type2__>
      void distance2(data__ &d, const point_type1__ &p1, const point_type2__ &p2)
      {
        Eigen::Matrix<typename point_type1__::Scalar, 1, Eigen::Dynamic> dif;

        dif=p1-p2;
        d=dif.dot(dif);
      }

      template<typename data__, typename point_type1__, typename point_type2__>
      void distance(data__ &d, const point_type1__ &p1, const point_type2__ &p2)
      {
        distance2(d, p1, p2);
        d=std::sqrt(d);
      }
    }
  }
}
#endif
