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

#ifndef eli_geom_intersect_minimum_distance_plane_hpp
#define eli_geom_intersect_minimum_distance_plane_hpp

#include <cmath>

#include "Eigen/Eigen"

#ifdef Success  // X11 #define collides with Eigen
#undef Success
#endif

#include "eli/geom/intersect/minimum_distance_line.hpp"

namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      template<typename Derived1__, typename Derived2__, typename Derived3__, typename Derived4__>
      typename Derived1__::Scalar minimum_distance(typename Derived1__::Scalar &u, typename Derived1__::Scalar &v,
                                                   const Eigen::MatrixBase<Derived1__> &a, const Eigen::MatrixBase<Derived2__> &b,
                                                   const Eigen::MatrixBase<Derived3__> &c, const Eigen::MatrixBase<Derived4__> &pt)
      {
        Eigen::Matrix<typename Derived1__::Scalar, 1, Eigen::Dynamic> pmc;
        typename Derived1__::Scalar aa, ab, bb, denom;

        aa=a.dot(a);
        ab=a.dot(b);
        bb=b.dot(b);
        pmc=pt-c;

        // degenerate point case
        if ((aa==0) && (bb==0))
        {
          u=0;
          v=0;
          return eli::geom::point::distance(c, pt);
        }

        // degenerate line case
        if (ab*ab==aa*bb)
        {
          v=0;
          return eli::geom::intersect::minimum_distance(u, c, a, pt);
        }

        // surface & point case
        denom=ab*ab-aa*bb;
        u=(ab*b-bb*a).dot(pmc)/denom;
        v=(ab*a-aa*b).dot(pmc)/denom;
        return (a*u+b*v-pmc).norm();
      }
    }
  }
}
#endif
