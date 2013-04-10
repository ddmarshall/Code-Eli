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

#ifndef eli_geom_intersect_minimum_distance_bounding_box_hpp
#define eli_geom_intersect_minimum_distance_bounding_box_hpp

#include "Eigen/Eigen"

#include "eli/geom/point/distance.hpp"
#include "eli/geom/general/bounding_box.hpp"


namespace eli
{
  namespace geom
  {
    namespace intersect
    {
      template<typename data__, unsigned short dim__, typename tol__>
      data__ minimum_distance(const eli::geom::general::bounding_box<data__, dim__, tol__> &bb,
                              const typename eli::geom::general::bounding_box<data__, dim__, tol__>::point_type &pt)
      {
        data__ dist2(0), len;
        typename eli::geom::general::bounding_box<data__, dim__, tol__>::index_type i;
        bool below_min, above_max;

        // for each dimension that is outside corresponding min/max add that to distance
        for (i=0; i<dim__; ++i)
        {
          below_min=pt(0, i) < bb.get_min()(0, i);
          above_max=pt(0, i) > bb.get_max()(0, i);

          if (below_min)
          {
            len=bb.get_min()(0,i)-pt(0,i);
          }
          else if (above_max)
          {
            len=pt(0,i)-bb.get_max()(0,i);
          }
          else
          {
            len=0;
          }
          dist2+=len*len;
        }

        return std::sqrt(dist2);
      }
    }
  }
}
#endif
