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

#ifdef Success  // X11 #define collides with Eigen
#undef Success
#endif

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

      template<typename data__, unsigned short dim__, typename tol__>
      data__ maximum_distance(const eli::geom::general::bounding_box<data__, dim__, tol__> &bb,
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
            len=bb.get_max()(0,i)-pt(0,i);
          }
          else if (above_max)
          {
            len=pt(0,i)-bb.get_min()(0,i);
          }
          else
          {
            data__ l1, l2;
            l1 = pt(0,i)-bb.get_min()(0,i);
            l2 = bb.get_max()(0,i)-pt(0,i);

            if(l1>=l2)
              len=l1;
            else
              len=l2;
          }
          dist2+=len*len;
        }

        return std::sqrt(dist2);
      }

      template<typename data__, unsigned short dim__, typename tol__>
      void minmax_distance(const eli::geom::general::bounding_box<data__, dim__, tol__> &bb,
                           const typename eli::geom::general::bounding_box<data__, dim__, tol__>::point_type &pt,
                           data__ &dmin,
                           data__ &dmax)
      {
        data__ maxdist2(0), mindist2(0), minlen, maxlen;
        typename eli::geom::general::bounding_box<data__, dim__, tol__>::index_type i;
        bool below_min, above_max;

        // for each dimension that is outside corresponding min/max add that to distance
        for (i=0; i<dim__; ++i)
        {
          below_min=pt(0, i) < bb.get_min()(0, i);
          above_max=pt(0, i) > bb.get_max()(0, i);

          if (below_min)
          {
            minlen=bb.get_min()(0,i)-pt(0,i);
            maxlen=bb.get_max()(0,i)-pt(0,i);
          }
          else if (above_max)
          {
            minlen=pt(0,i)-bb.get_max()(0,i);
            maxlen=pt(0,i)-bb.get_min()(0,i);
          }
          else
          {
            minlen=0;

            data__ l1, l2;
            l1 = pt(0,i)-bb.get_min()(0,i);
            l2 = bb.get_max()(0,i)-pt(0,i);

            if(l1>=l2)
              maxlen=l1;
            else
              maxlen=l2;
          }
          mindist2+=minlen*minlen;
          maxdist2+=maxlen*maxlen;
        }
        dmin=std::sqrt(mindist2);
        dmax=std::sqrt(maxdist2);
      }

    }
  }
}
#endif
