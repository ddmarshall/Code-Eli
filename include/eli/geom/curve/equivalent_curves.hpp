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

#ifndef eli_geom_curve_equivalent_curves_hpp
#define eli_geom_curve_equivalent_curves_hpp

#include "eli/code_eli.hpp"

#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      bool equivalent_curves(const bezier<data__, dim__, tol__> &c0, const bezier<data__, dim__, tol__> &c1)
      {
        // if are same degree then compare control points
        if (c0.degree()==c1.degree())
        {
          return c0.approximately_equal(c1);
        }
        // else need to degree promote lower degree curve
        else if (c0.degree()<c1.degree())
        {
          bezier<data__, dim__, tol__> c0e(c0);
          while(c0e.degree()<c1.degree())
          {
            c0e.degree_promote();
          }

          return c0e.approximately_equal(c1);
        }
        else
        {
          assert(c0.degree()>c1.degree());

          bezier<data__, dim__, tol__> c1e(c1);
          while(c1e.degree()<c0.degree())
          {
            c1e.degree_promote();
          }

          return c1e.approximately_equal(c0);
        }

        return false;
      }
    }
  }
}
#endif

