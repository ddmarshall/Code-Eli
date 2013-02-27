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

#ifndef eli_geom_curve_curvature_hpp
#define eli_geom_curve_curvature_hpp

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename curve__>
      void curvature(typename curve__::data_type &rho, const curve__ &c, const typename curve__::data_type &t)
      {
        // check to make sure have valid curve
        assert(c.degree()>0);

        typename curve__::point_type xp(c.fp(t)), xpp(c.fpp(t));

        if (xp.innerSize()==2)
        {
          typename curve__::data_type tmp1(std::abs(xp(0)*xpp(1)-xp(1)*xpp(0))), tmp2(xp.norm());

          rho=tmp1/(tmp2*tmp2*tmp2);
        }
        else
        {
          typename curve__::point_type tmp1;
          typename curve__::data_type tmp2(xp.norm());

          tmp1 << (xp(1)*xpp(2)-xp(2)*xpp(1)), (xp(2)*xpp(0)-xp(0)*xpp(2)), (xp(0)*xpp(1)-xp(1)*xpp(0));
          rho=tmp1.norm()/(tmp2*tmp2*tmp2);
        }
      }
    }
  }
}
#endif
