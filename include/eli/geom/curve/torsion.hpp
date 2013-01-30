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

#ifndef eli_geom_curve_torsion_hpp
#define eli_geom_curve_torsion_hpp

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      // FIX: (1) THIS NEEDS A UNIT TEST
      template<typename curve__>
      void torsion(typename curve__::data_type &tau, const curve__ &c, const typename curve__::data_type &t)
      {
        // check to make sure have valid curve
        assert(c.degree()>0);

        // check to make sure given valid parametric value
        assert((t>=0) && (t<=1));

        typename curve__::point_type xp(c.fp(t)), xpp(c.fpp(t)), xppp(c.fppp(t)), tmp1;
        typename curve__::data_type tmp2;

        tmp1=xp.cross(xpp);
        tmp2=tmp1.norm();
        tau=tmp1.dot(xppp)/(tmp2*tmp2);
      }
    }
  }
}
#endif
