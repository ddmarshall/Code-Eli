/*********************************************************************************
* Copyright (c) 2015 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_curve_utility_hpp
#define eli_geom_curve_utility_hpp

#include <cmath>

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace utility
      {
        template<typename data__>
        void calculate_circle(data__ &rad, data__ &x0, data__ &y0, const data__ &xa, const data__ &ya,
                              const data__ &xb, const data__ &yb, const data__ &xc, const data__ &yc)
        {
          data__ denom, xamb, yamb, xcma, ycma, xbmc, ybmc, alen2, blen2, clen2;

          xamb=(xa-xb);
          yamb=(ya-yb);
          xcma=(xc-xa);
          ycma=(yc-ya);
          xbmc=(xb-xc);
          ybmc=(yb-yc);
          alen2=xa*xa+ya*ya;
          blen2=xb*xb+yb*yb;
          clen2=xc*xc+yc*yc;

          denom=2*(xc*yamb+xb*ycma+xa*ybmc);
          rad=std::sqrt((xamb*xamb+yamb*yamb)*(xcma*xcma+ycma*ycma)*(xbmc*xbmc+ybmc*ybmc))/denom;
          x0= (clen2*yamb+blen2*ycma+alen2*ybmc)/denom;
          y0=-(clen2*xamb+blen2*xcma+alen2*xbmc)/denom;
        }
      }
    }
  }
}

#endif
