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

#ifndef eli_mutil_dm_binomial_coef_hpp
#define eli_mutil_dm_binomial_coef_hpp

#include "eli/code_eli.hpp"

namespace eli
{
  namespace mutil
  {
    namespace dm
    {
      template<typename data__, typename natural__>
      void n_choose_k(data__ &v, natural__ n, natural__ k)
      {
        natural__ i;

        if (n<k)
        {
          assert(false);
          v = -1;
          return;
        }
        if (n==k)
        {
          v = 1;
          return;
        }

        // use symmetry property to get smaller term
        if (k>n-k)
          k=n-k;

        v=1;
        for (i=1; i<=k; ++i)
        {
          // try and keep from overflowing integral data types
          if (static_cast<data__>(i)>v)
          {
            v*=n-k+i;
            v/=i;
          }
          else
          {
            v/=i;
            v*=n-k+i;
          }
        }
      }

      template<typename data__, typename natural__>
      void binomial_coefficient(data__ &v, const data__ &n, natural__ k)
      {
        natural__ i;

        if (n==k)
        {
          v = 1;
          return;
        }

        v=1;
        for (i=1; i<=k; ++i)
        {
          v*=n-k+i;
          v/=i;
        }
      }

      template<typename data__, typename natural__>
      void binomial_coefficient(data__ &v, natural__ n, natural__ k)
      {
        binomial_coefficient(v, static_cast<data__>(n), k);
      }
    }
  }
}

#endif
