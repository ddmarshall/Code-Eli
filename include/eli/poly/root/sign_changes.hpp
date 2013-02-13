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

#ifndef eli_poly_root_sign_changes_hpp
#define eli_poly_root_sign_changes_hpp

namespace eli
{
  namespace poly
  {
    namespace root
    {
      // this counts the number of sign changes in collection (excluding zeros)
      template<typename it__>
      int sign_changes(const it__ itb, const it__ ite)
      {
        if (itb==ite)
          return 0;

        int count(0);
        it__ it(itb);
        bool prev_sign=((*it)>0);

        for (++it; it!=ite; ++it)
        {
          if ( ((*it)!=0) && (prev_sign!=((*it)>0)) )
          {
            count++;
            prev_sign=!prev_sign;
          }
        }

        return count;
      }
    }
  }
}

#endif
