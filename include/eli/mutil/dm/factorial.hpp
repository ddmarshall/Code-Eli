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

#ifndef eli_mutil_dm_factorial_hpp
#define eli_mutil_dm_factorial_hpp

namespace eli
{
  namespace mutil
  {
    namespace dm
    {
      template<typename data__, typename natural__>
      void factorial(data__ &val, natural__ n)
      {
        val=static_cast<data__>(1);
        if (n<1)
        {
          return;
        }

        switch(n)
        {
          default:
          {
            natural__ i;
            for (i=n; i>1; --i)
            {
              val*=static_cast<data__>(i);
            }
            break;
          }
          case(10):
          {
            val=static_cast<data__>(3628800);
            break;
          }
          case(9):
          {
            val=static_cast<data__>(362880);
            break;
          }
          case(8):
          {
            val=static_cast<data__>(40320);
            break;
          }
          case(7):
          {
            val=static_cast<data__>(5040);
            break;
          }
          case(6):
          {
            val=static_cast<data__>(720);
            break;
          }
          case(5):
          {
            val=static_cast<data__>(120);
            break;
          }
          case(4):
          {
            val=static_cast<data__>(24);
            break;
          }
          case(3):
          {
            val=static_cast<data__>(6);
            break;
          }
          case(2):
          {
            val=static_cast<data__>(2);
            break;
          }
          case(1):
          {
            break;
          }
        }
      }
    }
  }
}

#endif

