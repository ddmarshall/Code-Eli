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

#ifndef eli_geom_general_continuity_hpp
#define eli_geom_general_continuity_hpp

namespace eli
{
  namespace geom
  {
    namespace general
    {
      enum continuity
      {
        NOT_CONNECTED = -1,
        C0            = 0,
        C1            = 1,
        C2            = 2,
        G1            = 101,
        G2            = 102
      };
    }
  }
}
#endif
