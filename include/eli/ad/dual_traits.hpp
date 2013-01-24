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

#ifndef eli_ad_dual_traits_hpp
#define eli_ad_dual_traits_hpp

#include "eli/ad/traits.hpp"

namespace eli
{
  namespace ad
  {
    // forward declare constant value type
    namespace dual_number
    {
      template <typename data__> class constant;
    }

    template <typename data__>
    class traits<dual_number::constant<data__> >
    {
      public:
        // how to refer to a constant reference
        typedef dual_number::constant<data__> const_expr_ref;

        // how to refer to a reference
        typedef dual_number::constant<data__> expr_ref;
    };

  }
}

#endif
