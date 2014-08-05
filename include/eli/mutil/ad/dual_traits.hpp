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

#ifndef eli_mutil_ad_dual_traits_hpp
#define eli_mutil_ad_dual_traits_hpp

#include "eli/code_eli.hpp"

#include "eli/util/traits.hpp"

namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      // forward declare constant value type
      namespace dual_number
      {
        template <typename data__> class constant;
      }
    }
  }
}

namespace eli
{
  namespace util
  {
    template <typename data__>
    class traits<eli::mutil::ad::dual_number::constant<data__> >
    {
      public:
        // how to refer to a constant reference
        typedef eli::mutil::ad::dual_number::constant<data__> const_expr_ref;

        // how to refer to a reference
        typedef eli::mutil::ad::dual_number::constant<data__> expr_ref;
    };
  }
}

#endif
