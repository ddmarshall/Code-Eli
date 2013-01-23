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

#ifndef mutil_poly_root_descartes_rule_hpp
#define mutil_poly_root_descartes_rule_hpp

#include <vector>

#include "eli/poly/root/sign_changes.hpp"
#include "eli/poly/polynomial.hpp"

namespace eli
{
  namespace poly
  {
    namespace root
    {
      template<typename data__>
      int descartes_rule(const polynomial<data__> &f, bool positive)
      {
        // catch special case
        if (f.degree()==0)
          return 0;

        // get the coefficients from the polynomial
        std::vector<data__> a(f.degree()+1);
        for (size_t i=0; i<a.size(); ++i)
        {
          a[i]=f.coefficient(i);

          // change the sign of odd powers if want the negative root count
          if (!positive && i%2==1)
            a[i]*=-1;
        }

        return eli::poly::root::sign_changes(a.begin(), a.end());
      }

      template<typename data__>
      int descartes_rule(const polynomial<data__> &f)
      {
        return descartes_rule(f, true)+descartes_rule(f, false);
      }
    }
  }
}

#endif
