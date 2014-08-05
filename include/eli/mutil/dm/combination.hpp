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

#ifndef eli_mutil_dm_combination_hpp
#define eli_mutil_dm_combination_hpp

#include <iterator>
#include <algorithm>
#include <functional>

#include "eli/code_eli.hpp"

namespace eli
{
  namespace mutil
  {
    namespace dm
    {
      /* This code is intended to be similar to the std::next_permutation() function
       * but instead returns the next next combination of std::distance(itb, itk) elements
       * in the collection of elements [itb, ite). It is inspired by Matthieu N. answer at
       * http://stackoverflow.com/questions/127704/algorithm-to-return-all-combinations-of-k-elements-from-n
       * The iterator used has to be bidirectional iterators.
       */
      template<typename it__, typename comp__>
      bool next_combination(const it__ itb, it__ itk, const it__ ite, comp__ comp)
      {
        // catch case where either: (1) length of collection is zero, (2) desired combinations are zero,
        // (3) desired combinations are the entire list, or (4) collection is only one element long
        if ( (itb==ite) || (itb==itk) || (ite==itk) || (std::distance(itb, ite)==1) )
          return false;

        // find the next combination
        it__ it1=itk;
        it__ it2=ite; --it2;
        while (it1!=itb)
        {
          --it1;
          if ( comp(*it1, *it2) )
          {
            it__ itj=itk;
            while (!comp(*it1, *itj))
              ++itj;

            std::iter_swap(it1, itj);
            ++it1;
            ++itj;
            it2=itk;
            std::rotate(it1, itj, ite);
            while (itj != ite)
            {
              ++itj;
              ++it2;
            }
            std::rotate(itk, it2, ite);
            return true;
          }
        }
        std::rotate(itb, itk, ite);
        return false;
      }

      template <typename it__>
      bool next_combination(const it__ itb, it__ itk, const it__ ite)
      {
        return next_combination(itb, itk, ite, std::less<typename it__::value_type>());
      }

      // TODO: There really should be a prev_combination() function as well
    }
  }
}
#endif
