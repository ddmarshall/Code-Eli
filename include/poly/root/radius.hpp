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

#ifndef eli_poly_root_radius_hpp
#define eli_poly_root_radius_hpp

#include <cmath>

#include "poly/polynomial.hpp"

namespace eli
{
  namespace poly
  {
    namespace root
    {
      template<typename data__>
      data__ max_radius(const polynomial<data__> &f)
      {
        data__ term, max_term(0), radius(1);
        typename polynomial<data__>::index_type i, deg(f.degree());
        typename polynomial<data__>::coefficient_type coef;

        // get the coefficient info
        f.get_coefficients(coef);

        // cycle through coefficients
        for (i=0; i<deg; ++i)
        {
          term=std::abs(coef(i)/coef(deg));
          if (term>max_term)
            max_term=term;
        }
        radius+=max_term;

        return radius;
      }

      template<typename data__>
      data__ at_least_radius(const polynomial<data__> &f)
      {
        data__ term1, term2, radius(1);
        typename polynomial<data__>::index_type i, deg(f.degree());
        typename polynomial<data__>::coefficient_type coef;

        // get the coefficient info
        f.get_coefficients(coef);

        // calculate the two terms
        term2=std::pow(std::abs(coef[0]/coef[deg]), 1/static_cast<data__>(deg));
        if (coef[1]==0)
          return term2;
        term1=deg*std::abs(coef[0]/coef[1]);
        return std::min(term1, term2);
      }
    }
  }
}

#endif
