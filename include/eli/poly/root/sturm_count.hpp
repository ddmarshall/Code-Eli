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

#ifndef eli_poly_root_sturm_count_hpp
#define eli_poly_root_sturm_count_hpp


#include <cmath>
#include <vector>

#include "eli/poly/root/sign_changes.hpp"
#include "eli/poly/polynomial.hpp"
#include "eli/poly/root/radius.hpp"

namespace eli
{
  namespace poly
  {
    namespace root
    {
      template<typename data__>
      void sturm_functions(std::vector< polynomial<data__> > &sturm_fun, const polynomial<data__> &f)
      {
        size_t i, n;

        // resize the Sturm function vector
        n=f.degree()+1;
        sturm_fun.resize(n);

        polynomial<data__> *ptemp(f.fp()), pjunk;

        // initialize the Sturm functions
        sturm_fun[0]=f;
        ptemp=f.fp();
        sturm_fun[1]=*ptemp;
        delete ptemp;

        // build the Sturm functions
        sturm_fun[0].divide(std::abs(sturm_fun[0].coefficient(sturm_fun[0].degree())));
        sturm_fun[1].divide(std::abs(sturm_fun[1].coefficient(sturm_fun[1].degree())));
        for (i=2; i<n; ++i)
        {
          if (sturm_fun[i-1].degree()>0)
          {
            data__ small_no, max_c(0);
            typename polynomial<data__>::coefficient_type c;

            sturm_fun[i-1].get_coefficients(c);
            for (typename polynomial<data__>::index_type j=0; j<=sturm_fun[i-1].degree(); ++j)
            {
              if (std::abs(c(j))>max_c)
                max_c=std::abs(c(j));
            }
            small_no=max_c*std::sqrt(std::numeric_limits<data__>::epsilon());

            pjunk.divide(sturm_fun[i], sturm_fun[i-2], sturm_fun[i-1]);
            sturm_fun[i].negative();
            if (std::abs(sturm_fun[i].coefficient(sturm_fun[i].degree()))>small_no)
              sturm_fun[i].divide(std::abs(sturm_fun[i].coefficient(sturm_fun[i].degree())));
            sturm_fun[i].adjust_zero(small_no);
          }
          else
          {
            sturm_fun[i]=0;
          }
        }
      }

      template<typename data__, typename data2__>
      int sturm_count(const std::vector< polynomial<data__> > &sturm_fun, const data2__ &xmin, const data2__ &xmax)
      {
        std::vector<data__> eval(sturm_fun.size());
        int count_min, count_max;

        // evaluate functions at end points and count sign changes
        for (size_t i=0; i<eval.size(); ++i)
          eval[i]=sturm_fun[i].f(xmin);
        count_min = eli::poly::root::sign_changes(eval.begin(), eval.end());
        for (size_t i=0; i<eval.size(); ++i)
          eval[i]=sturm_fun[i].f(xmax);
        count_max = eli::poly::root::sign_changes(eval.begin(), eval.end());

        // return the difference between sign changes for min value and sign changes for max value
        return count_min-count_max;
      }

      template<typename data__, typename data2__>
      int sturm_count(const polynomial<data__> &f, const data2__ &xmin, const data2__ &xmax)
      {
        // short circuit degenerate cases
        if (xmax<=xmin)
          return 0;

        // catch another degenerate case
        if (f.degree()==0)
          return 0;

        std::vector< polynomial<data__> > sturm_fun;

        // calculate the Sturm functions
        sturm_functions(sturm_fun, f);

        return sturm_count(sturm_fun, xmin, xmax);
      }

      template<typename data__>
      int sturm_count(const polynomial<data__> &f, bool positive)
      {
        typename polynomial<data__>::data_type xmin, xmax;

        // find appropriate xmin and xmax
        if (positive)
        {
          xmin=static_cast<typename polynomial<data__>::data_type>(0);
          xmax=eli::poly::root::max_radius(f);
        }
        else
        {
          xmin=-eli::poly::root::max_radius(f);
          xmax=static_cast<typename polynomial<data__>::data_type>(0);
        }

        return sturm_count(f, xmin, xmax);
      }

      template<typename data__>
      int sturm_count(const std::vector< polynomial<data__> > &sturm_fun, bool positive)
      {
        typename polynomial<data__>::data_type xmin, xmax;

        // find appropriate xmin and xmax
        if (positive)
        {
          xmin=static_cast<typename polynomial<data__>::data_type>(0);
          xmax=eli::poly::root::max_radius(sturm_fun[0]);
        }
        else
        {
          xmin=-eli::poly::root::max_radius(sturm_fun[0]);
          xmax=static_cast<typename polynomial<data__>::data_type>(0);
        }

        return sturm_count(sturm_fun, xmin, xmax);
      }

      template<typename data__>
      int sturm_count(const polynomial<data__> &f)
      {
        typename polynomial<data__>::data_type xmin, xmax;

        // find appropriate xmin and xmax
        xmax=eli::poly::root::max_radius(f);
        xmin=-xmax;

        return sturm_count(f, xmin, xmax);
      }

      template<typename data__>
      int sturm_count(const std::vector< polynomial<data__> > &sturm_fun)
      {
        typename polynomial<data__>::data_type xmin, xmax;

        // find appropriate xmin and xmax
        xmax=eli::poly::root::max_radius(sturm_fun[0]);
        xmin=-xmax;

        return sturm_count(sturm_fun, xmin, xmax);
      }
    }
  }
}

#endif
