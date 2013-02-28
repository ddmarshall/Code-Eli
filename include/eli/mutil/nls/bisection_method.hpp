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

#ifndef eli_mutil_nls_bisection_method_hpp
#define eli_mutil_nls_bisection_method_hpp

#include <cmath>

#include "eli/mutil/nls/iterative_root_base.hpp"

namespace eli
{
  namespace mutil
  {
    namespace nls
    {
      template<typename data__>
      class bisection_method : public iterative_root_base<data__>
      {
        private:
          data__ xmin, xmax;

        public:
          bisection_method() : iterative_root_base<data__>(), xmin(0), xmax(1)
          {
          }

          bisection_method(const bisection_method<data__> &bm)
            : iterative_root_base<data__>(bm), xmin(bm.xmin), xmax(bm.xmax)
          {
          }

          ~bisection_method()
          {
          }

          void set_low_bound(const data__ &xl)
          {
            xmin=xl;
          }

          const data__ & get_low_bound() const
          {
            return xmin;
          }

          void set_high_bound(const data__ &xh)
          {
            xmax=xh;
          }

          const data__ & get_high_bound() const
          {
            return xmax;
          }

          void set_bounds(const data__ &xl, const data__ &xh)
          {
            if (xh < xl)
            {
              assert(false);
              return;
            }

            set_low_bound(xl);
            set_high_bound(xh);
          }

          template<typename f__>
          typename iterative_root_base<data__>::status find_root(data__ &root, const f__ &fun, const data__ &f0) const
          {
            data__ xmn(xmin), xmx(xmax), two(2);
            data__ fmin, fmid, fmax, xmid, data_abs, data_abs2;
            typename iterative_root_base<data__>::iteration_type count;

            // calculate the function evaluated at the minimum location
            fmin=fun(xmn)-f0;
            data_abs=std::abs(fmin);
            if (this->test_converged(0, data_abs/f0, data_abs))
            {
              root=xmn;
              return iterative_root_base<data__>::converged;
            }

            // calculate the function evaluated at the maximum location
            fmax=fun(xmx)-f0;
            data_abs=std::abs(fmax);
            if (this->test_converged(0, data_abs/f0, data_abs))
            {
              root=xmx;
              return iterative_root_base<data__>::converged;
            }

            count=0;
            xmid=(xmx+xmn)/two;
            fmid=fun(xmid)-f0;
            data_abs=std::abs(xmx-xmn);
            data_abs2=std::abs(fmid);
            // this tests how close the min and max roots are to each other and tests how close the function evaluations are to zero
            while  (!this->test_converged(count, data_abs/xmid, data_abs) && !this->test_converged(count, data_abs2/f0, data_abs2))
            {
              // if middle point is opposite side of root as maximum point
              if (fmid*fmax<0)
              {
                fmin=fmid;
                xmn=xmid;
              }
              else if (fmid*fmin<0)
              {
                fmax=fmid;
                xmx=xmid;
              }
              else
              {
  #if (DEBUG>1)
                std::cerr << "bisection_method_def.hpp:" << " Iterative root not found! Bounds are [" << xmn << ", " << xmx << "]=>[" << fun(xmn) << ", " << fun(xmx) << "] with f0=" << f0 << std::endl;
  #endif
                root=(xmn+xmx)/2;
                return iterative_root_base<data__>::no_root_found;
              }

              xmid=(xmx+xmn)/two;
              fmid=fun(xmid)-f0;
              data_abs=std::abs(xmx-xmn);
              data_abs2=std::abs(fmid);
              ++count;
            }

            root=xmid;
            if (this->max_iteration_reached(count))
            {
  #if (DEBUG>1)
              std::cerr << "bisection_method_def.hpp:" << " Iterative root not converged! Bounds are [" << xmn << ", " << xmx << "]=>[" << fun(xmn) << ", " << fun(xmx) << "] with f0=" << f0 << std::endl;
  #endif
              return iterative_root_base<data__>::max_iteration; // could not converge
            }

            return iterative_root_base<data__>::converged;
          }
      };
    }
  }
}
#endif
