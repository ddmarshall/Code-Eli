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

#ifndef eli_mutil_nls_secant_method_hpp
#define eli_mutil_nls_secant_method_hpp

#include "eli/mutil/nls/iterative_root_base.hpp"

namespace eli
{
  namespace mutil
  {
    namespace nls
    {
      template<typename data__>
      class secant_method : public iterative_root_base<data__>
      {
        private:
          data__ x1, x2;

        public:
          secant_method() : iterative_root_base<data__>(), x1(0), x2(1)
          {
          }

          secant_method(const secant_method<data__> &sm)
          : iterative_root_base<data__>(sm), x1(sm.x1), x2(sm.x2)
          {
          }

          ~secant_method()
          {
          }

          void set_initial_guesses(const data__ &xg1, const data__ &xg2)
          {
            x1=xg1;
            x2=xg2;
          }

          void get_initial_guesses(data__ &xg1, data__ &xg2) const
          {
            xg1=x1;
            xg2=x2;
          }

          template<typename f__>
          int find_root(data__ &root, const f__ &fun, const data__ &f0) const
          {
            data__ x(x1), fx(fun(x1)), xm1(x2), fxm1(x2), eval, eval_abs, tmp;
            typename iterative_root_base<data__>::iteration_type count;

            // calculate the function evaluated at the initial location
            eval=fx-f0;
            eval_abs=std::abs(eval);
            if (this->test_converged(0, eval_abs/f0, eval_abs))
            {
              root=x;
              return this->converged;
            }

            count=0;
            while (!this->test_converged(count, eval_abs/f0, eval_abs))
            {
              if (fx==fxm1)
                return this->no_root_found;

              tmp=x;
              x-=eval*(x-xm1)/(fx-fxm1);
              xm1=tmp;
              fxm1=fx;
              fx=fun(x);
              eval=fx-f0;
              eval_abs=std::abs(eval);

              ++count;
            }

            root=x;
            if (this->max_iteration_reached(count))
              return this->max_iteration; // could not converge

            return this->converged;
          }
      };
    }
  }
}
#endif
