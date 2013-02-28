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

#ifndef eli_mutil_nls_newton_raphson_method_hpp
#define eli_mutil_nls_newton_raphson_method_hpp

#include "eli/mutil/nls/iterative_root_base.hpp"

namespace eli
{
  namespace mutil
  {
    namespace nls
    {
      template<typename data__>
      class newton_raphson_method : public iterative_root_base<data__>
      {
        private:
          data__ x0;

        public:
          newton_raphson_method() : iterative_root_base<data__>(), x0(0)
          {
          }

          newton_raphson_method(const newton_raphson_method<data__> &nrm)
          : iterative_root_base<data__>(nrm), x0(nrm.x0)
          {
          }

          ~newton_raphson_method()
          {
          }

          void set_initial_guess(const data__ &xg)
          {
            x0=xg;
          }

          const data__ & get_initial_guess() const
          {
            return x0;
          }

          template<typename f__, typename g__>
          typename iterative_root_base<data__>::status find_root(data__ &root, const f__ &fun, const g__ &fprime, const data__ &f0) const
          {
            data__ x(x0), fx(fun(x0)), fpx(fprime(x0)), eval, eval_abs;
            typename iterative_root_base<data__>::iteration_type count;

            // calculate the function evaluated at the initial location
            eval=fx-f0;
            eval_abs=std::abs(eval);
            if (this->test_converged(0, eval_abs/f0, eval_abs))
            {
              root=x;
              return iterative_root_base<data__>::converged;
            }

            count=0;
            while (!this->test_converged(count, eval_abs/f0, eval_abs))
            {
              if (fpx==0)
                return iterative_root_base<data__>::no_root_found;

              x-=eval/fpx;
              fx=fun(x);
              fpx=fprime(x);
              eval=fx-f0;
              eval_abs=std::abs(eval);

              ++count;
            }

            root=x;
            if (this->max_iteration_reached(count))
              return iterative_root_base<data__>::max_iteration; // could not converge

            return iterative_root_base<data__>::converged;
          }
      };
    }
  }
}
#endif
