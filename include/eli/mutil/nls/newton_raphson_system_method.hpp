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

#ifndef eli_mutil_nls_newton_raphson_system_method_hpp
#define eli_mutil_nls_newton_raphson_system_method_hpp

#include "eli/mutil/nls/iterative_system_root_base.hpp"

namespace eli
{
  namespace mutil
  {
    namespace nls
    {
      template<typename data__, size_t N__, size_t NSOL__=1>
      class newton_raphson_system_method : public iterative_system_root_base<data__, N__, NSOL__>
      {
        private:
          typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix x0;

        public:
          newton_raphson_system_method()
          : iterative_system_root_base<data__, N__, NSOL__>()
          {
            x0.setConstant(static_cast<data__>(0));
          }

          newton_raphson_system_method(const newton_raphson_system_method<data__, N__, NSOL__> &nrm)
          : iterative_system_root_base<data__, N__, NSOL__>(nrm), x0(nrm.x0)
          {
          }

          ~newton_raphson_system_method()
          {
          }

          void set_initial_guess(const typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix &xg)
          {
            x0=xg;
          }

          const typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix & get_initial_guess() const
          {
            return x0;
          }

          template<typename f__, typename g__>
          typename iterative_root_base<data__>::status find_root(typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix &root, const f__ &fun, const g__ &fprime, const typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix &f0) const
          {
            typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix dx, x(x0), fx(fun(x0)), eval1, eval2;
            typename iterative_system_root_base<data__, N__, NSOL__>::jacobian_matrix fpx(fprime(x0));
            data__ abs_tol_norm, rel_tol_norm;
            typename iterative_root_base<data__>::iteration_type count;

            // calculate the function evaluated at the initial location
            eval1=fx-f0;
            abs_tol_norm=this->calculate_norm(eval1);
            eval2=(fx-f0).array()/f0.array();
            eval2.setConstant(1);
            rel_tol_norm=this->calculate_norm(eval2);
            if (this->test_converged(0, rel_tol_norm, abs_tol_norm))
            {
              root=x;
              return iterative_root_base<data__>::converged;
            }

            count=0;
            while (!this->test_converged(count, rel_tol_norm, abs_tol_norm))
            {
  // Don't have any easy (efficient) way of determining if matrix in invertible
  //            if (fpx==0)
  //              return iterative_root_base<data__>::no_root_found;

              dx=fpx.lu().solve(eval1);
              x-=dx;
              fx=fun(x);
              fpx=fprime(x);
              eval1=fx-f0;
              abs_tol_norm=this->calculate_norm(eval1);
              bool nonzero(false);
              for (size_t i=0; i<N__; ++i)
              {
                if (std::abs(f0(i))<=std::numeric_limits<data__>::epsilon())
                  eval2(i)=std::numeric_limits<data__>::epsilon();
                else
                {
                  nonzero=true;
                  eval2(i)=eval1(i)/f0(i);
                }
              }
              if (nonzero)
                rel_tol_norm=this->calculate_norm(eval2);
              else
                rel_tol_norm=static_cast<data__>(0);

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
