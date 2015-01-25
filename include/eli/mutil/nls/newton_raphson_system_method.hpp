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

#include "eli/code_eli.hpp"

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
        public:
          static const int hit_constraint = 101;

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
          int find_root(typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix &root, const f__ &fun, const g__ &fprime, const typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix &f0) const
          {
            typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix dx, x(x0), fx(fun(x0)), eval1, eval2, eval3;
            typename iterative_system_root_base<data__, N__, NSOL__>::jacobian_matrix fpx(fprime(x0));
            data__ abs_tol_norm, rel_tol_norm, abs_x_norm, rel_x_norm;
            typename iterative_root_base<data__>::iteration_type count;

            // calculate the function evaluated at the initial location
            eval1=fx-f0;
            abs_tol_norm=this->calculate_norm(eval1);
            eval2=(fx-f0).array()/f0.array();
            eval2.setConstant(1);
            rel_tol_norm=this->calculate_norm(eval2);
            if (this->test_converged(0, rel_tol_norm, abs_tol_norm, 0, 0))
            {
              root=x;
              return this->converged;
            }

            bool all_zero(false);
            abs_x_norm=0;
            rel_x_norm=0;
            count=0;
            while (!this->test_converged(count, rel_tol_norm, abs_tol_norm, rel_x_norm, abs_x_norm) && !all_zero)
            {
  // FIX: Don't have any easy (efficient) way of determining if matrix in invertible
  //            if (fpx==0)
  //              return iterative_root_base<data__>::no_root_found;

              bool invertible;
              bool modified = false;

              if (N__ < 4)
              {
                typename iterative_system_root_base<data__, N__, NSOL__>::jacobian_matrix inverse;
                fpx.computeInverseWithCheck(inverse, invertible);

                std::vector<bool> zerodx(N__);

                if (!invertible)
                {
                  for (size_t i=0; i<N__; ++i)
                  {
                    if (std::abs(fpx(i,i)) < std::sqrt(std::numeric_limits<data__>::epsilon()))
                    {
                      zerodx[i] = true;
                      fpx(i,i) = 1.0;
                    }
                    else
                    {
                      zerodx[i] = false;
                    }
                  }
                  modified = true;
                  fpx.computeInverseWithCheck(inverse, invertible);
                  assert(invertible);
                }

                if (invertible)
                {
                  dx = - inverse * eval1;

                  if (modified)
                  {
                    for (size_t i=0; i<N__; ++i)
                    {
                      if (zerodx[i])
                      {
                        dx(i) = 0;
                      }
                    }
                  }
                }
                else
                {
                  dx.setZero();
                }

              }
              else
              {
                dx=-fpx.lu().solve(eval1);
              }

              dx=calculate_delta_factor(x, dx);
              x+=dx;
              fx=fun(x);
              fpx=fprime(x);
              eval1=fx-f0;
              abs_tol_norm=this->calculate_norm(eval1);
              abs_x_norm=this->calculate_norm(dx);
              bool nonzero(false);
              all_zero=true;
              for (size_t i=0; i<N__; ++i)
              {
                // check if stuck and cannot move x anymore
                if (std::abs(dx(i))<=std::numeric_limits<data__>::epsilon())
                {
                  eval3(i)=std::numeric_limits<data__>::epsilon();
                }
                else
                {
                  all_zero=false;
                  eval3(i)=dx(i)/x0(i);
                }

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

              if (all_zero)
                rel_x_norm=static_cast<data__>(0);
              else
                rel_x_norm=this->calculate_norm(eval3);

              ++count;
            }

            root=x;
            if (this->max_iteration_reached(count))
              return this->max_iteration; // could not converge
            if (all_zero)
              return this->hit_constraint; // constraints limited convergence

            return this->converged;
          }

        private:
          virtual typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix
                  calculate_delta_factor(const typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix &,
                                         const typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix &dx) const
          {
            return dx;
          }

        private:
          typename iterative_system_root_base<data__, N__, NSOL__>::solution_matrix x0;
      };
    }
  }
}
#endif
