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

#ifndef eli_mutil_nls_newton_raphson_shacham_system_method_hpp
#define eli_mutil_nls_newton_raphson_shacham_system_method_hpp

#include <limits>

#include "eli/mutil/nls/newton_raphson_system_method.hpp"

namespace eli
{
  namespace mutil
  {
    namespace nls
    {
      template<typename data__, size_t N__, size_t NSOL__=1>
      class newton_raphson_shacham_system_method : public newton_raphson_system_method<data__, N__, NSOL__>
      {
        public:
          typedef data__ data_type;
          enum end_condition_usage
          {
            NRS_NOT_USED  = -1,
            NRS_EXCLUSIVE =  0,
            NRS_INCLUSIVE =  1
          };

        public:
          newton_raphson_shacham_system_method()
            : newton_raphson_system_method<data_type, N__, NSOL__>()
          {
            for (size_t i=0; i<N__; ++i)
            {
              xmin[i]=0;
              xmax[i]=0;
              xmin_cond[i]=NRS_NOT_USED;
              xmax_cond[i]=NRS_NOT_USED;
            }
          }

          newton_raphson_shacham_system_method(const newton_raphson_shacham_system_method<data_type, N__, NSOL__> &nrm)
            : newton_raphson_system_method<data_type, N__, NSOL__>(nrm)
          {
            for (size_t i=0; i<N__; ++i)
            {
              xmin[i]=nrm.xmin[i];
              xmax[i]=nrm.xmax[i];
              xmin_cond[i]=nrm.xmin_cond[i];
              xmax_cond[i]=nrm.xmax_cond[i];
            }
          }

          ~newton_raphson_shacham_system_method()
          {
          }

          void unset_lower_condition()
          {
            for (size_t i=0; i<N__; ++i)
            {
              xmin_cond[i]=NRS_NOT_USED;
            }
          }
          void unset_lower_condition(size_t i)
          {
            if (i<N__)
            {
              xmin_cond[i]=NRS_NOT_USED;
            }
          }
          void set_lower_condition(size_t i, const data_type &d, end_condition_usage ec)
          {
            if (i<N__)
            {
              xmin[i]=d;
              xmin_cond[i]=ec;
            }
          }
          void get_lower_condition(data_type &d, end_condition_usage &ec, size_t i)
          {
            if (i<N__)
            {
            d=xmin[i];
            ec=xmin_cond[i];
            }
          }

          void unset_upper_condition()
          {
            for (size_t i=0; i<N__; ++i)
            {
              xmax_cond[i]=NRS_NOT_USED;
            }
          }
          void unset_upper_condition(size_t i)
          {
            if (i<N__)
            {
              xmax_cond[i]=NRS_NOT_USED;
            }
          }

          void set_upper_condition(size_t i, const data_type &d, end_condition_usage ec)
          {
            if (i<N__)
            {
              xmax[i]=d;
              xmax_cond[i]=ec;
            }
          }
          void get_upper_condition(data_type &d, end_condition_usage &ec, size_t i)
          {
            if (i<N__)
            {
              d=xmax[i];
              ec=xmax_cond[i];
            }
          }

        private:
          virtual typename iterative_system_root_base<data_type, N__, NSOL__>::solution_matrix
                  calculate_delta_factor(const typename iterative_system_root_base<data_type, N__, NSOL__>::solution_matrix &x,
                                         const typename iterative_system_root_base<data_type, N__, NSOL__>::solution_matrix &dx) const
          {
            data_type lambda(1.0);

            // calculate the minimum of all limits
            for (size_t i=0; i<N__; ++i)
            {
              data_type xinew(x(i)+lambda*dx(i));

              // check if min threshold is hit
              switch(xmin_cond[i])
              {
                case(NRS_EXCLUSIVE):
                {
                  if (xinew<xmin[i])
                  {
                    lambda=(xmin[i]-x(i))/dx(i);
                  }
                  break;
                }
                case(NRS_INCLUSIVE):
                {
                  if (xinew<=xmin[i])
                  {
                    lambda=(xmin[i]-x(i))*(1-std::numeric_limits<data_type>::epsilon())/dx(i);
                  }
                  break;
                }
                default:
                case(NRS_NOT_USED):
                {
                  break;
                }
              }

              // check if min threshold is hit
              switch(xmin_cond[i])
              {
                case(NRS_EXCLUSIVE):
                {
                  if (xinew>xmax[i])
                  {
                    lambda=(xmax[i]-x(i))/dx(i);
                  }
                  break;
                }
                case(NRS_INCLUSIVE):
                {
                  if (xinew>=xmax[i])
                  {
                    lambda=(xmax[i]-x(i))*(1-std::numeric_limits<data_type>::epsilon())/dx(i);
                  }
                  break;
                }
                default:
                case(NRS_NOT_USED):
                {
                  break;
                }
              }
            }

            return (lambda*dx).eval();
          }

        private:
          data_type xmin[N__], xmax[N__];
          end_condition_usage xmin_cond[N__], xmax_cond[N__];
      };
    }
  }
}
#endif
