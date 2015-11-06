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

#ifndef eli_mutil_nls_iterative_system_root_base_constrained_hpp
#define eli_mutil_nls_iterative_system_root_base_constrained_hpp

#include <limits>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/iterative_system_root_base.hpp"

namespace eli
{
  namespace mutil
  {
    namespace nls
    {
      template<typename data__, size_t N__, size_t NSOL__=1>
      class iterative_system_root_base_constrained : public iterative_system_root_base<data__, N__, NSOL__>
      {
        public:
          typedef data__ data_type;
          enum end_condition_usage
          {
            IRC_NOT_USED  = -1,
            IRC_EXCLUSIVE =  0,
            IRC_INCLUSIVE =  1,
            IRC_PERIODIC  =  2
          };

        public:
          iterative_system_root_base_constrained()
            : iterative_system_root_base<data_type, N__, NSOL__>()
          {
            for (size_t i=0; i<N__; ++i)
            {
              xmin[i]=0;
              xmax[i]=0;
              xmin_cond[i]=IRC_NOT_USED;
              xmax_cond[i]=IRC_NOT_USED;
            }
          }

          iterative_system_root_base_constrained(const iterative_system_root_base_constrained<data_type, N__, NSOL__> &nrm)
            : iterative_system_root_base<data_type, N__, NSOL__>(nrm)
          {
            for (size_t i=0; i<N__; ++i)
            {
              xmin[i]=nrm.xmin[i];
              xmax[i]=nrm.xmax[i];
              xmin_cond[i]=nrm.xmin_cond[i];
              xmax_cond[i]=nrm.xmax_cond[i];
            }
          }

          ~iterative_system_root_base_constrained()
          {
          }

          void unset_conditions()
          {
            for (size_t i=0; i<N__; ++i)
            {
              xmin_cond[i]=IRC_NOT_USED;
              xmax_cond[i]=IRC_NOT_USED;
            }
          }

          void set_periodic_condition(size_t i, const data_type &dmin, const data_type &dmax)
          {
            if (i<N__)
            {
              xmin[i]=dmin;
              xmax[i]=dmax;
              xmin_cond[i]=IRC_PERIODIC;
              xmax_cond[i]=IRC_PERIODIC;
            }
          }

          void unset_lower_condition()
          {
            for (size_t i=0; i<N__; ++i)
            {
              xmin_cond[i]=IRC_NOT_USED;
            }
          }
          void unset_lower_condition(size_t i)
          {
            if (i<N__)
            {
              xmin_cond[i]=IRC_NOT_USED;
            }
          }
          void set_lower_condition(size_t i, const data_type &d, end_condition_usage ec)
          {
            if (i<N__)
            {
              if ( (xmin_cond[i]==IRC_PERIODIC) && (ec!=IRC_PERIODIC) )
              {
                xmax_cond[i]=IRC_NOT_USED;
              }
              if (ec==IRC_PERIODIC)
              {
                xmax_cond[i]=IRC_PERIODIC;
              }

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
              xmax_cond[i]=IRC_NOT_USED;
            }
          }
          void unset_upper_condition(size_t i)
          {
            if (i<N__)
            {
              xmax_cond[i]=IRC_NOT_USED;
            }
          }

          void set_upper_condition(size_t i, const data_type &d, end_condition_usage ec)
          {
            if (i<N__)
            {
              if ( (xmax_cond[i]==IRC_PERIODIC) && (ec!=IRC_PERIODIC) )
              {
                xmin_cond[i]=IRC_NOT_USED;
              }
              if (ec==IRC_PERIODIC)
              {
                xmin_cond[i]=IRC_PERIODIC;
              }

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

        protected:
          virtual typename iterative_system_root_base<data_type, N__, NSOL__>::solution_matrix
                  calculate_delta_factor(const typename iterative_system_root_base<data_type, N__, NSOL__>::solution_matrix &x,
                                         const typename iterative_system_root_base<data_type, N__, NSOL__>::solution_matrix &dx) const
          {
            size_t i;
            typename iterative_system_root_base<data_type, N__, NSOL__>::solution_matrix dx_new(dx);

            // Limit constrained values to hit constraint, but allow other terms to remain.
            for (i=0; i<N__; ++i)
            {
              data_type xinew(x(i)+dx(i));

              // check if min threshold is hit
              switch(xmin_cond[i])
              {
                case(IRC_EXCLUSIVE):
                {
                  if (xinew<xmin[i])
                  {
                    dx_new(i)=xmin[i]-x(i);
                  }
                  break;
                }
                case(IRC_INCLUSIVE):
                {
                  if (xinew<=xmin[i])
                  {
                    dx_new(i)=(xmin[i]-x(i))*(1-std::numeric_limits<data_type>::epsilon());
                  }
                  break;
                }
                default:
                case(IRC_NOT_USED):
                {
                  break;
                }
              }

              // check if max threshold is hit
              switch(xmax_cond[i])
              {
                case(IRC_EXCLUSIVE):
                {
                  if (xinew>xmax[i])
                  {
                    dx_new(i)=xmax[i]-x(i);
                  }
                  break;
                }
                case(IRC_INCLUSIVE):
                {
                  if (xinew>=xmax[i])
                  {
                    dx_new(i)=(xmax[i]-x(i))*(1-std::numeric_limits<data_type>::epsilon());
                  }
                  break;
                }
                default:
                case(IRC_NOT_USED):
                {
                  break;
                }
              }
            }

            // Limit the new dx for periodic conditions.
            for (i=0; i<N__; ++i)
            {
              if (xmin_cond[i]==IRC_PERIODIC)
              {
                data_type xinew(x[i]+dx_new[i]), period(xmax[i]-xmin[i]);

                assert(xmax[i]>xmin[i]);
                assert(period>0);

                if (xinew<xmin[i])
                {
                  xinew-=period*std::floor((xinew-xmin[i])/period);

                  assert(xinew>=xmin[i]);
                  assert(xinew<=xmax[i]);

                  dx_new[i]=xinew-x[i];
                }
              }

              if (xmax_cond[i]==IRC_PERIODIC)
              {
                data_type xinew(x[i]+dx_new[i]), period(xmax[i]-xmin[i]);

                assert(xmax[i]>xmin[i]);
                assert(period>0);

                if (xinew>xmax[i])
                {
                  xinew-=period*std::ceil((xinew-xmax[i])/period);
                  dx_new[i]=fmod(xinew, period)-x[i];
                  assert(xinew>=xmin[i]);
                  assert(xinew<=xmax[i]);

                  dx_new[i]=xinew-x[i];
                }
              }
            }

            return dx_new;
          }

        private:
          data_type xmin[N__], xmax[N__];
          end_condition_usage xmin_cond[N__], xmax_cond[N__];
      };
    }
  }
}
#endif
