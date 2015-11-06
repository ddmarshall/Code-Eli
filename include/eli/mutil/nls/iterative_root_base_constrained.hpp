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

#ifndef eli_mutil_nls_iterative_root_base_constrained_hpp
#define eli_mutil_nls_iterative_root_base_constrained_hpp

#include <limits>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/nls/iterative_root_base.hpp"

namespace eli
{
  namespace mutil
  {
    namespace nls
    {
      template<typename data__>
      class iterative_root_base_constrained : public iterative_root_base<data__>
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
          iterative_root_base_constrained()
            : iterative_root_base<data_type>(), xmin(0), xmax(0), xmin_cond(IRC_NOT_USED), xmax_cond(IRC_NOT_USED)
          {
          }

          iterative_root_base_constrained(const iterative_root_base_constrained<data_type> &nrm)
            : iterative_root_base<data_type>(nrm), xmin(nrm.xmin), xmax(nrm.xmax), xmin_cond(nrm.xmin_cond), xmax_cond(nrm.xmax_cond)
          {
          }

          ~iterative_root_base_constrained()
          {
          }

          void set_periodic_condition(const data_type &dmin, const data_type &dmax)
          {
            xmin=dmin;
            xmax=dmax;
            xmin_cond=IRC_PERIODIC;
            xmax_cond=IRC_PERIODIC;
          }
          void unset_conditions()
          {
            xmin_cond=IRC_NOT_USED;
            xmax_cond=IRC_NOT_USED;
          }

          void unset_lower_condition() {xmin_cond=IRC_NOT_USED;}
          void set_lower_condition(const data_type &d, end_condition_usage ec)
          {
            if ( (xmin_cond==IRC_PERIODIC) && (ec!=IRC_PERIODIC) )
            {
              xmax_cond=IRC_NOT_USED;
            }
            if (ec==IRC_PERIODIC)
            {
              xmax_cond=IRC_PERIODIC;
            }

            xmin=d;
            xmin_cond=ec;
          }
          void get_lower_condition(data_type &d, end_condition_usage &ec)
          {
            d=xmin;
            ec=xmin_cond;
          }

          void unset_upper_condition() {xmax_cond=IRC_NOT_USED;}
          void set_upper_condition(const data_type &d, end_condition_usage ec)
          {
            if ( (xmax_cond==IRC_PERIODIC) && (ec!=IRC_PERIODIC) )
            {
              xmin_cond=IRC_NOT_USED;
            }
            if (ec==IRC_PERIODIC)
            {
              xmin_cond=IRC_PERIODIC;
            }

            xmax=d;
            xmax_cond=ec;
          }
          void get_upper_condition(data_type &d, end_condition_usage &ec)
          {
            d=xmax;
            ec=xmax_cond;
          }

        protected:
          virtual data_type calculate_delta_factor(const data_type &x, const data_type &dx) const
          {
            data_type xnew(x+dx);

            // check if min threshold is hit
            switch(xmin_cond)
            {
              case(IRC_EXCLUSIVE):
              {
                if (xnew<xmin)
                {
                  return (xmin-x);
                }
                break;
              }
              case(IRC_INCLUSIVE):
              {
                if (xnew<=xmin)
                {
                  return (xmin-x)*(1-std::numeric_limits<data_type>::epsilon());
                }
                break;
              }
              case(IRC_PERIODIC):
              {
                data_type period(xmax-xmin);

                assert(xmax>xmin);
                assert(period>0);

                if (xnew<xmin)
                {
                  xnew-=period*std::floor((xnew-xmin)/period);
                  assert(xnew>=xmin);
                  assert(xnew<=xmax);
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
            switch(xmax_cond)
            {
              case(IRC_EXCLUSIVE):
              {
                if (xnew>xmax)
                {
                  return (xmax-x);
                }
                break;
              }
              case(IRC_INCLUSIVE):
              {
                if (xnew>=xmax)
                {
                  return (xmax-x)*(1-std::numeric_limits<data_type>::epsilon());
                }
                break;
              }
              case(IRC_PERIODIC):
              {
                data_type period(xmax-xmin);

                assert(xmax>xmin);
                assert(period>0);

                if (xnew>xmax)
                {
                  xnew-=period*std::ceil((xnew-xmax)/period);
                  assert(xnew>=xmin);
                  assert(xnew<=xmax);
                }
                break;
              }
              default:
              case(IRC_NOT_USED):
              {
                break;
              }
            }

            // no threshold met
            return xnew-x;
          }

        private:
          data_type xmin, xmax;
          end_condition_usage xmin_cond, xmax_cond;
      };
    }
  }
}
#endif
