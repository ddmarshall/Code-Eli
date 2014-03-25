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

#ifndef eli_mutil_nls_newton_raphson_constrained_method_hpp
#define eli_mutil_nls_newton_raphson_constrained_method_hpp

#include <limits>
#include <algorithm>

#include "eli/mutil/nls/newton_raphson_method.hpp"

namespace eli
{
  namespace mutil
  {
    namespace nls
    {
      template<typename data__>
      class newton_raphson_constrained_method : public newton_raphson_method<data__>
      {
        public:
          typedef data__ data_type;

          enum end_condition_usage
          {
            NRC_NOT_USED  = -1,
            NRC_EXCLUSIVE =  0,
            NRC_INCLUSIVE =  1,
            NRC_PERIODIC  =  2
          };

        public:
          newton_raphson_constrained_method()
            : newton_raphson_method<data_type>(), xmin(0), xmax(0), xmin_cond(NRC_NOT_USED), xmax_cond(NRC_NOT_USED)
          {
          }

          newton_raphson_constrained_method(const newton_raphson_constrained_method<data_type> &nrm)
            : newton_raphson_method<data_type>(nrm), xmin(nrm.xmin), xmax(nrm.xmax), xmin_cond(nrm.xmin_cond), xmax_cond(nrm.xmax_cond)
          {
          }

          ~newton_raphson_constrained_method()
          {
          }

          void set_periodic_condition(const data_type &dmin, const data_type &dmax)
          {
            xmin=dmin;
            xmax=dmax;
            xmin_cond=NRC_PERIODIC;
            xmax_cond=NRC_PERIODIC;
          }
          void unset_conditions()
          {
            xmin_cond=NRC_NOT_USED;
            xmax_cond=NRC_NOT_USED;
          }

          void unset_lower_condition() {xmin_cond=NRC_NOT_USED;}
          void set_lower_condition(const data_type &d, end_condition_usage ec)
          {
            if ( (xmin_cond==NRC_PERIODIC) && (ec!=NRC_PERIODIC) )
            {
              xmax_cond=NRC_NOT_USED;
            }
            if (ec==NRC_PERIODIC)
            {
              xmax_cond=NRC_PERIODIC;
            }

            xmin=d;
            xmin_cond=ec;
          }
          void get_lower_condition(data_type &d, end_condition_usage &ec)
          {
            d=xmin;
            ec=xmin_cond;
          }

          void unset_upper_condition() {xmax_cond=NRC_NOT_USED;}
          void set_upper_condition(const data_type &d, end_condition_usage ec)
          {
            if ( (xmax_cond==NRC_PERIODIC) && (ec!=NRC_PERIODIC) )
            {
              xmin_cond=NRC_NOT_USED;
            }
            if (ec==NRC_PERIODIC)
            {
              xmin_cond=NRC_PERIODIC;
            }

            xmax=d;
            xmax_cond=ec;
          }
          void get_upper_condition(data_type &d, end_condition_usage &ec)
          {
            d=xmax;
            ec=xmax_cond;
          }

        private:
          virtual data_type calculate_delta_factor(const data_type &x, const data_type &dx) const
          {
            data_type xnew(x+dx);

            // check if min threshold is hit
            switch(xmin_cond)
            {
              case(NRC_EXCLUSIVE):
              {
                if (xnew<xmin)
                {
                  return (xmin-x);
                }
                break;
              }
              case(NRC_INCLUSIVE):
              {
                if (xnew<=xmin)
                {
                  return (xmin-x)*(1-std::numeric_limits<data_type>::epsilon());
                }
                break;
              }
              case(NRC_PERIODIC):
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
              case(NRC_NOT_USED):
              {
                break;
              }
            }

            // check if max threshold is hit
            switch(xmax_cond)
            {
              case(NRC_EXCLUSIVE):
              {
                if (xnew>xmax)
                {
                  return (xmax-x);
                }
                break;
              }
              case(NRC_INCLUSIVE):
              {
                if (xnew>=xmax)
                {
                  return (xmax-x)*(1-std::numeric_limits<data_type>::epsilon());
                }
                break;
              }
              case(NRC_PERIODIC):
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
              case(NRC_NOT_USED):
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
