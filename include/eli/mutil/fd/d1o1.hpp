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

#ifndef eli_mutil_fd_d1o1_hpp
#define eli_mutil_fd_d1o1_hpp

#include <vector>

namespace eli
{
  namespace mutil
  {
    namespace fd
    {
      template<typename data__>
      class d1o1
      {
        public:
          enum stencil
          {
            LEFT=0,
            RIGHT=2
          };

        private:
          const size_t nnodes;
          const int n_order;
          stencil st;

        protected:
          template<typename itc__, typename itphi__> data__ calculate_dot(itc__ a, itphi__ itphi) const
          {
            data__ d(static_cast<data__>(0));

            for (size_t i=0; i<number_nodes(); ++i, ++itphi)
              d+=a[i]*(*itphi);

            return d;
          }

        public:
          d1o1() : nnodes(2), n_order(1), st(LEFT)
          {
          }

          d1o1(const stencil &s) : nnodes(2), n_order(1), st(s)
          {
          }

          d1o1(const d1o1<data__> &d) : nnodes(2), n_order(1), st(d.st)
          {
          }

          ~d1o1()
          {
          }

          void set_stencil(const stencil &s)
          {
            st=s;
          }

          const stencil & get_stencil() const
          {
            return st;
          }

          int order(bool /*uniform*/) const
          {
            return n_order;
          }

          size_t number_nodes() const
          {
            return nnodes;
          }

          template<typename iti__> std::ptrdiff_t index(iti__ iti) const
          {
            size_t i0;

            switch(st)
            {
              case(LEFT):
              {
                i0=1;
                (*iti)=-1;++iti;
                (*iti)=0;
                break;
              }
              case(RIGHT):
              {
                i0=0;
                (*iti)=0;++iti;
                (*iti)=1;
                break;
              }
              default:
              {
                i0=3;
                assert(false);
                return i0;
              }
            }

            return i0;
          }

          template<typename itphi__> int evaluate(data__ &d, itphi__ itphi, const data__ &dx) const
          {
            std::vector<data__> a(number_nodes());
            int rtn;

            rtn=coefficients(a.begin(), dx);
            if (rtn==0)
              d=calculate_dot(a.begin(), itphi);

            return rtn;
          }

          template<typename itphi__, typename itx__> int evaluate(data__ &d, itphi__ itphi, itx__ itx) const
          {
            std::vector<data__> a(number_nodes());
            int rtn;

            rtn=coefficients(a.begin(), itx);
            if (rtn==0)
              d=calculate_dot(a.begin(), itphi);

            return rtn;
          }

          template<typename itc__> int coefficients(itc__ itc, const data__ &dx) const
          {
            switch(st)
            {
              case(RIGHT):
              case(LEFT):
              {
                (*itc)=static_cast<data__>(-1.0)/dx;++itc;
                (*itc)=static_cast<data__>( 1.0)/dx;

                break;
              }
              default:
              {
                assert(false);
                return -1;
              }
             }

            return 0;
          }

          template<typename itc__, typename itx__> int coefficients(itc__ itc, itx__ itx) const
          {
            std::vector<data__> x(number_nodes());
            data__ dx;

            // extract the x-locations
            x[0]=(*itx); ++itx;
            x[1]=(*itx);
            dx=x[1]-x[0];

            return coefficients(itc, dx);
          }

          int truncation_error(data__ &te, const data__ &phi2, const data__ &dx) const
          {
            switch(st)
            {
              case(LEFT):
              {
                te=phi2*dx/2;
                break;
              }
              case(RIGHT):
              {
                te=-phi2*dx/2;
                break;
              }
              default:
              {
                assert(false);
                return -1;
              }
             }

            return 0;
          }

          template<typename itx__> int truncation_error(data__ &te, const data__ &phi2, itx__ itx) const
          {
            std::vector<data__> x(number_nodes());
            data__ dx;

            // extract the x-locations
            x[0]=(*itx); ++itx;
            x[1]=(*itx);
            dx=x[1]-x[0];

            return truncation_error(te, phi2, dx);
          }
      };
    }
  }
}

#endif
