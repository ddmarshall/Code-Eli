/*********************************************************************************
* CopyRIGHT (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All RIGHTs reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_mutil_fd_d2o1_hpp
#define eli_mutil_fd_d2o1_hpp

#include <vector>

namespace eli
{
  namespace mutil
  {
    namespace fd
    {
      template<typename data__>
      class d2o1
      {
        public:
          enum stencil
          {
            LEFT=0,
            CENTER=1,
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
          d2o1() : nnodes(3), n_order(1), st(d2o1<data__>::CENTER)
          {
          }

          d2o1(const stencil &s) : nnodes(3), n_order(1), st(s)
          {
          }

          d2o1(const d2o1<data__> &d) : nnodes(3), n_order(1), st(d.st)
          {
          }

          ~d2o1()
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

          int order(bool uniform) const
          {
            if (uniform && (st==CENTER))
              return n_order+1;
            else
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
                i0=2;
                (*iti)=-2;++iti;
                (*iti)=-1;++iti;
                (*iti)=0;
                break;
              }
              case(CENTER):
              {
                i0=1;
                (*iti)=-1;++iti;
                (*iti)=0;++iti;
                (*iti)=1;
                break;
              }
              case(RIGHT):
              {
                i0=0;
                (*iti)=0;++iti;
                (*iti)=1;++iti;
                (*iti)=2;
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
              case(LEFT):
              {
                (*itc)=static_cast<data__>( 1)/dx/dx;++itc;
                (*itc)=static_cast<data__>(-2)/dx/dx;++itc;
                (*itc)=static_cast<data__>( 1)/dx/dx;

                break;
              }
              case(CENTER):
              {
                (*itc)=static_cast<data__>( 1)/dx/dx;++itc;
                (*itc)=static_cast<data__>(-2)/dx/dx;++itc;
                (*itc)=static_cast<data__>( 1)/dx/dx;

                break;
              }
              case(RIGHT):
              {
                (*itc)=static_cast<data__>( 1)/dx/dx;++itc;
                (*itc)=static_cast<data__>(-2)/dx/dx;++itc;
                (*itc)=static_cast<data__>( 1)/dx/dx;

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

            // extract the x-locations
            x[0]=(*itx); ++itx;
            x[1]=(*itx); ++itx;
            x[2]=(*itx);

            switch(st)
            {
              case(LEFT):
              {
                data__ alphaim1, betaim2;
                alphaim1=x[2]-x[1];
                betaim2=x[2]-x[0];

                (*itc)= 2/(betaim2*(betaim2-alphaim1));++itc;
                (*itc)=-2/(alphaim1*(betaim2-alphaim1));++itc;
                (*itc)= 2/(alphaim1*betaim2);

                break;
              }
              case(CENTER):
              {
                data__ alphai, alphaim1;
                alphaim1=x[1]-x[0];
                alphai=x[2]-x[1];

                (*itc)= 2/(alphaim1*(alphai+alphaim1));++itc;
                (*itc)=-2/(alphai*alphaim1);++itc;
                (*itc)= 2/(alphai*(alphai+alphaim1));

                break;
              }
              case(RIGHT):
              {
                data__ alphai, betai;
                alphai=x[1]-x[0];
                betai=x[2]-x[0];

                (*itc)= 2/(alphai*betai);++itc;
                (*itc)=-2/(alphai*(betai-alphai));++itc;
                (*itc)= 2/(betai*(betai-alphai));

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

          int truncation_error(data__ &te, const data__ &phi3, const data__ &dx) const
          {
            switch(st)
            {
              case(LEFT):
              {
                te=phi3*dx;
                break;
              }
              case(CENTER):
              {
                te=-dx*dx/12;
                break;
              }
              case(RIGHT):
              {
                te=-phi3*dx;
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

          template<typename itx__> int truncation_error(data__ &te, const data__ &phi3, itx__ itx) const
          {
            std::vector<data__> x(number_nodes());

            // extract the x-locations
            x[0]=(*itx); ++itx;
            x[1]=(*itx); ++itx;
            x[2]=(*itx);

            switch(st)
            {
              case(LEFT):
              {
                data__ alphaim1, betaim2;
                alphaim1=x[2]-x[1];
                betaim2=x[2]-x[0];

                te=phi3*(alphaim1+betaim2)/3;
                break;
              }
              case(CENTER):
              {
                data__ alphai, alphaim1;
                alphaim1=x[1]-x[0];
                alphai=x[2]-x[1];

                if (alphai==alphaim1)
                  te=-alphai*alphaim1/12;
                else
                  te=-phi3*(alphai-alphaim1)/3;
                break;
              }
              case(RIGHT):
              {
                data__ alphai, betai;
                alphai=x[1]-x[0];
                betai=x[2]-x[0];

                te=-phi3*(alphai+betai)/3;
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
      };
    }
  }
}

#endif
