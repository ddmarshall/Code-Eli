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

#ifndef eli_mutil_fd_d2o2_hpp
#define eli_mutil_fd_d2o2_hpp

#include <vector>
#include <assert.h>

namespace eli
{
  namespace mutil
  {
    namespace fd
    {
      template<typename data__>
      class d2o2
      {
        public:
          enum stencil
          {
            LEFT=0,
            LEFT_BIASED=1,
            RIGHT_BIASED=2,
            RIGHT=3
          };

        private:
          const size_t nnodes;
          const int n_order;
          stencil st;

        protected:
          template<typename __itc, typename __itphi> data__ calculate_dot(__itc a, __itphi itphi) const
          {
            data__ d(static_cast<data__>(0));

            for (size_t i=0; i<number_nodes(); ++i, ++itphi)
            {
              if (a[i]!=0)
                d+=a[i]*(*itphi);
            }

            return d;
          }

        public:
          d2o2() : nnodes(4), n_order(2), st(d2o2<data__>::RIGHT_BIASED)
          {
          }

          d2o2(const stencil &s) : nnodes(4), n_order(2), st(s)
          {
          }

          d2o2(const d2o2<data__> &d) : nnodes(4), n_order(2), st(d.st)
          {
          }

          ~d2o2()
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

          template<typename __iti> std::ptrdiff_t index(__iti iti) const
          {
            size_t i0;

            switch(st)
            {
              case(LEFT):
              {
                i0=3;
                (*iti)=-3;++iti;
                (*iti)=-2;++iti;
                (*iti)=-1;++iti;
                (*iti)= 0;
                break;
              }
              case(LEFT_BIASED):
              {
                i0=2;
                (*iti)=-2;++iti;
                (*iti)=-1;++iti;
                (*iti)= 0;++iti;
                (*iti)= 1;
                break;
              }
              case(RIGHT_BIASED):
              {
                i0=1;
                (*iti)=-1;++iti;
                (*iti)= 0;++iti;
                (*iti)= 1;++iti;
                (*iti)= 2;
                break;
              }
              case(RIGHT):
              {
                i0=0;
                (*iti)=0;++iti;
                (*iti)=1;++iti;
                (*iti)=2;++iti;
                (*iti)=3;
                break;
              }
              default:
              {
                i0=4;
                assert(false);
                return i0;
              }
            }

            return i0;
          }

          template<typename __itphi> int evaluate(data__ &d, __itphi itphi, const data__ &dx) const
          {
            std::vector<data__> a(number_nodes());
            int rtn;

            rtn=coefficients(a.begin(), dx);
            if (rtn==0)
              d=calculate_dot(a.begin(), itphi);

            return rtn;
          }

          template<typename __itphi, typename __itx> int evaluate(data__ &d, __itphi itphi, __itx itx) const
          {
            std::vector<data__> a(number_nodes());
            int rtn;

            rtn=coefficients(a.begin(), itx);
            if (rtn==0)
              d=calculate_dot(a.begin(), itphi);

            return rtn;
          }

          template<typename __itc> int coefficients(__itc itc, const data__ &dx) const
          {
            switch(st)
            {
              case(LEFT):
              {
                (*itc)=static_cast<data__>(-1)/dx/dx;++itc;
                (*itc)=static_cast<data__>( 4)/dx/dx;++itc;
                (*itc)=static_cast<data__>(-5)/dx/dx;++itc;
                (*itc)=static_cast<data__>( 2)/dx/dx;

                break;
              }
              case(LEFT_BIASED):
              {
                (*itc)=static_cast<data__>( 0)/dx/dx;++itc;
                (*itc)=static_cast<data__>( 1)/dx/dx;++itc;
                (*itc)=static_cast<data__>(-2)/dx/dx;++itc;
                (*itc)=static_cast<data__>( 1)/dx/dx;

                break;
              }
              case(RIGHT_BIASED):
              {
                (*itc)=static_cast<data__>( 1)/dx/dx;++itc;
                (*itc)=static_cast<data__>(-2)/dx/dx;++itc;
                (*itc)=static_cast<data__>( 1)/dx/dx;++itc;
                (*itc)=static_cast<data__>( 0)/dx/dx;

                break;
              }
              case(RIGHT):
              {
                (*itc)=static_cast<data__>( 2)/dx/dx;++itc;
                (*itc)=static_cast<data__>(-5)/dx/dx;++itc;
                (*itc)=static_cast<data__>( 4)/dx/dx;++itc;
                (*itc)=static_cast<data__>(-1)/dx/dx;

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

          template<typename __itc, typename __itx> int coefficients(__itc itc, __itx itx) const
          {
            std::vector<data__> x(number_nodes());

            // extract the x-locations
            for (size_t i=0; i<number_nodes(); ++i, ++itx)
              x[i]=(*itx);

            switch(st)
            {
              case(LEFT):
              {
                data__ alphaim1, betaim2, gammaim3;
                alphaim1=x[3]-x[2];
                betaim2=x[3]-x[1];
                gammaim3=x[3]-x[0];

                (*itc)=-2*(alphaim1+betaim2)/(gammaim3*(gammaim3-betaim2)*(gammaim3-alphaim1));++itc;
                (*itc)=+2*(alphaim1+gammaim3)/(betaim2*(gammaim3-betaim2)*(betaim2-alphaim1));++itc;
                (*itc)=-2*(gammaim3+betaim2)/(alphaim1*(betaim2-alphaim1)*(gammaim3-alphaim1));++itc;
                (*itc)=+2*(alphaim1+betaim2+gammaim3)/(alphaim1*betaim2*gammaim3);

                break;
              }
              case(LEFT_BIASED):
              {
                data__ alphai, alphaim1, betaim2;
                alphai=x[3]-x[2];
                alphaim1=x[2]-x[1];
                betaim2=x[2]-x[0];

                (*itc)=+2*(alphai-alphaim1)/(betaim2*(alphai+betaim2)*(betaim2-alphaim1));++itc;
                (*itc)=+2*(betaim2-alphai)/(alphaim1*(alphai+alphaim1)*(betaim2-alphaim1));++itc;
                (*itc)=-2*(alphaim1-alphai+betaim2)/(alphaim1*alphai*betaim2);++itc;
                (*itc)=+2*(alphaim1+betaim2)/(alphai*(alphai+alphaim1)*(alphai+betaim2));

                break;
              }
              case(RIGHT_BIASED):
              {
                data__ alphai, alphaim1, betai;
                alphaim1=x[1]-x[0];
                alphai=x[2]-x[1];
                betai=x[3]-x[1];

                (*itc)=+2*(alphai+betai)/(alphaim1*(alphaim1+betai)*(alphai+alphaim1));++itc;
                (*itc)=-2*(alphai-alphaim1+betai)/(alphai*alphaim1*betai);++itc;
                (*itc)=+2*(betai-alphaim1)/(alphai*(alphai+alphaim1)*(betai-alphai));++itc;
                (*itc)=-2*(alphai-alphaim1)/(betai*(alphaim1+betai)*(betai-alphai));

                break;
              }
              // FIX: Change for 3rd order
              case(RIGHT):
              {
                data__ alphai, betai, gammai;
                alphai=x[1]-x[0];
                betai=x[2]-x[0];
                gammai=x[3]-x[0];

                (*itc)=+2*(alphai+betai+gammai)/(alphai*betai*gammai);++itc;
                (*itc)=-2*(betai+gammai)/(alphai*(gammai-alphai)*(betai-alphai));++itc;
                (*itc)=+2*(alphai+gammai)/(betai*(gammai-betai)*(betai-alphai));++itc;
                (*itc)=-2*(alphai+betai)/(gammai*(gammai-alphai)*(gammai-betai));

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

          int truncation_error(data__ &te, const data__ &phi4, const data__ &dx) const
          {
            switch(st)
            {
              case(LEFT):
              {
                te=phi4*11*dx*dx/12;
                break;
              }
              case(LEFT_BIASED):
              {
                te=-phi4*dx*dx/12;
                break;
              }
              case(RIGHT_BIASED):
              {
                te=-phi4*dx*dx/12;
                break;
              }
              case(RIGHT):
              {
                te=phi4*11*dx*dx/12;
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

          template<typename __itx> int truncation_error(data__ &te, const data__ &phi4, __itx itx) const
          {
            std::vector<data__> x(number_nodes());

            // extract the x-locations
            for (size_t i=0; i<number_nodes(); ++i, ++itx)
              x[i]=(*itx);

            switch(st)
            {
              case(LEFT):
              {
                data__ alphaim1, betaim2, gammaim3;
                alphaim1=x[3]-x[2];
                betaim2=x[3]-x[1];
                gammaim3=x[3]-x[0];

                te=phi4*(alphaim1*gammaim3+betaim2*(alphaim1+gammaim3))/static_cast<data__>(12);
                break;
              }
              // FIX: These need to be fixed for d2o2
              case(LEFT_BIASED):
              {
                data__ alphai, alphaim1, betaim2;
                alphai=x[3]-x[2];
                alphaim1=x[2]-x[1];
                betaim2=x[2]-x[0];

                te=-phi4*(alphai*(alphaim1+betaim2)-alphaim1*betaim2)/static_cast<data__>(12);
                break;
              }
              case(RIGHT_BIASED):
              {
                data__ alphai, alphaim1, betai;
                alphaim1=x[1]-x[0];
                alphai=x[2]-x[1];
                betai=x[3]-x[1];

                te=-phi4*(alphaim1*(alphai+betai)-alphai*betai)/static_cast<data__>(12);
                break;
              }
              case(RIGHT):
              {
                data__ alphai, betai, gammai;
                alphai=x[1]-x[0];
                betai=x[2]-x[0];
                gammai=x[3]-x[0];

                te=phi4*(alphai*gammai+betai*(alphai+gammai))/static_cast<data__>(12);
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
