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

#ifndef eli_fd_d2o2_hpp
#define eli_fd_d2o2_hpp

namespace eli
{
  namespace fd
  {
    template<typename __data>
    class d2o2
    {
      public:
        enum stencil
        {
          left=0,
          left_biased=1,
          right_biased=2,
          right=3
        };

      private:
        const size_t nnodes;
        const int n_order;
        stencil st;

      protected:
        template<typename __itc, typename __itphi> __data calculate_dot(__itc a, __itphi itphi) const
        {
          __data d(static_cast<__data>(0));

          for (size_t i=0; i<number_nodes(); ++i, ++itphi)
            d+=a[i]*(*itphi);

          return d;
        }

      public:
        d2o2() : nnodes(4), n_order(2), st(d2o2<__data>::left_biased)
        {
        }

        d2o2(const stencil &s) : nnodes(4), n_order(2), st(s)
        {
        }

        d2o2(const d2o2<__data> &d) : nnodes(4), n_order(2), st(d.st)
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
            case(left):
            {
              i0=3;
              (*iti)=-3;++iti;
              (*iti)=-2;++iti;
              (*iti)=-1;++iti;
              (*iti)= 0;
              break;
            }
            case(left_biased):
            {
              i0=2;
              (*iti)=-2;++iti;
              (*iti)=-1;++iti;
              (*iti)= 0;++iti;
              (*iti)= 1;
              break;
            }
            case(right_biased):
            {
              i0=1;
              (*iti)=-1;++iti;
              (*iti)= 0;++iti;
              (*iti)= 1;++iti;
              (*iti)= 2;
              break;
            }
            case(right):
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

        template<typename __itphi> int evaluate(__data &d, __itphi itphi, const __data &dx) const
        {
          std::vector<__data> a(number_nodes());
          int rtn;

          rtn=coefficients(a.begin(), dx);
          if (rtn==0)
            d=calculate_dot(a.begin(), itphi);

          return rtn;
        }

        template<typename __itphi, typename __itx> int evaluate(__data &d, __itphi itphi, __itx itx) const
        {
          std::vector<__data> a(number_nodes());
          int rtn;

          rtn=coefficients(a.begin(), itx);
          if (rtn==0)
            d=calculate_dot(a.begin(), itphi);

          return rtn;
        }

        template<typename __itc> int coefficients(__itc itc, const __data &dx) const
        {
          switch(st)
          {
            case(left):
            {
              (*itc)=static_cast<__data>(-1)/dx/dx;++itc;
              (*itc)=static_cast<__data>( 4)/dx/dx;++itc;
              (*itc)=static_cast<__data>(-5)/dx/dx;++itc;
              (*itc)=static_cast<__data>( 2)/dx/dx;

              break;
            }
            case(left_biased):
            {
              (*itc)=static_cast<__data>( 0)/dx/dx;++itc;
              (*itc)=static_cast<__data>( 1)/dx/dx;++itc;
              (*itc)=static_cast<__data>(-2)/dx/dx;++itc;
              (*itc)=static_cast<__data>( 1)/dx/dx;

              break;
            }
            case(right_biased):
            {
              (*itc)=static_cast<__data>( 1)/dx/dx;++itc;
              (*itc)=static_cast<__data>(-2)/dx/dx;++itc;
              (*itc)=static_cast<__data>( 1)/dx/dx;++itc;
              (*itc)=static_cast<__data>( 0)/dx/dx;

              break;
            }
            case(right):
            {
              (*itc)=static_cast<__data>( 2)/dx/dx;++itc;
              (*itc)=static_cast<__data>(-5)/dx/dx;++itc;
              (*itc)=static_cast<__data>( 4)/dx/dx;++itc;
              (*itc)=static_cast<__data>(-1)/dx/dx;

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
          std::vector<__data> x(number_nodes());

          // extract the x-locations
          for (size_t i=0; i<number_nodes(); ++i, ++itx)
            x[i]=(*itx);

          switch(st)
          {
            case(left):
            {
              __data alphaim1, betaim2, gammaim3;
              alphaim1=x[3]-x[2];
              betaim2=x[3]-x[1];
              gammaim3=x[3]-x[0];

              (*itc)=-2*(alphaim1+betaim2)/(gammaim3*(gammaim3-betaim2)*(gammaim3-alphaim1));++itc;
              (*itc)=+2*(alphaim1+gammaim3)/(betaim2*(gammaim3-betaim2)*(betaim2-alphaim1));++itc;
              (*itc)=-2*(gammaim3+betaim2)/(alphaim1*(betaim2-alphaim1)*(gammaim3-alphaim1));++itc;
              (*itc)=+2*(alphaim1+betaim2+gammaim3)/(alphaim1*betaim2*gammaim3);

              break;
            }
            case(left_biased):
            {
              __data alphai, alphaim1, betaim2;
              alphai=x[3]-x[2];
              alphaim1=x[2]-x[1];
              betaim2=x[2]-x[0];

              (*itc)=+2*(alphai-alphaim1)/(betaim2*(alphai+betaim2)*(betaim2-alphaim1));++itc;
              (*itc)=+2*(betaim2-alphai)/(alphaim1*(alphai+alphaim1)*(betaim2-alphaim1));++itc;
              (*itc)=-2*(alphaim1-alphai+betaim2)/(alphaim1*alphai*betaim2);++itc;
              (*itc)=+2*(alphaim1+betaim2)/(alphai*(alphai+alphaim1)*(alphai+betaim2));

              break;
            }
            case(right_biased):
            {
              __data alphai, alphaim1, betai;
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
            case(right):
            {
              __data alphai, betai, gammai;
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

        int truncation_error(__data &te, const __data &phi4, const __data &dx) const
        {
          switch(st)
          {
            case(left):
            {
              te=phi4*11*dx*dx/12;
              break;
            }
            case(left_biased):
            {
              te=-phi4*dx*dx/12;
              break;
            }
            case(right_biased):
            {
              te=-phi4*dx*dx/12;
              break;
            }
            case(right):
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

        template<typename __itx> int truncation_error(__data &te, const __data &phi4, __itx itx) const
        {
          std::vector<__data> x(number_nodes());

          // extract the x-locations
          for (size_t i=0; i<number_nodes(); ++i, ++itx)
            x[i]=(*itx);

          switch(st)
          {
            case(left):
            {
              __data alphaim1, betaim2, gammaim3;
              alphaim1=x[3]-x[2];
              betaim2=x[3]-x[1];
              gammaim3=x[3]-x[0];

              te=phi4*(alphaim1*gammaim3+betaim2*(alphaim1+gammaim3))/static_cast<__data>(12);
              break;
            }
            // FIX: These need to be fixed for d2o2
            case(left_biased):
            {
              __data alphai, alphaim1, betaim2;
              alphai=x[3]-x[2];
              alphaim1=x[2]-x[1];
              betaim2=x[2]-x[0];

              te=-phi4*(alphai*(alphaim1+betaim2)-alphaim1*betaim2)/static_cast<__data>(12);
              break;
            }
            case(right_biased):
            {
              __data alphai, alphaim1, betai;
              alphaim1=x[1]-x[0];
              alphai=x[2]-x[1];
              betai=x[3]-x[1];

              te=-phi4*(alphaim1*(alphai+betai)-alphai*betai)/static_cast<__data>(12);
              break;
            }
            case(right):
            {
              __data alphai, betai, gammai;
              alphai=x[1]-x[0];
              betai=x[2]-x[0];
              gammai=x[3]-x[0];

              te=phi4*(alphai*gammai+betai*(alphai+gammai))/static_cast<__data>(12);
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

#endif
