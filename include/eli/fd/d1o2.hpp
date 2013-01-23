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

#ifndef eli_fd_d1o2_hpp
#define eli_fd_d1o2_hpp

namespace eli
{
  namespace fd
  {
    template<typename data__>
    class d1o2
    {
      public:
        enum stencil
        {
          left=0,
          center=1,
          right=2
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
        d1o2() : nnodes(3), n_order(2), st(d1o2<data__>::center)
        {
        }

        d1o2(const stencil &s) : nnodes(3), n_order(2), st(s)
        {
        }

        d1o2(const d1o2<data__> &d) : nnodes(3), n_order(2), st(d.st)
        {
        }

        ~d1o2()
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
            case(left):
            {
              i0=2;
              (*iti)=-2;++iti;
              (*iti)=-1;++iti;
              (*iti)=0;
              break;
            }
            case(center):
            {
              i0=1;
              (*iti)=-1;++iti;
              (*iti)=0;++iti;
              (*iti)=1;
              break;
            }
            case(right):
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
            case(left):
            {
              (*itc)=static_cast<data__>( 0.5)/dx;++itc;
              (*itc)=static_cast<data__>(-2.0)/dx;++itc;
              (*itc)=static_cast<data__>( 1.5)/dx;

              break;
            }
            case(center):
            {
              (*itc)=static_cast<data__>(-0.5)/dx;++itc;
              (*itc)=static_cast<data__>( 0.0)/dx;++itc;
              (*itc)=static_cast<data__>( 0.5)/dx;

              break;
            }
            case(right):
            {
              (*itc)=static_cast<data__>(-1.5)/dx;++itc;
              (*itc)=static_cast<data__>( 2.0)/dx;++itc;
              (*itc)=static_cast<data__>(-0.5)/dx;

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
            case(left):
            {
              data__ alphaim1, betaim2, alphaim2;
              alphaim1=x[2]-x[1];
              betaim2=x[2]-x[0];
              alphaim2=x[1]-x[0];

              (*itc)=alphaim1/(betaim2*alphaim2);++itc;
              (*itc)=-betaim2/(alphaim1*alphaim2);++itc;
              (*itc)=(alphaim1+betaim2)/(alphaim1*betaim2);

              break;
            }
            case(center):
            {
              data__ deli, delim1, gamip1;
              gamip1=x[2]-x[0];
              delim1=x[1]-x[0];
              deli=x[2]-x[1];

              (*itc)=-deli/(delim1*gamip1);++itc;
              (*itc)=(deli-delim1)/(deli*delim1);++itc;
              (*itc)=delim1/(deli*gamip1);

              break;
            }
            case(right):
            {
              data__ alphai, betai, alphaip1;
              alphai=x[1]-x[0];
              betai=x[2]-x[0];
              alphaip1=x[2]-x[1];

              (*itc)=-(alphai+betai)/(alphai*betai);++itc;
              (*itc)=betai/(alphai*alphaip1);++itc;
              (*itc)=-alphai/(betai*alphaip1);

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
            case(left):
            {
              te=phi3*dx*dx/3;
              break;
            }
            case(center):
            {
              te=-phi3*dx*dx/6;
              break;
            }
            case(right):
            {
              te=phi3*dx*dx/3;
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
            case(left):
            {
              data__ delim1, betaim2;

              delim1=x[2]-x[1];
              betaim2=x[2]-x[0];

              te=phi3*delim1*betaim2/6;
              break;
            }
            case(center):
            {
              data__ deli, delim1;

              deli=x[2]-x[1];
              delim1=x[1]-x[0];

              te=-phi3*deli*delim1/6;
              break;
            }
            case(right):
            {
              data__ alphai, betai;

              alphai=x[1]-x[0];
              betai=x[2]-x[0];

              te=phi3*alphai*betai/6;
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
