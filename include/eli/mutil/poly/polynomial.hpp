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

#ifndef eli_mutil_poly_polynomial_hpp
#define eli_mutil_poly_polynomial_hpp

#include <cmath>
#include <iterator>
#include <limits>

#include "eli/code_eli.hpp"

#include "eli/mutil/dm/combination.hpp"

namespace eli
{
  namespace mutil
  {
    namespace poly
    {
      template<typename data__>
      class polynomial
      {
        public:
          typedef data__ data_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> coefficient_type;
          typedef typename coefficient_type::Index index_type;

        private:
          coefficient_type a;

        public:
          polynomial() {}
          polynomial(const coefficient_type &c) : a(c) {}
          polynomial(const polynomial<data__> &p) : a(p.a) {}
          polynomial(const data_type &root) {set_roots(root);}
          polynomial(const data_type &root1, const data_type &root2) {set_roots(root1, root2);}
          polynomial(const data_type &root1, const data_type &root2, const data_type &root3) {set_roots(root1, root2, root3);}
          polynomial(const data_type &root1, const data_type &root2, const data_type &root3, const data_type &root4){set_roots(root1, root2, root3, root4);}
          template<typename itroot__>
          polynomial(itroot__ its, itroot__ ite) {set_roots(its, ite);}

          polynomial<data__> & operator=(const data_type &d)
          {
            coefficient_type a_new(1);

            a_new(0)=d;
            set_coefficients(a_new);

            return *this;
          }

          polynomial<data__> & operator=(const polynomial<data__> &p)
          {
            if (this!=&p)
            {
              a=p.a;
            }

            return *this;
          }

          index_type degree() const
          {
            return size()-1;
          }

          data_type coefficient(const index_type &i) const
          {
            return a(i);
          }
          void set_coefficients(const coefficient_type &ain)
          {
            a=ain;
          }

          void get_coefficients(coefficient_type &aout) const
          {
            aout=a;
          }

          void compress()
          {
            coefficient_type a_new;
            index_type i,j;
            for (i=this->degree(); i>0; --i)
            {
              if (a(i)!=static_cast<data_type>(0))
                break;
            }

            a_new.resize(i+1);
            for (j=0; j<=i; ++j)
              a_new(j)=a(j);

            set_coefficients(a_new);
          }

          void adjust_zero(const data_type &smval)
          {
            for (index_type i=0; i<this->size(); ++i)
            {
              if (std::abs(a(i))<smval)
                a(i)=static_cast<data_type>(0);
            }
            compress();
          }

          void set_roots(const data_type &root)
          {
            coefficient_type a_new(2);

            a_new(1)=1;
            a_new(0)=-root;
            set_coefficients(a_new);
          }

          void set_roots(const data_type &root1, const data_type &root2)
          {
            coefficient_type a_new(3);

            a_new(2)=1;
            a_new(1)=-(root1+root2);
            a_new(0)=root1*root2;
            set_coefficients(a_new);
          }

          void set_roots(const data_type &root1, const data_type &root2, const data_type &root3)
          {
            coefficient_type a_new(4);

            a_new(3)=1;
            a_new(2)=-(root1+root2+root3);
            a_new(1)=root1*root2+root1*root3+root2*root3;
            a_new(0)=-root1*root2*root3;
            set_coefficients(a_new);
          }

          void set_roots(const data_type &root1, const data_type &root2, const data_type &root3, const data_type &root4)
          {
            coefficient_type a_new(5);

            a_new(4)=1;
            a_new(3)=-(root1+root2+root3+root4);
            a_new(2)=root1*root2+root1*root3+root1*root4+root2*root3+root2*root4+root3*root4;
            a_new(1)=-(root1*root2*root3+root1*root2*root4+root1*root3*root4+root2*root3*root4);
            a_new(0)=root1*root2*root3*root4;
            set_coefficients(a_new);
          }

          /** This algorithm uses the fact that each coefficient is the sum of the multiplication
           *  of all combinations of a number of roots. See hand coded root cases to see
           *  the pattern.
           */
          template<typename itroot__>
          void set_roots(itroot__ its, itroot__ ite)
          {
            int deg=static_cast<int>(std::distance(its, ite)), n=deg+1;
            coefficient_type a_new(n);
            std::vector<data_type> root(deg);
            std::vector<int> index(deg);
            int i, j;

            // fill vector of roots
            i=0;
            for (itroot__ it=its; it!=ite; ++it, ++i)
            {
              index[i]=i;
              root[i]=(*it);
            }

            a_new[deg]=1;
            for (j=1; j<n; ++j)
            {
              a_new[deg-j]=static_cast<data_type>(0);
              do
              {
                data_type term(root[index[0]]);

                for (i=1; i<j; ++i)
                  term*=root[index[i]];
                a_new[deg-j]+=term;
              } while(eli::mutil::dm::next_combination(index.begin(), index.begin()+j, index.end()));
              if (j%2==1)
                a_new[deg-j]*=-1;
            }
            set_coefficients(a_new);
          }

          data_type f(const data_type &t) const
          {
            data_type rtn(0);

            // check to make sure have valid curve
            if (this->size()<=0)
            {
              assert(false);
              return rtn;
            }

            index_type i, n(this->degree());
            for (i=n; i>0; --i)
            {
              rtn=t*(a(i)+rtn);
            }
            i=0;
            rtn+=a(i);

            return rtn;
          }

          data_type fp(const data_type &t) const
          {
            data_type rtn(0);

            // check to make sure have valid curve
            if (this->size()<=0)
            {
              assert(false);
              return rtn;
            }

            index_type i, n;
            n=static_cast<index_type>(this->degree());
            for (i=n; i>1; --i)
            {
              rtn=t*(static_cast<data_type>(i)*a(i)+rtn);
            }
            if (n>0)
            {
              i=1;
              rtn+=static_cast<data_type>(i)*a(i);
            }

            return rtn;
          }

          data_type fpp(const data_type &t) const
          {
            data_type rtn(0);

            // check to make sure have valid curve
            if (this->size()<=0)
            {
              assert(false);
              return rtn;
            }

            index_type i, n;
            n=static_cast<index_type>(this->degree());
            for (i=n; i>2; --i)
            {
              rtn=t*(static_cast<data_type>(i)*static_cast<data_type>(i-1)*a(i)+rtn);
            }
            if (n>1)
            {
              i=2;
              rtn+=static_cast<data_type>(i)*static_cast<data_type>(i-1)*a(i);
            }

            return rtn;
          }

          data_type fppp(const data_type &t) const
          {
            data_type rtn(0);

            // check to make sure have valid curve
            if (this->size()<=0)
            {
              assert(false);
              return rtn;
            }

            index_type i, n;
            n=static_cast<index_type>(this->degree());
            for (i=n; i>3; --i)
            {
              rtn=t*(static_cast<data_type>(i)*static_cast<data_type>(i-1)*static_cast<data_type>(i-2)*a(i)+rtn);
            }
            if (n>2)
            {
              i=3;
              rtn+=static_cast<data_type>(i)*static_cast<data_type>(i-1)*static_cast<data_type>(i-2)*a(i);
            }

            return rtn;
          }

          polynomial<data_type> * f() const
          {
            if (this->size()<=0)
              return 0;

            return new polynomial<data_type>(a);
          }

          polynomial<data_type> * fp() const
          {
            if (this->size()<=0)
              return 0;

            index_type new_deg(this->degree()-1), deg(this->degree());
            coefficient_type a_new(new_deg+1);
            for (index_type i=1; i<=deg; ++i)
              a_new(i-1)=static_cast<data_type>(i)*a(i);

            return new polynomial<data_type>(a_new);
          }

          polynomial<data_type> * fpp() const
          {
            if (this->size()<=0)
              return 0;

            index_type new_deg(this->degree()-2), deg(this->degree());
            coefficient_type a_new(new_deg+1);
            for (index_type i=2; i<=deg; ++i)
              a_new(i-2)=static_cast<data_type>(i)*static_cast<data_type>(i-1)*a(i);

            return new polynomial<data_type>(a_new);
          }

          polynomial<data_type> * fppp() const
          {
            if (this->size()<=0)
              return 0;

            index_type new_deg(this->degree()-3), deg(this->degree());
            coefficient_type a_new(new_deg+1);
            for (index_type i=3; i<=deg; ++i)
              a_new(i-3)=static_cast<data_type>(i)*static_cast<data_type>(i-1)*static_cast<data_type>(i-2)*a(i);

            return new polynomial<data_type>(a_new);
          }

          void add(const polynomial<data_type> &p1, const polynomial<data_type> &p2)
          {
            // resize coefficients
            index_type deg1(p1.degree()), deg2(p2.degree()), deg, n, i;

            deg=std::max(deg1, deg2);
            n=deg+1;
            coefficient_type a_new(n);
            a_new.setZero();
            for (i=0; i<n; ++i)
            {
              if (i<=deg1)
                a_new(i)=p1.a(i);

              if (i<=deg2)
                a_new(i)+=p2.a(i);
            }

            // set new coefficients
            set_coefficients(a_new);
          }

          void subtract(const polynomial<data_type> &p1, const polynomial<data_type> &p2)
          {
            // resize coefficients
            index_type deg1(p1.degree()), deg2(p2.degree()), deg, n, i;

            deg=std::max(deg1, deg2);
            n=deg+1;
            coefficient_type a_new(n);
            a_new.setZero();
            for (i=0; i<n; ++i)
            {
              if (i<=deg1)
                a_new(i)=p1.a(i);

              if (i<=deg2)
                a_new(i)-=p2.a(i);
            }

            // set new coefficients
            set_coefficients(a_new);
          }

          void multiply(const polynomial<data_type> &p1, const polynomial<data_type> &p2)
          {
            // resize coefficients
            index_type deg1(p1.degree()), deg2(p2.degree()), deg, n, j1, j2;

            deg=deg1+deg2;
            n=deg+1;
            coefficient_type a_new(n);
            a_new.setZero();
            for (j1=0; j1<=deg1; ++j1)
            {
              for (j2=0; j2<=deg2; ++j2)
              {
                a_new(j1+j2)+=p1.a(j1)*p2.a(j2);
              }
            }

            // set new coefficients
            set_coefficients(a_new);
          }

          void multiply(const data_type &d)
          {
            for (index_type i=0; i<this->size(); ++i)
              a[i]*=d;
          }

          void divide(polynomial<data_type> &prem, const polynomial<data_type> &p1, const polynomial<data_type> &p2)
          {
            // resize coefficients
            index_type deg1(p1.degree()), deg2(p2.degree()), deg, n, j, j1, j2, nz;

            // catch case when deg2>deg1
            if (deg2>deg1)
            {
              prem=p1;
              coefficient_type a_new(1);
              a_new.setZero();
              set_coefficients(a_new);
              return;
            }

            // find the number of zero terms in the divisor to set the sizes
            for (j=0; j<=deg2; ++j)
            {
              if (p2.a(deg2-j)!=static_cast<data_type>(0))
                break;
            }
            nz=j;

            // if p2=0 then set results to nan and be done
            if(nz>deg2)
            {
              coefficient_type a_new(1);
              a_new(1)=std::numeric_limits<data__>::quiet_NaN();
              prem.set_coefficients(a_new);
              set_coefficients(a_new);
              return;
            }

            // set the size of quotient coefficients
            deg=deg1+nz-deg2;
            n=deg+1;
            coefficient_type a_new(n);
            a_new.setZero();

            // store current polynomial into remainder polynomial
            polynomial<data_type> p1_copy(p1);

            // perform the division
            deg2-=nz;

            for (j=0; j<=deg; ++j)
            {
              j1=deg-j;
              a_new(j1)=p1_copy.a(deg2+j1)/p2.a(deg2);
              for (j2=0; j2<=deg2; ++j2)
              {
                p1_copy.a(j2+j1)-=a_new(j1)*p2.a(j2);
              }
            }

            index_type nrem(deg2+nz);
            coefficient_type a_rem(nrem);
            for (j=0; j<nrem; ++j)
              a_rem(j)=p1_copy.a(j);

            // set new coefficients
            prem.set_coefficients(a_rem);
            set_coefficients(a_new);
          }

          void divide(const data_type &d)
          {
            for (index_type i=0; i<this->size(); ++i)
              a[i]/=d;
          }

          void negative()
          {
            for (index_type i=0; i<this->size(); ++i)
              a[i]*=static_cast<data_type>(-1);
          }

        private:
          index_type size() const
          {
            return a.rows();
          }

      };
    }
  }
}

#endif
