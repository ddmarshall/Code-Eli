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

#ifndef eli_geom_curve_utility_bezier_hpp
#define eli_geom_curve_utility_bezier_hpp

#include "Eigen/Eigen"

#include "eli/dm/binomial_coefficient.hpp"

namespace eli
{
  namespace geom
  {
    namespace utility
    {
      template<typename Derived1, typename Derived2>
      void de_casteljau(Eigen::MatrixBase<Derived1> &p, const Eigen::MatrixBase<Derived2> &cp, const typename Derived2::Scalar &t)
      {
        // do some checks on incoming matrix dimensions
        assert(p.cols()==cp.cols());

        Eigen::Matrix<typename Derived2::Scalar, Eigen::Dynamic, Eigen::Dynamic> Q(cp);
        typename Derived2::Scalar one(1);
        typename Derived2::Index i, k;

        for (k=1; k<Q.rows(); ++k)
        {
          for (i=0; i<Q.rows()-k; ++i)
          {
            Q.row(i)=(one-t)*Q.row(i)+t*Q.row(i+1);
          }
        }

        p=Q.row(0);
      }

      template<typename Derived1, typename Derived2>
      void bezier_p_control_point(Eigen::MatrixBase<Derived1> &cp_p, const Eigen::MatrixBase<Derived2> &cp)
      {
        // do some checks on incoming matrix dimensions
        assert(cp_p.rows()==cp.rows()-1);
        assert(cp_p.cols()==cp.cols());

        typename Derived2::Index i, n(cp.rows()-1);

        for (i=0; i<n; ++i)
        {
          cp_p.row(i)=n*(cp.row(i+1)-cp.row(i));
        }
      }

      template<typename Derived1, typename Derived2>
      void bezier_pp_control_point(Eigen::MatrixBase<Derived1> &cp_pp, const Eigen::MatrixBase<Derived2> &cp)
      {
        // do some checks on incoming matrix dimensions
        assert(cp_pp.rows()==cp.rows()-2);
        assert(cp_pp.cols()==cp.cols());

        typename Derived2::Index i, n(cp.rows()-1);

        for (i=0; i<n-1; ++i)
        {
          cp_pp.row(i)=n*(n-1)*(cp.row(i+2)-2*cp.row(i+1)+cp.row(i));
        }
      }

      template<typename Derived1, typename Derived2>
      void bezier_ppp_control_point(Eigen::MatrixBase<Derived1> &cp_ppp, const Eigen::MatrixBase<Derived2> &cp)
      {
        // do some checks on incoming matrix dimensions
        assert(cp_ppp.rows()==cp.rows()-3);
        assert(cp_ppp.cols()==cp.cols());

        typename Derived2::Index i, n(cp.rows()-1);

        for (i=0; i<n-2; ++i)
        {
          cp_ppp.row(i)=n*(n-1)*(n-2)*(cp.row(i+3)-3*cp.row(i+2)+3*cp.row(i+1)-cp.row(i));
        }
      }

      template<typename Derived1, typename Derived2>
      void bezier_promote_control_points(Eigen::MatrixBase<Derived1> &cp_out, const Eigen::MatrixBase<Derived2> &cp_in)
      {
        // do some dimension checks
        assert(cp_out.rows()==cp_in.rows()+1);
        assert(cp_out.cols()==cp_in.cols());

        typename Derived1::Index i, n(cp_out.rows()-1);

        cp_out.row(0)=cp_in.row(0);
        cp_out.row(n)=cp_in.row(n-1);
        for (i=1; i<n; ++i)
          cp_out.row(i)=(cp_in.row(i-1)-cp_in.row(i))*(static_cast<typename Derived1::Scalar>(i)/n)+cp_in.row(i);
      }

      // NOTE: Implements Eck's method: Matthias Eck. Least Squares Degree reduction of Bézier curves.
      //                                Computer-Aided Design. Volume 27, No. 11. (1995), 845-851
      // NOTE: This paper also provides an algorithm for creating a bezier spline order reduced
      //       curve given a bezier spline and a desired tolerance. It also provides a means of
      //       reducing more than one order at once.
      template<typename Derived1, typename Derived2>
      void bezier_demote_control_points(Eigen::MatrixBase<Derived1> &cp_out, const Eigen::MatrixBase<Derived2> &cp_in, int ncon)
      {
        typedef typename Derived1::Index index_type;
        typedef typename Derived1::Scalar data_type;

        // do some dimension checks
        assert(cp_out.rows()==cp_in.rows()-1);
        assert(cp_out.cols()==cp_in.cols());

        index_type n(cp_in.rows()-1), i;
        Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> B_I(n, cp_in.cols()), B_II(n, cp_in.cols());
        Eigen::Matrix<data_type, Eigen::Dynamic, 1> lambda(n);

        B_I.row(0)=cp_in.row(0);
        for (i=1; i<n; ++i)
          B_I.row(i)=(cp_in.row(i)*static_cast<data_type>(n)-B_I.row(i-1)*static_cast<data_type>(i))/static_cast<data_type>(n-i);
        B_II.row(n-1)=cp_in.row(n);
        for (i=n-1; i>=1; --i)
          B_II.row(i-1)=(cp_in.row(i)*static_cast<data_type>(n)-B_II.row(i)*static_cast<data_type>(n-i))/static_cast<data_type>(i);

        // calculate the blending terms
        if (ncon==0)
        {
          // no interpolation on end points
          // note: this comes from "Degree Reduction of Bezier Curves" by Dave Morgan
          data_type coef(1/static_cast<data_type>(1<<(2*n-1))), sum(1), tmp;
          lambda(0)=coef*sum;
          for (i=1; i<n; ++i)
          {
            eli::dm::n_choose_k(tmp, 2*n, 2*i);
            sum+=tmp;
            lambda(i)=coef*sum;
          }
        }
        else
        {
          // interpolating end points up-to alpha order derivative
          data_type tmp, coef, sum(0);
          index_type alpha(ncon/2);
          eli::dm::n_choose_k(tmp, 2*n, n+2*alpha);
          coef=1/tmp;
          lambda(0)=coef*sum;
          for (i=1; i<n; ++i)
          {
            if ( (i<alpha) || (i+alpha>n) )
              sum+=0;
            else
            {
              data_type tmp1, tmp2;
              eli::dm::n_choose_k(tmp1, n, i-alpha);
              eli::dm::n_choose_k(tmp2, n, i+alpha);
              sum+=tmp1*tmp2;
            }
            lambda(i)=coef*sum;
          }
        }

        // calculate new control points
        for (i=0; i<n; ++i)
          cp_out.row(i)=B_I.row(i)*(1-lambda(i))+B_II.row(i)*lambda(i);
      }

      // get the coefficients for the split curves by constructing a de Casteljau triangle scheme.
      // The lower side is the t<t0 section, and the higher side is the t>t0 section.
      template<typename Derived1, typename Derived2>
      void bezier_split_control_points(Eigen::MatrixBase<Derived1> &cp_lo, Eigen::MatrixBase<Derived1> &cp_hi,
                                       const Eigen::MatrixBase<Derived2> &cp_in, const typename Derived2::Scalar &t)
      {
        typename Derived2::Index i, j, n(cp_in.rows()-1);
        Eigen::Matrix<typename Derived2::Scalar, Eigen::Dynamic, Eigen::Dynamic> tri(cp_in);

        // do some dimensions check
        assert(cp_lo.rows()==cp_hi.rows());
        assert(cp_lo.cols()==cp_hi.cols());
        assert(cp_lo.rows()==cp_in.rows());
        assert(cp_lo.cols()==cp_in.cols());

        // set the control points using de Casteljau's algorithm
        for (i=0; i<=n; ++i)
        {
          cp_lo.row(i)=tri.row(0);
          cp_hi.row(n-i)=tri.row(n-i);
          for (j=0; j<n-i; ++j)
          {
            tri.row(j)=tri.row(j+1)*t+tri.row(j)*(1-t);
          }
        }
      }

      template<typename Derived>
      void bezier_N(Eigen::MatrixBase<Derived> &N, const typename Derived::Index &n)
      {
        typename Derived::Index i,j;
        typename Derived::Scalar bc1, bc2;

        // size N
        N.derived().resize(n+1,n+1);
        N.setZero();
        for (i=0; i<=n; ++i)
        {
          for (j=0; j<=n; ++j)
          {
            if (i+j<=n)
            {
              eli::dm::n_choose_k(bc1, n, j);
              eli::dm::n_choose_k(bc2, n-j, n-i-j);
              N(i,j)=bc1*bc2;
              if ((n-i-j)%2==1)
                N(i,j)*=-1;
            }
          }
        }
      }

      template<typename Derived>
      void bezier_T(Eigen::MatrixBase<Derived> &T, const typename Derived::Scalar &t, const typename Derived::Index &n)
      {
        // check to make sure have valid curve
        assert(n>=0);

        typename Derived::Index j, jj;

        T.derived().resize(n+1);
        T.fill(1);
        for (j=0; j<=n-1; ++j)
          for (jj=0; jj<j+1; ++jj)
            T(jj)*=t;
      }

      template<typename Derived>
      void bezier_T_p(Eigen::MatrixBase<Derived> &Tp, const typename Derived::Scalar &t, const typename Derived::Index &n)
      {
        // check to make sure have valid curve
        assert(n>=0);

        typename Derived::Index j, jj;

        Tp.derived().resize(n+1);
        Tp.setZero();
        for (j=0; j<=n-1; ++j)
          Tp(j)=static_cast<typename Derived::Scalar>(n-j);

        for (j=0; j<=n-2; ++j)
          for (jj=0; jj<j+1; ++jj)
            Tp(jj)*=t;
      }

      template<typename Derived>
      void bezier_T_pp(Eigen::MatrixBase<Derived> &Tpp, const typename Derived::Scalar &t, const typename Derived::Index &n)
      {
        // check to make sure have valid curve
        assert(n>=0);

        typename Derived::Index j, jj;

        Tpp.derived().resize(n+1);
        Tpp.setZero();
        for (j=0; j<=n-2; ++j)
          Tpp(j)=static_cast<typename Derived::Scalar>(n-j)*(n-j-1);

        for (j=0; j<=n-3; ++j)
          for (jj=0; jj<j+1; ++jj)
            Tpp(jj)*=t;
      }

      template<typename Derived>
      void bezier_T_ppp(Eigen::MatrixBase<Derived> &Tppp, const typename Derived::Scalar &t, const typename Derived::Index &n)
      {
        // check to make sure have valid curve
        assert(n>=0);

        typename Derived::Index j, jj;

        Tppp.derived().resize(n+1);
        Tppp.setZero();
        for (j=0; j<=n-3; ++j)
          Tppp(j)=static_cast<typename Derived::Scalar>(n-j)*(n-j-1)*(n-j-2);

        for (j=0; j<=n-4; ++j)
          for (jj=0; jj<j+1; ++jj)
            Tppp(jj)*=t;
      }

    }
  }
}

#endif
