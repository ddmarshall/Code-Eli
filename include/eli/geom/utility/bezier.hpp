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

#ifndef eli_geom_utility_bezier_hpp
#define eli_geom_utility_bezier_hpp

#include "eli/code_eli.hpp"

#include "eli/mutil/dm/binomial_coefficient.hpp"

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
        typename Derived2::Scalar oneminust(1-t);
        typename Derived2::Index i, k;

        for (k=1; k<Q.rows(); ++k)
        {
          for (i=0; i<Q.rows()-k; ++i)
          {
            Q.row(i)=oneminust*Q.row(i)+t*Q.row(i+1);
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
          cp_p.row(i)=static_cast<typename Derived2::Scalar>(n)*(cp.row(i+1)-cp.row(i));
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
          cp_pp.row(i)=static_cast<typename Derived2::Scalar>(n)*(n-1)*(cp.row(i+2)-2*cp.row(i+1)+cp.row(i));
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
          cp_ppp.row(i)=static_cast<typename Derived2::Scalar>(n)*(n-1)*(n-2)*(cp.row(i+3)-3*cp.row(i+2)+3*cp.row(i+1)-cp.row(i));
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

      template<typename Derived1, typename Derived2>
      void bezier_promote_control_points_to(Eigen::MatrixBase<Derived1> &cp_out, const Eigen::MatrixBase<Derived2> &cp_in)
      {
        typedef typename Derived1::Index index_type;
        typedef typename Derived1::Scalar data_type;

        // do some dimension checks
        assert(cp_out.rows()>=cp_in.rows());
        assert(cp_out.cols()==cp_in.cols());

        index_type i, ntarget(cp_out.rows()-1), nstart(cp_in.rows()-1);
        index_type n(nstart);

        // Make in-place copy of control points.
        for (i=0; i<n+1; ++i)
          cp_out.row(i)=cp_in.row(i);

        for (; n<ntarget; n++)
        {
          // Assign final value, n'th value no longer needed and can be replaced.
          cp_out.row(n+1)=cp_out.row(n);
          // Work backwards, calculating in-place.
          for (i=n; i>0; --i)
            cp_out.row(i)=(cp_out.row(i-1)-cp_out.row(i))*(static_cast<data_type>(i)/(n+1))+cp_out.row(i);
        }
      }

      // Calculate cubic 'equivalent' to bezier curve.  If input curve has degree less than cubic,
      // the curve is promoted to an exactly equivalent cubic curve.  If the input curve is cubic,
      // the output curve is copied from the input.  If the input curve is higher order, the output
      // curve is a cubic bezier curve that will match the input curve in endpoint position and
      // slopes.  This approximation is used because it is very simple and fast.  This curve may be
      // a terrible approximation of the original curve, but through successive subdivision, a
      // piecewise cubic approximation of any curve to any accuracy should be possible.
      template<typename Derived1, typename Derived2>
      void bezier_control_points_to_cubic(Eigen::MatrixBase<Derived1> &cp_out, const Eigen::MatrixBase<Derived2> &cp_in)
      {
        typedef typename Derived1::Index index_type;
        typedef typename Derived1::Scalar data_type;

        // do some dimension checks
        assert(cp_out.rows() == 4);
        assert(cp_out.cols() == cp_in.cols());

        if(cp_in.rows() < 4) // Promote
        {
          bezier_promote_control_points_to(cp_out, cp_in);
        }
        else if(cp_in.rows() == 4) // Do nothing
        {
          index_type i;
          for( i=0; i<4; ++i)
            cp_out.row(i) = cp_in.row(i);
        }
        else // C1 Cubic demote
        {
          index_type n(cp_in.rows()-1);
          data_type s = static_cast<data_type>(n)/static_cast<data_type>(3);

          cp_out.row(0) = cp_in.row(0);
          cp_out.row(1) = cp_in.row(0) + s * (cp_in.row(1)-cp_in.row(0));
          cp_out.row(2) = cp_in.row(n) + s * (cp_in.row(n-1)-cp_in.row(n));
          cp_out.row(3) = cp_in.row(n);
        }
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
          data_type coef(std::pow(2.0, 1-2*n)), sum(1), tmp;
          lambda(0)=coef*sum;
          for (i=1; i<n; ++i)
          {
            eli::mutil::dm::n_choose_k(tmp, 2*n, 2*i);
            sum+=tmp;
            lambda(i)=coef*sum;
          }
        }
        else
        {
          // interpolating end points up-to alpha order derivative
          data_type tmp, coef, sum(0);
          index_type alpha(ncon/2);
          eli::mutil::dm::n_choose_k(tmp, 2*n, n+2*alpha);
          coef=1/tmp;
          lambda(0)=coef*sum;
          for (i=1; i<n; ++i)
          {
            if ( (i<alpha) || (i+alpha>n) )
              sum+=0;
            else
            {
              data_type tmp1, tmp2;
              eli::mutil::dm::n_choose_k(tmp1, n, i-alpha);
              eli::mutil::dm::n_choose_k(tmp2, n, i+alpha);
              sum+=tmp1*tmp2;
            }
            lambda(i)=coef*sum;
          }
        }

        // calculate new control points
        for (i=0; i<n; ++i)
          cp_out.row(i)=B_I.row(i)*(1-lambda(i))+B_II.row(i)*lambda(i);
      }

      // Calculate upper bound of equiparametric distance between bezier curves.
      // This routine efficiently calculates an upper bound of the distance between
      // equal parameter points on two bezier curves.  While not useful for
      // calculating distances between arbitrary curves, this is just the ticket
      // for judging the quality of one curve meant to approximate another.
      // This algorithm requires that both curves have identical order, so the
      // lower order curve is first promoted to match the higher order curve.
      // Then, the control points of the curves are differenced, resulting in
      // the control points of the curve that is the equiparametric difference
      // between the curves.  Since the control points provide a convex hull of
      // a curve, they represent the worst-case difference between curves.  The
      // maximum distance between the origin and a control point of this curve
      // is an upper bound of the equipmarametric distance between the input
      // curves.
      template<typename Derived1, typename Derived2>
      void bezier_eqp_distance_bound(const Eigen::MatrixBase<Derived1> &cp_a, const Eigen::MatrixBase<Derived2> &cp_b, typename Derived1::Scalar &maxd)
      {
        typedef typename Derived1::Index index_type;
        typedef typename Derived1::Scalar data_type;

        // dimension check
        assert(cp_a.cols()==cp_b.cols());

        index_type na(cp_a.rows()-1), nb(cp_b.rows()-1);

        // make working copies
        Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> cp_A(cp_a);
        Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> cp_B(cp_b);

        if(na == nb)
        {
        }
        else if(na > nb)
        {
          cp_B=cp_A; // Copy to allocate proper size
          bezier_promote_control_points_to(cp_B, cp_b); // promote b to B
        }
        else
        {
          cp_A=cp_B; // Copy to allocate proper size
          bezier_promote_control_points_to(cp_A, cp_a); // promote a to A
        }

        // calculate control point differences and track farthest from origin
        maxd = (cp_B-cp_A).rowwise().norm().maxCoeff();
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
      void bezier_coefficient_factors(Eigen::MatrixBase<Derived> &co, const typename Derived::Scalar &t, const typename Derived::Index &n)
      {
        typename Derived::Scalar coef(1), k(1), tau(std::pow(1-t, n));
        typename Derived::Index j;

        co.derived().resize(n+1, 1);

        j=0;
        co(j)=coef*k*tau;
        for (j=1; j<=n; ++j)
        {
          k*=static_cast<typename Derived::Scalar>(n-j+1)/j;
          tau*=t/(1-t);
          co(j)=coef*k*tau;
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
              eli::mutil::dm::n_choose_k(bc1, n, j);
              eli::mutil::dm::n_choose_k(bc2, n-j, n-i-j);
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

      template<typename Derived1, typename Derived2>
      void monomial_to_bezier_control_points(Eigen::MatrixBase<Derived1> &cp, const Eigen::MatrixBase<Derived2> &a)
      {
        // do some checks on incoming matrix dimensions
        assert(cp.cols()==a.cols());
        assert(cp.rows()==a.rows());

        typename Derived1::Index i, j, deg(cp.rows()-1);
        typename Derived1::Scalar bc;

        cp.setZero();
        for (j=0; j<=deg; ++j)
        {
          for (i=0; i<=j; ++i)
          {
            eli::mutil::dm::binomial_coefficient(bc, deg-i, j-i);
            cp.row(j)+=bc*a.row(i);
          }
          eli::mutil::dm::binomial_coefficient(bc, deg, j);
          cp.row(j)/=bc;
        }
      }

      template<typename Derived1, typename Derived2>
      void bezier_control_points_to_monomial(Eigen::MatrixBase<Derived1> &a, const Eigen::MatrixBase<Derived2> &cp)
      {
        // do some checks on incoming matrix dimensions
        assert(cp.cols()==a.cols());
        assert(cp.rows()==a.rows());

        typename Derived1::Index i, j, deg(a.rows()-1);
        typename Derived1::Scalar bc1, bc2, sgn;

        a.setZero();
        for (j=0; j<=deg; ++j)
        {
          eli::mutil::dm::binomial_coefficient(bc1, deg, j);
          for (i=0; i<=j; ++i)
          {
            sgn=(((j-i)%2)==0)?(1):(-1);
            eli::mutil::dm::binomial_coefficient(bc2, j, i);
            a.row(j)+=bc1*bc2*sgn*cp.row(i);
          }
        }
      }

      template<typename Derived1, typename Derived2>
      void bezier_control_points_to_scaled_bezier(Eigen::MatrixBase<Derived1> &scp, const Eigen::MatrixBase<Derived2> &cp)
      {
        // dimension check
        assert(cp.cols()==scp.cols());
        assert(cp.rows()==scp.rows());

        typename Derived1::Index i, n(cp.rows()-1);
        typename Derived1::Scalar bc;

        for (i=0; i<=n; ++i)
        {
          eli::mutil::dm::binomial_coefficient(bc, n, i);
          scp.row(i) = cp.row(i) * bc;
        }
      }

      template<typename Derived1>
      void bezier_control_points_to_scaled_bezier(Eigen::MatrixBase<Derived1> &cp)
      {
        typename Derived1::Index i, n(cp.rows()-1);
        typename Derived1::Scalar bc;

        for (i=0; i<=n; ++i)
        {
          eli::mutil::dm::binomial_coefficient(bc, n, i);
          cp.row(i) = cp.row(i) * bc;
        }
      }

      template<typename Derived1, typename Derived2>
      void scaled_bezier_to_control_points_bezier(Eigen::MatrixBase<Derived1> &cp, const Eigen::MatrixBase<Derived2> &scp)
      {
        // dimension check
        assert(scp.cols()==cp.cols());
        assert(scp.rows()==cp.rows());

        typename Derived1::Index i, n(scp.rows()-1);
        typename Derived1::Scalar bc;

        for (i=0; i<=n; ++i)
        {
          eli::mutil::dm::binomial_coefficient(bc, n, i);
          cp.row(i) = scp.row(i) / bc;
        }
      }

      template<typename Derived1>
      void scaled_bezier_to_control_points_bezier(Eigen::MatrixBase<Derived1> &cp)
      {
        typename Derived1::Index i, n(cp.rows()-1);
        typename Derived1::Scalar bc;

        for (i=0; i<=n; ++i)
        {
          eli::mutil::dm::binomial_coefficient(bc, n, i);
          cp.row(i) = cp.row(i) / bc;
        }
      }

      template<typename Derived1, typename Derived2>
      void multiply_scaled_bezier(Eigen::MatrixBase<Derived1> &c, const Eigen::MatrixBase<Derived2> &a, const Eigen::MatrixBase<Derived2> &b)
      {
        typename Derived1::Index i, j, k, m( a.rows() - 1 ), n( b.rows() - 1 );

        // dimension check
        assert( a.cols() == b.cols() );
        assert( a.cols() == c.cols() );

        assert( c.rows() == m + n + 1 );

        for ( j = 0; j <= n; ++j )
        {
          for ( i = 0; i <= m; i++ )
          {
            k = i + j;
            c.row(k) = c.row(k) + a.row(i).cwiseProduct(b.row(j));
          }
        }
      }

      template<typename Derived1, typename Derived2, typename Derived3>
      void multiply_scaled_bezier1d(Eigen::MatrixBase<Derived1> &c, const Eigen::MatrixBase<Derived2> &a, const Eigen::MatrixBase<Derived3> &b)
      {
        typename Derived1::Index i, j, k, m( a.rows() - 1 ), n( b.rows() - 1 );

        // dimension check
        assert( b.cols() == 1 );
        assert( a.cols() == c.cols() );

        assert( c.rows() == m + n + 1 );

        for ( j = 0; j <= n; ++j )
        {
          for ( i = 0; i <= m; i++ )
          {
            k = i + j;
            c.row(k) = c.row(k) + b.row(j).col(0) * a.row(i);
          }
        }
      }

      template<typename Derived1, typename Derived2>
      void sqrt_scaled_bezier(Eigen::MatrixBase<Derived1> &c, const Eigen::MatrixBase<Derived2> &a)
      {
        typename Derived1::Index i, j, k, m( c.rows() - 1 );

        // dimension check
        assert( a.cols() == c.cols() );
        assert( a.rows() == m + m + 1 );

        c.row(0) = a.row(0).cwiseSqrt();
        for ( k = 1; k <= m; ++k )
        {
          c.row(k) = a.row(k); // Initialize to compute coeff in place.
          for ( j = 1; j < k; ++j ) // Coeffs before k are known.  Start subtraction at 1.
          {
            i = k - j; // Only interested in k = i + j boundary case.
            c.row(k) = c.row(k) - c.row(j).cwiseProduct(c.row(i));  // Subtract off known terms
          }

          c.row(k) = ( 0.5 * c.row(k) ).cwiseQuotient( c.row(0) );  // Divide off known zero coeff.
        }
      }

    }
  }
}

#endif
