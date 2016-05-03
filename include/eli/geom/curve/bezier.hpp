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

#ifndef eli_geom_curve_bezier_h
#define eli_geom_curve_bezier_h

#include <iostream>
#include <vector>

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/mutil/opt/least_squares.hpp"
#include "eli/mutil/dm/binomial_coefficient.hpp"

#include "eli/geom/utility/bezier.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/general/continuity.hpp"
#include "eli/geom/general/bounding_box.hpp"
#include "eli/geom/curve/fit_container.hpp"
#include "eli/geom/intersect/minimum_distance_curve.hpp"

// TODO: MOVE THESE TO utility namespace
namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace internal
      {
        template<typename Derived1, typename Derived2, typename PointType>
        void build_fit_Ab(Eigen::MatrixBase<Derived1> &A,
                          Eigen::MatrixBase<Derived2> &b,
                          std::vector<typename Derived1::Scalar> &t,
                          const std::vector<PointType, Eigen::aligned_allocator<PointType> > &pts,
                          const typename Derived1::Index &n, const size_t &dim)
        {
          typedef Eigen::Matrix<typename Derived1::Scalar, Eigen::Dynamic, Eigen::Dynamic> mat_type;
          typedef Eigen::Matrix<typename Derived1::Scalar, Eigen::Dynamic, 1> col_type;

          typename Derived1::Index i, sz(n+1), npts(pts.size()), nrows(std::max(npts, sz));
          mat_type N, T(nrows,sz);

          // resize the return matrices
          A.derived().resize(nrows, sz);
          b.derived().resize(nrows, dim);

          // build the least squares terms
          eli::geom::utility::bezier_N(N, n);
          for (i=0; i<npts; ++i)
          {
            col_type Tvec;
            eli::geom::utility::bezier_T(Tvec, t[i], n);
            T.row(i)=Tvec.transpose();
            b.row(i)=pts[i];
          }
          for (i=npts; i<nrows; ++i)
          {
            T.row(i).setZero();
            b.row(i).setZero();
          }

          A=T*N;
        }

        template<typename index_type1, typename index_type2, typename index_type3>
        index_type1 determine_n(const index_type1 &deg_in, const index_type2 &nconstrs, const index_type3 &npts)
        {
          index_type1 n;

          // if have an under-determined system then what should be done? For now, reducing
          // the order of resulting curve
          if (deg_in+1>static_cast<index_type1>(npts+nconstrs))
          {
            n=npts+nconstrs-1;
            std::cerr << "deg_in (" << deg_in << ") too low. Order should be " << n << std::endl;
            assert(false);
          }
          else
            n=deg_in;

          return n;
        }
      }
    }
  }
}

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      // TODO: Integrate the tol__ class into this class to replace open_flag and any other place
      //       where numerical error might affect an equivalence comparison
      template<typename data__, unsigned short dim__, typename tol__=eli::util::tolerance<data__> >
      class bezier
      {
        public:
          typedef unsigned short dimension_type;
          typedef data__ data_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef point_type control_point_type;
          typedef typename point_type::Index index_type;
          typedef geom::curve::fit_container<data_type, index_type, dim__, dim__> fit_container_type;
          typedef tol__ tolerance_type;
          typedef Eigen::Matrix<data_type, dim__, dim__> rotation_matrix_type;
          typedef eli::geom::general::bounding_box<data_type, dim__, tolerance_type> bounding_box_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> monomial_coefficient_type;

          typedef bezier<data_type, 1, tolerance_type> onedbezcurve;
          typedef bezier<data_type, 2, tolerance_type> twodbezcurve;
          typedef bezier<data_type, 3, tolerance_type> threedbezcurve;
          typedef bezier<data_type, 4, tolerance_type> fourdbezcurve;

        public:
          bezier() : B(1, dim__), deriv( NULL ) {}
          bezier(const index_type &n) : B((n<=0)?(1):(n+1), dim__), deriv( NULL ) {}

          bezier(const bezier<data_type, dim__, tolerance_type> &bc) : B(bc.B)
          {
            if (bc.deriv)
            {
              deriv = new bezier<data_type, dim__>( *(bc.deriv) );
            }
            else
            {
              deriv = NULL;
            }
          }

          ~bezier()
          {
            invalidate_deriv();
          }

          bezier & operator=(const bezier<data_type, dim__, tolerance_type> &bc)
          {
            if (this != &bc)
            {
              B=bc.B;
              if (bc.deriv)
              {
                deriv = new bezier<data_type, dim__>( *(bc.deriv) );
              }
              else
              {
                deriv = NULL;
              }
            }
            return *this;
          }

          bool operator==(const bezier<data_type, dim__, tolerance_type> &bc) const
          {
            if (this == &bc)
              return true;
            if ((B.rows()!=bc.B.rows()) || (B.cols()!=bc.B.cols()))
              return false;
            if (B!=bc.B)
              return false;
            return true;
          }

          bool operator!=(const bezier<data_type, dim__, tolerance_type> &bc) const
          {
            return !operator==(bc);
          }

          bool approximately_equal(const bezier<data_type, dim__, tolerance_type> &bc) const
          {
            tolerance_type tol;

            if (this==&bc)
              return true;

            if ((B.rows()!=bc.B.rows()) || (B.cols()!=bc.B.cols()))
              return false;

            for (index_type i=0; i<=degree(); ++i)
            {
              if (!tol.approximately_equal(get_control_point(i), bc.get_control_point(i)))
              {
                return false;
              }
            }

            return true;
          }

          bool abouteq(const bezier<data_type, dim__, tolerance_type> &bc, const data_type &ttol2 ) const
          {
            if (this==&bc)
              return true;

            if ((B.rows()!=bc.B.rows()) || (B.cols()!=bc.B.cols()))
              return false;

            for (index_type i=0; i<=degree(); ++i)
            {
              if ( eli::geom::point::distance2( get_control_point(i), bc.get_control_point(i) ) > ttol2 )
              {
                return false;
              }
            }

            return true;
          }

          data_type eqp_distance_bound(const bezier<data_type, dim__, tolerance_type> &bc) const
          {
            if (this==&bc)
              return 0;

            data_type d;
            eli::geom::utility::bezier_eqp_distance_bound(B, bc.B, d);
            return d;
          }

          void clear() {resize(0);}

          void resize(const index_type &t_dim)
          {
            B.resize(t_dim+1, dim__);
            invalidate_deriv();
          }

          index_type degree() const
          {
            return B.rows()-1;
          }

          data_type get_tmax() const {return 1;}
          data_type get_t0() const {return 0;}

          static dimension_type dimension() {return dim__;}

          void set_control_point(const control_point_type &cp, const index_type &i)
          {
            // make sure have valid index
            if (i>degree())
            {
              assert(false);
              return;
            }

            B.row(i)=cp;
            invalidate_deriv();
          }

          control_point_type get_control_point(const index_type &i) const
          {
            // make sure have valid index
            if (i>degree())
            {
              assert(false);
              return B.row(0);
            }

            return B.row(i);
          }

          void get_monomial_coefficients(monomial_coefficient_type &a) const
          {
            a.resize(degree()+1, dim__);
            eli::geom::utility::bezier_control_points_to_monomial(a, B);
          }

          void reflect_xy()
          {
            B.col(2)=-B.col(2);
            invalidate_deriv();
          }

          void reflect_xz()
          {
            B.col(1)=-B.col(1);
            invalidate_deriv();
          }

          void reflect_yz()
          {
            B.col(0)=-B.col(0);
            invalidate_deriv();
          }

          void reflect(const point_type &normal)
          {
            point_type n(normal);

            n.normalize();
            B=B-2*(B*n.transpose())*n;
            invalidate_deriv();
          }

          void reflect(const point_type &normal, const data_type &d)
          {
            point_type n(normal);

            n.normalize();
            B=B-2*(B*n.transpose()-d*Eigen::Matrix<data_type, Eigen::Dynamic, 1>::Ones(degree()+1, 1))*n;
            invalidate_deriv();
          }

          void reverse()
          {
            index_type i, n(degree());
            control_point_matrix_type B_new(n+1, dim__);

            for (i=0; i<=n; ++i)
            {
              B_new.row(n-i)=B.row(i);
            }

            // set the new control points
            B=B_new;
            invalidate_deriv();
          }

          void get_bounding_box(bounding_box_type &bb) const
          {
            index_type i, deg(degree());

            bb.clear();
            for (i=0; i<=deg; ++i)
            {
              bb.add(B.row(i));
            }
          }

          void rotate(const rotation_matrix_type &rmat)
          {
            B*=rmat.transpose();
            invalidate_deriv();
          }

          void rotate(const rotation_matrix_type &rmat, const point_type &rorig)
          {
            translate(-rorig);
            rotate(rmat);
            translate(rorig);
          }

          void translate(const point_type &trans)
          {
            index_type i, deg(degree());
            for (i=0; i<=deg; ++i)
            {
              B.row(i)+=trans;
            }
            invalidate_deriv();
          }

          void scale(const data_type &s)
          {
            index_type i, deg(degree());
            for (i=0; i<=deg; ++i)
            {
              B.row(i)*=s;
            }
            invalidate_deriv();
          }

          bool open() const {return !closed();}
          bool closed() const
          {
            tolerance_type tol;

            for (index_type i=0; i<dim__; ++i)
            {
              if (!tol.approximately_equal(B(0, i), B(degree(), i)))
                return false;
            }
            return true;
          }

          point_type f(const data_type &t) const
          {
            // check to make sure have valid curve
            assert(degree()>=0);

            // check to make sure given valid parametric value
            assert((t>=0) && (t<=1));

            // short circuit if degree not high enough
            if (degree()==0)
            {
              return B.row(0);
            }

            point_type rtn;
            eli::geom::utility::de_casteljau(rtn, B, t);

            return rtn;
          }

          point_type fp(const data_type &t) const
          {
            // check to make sure have valid curve
            assert(degree()>=0);

            // check to make sure given valid parametric value
            assert((t>=0) && (t<=1));

            // short circuit if degree not high enough
            if (degree()<1)
            {
              point_type rtn;
              rtn.setZero();
              return rtn;
            }

            validate_deriv();

            return deriv->f( t );
          }

          void fp(bezier<data_type, dim__> &bc_fp) const
          {
            // check to make sure have valid curve
            assert(degree()>=0);

            bc_fp.resize( degree()-1 );

            eli::geom::utility::bezier_p_control_point(bc_fp.B, B);
          }

          point_type fpp(const data_type &t) const
          {
            // check to make sure have valid curve
            assert(degree()>=0);

            // check to make sure given valid parametric value
            assert((t>=0) && (t<=1));

            point_type rtn;

            // short circuit if degree not high enough
            if (degree()<2)
            {
              rtn.setZero();
              return rtn;
            }

            validate_deriv();

            return deriv->fp( t );
          }

          point_type fppp(const data_type &t) const
          {
            // check to make sure have valid curve
            assert(degree()>=0);

            // check to make sure given valid parametric value
            assert((t>=0) && (t<=1));

            point_type rtn;

            // short circuit if degree not high enough
            if (degree()<3)
            {
              rtn.setZero();
              return rtn;
            }

            validate_deriv();

            return deriv->fpp( t );
          }

          point_type tangent(const data_type &t) const
          {
            point_type tgt(fp(t));

            tgt.normalize();
            return tgt;
          }

          void frenet_serret_frame(point_type &t, point_type &n, point_type &b, const data_type &t0)
          {
            t=tangent(t0);
            b=fp(t0).cross(fpp(t0));
            n=-t.cross(b)/b.norm();
            b.normalize();
            n.normalize();
          }

          void degree_promote()
          {
            // create vector of new control points
            control_point_matrix_type B_new(degree()+2, dim__);

            // build the new control points
            eli::geom::utility::bezier_promote_control_points(B_new, B);

            // set the new control points
            B=B_new;
            invalidate_deriv();
          }

          void degree_promote_to(const index_type target_degree)
          {
            // create vector of new control points
            control_point_matrix_type B_new(target_degree+1, dim__);

            // build the new control points
            eli::geom::utility::bezier_promote_control_points_to(B_new, B);

            // set the new control points
            B=B_new;
            invalidate_deriv();
          }

          bool degree_demote(const geom::general::continuity &continuity_degree=geom::general::C0)
          {
            // check if can demote
            int ncon(0);
            switch(continuity_degree)
            {
              case(eli::geom::general::NOT_CONNECTED):
                ncon=0;
                break;
              case(eli::geom::general::C0):
                ncon=2;
                break;
              case(eli::geom::general::C1):
                ncon=4;
                break;
              case(eli::geom::general::C2):
                ncon=6;
                break;
              default:
                ncon=-1;
            }
            if (ncon<0)
              return false;
            if (ncon>degree()-2)
              return false;

            // demote control points and set them
            control_point_matrix_type B_new(degree(), dim__);
            eli::geom::utility::bezier_demote_control_points(B_new, B, ncon);
            B=B_new;
            invalidate_deriv();

            return true;
          }

          void degree_to_cubic()
          {
              // allocate control points and set them
              control_point_matrix_type B_new(4, dim__);
              eli::geom::utility::bezier_control_points_to_cubic(B_new, B);
              B=B_new;
              invalidate_deriv();
          }

          void split(bezier<data_type, dim__> &bc_l, bezier<data_type, dim__> &bc_r, const data_type &t0) const
          {
            if ( (t0>1) || (t0<0) )
            {
              assert(false);
              return;
            }

            control_point_matrix_type bl(degree()+1, dim__), br(degree()+1, dim__);
            index_type n(degree());

            // resize the curves
            bc_l.resize(n);
            bc_r.resize(n);

            eli::geom::utility::bezier_split_control_points(bl, br, B, t0);

            // set the control points
            for (index_type i=0; i<=n; ++i)
            {
              bc_l.set_control_point(bl.row(i), i);
              bc_r.set_control_point(br.row(i), i);
            }
          }

          void fit(const fit_container_type &fcon, const index_type &deg_in)
          {
            std::vector<data_type> t;
            fit_only(t, fcon, deg_in);
          }

          void fit(std::vector<data_type> &t, const fit_container_type &fcon, const index_type &deg_in)
          {
            fit_with_error(t, fcon, deg_in);
          }

          data_type fit_with_error(const fit_container_type &fcon, const index_type &deg_in)
          {
            std::vector<data_type> t;
            return fit_with_error(t, fcon, deg_in);
          }

          data_type fit_with_error(std::vector<data_type> &t, const fit_container_type &fcon, const index_type &deg_in)
          {
            fit_only(t, fcon, deg_in);
            return est_fit_error(t, fcon);
          }

          void interpolate(const fit_container_type &fcon)
          {
            std::vector<data_type> t;
            interpolate(t, fcon);
          }

          void interpolate(std::vector<data_type> &t, const fit_container_type &fcon)
          {
            size_t i, npts(fcon.number_points()), n, nclosed, ai;
            std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(npts);

            // get the points from the container
            fcon.get_points(pts.begin());

            // account for the closed constraints
            switch(fcon.get_end_flag())
            {
              case(eli::geom::general::C2):
              {
                nclosed=3;
                break;
              }
              case(eli::geom::general::C1):
              {
                nclosed=2;
                break;
              }
              case(eli::geom::general::C0):
              {
                nclosed=1;
                break;
              }
              default:
              {
                nclosed=0;
                break;
              }
            }

            // calculate the t corresponding to input points via approximate arc-length
            determine_t(t, pts, fcon.closed());

            // determine the actual degree of curve
            n=fcon.number_constraints(false)+nclosed+npts-1;

            // build the fit terms from points
            mat_type A, x;
            row_pts_type b;
            internal::build_fit_Ab(A, b, t, pts, n, dim__);

            size_t nconpts=fcon.number_constraint_points();
            std::vector<typename fit_container_type::index_type> indexes(nconpts);

            // cycle through all of the constraints
            fcon.get_constraint_indexes(indexes.begin());
            ai=pts.size();
            for (size_t i=0; i<nconpts; ++i)
            {
              point_type pt;
              typename fit_container_type::constraint_info ci;
              typename fit_container_type::error_code ec;

              ec=fcon.get_constraint(indexes[i], ci);
              if (ec!=fit_container_type::NO_ERRORS)
              {
                assert(false);
              }
              else
              {
                // set the C1 constraint
                if (ci.using_fp()!=fit_container_type::constraint_info::NOT_USED)
                {
                  col_type Tp;
                  mat_type N;
                  eli::geom::utility::bezier_N(N, n);

                  eli::geom::utility::bezier_T_p(Tp, t[indexes[i]], n);
                  A.row(ai)=Tp.transpose()*N;
                  b.row(ai)=ci.get_fp();
                  ++ai;
                }

                // set the C2 constraint
                if (ci.using_fpp()!=fit_container_type::constraint_info::NOT_USED)
                {
                  col_type Tpp;
                  mat_type N;
                  eli::geom::utility::bezier_N(N, n);

                  eli::geom::utility::bezier_T_pp(Tpp, t[indexes[i]], n);
                  A.row(ai)=Tpp.transpose()*N;
                  b.row(ai)=ci.get_fpp();
                  ++ai;
                }
              }
            }

            // add the closed constraint if needed
            switch(fcon.get_end_flag())
            {
              case(eli::geom::general::C2):
              {
                A(ai,2)=1.0;
                A(ai,1)=-2.0;
                A(ai,0)=1.0;
                A(ai,n)=-1.0;
                A(ai,n-1)=2.0;
                A(ai,n-2)=-1.0;
                b.row(ai).setZero();
                ++ai;
              }
              case(eli::geom::general::C1):
              {
                A(ai,1)=1.0;
                A(ai,0)=-1.0;
                A(ai,n)=-1.0;
                A(ai,n-1)=1.0;
                ++ai;
              }
              case(eli::geom::general::C0):
              {
                A(ai,0)=1.0;
                A(ai,n)=-1.0;
                ++ai;
              }
              default:
              {
                // no need to do anything
                assert(ai==n+1);
              }
            }

            // solve for the control points
            x=A.lu().solve(b);

            // extract the control points and set them
            control_point_matrix_type ctrl(n+1, dim__);
            for (i=0; i<n+1; ++i)
              ctrl.row(i)=x.row(i);

            // ensure that the last control point and first are the same
            if (fcon.closed())
              ctrl.row(n)=ctrl.row(0);

            B=ctrl;
            invalidate_deriv();
          }

          void product( const bezier<data_type, dim__> &a, const bezier<data_type, dim__> &b)
          {
            assert( a.B.cols() == dim__ );
            assert( b.B.cols() == dim__ );

            index_type m( a.degree() ), n( b.degree() );
            control_point_matrix_type scaleda, scaledb, scaledc;

            scaleda.resize( m + 1, dim__ );
            eli::geom::utility::bezier_control_points_to_scaled_bezier( scaleda, a.B );

            scaledb.resize( n + 1, dim__ );
            eli::geom::utility::bezier_control_points_to_scaled_bezier( scaledb, b.B );

            scaledc.resize( m + n + 1, dim__ );
            scaledc.setZero();
            eli::geom::utility::multiply_scaled_bezier( scaledc, scaleda, scaledb );

            resize( m + n );
            eli::geom::utility::scaled_bezier_to_control_points_bezier( B, scaledc );
            invalidate_deriv();
          }

          void product1d( const bezier<data_type, dim__> &a, const bezier<data_type, 1> &b)
          {
            assert( a.B.cols() == dim__ );
            assert( b.B.cols() == 1 );

            index_type m( a.degree() ), n( b.degree() );
            control_point_matrix_type scaleda, scaledc;
            oned_control_point_matrix_type scaledb;

            scaleda.resize( m + 1, dim__ );
            eli::geom::utility::bezier_control_points_to_scaled_bezier( scaleda, a.B );

            scaledb.resize( n + 1, 1 );
            eli::geom::utility::bezier_control_points_to_scaled_bezier( scaledb, b.B );

            scaledc.resize( m + n + 1, dim__ );
            scaledc.setZero();
            eli::geom::utility::multiply_scaled_bezier1d( scaledc, scaleda, scaledb );

            resize( m + n );
            eli::geom::utility::scaled_bezier_to_control_points_bezier( B, scaledc );
            invalidate_deriv();
          }

          void square( const bezier<data_type, dim__> &a )
          {
            assert( a.B.cols() == dim__ );

            index_type m( a.degree() );
            control_point_matrix_type scaleda, scaledc;

            scaleda.resize( m + 1, dim__ );
            eli::geom::utility::bezier_control_points_to_scaled_bezier( scaleda, a.B );

            scaledc.resize( m + m + 1, dim__ );
            scaledc.setZero();
            eli::geom::utility::multiply_scaled_bezier( scaledc, scaleda, scaleda );

            resize( m + m );
            eli::geom::utility::scaled_bezier_to_control_points_bezier( B, scaledc );
            invalidate_deriv();
          }

          void sqrt( const bezier<data_type, dim__> &a )
          {
            typedef bezier<data_type, dim__> curve_type;

            curve_type ca( a );

            assert( ca.B.cols() == dim__ );

            index_type m( ca.degree() );

            if ( m % 2 ) // Promote once if odd
            {
              control_point_matrix_type a_new( m + 2, dim__);
              eli::geom::utility::bezier_promote_control_points(a_new, ca.B);
              ca.B=a_new;
              m++;
            }

            index_type n( m/2 );

            control_point_matrix_type scaleda, scaledc;

            scaleda.resize( m + 1, dim__ );
            eli::geom::utility::bezier_control_points_to_scaled_bezier( scaleda, ca.B );

            scaledc.resize( n + 1, dim__ );
            scaledc.setZero();
            eli::geom::utility::sqrt_scaled_bezier( scaledc, scaleda );

            resize( n );
            eli::geom::utility::scaled_bezier_to_control_points_bezier( B, scaledc );
            invalidate_deriv();
          }

          void sum( const bezier<data_type, dim__> &a, const bezier<data_type, dim__> &b)
          {
            typedef bezier<data_type, dim__> curve_type;

            curve_type ca( a );
            curve_type cb( b );

            index_type n;
            n = std::max( ca.degree(), cb.degree() );

            ca.degree_promote_to( n );
            cb.degree_promote_to( n );

            B = ca.B + cb.B;

            invalidate_deriv();
          }

          onedbezcurve sumcompcurve() const
          {
            onedbezcurve retcurve;

            index_type n(degree()), i, j;

            retcurve.resize(n);
            for (i=0; i<=n; ++i)
            {
              data_type d = 0;
              point_type p = get_control_point( i );
              for (j=0; j<dim__; ++j)
              {
                d += p(j);
              }
              typename onedbezcurve::point_type pd;
              pd(0) = d;

              retcurve.set_control_point( pd, i );
            }

            return retcurve;
          }

          onedbezcurve mindistcurve( const point_type & pt ) const
          {
            onedbezcurve retcurve;
            typedef bezier<data_type, dim__> curve_type;

            curve_type c(*this);

            c.translate( -pt );

            validate_deriv();

            curve_type prod;

            prod.product( *deriv, c );

            onedbezcurve dot;
            dot = prod.sumcompcurve();

            retcurve.product( dot, dot );

            return retcurve;
          }

          bool allpos( const data_type &smallpos ) const
          {
            index_type i, j;
            index_type n(degree());

            for (i=0; i<=n; ++i)
            {
              point_type p = get_control_point( i );

              for (j=0; j<dim__; ++j)
              {
                if ( p(j) <= smallpos )
                {
                  return false;
                }
              }
            }
            return true;
          }

        private:
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_point_matrix_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> oned_control_point_matrix_type;

          typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> row_pts_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> col_type;
          typedef Eigen::Matrix<data_type, 1, Eigen::Dynamic> row_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> mat_type;

        private:
          control_point_matrix_type B;      /** control polygon coordinates */
          mutable bezier<data_type, dim__> * deriv;

        private:
          void determine_t(std::vector<data_type> &t, const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts, bool closed) const
          {
            index_type i, npts(pts.size());
            data_type len, small_dist(std::numeric_limits<data_type>::epsilon());

            // calculate the corresponding t values for each point via approximate arc-length
            t.resize(npts);
            t[0]=0.0;
            for (i=1; i<npts; ++i)
            {
              data_type temp;

              temp=geom::point::distance(pts[i-1], pts[i]);

              // in case two adjacent points are given too close together
              if (temp<small_dist)
                temp=small_dist;

              t[i]=t[i-1]+temp;
            }
            len=t[npts-1];
            if (closed)
            {
              len+=geom::point::distance(pts[npts-1], pts[0]);
            }

            for (i=0; i<npts; ++i)
              t[i]/=len;
          }

          // This method is private because the t value returned in the first parameter is not the result of a nearest
          // neighbor fit as would be expected from the fit_with_error methods.
          void fit_only(std::vector<data_type> &t, const fit_container_type &fcon, const index_type &deg_in)
          {
            size_t i, npts(fcon.number_points()), n, nclosed;
            std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(npts);

            // get the points from the container
            fcon.get_points(pts.begin());

            // account for the closed constraints
            switch(fcon.get_end_flag())
            {
              case(eli::geom::general::C2):
              {
                nclosed=3;
                break;
              }
              case(eli::geom::general::C1):
              {
                nclosed=2;
                break;
              }
              case(eli::geom::general::C0):
              {
                nclosed=1;
                break;
              }
              default:
              {
                nclosed=0;
                break;
              }
            }

            // calculate the t corresponding to input points via approximate arc-length
            determine_t(t, pts, fcon.closed());

            // determine the actual degree of curve
            n=internal::determine_n(deg_in, fcon.number_constraints()+nclosed, npts);

            // build the fit terms from points
            mat_type A, x;
            row_pts_type b;
            internal::build_fit_Ab(A, b, t, pts, n, dim__);

            // handle special case of unconstrained optimization problem
            if ((fcon.number_constraints()==0) && (fcon.open()))
            {
              // determine the coefficients
              eli::mutil::opt::least_squares_uncon(x, A, b);
            }
            else
            {
              // now becomes a constrained least squares problem
              // the constraints come from closed flag and/or constraint collection
              size_t bi, ncon=fcon.number_constraints()+nclosed;

              mat_type B(ncon, n+1);
              row_pts_type d(ncon, dim__);

              // construct the system of constraints
              B.setZero();
              d.setZero();

              size_t nconpts=fcon.number_constraint_points();
              std::vector<typename fit_container_type::index_type> indexes(nconpts);

              // cycle through all of the constraints
              fcon.get_constraint_indexes(indexes.begin());
              bi=0;
              for (size_t i=0; i<nconpts; ++i)
              {
                point_type pt;
                typename fit_container_type::constraint_info ci;
                typename fit_container_type::error_code ec;

                ec=fcon.get_constraint(indexes[i], ci);
                if (ec!=fit_container_type::NO_ERRORS)
                {
                  assert(false);
                }
                else
                {
                  // set the C0 constraint
                  col_type T;
                  mat_type N;
                  eli::geom::utility::bezier_T(T, t[indexes[i]], n);
                  eli::geom::utility::bezier_N(N, n);
                  B.row(bi)=T.transpose()*N;
                  d.row(bi)=pts[indexes[i]];
                  ++bi;

                  // set the C1 constraint
                  if (ci.using_fp()!=fit_container_type::constraint_info::NOT_USED)
                  {
                    col_type Tp;

                    eli::geom::utility::bezier_T_p(Tp, t[indexes[i]], n);
                    B.row(bi)=Tp.transpose()*N;
                    d.row(bi)=ci.get_fp();
                    ++bi;
                  }

                  // set the C2 constraint
                  if (ci.using_fpp()!=fit_container_type::constraint_info::NOT_USED)
                  {
                    col_type Tpp;

                    eli::geom::utility::bezier_T_pp(Tpp, t[indexes[i]], n);
                    B.row(bi)=Tpp.transpose()*N;
                    d.row(bi)=ci.get_fpp();
                    ++bi;
                  }
                }
              }

              // add the closed constraint if needed
              switch(fcon.get_end_flag())
              {
                case(eli::geom::general::C2):
                {
                  B(bi,2)=1.0;
                  B(bi,1)=-2.0;
                  B(bi,0)=1.0;
                  B(bi,n)=-1.0;
                  B(bi,n-1)=2.0;
                  B(bi,n-2)=-1.0;
                  ++bi;
                }
                case(eli::geom::general::C1):
                {
                  B(bi,1)=1.0;
                  B(bi,0)=-1.0;
                  B(bi,n)=-1.0;
                  B(bi,n-1)=1.0;
                  ++bi;
                }
                case(eli::geom::general::C0):
                {
                  B(bi,0)=1.0;
                  B(bi,n)=-1.0;
                  ++bi;
                }
                default:
                {
                  // no need to do anything
                  assert(bi==ncon);
                }
              }

              // determine the coefficients
              eli::mutil::opt::least_squares_eqcon(x, A, b, B, d);
            }

            // extract the control points and set them
            control_point_matrix_type ctrl(n+1, dim__);
            for (i=0; i<n+1; ++i)
              ctrl.row(i)=x.row(i);

            // ensure that the last control point and first are the same
            if (fcon.closed())
              ctrl.row(n)=ctrl.row(0);

            B=ctrl;
            invalidate_deriv();
          }

          data_type est_fit_error(std::vector<data_type> &t, const fit_container_type &fcon)
          {
            size_t i, npts(fcon.number_points());
            std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(npts);

            // get the points from the container
            fcon.get_points(pts.begin());

            // calculate the error at the point
            data_type err(0);
            for (i=0; i<pts.size(); ++i)
            {
              err+=eli::geom::intersect::minimum_distance(t[i], *this, pts[i]);
            }

            return err;
          }

          void invalidate_deriv()
          {
            if ( deriv )
            {
              delete deriv;
              deriv = NULL;
            }
          }

          void validate_deriv() const
          {
            if ( !deriv )
            {
              deriv = new bezier<data_type, dim__>();
              fp( *deriv );
            }
          }

      };

      typedef bezier<float, 1> bezier1f;
      typedef bezier<float, 2> bezier2f;
      typedef bezier<float, 3> bezier3f;
      typedef bezier<double, 1> bezier1d;
      typedef bezier<double, 2> bezier2d;
      typedef bezier<double, 3> bezier3d;
      typedef bezier<long double, 1> bezier1ld;
      typedef bezier<long double, 2> bezier2ld;
      typedef bezier<long double, 3> bezier3ld;
    }
  }
}
#endif
