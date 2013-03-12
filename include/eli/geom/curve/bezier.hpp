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

#ifndef geom_curve_bezier_h
#define geom_curve_bezier_h

#include <cassert>
#include <iostream>

#include <Eigen/Eigen>

#include "eli/util/tolerance.hpp"

#include "eli/mutil/opt/least_squares.hpp"
#include "eli/mutil/dm/binomial_coefficient.hpp"

#include "eli/geom/utility/bezier.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/general/continuity.hpp"
#include "eli/geom/curve/fit_container.hpp"

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

        public:
          bezier() : B(1, dim__) {}
          bezier(const index_type &n) : B((n<=0)?(1):(n+1), dim__)
          {
          }
          bezier(const bezier<data_type, dim__> &bc) : B(bc.B) {}
          ~bezier() {}

          bezier & operator=(const bezier<data_type, dim__> &bc)
          {
            if (this != &bc)
            {
              B=bc.B;
            }
            return *this;
          }

          bool operator==(const bezier<data_type, dim__> &bc) const
          {
            if (this == &bc)
              return true;
            if ((B.rows()!=B.rows()) || (B.cols()!=B.cols()))
              return false;
            if (B!=bc.B)
              return false;
            return true;
          }

          bool operator!=(const bezier<data_type, dim__> &bc) const
          {
            return !operator==(bc);
          }

          void resize(const index_type &t_dim)
          {
            B.resize(t_dim+1, dim__);
          }

          index_type degree() const
          {
            return B.rows()-1;
          }

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
          }

          void get_bounding_box(point_type &pmin, point_type &pmax) const
          {
            index_type i, k, deg(degree());
            point_type tmp;

            pmin=B.row(0);
            pmax=pmin;
            for (i=0; i<=deg; ++i)
            {
              tmp=B.row(i);
              for (k=0; k<dim__; ++k)
              {
                if (tmp(k)<pmin(k))
                {
                  pmin(k)=tmp(k);
                }
                if (tmp(k)>pmax(k))
                {
                  pmax(k)=tmp(k);
                }
              }
            }
          }

          void rotate(const rotation_matrix_type &rmat)
          {
            B*=rmat.transpose();
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

            point_type rtn;

            // short circuit if degree not high enough
            if (degree()<1)
            {
              rtn.setZero();
              return rtn;
            }

            control_point_matrix_type B_p(degree()+1-1, dim__);

            eli::geom::utility::bezier_p_control_point(B_p, B);
            eli::geom::utility::de_casteljau(rtn, B_p, t);

            return rtn;
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

            control_point_matrix_type B_pp(degree()+1-2, dim__);

            eli::geom::utility::bezier_pp_control_point(B_pp, B);
            eli::geom::utility::de_casteljau(rtn, B_pp, t);

            return rtn;
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

            control_point_matrix_type B_ppp(degree()+1-3, dim__);

            eli::geom::utility::bezier_ppp_control_point(B_ppp, B);
            eli::geom::utility::de_casteljau(rtn, B_ppp, t);

            return rtn;
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

            return true;
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

          data_type fit(const fit_container_type &fcon, const index_type &deg_in)
          {
            std::vector<data_type> t;
            return fit(t, fcon, deg_in);
          }

          data_type fit(std::vector<data_type> &t, const fit_container_type &fcon, const index_type &deg_in)
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
                if (ec!=fit_container_type::NO_ERROR)
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

            // calculate the error at the point
            data_type err(0), d;
            for (i=0; i<pts.size(); ++i)
            {
              geom::point::distance(d, pts[i], f(t[i]));
              err+=d;
            }

            return err;
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
              if (ec!=fit_container_type::NO_ERROR)
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
          }

        private:
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_point_matrix_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> row_pts_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> col_type;
          typedef Eigen::Matrix<data_type, 1, Eigen::Dynamic> row_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> mat_type;

        private:
          control_point_matrix_type B;      /** control polygon coordinates */

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

              geom::point::distance(temp, pts[i-1], pts[i]);

              // in case two adjacent points are given too close together
              if (temp<small_dist)
                temp=small_dist;

              t[i]=t[i-1]+temp;
            }
            len=t[npts-1];
            if (closed)
            {
              data_type d;
              geom::point::distance(d, pts[npts-1], pts[0]);
              len+=d;
            }

            for (i=0; i<npts; ++i)
              t[i]/=len;
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
