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

#include <Eigen/Eigen>

#include "eli/opt/least_squares.hpp"
#include "eli/dm/binomial_coefficient.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/general/continuity.hpp"
#include "eli/geom/curve/fit_container.hpp"
#include "eli/geom/tolerance/simple.hpp"

// TODO: Integrate the tol__ class into this class to replace open_flag and any other place
//       where numerical error might affect an equivalence comparison
namespace eli
{
  namespace geom
  {
    namespace curve
    {
      // TODO: See if these could/should be put into an unnamed namespace
      namespace internal
      {
        template<typename Derived>
        void build_N(Eigen::MatrixBase<Derived> &N, const typename Derived::Index &n)
        {
          typename Derived::Index i,j;
          typename Derived::Scalar bc1, bc2;

          // size N
          N.derived().resize(n,n);
          N.fill(0);
          for (i=0; i<n; ++i)
          {
            for (j=0; j<n; ++j)
            {
              if (i+j<n)
              {
                eli::dm::n_choose_k(bc1, n-1, j);
                eli::dm::n_choose_k(bc2, n-j-1, n-i-j-1);
                N(i,j)=bc1*bc2;
                if ((n-i-j-1)%2==1)
                  N(i,j)*=-1;
              }
            }
          }
        }

        template<typename Derived>
        void build_T(Eigen::MatrixBase<Derived> &T, const typename Derived::Scalar &t, const typename Derived::Index &n)
        {
          // check to make sure have valid curve
          assert(n>0);

          typename Derived::Index j, jj;

          T.derived().resize(n);
          T.fill(1);
          for (j=0; j<n-1; ++j)
            for (jj=0; jj<j+1; ++jj)
              T(jj)*=t;
        }

        template<typename Derived>
        void build_Tp(Eigen::MatrixBase<Derived> &Tp, const typename Derived::Scalar &t, const typename Derived::Index &n)
        {
          // check to make sure have valid curve
          assert(n>0);

          typename Derived::Index j, jj;

          Tp.derived().resize(n);
          Tp.fill(0);
          for (j=0; j<n-1; ++j)
            Tp(j)=static_cast<typename Derived::Scalar>(n-j-1);

          for (j=0; j<n-2; ++j)
            for (jj=0; jj<j+1; ++jj)
              Tp(jj)*=t;
        }

        template<typename Derived>
        void build_Tpp(Eigen::MatrixBase<Derived> &Tpp, const typename Derived::Scalar &t, const typename Derived::Index &n)
        {
          // check to make sure have valid curve
          assert(n>0);

          typename Derived::Index j, jj;

          Tpp.derived().resize(n);
          Tpp.fill(0);
          for (j=0; j<n-2; ++j)
            Tpp(j)=static_cast<typename Derived::Scalar>(n-j-1)*(n-j-2);

          for (j=0; j<n-3; ++j)
            for (jj=0; jj<j+1; ++jj)
              Tpp(jj)*=t;
        }

        template<typename Derived>
        void build_Tppp(Eigen::MatrixBase<Derived> &Tppp, const typename Derived::Scalar &t, const typename Derived::Index &n)
        {
          // check to make sure have valid curve
          assert(n>0);

          typename Derived::Index j, jj;

          Tppp.derived().resize(n);
          Tppp.fill(0);
          for (j=0; j<n-3; ++j)
            Tppp(j)=static_cast<typename Derived::Scalar>(n-j-1)*(n-j-2)*(n-j-3);

          for (j=0; j<n-4; ++j)
            for (jj=0; jj<j+1; ++jj)
              Tppp(jj)*=t;
        }

        template<typename Derived1, typename Derived2>
        void build_De_Casteljau(Eigen::MatrixBase<Derived1> &beta, const typename Derived2::Scalar &t0, const Eigen::MatrixBase<Derived2> &B, const typename Derived2::Index &sz)
        {
          typename Derived2::Index i,j,n(sz-1);

          // resize the input matrix
          beta.derived().resize(sz,sz);

          // set the initial column
          j=0;
          for (i=0; i<=n; ++i)
            beta(i,j)=B.row(i);

          // fill in the rest of terms
          for (j=1; j<=n; ++j)
          {
            for (i=0; i<=n-j; ++i)
              beta(i,j)=beta(i,j-1)*(1-t0)+beta(i+1, j-1)*t0;
          }
        }

        template<typename Derived1, typename Derived2, typename PointType>
        void build_fit_Ab(Eigen::MatrixBase<Derived1> &A, Eigen::MatrixBase<Derived2> &b, std::vector<typename Derived1::Scalar> &t, const std::vector<PointType> &pts, const typename Derived1::Index &n, const size_t &dim)
        {
          typedef Eigen::Matrix<typename Derived1::Scalar, Eigen::Dynamic, Eigen::Dynamic> mat_type;
          typedef Eigen::Matrix<typename Derived1::Scalar, Eigen::Dynamic, 1> col_type;

          typename Derived1::Index i, sz(n+1), npts(pts.size()), nrows(std::max(npts, sz));
          mat_type N, T(nrows,sz);

          // resize the return matrices
          A.derived().resize(nrows, sz);
          b.derived().resize(nrows, dim);

          // build the least squares terms
          internal::build_N(N, sz);
          for (i=0; i<npts; ++i)
          {
            col_type Tvec;
            internal::build_T(Tvec, t[i], sz);
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

      template<typename data__, unsigned short dim__, typename tol__=geom::tolerance::simple<data__> >
      class bezier
      {
        public:
          typedef unsigned short dimension_type;
          typedef data__ data_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_point_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef typename control_point_type::Index index_type;
          typedef geom::curve::fit_container<data_type, index_type, dim__, dim__> fit_container_type;
          typedef tol__ tolerance_type;

        private:
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> row_pts_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> col_type;
          typedef Eigen::Matrix<data_type, 1, Eigen::Dynamic> row_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> mat_type;
          typedef Eigen::Matrix<point_type, Eigen::Dynamic, Eigen::Dynamic> mat_3d_type;

        private:
          bool open_flag;            /** flag representing whether curve is open or closed */
          control_point_type B;      /** control polygon coordinates */
          control_point_type Bstar;  /** convenience member to precompute terms that will not
                                       * change once control polygon is set */

        private:
          void build_N(mat_type &N) const
          {
            internal::build_N(N, this->size());
          }

          void build_T(col_type &T, const data_type &t) const
          {
            internal::build_T(T, t, this->size());
          }

          void build_Tp(col_type &Tp, const data_type &t) const
          {
            internal::build_Tp(Tp, t, this->size());
          }

          void build_Tpp(col_type &Tpp, const data_type &t) const
          {
            internal::build_Tpp(Tpp, t, this->size());
          }

          void build_Tppp(col_type &Tppp, const data_type &t) const
          {
            internal::build_Tppp(Tppp, t, this->size());
          }

          void build_De_Casteljau(mat_3d_type &beta, const data_type &t0) const
          {
            internal::build_De_Casteljau(beta, t0, B, this->size());
          }

          index_type size() const
          {
            return B.rows();
          }

          void determine_t(std::vector<data_type> &t, const std::vector<point_type> &pts, bool closed) const
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

        public:
          bezier() : open_flag(true) {}
          bezier(const control_point_type &bin) : open_flag(true)
          {
            set_control_points(bin);
          }
          bezier(const bezier<data_type, dim__> &bc) : open_flag(bc.open_flag), B(bc.B), Bstar(bc.Bstar) {}
          ~bezier() {}

          bezier & operator=(const bezier<data_type, dim__> &bc)
          {
            if (this != &bc)
            {
              open_flag=bc.open_flag;
              B=bc.B;
              Bstar=bc.Bstar;
            }
            return *this;
          }

          bool operator==(const bezier<data_type, dim__> &bc) const
          {
            if (this == &bc)
              return true;
            if (open_flag!=bc.open_flag)
              return false;
            if ((B.rows()!=B.rows()) || (B.cols()!=B.cols()))
              return false;
            if (B!=bc.B)
              return false;
            if ((Bstar.rows()!=Bstar.rows()) || (Bstar.cols()!=Bstar.cols()))
              return false;
            if (Bstar!=bc.Bstar)
              return false;
            return true;
          }

          bool operator!=(const bezier<data_type, dim__> &bc) const
          {
            return !operator==(bc);
          }

          static dimension_type dimension() {return dim__;}

          bool open() const {return open_flag;}
          bool closed() const {return !open_flag;}

          void degree_promote()
          {
            // create vector of new control points
            index_type i, n(this->size());
            control_point_type B_new(n+1, dim__);

            // build the new control points
            B_new.row(0)=B.row(0);
            B_new.row(n)=B.row(n-1);
            for (i=1; i<n; ++i)
              B_new.row(i)=(B.row(i-1)-B.row(i))*(static_cast<data_type>(i)/n)+B.row(i);

            // set the new control points
            set_control_points(B_new);
          }

          bool degree_demote(const geom::general::continuity &continuity_degree=geom::general::C0)
          {
            // before doing anything make sure can demote
            if (this->degree()==0)
              return false;

            // NOTE: Implements Eck's method: Matthias Eck. Least Squares Degree reduction of Bézier curves.
            //                                Computer-Aided Design. Volume 27, No. 11. (1995), 845-851
            // NOTE: This paper also provides an algorithm for creating a bezier spline order reduced
            //       curve given a bezier spline and a desired tolerance. It also provides a means of
            //       reducing more than one order at once.
            index_type n(this->size()-1), i;
            control_point_type B_I(n, dim__), B_II(n, dim__), B_new(n, dim__);
            col_type lambda(n);
            bool rtn;

            // calculate B_I and B_II (note: n is the index of the last element in b)
            B_I.row(0)=B.row(0);
            for (i=1; i<n; ++i)
              B_I.row(i)=(B.row(i)*static_cast<data_type>(n)-B_I.row(i-1)*static_cast<data_type>(i))/static_cast<data_type>(n-i);
            B_II.row(n-1)=B.row(n);
            for (i=n-1; i>=1; --i)
              B_II.row(i-1)=(B.row(i)*static_cast<data_type>(n)-B_II.row(i)*static_cast<data_type>(n-i))/static_cast<data_type>(i);

            // calculate the blending terms
            if (continuity_degree==geom::general::NOT_CONNECTED)
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

              rtn=true;
            }
            else
            {
              int ncd(static_cast<index_type>(continuity_degree));

              // store whether able to maintain desired continuity
              if (2*(ncd+1)<=n)
              {
                rtn=true;
              }
              else
              {
                rtn=false;
                ncd=n/2-1;
              }

              // interpolating end points up-to alpha order derivative
              data_type tmp, coef, sum(0);
              index_type alpha(ncd+1);
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
              B_new.row(i)=B_I.row(i)*(1-lambda(i))+B_II.row(i)*lambda(i);

            // set the new control points
            set_control_points(B_new);

            return rtn;
          }

          data_type fit(const fit_container_type &fcon, const index_type &deg_in)
          {
            std::vector<data_type> t;
            return fit(t, fcon, deg_in);
          }

          data_type fit(std::vector<data_type> &t, const fit_container_type &fcon, const index_type &deg_in)
          {
            size_t i, npts(fcon.number_points()), n, nclosed;
            std::vector<point_type> pts(npts);

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
              eli::opt::least_squares_uncon(x, A, b);
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
                  internal::build_T(T, t[indexes[i]], n+1);
                  internal::build_N(N, n+1);
                  B.row(bi)=T.transpose()*N;
                  d.row(bi)=pts[indexes[i]];
                  ++bi;

                  // set the C1 constraint
                  if (ci.using_fp()!=fit_container_type::constraint_info::NOT_USED)
                  {
                    col_type Tp;

                    internal::build_Tp(Tp, t[indexes[i]], n+1);
                    B.row(bi)=Tp.transpose()*N;
                    d.row(bi)=ci.get_fp();
                    ++bi;
                  }

                  // set the C2 constraint
                  if (ci.using_fpp()!=fit_container_type::constraint_info::NOT_USED)
                  {
                    col_type Tpp;

                    internal::build_Tpp(Tpp, t[indexes[i]], n+1);
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
              eli::opt::least_squares_eqcon(x, A, b, B, d);
            }

            // extract the control points and set them
            control_point_type ctrl(n+1, dim__);
            for (i=0; i<n+1; ++i)
              ctrl.row(i)=x.row(i);

            // ensure that the last control point and first are the same
            if (fcon.closed())
              ctrl.row(n)=ctrl.row(0);

            set_control_points(ctrl);

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
            std::vector<point_type> pts(npts);

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
                  internal::build_N(N, n+1);

                  internal::build_Tp(Tp, t[indexes[i]], n+1);
                  A.row(ai)=Tp.transpose()*N;
                  b.row(ai)=ci.get_fp();
                  ++ai;
                }

                // set the C2 constraint
                if (ci.using_fpp()!=fit_container_type::constraint_info::NOT_USED)
                {
                  col_type Tpp;
                  mat_type N;
                  internal::build_N(N, n+1);

                  internal::build_Tpp(Tpp, t[indexes[i]], n+1);
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
            control_point_type ctrl(n+1, dim__);
            for (i=0; i<n+1; ++i)
              ctrl.row(i)=x.row(i);

            // ensure that the last control point and first are the same
            if (fcon.closed())
              ctrl.row(n)=ctrl.row(0);

            set_control_points(ctrl);
          }

          index_type degree() const
          {
            return size()-1;
          }

          void split(bezier<data_type, dim__> &bc_l, bezier<data_type, dim__> &bc_r, const data_type &t0) const
          {
            if ( (t0>1) || (t0<0) )
            {
              assert(false);
              return;
            }

            // calculate De Casteljau matrix
            mat_3d_type beta;

            build_De_Casteljau(beta, t0);

            // extract out the sub-curve control points
            control_point_type bl(this->size(), dim__), br(this->size(), dim__);
            index_type k, n(this->degree());

            for (k=0; k<=n; ++k)
            {
              bl.row(k)=beta(0,k);
              br.row(k)=beta(k, n-k);
            }

            // set the control points for the two parts
            bc_l.set_control_points(bl);
            bc_r.set_control_points(br);
          }

          void set_control_points(const control_point_type &Bin)
          {
            // set the b and size
            B=Bin;

            // compute the Bstar values
            mat_type N;

            build_N(N);
            Bstar=N*Bin;

            // check whether first and last control point are coincident
            point_type b0, bn;

            b0=B.row(0);
            bn=B.row(degree());

            open_flag=(b0!=bn);
          }

          void get_control_points(control_point_type &Bout) const
          {
            // check to make sure have valid curve
            assert(this->size()>0);

            Bout=B;
          }

          point_type f(const data_type &t) const
          {
            // check to make sure have valid curve
            assert(this->size()>0);

            // check to make sure given valid parametric value
            assert((t>=0) && (t<=1));

            point_type rtn;
            col_type T;

            build_T(T, t);
            rtn = Bstar.transpose()*T;

            return rtn;
          }

          point_type fp(const data_type &t) const
          {
            // check to make sure have valid curve
            assert(this->size()>0);

            // check to make sure given valid parametric value
            assert((t>=0) && (t<=1));

            point_type rtn;
            col_type Tp;

            build_Tp(Tp, t);
            rtn = Bstar.transpose()*Tp;
            return rtn;
          }

          point_type fpp(const data_type &t) const
          {
            // check to make sure have valid curve
            assert(this->size()>0);

            // check to make sure given valid parametric value
            assert((t>=0) && (t<=1));

            point_type rtn;
            col_type Tpp;

            build_Tpp(Tpp, t);
            rtn = Bstar.transpose()*Tpp;
            return rtn;
          }

          point_type fppp(const data_type &t) const
          {
            // check to make sure have valid curve
            assert(this->size()>0);

            // check to make sure given valid parametric value
            assert((t>=0) && (t<=1));

            point_type rtn;
            col_type Tppp;

            build_Tppp(Tppp, t);
            rtn = Bstar.transpose()*Tppp;
            return rtn;
          }

          // TODO: Implement additional functionality
          // derivative of the curve with respect to the ith control point
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
