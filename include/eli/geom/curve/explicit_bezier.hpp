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

#ifndef eli_geom_curve_explicit_bezier_hpp
#define eli_geom_curve_explicit_bezier_hpp

#include <Eigen/Eigen>

#include "eli/util/tolerance.hpp"

#include "eli/geom/point/distance.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/general/continuity.hpp"
#include "eli/geom/curve/fit_container.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, typename tol__=eli::util::tolerance<data__> >
      class explicit_bezier
      {
        private:
          typedef bezier<data__, 1, tol__> curve_type;

        public:
          typedef typename curve_type::dimension_type dimension_type;
          typedef typename curve_type::data_type data_type;
          typedef Eigen::Matrix<data_type, 1, 2> point_type;
          typedef typename curve_type::point_type control_point_type;
          typedef typename curve_type::index_type index_type;
          typedef geom::curve::fit_container<data_type, index_type, 2, 1> fit_container_type;
          typedef typename curve_type::tolerance_type tolerance_type;

        public:
          explicit_bezier() {}
          explicit_bezier(const index_type &n) : y_curve(n) {}
          explicit_bezier(const explicit_bezier<data_type, tolerance_type> &eb) : y_curve(eb.y_curve) {}
          ~explicit_bezier() {}

          bool operator==(const explicit_bezier<data_type, tolerance_type> &eb) const
          {
            if (this==&eb)
              return true;
            return (y_curve==eb.y_curve);
          }

          bool operator!=(const explicit_bezier<data_type, tolerance_type> &eb) const
          {
            return !operator==(eb);
          }

          static dimension_type dimension() {return 2;}

          void resize(const index_type &t_dim)
          {
            y_curve.resize(t_dim);
          }

          bool open() const
          {
            return true;
          }
          bool closed() const
          {
            return false;
          }

          index_type degree() const
          {
            return y_curve.degree();
          }

          void set_control_point(const control_point_type &cp_in, const index_type &i)
          {
            y_curve.set_control_point(cp_in, i);
          }

          control_point_type get_control_point(const index_type &i) const
          {
            return y_curve.get_control_point(i);
          }

          void reverse()
          {
            y_curve.reverse();
          }

          point_type f(const data_type &t) const
          {
            point_type rtn;
            rtn << t, y_curve.f(t);
            return rtn;
          }

          point_type fp(const data_type &t) const
          {
            point_type rtn;
            rtn << 1, y_curve.fp(t);
            return rtn;
          }

          point_type fpp(const data_type &t) const
          {
            point_type rtn;
            rtn << 0, y_curve.fpp(t);
            return rtn;
          }

          point_type fppp(const data_type &t) const
          {
            point_type rtn;
            rtn << 0, y_curve.fppp(t);
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
            n(0)=-t(1);
            n(1)=t(0);
            b.setZero();
          }

          void degree_promote()
          {
            y_curve.degree_promote();
          }

          bool degree_demote(const geom::general::continuity &continuity_degree=geom::general::C0)
          {
            return y_curve.degree_demote(continuity_degree);
          }

          data_type fit(const fit_container_type &fcon, const index_type &deg_in)
          {
            std::vector<data_type> t;
            return fit(t, fcon, deg_in);
          }

          data_type fit(std::vector<data_type> &t, const fit_container_type &fcon, const index_type &deg_in)
          {
            size_t i, npts(fcon.number_points()), n;
            std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(npts);
            std::vector<Eigen::Matrix<data_type, 1, 1>, Eigen::aligned_allocator<Eigen::Matrix<data_type, 1, 1> > > ypts(npts);

            // get the points from the container
            fcon.get_points(pts.begin());

            // cannot have closed explicit curve
            if (fcon.closed())
            {
              assert(false);
              return -1;
            }

            // the t values are the x-coordinates of fit points
            t.resize(npts);
            for (i=0; i<npts; ++i)
            {
              t[i]=pts[i](0);
              ypts[i](0)=pts[i](1);
            }

            // determine the actual degree of curve
            n=internal::determine_n(deg_in, fcon.number_constraints(), npts);

            // build the fit terms from points
            mat_type A, x;
            row_pts_type b;
            internal::build_fit_Ab(A, b, t, ypts, n, 1);

            // handle special case of unconstrained optimization problem
            if (fcon.number_constraints()==0)
            {
              // determine the coefficients
              eli::mutil::opt::least_squares_uncon(x, A, b);
            }
            else
            {
              // now becomes a constrained least squares problem
              // the constraints come from closed flag and/or constraint collection
              size_t bi, ncon=fcon.number_constraints();

              mat_type B(ncon, n+1);
              row_pts_type d(ncon, 1);

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
                  d.row(bi)=ypts[indexes[i]];
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

              // determine the coefficients
              eli::mutil::opt::least_squares_eqcon(x, A, b, B, d);
            }

            // extract the control points and set them
            resize(n);
            for (i=0; i<=n; ++i)
            {
              set_control_point(x.row(i), i);
            }

            // calculate the error at the point
            data_type err(0);
            for (i=0; i<pts.size(); ++i)
            {
              err+=eli::geom::point::distance(pts[i], f(t[i]));
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
            size_t i, npts(fcon.number_points()), n, ai;
            std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(npts);
            std::vector<Eigen::Matrix<data_type, 1, 1>, Eigen::aligned_allocator<Eigen::Matrix<data_type, 1, 1> > > ypts(npts);

            // get the points from the container
            fcon.get_points(pts.begin());

            // cannot have closed explicit curve
            if (fcon.closed())
            {
              assert(false);
              return;
            }

            // the t values are the x-coordinates of fit points
            t.resize(npts);
            for (i=0; i<npts; ++i)
            {
              t[i]=pts[i](0);
              ypts[i](0)=pts[i](1);
            }

            // determine the actual degree of curve
            n=fcon.number_constraints(false)+npts-1;

            // build the fit terms from points
            mat_type A, x;
            row_pts_type b;
            internal::build_fit_Ab(A, b, t, ypts, n, 1);

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

            // solve for the control points
            x=A.lu().solve(b);

            // extract the control points and set them
            resize(n);
            for (i=0; i<=n; ++i)
            {
              set_control_point(x.row(i), i);
            }
          }

        private:
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> row_pts_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> col_type;
          typedef Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> mat_type;

        private:
          curve_type y_curve;
      };
    }
  }
}
#endif
