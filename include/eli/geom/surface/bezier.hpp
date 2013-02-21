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

#ifndef eli_geom_surface_bezier_hpp
#define eli_gome_surface_bezier_hpp

#include <cstddef>  // nullptr

#include <vector>

#include "Eigen/Eigen"

#include "eli/geom/general/continuity.hpp"
#include "eli/geom/tolerance/simple.hpp"
#include "eli/dm/binomial_coefficient.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
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
            index_type ncd(ncon/2);

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
            cp_out.row(i)=B_I.row(i)*(1-lambda(i))+B_II.row(i)*lambda(i);
        }

        // get the coefficients for the split curves by constructing a de Casteljau triangle scheme.
        // The upper side is the t<t0 section, and the lower side is the t>t0 section.
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
      }
    }
  }
}

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      template<typename data__, unsigned short dim__, typename tol__=geom::tolerance::simple<data__> >
      class bezier
      {
        public:
          typedef unsigned short dimension_type;
          typedef data__ data_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef Eigen::Map<Eigen::Matrix<data_type, Eigen::Dynamic, dim__>,
                             Eigen::Unaligned,
                             Eigen::Stride<1, dim__> > control_point_type;
          typedef typename control_point_type::Index index_type;
          typedef Eigen::Map<Eigen::Matrix<data_type, Eigen::Dynamic, dim__>,
                             Eigen::Unaligned,
                             Eigen::Stride<1, Eigen::Dynamic> > v_dir_control_point_type;
//           typedef geom::curve::fit_container<data_type, index_type, dim__, dim__> fit_container_type;
          typedef tol__ tolerance_type;

        private:
//           typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> row_pts_type;
//           typedef Eigen::Matrix<data_type, Eigen::Dynamic, 1> col_type;
//           typedef Eigen::Matrix<data_type, 1, Eigen::Dynamic> row_type;
//           typedef Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> mat_type;
//           typedef Eigen::Matrix<point_type, Eigen::Dynamic, Eigen::Dynamic> mat_3d_type;

        private:
          // TODO: Should the curve class have similar setup?
          // TODO: Should the curve class be changed to use deCasteljau algorithm for evaluations?
          std::vector<data_type> point_data; /** raw control points stored in vector. The order is {x,y...}_(0,0) {x,y...}_(1,0) ... {x,y...}_(n,0) {x,y...}_(0,1) {x,y...}_(1,1) ... {x,y...}_(n,1) ... {x,y...}_(n,m) */
          std::vector<control_point_type> B_u; /** vector of u-direction, i.e. direction where v is constant, curve control points in point_data */
          std::vector<v_dir_control_point_type> B_v; /** vector of v-direction, i.e. direction where u is constant, curve control points in point_data */

        public:
          bezier() : point_data(3, 0)
          {
            // set the B_u and B_v maps
            set_Bs(1, 1);
          }
          bezier(const index_type &u_dim, const index_type &v_dim) : point_data(dim__*(u_dim+1)*(v_dim+1))
          {
            // set the B_u and B_v maps
            set_Bs(u_dim, v_dim);
          }
          bezier(const bezier<data_type, dim__, tol__> &bs) : point_data(bs.point_data)
          {
            // set the B_u and B_v maps
            set_Bs(bs.degree_u(), bs.degree_v());
          }
          ~bezier() {}

          bezier & operator=(const bezier<data_type, dim__, tol__> &bs)
          {
            if (this!=&bs)
            {
              point_data=bs.point_data;
              set_Bs(bs.degree_u(), bs.degree_v());
            }

            return (*this);
          }

          bool operator==(const bezier<data_type, dim__, tol__> &bs) const
          {
            if (this==&bs)
              return true;

            if (point_data!=bs.point_data)
              return false;

            if (B_u.size()!=bs.B_u.size())
              return false;

            if (B_v.size()!=bs.B_v.size())
              return false;

            return true;
          }

          bool operator!=(const bezier<data_type, dim__, tol__> &bs) const
          {
            return !operator==(bs);
          }

          static dimension_type dimension() {return dim__;}

          index_type degree_u() const {return static_cast<index_type>(B_v.size())-1;}
          index_type degree_v() const {return static_cast<index_type>(B_u.size())-1;}

          void resize(const index_type &u_dim, const index_type &v_dim)
          {
            // allocate the control points
            point_data.resize(dim__*(u_dim+1)*(v_dim+1));

            // set the B_u and B_v maps
            set_Bs(u_dim, v_dim);
          }

          point_type get_control_point(const index_type &i, const index_type &j) const
          {
            // make sure have valid indexes
            if ( (i>degree_u()) || (j>degree_v()) )
            {
              assert(false);
              return B_u[0].row(0);
            }

            // make sure that the two maps point to same items
            assert(B_u[j].row(i)==B_v[i].row(j));

            return B_u[j].row(i);
          }

          void set_control_point(const point_type &cp, const index_type &i, const index_type &j)
          {
            // make sure have valid indexes
            if ( (i>degree_u()) || (j>degree_v()) )
            {
              assert(false);
              return;
            }

            B_u[j].row(i) = cp;
          }

// NOTE: These would need to be done for u-direction and v-direction
//           bool open() const {return open_flag;}
//           bool closed() const {return !open_flag;}

          point_type f(const data_type &u, const data_type &v) const
          {
            point_type ans, tmp;
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp;
            index_type i, n(degree_u()), m(degree_v());

            // check to make sure have valid curve
            assert(this->degree_u()>=0);
            assert(this->degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (n<=m)
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                eli::geom::curve::utility::de_casteljau(tmp, B_u[i], u);
                temp_cp.row(i)=tmp;
              }
              eli::geom::curve::utility::de_casteljau(ans, temp_cp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                eli::geom::curve::utility::de_casteljau(tmp, B_v[i], v);
                temp_cp.row(i)=tmp;
              }
              eli::geom::curve::utility::de_casteljau(ans, temp_cp, u);
            }

            return ans;
          }

          point_type f_u(const data_type &u, const data_type &v) const
          {
            point_type ans, tmp;
            index_type i, n(degree_u()), m(degree_v());
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp, B_up(n+1-1, dim__);

            // check to make sure have valid curve
            assert(this->degree_u()>=0);
            assert(this->degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_u()<1)
            {
              ans.setZero();
              return ans;
            }

            if ((n-1)<=m)
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                // create the control points for first derivative curve
                eli::geom::curve::utility::bezier_p_control_point(B_up, B_u[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_up, u);
                temp_cp.row(i)=tmp;
              }
              eli::geom::curve::utility::de_casteljau(ans, temp_cp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                eli::geom::curve::utility::de_casteljau(tmp, B_v[i], v);
                temp_cp.row(i)=tmp;
              }
              // create the control points for first derivative curve
              eli::geom::curve::utility::bezier_p_control_point(B_up, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_up, u);
            }

            return ans;
          }

          point_type f_v(const data_type &u, const data_type &v) const
          {
            point_type ans, tmp;
            index_type i, n(degree_u()), m(degree_v());
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp, B_vp(m+1-1, dim__);

            // check to make sure have valid curve
            assert(this->degree_u()>=0);
            assert(this->degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_v()<1)
            {
              ans.setZero();
              return ans;
            }

            if (n<=(m-1))
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                eli::geom::curve::utility::de_casteljau(tmp, B_u[i], u);
                temp_cp.row(i)=tmp;
              }
              // create the control points for first derivative curve
              eli::geom::curve::utility::bezier_p_control_point(B_vp, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_vp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                // create the control points for first derivative curve
                eli::geom::curve::utility::bezier_p_control_point(B_vp, B_v[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_vp, v);
                temp_cp.row(i)=tmp;
              }
              eli::geom::curve::utility::de_casteljau(ans, temp_cp, u);
            }

            return ans;
          }

          point_type f_uu(const data_type &u, const data_type &v) const
          {
            point_type ans, tmp;
            index_type i, n(degree_u()), m(degree_v());
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp, B_upp(n+1-2, dim__);

            // check to make sure have valid curve
            assert(this->degree_u()>=0);
            assert(this->degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_u()<2)
            {
              ans.setZero();
              return ans;
            }

            if ((n-2)<=m)
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                // create the control points for second derivative curve
                eli::geom::curve::utility::bezier_pp_control_point(B_upp, B_u[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_upp, u);
                temp_cp.row(i)=tmp;
              }
              eli::geom::curve::utility::de_casteljau(ans, temp_cp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                eli::geom::curve::utility::de_casteljau(tmp, B_v[i], v);
                temp_cp.row(i)=tmp;
              }

              // create the control points for second derivative curve
              eli::geom::curve::utility::bezier_pp_control_point(B_upp, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_upp, u);
            }

            return ans;
          }

          point_type f_uv(const data_type &u, const data_type &v) const
          {
            point_type ans, tmp;
            index_type i, n(degree_u()), m(degree_v());
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp, B_up(n+1-1, dim__), B_vp(m+1-1, dim__);

            // check to make sure have valid curve
            assert(this->degree_u()>=0);
            assert(this->degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if ( (this->degree_u()<1) || (this->degree_v()<1) )
            {
              ans.setZero();
              return ans;
            }

            if (n<=m)
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                // create the control points for first u-derivative curve
                eli::geom::curve::utility::bezier_p_control_point(B_up, B_u[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_up, u);
                temp_cp.row(i)=tmp;
              }
              // create the control points for first v-derivative curve
              eli::geom::curve::utility::bezier_p_control_point(B_vp, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_vp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                // create the control points for first v-derivative curve
                eli::geom::curve::utility::bezier_p_control_point(B_vp, B_v[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_vp, v);
                temp_cp.row(i)=tmp;
              }
              // create the control points for first u-derivative curve
              eli::geom::curve::utility::bezier_p_control_point(B_up, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_up, u);
            }

            return ans;
          }

          point_type f_vv(const data_type &u, const data_type &v) const
          {
            point_type ans, tmp;
            index_type i, n(degree_u()), m(degree_v());
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp, B_vpp(m+1-2, dim__);

            // check to make sure have valid curve
            assert(this->degree_u()>=0);
            assert(this->degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_v()<2)
            {
              ans.setZero();
              return ans;
            }

            if (n<=(m-2))
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                eli::geom::curve::utility::de_casteljau(tmp, B_u[i], u);
                temp_cp.row(i)=tmp;
              }
              // create the control points for second derivative curve
              eli::geom::curve::utility::bezier_pp_control_point(B_vpp, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_vpp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                // create the control points for second derivative curve
                eli::geom::curve::utility::bezier_pp_control_point(B_vpp, B_v[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_vpp, v);
                temp_cp.row(i)=tmp;
              }
              eli::geom::curve::utility::de_casteljau(ans, temp_cp, u);
            }

            return ans;
          }

          point_type f_uuu(const data_type &u, const data_type &v) const
          {

            point_type ans, tmp;
            index_type i, n(degree_u()), m(degree_v());
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp, B_uppp(n+1-3, dim__);

            // check to make sure have valid curve
            assert(this->degree_u()>=0);
            assert(this->degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_u()<3)
            {
              ans.setZero();
              return ans;
            }

            if ((n-3)<=m)
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                // create the control points for third derivative curve
                eli::geom::curve::utility::bezier_ppp_control_point(B_uppp, B_u[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_uppp, u);
                temp_cp.row(i)=tmp;
              }
              eli::geom::curve::utility::de_casteljau(ans, temp_cp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                eli::geom::curve::utility::de_casteljau(tmp, B_v[i], v);
                temp_cp.row(i)=tmp;
              }

              // create the control points for third derivative curve
              eli::geom::curve::utility::bezier_pp_control_point(B_uppp, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_uppp, u);
            }

            return ans;
        }

          point_type f_uuv(const data_type &u, const data_type &v) const
          {
            point_type ans, tmp;
            index_type i, n(degree_u()), m(degree_v());
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp, B_upp(n+1-2, dim__), B_vp(m+1-1, dim__);

            // check to make sure have valid curve
            assert(this->degree_u()>=0);
            assert(this->degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if ( (this->degree_u()<2) || (this->degree_v()<1) )
            {
              ans.setZero();
              return ans;
            }

            if ((n-1)<=m)
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                // create the control points for second u-derivative curve
                eli::geom::curve::utility::bezier_pp_control_point(B_upp, B_u[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_upp, u);
                temp_cp.row(i)=tmp;
              }
              // create the control points for first v-derivative curve
              eli::geom::curve::utility::bezier_p_control_point(B_vp, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_vp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                // create the control points for first v-derivative curve
                eli::geom::curve::utility::bezier_p_control_point(B_vp, B_v[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_vp, v);
                temp_cp.row(i)=tmp;
              }
              // create the control points for second u-derivative curve
              eli::geom::curve::utility::bezier_pp_control_point(B_upp, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_upp, u);
            }

            return ans;
          }

          point_type f_uvv(const data_type &u, const data_type &v) const
          {
            point_type ans, tmp;
            index_type i, n(degree_u()), m(degree_v());
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp, B_up(n+1-1, dim__), B_vpp(m+1-2, dim__);

            // check to make sure have valid curve
            assert(this->degree_u()>=0);
            assert(this->degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if ( (this->degree_u()<1) || (this->degree_v()<2) )
            {
              ans.setZero();
              return ans;
            }

            if (n<=(m-1))
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                // create the control points for first u-derivative curve
                eli::geom::curve::utility::bezier_p_control_point(B_up, B_u[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_up, u);
                temp_cp.row(i)=tmp;
              }
              // create the control points for second v-derivative curve
              eli::geom::curve::utility::bezier_pp_control_point(B_vpp, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_vpp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                // create the control points for second v-derivative curve
                eli::geom::curve::utility::bezier_pp_control_point(B_vpp, B_v[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_vpp, v);
                temp_cp.row(i)=tmp;
              }
              // create the control points for first u-derivative curve
              eli::geom::curve::utility::bezier_p_control_point(B_up, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_up, u);
            }

            return ans;
          }

          point_type f_vvv(const data_type &u, const data_type &v) const
          {
            point_type ans, tmp;
            index_type i, n(degree_u()), m(degree_v());
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp, B_vppp(m+1-3, dim__);

            // check to make sure have valid curve
            assert(this->degree_u()>=0);
            assert(this->degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_v()<3)
            {
              ans.setZero();
              return ans;
            }

            if (n<=(m-2))
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                eli::geom::curve::utility::de_casteljau(tmp, B_u[i], u);
                temp_cp.row(i)=tmp;
              }
              // create the control points for third derivative curve
              eli::geom::curve::utility::bezier_ppp_control_point(B_vppp, temp_cp);
              eli::geom::curve::utility::de_casteljau(ans, B_vppp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                // create the control points for third derivative curve
                eli::geom::curve::utility::bezier_ppp_control_point(B_vppp, B_v[i]);
                eli::geom::curve::utility::de_casteljau(tmp, B_vppp, v);
                temp_cp.row(i)=tmp;
              }
              eli::geom::curve::utility::de_casteljau(ans, temp_cp, u);
            }

            return ans;
          }

          void promote_u()
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_row_type;
            typedef std::vector<control_row_type, Eigen::aligned_allocator<control_row_type> > control_row_collection_type;

            index_type i, n(degree_u()), m(degree_v());
            control_row_collection_type current_row(m+1, control_row_type(n+1, dim__));

            // copy the control rows
            for (i=0; i<=m; ++i)
              current_row[i]=B_u[i];

            // resize current surface
            resize(n+1, m);

            // set the new control points
            control_row_type tmp_cp(n+2, dim__);
            for (i=0; i<=m; ++i)
            {
              eli::geom::curve::utility::bezier_promote_control_points(tmp_cp, current_row[i]);
              B_u[i]=tmp_cp;
            }
          }

          void promote_v()
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_col_type;
            typedef std::vector<control_col_type, Eigen::aligned_allocator<control_col_type> > control_col_collection_type;

            index_type i, n(degree_u()), m(degree_v());
            control_col_collection_type current_col(n+1, control_col_type(m+1, dim__));

            // copy the control cols
            for (i=0; i<=n; ++i)
              current_col[i]=B_v[i];

            // resize current surface
            resize(n, m+1);

            // set the new control points
            control_col_type tmp_cp(m+2, dim__);
            for (i=0; i<=n; ++i)
            {
              eli::geom::curve::utility::bezier_promote_control_points(tmp_cp, current_col[i]);
              B_v[i]=tmp_cp;
            }
          }

          bool demote_u(const geom::general::continuity &u_continuity_degree=geom::general::C0)
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_row_type;
            typedef std::vector<control_row_type, Eigen::aligned_allocator<control_row_type> > control_row_collection_type;

            index_type i, n(degree_u()), m(degree_v());
            control_row_collection_type current_row(m+1, control_row_type(n+1, dim__));

            // check if can demote
            int ncon(0);
            switch(u_continuity_degree)
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
            if (ncon>n-2)
              return false;

            // copy the control rows
            for (i=0; i<=m; ++i)
              current_row[i]=B_u[i];

            // resize current surface
            resize(n-1, m);

            // set the new control points
            control_row_type tmp_cp(n, dim__);
            for (i=0; i<=m; ++i)
            {
              eli::geom::curve::utility::bezier_demote_control_points(tmp_cp, current_row[i], ncon);
              B_u[i]=tmp_cp;
            }

            return true;
          }

          bool demote_v(const geom::general::continuity &v_continuity_degree=geom::general::C0)
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_col_type;
            typedef std::vector<control_col_type, Eigen::aligned_allocator<control_col_type> > control_col_collection_type;

            index_type i, n(degree_u()), m(degree_v());
            control_col_collection_type current_col(n+1, control_col_type(m+1, dim__));

            // check if can demote
            int ncon(0);
            switch(v_continuity_degree)
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
            if (ncon>m-2)
              return false;

            // copy the control rows
            for (i=0; i<=n; ++i)
              current_col[i]=B_v[i];

            // resize current surface
            resize(n, m-1);

            // set the new control points
            control_col_type tmp_cp(m, dim__);
            for (i=0; i<=n; ++i)
            {
              eli::geom::curve::utility::bezier_demote_control_points(tmp_cp, current_col[i], ncon);
              B_v[i]=tmp_cp;
            }

            return true;
          }

          void split_u(bezier<data_type, dim__, tol__> &bs_lo, bezier<data_type, dim__, tol__> &bs_hi, const data_type &u0) const
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_row_type;

            index_type i, j, n(degree_u()), m(degree_v());
            control_row_type cp_lo(n+1, dim__), cp_hi(n+1, dim__);

            // make sure have valid index
            assert((u0>=0) && (u0<=1));

            // resize the surfaces
            bs_lo.resize(n, m);
            bs_hi.resize(n, m);

            // cycle through each row and split each it
            for (j=0; j<=m; ++j)
            {
              eli::geom::curve::utility::bezier_split_control_points(cp_lo, cp_hi, B_u[j], u0);
              for (i=0; i<=n; ++i)
              {
                bs_lo.set_control_point(cp_lo.row(i), i, j);
                bs_hi.set_control_point(cp_hi.row(i), i, j);
              }
            }
          }

          void split_v(bezier<data_type, dim__, tol__> &bs_lo, bezier<data_type, dim__, tol__> &bs_hi, const data_type &v0) const
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_col_type;

            index_type i, j, n(degree_u()), m(degree_v());
            control_col_type cp_lo(m+1, dim__), cp_hi(m+1, dim__);

            // make sure have valid index
            assert((v0>=0) && (v0<=1));

            // resize the surfaces
            bs_lo.resize(n, m);
            bs_hi.resize(n, m);

            // cycle through each col and split each it
            for (i=0; i<=n; ++i)
            {
              eli::geom::curve::utility::bezier_split_control_points(cp_lo, cp_hi, B_v[i], v0);
              for (j=0; j<=m; ++j)
              {
                bs_lo.set_control_point(cp_lo.row(j), i, j);
                bs_hi.set_control_point(cp_hi.row(j), i, j);
              }
            }
          }

        private:
          void set_Bs(index_type n, index_type m)
          {
            // allocate vectors of control point maps
            B_u.resize(m+1, control_point_type(nullptr, m+1, dim__, Eigen::Stride<1, dim__>()));
            for (index_type j=0; j<=m; ++j)
            {
              new (B_u.data()+j) control_point_type(point_data.data()+j*(n+1)*dim__, n+1, dim__, Eigen::Stride<1, dim__>());
            }

            B_v.resize(n+1, v_dir_control_point_type(nullptr, n+1, dim__, Eigen::Stride<1, Eigen::Dynamic>(1, (n+1)*dim__)));
            for (index_type i=0; i<=n; ++i)
            {
              new (B_v.data()+i) v_dir_control_point_type(point_data.data()+i*dim__, m+1, dim__, Eigen::Stride<1, Eigen::Dynamic>(1, (n+1)*dim__));
            }
          }
      };

      typedef bezier<float, 3> bezier3f;
      typedef bezier<double, 2> bezier2d;
      typedef bezier<double, 3> bezier3d;
      typedef bezier<long double, 2> bezier2ld;
      typedef bezier<long double, 3> bezier3ld;
    }
  }
}

#endif
