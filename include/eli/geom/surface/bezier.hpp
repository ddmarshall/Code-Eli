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
#define eli_geom_surface_bezier_hpp

#include <cstddef>  // nullptr

#include <vector>

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/geom/utility/bezier.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/equivalent_curves.hpp"
#include "eli/geom/general/continuity.hpp"
#include "eli/geom/general/bounding_box.hpp"
#include "eli/geom/point/distance.hpp"

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      template<typename data__, unsigned short dim__, typename tol__=eli::util::tolerance<data__> >
      class bezier
      {
        public:
          typedef unsigned short dimension_type;
          typedef data__ data_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef point_type control_point_type;
          typedef typename control_point_type::Index index_type;
          typedef tol__ tolerance_type;
          typedef Eigen::Matrix<data_type, dim__, dim__> rotation_matrix_type;
          typedef eli::geom::general::bounding_box<data_type, dim__, tolerance_type> bounding_box_type;
          typedef eli::geom::curve::bezier<data_type, dim__, tolerance_type> curve_type;

          typedef bezier<data_type, 1, tolerance_type> onedbezsurf;

        private:
          typedef Eigen::Map<Eigen::Matrix<data_type, Eigen::Dynamic, dim__>,
                             Eigen::Unaligned,
                             Eigen::Stride<1, dim__> > control_point_matrix_type;
          typedef Eigen::Map<Eigen::Matrix<data_type, Eigen::Dynamic, dim__>,
                             Eigen::Unaligned,
                             Eigen::Stride<1, Eigen::Dynamic> > v_dir_control_point_matrix_type;
          typedef std::vector<data_type> control_point_container;
          typedef std::vector<control_point_matrix_type> u_control_point_matrix_container;
          typedef std::vector<v_dir_control_point_matrix_type> v_control_point_matrix_container;

        private:
          control_point_container point_data; /** raw control points stored in vector. The order is {x,y...}_(0,0) {x,y...}_(1,0) ... {x,y...}_(n,0) {x,y...}_(0,1) {x,y...}_(1,1) ... {x,y...}_(n,1) ... {x,y...}_(n,m) */
          u_control_point_matrix_container B_u; /** vector of u-direction, i.e. direction where v is constant, curve control points in point_data */
          v_control_point_matrix_container B_v; /** vector of v-direction, i.e. direction where u is constant, curve control points in point_data */

          mutable bezier<data_type, dim__, tol__> * deriv_u;
          mutable bezier<data_type, dim__, tol__> * deriv_v;

        public:
          bezier() : point_data(3, 0), deriv_u( NULL ), deriv_v( NULL )
          {
            // set the B_u and B_v maps
            set_Bs(1, 1);
          }
          bezier(const index_type &u_dim, const index_type &v_dim) : point_data(dim__*(u_dim+1)*(v_dim+1)), deriv_u( NULL ), deriv_v( NULL )
          {
            // set the B_u and B_v maps
            set_Bs(u_dim, v_dim);
          }
          bezier(const bezier<data_type, dim__, tol__> &bs) : point_data(bs.point_data)
          {
            // set the B_u and B_v maps
            set_Bs(bs.degree_u(), bs.degree_v());

            if (bs.deriv_u)
            {
              deriv_u = new bezier<data_type, dim__>( *(bs.deriv_u) );
            }
            else
            {
              deriv_u = NULL;
            }

            if (bs.deriv_v)
            {
              deriv_v = new bezier<data_type, dim__>( *(bs.deriv_v) );
            }
            else
            {
              deriv_v = NULL;
            }
          }

          ~bezier()
          {
            invalidate_deriv();
          }

          bezier & operator=(const bezier<data_type, dim__, tol__> &bs)
          {
            if (this!=&bs)
            {
              point_data=bs.point_data;
              set_Bs(bs.degree_u(), bs.degree_v());

              if (bs.deriv_u)
              {
                deriv_u = new bezier<data_type, dim__>( *(bs.deriv_u) );
              }
              else
              {
                deriv_u = NULL;
              }

              if (bs.deriv_v)
              {
                deriv_v = new bezier<data_type, dim__>( *(bs.deriv_v) );
              }
              else
              {
                deriv_v = NULL;
              }

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

          bool abouteq(const bezier<data_type, dim__, tol__> &bs, const data_type &ttol2 ) const
          {
            if (this==&bs)
              return true;

            if (B_u.size()!=bs.B_u.size())
              return false;

            if (B_v.size()!=bs.B_v.size())
              return false;

            index_type i, j, degu(degree_u()), degv(degree_v());
            for (j=0; j<=degv; ++j)
            {
              for (i=0; i<=degu; ++i)
              {
                if ( eli::geom::point::distance2( B_u[j].row(i), bs.B_u[j].row(i) ) > ttol2 )
                {
                  return false;
                }
              }
            }

            return true;
          }

          bool operator!=(const bezier<data_type, dim__, tol__> &bs) const
          {
            return !operator==(bs);
          }

          static dimension_type dimension() {return dim__;}

          index_type degree_u() const {return static_cast<index_type>(B_v.size())-1;}
          index_type degree_v() const {return static_cast<index_type>(B_u.size())-1;}

          void get_parameter_min(data_type &umin, data_type &vmin) const
          {
            umin=0;
            vmin=0;
          }

          void get_parameter_max(data_type &umax, data_type &vmax) const
          {
            umax=1;
            vmax=1;
          }

          data_type get_umin() const {return static_cast<data_type>(0);}
          data_type get_vmin() const {return static_cast<data_type>(0);}
          data_type get_umax() const {return static_cast<data_type>(1);}
          data_type get_vmax() const {return static_cast<data_type>(1);}

          void resize(const index_type &u_dim, const index_type &v_dim)
          {
            resize( B_u, B_v, point_data, u_dim, v_dim );
            invalidate_deriv();
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

          void get_bounding_box(bounding_box_type &bb) const
          {
            index_type i, j, degu(degree_u()), degv(degree_v());

            bb.clear();
            for (i=0; i<=degu; ++i)
            {
              for (j=0; j<=degv; ++j)
              {
                bb.add(B_u[j].row(i));
              }
            }
          }

          void rotate(const rotation_matrix_type &rmat)
          {
            index_type j, degv(degree_v());
            for (j=0; j<=degv; ++j)
            {
              B_u[j]*=rmat.transpose();
            }
            invalidate_deriv();
          }

          void rotate(const rotation_matrix_type &rmat, const point_type &rorig)
          {
            translate(-rorig);
            rotate(rmat);
            translate(rorig);
            invalidate_deriv();
          }

          void translate(const point_type &trans)
          {
            index_type i, j, degu(degree_u()), degv(degree_v());
            for (j=0; j<=degv; ++j)
            {
              for (i=0; i<=degu; ++i)
              {
                B_u[j].row(i)+=trans;
              }
            }
            invalidate_deriv();
          }

          void scale(const data_type &s)
          {
            index_type i, j, degu(degree_u()), degv(degree_v());
            for (j=0; j<=degv; ++j)
            {
              for (i=0; i<=degu; ++i)
              {
                B_u[j].row(i)*=s;
              }
            }
            invalidate_deriv();
          }

          bool open_u() const {return !closed_u();}
          bool closed_u() const
          {
            curve_type bc0, bc1;

            get_uconst_curve(bc0, 0);
            get_uconst_curve(bc1, 1);

            return eli::geom::curve::equivalent_curves(bc0, bc1);
          }
          bool open_v() const {return !closed_v();}
          bool closed_v() const
          {
            curve_type bc0, bc1;

            get_vconst_curve(bc0, 0);
            get_vconst_curve(bc1, 1);

            return eli::geom::curve::equivalent_curves(bc0, bc1);
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
            invalidate_deriv();
          }

          void reverse_u()
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_col_type;
            typedef std::vector<control_col_type, Eigen::aligned_allocator<control_col_type> > control_col_collection_type;

            index_type i, n(degree_u()), m(degree_v());
            control_col_collection_type current_col(n+1, control_col_type(m+1, dim__));

            // copy the control cols
            for (i=0; i<=n; ++i)
              current_col[i]=B_v[i];

            // set the new control points
            for (i=0; i<=n; ++i)
            {
              B_v[i]=current_col[n-i];
            }
            invalidate_deriv();
          }

          void reverse_v()
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_row_type;
            typedef std::vector<control_row_type, Eigen::aligned_allocator<control_row_type> > control_row_collection_type;

            index_type i, n(degree_u()), m(degree_v());
            control_row_collection_type current_row(m+1, control_row_type(n+1, dim__));

            // copy the control rows
            for (i=0; i<=m; ++i)
              current_row[i]=B_u[i];

            // set the new control points
            for (i=0; i<=m; ++i)
            {
              B_u[i]=current_row[m-i];
            }
            invalidate_deriv();
          }

          void swap_uv()
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_row_type;
            typedef std::vector<control_row_type, Eigen::aligned_allocator<control_row_type> > control_row_collection_type;

            index_type i, j, n(degree_u()), m(degree_v());
            control_row_collection_type current_row(m+1, control_row_type(n+1, dim__));

            // copy the control rows
            for (i=0; i<=m; ++i)
              current_row[i]=B_u[i];

            // resize current surface
            resize(m, n);

            // set the new control points
            for (i=0; i<=m; ++i)
            {
              for (j=0; j<=n; ++j)
              {
                set_control_point(current_row[i].row(j), i, j);
              }
            }
            invalidate_deriv();
          }

          void get_uconst_curve(curve_type &bc, const data_type &u) const
          {
            index_type j, m(degree_v());

            // check to make sure have valid curve
            assert(m>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));

            // build curve
            point_type cp;

            bc.resize(m);
            for (j=0; j<=m; ++j)
            {
              eli::geom::utility::de_casteljau(cp, B_u[j], u);
              bc.set_control_point(cp, j);
            }
          }

          void get_umin_bndy_curve( curve_type &bc ) const
          {
            index_type j, m(degree_v());

            // check to make sure have valid curve
            assert(m>=0);

            // build curve
            point_type cp;

            bc.resize(m);
            for (j=0; j<=m; ++j)
            {
              bc.set_control_point(B_v[0].row(j), j);
            }
          }

          void get_umax_bndy_curve( curve_type &bc ) const
          {
            index_type j, m(degree_v()), n(degree_u());

            // check to make sure have valid curve
            assert(m>=0);

            // build curve
            point_type cp;

            bc.resize(m);
            for (j=0; j<=m; ++j)
            {
              bc.set_control_point(B_v[n].row(j), j);
            }
          }

          void get_vmin_bndy_curve( curve_type &bc ) const
          {
            index_type i, n(degree_u());

            // check to make sure have valid curve
            assert(n>=0);

            // build curve
            point_type cp;

            bc.resize(n);
            for (i=0; i<=n; ++i)
            {
              bc.set_control_point(B_u[0].row(i), i);
            }
          }

          void get_vmax_bndy_curve( curve_type &bc ) const
          {
            index_type i, n(degree_u()), m(degree_v());

            // check to make sure have valid curve
            assert(n>=0);

            // build curve
            point_type cp;

            bc.resize(n);
            for (i=0; i<=n; ++i)
            {
              bc.set_control_point(B_u[m].row(i), i);
            }
          }

          void get_umin_ndelta_pcurve( curve_type &bc ) const
          {
            index_type j, m(degree_v());

            tolerance_type tol;

            // check to make sure have valid curve
            assert(m>=0);

            // build curve
            point_type cp;

            bc.resize(m);
            for (j=0; j<=m; ++j)
            {
              point_type p = B_v[0].row(j) - B_v[1].row(j);
              data_type len = p.norm();
              if ( tol.approximately_equal(len, 0) )
              {
                p << 0, 0, 0;
              }
              else
              {
                p = p / len;
              }
              bc.set_control_point( p, j);
            }
          }

          void get_umax_ndelta_pcurve( curve_type &bc ) const
          {
            index_type j, m(degree_v()), n(degree_u());

            tolerance_type tol;

            // check to make sure have valid curve
            assert(m>=0);

            // build curve
            point_type cp;

            bc.resize(m);
            for (j=0; j<=m; ++j)
            {
              point_type p = B_v[n].row(j) - B_v[n-1].row(j);
              data_type len = p.norm();
              if ( tol.approximately_equal(len, 0) )
              {
                p << 0, 0, 0;
              }
              else
              {
                p = p / len;
              }
              bc.set_control_point( p, j);
            }
          }

          void get_vmin_ndelta_pcurve( curve_type &bc ) const
          {
            index_type i, n(degree_u());

            tolerance_type tol;

            // check to make sure have valid curve
            assert(n>=0);

            // build curve
            point_type cp;

            bc.resize(n);
            for (i=0; i<=n; ++i)
            {
              point_type p = B_u[0].row(i) - B_u[1].row(i);
              data_type len = p.norm();
              if ( tol.approximately_equal(len, 0) )
              {
                p << 0, 0, 0;
              }
              else
              {
                p = p / len;
              }
              bc.set_control_point( p, i);
            }
          }

          void get_vmax_ndelta_pcurve( curve_type &bc ) const
          {
            index_type i, n(degree_u()), m(degree_v());

            tolerance_type tol;

            // check to make sure have valid curve
            assert(n>=0);

            // build curve
            point_type cp;

            bc.resize(n);
            for (i=0; i<=n; ++i)
            {
              point_type p = B_u[m].row(i) - B_u[m-1].row(i);
              data_type len = p.norm();
              if ( tol.approximately_equal(len, 0) )
              {
                p << 0, 0, 0;
              }
              else
              {
                p = p / len;
              }
              bc.set_control_point( p, i);
            }
          }

          void get_uconst_f_u_curve(curve_type &bc, const data_type &u) const
          {
            validate_u();

            return deriv_u->get_uconst_curve( bc, u );
          }

          void get_uconst_f_v_curve(curve_type &bc, const data_type &u) const
          {
            validate_v();

            return deriv_v->get_uconst_curve( bc, u );
          }

          void get_vconst_f_u_curve(curve_type &bc, const data_type &v) const
          {
            validate_u();

            return deriv_u->get_vconst_curve( bc, v );
          }

          void get_vconst_f_v_curve(curve_type &bc, const data_type &v) const
          {
            validate_v();

            return deriv_v->get_vconst_curve( bc, v );
          }

          void get_vconst_curve(curve_type &bc, const data_type &v) const
          {
            index_type i, n(degree_u());

            // check to make sure have valid curve
            assert(n>=0);

            // check to make sure given valid parametric value
            assert((v>=0) && (v<=1));

            // build curve
            point_type cp;

            bc.resize(n);
            for (i=0; i<=n; ++i)
            {
              eli::geom::utility::de_casteljau(cp, B_v[i], v);
              bc.set_control_point(cp, i);
            }
          }

          point_type f(const data_type &u, const data_type &v) const
          {
            point_type ans, tmp;
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> temp_cp;
            index_type i, n(degree_u()), m(degree_v());

            // check to make sure have valid curve
            assert(n>=0);
            assert(m>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (n<=m)
            {
              temp_cp.resize(m+1, dim__);
              // build the temporary control points
              for (i=0; i<=m; ++i)
              {
                eli::geom::utility::de_casteljau(tmp, B_u[i], u);
                temp_cp.row(i)=tmp;
              }
              eli::geom::utility::de_casteljau(ans, temp_cp, v);
            }
            else
            {
              temp_cp.resize(n+1, dim__);
              // build the temporary control points
              for (i=0; i<=n; ++i)
              {
                eli::geom::utility::de_casteljau(tmp, B_v[i], v);
                temp_cp.row(i)=tmp;
              }
              eli::geom::utility::de_casteljau(ans, temp_cp, u);
            }

            return ans;
          }

          point_type f_u(const data_type &u, const data_type &v) const
          {
            // check to make sure have valid curve
            assert(degree_u()>=0);
            assert(degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_u()<1)
            {
              point_type ans;
              ans.setZero();
              return ans;
            }

            validate_u();

            return deriv_u->f( u, v );
          }

          void f_u( bezier<data_type, dim__, tol__> &bs_fu ) const
          {
            index_type i, n(degree_u()), m(degree_v());

            // check to make sure have valid curve
            assert(n>=0);
            assert(m>=0);

            // resize the surfaces
            bs_fu.resize(n-1, m);

            for (i=0; i<=m; ++i)
            {
              // create the control points for first derivative curve
              eli::geom::utility::bezier_p_control_point(bs_fu.B_u[i], B_u[i]);
            }
          }

          point_type f_v(const data_type &u, const data_type &v) const
          {
            // check to make sure have valid curve
            assert(degree_u()>=0);
            assert(degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_v()<1)
            {
              point_type ans;
              ans.setZero();
              return ans;
            }

            validate_v();

            return deriv_v->f( u, v );
          }

          void f_v( bezier<data_type, dim__, tol__> &bs_fv ) const
          {
            index_type i, n(degree_u()), m(degree_v());

            // check to make sure have valid curve
            assert(n>=0);
            assert(m>=0);

            // resize the surfaces
            bs_fv.resize(n, m-1);

            for (i=0; i<=n; ++i)
            {
              // create the control points for first derivative curve
              eli::geom::utility::bezier_p_control_point(bs_fv.B_v[i], B_v[i]);
            }
          }

          point_type f_uu(const data_type &u, const data_type &v) const
          {
            // check to make sure have valid curve
            assert(degree_u()>=0);
            assert(degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_u()<2)
            {
              point_type ans;
              ans.setZero();
              return ans;
            }

            validate_u();

            return deriv_u->f_u( u, v );
          }

          point_type f_uv(const data_type &u, const data_type &v) const
          {
            // check to make sure have valid curve
            assert(degree_u()>=0);
            assert(degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if ( (this->degree_u()<1) || (this->degree_v()<1) )
            {
              point_type ans;
              ans.setZero();
              return ans;
            }

            validate_u();

            return deriv_u->f_v( u, v );
          }

          point_type f_vv(const data_type &u, const data_type &v) const
          {
            // check to make sure have valid curve
            assert(degree_u()>=0);
            assert(degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_v()<2)
            {
              point_type ans;
              ans.setZero();
              return ans;
            }

            validate_v();

            return deriv_v->f_v( u, v );
          }

          point_type f_uuu(const data_type &u, const data_type &v) const
          {
            // check to make sure have valid curve
            assert(degree_u()>=0);
            assert(degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_u()<3)
            {
              point_type ans;
              ans.setZero();
              return ans;
            }

            validate_u();

            return deriv_u->f_uu( u, v );
        }

          point_type f_uuv(const data_type &u, const data_type &v) const
          {
            // check to make sure have valid curve
            assert(degree_u()>=0);
            assert(degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if ( (this->degree_u()<2) || (this->degree_v()<1) )
            {
              point_type ans;
              ans.setZero();
              return ans;
            }

            validate_u();

            return deriv_u->f_uv( u, v );
          }

          point_type f_uvv(const data_type &u, const data_type &v) const
          {
            // check to make sure have valid curve
            assert(degree_u()>=0);
            assert(degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if ( (this->degree_u()<1) || (this->degree_v()<2) )
            {
              point_type ans;
              ans.setZero();
              return ans;
            }

            validate_u();

            return deriv_u->f_vv( u, v );
          }

          point_type f_vvv(const data_type &u, const data_type &v) const
          {
            // check to make sure have valid curve
            assert(degree_u()>=0);
            assert(degree_v()>=0);

            // check to make sure given valid parametric value
            assert((u>=0) && (u<=1));
            assert((v>=0) && (v<=1));

            if (this->degree_v()<3)
            {
              point_type ans;
              ans.setZero();
              return ans;
            }

            validate_v();

            return deriv_v->f_vv( u, v );
          }

          point_type normal(const data_type &u, const data_type &v) const
          {
            point_type n=f_u(u, v).cross(f_v(u, v));
            data_type nlen(n.norm());
            tolerance_type tol;

            // If have degenerate surface try higher order terms.
            // Since want direction and don't care about magnitude,
            // can use du=dv=1  in Taylor serier expansion:
            // N(du, dv)=N+N_u*du+N_v*dv+1/2*(N_uu*du^2+2*N_uv*du*dv+N_vv*dv^2)+H.O.T.
            // see "Bezier Normal Vector Surface and Its Applications" by Yamaguchi
            if (tol.approximately_equal(nlen, 0))
            {
              point_type S_u, S_v, S_uu, S_uv, S_vv, N_u, N_v;
              data_type du(1), dv(1);

              // calculate Taylor series second term
              S_u=f_u(u, v);
              S_v=f_v(u, v);
              S_uu=f_uu(u, v);
              S_uv=f_uv(u, v);
              S_vv=f_vv(u, v);
              N_u=S_uu.cross(S_v)+S_u.cross(S_uv);
              N_v=S_uv.cross(S_v)+S_u.cross(S_vv);
              n=N_u*du+N_v*dv;
              nlen=n.norm();

              // if still have zero normal then calculate Taylor series third term
              if (tol.approximately_equal(nlen, 0))
              {
                point_type S_uuu, S_uuv, S_uvv, S_vvv, N_uu, N_uv, N_vv;

                S_uuu=f_uuu(u, v);
                S_uuv=f_uuv(u, v);
                S_uvv=f_uvv(u, v);
                S_vvv=f_vvv(u, v);
                N_uu=S_uuu.cross(S_v)+2*S_uu.cross(S_uv)+S_u.cross(S_uuv);
                N_uv=S_uuv.cross(S_v)+S_uu.cross(S_vv)+S_uv.cross(S_uv)+S_u.cross(S_uvv);
                N_vv=S_uvv.cross(S_v)+2*S_uv.cross(S_vv)+S_u.cross(S_vvv);
                n=0.5*(N_uu*du*du+2*N_uv*du*dv+N_vv*dv*dv);
                nlen=n.norm();

                // if still get zero normal then give up
                if (tol.approximately_equal(nlen, 0))
                {
                  n.setZero();
                  nlen=1;
                }
              }
            }

            n/=nlen;
            return n;
          }

          void f_pt_normal(const data_type &u, const data_type &v, point_type &pt, point_type &norm ) const
          {
            pt = f( u, v );
            norm = normal( u, v );
          }

          void f_pt_derivs(const data_type &u, const data_type &v, point_type &pt, point_type &pt_u, point_type &pt_v ) const
          {
            pt = f( u, v );
            pt_u = f_u( u, v );
            pt_v = f_v( u, v );
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
              eli::geom::utility::bezier_promote_control_points(tmp_cp, current_row[i]);
              B_u[i]=tmp_cp;
            }
            invalidate_deriv();
          }

          void promote_u_to(index_type target_degree)
          {
              typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_row_type;
              typedef std::vector<control_row_type, Eigen::aligned_allocator<control_row_type> > control_row_collection_type;

              index_type i, n(degree_u()), m(degree_v());
              control_row_collection_type current_row(m+1, control_row_type(n+1, dim__));

              // copy the control rows
              for (i=0; i<=m; ++i)
              {
                current_row[i]=B_u[i];
              }

              // resize current surface
              resize(target_degree, m);

              // set the new control points
              control_row_type tmp_cp(target_degree+1, dim__);
              for (i=0; i<=m; ++i)
              {
                eli::geom::utility::bezier_promote_control_points_to(tmp_cp, current_row[i]);
                B_u[i]=tmp_cp;
              }
              invalidate_deriv();
          }

          void promote_v()
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_col_type;
            typedef std::vector<control_col_type, Eigen::aligned_allocator<control_col_type> > control_col_collection_type;

            index_type i, n(degree_u()), m(degree_v());
            control_col_collection_type current_col(n+1, control_col_type(m+1, dim__));

            // copy the control cols
            for (i=0; i<=n; ++i)
            {
              current_col[i]=B_v[i];
            }

            // resize current surface
            resize(n, m+1);

            // set the new control points
            control_col_type tmp_cp(m+2, dim__);
            for (i=0; i<=n; ++i)
            {
              eli::geom::utility::bezier_promote_control_points(tmp_cp, current_col[i]);
              B_v[i]=tmp_cp;
            }
            invalidate_deriv();
          }

          void promote_v_to(index_type target_degree)
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_col_type;
            typedef std::vector<control_col_type, Eigen::aligned_allocator<control_col_type> > control_col_collection_type;

            index_type i, n(degree_u()), m(degree_v());
            control_col_collection_type current_col(n+1, control_col_type(m+1, dim__));

            // copy the control cols
            for (i=0; i<=n; ++i)
            {
              current_col[i]=B_v[i];
            }

            // resize current surface
            resize(n, target_degree);

            // set the new control points
            control_col_type tmp_cp(target_degree+1, dim__);
            for (i=0; i<=n; ++i)
            {
              eli::geom::utility::bezier_promote_control_points_to(tmp_cp, current_col[i]);
              B_v[i]=tmp_cp;
            }
            invalidate_deriv();
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
              eli::geom::utility::bezier_demote_control_points(tmp_cp, current_row[i], ncon);
              B_u[i]=tmp_cp;
            }
            invalidate_deriv();
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
              eli::geom::utility::bezier_demote_control_points(tmp_cp, current_col[i], ncon);
              B_v[i]=tmp_cp;
            }
            invalidate_deriv();
            return true;
          }

          void to_cubic_u()
          {
              typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_row_type;
              typedef std::vector<control_row_type, Eigen::aligned_allocator<control_row_type> > control_row_collection_type;

              index_type i, n(degree_u()), m(degree_v());
              control_row_collection_type current_row(m+1, control_row_type(n+1, dim__));

              // copy the control rows
              for (i=0; i<=m; ++i)
              {
                current_row[i]=B_u[i];
              }

              // resize current surface
              resize(3, m);

              // set the new control points
              control_row_type tmp_cp(4, dim__);
              for (i=0; i<=m; ++i)
              {
                eli::geom::utility::bezier_control_points_to_cubic(tmp_cp, current_row[i]);
                B_u[i]=tmp_cp;
              }
              invalidate_deriv();
          }

          void to_cubic_v()
          {
            typedef Eigen::Matrix<data_type, Eigen::Dynamic, dim__> control_col_type;
            typedef std::vector<control_col_type, Eigen::aligned_allocator<control_col_type> > control_col_collection_type;

            index_type i, n(degree_u()), m(degree_v());
            control_col_collection_type current_col(n+1, control_col_type(m+1, dim__));

            // copy the control rows
            for (i=0; i<=n; ++i)
            {
              current_col[i]=B_v[i];
            }

            // resize current surface
            resize(n, 3);

            // set the new control points
            control_col_type tmp_cp(4, dim__);
            for (i=0; i<=n; ++i)
            {
              eli::geom::utility::bezier_control_points_to_cubic(tmp_cp, current_col[i]);
              B_v[i]=tmp_cp;
            }
            invalidate_deriv();
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
              eli::geom::utility::bezier_split_control_points(cp_lo, cp_hi, B_u[j], u0);
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
              eli::geom::utility::bezier_split_control_points(cp_lo, cp_hi, B_v[i], v0);
              for (j=0; j<=m; ++j)
              {
                bs_lo.set_control_point(cp_lo.row(j), i, j);
                bs_hi.set_control_point(cp_hi.row(j), i, j);
              }
            }
          }

          data_type eqp_distance_bound(const bezier<data_type, dim__, tol__> &bs) const
          {
            typedef bezier<data_type, dim__, tol__> surf_type;

            // Make working copies of surfaces.
            surf_type bsa(*this);
            surf_type bsb(bs);

            // Find maximum common order.
            index_type n(bsa.degree_u()), m(bsa.degree_v());
            if(bsb.degree_u() > n)
            {
              n = bsb.degree_u();
            }

            if(bsb.degree_v() > m)
            {
              m = bsb.degree_v();
            }

            // Promote both to max common order.
            bsa.promote_u_to(n);
            bsa.promote_v_to(m);

            bsb.promote_u_to(n);
            bsb.promote_v_to(m);

            // Find maximum distance between control points.
            index_type i, j;
            data_type d, maxd(0);
            for (i=0; i<=n; ++i)
            {
              for (j=0; j<=m; ++j)
              {
                d = (bsa.get_control_point(i, j) - bsb.get_control_point(i, j)).norm();
                if(d > maxd)
                {
                  maxd=d;
                }
              }
            }

            return maxd;
          }


          void product( const bezier<data_type, dim__, tol__> &a, const bezier<data_type, dim__, tol__> &b)
          {
            index_type i, mu(a.degree_u()), mv(a.degree_v()), nu(b.degree_u()), nv(b.degree_v());

            u_control_point_matrix_container a_u, b_u;
            v_control_point_matrix_container a_v, b_v;
            control_point_container a_data, b_data;

            a_data = a.point_data;
            set_Bs(a_u, a_v, a_data, mu, mv);
            b_data = b.point_data;
            set_Bs( b_u, b_v, b_data, nu, nv );

            resize( mu + nu, mv + nv );

            for ( size_t i=1; i<point_data.size(); i++)
            {
              point_data[i] = 0.0;
            }

            for (i=0; i<=mu; ++i)
            {
              eli::geom::utility::bezier_control_points_to_scaled_bezier( a_v[i] );
            }

            for (i=0; i<=mv; ++i)
            {
              eli::geom::utility::bezier_control_points_to_scaled_bezier( a_u[i] );
            }

            for (i=0; i<=nu; ++i)
            {
              eli::geom::utility::bezier_control_points_to_scaled_bezier( b_v[i] );
            }

            for (i=0; i<=nv; ++i)
            {
              eli::geom::utility::bezier_control_points_to_scaled_bezier( b_u[i] );
            }

            index_type ia, ib, ic;
            index_type ja, jb, jc;

            for (ia = 0; ia <= mu; ia++ )
            {
              for (ib = 0; ib <= nu; ib++ )
              {
                ic = ia + ib;

                for (ja = 0; ja <= mv; ja++ )
                {
                  for (jb = 0; jb <= nv; jb++ )
                  {
                    jc = ja + jb;

                    B_u[jc].row(ic) = B_u[jc].row(ic) + a_u[ja].row(ia).cwiseProduct( b_u[jb].row(ib) );
                  }
                }
              }
            }

            for (i=0; i<=mu+nu; ++i)
            {
              eli::geom::utility::scaled_bezier_to_control_points_bezier( B_v[i] );
            }

            for (i=0; i<=mv+nv; ++i)
            {
              eli::geom::utility::scaled_bezier_to_control_points_bezier( B_u[i] );
            }
            invalidate_deriv();
          }

          void sum( const bezier<data_type, dim__, tol__> &a, const bezier<data_type, dim__, tol__> &b)
          {
            typedef bezier<data_type, dim__, tol__> surf_type;

            // Make working copies of surfaces.
            surf_type bsa(a);
            surf_type bsb(b);

            // Find maximum common order.
            index_type n(bsa.degree_u()), m(bsa.degree_v());
            if(bsb.degree_u() > n)
            {
              n = bsb.degree_u();
            }

            if(bsb.degree_v() > m)
            {
              m = bsb.degree_v();
            }

            // Promote both to max common order.
            bsa.promote_u_to(n);
            bsa.promote_v_to(m);

            bsb.promote_u_to(n);
            bsb.promote_v_to(m);

            resize( n, m );

            // Find maximum distance between control points.
            index_type i, j;
            for (i=0; i<=n; ++i)
            {
              for (j=0; j<=m; ++j)
              {
                set_control_point( bsa.get_control_point(i, j) + bsb.get_control_point(i, j), i, j );
              }
            }
          }

          onedbezsurf sumcompsurf() const
          {
            onedbezsurf retsurf;
            index_type n(degree_u()), m(degree_v());

            retsurf.resize( n, m );

            index_type i, j, k;
            for (i=0; i<=n; ++i)
            {
              for (j=0; j<=m; ++j)
              {
                data_type d = 0;
                point_type p = get_control_point( i, j );

                for (k=0; k<dim__; ++k)
                {
                  d += p(k);
                }

                typename onedbezsurf::point_type pd;
                pd(0) = d;

                retsurf.set_control_point( pd, i, j );
              }
            }

            return retsurf;
          }

          onedbezsurf mindistsurf( const point_type & pt ) const
          {
            onedbezsurf retsurf;
            typedef bezier<data_type, dim__, tol__> surf_type;

            surf_type d(*this);

            d.translate( -pt );

            validate_u();

            validate_v();

            surf_type prod1, prod2;

            prod1.product( *deriv_u, d );
            prod2.product( *deriv_v, d );

            onedbezsurf dot1, dot2;
            dot1 = prod1.sumcompsurf();
            dot2 = prod2.sumcompsurf();

            onedbezsurf sq1, sq2;

            sq1.product( dot1, dot1 );
            sq2.product( dot2, dot2 );

            retsurf.sum( sq1, sq2 );

            return retsurf;
          }

          bool allpos( const data_type &smallpos ) const
          {
            index_type i, j, k;
            index_type n(degree_u()), m(degree_v());

            for (i=0; i<=n; ++i)
            {
              for (j=0; j<=m; ++j)
              {
                point_type p = get_control_point( i, j );

                for (k=0; k<dim__; ++k)
                {
                  if ( p(k) <= smallpos )
                  {
                    return false;
                  }
                }
              }
            }
            return true;
          }

        private:

          static void resize( u_control_point_matrix_container &uvec, v_control_point_matrix_container &vvec, control_point_container &data, const index_type &u_dim, const index_type &v_dim)
          {
            // allocate the control points
            data.resize(dim__*(u_dim+1)*(v_dim+1));

            // set the B_u and B_v maps
            set_Bs( uvec, vvec, data, u_dim, v_dim);
          }

          void set_Bs(index_type n, index_type m)
          {
            set_Bs( B_u, B_v, point_data, n, m );
          }

          static void set_Bs( u_control_point_matrix_container &uvec, v_control_point_matrix_container &vvec, control_point_container &data, index_type n, index_type m)
          {
            // allocate vectors of control point maps
            uvec.resize(m+1, control_point_matrix_type(nullptr, m+1, dim__, Eigen::Stride<1, dim__>()));
            for (index_type j=0; j<=m; ++j)
            {
#ifdef ELI_NO_VECTOR_DATA
              new (&(uvec.at(0))+j) control_point_matrix_type(&(data.at(0))+j*(n+1)*dim__, n+1, dim__, Eigen::Stride<1, dim__>());
#else
              new (uvec.data()+j) control_point_matrix_type(data.data()+j*(n+1)*dim__, n+1, dim__, Eigen::Stride<1, dim__>());
#endif
            }

            vvec.resize(n+1, v_dir_control_point_matrix_type(nullptr, n+1, dim__, Eigen::Stride<1, Eigen::Dynamic>(1, (n+1)*dim__)));
            for (index_type i=0; i<=n; ++i)
            {
#ifdef ELI_NO_VECTOR_DATA
              new (&(vvec.at(0))+i) v_dir_control_point_matrix_type(&(data.at(0))+i*dim__, m+1, dim__, Eigen::Stride<1, Eigen::Dynamic>(1, (n+1)*dim__));
#else
              new (vvec.data()+i) v_dir_control_point_matrix_type(data.data()+i*dim__, m+1, dim__, Eigen::Stride<1, Eigen::Dynamic>(1, (n+1)*dim__));
#endif
            }
          }

          void invalidate_deriv()
          {
            if ( deriv_u )
            {
              delete deriv_u;
              deriv_u = NULL;
            }

            if ( deriv_v )
            {
              delete deriv_v;
              deriv_v = NULL;
            }
          }

          void validate_u() const
          {
            if ( !deriv_u )
            {
              deriv_u = new bezier<data_type, dim__>();
              f_u( *deriv_u );
            }
          }

          void validate_v() const
          {
            if ( !deriv_v )
            {
              deriv_v = new bezier<data_type, dim__>();
              f_v( *deriv_v );
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
