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

#ifndef eli_geom_curve_piecewise_circle_creator_hpp
#define eli_geom_curve_piecewise_circle_creator_hpp

#include <iterator>

#include "eli/geom/point/distance.hpp"

#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace utility
      {
        template<typename data__>
        void calculate_circle(data__ &rad, data__ &x0, data__ &y0, const data__ &xa, const data__ &ya,
                              const data__ &xb, const data__ &yb, const data__ &xc, const data__ &yc)
        {
          data__ denom, xamb, yamb, xcma, ycma, xbmc, ybmc, alen2, blen2, clen2;

          xamb=(xa-xb);
          yamb=(ya-yb);
          xcma=(xc-xa);
          ycma=(yc-ya);
          xbmc=(xb-xc);
          ybmc=(yb-yc);
          alen2=xa*xa+ya*ya;
          blen2=xb*xb+yb*yb;
          clen2=xc*xc+yc*yc;

          denom=2*(xc*yamb+xb*ycma+xa*ybmc);
          rad=std::sqrt((xamb*xamb+yamb*yamb)*(xcma*xcma+ycma*ycma)*(xbmc*xbmc+ybmc*ybmc))/denom;
          x0= (clen2*yamb+blen2*ycma+alen2*ybmc)/denom;
          y0=-(clen2*xamb+blen2*xcma+alen2*xbmc)/denom;
        }
      }

      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_circle_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef data__  data_type;
          typedef int index_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef tol__ tolerance_type;

          piecewise_circle_creator() : piecewise_creator_base<data_type, dim__, tolerance_type>(4, 0), radius(1) {}
          piecewise_circle_creator(const index_type &ns) : piecewise_creator_base<data_type, dim__, tolerance_type>(ns, 0), radius(1) {}
          piecewise_circle_creator(const piecewise_circle_creator<data_type, dim__, tolerance_type> &pcc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(pcc), radius(pcc.radius), x(pcc.x), y(pcc.y) {}

          void set_radius(const data_type &r)
          {
            tolerance_type tol;

            if (tol.approximately_equal(r, 0))
            {
              radius=0;
            }
            else if (r<0)
            {
              assert(false);
            }
            else
            {
              radius=r;
            }
          }
          data_type get_radius() const {return radius;}

          void set_origin(const point_type &orig) {origin=orig;}
          point_type get_origin() const {return origin;}

          void set_xy_directions(const point_type &xdir, const point_type &ydir)
          {
            tolerance_type tol;

            if (tol.approximately_equal(xdir.dot(ydir), 0))
            {
              x=xdir;
              y=ydir;
            }
            else
            {
              assert(false);
            }
          }
          void get_xy_directions(point_type &xdir, point_type &ydir) const
          {
            xdir=x;
            ydir=y;
          }

          void set(const point_type &start, const point_type &orig)
          {
            // need to be on the same z-plane
            if (dim__!=2)
            {
              if (start.col(2)!=orig.col(2))
              {
                assert(false);
                return;
              }
            }

            // set radius & origin
            data_type r;
            eli::geom::point::distance(r, start, orig);
            set_radius(r);
            origin=orig;

            if (radius==0)
            {
              x.setZero();
              x(0)=1;
              y.setZero();
              y(1)=1;
            }
            else
            {
              x=start-orig;
              x.normalize();
              y=origin;
              y.x()=-x.y();
              y.y()=x.x();
            }
          }

          void set(const point_type &start, const point_type &orig, const point_type &normal)
          {
            if (dim__==2)
            {
              set(start, origin);
              return;
            }

            // set radius & origin
            data_type r;
            eli::geom::point::distance(r, start, orig);
            set_radius(r);
            origin=orig;

            if (radius==0)
            {
              x.setZero();
              x(0)=1;
              y.setZero();
              y(1)=1;
            }
            else
            {
              x=start-origin;
              x.normalize();
              y << normal(1)*x(2)-normal(2)*x(1),
                   normal(2)*x(0)-normal(0)*x(2),
                   normal(0)*x(1)-normal(1)*x(0);
              y.normalize();

#ifdef DEBUG
              point_type n(normal);
              tolerance_type tol;
              n.normalize();
              assert(tol.approximately_equal(x.dot(n), 0));
#endif
            }
          }

          void set_3pt(const point_type &start, const point_type &middle, const point_type &end)
          {
            tolerance_type tol;

            // cannot have points in same location
            if ( tol.approximately_equal(start, end)
              || tol.approximately_equal(start, middle)
              || tol.approximately_equal(middle, end) )
            {
              assert(false);
              return;
            }

            if (dim__==2)
            {
              point_type orig;
              data_type r;

              eli::geom::curve::utility::calculate_circle(r, orig.x(), orig.y(), start.x(), start.y(),
                                                          middle.x(), middle.y(), end.x(), end.y());
              set(start, orig);
            }
            else
            {
              point_type orig, normal, xtmp, ytmp;

              // set the temporary x-direction as from start to middle
              xtmp=middle-start;
              xtmp.normalize();

              // find normal
              ytmp=end-start;
              normal << xtmp(1)*ytmp(2)-xtmp(2)*ytmp(1),
                        xtmp(2)*ytmp(0)-xtmp(0)*ytmp(2),
                        xtmp(0)*ytmp(1)-xtmp(1)*ytmp(0);
              normal.normalize();

              // get the temporary y-direction
              ytmp << normal(1)*xtmp(2)-normal(2)*xtmp(1),
                      normal(2)*xtmp(0)-normal(0)*xtmp(2),
                      normal(0)*xtmp(1)-normal(1)*xtmp(0);
              ytmp.normalize();

              // transform points into temp x-y coordinates
              data_type r, x0, y0, xa, ya, xb, yb, xc, yc;
              point_type tmp;

              xa=0;
              ya=0;
              tmp=middle-start;
              xb=tmp.dot(xtmp);
              yb=tmp.dot(ytmp);
              tmp=end-start;
              xc=tmp.dot(xtmp);
              yc=tmp.dot(ytmp);

              // use formula to calculate radius and origin (will be in temp x-y coordinates)
              eli::geom::curve::utility::calculate_circle(r, x0, y0, xa, ya, xb, yb, xc, yc);

              // calculate the actual coordinates of the origin and set values
              orig=x0*xtmp+y0*ytmp+start;
              set(start, orig, normal);
            }
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;
            typedef typename curve_type::control_point_type control_point_type;

            curve_type c(3);
            control_point_type cp[4];
            error_code err;
            data_type k;
            index_type i;

            // can only handle 4 segments for now
            if (this->get_number_segments()!=4)
            {
              assert(false);
              return false;
            }

            // set the start parameter
            pc.set_t0(this->get_t0());

            // set up for curve creation
            pc.clear();
            k=4*(eli::constants::math<data_type>::sqrt_two()-1)/3;

            // set 1st quadrant curve
            cp[0]=radius*x+origin;
            cp[1]=radius*(x+k*y)+origin;
            cp[2]=radius*(k*x+y)+origin;
            cp[3]=radius*y+origin;
            for (i=0; i<4; ++i)
            {
              c.set_control_point(cp[i], i);
            }
            err=pc.push_back(c, this->get_segment_dt(0));
            if (err!=piecewise_curve_type::NO_ERROR)
            {
              pc.clear();
              pc.set_t0(0);
              return false;
            }

            // set 2nd quadrant curve
            cp[0]=radius*y+origin;
            cp[1]=radius*(y-k*x)+origin;
            cp[2]=radius*(k*y-x)+origin;
            cp[3]=-radius*x+origin;
            for (i=0; i<4; ++i)
            {
              c.set_control_point(cp[i], i);
            }
            err=pc.push_back(c, this->get_segment_dt(1));
            if (err!=piecewise_curve_type::NO_ERROR)
            {
              pc.clear();
              pc.set_t0(0);
              return false;
            }

            // set 3rd quadrant curve
            cp[0]=-radius*x+origin;
            cp[1]=-radius*(x+k*y)+origin;
            cp[2]=-radius*(k*x+y)+origin;
            cp[3]=-radius*y+origin;
            for (i=0; i<4; ++i)
            {
              c.set_control_point(cp[i], i);
            }
            err=pc.push_back(c, this->get_segment_dt(2));
            if (err!=piecewise_curve_type::NO_ERROR)
            {
              pc.clear();
              pc.set_t0(0);
              return false;
            }

            // set 4th quadrant curve
            cp[0]=-radius*y+origin;
            cp[1]=-radius*(y-k*x)+origin;
            cp[2]=-radius*(k*y-x)+origin;
            cp[3]=radius*x+origin;
            for (i=0; i<4; ++i)
            {
              c.set_control_point(cp[i], i);
            }
            err=pc.push_back(c, this->get_segment_dt(3));
            if (err!=piecewise_curve_type::NO_ERROR)
            {
              pc.clear();
              pc.set_t0(0);
              return false;
            }

            return true;
          }

        private:
          point_type origin, x, y;
          data_type radius;
      };
    }
  }
}
#endif
