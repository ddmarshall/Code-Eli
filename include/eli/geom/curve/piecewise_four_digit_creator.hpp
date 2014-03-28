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

#ifndef eli_geom_curve_piecewise_four_digit_creator_hpp
#define eli_geom_curve_piecewise_four_digit_creator_hpp

#ifndef eli_geom_curve_piecewise_airfoil_creator_base_hpp
#define eli_geom_curve_piecewise_airfoil_creator_base_hpp

#include <vector>

#ifdef Success // X11 #define collides with Eigen
#undef Success
#endif

#include "Eigen/Eigen"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/pseudo/four_digit.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_airfoil_creator_base
      {
      public:
        typedef data__  data_type;
        typedef Eigen::Matrix<data_type, 1, dim__> point_type;
        typedef typename point_type::Index index_type;
        typedef tol__ tolerance_type;
        
      public:
        piecewise_airfoil_creator_base(const data_type &tt0) : dt(4), t0(tt0)
        {
          for (index_type i=0; i<static_cast<index_type>(dt.size()); ++i)
            dt[i]=1;
        }
        piecewise_airfoil_creator_base(const piecewise_airfoil_creator_base<data_type, dim__, tolerance_type> &pac) : dt(pac.dt), t0(pac.t0) {}
        virtual ~piecewise_airfoil_creator_base() {}
        
        index_type get_number_segments() const
        {
          return static_cast<index_type>(dt.size());
        }
        
        void set_t0(const data_type &tt0) {t0=tt0;}
        data_type get_t0() const {return t0;}
        
        void set_segment_dt(const data_type &dtt, const index_type &i)
        {
          if ((dtt>0) && (i>=0) && (i<static_cast<index_type>(dt.size())))
            dt[i]=dtt;
          else
            assert(false);
        }
        data_type get_segment_dt(const index_type &i) const
        {
          if ((i<0) || (i>=static_cast<index_type>(dt.size())))
          {
            assert(false);
            return static_cast<data_type>(-1);
          }
          
          return dt[i];
        }
        
        virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const = 0;

      private:
        std::vector<data_type> dt;
        data_type t0;
      };
    }
  }
}
#endif

#include "Eigen/Eigen"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_four_digit_creator : public piecewise_airfoil_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_airfoil_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef eli::geom::curve::pseudo::four_digit<data_type> airfoil_type;
          typedef typename airfoil_type::point_type af_point_type;

          piecewise_four_digit_creator() : piecewise_airfoil_creator_base<data_type, dim__, tolerance_type>(0) {}
          piecewise_four_digit_creator(const piecewise_four_digit_creator<data_type, dim__, tolerance_type> &ppc)
            : piecewise_airfoil_creator_base<data_type, dim__, tolerance_type>(ppc) {}
          ~piecewise_four_digit_creator() {}

          void set_sharp_trailing_edge(bool fl)
          {
            af.set_sharp_trailing_edge(fl);
          }
          bool sharp_trailing_edge() const
          {
            return af.sharp_trailing_edge();
          }

          bool set_thickness(const data_type &t)
          {
            return af.set_thickness(t);
          }

          data_type get_thickness() const
          {
            return af.get_thickness();
          }

          bool set_camber(const data_type &cam, const data_type &cam_loc)
          {
            return af.set_camber(cam, cam_loc);
          }

          data_type get_maximum_camber() const
          {
            return af.get_maximum_camber();
          }

          data_type get_maximum_camber_location() const
          {
            return af.get_maximum_camber_location();
          }

          bool set_name(const std::string &name)
          {
            return af.set_name(name);
          }

          std::string get_name() const
          {
            return af.get_name();
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;
            typedef typename curve_type::fit_container_type fit_container_type;
            typedef typename curve_type::dimension_type dimension_type;

            af_point_type temppt;

            std::vector<point_type, Eigen::aligned_allocator<point_type> > pts;


            index_type i;
            index_type nseg(this->get_number_segments());

            // Number of sample points per segment.
            index_type nref=25;

            // Number of points evaluated around airfoil.
            index_type npt = nref*nseg+1;
            pts.resize(npt);
            // Leading edge point index.
            index_type ile(1+(npt-1)/2);

            // Set up initial parameter and parameter step.
            data_type xi(af.get_u_min());
            data_type dxi((af.get_u_max()-af.get_u_min())/(npt-1));
            std::vector< data_type > xis;
            xis.resize(npt);

            // Evaluate airfoil.
            for( i = 0; i < npt; i++ )
            {
              if(i==ile)
              {
                xi=0;  // Force exact floating point value for le.
              }
              else if(i==npt-1)
              {
                xi=af.get_u_max();  // Force exact floating point value for te.
              }

              temppt = af.f(xi);
              pts[i] = point_type(temppt.x(), temppt.y(), 0);
              xis[i] = xi;
              xi += dxi;
            }

            pc.clear();
            pc.set_t0(this->get_t0());

            index_type istart(0), iend(nref);


            for( i = 0; i < nseg; i++ )
            {

              // set up fit container
              fit_container_type fcon;

              fcon.set_points(pts.begin()+istart, pts.begin()+iend+1);
              fcon.add_start_C0_constraint();
              fcon.add_end_C0_constraint();


              // do fit
              dimension_type dim(10);
              curve_type bez;

              bez.fit(fcon, dim);

              // Push back to piecewise curve.
              error_code err = pc.push_back(bez, this->get_segment_dt(i));
              assert(err==piecewise_curve_type::NO_ERRORS);

              istart = iend;

              iend = iend + nref;
            }

            return false;
          }

        private:
          airfoil_type af;
      };
    }
  }
}
#endif
