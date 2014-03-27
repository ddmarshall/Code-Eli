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

          piecewise_four_digit_creator() : piecewise_airfoil_creator_base<data_type, dim__, tolerance_type>(0) {}
          piecewise_four_digit_creator(const piecewise_four_digit_creator<data_type, dim__, tolerance_type> &ppc)
            : piecewise_airfoil_creator_base<data_type, dim__, tolerance_type>(ppc) {}
          ~piecewise_four_digit_creator() {}

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
#if 0
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;

            pc.clear();

            curve_type c(1);
            error_code err;
            index_type nsegs(this->get_number_segments());

            // do sanity check
            if (corner.size()!=static_cast<size_t>(nsegs+1))
            {
              assert(false);
              return false;
            }

            // set the start parameter
            pc.set_t0(this->get_t0());

            // set the first n edges
            for (index_type i=0; i<nsegs; ++i)
            {
              c.set_control_point(corner[i], 0);
              c.set_control_point(corner[i+1], 1);
              err=pc.push_back(c, this->get_segment_dt(i));
              if (err!=piecewise_curve_type::NO_ERRORS)
              {
                pc.clear();
                pc.set_t0(0);
                return false;
              }
            }
#endif
            
            return false;
          }

        private:
      };
    }
  }
}
#endif
