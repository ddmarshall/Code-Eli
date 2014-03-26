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

#ifndef eli_geom_curve_piecewise_general_creator_hpp
#define eli_geom_curve_piecewise_general_creator_hpp

#include <vector>

#include "eli/geom/general/continuity.hpp"

#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {

      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_general_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;

          enum joint_conditions {NOT_SET=0,
                                 VALUE_SET=1/*,
                                 DIRECTION_SET=2*/};

          enum joint_continuity {C0=general::C0,
                                 C1=general::C1,
                                 C2=general::C2};

          class joint_data
          {
            public:
              void set_f(const point_type &p) {f=p;}
              point_type get_f() const {return f;}

              void unset_fp() {use_fp=NOT_SET;}
              void set_fp(const point_type &p)
              {
                use_fp=VALUE_SET;
                fp=p;
              }
              joint_conditions fp_state() const {return use_fp;}
              point_type get_fp() const {return fp;}

              void unset_fpp() {use_fpp=NOT_SET;}
              void set_fpp(const point_type &p)
              {
                use_fpp=VALUE_SET;
                fpp=p;
              }
              joint_conditions fpp_state() const {return use_fpp;}
              point_type get_fpp() const {return fpp;}

            private:
              point_type f;
              joint_conditions use_fp;
              point_type fp;
              joint_conditions use_fpp;
              point_type fpp;
          };

        public:
          piecewise_general_creator()
            : piecewise_creator_base<data_type, dim__, tolerance_type>(1, 0), joints(2),
              joint_cont(1), max_degree(1), closed(false)
          {
          }
          piecewise_general_creator(const index_type &ns)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(ns, 0), joints(ns+1),
              joint_cont(ns), max_degree(ns), closed_cont(C0), closed(false)
          {
          }
          piecewise_general_creator(const piecewise_general_creator<data_type, dim__, tolerance_type> &pcc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(pcc), joints(pcc.joints),
              joint_cont(pcc.joint_cont), max_degree(pcc.max_degree()), closed(pcc.closed)
          {
          }
          virtual ~piecewise_general_creator()
          {
          };

          void set_closed_continuity(const joint_continuity &cnt)
          {
            closed=true;
            closed_cont=cnt;
          }
          void set_open() {closed=false;}
          bool is_closed() const {return closed;}
          bool is_open() const {return !closed;}


          bool set_conditions(const std::vector<joint_data> &jnts, const std::vector<joint_continuity> &cnt,
                              const std::vector<index_type> &maxd, bool cl=false)
          {
            index_type nsegs(static_cast<index_type>(maxd.size()));

            // ensure input vectors are correct size
            if (jnts.size()!=(nsegs+1))
              return false;
            if (nsegs!=(cnt.size()+1))
              return false;

            // reset the number of segments
            set_number_segments(nsegs);

            joints=jnts;
            joint_cont=cnt;
            closed=cl;

            return true;
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;
            typedef typename curve_type::control_point_type control_point_type;

            pc.clear();


            return false;
          }

        private:
          std::vector<joint_data> joints;
          std::vector<joint_continuity> joint_cont;
          std::vector<index_type> max_degree;
          joint_continuity closed_cont;
          bool closed;
      };
    }
  }
}
#endif
