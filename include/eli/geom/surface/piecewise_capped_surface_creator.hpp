/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_surface_piecewise_body_of_revolution_creator_hpp
#define eli_geom_surface_piecewise_body_of_revolution_creator_hpp

#include <list>
#include <iterator>

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise.hpp"

#include "eli/geom/surface/bezier.hpp"
#include "eli/geom/surface/piecewise.hpp"

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class capped_surface_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          enum edge_split_identifier
          {
            SPLIT_NONE,
            SPLIT_UMIN,
            SPLIT_UMAX,
            SPLIT_VMIN,
            SPLIT_VMAX
          };

          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef typename base_class_type::piecewise_surface_type piecewise_surface_type;

          capped_surface_creator()
            : piecewise_creator_base<data__, dim__, tol__>(0, 0), split_param(0), edge_to_split(SPLIT_NONE)
          {
          }
          capped_surface_creator(const data_type &uu0, const data_type &sp, edge_split_identifier esi)
            : piecewise_creator_base<data__, dim__, tol__>(0, 0), edge_curve
          {
          }
          capped_surface_creator(const general_skinning_surface_creator<data_type, dim__, tolerance_type> & gs)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(gs), ribs(gs.ribs),
              max_degree(gs.max_degree), closed(gs.closed)
          {
          }
          virtual ~capped_surface_creator()
          {
          }

          bool set_conditions(const piecewise_curve_type &rbs, const std::vector<index_type> &maxd, bool cl=false)
          {
          }

          virtual bool create(piecewise_surface_type &ps) const
          {
            typedef typename eli::geom::curve::piecewise<eli::geom::curve::bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            // need to implement
            assert(false);

            return false;
          }
#if 0
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_surface_type;
            typedef typename piecewise_surface_type::surface_type surface_type;

            index_type nribs(this->get_number_u_segments()+1), i, j;
            std::vector<index_type> seg_degree(nribs-1);
            std::vector<rib_data_type> rib_states(ribs);
            tolerance_type tol;

            // FIX: Should be able to handle closed surfaces
            assert(!closed);

            // FIX: Need to be able to handle v-direction discontinuous fu and fuu specifications

            // reset the incoming piecewise surface
            ps.clear();

            // split ribs so have same number of curves (with same joint parameters) for all ribs and get degree
            index_type njoints(this->get_number_v_segments()+1);
            std::vector<data_type> joints(njoints);
            std::vector<index_type> max_jdegs(njoints-1,0);

            joints[0]=this->get_v0();
            for (j=0; j<(njoints-1); ++j)
            {
              joints[j+1]=joints[j]+this->get_segment_dv(j);
            }

            for (i=0; i<nribs; ++i)
            {
              std::vector<index_type> jdegs;
              rib_states[i].split(joints.begin(), joints.end(), std::back_inserter(jdegs));
              for (j=0; j<(njoints-1); ++j)
              {
                if (jdegs[j]>max_jdegs[j])
                {
                  max_jdegs[j]=jdegs[j];
                }
              }
            }

            // set degree in u-direction for each rib segment strip
            for (i=0; i<nribs; ++i)
            {
              rib_states[i].promote(max_jdegs.begin(), max_jdegs.end());
            }

            // resize the piecewise surface
            index_type u, v, nu(nribs-1), nv(njoints-1);

            ps.init_uv(this->du_begin(), this->du_end(), this->dv_begin(), this->dv_end(), this->get_u0(), this->get_v0());

            // build segments based on rib information
            // here should have everything to make an nribs x njoints piecewise surface with all
            // of the j-degrees matching in the u-direction so that can use general curve creator
            // techniques to create control points
            for (v=0; v<nv; ++v)
            {
              typedef eli::geom::curve::piecewise_general_creator<data_type, dim__, tolerance_type> piecewise_curve_creator_type;
              typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
              typedef typename piecewise_curve_type::curve_type curve_type;

              std::vector<typename piecewise_curve_creator_type::joint_data> joints(nu+1);
              piecewise_curve_creator_type gc;
              piecewise_curve_type c;

              std::vector<surface_type> surfs(nu);

              for (j=0; j<=max_jdegs[v]; ++j)
              {
                // cycle through each rib to set corresponding joint info
                for (u=0; u<=nu; ++u)
                {
                  curve_type jcrv;

                  joints[u].set_continuity(static_cast<typename piecewise_curve_creator_type::joint_continuity>(rib_states[u].get_continuity()));

                  rib_states[u].get_f().get(jcrv, v);
                  joints[u].set_f(jcrv.get_control_point(j));
                  if (rib_states[u].use_left_fp())
                  {
                    rib_states[u].get_left_fp().get(jcrv, v);
                    joints[u].set_left_fp(jcrv.get_control_point(j));
                  }
                  if (rib_states[u].use_right_fp())
                  {
                    rib_states[u].get_right_fp().get(jcrv, v);
                    joints[u].set_right_fp(jcrv.get_control_point(j));
                  }
                  if (rib_states[u].use_left_fpp())
                  {
                    rib_states[u].get_left_fpp().get(jcrv, v);
                    joints[u].set_left_fpp(jcrv.get_control_point(j));
                  }
                  if (rib_states[u].use_right_fpp())
                  {
                    rib_states[u].get_right_fpp().get(jcrv, v);
                    joints[u].set_right_fpp(jcrv.get_control_point(j));
                  }
                }

                // set the conditions for the curve creator
                bool rtn_flag(gc.set_conditions(joints, max_degree, closed));
                if (!rtn_flag)
                {
                  return false;
                }

                // set the parameterizations and create curve
                gc.set_t0(this->get_u0());
                for (u=0; u<nu; ++u)
                {
                  gc.set_segment_dt(this->get_segment_du(u), u);
                }
                rtn_flag=gc.create(c);
                if (!rtn_flag)
                {
                  return false;
                }

                // extract the control points from piecewise curve and set the surface control points
                for (u=0; u<nu; ++u)
                {
                  curve_type crv;

                  c.get(crv, u);

                  // resize the temp surface
                  if (j==0)
                  {
                    surfs[u].resize(crv.degree(), max_jdegs[v]);
                  }

                  for (i=0; i<=crv.degree(); ++i)
                  {
                    surfs[u].set_control_point(crv.get_control_point(i), i, j);
                  }
                }
              }

              // put these surfaces into piecewise surface
              typename piecewise_surface_type::error_code ec;
              for (u=0; u<nu; ++u)
              {
                ec=ps.set(surfs[u], u, v);
                if (ec!=piecewise_surface_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
            }

            return true;
          }
#endif

        private:
          piecewise_surface_type orig_surface;
          data_type split_param;
          edge_split_identifier edge_to_split;
      };
    }
  }
}
#endif
