/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - piecewise general surface creator
*    Rob McDonald - Modified to piecewise cubic spline surface creator
********************************************************************************/

#ifndef eli_geom_surface_piecewise_cubic_spline_skinning_surface_creator_hpp
#define eli_geom_surface_piecewise_cubic_spline_skinning_surface_creator_hpp

#include <list>
#include <vector>
#include <iterator>

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_cubic_spline_creator.hpp"

#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/bezier.hpp"
#include "eli/geom/surface/piecewise_connection_data.hpp"
#include "eli/geom/surface/piecewise_creator_base.hpp"

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_cubic_spline_skinning_surface_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef typename base_class_type::piecewise_surface_type piecewise_surface_type;

          typedef connection_data<data_type, dim__, tolerance_type> rib_data_type;

          piecewise_cubic_spline_skinning_surface_creator()
            : piecewise_creator_base<data__, dim__, tol__>(0, 0), ribs(2),
              max_degree(1), closed(false)
          {
          }
          piecewise_cubic_spline_skinning_surface_creator(const data_type &uu0, const data_type &vv0)
            : piecewise_creator_base<data__, dim__, tol__>(uu0, vv0), ribs(2),
              max_degree(1), closed(false)
          {
          }
          piecewise_cubic_spline_skinning_surface_creator(const piecewise_general_skinning_surface_creator<data_type, dim__, tolerance_type> & gs)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(gs), ribs(gs.ribs),
              max_degree(gs.max_degree), closed(gs.closed)
          {
          }
          virtual ~piecewise_cubic_spline_skinning_surface_creator()
          {
          }

          void set_closed() {closed=true;}
          void set_open() {closed=false;}
          bool is_closed() const {return closed;}
          bool is_open() const {return !closed;}

          void set_u0(const data_type &uu0) {this->set_initial_u(uu0);}

          void set_segment_du(const data_type &duu, const index_type &i)
          {
            this->set_du(duu, i);
          }

          void set_tdisc( const std::vector<data_type> &td )
          {
            tdisc = td;
          }

          bool set_conditions(const std::vector<rib_data_type> &rbs, const std::vector<index_type> &maxd, bool cl=false)
          {
            index_type i, j, nsegs(static_cast<index_type>(maxd.size())), nribs(rbs.size());

            // ensure input vectors are correct size
            if (!cl && (nribs!=(nsegs+1)))
              return false;
            if (cl && (nribs!=nsegs))
              return false;

            // check to make sure have valid end conditions
            if (!cl)
            {
              if (rbs[0].use_left_fp() || rbs[0].use_left_fpp() || rbs[0].get_continuity()!=rib_data_type::C0)
              {
                return false;
              }
              if (rbs[nsegs].use_right_fp() || rbs[nsegs].use_right_fpp() || rbs[nsegs].get_continuity()!=rib_data_type::C0)
              {
                return false;
              }
            }

            // make sure ribs are in valid state
            data_type v_start(rbs[0].get_t0()), v_end(rbs[0].get_tmax());
            tolerance_type tol;
            for (i=0; i<nribs; ++i)
            {
              if (!rbs[i].check_state())
                return false;
              if (!tol.approximately_equal(rbs[i].get_t0(), v_start) || !tol.approximately_equal(rbs[i].get_tmax(), v_end))
                return false;
            }

            // find all unique v-coordinates on joints for each rib
            auto comp = [&tol](const data_type &x1, const data_type &x2)->bool
            {
              return tol.approximately_less_than(x1, x2);
            };

            std::vector<data_type> joints;
            data_type t0(rbs[0].get_t0()), tmax(rbs[0].get_tmax());

            rbs[0].get_joints(std::back_inserter(joints));
            for (i=1; i<nribs; ++i)
            {
              // test to make sure this rib's parameterization matches rest
              if (!tol.approximately_equal(rbs[i].get_t0(), t0) || !tol.approximately_equal(rbs[i].get_tmax(), tmax))
              {
                return false;
              }

              // get the joints on the current rib
              std::vector<data_type> rjoints, jts_out;
              rbs[i].get_joints(std::back_inserter(rjoints));

              // merge these joints with current list of joints
              std::set_union(joints.begin(), joints.end(), rjoints.begin(), rjoints.end(), std::back_inserter(jts_out), comp);
              std::swap(joints, jts_out);
            }

            // record where the joints need to be for create()
            index_type njoints(static_cast<index_type>(joints.size()));

            // set the v-parameterization
            this->set_number_v_segments(njoints-1);
            this->set_initial_v(joints[0]);
            for (j=0; j<(njoints-1); ++j)
            {
              this->set_dv(joints[j+1]-joints[j], j);
            }

            // reset the number of u-segments
            this->set_number_u_segments(nsegs);

            ribs=rbs;
            max_degree=maxd;
            closed=cl;

            return true;
          }

          virtual bool create(piecewise_surface_type &ps) const
          {
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

            std::vector<index_type> idisc( tdisc.size() );
            std::vector<data_type> uvalvec( nu + 1 );
            uvalvec[0] = this->get_u0();
            for ( u = 0; u < nu; ++u)
            {
              uvalvec[u+1] = uvalvec[u] + this->get_segment_du( u );
            }
            uvalvec.back()=tdisc.back();

            u = 0;
            for ( i = 0; i < tdisc.size(); i++ )
            {
              while ( u <= nu && uvalvec[u] < tdisc[i] )
              {
                u++;
              }
              idisc[i] = u;
              u++;
            }


            // build segments based on rib information
            // here should have everything to make an nribs x njoints piecewise surface with all
            // of the j-degrees matching in the u-direction so that can use general curve creator
            // techniques to create control points
            for (v=0; v<nv; ++v)
            {
              typedef eli::geom::curve::piecewise_cubic_spline_creator<data_type, dim__, tolerance_type> piecewise_curve_creator_type;
              typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
              typedef typename piecewise_curve_type::curve_type curve_type;

              std::vector < point_type > pts( nu+1 );
              piecewise_curve_creator_type gc;
              piecewise_curve_type c, cseg;

              std::vector<surface_type> surfs(nu);

              for (j=0; j<=max_jdegs[v]; ++j)
              {
                // cycle through each rib to set corresponding joint info
                for (u=0; u<=nu; ++u)
                {
                  curve_type jcrv;

                  rib_states[u].get_f().get(jcrv, v);

                  pts[u] = jcrv.get_control_point(j);
                }

                c.clear();
                c.set_t0(this->get_u0());

                for ( i = 0; i < idisc.size()-1; i++ )
                {
                  int nseg( idisc[i+1] - idisc[i] );

                  std::vector < point_type > subpts( nseg+1 );

                  gc.set_number_segments( nseg );
                  gc.set_t0( tdisc[i] );

                  // set the parameterizations and create curve
                  index_type iseg = 0;
                  for ( u = idisc[i]; u < idisc[i+1]; ++u)
                  {
                    gc.set_segment_dt( this->get_segment_du( u ), iseg);
                    subpts[iseg] = pts[u];
                    iseg++;
                  }
                  subpts[iseg] = pts[u];


                  gc.set_chip( subpts.begin(), eli::geom::general::NOT_CONNECTED );

                  bool rtn_flag=gc.create(cseg);
                  if (!rtn_flag)
                  {
                    return false;
                  }

                  c.push_back( cseg );
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

        private:
          std::vector<rib_data_type> ribs;
          std::vector<index_type> max_degree;
          std::vector<data_type> tdisc;
          bool closed;
      };
    }
  }
}
#endif
