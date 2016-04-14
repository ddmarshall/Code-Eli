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

#ifndef eli_geom_surface_piecewise_capped_surface_creator_hpp
#define eli_geom_surface_piecewise_capped_surface_creator_hpp

#include <list>
#include <iterator>

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise.hpp"

#include "eli/geom/surface/bezier.hpp"
#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/piecewise_connection_data.hpp"
#include "eli/geom/surface/piecewise_general_skinning_surface_creator.hpp"

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_capped_surface_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          enum edge_cap_identifier
          {
            CAP_NONE,
            CAP_UMIN,
            CAP_UMAX,
            CAP_VMIN,
            CAP_VMAX
          };

          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef typename base_class_type::piecewise_surface_type piecewise_surface_type;

          piecewise_capped_surface_creator()
            : piecewise_creator_base<data__, dim__, tol__>(0, 0), delta_param(1), edge_to_cap(CAP_NONE)
          {
          }
          piecewise_capped_surface_creator(const piecewise_surface_type &os, const data_type &dp, edge_cap_identifier ec)
            : piecewise_creator_base<data__, dim__, tol__>(0, 0), orig_surface(os), delta_param(dp), edge_to_cap(ec)
          {
          }
          piecewise_capped_surface_creator(const piecewise_capped_surface_creator<data_type, dim__, tolerance_type> & gs)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(gs), orig_surface(gs.orig_surface),
              delta_param(gs.delta_param), edge_to_cap(gs.edge_to_cap)
          {
          }
          virtual ~piecewise_capped_surface_creator()
          {
          }

          bool set_conditions(const piecewise_surface_type &os, const data_type &dp, edge_cap_identifier ec)
          {
            // all deltas are positive, even on the min edges
            if (dp<=0)
              return false;

            orig_surface = os;
            edge_to_cap = ec;
            delta_param = dp;

            return true;
          }

          virtual bool create(piecewise_surface_type &ps) const
          {
            typedef typename eli::geom::curve::piecewise<eli::geom::curve::bezier, data_type, dim__, tolerance_type> piecewise_curve_type;

            tolerance_type tol;

            // extract the curve corresonding to the edge of surface wanting to cap
            piecewise_curve_type edge;
            switch (edge_to_cap)
            {
              case (CAP_UMIN):
              {
                orig_surface.get_umin_bndy_curve( edge );
                break;
              }
              case (CAP_UMAX):
              {
                orig_surface.get_umax_bndy_curve( edge );
                break;
              }
              case (CAP_VMIN):
              {
                orig_surface.get_vmin_bndy_curve( edge );
                break;
              }
              case (CAP_VMAX):
              {
                orig_surface.get_vmax_bndy_curve( edge );
                break;
              }
              default:
              {
                // must not want any edge capped
                assert(edge_to_cap==CAP_NONE);
                return false;
              }
            }

            // make sure extracted curve is closed
            if (edge.open())
            {
              return false;
            }

            // make sure that the split location is different than the start/end location otherwise
            // the edge is a point and no need to cap
            data_type tmin(edge.get_t0()), tmax(edge.get_tmax()), tsplit((tmin+tmax)/2);
            if (tol.approximately_equal(edge.f(tmin), edge.f(tsplit)))
            {
              return false;
            }

            // split the curve at the mid-parameter location and reverse the second piece's parameterization
            piecewise_curve_type first_half, second_half;
            edge.split(first_half, second_half, tsplit);
            second_half.reverse();
            second_half.set_t0(first_half.get_t0());

            // assert(first_half.get_t0()==second_half.get_t0());
            // assert(first_half.get_tmax()==second_half.get_tmax());

            // create surface connecting two edges with param spacing of 2*delta_param and points for other two edges
            // cap u-direction goes from second_half curve to first_half curve
            // cap v-direction follows first_half direction
            piecewise_surface_type cap;
            {
              typedef typename eli::geom::surface::connection_data<data__, 3, tolerance_type> rib_data_type;
              typedef typename eli::geom::surface::piecewise_general_skinning_surface_creator<data__, 3, tolerance_type> general_creator_type;

              std::vector<rib_data_type> ribs(2);
              std::vector<typename general_creator_type::index_type> max_degree(1);
              general_creator_type gc;
              bool rtn_flag;

              // set the rib data
              ribs[0].set_f(second_half);
              ribs[1].set_f(first_half);

              // set the maximum degrees of each segment
              max_degree[0]=0;

              // create surface
              rtn_flag=gc.set_conditions(ribs, max_degree, false);
              if (!rtn_flag)
              {
                assert(false);
                return false;
              }
              gc.set_u0(orig_surface.get_u0()-2*delta_param);
              gc.set_segment_du(2*delta_param, 0);
              rtn_flag=gc.create(cap);
              if (!rtn_flag)
              {
                assert(false);
                return false;
              }
            }

            // split the cap surface at mid point
            cap.split_u(orig_surface.get_u0()-delta_param);

            // resize output surface to needed size
            data_type u0, v0;
            std::vector<data_type> ucap_param, vcap_param, uparam, vparam, du, dv;
            index_type i, j, icap_mid, umin_offset(0), vmin_offset(0);

            orig_surface.get_pmap_u(uparam);
            cap.get_pmap_uv(ucap_param, vcap_param);
            // NOTE: This assumes that there are same number of patches on both sides of split
            icap_mid=ucap_param.size()/2;
            switch (edge_to_cap)
            {
              case (CAP_UMIN):
              {
                // establish the new u and v parameterization of surface
                const index_type nucap_patch(icap_mid), nvcap_patch(cap.number_v_patches());
                u0=ucap_param[icap_mid];
                du.resize(uparam.size()-1+nucap_patch);
                for (i=icap_mid; i<static_cast<index_type>(ucap_param.size())-1; ++i)
                {
                  du[i-icap_mid]=ucap_param[i+1]-ucap_param[i];
                }
                for (i=0; i<static_cast<index_type>(uparam.size())-1; ++i)
                {
                  du[i+nucap_patch]=uparam[i+1]-uparam[i];
                }

                // set the new v-parameterization
                vparam.resize(2*vcap_param.size()-1);
                data_type vmid(vcap_param[vcap_param.size()-1]);
                for (j=0; j<static_cast<index_type>(vcap_param.size()); ++j)
                {
                  vparam[j]=vcap_param[j];
                  vparam[vparam.size()-1-j]=2*vmid-vparam[j];
                }

                // set the v-parameters
                v0=vparam[0];
                dv.resize(vparam.size()-1);
                for (j=0; j<static_cast<index_type>(dv.size()); ++j)
                {
                  dv[j]=vparam[j+1]-vparam[j];
                }

                // split a copy original surface at any new v-coordinate locations
                piecewise_surface_type orig_copy(orig_surface);
                typename piecewise_surface_type::surface_type patch;

                for (j=0; j<static_cast<index_type>(vparam.size()); ++j)
                {
                  orig_copy.split_v(vparam[j]);
                }
                assert(orig_copy.number_v_patches()==(static_cast<index_type>(vparam.size())-1));

                // resize the output surface
                ps.init_uv(du.begin(), du.end(), dv.begin(), dv.end(), u0, v0);
                umin_offset=nucap_patch;

                // add first half of cap surfaces
                for (i=0; i<nucap_patch; ++i)
                {
                  for (j=0; j<nvcap_patch; ++j)
                  {
                    cap.get(patch, i + nucap_patch, j);
                    ps.set(patch, i, j);
                  }
                }

                // reverse the cap so that the second half
                cap.reverse_u();
                cap.reverse_v();

                // add second half of cap surfaces
                for (i=0; i<nucap_patch; ++i)
                {
                  for (j=nvcap_patch; j<static_cast<index_type>(dv.size()); ++j)
                  {
                    cap.get(patch, i + nucap_patch, j-nvcap_patch);
                    ps.set(patch, i, j);
                  }
                }

                // add non-cap surfaces
                for (i=0; i<orig_copy.number_u_patches(); ++i)
                {
                  for (j=0; j<orig_copy.number_v_patches(); ++j)
                  {
                    orig_copy.get(patch, i, j);
                    ps.set(patch, i+umin_offset, j+vmin_offset);
                  }
                }

                break;
              }
              case (CAP_UMAX):
              {
                // need to reverse u-direction
                cap.reverse_u();

                // establish the new u and v parameterization of surface
                const index_type nucap_patch(icap_mid), nu_patch(uparam.size()-1), nvcap_patch(cap.number_v_patches());
                u0=uparam[0];
                du.resize(uparam.size()-1+nucap_patch);
                for (i=0; i<static_cast<index_type>(uparam.size())-1; ++i)
                {
                  du[i]=uparam[i+1]-uparam[i];
                }
                for (i=nu_patch; i<static_cast<index_type>(du.size()); ++i)
                {
                  du[i]=ucap_param[(i-nu_patch)+1]-ucap_param[i-nu_patch];
                }

                // set the new v-parameterization
                vparam.resize(2*vcap_param.size()-1);
                data_type vmid(vcap_param[vcap_param.size()-1]);
                for (j=0; j<static_cast<index_type>(vcap_param.size()); ++j)
                {
                  vparam[j]=vcap_param[j];
                  vparam[vparam.size()-1-j]=2*vmid-vparam[j];
                }

                // set the v-parameters
                v0=vparam[0];
                dv.resize(vparam.size()-1);
                for (j=0; j<static_cast<index_type>(dv.size()); ++j)
                {
                  dv[j]=vparam[j+1]-vparam[j];
                }

                // split a copy original surface at any new v-coordinate locations
                piecewise_surface_type orig_copy(orig_surface);
                typename piecewise_surface_type::surface_type patch;

                for (j=0; j<static_cast<index_type>(vparam.size()); ++j)
                {
                  orig_copy.split_v(vparam[j]);
                }
                assert(orig_copy.number_v_patches()==(static_cast<index_type>(vparam.size())-1));

                // resize the output surface
                ps.init_uv(du.begin(), du.end(), dv.begin(), dv.end(), u0, v0);

                // add first half of cap surfaces
                for (i=0; i<nucap_patch; ++i)
                {
                  for (j=0; j<nvcap_patch; ++j)
                  {
                    cap.get(patch, i, j);
                    ps.set(patch, i + nu_patch, j);
                  }
                }

                // reverse the cap so that the second half
                cap.reverse_u();
                cap.reverse_v();

                // add second half of cap surfaces
                for (i=0; i<nucap_patch; ++i)
                {
                  for (j=nvcap_patch; j<static_cast<index_type>(dv.size()); ++j)
                  {
                    cap.get(patch, i, j-nvcap_patch);
                    ps.set(patch, i + nu_patch, j);
                  }
                }

                // add non-cap surfaces
                for (i=0; i<orig_copy.number_u_patches(); ++i)
                {
                  for (j=0; j<orig_copy.number_v_patches(); ++j)
                  {
                    orig_copy.get(patch, i, j);
                    ps.set(patch, i+umin_offset, j+vmin_offset);
                  }
                }

                break;
              }
              case (CAP_VMIN):
              {
                assert(false);
                return false;
                break;
              }
              case (CAP_VMAX):
              {
                assert(false);
                return false;
                break;
              }
              default:
              {
                // must not want any edge capped
                assert(edge_to_cap==CAP_NONE);
                return false;
              }
            }

            return true;
          }

        private:
          piecewise_surface_type orig_surface;
          data_type delta_param;
          edge_cap_identifier edge_to_cap;
      };
    }
  }
}
#endif
