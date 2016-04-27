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

#ifndef eli_geom_surface_piecewise_hpp
#define eli_geom_surface_piecewise_hpp

#include <vector>
#include <iterator>
#include <utility>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"
#include "eli/geom/curve/piecewise.hpp"

#include "eli/geom/general/continuity.hpp"
#include "eli/geom/curve/equivalent_curves.hpp"

namespace eli
{
  namespace geom
  {

    namespace surface
    {
      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__ >
      class piecewise;
    }

    namespace intersect
    {
      template<template<typename, unsigned short, typename> class surface1__, typename data1__, unsigned short dim1__, typename tol1__ >
      typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::data_type
        minimum_distance(
          typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::data_type &u,
	        typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::data_type &v,
          const surface::piecewise<surface1__, data1__, dim1__, tol1__> &ps,
          const typename surface::piecewise<surface1__, data1__, dim1__, tol1__>::point_type &pt);
    }

    namespace surface
    {
      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__=eli::util::tolerance<data__> >
      class piecewise
      {
        public:
          typedef surface__<data__, dim__, tol__> surface_type;
          typedef typename surface_type::index_type index_type;
          typedef typename surface_type::point_type point_type;
          typedef typename surface_type::control_point_type control_point_type;
          typedef typename surface_type::rotation_matrix_type rotation_matrix_type;
          typedef typename surface_type::bounding_box_type bounding_box_type;
          typedef typename surface_type::curve_type curve_type;
          typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, dim__, tol__> piecewise_curve_type;

          typedef std::pair< index_type, index_type > patch_boundary_code_type;

          typedef data__ data_type;
          typedef unsigned short dimension_type;
          typedef tol__ tolerance_type;
          enum error_code
          {
            NO_ERRORS=0,
            INVALID_INDEX=1,
            INDEX_NOT_FOUND=2,
            INVALID_PARAM=50,
            INVALID_PARAM_DIFFERENCE=51,
            PATCH_NOT_CONNECTED=100,
            UNKNOWN_ERROR=999
          };

        public:
          piecewise() : nu(0), nv(0) {}
          piecewise(const piecewise<surface__, data_type, dim__, tol__> &p)
            : patches(p.patches), ukey(p.ukey), vkey(p.vkey), nu(p.nu), nv(p.nv) {}
          ~piecewise() {}

          piecewise & operator=(const piecewise<surface__, data_type, dim__> &p)
          {
            if (this==&p)
              return (*this);

            patches=p.patches;
            ukey=p.ukey;
            vkey=p.vkey;
            nu=p.nu;
            nv=p.nv;

            return (*this);
          }

          bool operator==(const piecewise<surface__, data_type, dim__> &p) const
          {
            if (this==&p)
              return true;
            if (nu!=p.nu)
              return false;
            if (nv!=p.nv)
              return false;
            if (ukey!=p.ukey)
              return false;
            if (vkey!=p.vkey)
              return false;
            if (number_u_patches()!=p.number_u_patches())
              return false;
            if (number_v_patches()!=p.number_v_patches())
              return false;
            typename patch_collection_type::const_iterator scit, it;
            for (scit=patches.begin(), it=p.patches.begin(); scit!=patches.end(); ++scit, ++it)
            {
              if ((*it)!=(*scit))
                return false;
            }

            return true;
          }

          bool operator!=(const piecewise<surface__, data_type, dim__> &p) const
          {
            return !operator==(p);
          }

          static dimension_type dimension() {return dim__;}

          data_type get_u0() const {return ukey.get_pmin();}
          void set_u0(const data_type &u0_in) {ukey.set_pmin(u0_in);}

          data_type get_v0() const {return vkey.get_pmin();}
          void set_v0(const data_type &v0_in) {vkey.set_pmin(v0_in);}

          data_type get_umax() const {return ukey.get_pmax();}
          data_type get_vmax() const {return vkey.get_pmax();}

          index_type number_u_patches() const {return nu;}
          index_type number_v_patches() const {return nv;}

          surface_type * get_patch( const index_type &ui, const index_type &vi)
          {
              index_type uk, vk;
              find_patch(uk, vk, ui, vi);
              return &patches[uk][vk];
          }

          const surface_type * get_patch( const index_type &ui, const index_type &vi) const
          {
              index_type uk, vk;
              find_patch(uk, vk, ui, vi);
              return &patches[uk][vk];
          }

          surface_type * get_patch( const index_type &ui, const index_type &vi, double &ustart, double &du, double &vstart, double &dv)
          {
              index_type uk, vk;
              typename keymap_type::const_iterator uit, vit;

              find_patch( uk, vk, uit, vit, ui, vi );

              ustart = uit->first;
              du = ukey.get_delta_parm( uit );

              vstart = vit->first;
              dv = vkey.get_delta_parm( vit );

              return &patches[uk][vk];
          }

          const surface_type * get_patch( const index_type &ui, const index_type &vi, double &ustart, double &du, double &vstart, double &dv) const
          {
              index_type uk, vk;
              typename keymap_type::const_iterator uit, vit;

              find_patch( uk, vk, uit, vit, ui, vi );

              ustart = uit->first;
              du = ukey.get_delta_parm( uit );

              vstart = vit->first;
              dv = vkey.get_delta_parm( vit );

              return &patches[uk][vk];
          }

          surface_type * get_patch_unordered( const index_type &uk, const index_type &vk)
          {
              return &patches[uk][vk];
          }

          const surface_type * get_patch_unordered( const index_type &uk, const index_type &vk) const
          {
              return &patches[uk][vk];
          }

          void get_parameter_min(data_type &umin, data_type &vmin) const
          {
            umin=ukey.get_pmin();
            vmin=vkey.get_pmin();
          }

          void get_parameter_max(data_type &umax, data_type &vmax) const
          {
            umax=ukey.get_pmax();
            vmax=vkey.get_pmax();
          }

          void parameter_report() const
          {
            printf("U parameter:\n");
            ukey.parameter_report();
            printf("V parameter:\n");
            vkey.parameter_report();
          }

          void octave_print(int figno ) const
          {
            index_type i, j, pp, qq, nup, nvp;
            data_type umin, vmin, umax, vmax;

            nup = number_u_patches();
            nvp = number_v_patches();
            get_parameter_min(umin, vmin);
            get_parameter_max(umax, vmax);

            std::cout << "figure(" << figno << ");" << std::endl;

            // initialize the u & v parameters
            std::vector<data__> u(31), v(31);
            for (i=0; i<static_cast<index_type>(u.size()); ++i)
            {
              u[i]=umin+(umax-umin)*static_cast<data__>(i)/(u.size()-1);
            }
            for (j=0; j<static_cast<index_type>(v.size()); ++j)
            {
              v[j]=vmin+(vmax-vmin)*static_cast<data__>(j)/(v.size()-1);
            }

            // set the surface points
            std::cout << "surf_x=[";
            for (i=0; i<static_cast<index_type>(u.size()); ++i)
            {
              std::cout << this->f(u[i], v[0]).x();
              for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
              {
                std::cout << ", " << this->f(u[i], v[j]).x();
              }
              j=static_cast<index_type>(v.size()-1);
              std::cout << ", " << this->f(u[i], v[j]).x();
              if (i<static_cast<index_type>(u.size()-1))
                std::cout << "; " << std::endl;
            }
            std::cout << "];" << std::endl;

            std::cout << "surf_y=[";
            for (i=0; i<static_cast<index_type>(u.size()); ++i)
            {
              std::cout << f(u[i], v[0]).y();
              for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
              {
                std::cout << ", " << f(u[i], v[j]).y();
              }
              j=static_cast<index_type>(v.size()-1);
              std::cout << ", " << f(u[i], v[j]).y();
              if (i<static_cast<index_type>(u.size()-1))
                std::cout << "; " << std::endl;
            }
            std::cout << "];" << std::endl;

            std::cout << "surf_z=[";
            for (i=0; i<static_cast<index_type>(u.size()); ++i)
            {
              std::cout << f(u[i], v[0]).z();
              for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
              {
                std::cout << ", " << f(u[i], v[j]).z();
              }
              j=static_cast<index_type>(v.size()-1);
              std::cout << ", " << f(u[i], v[j]).z();
              if (i<static_cast<index_type>(u.size()-1))
                std::cout << "; " << std::endl;
            }
            std::cout << "];" << std::endl;

            std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
            std::cout << "mesh(surf_x, surf_y, surf_z, zeros(size(surf_z)), 'EdgeColor', [0 0 0]);" << std::endl;
            std::cout << "axis equal" << std::endl;
            std::cout << "axis off" << std::endl;
          }

          void get_pmap_u(std::vector<data_type> &pmap) const
          {
            ukey.get_pmap(pmap);
          }

          void get_pmap_v(std::vector<data_type> &pmap) const
          {
            vkey.get_pmap(pmap);
          }

          void get_pmap_uv(std::vector<data_type> &upmap, std::vector<data_type> &vpmap) const
          {
            ukey.get_pmap(upmap);
            vkey.get_pmap(vpmap);
          }

          void init_u(const index_type &nsegu, const data_type &du = 1, const data_type &u0 = 0)
          {
            patches.clear();
            resize_store(nsegu, nv);
            ukey.init(nsegu, du, u0);
          }

          void init_v(const index_type &nsegv, const data_type &dv = 1, const data_type &v0 = 0)
          {
            patches.clear();
            resize_store(nu, nsegv);
            vkey.init(nsegv, dv, v0);
          }

          void init_uv(const index_type &nsegu, const index_type &nsegv, const data_type &du = 1, const data_type &dv = 1, const data_type &u0 = 0, const data_type &v0 = 0)
          {
            patches.clear();
            resize_store(nsegu, nsegv);
            ukey.init(nsegu, du, u0);
            vkey.init(nsegv, dv, v0);
          }

          template<typename it__>
          void init_u(const it__ &dus, const it__ &due, const data_type &u0 = 0)
          {
            patches.clear();
            ukey.init(dus, due, u0);
            resize_store(ukey.key.size(), nv);
          }

          template<typename it__>
          void init_v(const it__ &dvs, const it__ &dve, const data_type &v0 = 0)
          {
            patches.clear();
            vkey.init(dvs, dve, v0);
            resize_store(nu, vkey.key.size());
          }

          template<typename it__>
          void init_uv(const it__ &dus, const it__ &due, const it__ &dvs, const it__ &dve, const data_type &u0 = 0, const data_type &v0 = 0)
          {
            patches.clear();
            ukey.init(dus, due, u0);
            vkey.init(dvs, dve, v0);
            resize_store(ukey.key.size(), vkey.key.size());
          }

          void degree_u(index_type &mind, index_type &maxd)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            uit=patches.begin();
            vit=(*uit).begin();

            index_type d = vit->degree_u();
            mind = d;
            maxd = d;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                d = vit->degree_u();

                if(d<mind)
                {
                  mind = d;
                }
                if(d>maxd)
                {
                  maxd=d;
                }
              }
            }
          }

          void degree_u( std::vector<index_type> &deg)
          {
            index_type vk, j;
            typename keymap_type::const_iterator vit;

            deg.resize( nv );

            for ( j = 0, vit = vkey.key.begin(); vit != vkey.key.end(); ++vit, ++j )
            {
              vk = vit->second;

              deg[j] = patches[0][vk].degree_u();
            }
          }

          void degree_v(index_type &mind, index_type &maxd)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            uit=patches.begin();
            vit=(*uit).begin();

            index_type d = vit->degree_v();
            mind = d;
            maxd = d;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                d = vit->degree_v();

                if(d<mind)
                {
                  mind = d;
                }
                if(d>maxd)
                {
                  maxd=d;
                }
              }
            }
          }

          void degree_v( std::vector<index_type> &deg)
          {
            index_type uk, i;
            typename keymap_type::const_iterator uit;

            deg.resize( nu );

            for ( i = 0, uit = ukey.key.begin(); uit != ukey.key.end(); ++uit, ++i )
            {
              uk = uit->second;

              deg[i] = patches[uk][0].degree_v();
            }
          }

          bool open_u() const
          {
            return !closed_u();
          }
          bool closed_u() const
          {
            index_type ifirst, ilast, j;
            typename surface_type::curve_type bc0, bc1;

            ifirst = ukey.key.begin()->second;
            ilast = ukey.key.rbegin()->second;

            for (j=0; j<nv; ++j)
            {
              patches[ifirst][j].get_uconst_curve(bc0, 0);
              patches[ilast][j].get_uconst_curve(bc1, 1);
              if (!eli::geom::curve::equivalent_curves(bc0, bc1))
                return false;
            }

            return true;
          }

          bool open_v() const
          {
            return !closed_v();
          }
          bool closed_v() const
          {
            index_type i, jfirst, jlast;
            typename surface_type::curve_type bc0, bc1;

            jfirst = vkey.key.begin()->second;
            jlast = vkey.key.rbegin()->second;

            for (i=0; i<nu; ++i)
            {
              patches[i][jfirst].get_vconst_curve(bc0, 0);
              patches[i][jlast].get_vconst_curve(bc1, 1);
              if (!eli::geom::curve::equivalent_curves(bc0, bc1))
              {
                return false;
              }
            }

            return true;
          }

          void get_bounding_box(bounding_box_type &bb) const
          {
            typename patch_collection_type::const_iterator uit;
            typename patch_strip_type::const_iterator vit;
            bounding_box_type bb_local;

            bb.clear();

            // cycle through all patches to get each bounding box to compare
            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->get_bounding_box(bb_local);
                bb.add(bb_local);
              }
            }
          }

          void rotate(const rotation_matrix_type &rmat)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->rotate(rmat);
              }
            }
          }

          void rotate(const rotation_matrix_type &rmat, const point_type &rorig)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->rotate(rmat, rorig);
              }
            }
          }

          void translate(const point_type &trans)
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->translate(trans);
              }
            }
          }

          void reverse_u()
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->reverse_u();
              }
            }
            ukey.reverse_keymap();
          }

          void reverse_v()
          {
            typename patch_collection_type::iterator uit;
            typename patch_strip_type::iterator vit;

            for (uit=patches.begin(); uit!=patches.end(); ++uit)
            {
              for (vit=(*uit).begin(); vit!=(*uit).end(); ++vit)
              {
                vit->reverse_v();
              }
            }
            vkey.reverse_keymap();
          }

          void swap_uv()
          {
            patch_collection_type old_patches;
            old_patches.swap(patches);

            index_type nu_old(nu), nv_old(nv);

            // Resizes patches and also assigns nu, nv.
            resize_store(nv_old, nu_old);

            for (index_type i=0; i<nu; ++i)
            {
              for (index_type j=0; j<nv; ++j)
              {
                patches[i][j]=old_patches[j][i];
                patches[i][j].swap_uv();
              }
            }

            data_type pmaxtmp;
            pmaxtmp = ukey.pmax;
            ukey.pmax = vkey.pmax;
            vkey.pmax = pmaxtmp;

            swap(ukey.key, vkey.key);
          }

          void clear()
          {
            nu=0;
            nv=0;
            patches.clear();
            ukey.clear();
            vkey.clear();
          }

          error_code get(surface_type &surf, const index_type &ui, const index_type &vi) const
          {
            data_type du, dv;
            return get(surf, du, dv, ui, vi);
          }

          error_code get(surface_type &surf, data_type &du, data_type &dv, const index_type &ui, const index_type &vi) const
          {
            if ((ui>=number_u_patches()) || (vi>=number_v_patches()))
              return INVALID_INDEX;

            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            find_patch(uk, vk, uit, vit, ui, vi);

            du = ukey.get_delta_parm(uit);
            dv = vkey.get_delta_parm(vit);
            surf = patches[uk][vk];

            return NO_ERRORS;
          }

          error_code set(const surface_type &surf, const index_type &ui, const index_type &vi)
          {
            if ((ui>=number_u_patches()) || (vi>=number_v_patches()))
              return INVALID_INDEX;

            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            find_patch(uk, vk, uit, vit, ui, vi);

            // set the new surf
            patches[uk][vk]=surf;

            return NO_ERRORS;
          }

          error_code replace(const surface_type &surf, const index_type &ui, const index_type &vi)
          {
            if ((ui>=number_u_patches()) || (vi>=number_v_patches()))
              return INVALID_INDEX;

            // advance to desired index
            typename surface_type::curve_type bc0, bc1;
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            find_patch(uk, vk, uit, vit, ui, vi);

            surface_type s = patches[uk][vk];

            if (ui>0)
            {
              surf.get_uconst_curve(bc0, 0);
              s.get_uconst_curve(bc1, 0);
              if (!eli::geom::curve::equivalent_curves(bc0, bc1))
              {
                return PATCH_NOT_CONNECTED;
              }
            }
            if ((ui+1)<number_u_patches())
            {
              surf.get_uconst_curve(bc0, 1);
              s.get_uconst_curve(bc1, 1);
              if (!eli::geom::curve::equivalent_curves(bc0, bc1))
              {
                return PATCH_NOT_CONNECTED;
              }
            }
            if (vi>0)
            {
              surf.get_vconst_curve(bc0, 0);
              s.get_vconst_curve(bc1, 0);
              if (!eli::geom::curve::equivalent_curves(bc0, bc1))
              {
                return PATCH_NOT_CONNECTED;
              }
            }
            if ((vi+1)<number_v_patches())
            {
              surf.get_vconst_curve(bc0, 1);
              s.get_vconst_curve(bc1, 1);
              if (!eli::geom::curve::equivalent_curves(bc0, bc1))
              {
                return PATCH_NOT_CONNECTED;
              }
            }

            // set the new surf
            patches[uk][vk]=surf;

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERRORS;
          }

          error_code split_u(const data_type &u_in)
          {
            index_type uk, vk;
            typename keymap_type::iterator uit, vit;
            data_type uu(0), vv(0);
            data_type vmin = vkey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, u_in, vmin);

            // check for out of range input
            if ((uk == -1) || (vk == -1))
              return INVALID_PARAM;

            // check if no need to split
            tolerance_type tol;
            if (tol.approximately_equal(uu, 0))
              return NO_ERRORS;
            if (tol.approximately_equal(uu, 1))
              return NO_ERRORS;

            return split_u(uk, uit, u_in, uu);
          }

          error_code split_u(piecewise<surface__, data_type, dim__, tol__> &before, piecewise<surface__, data_type, dim__, tol__> &after, const data_type &u_in) const
          {
            before.clear();
            after.clear();

            if (u_in < ukey.get_pmin())
            {
              after=(*this);
              return NO_ERRORS;
            }

            if (u_in > ukey.get_pmax())
            {
              before=(*this);
              return NO_ERRORS;
            }

            piecewise<surface__, data_type, dim__, tol__> s(*this);

            error_code spliterr = s.split_u( u_in );

            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type vmin = vkey.get_pmin();

            s.find_patch(uit, vit, uu, vv, u_in, vmin);

            s.subsurf(before, s.ukey.key.begin(), uit, s.vkey.key.begin(), s.vkey.key.end());
            s.subsurf(after, uit, s.ukey.key.end(), s.vkey.key.begin(), s.vkey.key.end());

            return spliterr;
          }

          error_code split_v(const data_type &v_in)
          {
            index_type uk, vk;
            typename keymap_type::iterator uit, vit;
            data_type uu(0), vv(0);
            data_type umin = ukey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, umin, v_in);

            // check for out of range input
            if ((uk == -1) || (vk == -1))
              return INVALID_PARAM;

            // check if no need to split
            tolerance_type tol;
            if (tol.approximately_equal(vv, 0))
              return NO_ERRORS;
            if (tol.approximately_equal(vv, 1))
              return NO_ERRORS;

            return split_v(vk, vit, v_in, vv);
          }

          error_code split_v(piecewise<surface__, data_type, dim__, tol__> &before, piecewise<surface__, data_type, dim__, tol__> &after, const data_type &v_in) const
          {
            before.clear();
            after.clear();

            if (v_in < vkey.get_pmin())
            {
              after=(*this);
              return NO_ERRORS;
            }

            if (v_in > vkey.get_pmax())
            {
              before=(*this);
              return NO_ERRORS;
            }

            piecewise<surface__, data_type, dim__, tol__> s(*this);

            error_code spliterr = s.split_v( v_in );

            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type umin = ukey.get_pmin();

            s.find_patch(uit, vit, uu, vv, umin, v_in);

            s.subsurf(before, s.ukey.key.begin(), s.ukey.key.end(), s.vkey.key.begin(), vit);
            s.subsurf(after, s.ukey.key.begin(), s.ukey.key.end(), vit, s.vkey.key.end());

            return spliterr;
          }

          void to_cubic_u(const data_type &ttol)
          {
            typename keymap_type::iterator uit, vit;

            // First pass to split patches until cubic approximation is within tolerance.
            for(uit = ukey.key.begin(); uit != ukey.key.end(); ++uit)
            {
              for(vit = vkey.key.begin(); vit != vkey.key.end(); ++vit)
              {
                index_type uk = uit->second;
                index_type vk = vit->second;

                surface_type s = patches[uk][vk];
                surface_type sc(s);

                sc.to_cubic_u();

                data_type d = s.eqp_distance_bound(sc);

                while(d > ttol)
                {
                  data_type delta_u = ukey.get_delta_parm(uit);
                  data_type u_in = uit->first + static_cast<data_type>(0.5) * delta_u;

                  split_u(uk, uit, u_in, 0.5);

                  s = patches[uk][vk];
                  sc = s;

                  sc.to_cubic_u();

                  d = s.eqp_distance_bound(sc);
                }
              }
            }

            // Second pass to convert all patches to cubic.
            for (index_type uk=0; uk<nu; ++uk)
            {
              for (index_type vk=0; vk<nv; ++vk)
              {
                patches[uk][vk].to_cubic_u();
              }
            }
          }

          void to_cubic_v(const data_type &ttol)
          {
            typename keymap_type::iterator uit, vit;

            // First pass to split patches until cubic approximation is within tolerance.
            for(uit = ukey.key.begin(); uit != ukey.key.end(); ++uit)
            {
              for(vit = vkey.key.begin(); vit != vkey.key.end(); ++vit)
              {
                index_type uk = uit->second;
                index_type vk = vit->second;

                surface_type s = patches[uk][vk];
                surface_type sc(s);

                sc.to_cubic_v();

                data_type d = s.eqp_distance_bound(sc);

                while(d > ttol)
                {
                  data_type delta_v = vkey.get_delta_parm(vit);
                  data_type v_in = vit->first + static_cast<data_type>(0.5) * delta_v;

                  split_v(vk, vit, v_in, 0.5);

                  s = patches[uk][vk];
                  sc = s;

                  sc.to_cubic_v();

                  d = s.eqp_distance_bound(sc);
                }
              }
            }

            // Second pass to convert all patches to cubic.
            for (index_type uk=0; uk<nu; ++uk)
            {
              for (index_type vk=0; vk<nv; ++vk)
              {
                patches[uk][vk].to_cubic_v();
              }
            }
          }

          void to_cubic(const data_type &ttol)
          {
            to_cubic_u(ttol);
            to_cubic_v(ttol);
          }

          void promote_u_to( const std::vector< index_type > ord )
          {
            assert ( ord.size() == nu );

            index_type uk;
            typename keymap_type::const_iterator uit;

            index_type i;

            for ( i = 0, uit = ukey.key.begin(); uit != ukey.key.end(); ++uit, ++i )
            {
              uk = uit->second;

              for (index_type j=0; j<nv; ++j)
              {
                surface_type s=patches[uk][j];

                s.promote_u_to( ord[i] );
              }
            }
          }

          void promote_v_to( const std::vector< index_type > ord )
          {
            assert ( ord.size() == nv );

            index_type vk;
            typename keymap_type::const_iterator vit;

            index_type j;

            for ( j = 0, vit = vkey.key.begin(); vit != vkey.key.end(); ++vit, ++j )
            {
              vk = vit->second;

              for (index_type i=0; i<nu; ++i)
              {
                surface_type s=patches[i][vk];

                s.promote_v_to( ord[j] );
              }
            }
          }

          void get_uconst_curve(piecewise_curve_type &pwc, const data_type &u) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type vmin = vkey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, u, vmin);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_uconst_curve(c, uu);

              pwc.push_back(c,dv);
            }
          }

          void get_umin_bndy_curve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator vit;
            data_type vmin = vkey.get_pmin();

            uk = ukey.key.begin()->second;

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_umin_bndy_curve(c);

              pwc.push_back(c,dv);
            }
          }

          void get_umax_bndy_curve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator vit;
            data_type vmin = vkey.get_pmin();

            uk = ukey.key.rbegin()->second;

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_umax_bndy_curve(c);

              pwc.push_back(c,dv);
            }
          }

          void get_vmin_bndy_curve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit;
            data_type umin = ukey.get_pmin();

            vk = vkey.key.begin()->second;

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_vmin_bndy_curve(c);

              pwc.push_back(c,du);
            }
          }

          void get_vmax_bndy_curve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit;
            data_type umin = ukey.get_pmin();

            vk = vkey.key.rbegin()->second;

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_vmax_bndy_curve(c);

              pwc.push_back(c,du);
            }
          }

          void get_umin_ndelta_pcurve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator vit;
            data_type vmin = vkey.get_pmin();

            uk = ukey.key.begin()->second;

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_umin_ndelta_pcurve(c);

              pwc.push_back(c,dv);
            }
          }

          void get_umax_ndelta_pcurve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator vit;
            data_type vmin = vkey.get_pmin();

            uk = ukey.key.rbegin()->second;

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_umax_ndelta_pcurve(c);

              pwc.push_back(c,dv);
            }
          }

          void get_vmin_ndelta_pcurve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit;
            data_type umin = ukey.get_pmin();

            vk = vkey.key.begin()->second;

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_vmin_ndelta_pcurve(c);

              pwc.push_back(c,du);
            }
          }

          void get_vmax_ndelta_pcurve( piecewise_curve_type &pwc ) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit;
            data_type umin = ukey.get_pmin();

            vk = vkey.key.rbegin()->second;

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_vmax_ndelta_pcurve(c);

              pwc.push_back(c,du);
            }
          }

          void get_uconst_f_u_curve(piecewise_curve_type &pwc, const data_type &u) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type vmin = vkey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, u, vmin);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_uconst_f_u_curve(c, uu);

              pwc.push_back(c,dv);
            }
          }

          void get_uconst_f_v_curve(piecewise_curve_type &pwc, const data_type &u) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type vmin = vkey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, u, vmin);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(vmin);

            for ( vit = vkey.key.begin(); vit != vkey.key.end(); ++vit )
            {
              vk = vit->second;

              data_type dv=vkey.get_delta_parm(vit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_uconst_f_v_curve(c, uu);

              pwc.push_back(c,dv);
            }
          }

          void get_vconst_curve(piecewise_curve_type &pwc, const data_type &v) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type umin = ukey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, umin, v);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_vconst_curve(c, vv);

              pwc.push_back(c,du);
            }
          }

          void get_vconst_f_u_curve(piecewise_curve_type &pwc, const data_type &v) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type umin = ukey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, umin, v);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_vconst_f_u_curve(c, vv);

              pwc.push_back(c,du);
            }
          }

          void get_vconst_f_v_curve(piecewise_curve_type &pwc, const data_type &v) const
          {
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);
            data_type umin = ukey.get_pmin();

            find_patch(uk, vk, uit, vit, uu, vv, umin, v);

            assert ((uk != -1) && (vk != -1));

            pwc.clear();
            pwc.set_t0(umin);

            for ( uit = ukey.key.begin(); uit != ukey.key.end(); ++uit )
            {
              uk = uit->second;

              data_type du=ukey.get_delta_parm(uit);

              surface_type s=patches[uk][vk];

              curve_type c;

              s.get_vconst_f_v_curve(c, vv);

              pwc.push_back(c,du);
            }
          }
          void find_interior_feature_edges(std::vector<data_type> &uconst, std::vector<data_type> &vconst, const data_type &angle_tol) const
          {
            index_type nu, nv, iu, iv;

            piecewise_curve_type c;
            std::vector<data_type> pmap, ldis, ldis_out;
            tolerance_type tol;

            // initialize the output
            uconst.clear();
            vconst.clear();

            // unnamed function to test if two parameters are close enough
            auto comp = [&tol](const data_type &x1, const data_type &x2)->bool
            {
              return tol.approximately_less_than(x1, x2);
            };

            // extract each v-const curve that is an edge of one of the patches to find u-parameters
            // that contain C0 only edges
            get_pmap_v(pmap);
            nv=pmap.size();
            assert(nv-1==number_v_patches());
            for (iv=0; iv<nv; ++iv)
            {
              get_vconst_curve(c, pmap[iv]);
              c.find_discontinuities(angle_tol, ldis);

              // merge these parameters with current list
              ldis_out.clear();
              std::set_union(uconst.begin(), uconst.end(), ldis.begin(), ldis.end(), std::back_inserter(ldis_out), comp);
              std::swap(uconst, ldis_out);
            }

            // extract each u-const curve that is an edge of one of the patches to find v-parameters
            // that contain C0 only edges
            pmap.clear();
            get_pmap_u(pmap);
            nu=pmap.size();
            assert(nu-1==number_u_patches());
            for (iu=0; iu<nu; ++iu)
            {
              get_uconst_curve(c, pmap[iu]);
              c.find_discontinuities(angle_tol, ldis);

              // merge these parameters with current list
              ldis_out.clear();
              std::set_union(vconst.begin(), vconst.end(), ldis.begin(), ldis.end(), std::back_inserter(ldis_out), comp);
              std::swap(vconst, ldis_out);
            }

            // TODO: Need to compare actual control points next to edges to catch cases where the
            //       patch corners satisfy the constraints but internally the constraints are not.

          }

          void find_interior_C0_edges(std::vector<data_type> &uconst, std::vector<data_type> &vconst) const
          {
            index_type nu, nv, iu, iv;

            piecewise_curve_type c;
            std::vector<data_type> pmap, ldis, ldis_out;
            tolerance_type tol;

            // initialize the output
            uconst.clear();
            vconst.clear();

            // unnamed function to test if two parameters are close enough
            auto comp = [&tol](const data_type &x1, const data_type &x2)->bool
            {
              return tol.approximately_less_than(x1, x2);
            };

            // extract each v-const curve that is an edge of one of the patches to find u-parameters
            // that contain C0 only edges
            get_pmap_v(pmap);
            nv=pmap.size();
            assert(nv-1==number_v_patches());
            for (iv=0; iv<nv; ++iv)
            {
              get_vconst_curve(c, pmap[iv]);
              c.find_discontinuities(eli::geom::general::G1, ldis);

              // merge these parameters with current list
              ldis_out.clear();
              std::set_union(uconst.begin(), uconst.end(), ldis.begin(), ldis.end(), std::back_inserter(ldis_out), comp);
              std::swap(uconst, ldis_out);
            }

            // extract each u-const curve that is an edge of one of the patches to find v-parameters
            // that contain C0 only edges
            pmap.clear();
            get_pmap_u(pmap);
            nu=pmap.size();
            assert(nu-1==number_u_patches());
            for (iu=0; iu<nu; ++iu)
            {
              get_uconst_curve(c, pmap[iu]);
              c.find_discontinuities(eli::geom::general::G1, ldis);

              // merge these parameters with current list
              ldis_out.clear();
              std::set_union(vconst.begin(), vconst.end(), ldis.begin(), ldis.end(), std::back_inserter(ldis_out), comp);
              std::swap(vconst, ldis_out);
            }

            // TODO: Need to compare actual control points next to edges to catch cases where the
            //       patch corners are continuous but internally it is not.
          }

          point_type f(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            return patches[uk][vk].f(uu, vv);
          }

          point_type f_u(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);

            return patches[uk][vk].f_u(uu, vv)/delta_u;
          }

          point_type f_v(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_v(uu, vv)/delta_v;
          }

          point_type f_uu(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);

            return patches[uk][vk].f_uu(uu, vv)/(delta_u*delta_u);
          }

          point_type f_uv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);
            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_uv(uu, vv)/(delta_u*delta_v);
          }

          point_type f_vv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_vv(uu, vv)/(delta_v*delta_v);
          }

          point_type f_uuu(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);

            return patches[uk][vk].f_uuu(uu, vv)/(delta_u*delta_u*delta_u);
          }

          point_type f_uuv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);
            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_uuv(uu, vv)/(delta_u*delta_u*delta_v);
          }

          point_type f_uvv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);
            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_uvv(uu, vv)/(delta_u*delta_v*delta_v);
          }

          point_type f_vvv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            data_type delta_u = ukey.get_delta_parm(uit);
            data_type delta_v = vkey.get_delta_parm(vit);

            return patches[uk][vk].f_vvv(uu, vv)/(delta_v*delta_v*delta_v);
          }

          point_type normal(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            return patches[uk][vk].normal(uu, vv);
          }

          void f_pt_normal(const data_type &u, const data_type &v, point_type &pt, point_type &norm, const index_type &utie = 0, const index_type &vtie = 0 ) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            data_type uu(0), vv(0);
            typename keymap_type::const_iterator uit, vit;

            patch_boundary_code_type bcode = find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            // Handle u boundary code
            if ( bcode.first == -1 && utie == 1 ) // u is at start of interval, uu = 0.
            {
              if ( uit != ukey.key.begin() )
              {
                uit--;
                uk = uit->second;
                uu = 1.0;
              }
            }
            else if ( bcode.first == 1 && utie == -1 ) // u is at end of interval uu = 1;
            {
              typename keymap_type::const_iterator uitnext = uit;
              uitnext++;

              if ( uitnext != ukey.key.end() )
              {
                uit = uitnext;
                uk = uit->second;
                uu = 0.0;
              }
            }

            // Handle v boundary code
            if ( bcode.second == -1 && vtie == 1 ) // v is at start of interval, vv = 0.
            {
              if ( vit != vkey.key.begin() )
              {
                vit--;
                vk = vit->second;
                vv = 1.0;
              }
            }
            else if ( bcode.second == 1 && vtie == -1) // v is at end of interval vv = 1;
            {
              typename keymap_type::const_iterator vitnext = vit;
              vitnext++;

              if ( vitnext != vkey.key.end() )
              {
                vit = vitnext;
                vk = vit->second;
                vv = 0.0;
              }
            }

            pt = patches[uk][vk].f(uu, vv);
            norm = patches[uk][vk].normal(uu, vv);
          }

          void f_pt_normal_grid(const std::vector < data_type > &uvec, const std::vector < data_type > &vvec, std::vector < std::vector < point_type > > &ptmat, std::vector < std::vector < point_type > > &normmat ) const
          {
            typedef std::vector < data_type > pvec_type;

            pvec_type uuvec( uvec.size() );
            std::vector < index_type > ukvec( uvec.size() );
            pvec_type vvvec( vvec.size() );
            std::vector < index_type > vkvec( vvec.size() );

            int i, j;

            for ( i = 0; i < uvec.size(); i++ )
            {
              typename keymap_type::const_iterator uit;
              index_type uk;
              data_type uu;
              index_type code = ukey.find_segment( uk, uit, uu, uvec[i] );

              if ( code == -1 && i == ( uvec.size() - 1 ) ) // u is at start of interval, uu = 0.
              {
                if ( uit != ukey.key.begin() )
                {
                  uit--;
                  uk = uit->second;
                  uu = 1.0;
                }
              }
              else if ( code == 1 && i == 0 ) // u is at end of interval uu = 1;
              {
                typename keymap_type::const_iterator uitnext = uit;
                uitnext++;

                if ( uitnext != ukey.key.end() )
                {
                  uit = uitnext;
                  uk = uit->second;
                  uu = 0.0;
                }
              }

              uuvec[i] = uu;
              ukvec[i] = uk;
            }

            for ( i = 0; i < vvec.size(); i++ )
            {
              typename keymap_type::const_iterator vit;
              index_type vk;
              data_type vv;
              index_type code = vkey.find_segment( vk, vit, vv, vvec[i] );

              if ( code == -1 && i == ( vvec.size() - 1 ) ) // v is at start of interval, vv = 0.
              {
                if ( vit != vkey.key.begin() )
                {
                  vit--;
                  vk = vit->second;
                  vv = 1.0;
                }
              }
              else if ( code == 1 && i == 0 ) // v is at end of interval vv = 1;
              {
                typename keymap_type::const_iterator vitnext = vit;
                vitnext++;

                if ( vitnext != vkey.key.end() )
                {
                  vit = vitnext;
                  vk = vit->second;
                  vv = 0.0;
                }
              }

              vvvec[i] = vv;
              vkvec[i] = vk;
            }

            ptmat.resize( uvec.size() );
            normmat.resize( uvec.size() );
            for ( i = 0; i < uvec.size(); i++ )
            {
              ptmat[i].resize( vvec.size() );
              normmat[i].resize( vvec.size() );

              for ( j = 0; j < vvec.size(); j++ )
              {
                ptmat[i][j] = patches[ukvec[i]][vkvec[j]].f(uuvec[i], vvvec[j]);
                normmat[i][j] = patches[ukvec[i]][vkvec[j]].normal(uuvec[i], vvvec[j]);
              }
            }
          }

          void f_pt_derivs(const data_type &u, const data_type &v, point_type &pt, point_type &pt_u, point_type &pt_v ) const
          {
            // find patch that corresponds to given u & v
            index_type uk, vk;
            typename keymap_type::const_iterator uit, vit;
            data_type uu(0), vv(0);

            find_patch(uk, vk, uit, vit, uu, vv, u, v);

            assert((uk != -1) && (vk != -1));

            pt = patches[uk][vk].f(uu, vv);

            data_type delta_u = ukey.get_delta_parm(uit);

            pt_u = patches[uk][vk].f_u(uu, vv)/delta_u;

            data_type delta_v = vkey.get_delta_parm(vit);

            pt_v = patches[uk][vk].f_v(uu, vv)/delta_v;
          }

          // TODO: NEED TO IMPLEMENT
          //       * fit
          //       * interpolate

        private:
//           template<template<typename, unsigned short, typename> class surf1__,
//                    typename data1__, unsigned short dim1__, typename tol1__>
//           friend void area(typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &len,
//                            const piecewise<surf1__, data1__, dim1__, tol1__> &pc,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &tol);
//           template<template<typename, unsigned short, typename> class surf1__,
//                             typename data1__, unsigned short dim1__, typename tol1__>
//           friend void area(typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &len,
//                            const piecewise<surf1__, data1__, dim1__, tol1__> &pc,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &t0,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &t1,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &tol);

          template<template<typename, unsigned short, typename> class surface1__, typename data1__, unsigned short dim1__, typename tol1__ >
          friend typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type
            eli::geom::intersect::minimum_distance(
              typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type &u,
              typename piecewise<surface1__, data1__, dim1__, tol1__>::data_type &v,
              const piecewise<surface1__, data1__, dim1__, tol1__> &ps,
              const typename piecewise<surface1__, data1__, dim1__, tol1__>::point_type &pt);

          typedef std::map< data_type, index_type > keymap_type;

          struct parameter_key
          {
            keymap_type key;
            data_type pmax;

            parameter_key() : pmax(0) {}
            parameter_key(const parameter_key &pk) : key(pk.key), pmax(pk.pmax) {}
            ~parameter_key() {}

            bool operator==(const parameter_key &pk) const
            {
              if (this==&pk)
                return true;
              if (pmax!=pk.pmax)
                return false;
              if (key!=pk.key)
                return false;

              return true;
            }

            bool operator!=(const parameter_key &pk) const
            {
              return !operator==(pk);
            }

            void clear()
            {
              pmax=0;
              key.clear();
            }

            data_type get_pmax() const
            {
              return pmax;
            }

            data_type get_pmin() const
            {
              if(!key.empty())
                return key.begin()->first;
              else
                return pmax;
            }

            void set_pmax(const data_type &pmax_in)
            {
              pmax = pmax_in;
            }

            void set_pmin(const data_type &pmin_in)
            {
              if(!key.empty())
              {
                if(pmin_in != key.begin()->first)
                {
                  data_type p = pmin_in;
                  keymap_type shiftkey;
                  for (typename keymap_type::iterator it=key.begin(); it!=key.end(); ++it)
                  {
                    data_type delta_p = get_delta_parm(it);

                    shiftkey.insert(shiftkey.end(), std::make_pair(p, it->second));

                    p+=delta_p;
                  }
                  key.swap(shiftkey);
                  pmax = p;
                }
              }
              else
              {
                pmax=pmin_in;
              }
            }

            void init(const index_type &nseg, const data_type &dp = 1, const data_type &p0 = 0)
            {
              key.clear();
              pmax = p0;
              append(nseg, dp);
            }

            void append(const data_type &dp = 1)
            {
              typename keymap_type::iterator itguess = key.end();
              itguess = key.insert(itguess, std::make_pair(pmax, key.size()));
              pmax += dp;
            }

            void append(const index_type &nseg, const data_type &dp = 1)
            {
              typename keymap_type::iterator itguess = key.end();
              index_type j = key.size();
              data_type p = pmax;
              for(index_type i = 0; i < nseg; ++i)
              {
                itguess = key.insert(itguess, std::make_pair(p, j));
                p += dp;
                ++j;
              }
              pmax = p;
            }

            template<typename it__>
            void init(const it__ &dps, const it__ &dpe, const data_type &p0 = 0)
            {
              key.clear();
              pmax = p0;
              append(dps, dpe);
            }

            template<typename it__>
            void append(const it__ &dps, const it__ &dpe)
            {
              typename keymap_type::iterator itguess = key.end();
              index_type j = key.size();
              data_type p = pmax;
              for(it__ dp = dps; dp != dpe; ++dp)
              {
                itguess = key.insert(itguess, std::make_pair(p, j));
                p += (*dp);
                ++j;
              }
              pmax = p;
            }

            void parameter_report() const
            {
              printf("Parameter report:\n");
              typename keymap_type::const_iterator it;

              int i = 0;
              // cycle through all segments to get each bounding box to add
              for (it=key.begin(); it!=key.end(); ++it)
              {
                printf(" seg: %d \t p: %f \t pk %d\n", i, it->first, it->second);
                ++i;
              }
              printf(" pmax: %f\n", pmax);
              printf("End report\n");
            }

            void get_pmap(std::vector<data_type> &pmap) const
            {
              pmap.clear();

              typename keymap_type::const_iterator it;
              for (it=key.cbegin(); it!=key.cend(); ++it)
              {
                pmap.push_back( it->first );
              }
              pmap.push_back( pmax );
            }

            void reverse_keymap()
            {
              keymap_type rkey;
              typename keymap_type::iterator itr;
              typename keymap_type::iterator itrguess = rkey.begin();

              data_type p = get_pmin();

              for (typename keymap_type::reverse_iterator it=key.rbegin(); it!=key.rend(); ++it)
              {
                itr = rkey.insert(itrguess, std::make_pair(p, it->second));

                data_type delta_p = get_delta_parm(it);
                p += delta_p;

                itrguess = itr;
              }
              key.swap(rkey);

              // Parametric length should stay the same.
              tolerance_type tol;
              assert(tol.approximately_equal(p, pmax));
            }


            data_type get_delta_parm(const typename keymap_type::iterator &it) const
            {
              assert (it != key.end());

              typename keymap_type::iterator itnext = it;
              itnext++;

              data_type delta_p;

              if(itnext != key.end())
                delta_p = itnext->first - it->first;
              else
                delta_p = pmax - it->first;

              return delta_p;
            }

            data_type get_delta_parm(const typename keymap_type::const_iterator &it) const
            {
              assert (it != key.end());

              typename keymap_type::const_iterator itnext = it;
              itnext++;

              data_type delta_p;

              if(itnext != key.end())
                delta_p = itnext->first - it->first;
              else
                delta_p = pmax - it->first;

              return delta_p;
            }

            data_type get_delta_parm(const typename keymap_type::reverse_iterator &it) const
            {
              assert (it != key.rend());

              data_type delta_p;

              if(it != key.rbegin())
              {
                typename keymap_type::reverse_iterator itprev = it;
                itprev--;
                delta_p = itprev->first - it->first;
              }
              else
              {
                delta_p = pmax - it->first;
              }

              return delta_p;
            }

            data_type get_delta_parm(const typename keymap_type::const_reverse_iterator &it) const
            {
              assert (it != key.rend());

              data_type delta_p;

              if(it != key.rbegin())
              {
                typename keymap_type::const_reverse_iterator itprev = it;
                itprev--;
                delta_p = itprev->first - it->first;
              }
              else
              {
                delta_p = pmax - it->first;
              }

              return delta_p;
            }

            void find_segment(index_type &ikey, typename keymap_type::const_iterator &it, const index_type &index) const
            {
              if(index >= (int) key.size() || index < 0)
              {
                it=key.end();
                ikey=-1;
                return;
              }

              // advance to desired index
              index_type i;
              for (i=0, it=key.begin(); i<index; ++i, ++it) {}

              ikey=it->second;
            }

            void find_segment(index_type &ikey, typename keymap_type::iterator &it, const index_type &index) const
            {
              if(index >= (int) key.size() || index < 0)
              {
                it=key.end();
                ikey=-1;
                return;
              }

              // advance to desired index
              index_type i;
              for (i=0, it=key.begin(); i<index; ++i, ++it) {}

              ikey=it->second;
            }

            index_type find_segment(index_type &ikey, typename keymap_type::iterator &it, data_type &pp, const data_type &p_in)
            {
              tol__ tol;

              if(p_in>pmax)
              {
                it=key.end();
                ikey = -1;
                return 0;
              }

              data_type pmin = get_pmin();

              if(p_in<pmin)
              {
                it=key.end();
                ikey = -1;
                return 0;
              }

              // Use map::upper_bound for fast lookup of segment after p_in
              it=key.upper_bound(p_in);

              // Decrement to segment containing p_in
              if(it != key.begin())
                it--;

              ikey = it->second;

              // At start of segment
              if(tol.approximately_equal(p_in, it->first))
              {
                pp=static_cast<data_type>(0);
                return -1;
              }

              data_type delta_p = get_delta_parm(it);

              // At end of segment
              if(tol.approximately_equal(p_in, it->first + delta_p))
              {
                pp=static_cast<data_type>(1);
                return 1;
              }

              // Typical case
              pp=(p_in-it->first)/delta_p;

              // Super careful checks
              if (pp>static_cast<data_type>(1))
                pp=static_cast<data_type>(1);
              if (pp<static_cast<data_type>(0))
                pp=static_cast<data_type>(0);

              return 0;
            }

            index_type find_segment(index_type &ikey, typename keymap_type::const_iterator &it, data_type &pp, const data_type &p_in) const
            {
              tol__ tol;

              if(p_in>pmax)
              {
                it=key.end();
                ikey = -1;
                return 0;
              }

              data_type pmin = get_pmin();

              if(p_in<pmin)
              {
                it=key.end();
                ikey = -1;
                return 0;
              }

              // Use map::upper_bound for fast lookup of segment after p_in
              it=key.upper_bound(p_in);

              // Decrement to segment containing p_in
              if(it != key.begin())
                it--;

              ikey = it->second;

              // At start of segment
              if(tol.approximately_equal(p_in, it->first))
              {
                pp=static_cast<data_type>(0);
                return -1;
              }

              data_type delta_p = get_delta_parm(it);

              // At end of segment
              if(tol.approximately_equal(p_in, it->first + delta_p))
              {
                pp=static_cast<data_type>(1);
                return 1;
              }

              // Typical case
              pp=(p_in-it->first)/delta_p;

              // Super careful checks
              if (pp>static_cast<data_type>(1))
                pp=static_cast<data_type>(1);
              if (pp<static_cast<data_type>(0))
                pp=static_cast<data_type>(0);

              return 0;
            }
          };


          typedef std::vector< surface_type > patch_strip_type;
          typedef std::vector< patch_strip_type > patch_collection_type;

          patch_collection_type patches;
          // By convention, patches[uk][vk]

          parameter_key ukey, vkey;
          index_type nu, nv;

        protected:
          bool check_continuity(const eli::geom::general::continuity &/*cont*/) const
          {
            // TODO: Need to implement this
            return true;
          }

          bool check_u_continuity(const surface_type &/*s1*/, const surface_type &/*s2*/, const eli::geom::general::continuity &/*cont*/) const
          {
            // TODO: Need to implement this
            return true;
          }

          bool check_v_continuity(const surface_type &/*s1*/, const surface_type &/*s2*/, const eli::geom::general::continuity &/*cont*/) const
          {
            // TODO: Need to implement this
            return true;
          }

        private:

          void resize_store(const index_type &nu_in, const index_type &nv_in)
          {
            if ((nu_in<=0) || (nv_in<=0))
              return;

            patches.resize(nu_in);
            nu = nu_in;

            // Unconditionally do this to make sure newly added rows are properly sized.
            for(index_type i = 0; i < nu_in; i++)
              patches[i].resize(nv_in);

            nv = nv_in;
          }

          error_code subsurf(piecewise<surface__, data_type, dim__, tol__> &surf, const typename keymap_type::const_iterator &ustart, const typename keymap_type::const_iterator &uend, const typename keymap_type::const_iterator &vstart, const typename keymap_type::const_iterator &vend ) const
          {
            surf.clear();

            surf.set_u0( ustart->first );
            surf.set_v0( vstart->first );

            typename keymap_type::const_iterator uit, vit;

            index_type nusub = 0;
            for ( uit = ustart; uit != uend; uit++ )
            {
              data_type du = ukey.get_delta_parm( uit );
              surf.ukey.append( du );
              nusub++;
            }

            index_type nvsub = 0;
            for ( vit = vstart; vit != vend; vit++ )
            {
              data_type dv = vkey.get_delta_parm( vit );
              surf.vkey.append( dv );
                nvsub++;
            }

            surf.resize_store( nusub, nvsub );

            index_type ikstore, jkstore;
            ikstore = 0;
            for ( uit = ustart; uit != uend; uit++ )
            {
              jkstore = 0;
              for ( vit = vstart; vit != vend; vit++ )
              {
                surf.patches[ikstore][jkstore] = patches[(*uit).second][(*vit).second];
                jkstore++;
              }
              ikstore++;
            }

            return NO_ERRORS;
          }

          error_code split_u(const index_type &uk, const typename keymap_type::iterator &uit, const data_type &u_in, const data_type &uu)
          {
            tolerance_type tol;
            assert(!tol.approximately_equal(uu, 0));
            assert(!tol.approximately_equal(uu, 1));

            index_type ukr, vk;
            // Right half will be added at end of patch matrix.
            ukr=nu;
            ukey.key.insert(uit, std::make_pair(u_in, ukr));

            // Increase matrix size.
            resize_store(nu+1, nv);

            for (vk=0; vk<nv; ++vk)
            {
              surface_type s = patches[uk][vk];
              s.split_u(patches[uk][vk], patches[ukr][vk], uu);
            }

            return NO_ERRORS;
          }

          error_code split_v(const index_type &vk, const typename keymap_type::iterator &vit, const data_type &v_in, const data_type &vv)
          {
            tolerance_type tol;
            assert(!tol.approximately_equal(vv, 0));
            assert(!tol.approximately_equal(vv, 1));

            index_type uk, vkr;
            // Right half will be added at end of patch matrix.
            vkr=nv;
            vkey.key.insert(vit, std::make_pair(v_in, vkr));

            // Increase matrix size.
            resize_store(nu, nv+1);

            for (uk=0; uk<nu; ++uk)
            {
              surface_type s = patches[uk][vk];
              s.split_v(patches[uk][vk], patches[uk][vkr], vv);
            }

            return NO_ERRORS;
          }

          // Lookup based on i,j
          void find_patch(index_type &uk, index_type &vk,
                          typename keymap_type::iterator &uit, typename keymap_type::iterator &vit,
                          const index_type & ui, const index_type &vi)
          {
            ukey.find_segment(uk, uit, ui);
            vkey.find_segment(vk, vit, vi);
          }

          void find_patch(typename keymap_type::iterator &uit, typename keymap_type::iterator &vit,
                          const index_type & ui, const index_type &vi)
          {
            index_type uk, vk;
            find_patch(uk, vk, uit, vit, ui, vi);
          }

          void find_patch(index_type &uk, index_type &vk,
                          typename keymap_type::const_iterator &uit, typename keymap_type::const_iterator &vit,
                          const index_type & ui, const index_type &vi) const
          {
            ukey.find_segment(uk, uit, ui);
            vkey.find_segment(vk, vit, vi);
          }

          void find_patch(typename keymap_type::const_iterator &uit, typename keymap_type::const_iterator &vit,
                          const index_type & ui, const index_type &vi) const
          {
            index_type uk, vk;
            find_patch(uk, vk, uit, vit, ui, vi);
          }

          void find_patch(index_type &uk, index_type &vk,
                          const index_type & ui, const index_type &vi) const
          {
            typename keymap_type::const_iterator uit, vit;
            find_patch(uk, vk, uit, vit, ui, vi);
          }

          // Lookup based on u_in, v_in.
          patch_boundary_code_type find_patch(index_type &uk, index_type &vk,
                          typename keymap_type::iterator &uit, typename keymap_type::iterator &vit,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in)
          {
            index_type ucode, vcode;
            ucode = ukey.find_segment(uk, uit, uu, u_in);
            vcode = vkey.find_segment(vk, vit, vv, v_in);
            return std::make_pair( ucode, vcode );
          }

          patch_boundary_code_type find_patch(typename keymap_type::iterator &uit, typename keymap_type::iterator &vit,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in)
          {
            index_type uk, vk;
            return find_patch(uk, vk, uit, vit, uu, vv, u_in, v_in);
          }

          patch_boundary_code_type find_patch(index_type &uk, index_type &vk,
                          typename keymap_type::const_iterator &uit, typename keymap_type::const_iterator &vit,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in) const
          {
            index_type ucode, vcode;
            ucode = ukey.find_segment(uk, uit, uu, u_in);
            vcode = vkey.find_segment(vk, vit, vv, v_in);
            return std::make_pair( ucode, vcode );
          }

          patch_boundary_code_type find_patch(typename keymap_type::const_iterator &uit, typename keymap_type::const_iterator &vit,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in) const
          {
            index_type uk, vk;
            return find_patch(uk, vk, uit, vit, uu, vv, u_in, v_in);
          }

          patch_boundary_code_type find_patch(index_type &uk, index_type &vk,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in) const
          {
            typename keymap_type::const_iterator uit, vit;
            return find_patch(uk, vk, uit, vit, uu, vv, u_in, v_in);
          }

      };
    }
  }
}
#endif
