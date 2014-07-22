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

#ifndef eli_geom_surface_piecewise_general_skinning_surface_creator_hpp
#define eli_geom_surface_piecewise_general_skinning_surface_creator_hpp

#include <list>
#include <vector>
#include <iterator>

#include "eli/util/tolerance.hpp"

#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_general_creator.hpp"

#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/bezier.hpp"
#include "eli/geom/surface/piecewise_connection_data.hpp"
#include "eli/geom/surface/piecewise_creator_base.hpp"

#include <string>

namespace eli
{
  inline
  std::string random_string( size_t length )
  {
      auto randchar = []() -> char
      {
          const char charset[] =
          "0123456789"
          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
          "abcdefghijklmnopqrstuvwxyz";
          const size_t max_index = (sizeof(charset) - 1);
          return charset[ rand() % max_index ];
      };
      std::string str(length,0);
      std::generate_n( str.begin(), length, randchar );
      return str;
  }

  inline
  void octave_start(int figno)
  {
    std::cout << "figure(" << figno << ");" << std::endl;
    std::cout << "clf(" << figno << ", 'reset');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "grid on;" << std::endl;
  }

  inline
  void octave_finish(int figno)
  {
    std::cout << "figure(" << figno << ");" << std::endl;
    std::cout << "hold off;" << std::endl;
    std::cout << "rotate3d on;" << std::endl;
    std::cout << "xlabel('x');" << std::endl;
    std::cout << "ylabel('y');" << std::endl;
    std::cout << "zlabel('z');" << std::endl;
  }

  template<typename data__>
  void octave_print(int figno, const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &pc,
                    const std::string &name="", bool show_control_points=true)
  {
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;

    std::string nm, cpxbuf, cpybuf, cpzbuf, cxbuf, cybuf, czbuf;

    index_type i, pp, ns;
    data_type tmin, tmax;

    ns=pc.number_segments();
    pc.get_parameter_min(tmin);
    pc.get_parameter_max(tmax);

    // build name
    if (name=="")
    {
      nm=random_string(5);
    }
    else
    {
      nm=name;
    }

    // set control points
    if (show_control_points)
    {
      cpxbuf=nm+"_curv_cp_x=[";
      cpybuf=nm+"_curv_cp_y=[";
      cpzbuf=nm+"_curv_cp_z=[";
      for (pp=0; pp<ns; ++pp)
      {
        index_type bez_deg;
        curve_type bez;
        pc.get(bez, pp);
        bez_deg=bez.degree();
        for (i=0; i<=bez_deg; ++i)
        {
          cpxbuf+=std::to_string(bez.get_control_point(i).x());
          cpybuf+=std::to_string(bez.get_control_point(i).y());
          cpzbuf+=std::to_string(bez.get_control_point(i).z());
          if ((pp<(ns-1)) || ((pp==(ns-1)) && (i<bez_deg)))
          {
            cpxbuf+=", ";
            cpybuf+=", ";
            cpzbuf+=", ";
          }
        }
      }
      cpxbuf+="];";
      cpybuf+="];";
      cpzbuf+="];";
    }

    // initialize the t parameters
    index_type nt(129);
    std::vector<data__> t(nt);
    for (i=0; i<nt; ++i)
    {
      t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(nt-1);
    }

    // set the curve points
    cxbuf=nm+"_curv_x=[";
    cybuf=nm+"_curv_y=[";
    czbuf=nm+"_curv_z=[";
    for (i=0; i<nt; ++i)
    {
      cxbuf+=std::to_string(pc.f(t[i]).x());
      cybuf+=std::to_string(pc.f(t[i]).y());
      czbuf+=std::to_string(pc.f(t[i]).z());
      if (i<nt)
      {
        cxbuf+=", ";
        cybuf+=", ";
        czbuf+=", ";
      }
    }
    cxbuf+="];";
    cybuf+="];";
    czbuf+="];";

    std::cout << "% curve: " << nm << std::endl;
    std::cout << "figure(" << figno << ");" << std::endl;
    std::cout << cxbuf << std::endl;
    std::cout << cybuf << std::endl;
    std::cout << czbuf << std::endl;
    if (show_control_points)
    {
      std::cout << cpxbuf << std::endl;
      std::cout << cpybuf << std::endl;
      std::cout << cpzbuf << std::endl;
    }
    std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
    std::cout << "plot3(" << nm << "_curv_x, "
                          << nm << "_curv_y, "
                          << nm << "_curv_z, '-g');" << std::endl;
    if (show_control_points)
    {
      std::cout << "plot3(" << nm << "_curv_cp_x', "
                            << nm << "_curv_cp_y', "
                            << nm << "_curv_cp_z', '-o', 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0]);" << std::endl;
    }
  }

  template<typename data__>
  void octave_print(int figno, const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &pc,
                    const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &vec,
                    const std::string &name="")
  {
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;

    std::string nm(random_string(5)), vecxbuf, vecybuf, veczbuf, cxbuf, cybuf, czbuf;

    index_type i, pp, ns;
    data_type tmin, tmax;

    tolerance_type tol;

    ns=pc.number_segments();
    pc.get_parameter_min(tmin);
    pc.get_parameter_max(tmax);

    // check parameterization of vec curve
    if (!tol.approximately_equal(vec.get_t0(), tmin))
    {
      return;
    }
    if (!tol.approximately_equal(vec.get_tmax(), tmax))
    {
      return;
    }

    // build name
    if (name=="")
    {
      nm=random_string(5);
    }
    else
    {
      nm=name;
    }

    // initialize the t parameters
    index_type nt(11);
    std::vector<data__> t(nt);
    for (i=0; i<nt; ++i)
    {
      t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(nt-1);
    }

    // set the surface points
    cxbuf=nm+"_loc_x=[";
    cybuf=nm+"_loc_y=[";
    czbuf=nm+"_loc_z=[";
    vecxbuf=nm+"_vec_x=[";
    vecybuf=nm+"_vec_y=[";
    veczbuf=nm+"_vec_z=[";
    for (i=0; i<nt; ++i)
    {
      cxbuf+=std::to_string(pc.f(t[i]).x());
      cybuf+=std::to_string(pc.f(t[i]).y());
      czbuf+=std::to_string(pc.f(t[i]).z());
      vecxbuf+=std::to_string(vec.f(t[i]).x());
      vecybuf+=std::to_string(vec.f(t[i]).y());
      veczbuf+=std::to_string(vec.f(t[i]).z());
      if (i<nt)
      {
        cxbuf+=", ";
        cybuf+=", ";
        czbuf+=", ";
        vecxbuf+=", ";
        vecybuf+=", ";
        veczbuf+=", ";
      }
    }
    cxbuf+="];";
    cybuf+="];";
    czbuf+="];";
    vecxbuf+="];";
    vecybuf+="];";
    veczbuf+="];";

    std::cout << "% curve: " << nm << std::endl;
    std::cout << "figure(" << figno << ");" << std::endl;
    std::cout << cxbuf << std::endl;
    std::cout << cybuf << std::endl;
    std::cout << czbuf << std::endl;
    std::cout << vecxbuf << std::endl;
    std::cout << vecybuf << std::endl;
    std::cout << veczbuf << std::endl;
    std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
    std::cout << "quiver3(" << nm << "_loc_x, "
                            << nm << "_loc_y, "
                            << nm << "_loc_z, "
                            << nm << "_vec_x, "
                            << nm << "_vec_y, "
                            << nm << "_vec_z, "
                            << "'r');" << std::endl;
  }

  template<typename data__>
  void octave_print(int figno, const eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> &ps,
                    const std::string &name="", bool show_control_points=true)
  {
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;

    std::string nm, sxbuf, sybuf, szbuf;
    std::vector<std::string> cpxbuf, cpybuf, cpzbuf;

    index_type i, j, pp, qq, nup, nvp;
    data_type umin, vmin, umax, vmax;

    nup=ps.number_u_patches();
    nvp=ps.number_v_patches();
    ps.get_parameter_min(umin, vmin);
    ps.get_parameter_max(umax, vmax);


    // build name
    if (name=="")
    {
      nm=random_string(5);
    }
    else
    {
      nm=name;
    }

    // set control points
    if (show_control_points)
    {
      cpxbuf.resize(nup*nvp);
      cpybuf.resize(nup*nvp);
      cpzbuf.resize(nup*nvp);

      for (pp=0; pp<nup; ++pp)
      {
        for (qq=0; qq<nvp; ++qq)
        {
          surface_type bez;
          ps.get(bez, pp, qq);
          index_type ppqq;

          ppqq=qq*nup+pp;
          cpxbuf[ppqq]=nm+"_surf_cp"+std::to_string(pp)+std::to_string(qq)+"_x=[";
          cpybuf[ppqq]=nm+"_surf_cp"+std::to_string(pp)+std::to_string(qq)+"_y=[";
          cpzbuf[ppqq]=nm+"_surf_cp"+std::to_string(pp)+std::to_string(qq)+"_z=[";

          for (i=0; i<=bez.degree_u(); ++i)
          {
            for (j=0; j<=bez.degree_v(); ++j)
            {
              cpxbuf[ppqq]+=std::to_string(bez.get_control_point(i, j).x());
              cpybuf[ppqq]+=std::to_string(bez.get_control_point(i, j).y());
              cpzbuf[ppqq]+=std::to_string(bez.get_control_point(i, j).z());

              if ((i==bez.degree_u()) && (j==bez.degree_v()))
              {
              }
              else if ((j==bez.degree_v()) && (i<bez.degree_u()))
              {
                cpxbuf[ppqq]+="; ";
                cpybuf[ppqq]+="; ";
                cpzbuf[ppqq]+="; ";
              }
              else
              {
                cpxbuf[ppqq]+=", ";
                cpybuf[ppqq]+=", ";
                cpzbuf[ppqq]+=", ";
              }
            }
          }
          cpxbuf[ppqq]+="];";
          cpybuf[ppqq]+="];";
          cpzbuf[ppqq]+="];";
        }
      }
    }

    // initialize the u & v parameters
    index_type nu(10*nup+1), nv(10*nvp+1);
    std::vector<data__> u(nu), v(nv);
    for (i=0; i<static_cast<index_type>(u.size()); ++i)
    {
      u[i]=umin+(umax-umin)*static_cast<data__>(i)/(u.size()-1);
    }
    for (j=0; j<static_cast<index_type>(v.size()); ++j)
    {
      v[j]=vmin+(vmax-vmin)*static_cast<data__>(j)/(v.size()-1);
    }

    // set the surface points
    sxbuf=nm+"_surf_x=[";
    sybuf=nm+"_surf_y=[";
    szbuf=nm+"_surf_z=[";
    for (i=0; i<nu; ++i)
    {
      for (j=0; j<nv; ++j)
      {
        sxbuf+=std::to_string(ps.f(u[i], v[j]).x());
        sybuf+=std::to_string(ps.f(u[i], v[j]).y());
        szbuf+=std::to_string(ps.f(u[i], v[j]).z());
        if (j==(nv-1))
        {
          if (i<(nu-1))
          {
            sxbuf+=";\n";
            sybuf+=";\n";
            szbuf+=";\n";
          }
        }
        else
        {
          sxbuf+=", ";
          sybuf+=", ";
          szbuf+=", ";
        }
      }
    }
    sxbuf+="];";
    sybuf+="];";
    szbuf+="];";

    std::cout << "% surface: " << nm << std::endl;
    std::cout << "figure(" << figno << ");" << std::endl;
    std::cout << sxbuf << std::endl;
    std::cout << sybuf << std::endl;
    std::cout << szbuf << std::endl;
    std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
    std::cout << "mesh(" << nm << "_surf_x, "
                         << nm << "_surf_y, "
                         << nm << "_surf_z, "
                         << "'EdgeColor', [0 0 0]);" << std::endl;
    if (show_control_points)
    {
      for (pp=0; pp<nup; ++pp)
      {
        for (qq=0; qq<nvp; ++qq)
        {
          index_type ppqq;

          ppqq=qq*nup+pp;
          std::cout << cpxbuf[ppqq] << std::endl;
          std::cout << cpybuf[ppqq] << std::endl;
          std::cout << cpzbuf[ppqq] << std::endl;

          std::cout << "mesh(" << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_x', "
                               << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_y', "
                               << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_z', "
                               << "'EdgeColor', [0.5 0.5 0.5], 'FaceColor', 'none');" << std::endl;
          std::cout << "plot3(" << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_x', "
                                << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_y', "
                                << nm << "_surf_cp" << std::to_string(pp) << std::to_string(qq) <<  "_z', "
                                << "'o', 'MarkerEdgeColor', [0.5 0.5 0.5],'MarkerFaceColor', [0.5 0.5 0.5]);" << std::endl;
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
      template<typename data__, unsigned short dim__, typename tol__>
      class general_skinning_surface_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;

          typedef connection_data<data_type, dim__, tolerance_type> rib_data_type;

          general_skinning_surface_creator()
            : piecewise_creator_base<data__, dim__, tol__>(0, 0), ribs(2),
              max_degree(1), closed(false)
          {
          }
          general_skinning_surface_creator(const data_type &uu0, const data_type &vv0)
            : piecewise_creator_base<data__, dim__, tol__>(uu0, vv0), ribs(2),
              max_degree(1), closed(false)
          {
          }
          general_skinning_surface_creator(const general_skinning_surface_creator<data_type, dim__, tolerance_type> & gs)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(gs), ribs(gs.ribs),
              max_degree(gs.max_degree), closed(gs.closed)
          {
          }
          virtual ~general_skinning_surface_creator()
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

            // split ribs so have same number of curves (with same joint parameters) for all ribs and get degree
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

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &ps) const
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

        private:
          std::vector<rib_data_type> ribs;
          std::vector<index_type> max_degree;
          bool closed;
      };
    }
  }
}
#endif
