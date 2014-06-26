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

#ifndef eli_geom_surface_piecewise_connection_data_hpp
#define eli_geom_surface_piecewise_connection_data_hpp

#include "eli/util/tolerance.hpp"

#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise.hpp"

#include <iterator>

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class connection_data
      {
        public:
          enum
          {
            CONNECTION_SET=0x000001, // must be set for valid joint
            LEFT_FP_SET   =0x000010,
            RIGHT_FP_SET  =0x000100,
            LEFT_FPP_SET  =0x001000,
            RIGHT_FPP_SET =0x010000
          };

          enum connection_continuity
          {
            C0=general::C0,
            C1=general::C1,
            C2=general::C2
          };

          typedef data__ data_type;
          typedef tol__ tolerance_type;
          typedef typename eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, dim__, tol__> curve_type;
          typedef typename curve_type::index_type index_type;

        public:
          connection_data() : conditions(0), continuity(C0)
          {
          }

          connection_data(const connection_data &cd)
            : f(cd.f), fp_left(cd.fp_left), fp_right(cd.fp_right), fpp_left(cd.fpp_left),
              fpp_right(cd.fpp_right), conditions(cd.conditions), continuity(cd.continuity)
          {
          }

          ~connection_data() {}

          const connection_data & operator=(const connection_data &cd)
          {
            if (this!=&cd)
            {
              f=cd.f;
              fp_left=cd.fp_left;
              fp_right=cd.fp_right;
              fpp_left=cd.fpp_left;
              fpp_right=cd.fpp_right;
              conditions=cd.conditions;
              continuity=cd.continuity;
            }

            return (*this);
          }

          bool operator==(const connection_data &cd) const
          {
            tolerance_type tol;

            if (conditions!=cd.conditions)
              return false;
            if (continuity!=cd.continuity)
              return false;
            if (cd.use_f())
            {
              if (f!=cd.f)
                return false;
            }
            if (cd.use_left_fp())
            {
              if (fp_left!=cd.fp_left)
                return false;
            }
            if (cd.use_right_fp())
            {
              if (fp_right!=cd.fp_right)
                return false;
            }
            if (cd.use_left_fpp())
            {
              if (fpp_left!=cd.fpp_left)
                return false;
            }
            if (cd.use_right_fpp())
            {
              if (fpp_right!=cd.fpp_right)
                return false;
            }

            return true;
          }

          bool operator!=(const connection_data &cd) const {return !operator==(cd);}

          data_type get_t0() const
          {
            return f.get_t0();
          }
          data_type get_tmax() const
          {
            return f.get_tmax();
          }

          // get joints for all active curves
          template<typename output_it__>
          void get_joints(output_it__ it_in) const
          {
            tolerance_type tol;
            auto comp = [&tol](const data_type &x1, const data_type &x2)->bool
            {
              return tol.approximately_less_than(x1, x2);
            };

            std::vector<data_type> jts1, jts, jts_out;
            if (use_f())
            {
              f.get_parameters(std::back_inserter(jts1));
              if (jts.empty())
              {
                std::swap(jts, jts1);
              }
              else
              {
                std::set_union(jts.begin(), jts.end(), jts1.begin(), jts1.end(), std::back_inserter(jts_out), comp);
                std::swap(jts, jts_out);
                jts1.clear();
                jts_out.clear();
              }
            }
            if (use_left_fp())
            {
              fp_left.get_parameters(std::back_inserter(jts1));
              if (jts.empty())
              {
                std::swap(jts, jts1);
              }
              else
              {
                std::set_union(jts.begin(), jts.end(), jts1.begin(), jts1.end(), std::back_inserter(jts_out), comp);
                std::swap(jts, jts_out);
                jts1.clear();
                jts_out.clear();
              }
            }
            if (use_right_fp())
            {
              fp_right.get_parameters(std::back_inserter(jts1));
              if (jts.empty())
              {
                std::swap(jts, jts1);
              }
              else
              {
                std::set_union(jts.begin(), jts.end(), jts1.begin(), jts1.end(), std::back_inserter(jts_out), comp);
                std::swap(jts, jts_out);
                jts1.clear();
                jts_out.clear();
              }
            }
            if (use_left_fpp())
            {
              fpp_left.get_parameters(std::back_inserter(jts1));
              if (jts.empty())
              {
                std::swap(jts, jts1);
              }
              else
              {
                std::set_union(jts.begin(), jts.end(), jts1.begin(), jts1.end(), std::back_inserter(jts_out), comp);
                std::swap(jts, jts_out);
                jts1.clear();
                jts_out.clear();
              }
            }
            if (use_right_fpp())
            {
              fpp_right.get_parameters(std::back_inserter(jts1));
              if (jts.empty())
              {
                std::swap(jts, jts1);
              }
              else
              {
                std::set_union(jts.begin(), jts.end(), jts1.begin(), jts1.end(), std::back_inserter(jts_out), comp);
                std::swap(jts, jts_out);
                jts1.clear();
                jts_out.clear();
              }
            }

            std::copy(jts.begin(), jts.end(), it_in);
          }

          // split all curves
          template<typename it1__, typename it2__>
          void split(it1__ itb, it1__ ite, it2__ itd)
          {
            index_type i, njoints(0);
            for (it1__ it=itb; it!=ite; ++it, ++njoints)
            {
              if (use_f())
              {
                f.split(*it);
              }
              if (use_left_fp())
              {
                fp_left.split(*it);
              }
              if (use_right_fp())
              {
                fp_right.split(*it);
              }
              if (use_left_fpp())
              {
                fpp_left.split(*it);
              }
              if (use_right_fpp())
              {
                fpp_right.split(*it);
              }
            }

            index_type nsegs(njoints-1);
            std::vector<index_type> deg(nsegs), max_deg(nsegs, 0);
            if (use_f())
            {
              f.degrees(deg.begin());
              for (i=0; i<nsegs; ++i)
              {
                if (deg[i]>max_deg[i])
                {
                  max_deg[i]=deg[i];
                }
              }
            }
            if (use_left_fp())
            {
              fp_left.degrees(deg.begin());
              for (i=0; i<nsegs; ++i)
              {
                if (deg[i]>max_deg[i])
                {
                  max_deg[i]=deg[i];
                }
              }
            }
            if (use_right_fp())
            {
              fp_right.degrees(deg.begin());
              for (i=0; i<nsegs; ++i)
              {
                if (deg[i]>max_deg[i])
                {
                  max_deg[i]=deg[i];
                }
              }
            }
            if (use_left_fpp())
            {
              fpp_left.degrees(deg.begin());
              for (i=0; i<nsegs; ++i)
              {
                if (deg[i]>max_deg[i])
                {
                  max_deg[i]=deg[i];
                }
              }
            }
            if (use_right_fpp())
            {
              fpp_right.degrees(deg.begin());
              for (i=0; i<nsegs; ++i)
              {
                if (deg[i]>max_deg[i])
                {
                  max_deg[i]=deg[i];
                }
              }
            }

            for (i=0; i<nsegs; ++i, ++itd)
            {
              (*itd)=max_deg[i];
            }
          }

          template<typename it__>
          bool promote(it__ itb, it__ ite)
          {
            // NOTE: the piecewise curve class will return error if degree vector is
            //       larger than number of segments
            typename curve_type::error_code ec;
            it__ it;
            index_type i;
            for (i=0, it=itb; it!=ite; ++it, ++i)
            {
              if (use_f())
              {
                ec=f.degree_promote_to(i, (*it));
                if (ec!=curve_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              if (use_left_fp())
              {
                ec=fp_left.degree_promote_to(i, (*it));
                if (ec!=curve_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              if (use_right_fp())
              {
                ec=fp_right.degree_promote_to(i, (*it));
                if (ec!=curve_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              if (use_left_fpp())
              {
                ec=fpp_left.degree_promote_to(i, (*it));
                if (ec!=curve_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              if (use_right_fpp())
              {
                ec=fpp_right.degree_promote_to(i, (*it));
                if (ec!=curve_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
            }

            return true;
          }
          
          // connection interface
          bool set_f(const curve_type &ff)
          {
            f=ff;
            conditions|=CONNECTION_SET;
            return check_state();
          }
          const curve_type & get_f() const {return f;}
          bool unset_f()
          {
            conditions&=~CONNECTION_SET;
            return check_state();
          }
          bool use_f() const
          {
            return ((conditions & CONNECTION_SET) == CONNECTION_SET);
          }

          // first derivative interface
          bool set_left_fp(const curve_type &fpl)
          {
            fp_left=fpl;
            conditions|=LEFT_FP_SET;

            // if already wanting C1 continuous then set right as well
            if (get_continuity()>C0)
            {
              fp_right=fpl;
              conditions|=RIGHT_FP_SET;
            }

            return check_state();
          }
          bool set_right_fp(const curve_type &fpr)
          {
            fp_right=fpr;
            conditions|=RIGHT_FP_SET;

            // if already wanting C1 continuous then set left as well
            if (get_continuity()>C0)
            {
              fp_left=fpr;
              conditions|=LEFT_FP_SET;
            }

            return check_state();
          }
          bool set_fp(const curve_type &p)
          {
            conditions|=(LEFT_FP_SET | RIGHT_FP_SET);
            fp_left=p;
            fp_right=p;
            return check_state();
          }

          const curve_type & get_left_fp() const
          {
            return fp_left;
          }
          const curve_type & get_right_fp() const
          {
            return fp_right;
          }
          void get_fp(curve_type &fpl, curve_type &fpr) const
          {
            fpl=fp_left;
            fpr=fp_right;
          }

          bool unset_fp()
          {
            conditions&=~(LEFT_FP_SET | RIGHT_FP_SET);
            return check_state();
          }
          bool unset_left_fp()
          {
            conditions&=~LEFT_FP_SET;

            return check_state();
          }
          bool unset_right_fp()
          {
            conditions&=~RIGHT_FP_SET;

            return check_state();
          }

          bool use_left_fp() const
          {
            return ((conditions & LEFT_FP_SET) == LEFT_FP_SET);
          }
          bool use_right_fp() const
          {
            return ((conditions & RIGHT_FP_SET) == RIGHT_FP_SET);
          }

          // second derivative interface
          bool set_left_fpp(const curve_type &fppl)
          {
            fpp_left=fppl;
            conditions|=LEFT_FPP_SET;

            // if already wanting C2 continuous then set right as well
            if (get_continuity()>C1)
            {
              fpp_right=fppl;
              conditions|=RIGHT_FPP_SET;
            }

            return check_state();
          }
          bool set_right_fpp(const curve_type &fppr)
          {
            fpp_right=fppr;
            conditions|=RIGHT_FPP_SET;

            // if already wanting C2 continuous then set left as well
            if (get_continuity()>C1)
            {
              fpp_left=fppr;
              conditions|=LEFT_FPP_SET;
            }

            return check_state();
          }
          bool set_fpp(const curve_type &p)
          {
            conditions|=(LEFT_FPP_SET | RIGHT_FPP_SET);
            fpp_left=p;
            fpp_right=p;
            return check_state();
          }

          const curve_type & get_left_fpp() const
          {
            return fpp_left;
          }
          const curve_type & get_right_fpp() const
          {
            return fpp_right;
          }
          void get_fpp(curve_type &fppl, curve_type &fppr) const
          {
            fppl=fpp_left;
            fppr=fpp_right;
          }

          bool unset_fpp()
          {
            conditions&=~(LEFT_FPP_SET | RIGHT_FPP_SET);
            return check_state();
          }
          bool unset_left_fpp()
          {
            conditions&=~LEFT_FPP_SET;

            return check_state();
          }
          bool unset_right_fpp()
          {
            conditions&=~RIGHT_FPP_SET;

            return check_state();
          }

          bool use_left_fpp() const
          {
            return ((conditions & LEFT_FPP_SET) == LEFT_FPP_SET);
          }
          bool use_right_fpp() const
          {
            return ((conditions & RIGHT_FPP_SET) == RIGHT_FPP_SET);
          }

          bool set_continuity(connection_continuity jc)
          {
            continuity=jc;
            return check_state();
          }
          connection_continuity get_continuity() const
          {
            return continuity;
          }

          bool check_state() const
          {
            tolerance_type tol;

            // check if point set
            if ((conditions & CONNECTION_SET)==0)
              return false;

            // if highest continuity is C0 then done
            if (continuity==C0)
              return true;

            // check first derivatives match on both sides
            if ((conditions & (LEFT_FP_SET | RIGHT_FP_SET)) == (LEFT_FP_SET | RIGHT_FP_SET))
            {
              // TODO: make this comparison an apprimately_equal() call
//              if (!tol.approximately_equal(fp_left, fp_right))
              if (fp_left!=fp_right)
                return false;
            }
            // check neither side first derivatives set set
            else if ((conditions & (LEFT_FP_SET | RIGHT_FP_SET)) != 0)
            {
              return false;
            }

            // check that the parameterization of fp is same as f
            data_type t0(get_t0()), tmax(get_tmax());
            if (use_left_fp())
            {
              if (!tol.approximately_equal(t0, fp_left.get_t0()) || !tol.approximately_equal(tmax, fp_left.get_tmax()))
                return false;
            }
            if (use_right_fp())
            {
              if (!tol.approximately_equal(t0, fp_right.get_t0()) || !tol.approximately_equal(tmax, fp_right.get_tmax()))
                return false;
            }

            // if highest continuity is C1 then done
            if (continuity==C1)
              return true;

            // check second derivatives match on both sides
            if ((conditions & (LEFT_FPP_SET | RIGHT_FPP_SET)) == (LEFT_FPP_SET | RIGHT_FPP_SET))
            {
              // TODO: make this comparison an apprimately_equal() call
//              if (!tol.approximately_equal(fpp_left, fpp_right))
              if (fpp_left!=fpp_right)
                return false;
            }
            // check neither side second derivatives set set
            else if ((conditions & (LEFT_FPP_SET | RIGHT_FPP_SET)) != 0)
            {
              return false;
            }

            // check that the parameterization of fpp is same as f
            if (use_left_fpp())
            {
              if (!tol.approximately_equal(t0, fpp_left.get_t0()) || !tol.approximately_equal(tmax, fpp_left.get_tmax()))
                return false;
            }
            if (use_right_fpp())
            {
              if (!tol.approximately_equal(t0, fpp_right.get_t0()) || !tol.approximately_equal(tmax, fpp_right.get_tmax()))
                return false;
            }

            // since highest continuity can only be C2 then done
            return true;
          }

        private:
          curve_type f;
          curve_type fp_left, fp_right;
          curve_type fpp_left, fpp_right;
          unsigned int conditions;
          connection_continuity continuity;
      };
    }
  }
}
#endif

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
#include "eli/geom/surface/piecewise_creator_base.hpp"

#include <string>

namespace eli
{
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

  void octave_start(int figno)
  {
    std::cout << "clf(" << figno << ", 'reset');" << std::endl;
    std::cout << "figure(" << figno << ");" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "grid on;" << std::endl;
  }

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
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;

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
          if ((pp<(ns-1)) || (pp==(ns-1)) && (i<bez_deg))
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
                            << nm << "_curv_cp_z', '-og', 'MarkerFaceColor', [0 1 0]);" << std::endl;
    }
  }

  template<typename data__>
  void octave_print(int figno, const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &pc,
                    const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &vec,
                    const std::string &name="")
  {
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
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
  void octave_print(int figno, const const eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> &ps,
                    const std::string &name="", bool show_control_points=true)
  {
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;
    typedef typename piecewise_surface_type::tolerance_type tolerance_type;

    std::string nm, cpxbuf, cpybuf, cpzbuf, sxbuf, sybuf, szbuf;

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
      assert(false); // Need to create separate control point surfaces for each patch
      cpxbuf=nm+"_surf_cp_x=[";
      cpybuf=nm+"_surf_cp_y=[";
      cpzbuf=nm+"_surf_cp_z=[";
      for (pp=0; pp<nup; ++pp)
      {
        for (qq=0; qq<nvp; ++qq)
        {
          surface_type bez;
          ps.get(bez, pp, qq);
          for (i=0; i<=bez.degree_u(); ++i)
          {
            for (j=0; j<=bez.degree_v(); ++j)
            {
              cpxbuf+=std::to_string(bez.get_control_point(i, j).x());
              cpybuf+=std::to_string(bez.get_control_point(i, j).y());
              cpzbuf+=std::to_string(bez.get_control_point(i, j).z());

              if ((pp<(nup-1)) || (qq<(nvp-1)))
              {
                cpxbuf+=", ";
                cpybuf+=", ";
                cpzbuf+=", ";
              }
              else
              {
                if ((i<bez.degree_u()) || (j<bez.degree_v()))
                {
                  cpxbuf+=", ";
                  cpybuf+=", ";
                  cpzbuf+=", ";
                }
              }
            }
          }
        }
      }
      cpxbuf+="];";
      cpybuf+="];";
      cpzbuf+="];";
    }

    // initialize the u & v parameters
    index_type nu(21), nv(21);
    std::vector<data__> u(nu), v(nv);
    for (i=0; i<static_cast<index_type>(u.size()); ++i)
    {
      u[i]=umin+(umax-umin)*static_cast<data__>(i)/(u.size()-1);
    }
    for (j=0; j<static_cast<index_type>(v.size()); ++j)
    {
      v[j]=vmin+(vmax-vmin)*static_cast<data__>(j)/(v.size()-1);
    }

    // set the curve points
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
    if (show_control_points)
    {
      std::cout << cpxbuf << std::endl;
      std::cout << cpybuf << std::endl;
      std::cout << cpzbuf << std::endl;
    }
    std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
    std::cout << "mesh(" << nm << "_surf_x, "
                         << nm << "_surf_y, "
                         << nm << "_surf_z, "
                         << "zeros(size(surf_z)), 'EdgeColor', [0 0 0]);" << std::endl;
    if (show_control_points)
    {
      std::cout << "mesh(" << nm << "_surf_cp_x', "
                           << nm << "_surf_cp_y', "
                           << nm << "_surf_cp_z', 'EdgeColor', [0 0 0]);" << std::endl;
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
            data_type v_start(ribs[0].get_t0()), v_end(ribs[0].get_tmax());
            tolerance_type tol;
            for (i=0; i<nribs; ++i)
            {
              if (!rbs[i].check_state())
                return false;
              if (!tol.approximately_equal(ribs[i].get_t0(), v_start) || !tol.approximately_equal(ribs[i].get_tmax(), v_end))
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
            this->set_v0(joints[0]);
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
            typedef typename rib_data_type::curve_type curve_type;
            typedef typename rib_data_type::tolerance_type;

            index_type nribs(this->get_number_u_segments()+1), i, j;
            std::vector<index_type> seg_degree(nribs-1);
            std::vector<rib_data_type> rib_states(ribs);
            tolerance_type tol;

            // FIX: Should be able to handle closed surfaces
            assert(!closed);

            // FIX: Need way of handling maximum u-degree specifi0]);ation for each surface patch
            for (i=0; i<(nribs-1); ++i)
            {
              assert(max_degree[i]==0);
            }

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
                if (jdegs[i]>max_jdegs[i])
                {
                  max_jdegs[i]=jdegs[i];
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

                  ribs[u].get_f().get(jcrv, v);
                  joints[u].set_f(jcrv.get_control_point(j));
                  if (ribs[u].use_left_fp())
                  {
                    ribs[u].get_left_fp().get(jcrv, v);
                    joints[u].set_left_fp(jcrv.get_control_point(j));
                  }
                  if (ribs[u].use_right_fp())
                  {
                    ribs[u].get_right_fp().get(jcrv, v);
                    joints[u].set_right_fp(jcrv.get_control_point(j));
                  }
                  if (ribs[u].use_left_fpp())
                  {
                    ribs[u].get_left_fpp().get(jcrv, v);
                    joints[u].set_left_fpp(jcrv.get_control_point(j));
                  }
                  if (ribs[u].use_right_fpp())
                  {
                    ribs[u].get_right_fpp().get(jcrv, v);
                    joints[u].set_right_fpp(jcrv.get_control_point(j));
                  }
                }

                // set the conditions for the curve creator
                bool rtn_flag(gc.set_conditions(joints, max_degree, closed));
                if (!rtn_flag)
                {
                  assert(false);
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
                  assert(false);
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
#if 0
                                   (piecewise<bezier, data__, dim__, tol__> &ps,
                                   const std::vector<eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, dim__, tol__> > &pc,
                                   bool outward_normal)
        typedef piecewise<bezier, data__, dim__, tol__> piecewise_surface_type;
        typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, dim__, tol__> piecewise_curve_type;
        typedef eli::geom::curve::piecewise_circle_creator<data__, dim__, tol__> circle_creator_type;
        typedef typename piecewise_curve_type::curve_type curve_type;
        typedef typename piecewise_surface_type::surface_type surface_type;
        typedef typename curve_type::point_type point_type;
        typedef typename curve_type::data_type data_type;
        typedef typename curve_type::index_type index_type;


        // not implemented
        assert(false);

        // make sure each curve is parameterized the same way

        // for each strip of surfaces need to sub-divide edge curves so that all edge curves have
        // the same number of segments with joints at same parameter values

        // 
        return false;

        curve_type c, arc;
        surface_type s[4];
        piecewise_curve_type pc_circle;
        point_type origin, normal, start;
        data_type du, dv;
        index_type i, j, pp, qq, nu=pc.number_segments(), nv=4, udim, vdim;

        // resize the surface
        ps.init_uv(nu, nv);

        // set the axis of rotation
        switch(axis)
        {
          case(0):
          {
            normal << 1, 0, 0;
            break;
          }
          case(1):
          {
            normal << 0, 1, 0;
            break;
          }
          case(2):
          {
            normal << 0, 0, 1;
            break;
          }
          default:
          {
            assert(false);
            return false;
          }
        }
        if (outward_normal)
        {
          normal*=-1;
        }

        // cycle through each curve segment
        circle_creator_type circle_creator;
        for (pp=0; pp<nu; ++pp)
        {
          // resize the surface patch
          udim=c.dimension();
          vdim=3;
          pc.get(c, du, pp);
          s[0].resize(udim, vdim);
          s[1].resize(udim, vdim);
          s[2].resize(udim, vdim);
          s[3].resize(udim, vdim);

          // cycle through each control point in current curve and create patch control points
          for (i=0; i<=udim; ++i)
          {
            // set up the vectors for this circle
            start=c.get_control_point(i);
            origin=normal.dot(start)*normal;

            // get the circle
            circle_creator.set(start, origin, normal);
            if (!circle_creator.create(pc_circle))
            {
              assert(false);
              return false;
            }

            // for each segment of circle set the control points
            for (qq=0; qq<4; ++qq)
            {
              pc_circle.get(arc, dv, qq);
              for (j=0; j<4; ++j)
              {
                s[qq].set_control_point(arc.get_control_point(j), i, j);
              }
            }
          }

          // set the surface patches
          for (qq=0; qq<4; ++qq)
          {
            ps.set(s[qq], pp, qq);
          }
        }

        return true;
      }
#endif
    }
  }
}
#endif
