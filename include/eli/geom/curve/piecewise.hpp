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

#ifndef eli_geom_curve_piecewise_hpp
#define eli_geom_curve_piecewise_hpp

#include <list>

#include "eli/util/tolerance.hpp"

#include "eli/geom/general/continuity.hpp"

namespace eli
{
  namespace geom
  {
    namespace utility
    {
      // NOTE: Create check_equivalence function that demotes higher order curve to see
      //       if get (nearly) the same control points for the two curves. This could
      //       work for all curve types that have control points and support demotion.

      template<typename curve1__, typename curve2__, typename tol__>
      static bool check_joint_continuity(const curve1__ &curve1, const curve2__ &curve2, const eli::geom::general::continuity &cont, const tol__ &tol)
      {
        switch(cont)
        {
          case(eli::geom::general::G2):
          {
          }
          case(eli::geom::general::C2):
          {
            if (!tol.approximately_equal(curve1.fpp(1), curve2.fpp(0)))
              return false;
          }
          case(eli::geom::general::C1):
          {
            if (!tol.approximately_equal(curve1.fp(1), curve2.fp(0)))
              return false;
          }
          case(eli::geom::general::C0):
          {
            return tol.approximately_equal(curve1.f(1), curve2.f(0));
            break;
          }
          case(eli::geom::general::NOT_CONNECTED):
          {
            return !tol.approximately_equal(curve1.f(1), curve2.f(0));
            break;
          }
          default:
          {
            return false;
            break;
          }
        }

        return false;
      }
    }
  }
}

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__=eli::util::tolerance<data__> >
      class piecewise
      {
        public:
          typedef curve__<data__, dim__, tol__> curve_type;
          typedef typename curve_type::index_type index_type;
          typedef typename curve_type::point_type point_type;
          typedef typename curve_type::control_point_type control_point_type;
          typedef typename curve_type::rotation_matrix_type rotation_matrix_type;
          typedef data__ data_type;
          typedef unsigned short dimension_type;
          typedef tol__ tolerance_type;
          enum error_code
          {
            NO_ERROR=0,
            INVALID_INDEX=1,
            INDEX_NOT_FOUND=2,
            INVALID_PARAM=50,
            INVALID_PARAM_DIFFERENCE=51,
            SEGMENT_NOT_CONNECTED=100,
            UNKNOWN_ERROR=999
          };

        public:
          piecewise() : t0(0) {}
          piecewise(const piecewise<curve__, data_type, dim__, tol__> &p) : segments(p.segments), t0(p.t0), tol(p.tol) {}
          ~piecewise() {}

          bool operator==(const piecewise<curve__, data_type, dim__> &p) const
          {
            if (this==&p)
              return true;
            if (t0!=p.t0)
              return false;
            if (tol!=p.tol)
              return false;
            if (number_segments()!=p.number_segments())
              return false;
            typename segment_collection_type::const_iterator scit, it;
            for (scit=segments.begin(), it=p.segments.begin(); scit!=segments.end(); ++scit, ++it)
            {
              if ((*it)!=(*scit))
                return false;
            }

            return true;
          }

          bool operator!=(const piecewise<curve__, data_type, dim__> &p) const
          {
            return !operator==(p);
          }

          static dimension_type dimension() {return dim__;}

          const data_type & get_t0() const {return t0;}
          void set_t0(const data_type &t0_in) {t0=t0_in;}

          void get_parameter_min(data_type &tmin) const
          {
            tmin=t0;
          }

          void get_parameter_max(data_type &tmax) const
          {
            typename segment_collection_type::const_iterator psit;

            tmax=t0;

            for (psit=segments.begin(); psit!=segments.end(); ++psit)
            {
              tmax+=psit->delta_t;
            }
          }

          index_type number_segments() const {return static_cast<index_type>(segments.size());}

          void get_bounding_box(point_type &pmin, point_type &pmax) const
          {
            typename segment_collection_type::const_iterator it=segments.begin();
            point_type pmintmp, pmaxtmp;

            // cycle through all segments to get each bounding box to compare
            it->c.get_bounding_box(pmin, pmax);
            for (++it; it!=segments.end(); ++it)
            {
              it->c.get_bounding_box(pmintmp, pmaxtmp);
              for (index_type i=0; i<dim__; ++i)
              {
                if (pmintmp(i)<pmin(i))
                {
                  pmin(i)=pmintmp(i);
                }
                if (pmaxtmp(i)>pmax(i))
                {
                  pmax(i)=pmaxtmp(i);
                }
              }
            }
          }

          void rotate(const rotation_matrix_type &rmat)
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->c.rotate(rmat);
            }
          }

          void rotate(const rotation_matrix_type &rmat, const point_type &rorig)
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->c.rotate(rmat, rorig);
            }
          }

          void translate(const point_type &trans)
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->c.translate(trans);
            }
          }

          bool closed() const
          {
            return eli::geom::utility::check_joint_continuity(segments.rbegin()->c, segments.begin()->c, eli::geom::general::C0, tol);
          }
          bool open() const
          {
            return !closed();
          }

          bool continuous(const eli::geom::general::continuity &cont)
          {
            switch (cont)
            {
              case(eli::geom::general::NOT_CONNECTED):
              {
                return check_continuity(cont);
              }
              default:
              {
                assert(false);
                return false;
              }
            }

            return true;
          }

          void reverse()
          {
            // reverse order of segments
            segments.reverse();

            // reverse each segment
            for (typename segment_collection_type::iterator it=segments.begin(); it!=segments.end(); ++it)
            {
              it->c.reverse();
            }

            // check if still connected
            assert(check_continuity(eli::geom::general::C0));
          }

          void clear() {segments.clear();}

          template<typename it__>
          error_code set(it__ itb, it__ ite)
          {
            segments.clear();

            segment_info si;
            index_type i;
            it__ it;
            for (i=0, it=itb; it!=ite; ++i, ++it)
            {
              error_code err=push_back(*it);
              if (err!=NO_ERROR)
              {
                segments.clear();
                return err;
              }
            }

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

          template<typename it__, typename itd__>
          error_code set(it__ itb, it__ ite, itd__ itd)
          {
            segments.clear();

            segment_info si;
            index_type i;
            it__ it;
            itd__ itdt;
            for (i=0, it=itb, itdt=itd; it!=ite; ++i, ++it, ++itdt)
            {
              error_code err=push_back(*it, *itdt);
              if (err!=NO_ERROR)
              {
                segments.clear();
                return err;
              }
            }

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

          error_code push_front(const curve_type &curve, const data_type &dt=1.0)
          {
            segment_info si;

            if (dt<=0)
              return INVALID_PARAM_DIFFERENCE;

            // check to make sure have valid segments
            if (!segments.empty())
            {
              if (!eli::geom::utility::check_joint_continuity(curve, segments.begin()->c, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // add segment
            si.c=curve;
            si.delta_t=dt;
            segments.push_front(si);

            // adjust start so that existing parameter values match for existing segments
            t0-=dt;

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

          error_code push_back(const curve_type &curve, const data_type &dt=1.0)
          {
            segment_info si;

            if (dt<=0)
              return INVALID_PARAM_DIFFERENCE;

            // check to make sure have valid segments
            if (!segments.empty())
            {
              if (!eli::geom::utility::check_joint_continuity(segments.rbegin()->c, curve, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            si.c=curve;
            si.delta_t=dt;
            segments.push_back(si);

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

          error_code get(curve_type &curve, const index_type &index) const
          {
            data_type dt;
            return get(curve, dt, index);
          }

          error_code get(curve_type &curve, data_type &dt, const index_type &index) const
          {
            if (index>=number_segments())
              return INVALID_INDEX;

            // advance to desired index
            index_type i;
            typename segment_collection_type::const_iterator scit;
            for (i=0, scit=segments.begin(); i<index; ++i, ++scit) {}

            curve=scit->c;
            dt=scit->delta_t;
            return NO_ERROR;
          }

          error_code replace(const curve_type &curve, const index_type &index)
          {
            if (index>=number_segments())
              return INVALID_INDEX;

            // advance to desired index
            index_type i;
            typename segment_collection_type::iterator scit, scito;
            for (i=0, scit=segments.begin(); i<index; ++i, ++scit) {}

            // check the connectivity on adjacent nodes (if available)
            if (index>0)
            {
              scito=scit;
              --scito;
              if (!eli::geom::utility::check_joint_continuity(scito->c, curve, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if ((index+1)<number_segments())
            {
              scito=scit;
              ++scito;
              if (!eli::geom::utility::check_joint_continuity(curve, scito->c, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // set the new curve and delta t
            scit->c=curve;

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

          error_code replace(const curve_type &curve, const data_type &dt, const index_type &index)
          {
            if (index>=number_segments())
              return INVALID_INDEX;

            if (dt<=0)
              return INVALID_PARAM_DIFFERENCE;

            // advance to desired index
            index_type i;
            typename segment_collection_type::iterator scit, scito;
            for (i=0, scit=segments.begin(); i<index; ++i, ++scit) {}

            // check the connectivity on adjacent nodes (if available)
            if (index>0)
            {
              scito=scit;
              --scito;
              if (!eli::geom::utility::check_joint_continuity(scito->c, curve, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if ((index+1)<number_segments())
            {
              scito=scit;
              ++scito;
              if (!eli::geom::utility::check_joint_continuity(curve, scito->c, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // set the new curve and delta t
            scit->c=curve;
            scit->delta_t=dt;

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

          error_code replace(const curve_type &curve, const data_type &dt, const index_type &index0, const index_type &index1)
          {
            if (index0>=number_segments())
              return INVALID_INDEX;
            if (index1>=number_segments())
              return INVALID_INDEX;
            if (index0>=index1)
              return INVALID_INDEX;

            if (dt<=0)
              return INVALID_PARAM_DIFFERENCE;

            // advance to desired index
            index_type i;
            typename segment_collection_type::iterator scit0, scit1, scito;
            for (i=0, scit0=segments.begin(); i<index0; ++i, ++scit0) {}
            for (scit1=scit0; i<index1; ++i, ++scit1) {}

            // check the connectivity on adjacent nodes (if available)
            if (index0>0)
            {
              scito=scit0;
              --scito;
              if (!eli::geom::utility::check_joint_continuity(scito->c, curve, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (index1<number_segments())
            {
              scito=scit1;
              if (!eli::geom::utility::check_joint_continuity(curve, scito->c, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // set the new curve and delta t
            scit0->c=curve;
            scit0->delta_t=dt;
            ++scit0;
            segments.erase(scit0, scit1);

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

          error_code replace(const piecewise<curve__, data_type, dim__> &p, const index_type &index)
          {
            if (index>=number_segments())
              return INVALID_INDEX;

            // advance to desired index
            index_type i;
            typename segment_collection_type::iterator scit, scito;
            for (i=0, scit=segments.begin(); i<index; ++i, ++scit) {}

            // get the first and last curve
            curve_type cs, ce;
            data_type dt;
            p.get(cs, dt, 0);
            p.get(ce, dt, p.number_segments()-1);

            // check the connectivity on adjacent nodes (if available)
            if (index>0)
            {
              scito=scit;
              --scito;
              if (!eli::geom::utility::check_joint_continuity(scito->c, cs, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if ((index+1)<number_segments())
            {
              scito=scit;
              ++scito;
              if (!eli::geom::utility::check_joint_continuity(ce, scito->c, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // replace element
            typename segment_collection_type::const_iterator it=p.segments.begin();
            (*scit)=(*it);

            // insert the rest
            ++scit; ++it;
            segments.insert(scit, it, p.segments.end());

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

          error_code replace(const piecewise<curve__, data_type, dim__> &p, const index_type &index0, const index_type &index1)
          {
            if (index0>=number_segments())
              return INVALID_INDEX;
            if (index1>=number_segments())
              return INVALID_INDEX;
            if (index0>=index1)
              return INVALID_INDEX;

            // advance to desired index
            index_type i;
            typename segment_collection_type::iterator scit0, scit1, scito;
            for (i=0, scit0=segments.begin(); i<index0; ++i, ++scit0) {}
            for (scit1=scit0; i<index1; ++i, ++scit1) {}

            // get the first and last curve
            curve_type cs, ce;
            data_type dt;
            p.get(cs, dt, 0);
            p.get(ce, dt, p.number_segments()-1);

            // check the connectivity on adjacent nodes (if available)
            if (index0>0)
            {
              scito=scit0;
              --scito;
              if (!eli::geom::utility::check_joint_continuity(scito->c, cs, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (index1<number_segments())
            {
              scito=scit1;
              if (!eli::geom::utility::check_joint_continuity(ce, scito->c, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // erase old segments
            segments.erase(scit0, scit1);

            // cycle through piecewise to get each segment and insert
            segments.insert(scit1, p.segments.begin(), p.segments.end());

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

          error_code split(const data_type &t)
          {
            // find segment that corresponds to given t
            typename segment_collection_type::iterator it;
            data_type tt;
            find_segment(it, tt, t);

            if (it==segments.end())
              return INVALID_PARAM;

            // split the segment and replace
            curve_type cl, cr;
            segment_info stl, str;
            it->c.split(cl, cr, tt);
            stl.c=cl;
            stl.delta_t=it->delta_t*tt;
            str.c=cr;
            str.delta_t=it->delta_t*(1-tt);
            (*it)=str;
            segments.insert(it, stl);

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

          point_type f(const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it;
            data_type tt(0);
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              --it;
            }

            return it->c.f(tt);
          }

          point_type fp(const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it;
            data_type tt;
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              --it;
            }

            return it->c.fp(tt)/it->delta_t;
          }

          point_type fpp(const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it;
            data_type tt;
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              --it;
            }

            return it->c.fpp(tt)/(it->delta_t*it->delta_t);
          }

          point_type fppp(const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it;
            data_type tt;
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              --it;
            }

            return it->c.fppp(tt)/(it->delta_t*it->delta_t*it->delta_t);
          }

          point_type tanget(const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it;
            data_type tt(0);
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              --it;
            }

            return it->c.tangent(tt);
          }

          void frenet_serret_frame(point_type &t, point_type &n, point_type &b, const data_type &t0)
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it;
            data_type tt(0);
            find_segment(it, tt, t0);

            if (it==segments.end())
            {
              assert(false);
              --it;
            }

            it->c.frenet_serret_frame(t, n, b, tt);
          }

          // TODO: NEED TO IMPLEMENT
          //       * fit
          //       * interpolate

        private:
          template<template<typename, unsigned short, typename> class curve1__,
                   typename data1__, unsigned short dim1__, typename tol1__>
          friend void length(typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &len,
                             const piecewise<curve1__, data1__, dim1__, tol1__> &pc,
                             const typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &tol);
          template<template<typename, unsigned short, typename> class curve1__,
                            typename data1__, unsigned short dim1__, typename tol1__>
          friend void length(typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &len,
                             const piecewise<curve1__, data1__, dim1__, tol1__> &pc,
                             const typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &t0,
                             const typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &t1,
                             const typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &tol);

          struct segment_info
          {
            curve_type c;
            data_type delta_t;

            segment_info() : delta_t(1) {}
            segment_info(const segment_info &si) : c(si.c), delta_t(si.delta_t) {}
            ~segment_info() {}

            bool operator==(const segment_info &si) const
            {
              if (this==&si)
                return true;
              if (delta_t!=si.delta_t)
                return false;
              if (c!=si.c)
                return false;

              return true;
            }

            bool operator!=(const segment_info &si) const
            {
              return !operator==(si);
            }
          };
          typedef std::list<segment_info> segment_collection_type;

          segment_collection_type segments;
          data_type t0;
          tolerance_type tol;

        private:
          bool check_continuity(const eli::geom::general::continuity &cont) const
          {
            typename segment_collection_type::const_iterator it(segments.begin()), itp(it);

            for (++it; it!=segments.end(); ++it, ++itp)
            {
              if (!eli::geom::utility::check_joint_continuity(itp->c, it->c, cont, tol))
              {
                return false;
              }
            }

            return true;
          }

          void find_segment(typename segment_collection_type::const_iterator &it, data_type &tt, const data_type &t_in) const
          {
            data_type t(t0);

            // check to see if have invalid t_in
            if (t_in<t0)
            {
              it=segments.end();
              return;
            }

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              if (t_in<=t+it->delta_t)
              {
                tt=(t_in-t)/it->delta_t;
                return;
              }
              t+=it->delta_t;
            }
          }

          void find_segment(typename segment_collection_type::iterator &it, data_type &tt, const data_type &t_in)
          {
            data_type t(t0);

            // check to see if have invalid t_in
            if (t_in<t0)
            {
              it=segments.end();
              return;
            }

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              if (t_in<=t+it->delta_t)
              {
                tt=(t_in-t)/it->delta_t;
                return;
              }
              t+=it->delta_t;
            }
          }
      };
    }
  }
}
#endif
