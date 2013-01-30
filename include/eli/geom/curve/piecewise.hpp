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

#include "eli/geom/tolerance/simple.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__=geom::tolerance::simple<data__> >
      class piecewise
      {
        public:
          typedef curve__<data__, dim__, tol__> curve_type;
          typedef typename curve_type::index_type index_type;
          typedef typename curve_type::point_type point_type;
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

        private:
          template<template<typename, unsigned short, typename> class curve1__, typename data1__, unsigned short dim1__, typename tol1__>
          friend void length(typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &len, const piecewise<curve1__, data1__, dim1__, tol1__> &pc);
          template<template<typename, unsigned short, typename> class curve1__, typename data1__, unsigned short dim1__, typename tol1__>
          friend void length(typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &len, const piecewise<curve1__, data1__, dim1__, tol1__> &pc,
                             const typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &t0, const typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &t1);

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

        private:
          static bool check_continuity(const curve_type &curve1, const curve_type &curve2, const eli::geom::general::continuity &cont)
          {
            tolerance_type tol;

            switch(cont)
            {
              case(geom::general::C2):
              {
                if (!tol(curve1.fpp(1), curve2.fpp(0)))
                  return false;
              }
              case(geom::general::C1):
              {
                if (!tol(curve1.fp(1), curve2.fp(0)))
                  return false;
              }
              case(geom::general::C0):
              {
                return tol(curve1.f(1), curve2.f(0));
                break;
              }
              case(geom::general::NOT_CONNECTED):
              {
                return !tol(curve1.f(1), curve2.f(0));
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

        public:
          piecewise() : t0(0) {}
          piecewise(const piecewise<curve__, data_type, dim__, tol__> &p) : segments(p.segments), t0(p.t0) {}
          ~piecewise() {}

          bool operator==(const piecewise<curve__, data_type, dim__> &p) const
          {
            if (this==&p)
              return true;
            if (t0!=p.t0)
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

          size_t number_segments() const {return segments.size();}

          bool open() const
          {
            return check_continuity(segments.rbegin()->c, segments.begin()->c, eli::geom::general::C0);
          }
          bool closed() const
          {
            return !open();
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
              if (!check_continuity(curve, segments.begin()->c, eli::geom::general::C0))
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
              if (!check_continuity(segments.rbegin()->c, curve, eli::geom::general::C0))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            si.c=curve;
            si.delta_t=dt;
            segments.push_back(si);

            return NO_ERROR;
          }

          error_code get(curve_type &curve, data_type &dt, const index_type &index) const
          {
            if (static_cast<size_t>(index)>=number_segments())
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
            if (static_cast<size_t>(index)>=number_segments())
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
              if (!check_continuity(scito->c, curve, eli::geom::general::C0))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (static_cast<size_t>(index+1)<number_segments())
            {
              scito=scit;
              ++scito;
              if (!check_continuity(curve, scito->c, eli::geom::general::C0))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // set the new curve and delta t
            scit->c=curve;

            return NO_ERROR;
          }

          error_code replace(const curve_type &curve, const data_type &dt, const index_type &index)
          {
            if (static_cast<size_t>(index)>=number_segments())
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
              if (!check_continuity(scito->c, curve, eli::geom::general::C0))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (static_cast<size_t>(index+1)<number_segments())
            {
              scito=scit;
              ++scito;
              if (!check_continuity(curve, scito->c, eli::geom::general::C0))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // set the new curve and delta t
            scit->c=curve;
            scit->delta_t=dt;

            return NO_ERROR;
          }

          error_code replace(const curve_type &curve, const data_type &dt, const index_type &index0, const index_type &index1)
          {
            if (static_cast<size_t>(index0)>=number_segments())
              return INVALID_INDEX;
            if (static_cast<size_t>(index1)>=number_segments())
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
              if (!check_continuity(scito->c, curve, eli::geom::general::C0))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (static_cast<size_t>(index1)<number_segments())
            {
              scito=scit1;
              if (!check_continuity(curve, scito->c, eli::geom::general::C0))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // set the new curve and delta t
            scit0->c=curve;
            scit0->delta_t=dt;
            ++scit0;
            segments.erase(scit0, scit1);

            return NO_ERROR;
          }

          error_code replace(const piecewise<curve__, data_type, dim__> &p, const index_type &index)
          {
            if (static_cast<size_t>(index)>=number_segments())
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
              if (!check_continuity(scito->c, cs, eli::geom::general::C0))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (static_cast<size_t>(index+1)<number_segments())
            {
              scito=scit;
              ++scito;
              if (!check_continuity(ce, scito->c, eli::geom::general::C0))
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

            return NO_ERROR;
          }

          error_code replace(const piecewise<curve__, data_type, dim__> &p, const index_type &index0, const index_type &index1)
          {
            if (static_cast<size_t>(index0)>=number_segments())
              return INVALID_INDEX;
            if (static_cast<size_t>(index1)>=number_segments())
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
              if (!check_continuity(scito->c, cs, eli::geom::general::C0))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (static_cast<size_t>(index1)<number_segments())
            {
              scito=scit1;
              if (!check_continuity(ce, scito->c, eli::geom::general::C0))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // erase old segments
            segments.erase(scit0, scit1);

            // cycle through piecewise to get each segment and insert
            segments.insert(scit1, p.segments.begin(), p.segments.end());

            return NO_ERROR;
          }

          // TODO: NEED TO IMPLEMENT
          //       * fit
          //       * interpolate

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

            return NO_ERROR;
          }

          point_type f(const data_type &t) const
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
      };
    }
  }
}
#endif
