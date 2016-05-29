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

#include <map>
#include <iterator>

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/util/tolerance.hpp"

#include "eli/geom/general/continuity.hpp"

namespace eli
{
  namespace geom
  {
    namespace utility
    {
      template<typename curve1__, typename curve2__, typename tol__>
      bool check_point_continuity(const curve1__ &curve1, const typename curve1__::data_type &dt1,
                                  const curve2__ &curve2, const typename curve2__::data_type &dt2,
                                  const eli::geom::general::continuity &cont, const tol__ &tol)
      {
        switch(cont)
        {
          case(eli::geom::general::G3):
          {
            typename curve1__::point_type fppp1(curve1.fppp(1)); fppp1.normalize();
            typename curve2__::point_type fppp2(curve2.fppp(0)); fppp2.normalize();

            if (!tol.approximately_equal(fppp1, fppp2))
              return false;
            else
              return check_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::G2, tol);
            break;
          }
          case(eli::geom::general::C3):
          {
            if (!tol.approximately_equal(curve1.fppp(1)/dt1/dt1/dt1, curve2.fppp(0)/dt2/dt2/dt2))
              return false;
            else
              return check_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::C2, tol);
            break;
          }
          case(eli::geom::general::G2):
          {
            typename curve1__::point_type fpp1(curve1.fpp(1)); fpp1.normalize();
            typename curve2__::point_type fpp2(curve2.fpp(0)); fpp2.normalize();

            if (!tol.approximately_equal(fpp1, fpp2))
              return false;
            else
              return check_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::G1, tol);
            break;
          }
          case(eli::geom::general::C2):
          {
            if (!tol.approximately_equal(curve1.fpp(1)/dt1/dt1, curve2.fpp(0)/dt2/dt2))
              return false;
            else
              return check_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::C1, tol);
            break;
          }
          case(eli::geom::general::G1):
          {
            typename curve1__::point_type fp1(curve1.fp(1)); fp1.normalize();
            typename curve2__::point_type fp2(curve2.fp(0)); fp2.normalize();

            if (!tol.approximately_equal(fp1, fp2))
              return false;
            else
              return check_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::G0, tol);
            break;
          }
          case(eli::geom::general::C1):
          {
            if (!tol.approximately_equal(curve1.fp(1)/dt1, curve2.fp(0)/dt2))
              return false;
            else
              return check_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::C0, tol);
            break;
          }
          case(eli::geom::general::G0):
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
            // shouldn't get here
            assert(false);
            return false;
            break;
          }
        }

        // shouldn't get here
        assert(false);
        return false;
      }

      namespace internal
      {
        template<typename curve1__, typename curve2__, typename tol__>
        eli::geom::general::continuity report_point_continuity(const curve1__ &curve1, const typename curve1__::data_type &dt1,
                                                               const curve2__ &curve2, const typename curve2__::data_type &dt2,
                                                               const eli::geom::general::continuity &cont, const tol__ &tol)
        {
          typename curve1__::point_type v1;
          typename curve2__::point_type v2;

          switch(cont)
          {
            case(eli::geom::general::NOT_CONNECTED):
            {
              v1=curve1.f(1);
              v2=curve2.f(0);

              if (tol.approximately_equal(v1, v2))
                return report_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::C0, tol);
              else
                return cont;

              break;
            }
            case(eli::geom::general::C0):
            {
              v1=curve1.fp(1)/dt1;
              v2=curve2.fp(0)/dt2;

              if (tol.approximately_equal(v1, v2))
                return report_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::C1, tol);
            }
            case(eli::geom::general::G0):
            {
              v1=curve1.fp(1).normalized();
              v2=curve2.fp(0).normalized();

              if (tol.approximately_equal(v1, v2))
                return report_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::G1, tol);
              else
                return cont;

              break;
            }
            case(eli::geom::general::C1):
            {
              v1=curve1.fpp(1)/dt1/dt1;
              v2=curve2.fpp(0)/dt2/dt2;

              if (tol.approximately_equal(v1, v2))
                return report_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::C2, tol);
            }
            case(eli::geom::general::G1):
            {
              v1=curve1.fpp(1).normalized();
              v2=curve2.fpp(0).normalized();

              if (tol.approximately_equal(v1, v2))
                return report_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::G2, tol);
              else
                  return cont;

              break;
            }
            case(eli::geom::general::C2):
            {
              v1=curve1.fppp(1)/dt1/dt1/dt1;
              v2=curve2.fppp(0)/dt2/dt2/dt2;

              if (tol.approximately_equal(v1, v2))
                return report_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::C3, tol);
            }
            case(eli::geom::general::G2):
            {
              v1=curve1.fppp(1).normalized();
              v2=curve2.fppp(0).normalized();

              if (tol.approximately_equal(v1, v2))
                return report_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::G3, tol);
              else
                  return cont;

              break;
            }
            case(eli::geom::general::C3):
            case(eli::geom::general::G3):
            {
              return cont;
            }
            default:
            {
              // shouldn't get here
              assert(false);
              return eli::geom::general::NOT_CONNECTED;
              break;
            }
          }

          // shouldn't get here
          assert(false);
          return eli::geom::general::NOT_CONNECTED;
        }
      }

      template<typename curve1__, typename curve2__, typename tol__>
      eli::geom::general::continuity report_point_continuity(const curve1__ &curve1, const typename curve1__::data_type &dt1,
                                                             const curve2__ &curve2, const typename curve2__::data_type &dt2,
                                                             const tol__ &tol)
      {
        return internal::report_point_continuity(curve1, dt1, curve2, dt2, eli::geom::general::NOT_CONNECTED, tol);
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
      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__ >
      class piecewise;
    }

    namespace intersect
    {
      template<template<typename, unsigned short, typename> class curve1__, typename data1__, unsigned short dim1__, typename tol1__>
      typename curve::piecewise<curve1__, data1__, dim1__, tol1__>::data_type
                      minimum_distance(
                              typename curve::piecewise<curve1__, data1__, dim1__, tol1__>::data_type &t,
                              const curve::piecewise<curve1__, data1__, dim1__, tol1__> &pc,
                              const typename curve::piecewise<curve1__, data1__, dim1__, tol1__>::point_type &pt);
    }

    namespace curve
    {
      // forward declaration of length function used in methods below. The length function
      // includes this header.
      template<typename curve__>
      void length(typename curve__::data_type &len, const curve__ &c, const typename curve__::data_type &tol);
      template<typename curve__>
      void length(typename curve__::data_type &len, const curve__ &c, const typename curve__::data_type &t0, const typename curve__::data_type &t1, const typename curve__::data_type &tol);

      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__=eli::util::tolerance<data__> >
      class piecewise
      {
        public:
          typedef curve__<data__, dim__, tol__> curve_type;
          typedef typename curve_type::index_type index_type;
          typedef typename curve_type::point_type point_type;
          typedef typename curve_type::control_point_type control_point_type;
          typedef typename curve_type::rotation_matrix_type rotation_matrix_type;
          typedef typename curve_type::bounding_box_type bounding_box_type;
          typedef data__ data_type;
          typedef unsigned short dimension_type;
          typedef tol__ tolerance_type;

          typedef piecewise<curve__, data_type, dim__> piecewise_curve_type;

          typedef typename curve_type::onedbezcurve onedbezcurve;
          typedef piecewise<curve__, data_type, 1, tol__> onedpiecewisecurve;
          typedef piecewise<curve__, data_type, 2, tol__> twodpiecewisecurve;
          typedef piecewise<curve__, data_type, 3, tol__> threedpiecewisecurve;
          typedef piecewise<curve__, data_type, 4, tol__> fourdpiecewisecurve;

          enum error_code
          {
            NO_ERRORS=0,
            INVALID_INDEX=1,
            INDEX_NOT_FOUND=2,
            INVALID_PARAM=50,
            INVALID_PARAM_DIFFERENCE=51,
            SEGMENT_NOT_CONNECTED=100,
            UNKNOWN_ERROR=999
          };

          enum mod_type
          {
            FLAT,
            ROUND,
            EDGE,
            SHARP
          };

        public:
          piecewise() : tmax(0) {}
          piecewise(const piecewise<curve__, data_type, dim__, tol__> &p) : segments(p.segments), tmax(p.tmax), tol(p.tol) {}
          ~piecewise() {}

          bool operator==(const piecewise<curve__, data_type, dim__> &p) const
          {
            if (this==&p)
              return true;
            if (tmax!=p.tmax)
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

          piecewise & operator=(const piecewise<curve__, data_type, dim__> &p)
          {
            if (this==&p)
              return (*this);
            segments=p.segments;
            tmax=p.tmax;
            tol=p.tol;

            return (*this);
          }

          bool operator!=(const piecewise<curve__, data_type, dim__> &p) const
          {
            return !operator==(p);
          }

          bool abouteq(const piecewise<curve__, data_type, dim__> &p, const data_type ttol2 ) const
          {
            if (this==&p)
              return true;
            if (tmax!=p.tmax)
              return false;
            if (tol!=p.tol)
              return false;
            if (number_segments()!=p.number_segments())
              return false;
            typename segment_collection_type::const_iterator scit, it;
            for (scit=segments.begin(), it=p.segments.begin(); scit!=segments.end(); ++scit, ++it)
            {
              if ( ! it->second.abouteq( scit->second, ttol2 ) )
              return false;
            }

            return true;
          }

          static dimension_type dimension() {return dim__;}

          data_type get_tmax() const {return tmax;}

          data_type get_t0() const
          {
            return get_parameter_min();
          }

          void set_tmax(const data_type &tmax_in)
          {
            tmax = tmax_in;
          }

          void set_t0(const data_type &t0_in)
          {
            if(!segments.empty())
            {
              if(t0_in != segments.begin()->first)
              {
                data_type t = t0_in;
                segment_collection_type shiftseg;
                for (typename segment_collection_type::iterator it=segments.begin(); it!=segments.end(); ++it)
                {
                  data_type delta_t = get_delta_t(it);

                  shiftseg.insert(shiftseg.end(), std::make_pair(t, it->second));

                  t+=delta_t;
                }
                segments.swap(shiftseg);
                tmax = t;
              }
            }
            else
            {
              tmax=t0_in;
            }
          }

          void set_t( const data_type &t_old, const data_type &t_new )
          {
            tol__ tol;

            if(!segments.empty())
            {
              segment_collection_type shiftseg;
              for (typename segment_collection_type::iterator it=segments.begin(); it!=segments.end(); ++it)
              {

                if( tol.approximately_equal(t_old, it->first) )
                {
                  shiftseg.insert( shiftseg.end(), std::make_pair( t_new, it->second) );
                }
                else
                {
                  shiftseg.insert( shiftseg.end(), std::make_pair( it->first, it->second) );
                }
              }
              segments.swap(shiftseg);
            }
          }

          void scale_t( const data_type &t_min_new, const data_type &t_max_new )
          {
            if(segments.empty())
            {
              tmax = t_min_new;
            }
            else
            {
              data_type t_min = segments.begin()->first;
              data_type t_max = tmax;
              data_type t_scale = (t_max_new - t_min_new) / (t_max - t_min);

              segment_collection_type scaleseg;
              for (typename segment_collection_type::iterator it=segments.begin(); it!=segments.end(); ++it)
              {
                data_type t = t_min_new + t_scale * (it->first - t_min);
                scaleseg.insert( scaleseg.end(), std::make_pair( t, it->second) );
              }
              segments.swap(scaleseg);
              tmax = t_max_new;
            }
          }

          data_type get_parameter_min() const
          {
            if(!segments.empty())
              return segments.begin()->first;
            else
              return tmax;
          }

          data_type get_parameter_max() const
          {
            return tmax;
          }

          template<typename it__>
          void get_parameters(it__ itt) const
          {
            typename segment_collection_type::const_iterator its;

            for (its=segments.begin(); its!=segments.end(); ++its)
            {
              (*itt)=its->first;++itt;
            }
            (*itt)=tmax;
          }

          void parameter_report() const
          {
            printf("Parameter report:\n");
            typename segment_collection_type::const_iterator it;

            int i = 0;
            // cycle through all segments to get each bounding box to add
            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              printf(" seg: %d \t t: %f\n", i, it->first);
              ++i;
            }
            printf(" tmax: %f\n", tmax);
            printf("End report\n");
          }

          void get_pmap( std::vector < data_type > &pmap )
          {
            pmap.clear();

            typename segment_collection_type::const_iterator it;
            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              pmap.push_back( it->first );
            }
            pmap.push_back( tmax );
          }

          index_type number_segments() const {return static_cast<index_type>(segments.size());}

          void degree(index_type &mind, index_type &maxd) const
          {
            typename segment_collection_type::const_iterator it;

            it=segments.begin();

            index_type d = it->second.degree();
            mind = d;
            maxd = d;
            ++it;

            // cycle through all segments to get each bounding box to add
            for (; it!=segments.end(); ++it)
            {
              index_type d = it->second.degree();

              if(d<mind)
                mind=d;
              if(d>maxd)
                maxd=d;
            }
          }

          template<typename it__>
          void degrees(it__ itd)
          {
            typename segment_collection_type::const_iterator it;

            // cycle through all segments to get each degree
            for (it=segments.begin(); it!=segments.end(); ++it, ++itd)
            {
              (*itd) = it->second.degree();
            }
          }

          void get_bounding_box(bounding_box_type &bb) const
          {
            typename segment_collection_type::const_iterator it;
            bounding_box_type bb_local;

            bb.clear();

            // cycle through all segments to get each bounding box to add
            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->second.get_bounding_box(bb_local);
              bb.add(bb_local);
            }
          }

          void scale(const data_type &s)
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->second.scale(s);
            }
          }

          void rotate(const rotation_matrix_type &rmat)
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->second.rotate(rmat);
            }
          }

          void rotate(const rotation_matrix_type &rmat, const point_type &rorig)
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->second.rotate(rmat, rorig);
            }
          }

          void translate(const point_type &trans)
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->second.translate(trans);
            }
          }

          bool closed() const
          {
            typename segment_collection_type::const_iterator itlast, itfirst;

            itlast=segments.end(); --itlast;
            itfirst=segments.begin();

            data_type dtlast = get_delta_t(itlast);
            data_type dtfirst = get_delta_t(itfirst);

            return eli::geom::utility::check_point_continuity(itlast->second, dtlast, itfirst->second, dtfirst, eli::geom::general::C0, tol);
          }

          bool open() const
          {
            return !closed();
          }

          void reverse()
          {
            segment_collection_type rseg;
            typename segment_collection_type::iterator itr;
            typename segment_collection_type::iterator itrguess = rseg.begin();

            data_type t(get_parameter_min());

            for (typename segment_collection_type::reverse_iterator it=segments.rbegin(); it!=segments.rend(); ++it)
            {
              itr = rseg.insert(itrguess, std::make_pair(t, it->second));

              // reverse each segment
              itr->second.reverse();

              data_type delta_t = get_delta_t(it);
              t += delta_t;

              itrguess = itr;
            }
            segments.swap(rseg);

            // Parametric length should stay the same.
            // assert(t == tmax);

            // check if still connected
            assert(check_continuity(eli::geom::general::C0));
          }


          void reflect_xy()
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->second.reflect_xy();
            }
          }

          void reflect_xz()
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->second.reflect_xz();
            }
          }

          void reflect_yz()
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->second.reflect_yz();
            }
          }

          void reflect(const point_type &normal)
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->second.reflect(normal);
            }
          }

          void reflect(const point_type &normal, const data_type &d)
          {
            typename segment_collection_type::iterator it;

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->second.reflect(normal, d);
            }
          }

          void clear() {segments.clear();}

          template<typename it__>
          error_code set(it__ itb, it__ ite)
          {
            segments.clear();

            index_type i;
            it__ it;
            for (i=0, it=itb; it!=ite; ++i, ++it)
            {
              error_code err=push_back(*it);
              if (err!=NO_ERRORS)
              {
                segments.clear();
                return err;
              }
            }

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERRORS;
          }

          template<typename it__, typename itd__>
          error_code set(it__ itb, it__ ite, itd__ itd)
          {
            segments.clear();

            index_type i;
            it__ it;
            itd__ itdt;
            for (i=0, it=itb, itdt=itd; it!=ite; ++i, ++it, ++itdt)
            {
              error_code err=push_back(*it, *itdt);
              if (err!=NO_ERRORS)
              {
                segments.clear();
                return err;
              }
            }

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERRORS;
          }

          error_code push_front(const curve_type &curve, const data_type &dt=1.0)
          {
            if (dt<=0)
              return INVALID_PARAM_DIFFERENCE;

            // check to make sure have valid segments
            if (!segments.empty())
            {
              data_type start_dt = get_delta_t(segments.begin());
              if (!eli::geom::utility::check_point_continuity(curve, dt, segments.begin()->second, start_dt, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // add segment
            data_type t0(get_parameter_min());
            t0 -= dt;
            segments.insert(segments.begin(), std::make_pair(t0, curve));

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERRORS;
          }

          error_code push_front(const piecewise<curve__, data_type, dim__, tol__> &p)
          {
            typename segment_collection_type::const_reverse_iterator itp;
            error_code err;

            for (itp=p.segments.rbegin(); itp!=p.segments.rend(); ++itp)
            {
              err=push_front(itp->second, p.get_delta_t(itp));
              if (err!=NO_ERRORS)
              {
                return err;
              }
            }

            return NO_ERRORS;
          }

          error_code push_back(const curve_type &curve, const data_type &dt=1.0)
          {
            if (dt<=0)
              return INVALID_PARAM_DIFFERENCE;

            // check to make sure have valid segments
            if (!segments.empty())
            {
              data_type end_dt = get_delta_t(segments.rbegin());
              if (!eli::geom::utility::check_point_continuity(segments.rbegin()->second, end_dt, curve, dt, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // add segment
            segments.insert(segments.end(), std::make_pair(tmax, curve));
            tmax+=dt;

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERRORS;
          }

          error_code push_back(const piecewise<curve__, data_type, dim__, tol__> &p)
          {
            typename segment_collection_type::const_iterator itp;
            error_code err;

            for (itp=p.segments.begin(); itp!=p.segments.end(); ++itp)
            {
              err=push_back(itp->second, p.get_delta_t(itp));
              if (err!=NO_ERRORS)
              {
                return err;
              }
            }

            return NO_ERRORS;
          }

          error_code get(curve_type &curve, const index_type &index) const
          {
            data_type dt;
            return get(curve, dt, index);
          }

          error_code degree_promote()
          {
            typename segment_collection_type::const_iterator scit;
            for (scit=segments.begin(); scit!=segments.end(); ++scit)
            {
              scit->second.degree_promote();
            }

            return NO_ERRORS;
          }

          error_code degree_promote(const index_type &index)
          {
            if (index>=number_segments())
              return INVALID_INDEX;

            typename segment_collection_type::iterator scit;
            find_segment(scit, index);

            scit->second.degree_promote();
            return NO_ERRORS;
          }

          error_code degree_promote_to(const index_type &deg)
          {
            typename segment_collection_type::const_iterator scit;
            for (scit=segments.begin(); scit!=segments.end(); ++scit)
            {
              scit->second.degree_promote_to(deg);
            }

            return NO_ERRORS;
          }

          error_code degree_promote_to(const index_type &index, const index_type &deg)
          {
            if (index>=number_segments())
              return INVALID_INDEX;

            typename segment_collection_type::iterator scit;
            find_segment(scit, index);

            scit->second.degree_promote_to(deg);
            return NO_ERRORS;
          }

          error_code get(curve_type &curve, data_type &dt, const index_type &index) const
          {
            if (index>=number_segments())
              return INVALID_INDEX;

            // advance to desired index
            index_type i;
            typename segment_collection_type::const_iterator scit;
            for (i=0, scit=segments.begin(); i<index; ++i, ++scit) {}

            curve=scit->second;
            dt=get_delta_t(scit);
            return NO_ERRORS;
          }

          error_code replace(const curve_type &curve, const index_type &index)
          {
            if (index>=number_segments())
              return INVALID_INDEX;

            typename segment_collection_type::iterator scit;
            data_type dt, dto;

            find_segment(scit, index);

            dt = get_delta_t(scit);

            typename segment_collection_type::iterator scito;

            // check the connectivity on adjacent nodes (if available)
            if (index>0)
            {
              scito=scit;
              --scito;
              dto = get_delta_t(scito);
              if (!eli::geom::utility::check_point_continuity(scito->second, dto, curve, dt, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if ((index+1)<number_segments())
            {
              scito=scit;
              ++scito;
              dto = get_delta_t(scito);
              if (!eli::geom::utility::check_point_continuity(curve, dt, scito->second, dto, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            // set the new curve
            scit->second=curve;

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERRORS;
          }

          error_code replace(const curve_type &curve, const index_type &index0, const index_type &index1)
          {
            if (index0>=number_segments())
              return INVALID_INDEX;
            if (index1>number_segments())
              return INVALID_INDEX;
            if (index0>=index1)
              return INVALID_INDEX;

            typename segment_collection_type::iterator scit0, scit1, scito;
            find_segment(scit0, index0);
            find_segment(scit1, index1);

            data_type dt0, dt1, dto;
            dt0 = get_delta_t(scit0);
            dt1 = get_delta_t(scit1);

            // check the connectivity on adjacent nodes (if available)
            if (index0>0)
            {
              scito=scit0;
              --scito;
              dto = get_delta_t(scito);
              if (!eli::geom::utility::check_point_continuity(scito->second, dto, curve, dt0, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if ((index1+1)<number_segments())
            {
              scito=scit1;
              dto = get_delta_t(scito);
              if (!eli::geom::utility::check_point_continuity(curve, dt1, scito->second, dto, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            data_type t = scit0->first;

            typename segment_collection_type::iterator itguess;

            // erase old segments
            itguess = segments.erase(scit0, scit1);

            segments.insert(itguess, std::make_pair(t, curve));

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERRORS;
          }

          error_code replace(const piecewise<curve__, data_type, dim__> &p, const index_type &index)
          {
            if (index>=number_segments())
              return INVALID_INDEX;

            typename segment_collection_type::iterator scit;
            find_segment(scit, index);

            return replace_it( p, scit );
          }

          error_code replace_t(const piecewise<curve__, data_type, dim__> &p, const data_type &t0)
          {
            data_type tt;
            typename segment_collection_type::iterator scit;
            find_segment(scit, tt, t0 );

            return replace_it( p, scit );
          }





          /** index0 - start index of segments to replace
            * index1 - index of the first segment after the ones to be replaced
            */

          error_code replace(const piecewise<curve__, data_type, dim__> &p, const index_type &index0, const index_type &index1)
          {
            if (index0>=number_segments())
              return INVALID_INDEX;
            if (index1>number_segments())
              return INVALID_INDEX;
            if (index0>=index1)
              return INVALID_INDEX;

            typename segment_collection_type::iterator scit0, scit1, scito;
            find_segment(scit0, index0);
            find_segment(scit1, index1);

            // Find parameter span to insert
            data_type pt0(p.get_parameter_min()), ptmax(p.get_parameter_max()), ptspan(ptmax - pt0);

            // get the first and last curve to insert
            typename segment_collection_type::const_iterator itps, itpe;
            curve_type cs, ce;
            data_type dts, dte;

            itps = p.segments.begin();
            cs = itps->second;
            dts = p.get_delta_t(itps);

            itpe = p.segments.end(); --itpe;
            ce = itpe->second;
            dte = p.get_delta_t(itpe);

            // Find parameter span of segment to replace
            data_type dti, dto;

            if(scit1 != segments.end())
              dti = scit1->first - scit0->first;
            else
              dti = tmax - scit0->first;

            // Ratio of parameter lengths, replace/insert
            data_type pratio = dti/ptspan;

            // check the connectivity on adjacent nodes (if available)
            if (index0>0)
            {
              scito=scit0;
              --scito;
              dto = get_delta_t(scito);
              if (!eli::geom::utility::check_point_continuity(scito->second, dto, cs, dts*pratio, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (index1<number_segments())
            {
              scito=scit1;
              dto = get_delta_t(scito);
              if (!eli::geom::utility::check_point_continuity(ce, dte*pratio, scito->second, dto, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }

            data_type t = scit0->first;
            data_type dtp;

            typename segment_collection_type::const_iterator it;
            typename segment_collection_type::iterator itguess;

            // erase old segments
            itguess = segments.erase(scit0, scit1);

            for(it=p.segments.begin(); it != p.segments.end(); ++it)
            {
              itguess = segments.insert(itguess, std::make_pair(t, it->second));
              dtp = p.get_delta_t(it);
              t += dtp*pratio;
            }

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERRORS;
          }

          error_code split(const data_type &t)
          {
            // find segment that corresponds to given t
            typename segment_collection_type::iterator it;
            data_type tt;
            find_segment(it, tt, t);

            // do some checking to see if even need to split
            if (tol.approximately_equal(tt, 0))
            	return NO_ERRORS;
            if (tol.approximately_equal(tt, 1))
            	return NO_ERRORS;
            if (it==segments.end())
              return INVALID_PARAM;

            // split the segment and replace
            return split_seg(it, tt);
          }

          error_code split(piecewise<curve__, data_type, dim__> &before, piecewise<curve__, data_type, dim__> &after, const data_type &tsplit) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it, iit;
            data_type tt;
            error_code ec;

            find_segment(it, tt, tsplit);

            // do some checking to see if even need to split
            if (it==segments.end())
              return INVALID_PARAM;

            // if found end of segment adjust to make beginning of next
            if (tol.approximately_equal(tt, 1))
            {
              ++it;
              tt=0;

              // catch case that wanted to split at end
              if (it==segments.end())
              {
                before=(*this);
                after.clear();
                return NO_ERRORS;
              }
            }

            // catch case that wanted to split at beginning
            if ((tol.approximately_equal(tt, 0)) && (it==segments.begin()))
            {
              before.clear();
              after=(*this);
              return NO_ERRORS;
            }

            // fill first half
            before.clear();
            after.clear();
            before.set_t0(get_t0());
            for (iit=segments.cbegin(); iit!=it; ++iit)
            {
              ec=before.push_back(iit->second, get_delta_t(iit));
              if (ec!=NO_ERRORS)
              {
                before.clear();
                assert(false);
                return ec;
              }
            }

            // split the segment and replace (if needed)
            if (!tol.approximately_equal(tt, 0))
            {
              data_type delta_t(get_delta_t(it));
              curve_type cl, cr;

              // split underlying curve and add left to before and right to after
              it->second.split(cl, cr, tt);
              ec=before.push_back(cl, tt*delta_t);
              if (ec!=NO_ERRORS)
              {
                before.clear();
                assert(false);
                return ec;
              }

              after.set_t0(it->first+tt*delta_t);
              after.push_back(cr, (1-tt)*delta_t);
            }
            else
            {
              after.set_t0(it->first);
              after.push_back(it->second, get_delta_t(it));
            }

            // fill second half
            iit=it;
            for (++iit; iit!=segments.cend(); ++iit)
            {
              ec=after.push_back(iit->second, get_delta_t(iit));
              if (ec!=NO_ERRORS)
              {
                before.clear();
                after.clear();
                assert(false);
                return ec;
              }
            }

            return NO_ERRORS;
          }

          void match_pmap( piecewise<curve__, data_type, dim__> &other )
          {
            std::vector < data_type > pmap, omap, cmap;

            get_pmap( pmap );
            other.get_pmap( omap );


            tolerance_type ttol(tol);

            // Comparison function for set_union.
            auto comp = [&ttol](const data_type &x1, const data_type &x2)->bool
            {
              return ttol.approximately_less_than(x1, x2);
            };

            // Place union of pmap and omap into cmap.
            std::set_union( pmap.begin(), pmap.end(), omap.begin(), omap.end(), std::back_inserter(cmap), comp );

            for ( int i = 0; i < cmap.size(); i++ )
            {
              split( cmap[i] );
              other.split( cmap[i] );
            }
          }


          void to_cubic(const data_type &ttol)
          {
            typename segment_collection_type::iterator it;
            for (it=segments.begin(); it!=segments.end(); ++it)
              segment_to_cubic(it, ttol);
          }

          void round(const data_type &rad)
          {
            // catch special case of no rounding
            if (rad<=0)
            {
              return;
            }

            // Note: this could be implemented more efficiently if wanted to track the
            //       segment container iterators, but that would require more code
            //       duplication since the round(rad, i) method calls other methods with
            //       useful error checking.

            // Note: need to keep calling number_segments() because a call to round(rad, i)
            //       might increase the number of sections and need to make sure that this
            //       loop gets to the last segment
            for (index_type i=0; i<=number_segments(); ++i)
            {
              round(rad, i);
            }
          }

          bool round(const data_type &rad, const index_type &joint)
          {
            // if joint doesn't exist then return
            if ((joint<0) || (joint>number_segments()))
            {
              assert(false);
              return false;
            }

            // catch special case of no rounding
            if (rad<=0)
            {
              return false;
            }

            // catch special case of rounding first or last joint of open curve
            if (((joint==0) || (joint==number_segments())) && open())
            {
              return false;
            }

            index_type im1, i;
            bool rounding_end(false);

            if ((joint==0) || (joint==number_segments()))
            {
              i=0;
              im1=number_segments()-1;
              rounding_end=true;
            }
            else
            {
              i=joint;
              im1=joint-1;
            }

            curve_type cim1, ci, arc, c, ctrim;
            data_type dtim1(0), dti(0), r, lenim1, leni, tim1_split(-1), ti_split(-1);
            point_type fim1, fi, fpim1, fpi;
            control_point_type cp[4];

            // get the two curve segments
            get(cim1, dtim1, im1);
            get(ci, dti, i);

            // check to see if joint needs to be rounded
            fpim1=cim1.fp(1); fpim1.normalize();
            fpi=ci.fp(0); fpi.normalize();
            if (tol.approximately_equal(fpi.dot(fpim1), 1))
            {
              return false;
            }

            // determine what the actual radius will be
            r=rad;
            eli::geom::curve::length(lenim1, cim1, tol.get_absolute_tolerance());
            eli::geom::curve::length(leni, ci, tol.get_absolute_tolerance());
            if (lenim1<r)
            {
              tim1_split=0;
              r=lenim1;
            }
            if (leni<r)
            {
              tim1_split=-1;
              ti_split=1;
              r=leni;
            }

            // find coordinate that corresponds to location of radius on each curve
            double small_t(1e-3);
            if (tim1_split<0)
            {
              // FIX: this is exact for straight lines and approximate for other curves
              tim1_split=1-r/lenim1;

              // if resulting segment is too small then make entire edge the round
              if (tim1_split<small_t)
              {
            	tim1_split=0;
              }
            }
            if (ti_split<0)
            {
              // FIX: this is exact for straight lines and approximate for other curves
              ti_split=r/leni;

              // if resulting segment is too small then make entire edge the round
              if (ti_split>(1-small_t))
              {
            	ti_split=1;
              }
            }

            // if resulting round is too small compared then don't round
            if ((tim1_split>(1-small_t)) || (ti_split<small_t))
            	return false;

            // calculate the points & slopes for end of round
            fim1=cim1.f(tim1_split);
            fpim1=cim1.fp(tim1_split); fpim1.normalize();
            fi=ci.f(ti_split);
            fpi=ci.fp(ti_split); fpi.normalize();

            // build round curve
            data_type k=static_cast<data_type>(4)*(eli::constants::math<data_type>::sqrt_two()-static_cast<data_type>(1))/static_cast<data_type>(3); // use this value so that 90 deg. corners will have close circle
            arc.resize(3);
            cp[0]=fim1;
            cp[1]=fim1+fpim1*(k*r);
            cp[2]=fi-fpi*(k*r);
            cp[3]=fi;
            arc.set_control_point(cp[0], 0);
            arc.set_control_point(cp[1], 1);
            arc.set_control_point(cp[2], 2);
            arc.set_control_point(cp[3], 3);

            // split the two curves
            if (tim1_split>0)
            {
              cim1.split(c, ctrim, tim1_split); cim1=c;
            }
            if (ti_split<1)
            {
              ci.split(ctrim, c, ti_split); ci=c;
            }

            // replace/add curves
            error_code ec;
            if (rounding_end)
            {
              curve_type arc1, arc2;
              piecewise<curve__, data_type, dim__> pct0, pct1;

              // split the arc
              arc.split(arc1, arc2, static_cast<data_type>(0.5));

              // put the ith segment and arc2 onto the end of the list of curves
              if ( tim1_split > 0 )
              {
                ec=pct0.push_back(cim1, dtim1*tim1_split);
                if (ec!=NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              if ( tim1_split < 1 )
              {
                ec=pct0.push_back(arc1, dtim1*(static_cast<data_type>(1)-tim1_split));
                if (ec!=NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }

              ec=replace(pct0, number_segments()-1);
              if (ec!=NO_ERRORS)
              {
                assert(false);
                return false;
              }

              // put arc1 and the (i-1)st segment onto the front of the list of curves
              if ( ti_split < 1 )
              {
                ec=pct1.push_front(ci, dti*(1-ti_split));
                if (ec!=NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              if ( ti_split > 0 )
              {
                ec=pct1.push_front(arc2, dti*ti_split);
                if (ec!=NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }

              ec=replace(pct1, 0);
              if (ec!=NO_ERRORS)
              {
                assert(false);
                return false;
              }
            }
            else
            {
              piecewise<curve__, data_type, dim__> pct;

              if (tim1_split>0)
              {
            	ec=pct.push_back(cim1, dtim1*tim1_split);
                if (ec!=NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              ec=pct.push_back(arc, dtim1*(1-tim1_split)+dti*ti_split);
              if (ec!=NO_ERRORS)
              {
                assert(false);
                return false;
              }
              if (ti_split<1)
              {
            	ec=pct.push_back(ci, dti*(1-ti_split));
                if (ec!=NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }

              // replace the two segments with the piecewise curve
              ec=replace(pct, i-1, i+1);
              if (ec!=NO_ERRORS)
              {
                assert(false);
                return false;
              }
            }

            return true;
          }

          // Modify a curve around parameter t.  Remove +-dt and replace with generated curve segment.
          void modify( const index_type &mod, const data_type &t, const data_type &dt, const data_type &len_factor, const data_type &off_factor, const data_type &str_factor )
          {
            data_type tmin = get_parameter_min();
            data_type tmax_save = tmax;

            data_type tstart, tend;

            bool mod_ends = false;

            if ( tol.approximately_equal( tmin, t ) || tol.approximately_equal( tmax, t ) )
            {
              tstart = tmax - dt;
              tend = tmin + dt;
              mod_ends = true;
            }
            else
            {
              tstart = t - dt;
              tend = t + dt;
              mod_ends = false;
            }

            point_type pstart, pend, tanstart, tanend, tjunk;

            // Evaluate point and tangents as needed.
            pstart = f( tstart );
            tangents( tstart, tanstart, tjunk );
            pend = f ( tend );
            tangents( tend, tjunk, tanend );
            tanend = -1.0 * tanend;

            // Variable to store new curve segments.
            curve_type cstart, cend;

            // Distance to close.
            data_type d = dist( pstart, pend );

            if ( tol.approximately_equal( d, 0 ) )
            {
              return;
            }

            switch ( mod )
            {
            case ROUND:
              {
                point_type disp = pstart - pend;
                disp.normalize();

                // Find angles between tangents and displacement
                data_type dottanstart = disp.dot( tanstart );
                data_type thetastart = acos( dottanstart );

                data_type dottanend = disp.dot( tanend );
                data_type thetaend = acos( dottanend );

                // Find circular arc to include
                data_type theta = eli::constants::math<data__>::pi() + thetaend - thetastart;

                if ( tol.approximately_equal( theta, 0 ) )
                {
                  return;
                }

                // Find radius from chord length.
                data_type radius = d / ( 2.0 * sin( theta / 2.0 ) );

                // Distance to quad Bezier construction point.
                data_type b = radius * tan( theta / 4.0 );

                // Extrapolated tip up/down points.
                point_type ptstart, ptend;
                ptstart = pstart + tanstart * b * len_factor;
                ptend = pend + tanend * b * len_factor;

                // Displacement vector at extrapolated tip
                point_type dispt = ptstart - ptend;

                // Point on center of circle.
                point_type psplit = ( ptstart + ptend ) * 0.5 + dispt * off_factor;

                // Fraction of radius to place cubic control point
                data_type f = (4.0/3.0) * tan(theta / 8.0);

                cstart.resize(3);
                cend.resize(3);

                cstart.set_control_point( pstart, 0 );
                cstart.set_control_point( pstart + tanstart * radius * f * len_factor, 1 );
                cstart.set_control_point( psplit + disp * radius * f, 2 );
                cstart.set_control_point( psplit, 3 );

                cend.set_control_point( psplit, 0 );
                cend.set_control_point( psplit - disp * radius * f, 1 );
                cend.set_control_point( pend + tanend * radius * f * len_factor, 2 );
                cend.set_control_point( pend, 3 );
              }
              break;
            case EDGE:
              {
                point_type pt1 = pstart + d * tanstart * len_factor * 0.5;
                point_type pt2 = pend + d * tanend * len_factor * 0.5;

                point_type disp = pstart - pend;
                disp.normalize();

                point_type dispt = pt1 - pt2;
                data_type tipht = dispt.norm();

                point_type psplit = ( pt1 + pt2 ) * 0.5 + disp * tipht * off_factor;

                cstart.resize(1);
                cend.resize(1);

                cstart.set_control_point( pstart, 0 );
                cstart.set_control_point( psplit, 1 );

                cend.set_control_point( psplit, 0 );
                cend.set_control_point( pend, 1 );
              }
              break;
            case SHARP:
              {
                point_type pt1 = pstart + d * tanstart * len_factor * 0.5;
                point_type pt2 = pend + d * tanend * len_factor * 0.5;

                point_type disp = pstart - pend;
                disp.normalize();

                point_type dispt = pt1 - pt2;
                data_type tipht = dispt.norm();

                point_type psplit = ( pt1 + pt2 ) * 0.5 + disp * tipht * off_factor;

                point_type cp1 = pstart + d * tanstart * len_factor * 0.5 * str_factor;
                point_type cp2 = pend + d * tanend * len_factor * 0.5 * str_factor;

                cstart.resize(2);
                cend.resize(2);

                cstart.set_control_point( pstart, 0 );
                cstart.set_control_point( cp1, 1 );
                cstart.set_control_point( psplit, 2 );

                cend.set_control_point( psplit, 0 );
                cend.set_control_point( cp2, 1 );
                cend.set_control_point( pend, 2 );
              }
              break;
            case FLAT:
            default:
              {
                point_type padd = ( pstart + pend ) * 0.5;

                cstart.resize(1);
                cend.resize(1);

                cstart.set_control_point( pstart, 0 );
                cstart.set_control_point( padd, 1 );

                cend.set_control_point( padd, 0 );
                cend.set_control_point( pend, 1 );
                break;
              }
            }

            if ( mod_ends )
            {
              // Make sure curves are split at modification points.
              split( tstart );
              split( tend );

              replace( cend, 0 );
              replace( cstart, number_segments() - 1 );
            }
            else
            {
              piecewise_curve_type cbefore, cjunk, cafter;

              split( cbefore, cjunk, tstart );
              split( cjunk, cafter, tend );

              cbefore.push_back( cstart, dt );
              cbefore.push_back( cend, dt );
              cbefore.push_back( cafter );

              clear();
              set_t0( tmin );
              push_back( cbefore );

              set_tmax( tmax_save ); // Just to be sure.

            }
          }


          bool continuous(eli::geom::general::continuity cont, const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it, itfirst, itsecond;
            data_type tt(0);
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              return false;
            }

            if (tt==0)
            {
              if (it==segments.begin())
              {
                if (open())
                {
                  return (cont==eli::geom::general::NOT_CONNECTED);
                }
                else
                {
                  typename segment_collection_type::const_iterator itlast(segments.end());

                  --itlast;
                  itfirst=itlast;
                  itsecond=it;
                }
              }
              else
              {
                  itsecond=it;
                  itfirst=it; --itfirst;
              }
            }
            else if (tt==1)
            {
              typename segment_collection_type::const_iterator itlast(segments.end());

              --itlast;
              if (it==itlast)
              {
                if (open())
                {
                  return (cont==eli::geom::general::NOT_CONNECTED);
                }
                else
                {
                  itfirst=it;
                  itsecond=segments.begin();
                }
              }
              else
              {
                itfirst=it;
                itsecond=it; ++itsecond;
              }
            }
            else
            {
              return (cont!=eli::geom::general::NOT_CONNECTED);
            }

            data_type dtfirst = get_delta_t(itfirst);
            data_type dtsecond = get_delta_t(itsecond);

            // check the continuity of the two sections
            return eli::geom::utility::check_point_continuity(itfirst->second, dtfirst, itsecond->second, dtsecond, cont, tol);
          }

          eli::geom::general::continuity continuity(const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it, itfirst, itsecond;
            data_type tt(0);
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              return eli::geom::general::NOT_CONNECTED;
            }

            if (tt==0)
            {
              if (it==segments.begin())
              {
                if (open())
                {
                  return eli::geom::general::NOT_CONNECTED;
                }
                else
                {
                  typename segment_collection_type::const_iterator itlast(segments.end());

                  --itlast;
                  itfirst=itlast;
                  itsecond=it;
                }
              }
              else
              {
                itsecond=it;
                itfirst=it; --itfirst;
              }
            }
            else if (tt==1)
            {
              typename segment_collection_type::const_iterator itlast(segments.end());

              --itlast;
              if (it==itlast)
              {
                if (open())
                {
                  return eli::geom::general::NOT_CONNECTED;
                }
                else
                {
                  itfirst=it;
                  itsecond=segments.begin();
                }
              }
              else
              {
                itfirst=it;
                itsecond=it; ++itsecond;
              }
            }
            else
            {
              return eli::geom::general::CINFINITY;
            }

            data_type dtfirst = get_delta_t(itfirst);
            data_type dtsecond = get_delta_t(itsecond);

            // check the continuity of the two sections
            return eli::geom::utility::report_point_continuity(itfirst->second, dtfirst, itsecond->second, dtsecond, tol);
          }

          bool smooth(const data_type &angle_tol, const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it, itfirst, itsecond;
            data_type tt(0);
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              return false;
            }

            if (tt==0)
            {
              if (it==segments.begin())
              {
                if (open())
                {
                  return false;
                }
                else
                {
                  typename segment_collection_type::const_iterator itlast(segments.end());

                  --itlast;
                  itfirst=itlast;
                  itsecond=it;
                }
              }
              else
              {
                  itsecond=it;
                  itfirst=it; --itfirst;
              }
            }
            else if (tt==1)
            {
              typename segment_collection_type::const_iterator itlast(segments.end());

              --itlast;
              if (it==itlast)
              {
                if (open())
                {
                  return false;
                }
                else
                {
                  itfirst=it;
                  itsecond=segments.begin();
                }
              }
              else
              {
                itfirst=it;
                itsecond=it; ++itsecond;
              }
            }
            else
            {
              return true;
            }

            data_type dtfirst = get_delta_t(itfirst);
            data_type dtsecond = get_delta_t(itsecond);

            // check if coincident points
            if (eli::geom::utility::check_point_continuity(itfirst->second, dtfirst, itsecond->second, dtsecond, eli::geom::general::C0, tol))
            {
              // check if angles between fp's are less than angle
              point_type fp1(itfirst->second.fp(1));
              point_type fp2(itsecond->second.fp(0));
              data_type val;

              val=fp1.norm();
              // NOTE: This is a case where the geometric tolerance object should catch this
              if (val>tol.get_relative_tolerance())
              {
                fp1/=val;
              }
              else
              {
                fp1.setZero();
                // catch case where both vectors are nearly zero
                if (tol.approximately_equal(fp1, fp2))
                {
                  return true;
                }
              }
              val=fp2.norm();
              // NOTE: This is a case where the geometric tolerance object should catch this
              if (val>tol.get_relative_tolerance())
              {
                fp2/=val;
              }
              else
              {
                fp2.setZero();
                // since other vector was non zero (otherwise previous else case would have triggered)
                if (tol.approximately_equal(fp1, fp2))
                {
                  return true;
                }
                // know the vectors are not in same direction
                else
                {
                  return false;
                }
              }

              val=std::abs(1-fp1.dot(fp2));
              return (val<=angle_tol);
            }

            // check the continuity of the two sections
            return false;
          }

          void find_discontinuities(eli::geom::general::continuity cont, std::vector<data_type> &tdisc) const
          {
            // clear input vector
            tdisc.clear();

            index_type i, istart(0), njoints(number_segments()+1);
            std::vector<data_type> joints;

            // get all of the joints
            get_parameters(std::back_inserter(joints));

            // if curve is open then don't check last joint
            if (open())
            {
              ++istart;
              --njoints;
            }

            // check each joint (after starting joint) if it is continuous
            for (i=istart; i<njoints; ++i)
            {
              if (!continuous(cont, joints[i]))
              {
                tdisc.push_back(joints[i]);
              }
            }
          }

          void find_discontinuities(const data_type &angle_tol, std::vector<data_type> &tdisc) const
          {
            // clear input vector
            tdisc.clear();

            index_type i, istart(0), njoints(number_segments()+1);
            std::vector<data_type> joints;

            // get all of the joints
            get_parameters(std::back_inserter(joints));

            // if curve is open then don't check last joint
            if (open())
            {
              ++istart;
              --njoints;
            }

            // check each joint (after starting joint) if it is continuous
            for (i=istart; i<njoints; ++i)
            {
              if (!smooth(angle_tol, joints[i]))
              {
                tdisc.push_back(joints[i]);
              }
            }
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

            return it->second.f(tt);
          }

          point_type fp(const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it;
            data_type tt, delta_t;
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              --it;
            }

            delta_t = get_delta_t(it);
            return it->second.fp(tt)/delta_t;
          }

          void fps(const data_type &t, point_type &fp1, point_type &fp2) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it;
            data_type tt(0), delta_t;
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              --it;
            }

            if ( tol.approximately_equal( tt, 0 ) )
            {
              delta_t = get_delta_t(it);
              fp2 = it->second.fp( tt )/delta_t;

              if ( it == segments.begin() )
              {
                fp1 = fp2;
              }
              else
              {
                --it;
                delta_t = get_delta_t(it);
                fp1 = it->second.fp( 1 )/delta_t;
              }
            }
            else if ( tol.approximately_equal( tt, 1 ) )
            {
              delta_t = get_delta_t(it);
              fp1 = it->second.fp( tt )/delta_t;
              ++it;
              if ( it == segments.end() )
              {
                fp2 = fp1;
              }
              else
              {
                delta_t = get_delta_t(it);
                fp2 = it->second.fp( 0 )/delta_t;
              }
            }
            else
            {
              delta_t = get_delta_t(it);
              fp1 = it->second.fp(tt)/delta_t;
              fp2 = fp1;
            }
          }

          point_type fpp(const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it;
            data_type tt, delta_t;
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              --it;
            }

            delta_t = get_delta_t(it);
            return it->second.fpp(tt)/(delta_t*delta_t);
          }

          point_type fppp(const data_type &t) const
          {
            // find segment that corresponds to given t
            typename segment_collection_type::const_iterator it;
            data_type tt, delta_t;
            find_segment(it, tt, t);

            if (it==segments.end())
            {
              assert(false);
              --it;
            }

            delta_t = get_delta_t(it);
            return it->second.fppp(tt)/(delta_t*delta_t*delta_t);
          }

          point_type tangent(const data_type &t) const
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

            return it->second.tangent(tt);
          }

          void tangents(const data_type &t, point_type &t1, point_type &t2) const
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

            if ( tol.approximately_equal( tt, 0 ) )
            {
              t2 = it->second.tangent( tt );

              if ( it == segments.begin() )
              {
                t1 = t2;
              }
              else
              {
                --it;
                t1 = it->second.tangent( 1 );
              }
            }
            else if ( tol.approximately_equal( tt, 1 ) )
            {
              t1 = it->second.tangent( tt );
              ++it;
              if ( it == segments.end() )
              {
                t2 = t1;
              }
              else
              {
                t2 = it->second.tangent( 0 );
              }
            }
            else
            {
              t1 = it->second.tangent(tt);
              t2 = t1;
            }
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

            it->second.frenet_serret_frame(t, n, b, tt);
          }

          void product(const piecewise<curve__, data_type, dim__> &a, const piecewise<curve__, data_type, dim__> &b)
          {
            set_t0( a.get_t0() );

            typename segment_collection_type::const_iterator scita, scitb;
            scita=a.segments.begin();
            scitb=b.segments.begin();
            for ( ; scita!=a.segments.end(); ++scita, ++scitb)
            {
              curve_type c;

              c.product( scita->second, scitb->second );

              push_back( c, a.get_delta_t(scita) );
            }
          }

          void product1d(const piecewise<curve__, data_type, dim__> &a, const piecewise<curve__, data_type, 1> &b)
          {
            set_t0( a.get_t0() );

            typename segment_collection_type::const_iterator scita;
            typename onedpiecewisecurve::segment_collection_type::const_iterator scitb;
            scita=a.segments.begin();
            scitb=b.segments.begin();
            for ( ; scita!=a.segments.end(); ++scita, ++scitb)
            {
              curve_type c;

              c.product1d( scita->second, scitb->second );

              push_back( c, a.get_delta_t(scita) );
            }
          }

          void square(const piecewise<curve__, data_type, dim__> &a)
          {
            set_t0( a.get_t0() );

            typename segment_collection_type::const_iterator scita;
            for ( scita=a.segments.begin(); scita!=a.segments.end(); ++scita)
            {
              curve_type c;

              c.square( scita->second );

              push_back( c, a.get_delta_t(scita) );
            }
          }

          void sqrt(const piecewise<curve__, data_type, dim__> &a)
          {
            set_t0( a.get_t0() );

            typename segment_collection_type::const_iterator scita;
            for ( scita=a.segments.begin(); scita!=a.segments.end(); ++scita)
            {
              curve_type c, b( scita->second );

              bool flip = false;

              if ( b.f(0).x() < 1e-12 )
              {
                b.reverse();
                flip = true;
              }

              c.sqrt( b );

              if ( flip )
              {
                c.reverse();
              }

              push_back( c, a.get_delta_t(scita) );
            }
          }

          void sum(const piecewise<curve__, data_type, dim__> &a, const piecewise<curve__, data_type, dim__> &b)
          {
            set_t0( a.get_t0() );

            typename segment_collection_type::const_iterator scita, scitb;
            scita=a.segments.begin();
            scitb=b.segments.begin();
            for ( ; scita!=a.segments.end(); ++scita, ++scitb)
            {
              curve_type c;

              c.sum( scita->second, scitb->second );

              push_back( c, a.get_delta_t(scita) );
            }
          }

          onedpiecewisecurve sumcompcurve()
          {
            onedpiecewisecurve retcurve;

            retcurve.set_t0( get_t0() );

            typename segment_collection_type::const_iterator scit;
            for ( scit = segments.begin(); scit!=segments.end(); ++scit)
            {
              typename curve_type::onedbezcurve c;

              c = scit->second.sumcompcurve();

              retcurve.push_back( c, get_delta_t(scit) );
            }

            return retcurve;
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

          template<template<typename, unsigned short, typename> class curve1__, typename data1__, unsigned short dim1__, typename tol1__>
          friend typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type
                          eli::geom::intersect::minimum_distance(
                                  typename piecewise<curve1__, data1__, dim1__, tol1__>::data_type &t,
                                  const piecewise<curve1__, data1__, dim1__, tol1__> &pc,
                                  const typename piecewise<curve1__, data1__, dim1__, tol1__>::point_type &pt);

          typedef std::map<data_type, curve_type> segment_collection_type;

          segment_collection_type segments;
          data_type tmax;
          tolerance_type tol;

        private:
          template<typename it__>
          error_code split_seg(it__ it, const data_type &tt)
          {
            it__ itinsert;
            return split_seg(it, itinsert, tt);
          }

          template<typename it__>
          error_code split_seg(it__ it, it__ &itinsert, const data_type &tt)
          {
            // split the segment and replace
            data_type tr;
            curve_type cl, cr;

            it__ itnext = it;
            itnext++;

            it->second.split(cl, cr, tt);
            it->second = cl;

            data_type delta_t = get_delta_t(it);

            tr = it->first + delta_t*tt;

            itinsert = segments.insert(itnext, std::make_pair(tr, cr));

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERRORS;
          }

          template<typename it__>
          void segment_to_cubic(it__ it, const data_type &ttol)
          {
            curve_type c = it->second;
            curve_type cc(c);
            cc.degree_to_cubic();

            data_type d = c.eqp_distance_bound(cc);

            if(d<ttol)
            {
              // replace c with cc.
              it->second = cc;
            }
            else
            {
              // split and recurse
              it__ itl;
              split_seg(it, itl, static_cast<data_type>(0.5));

              segment_to_cubic(itl, ttol);
              segment_to_cubic(it, ttol);
            }
          }

          bool check_continuity(const eli::geom::general::continuity &cont) const
          {
            assert(!segments.empty());

            typename segment_collection_type::const_iterator it(segments.begin()), itp(it);

            for (++it; it!=segments.end(); ++it, ++itp)
            {
              data_type dt, dtp;
              dt = get_delta_t(it);
              dtp = get_delta_t(itp);
              if (!eli::geom::utility::check_point_continuity(itp->second, dtp, it->second, dt, cont, tol))
              {
                return false;
              }
            }

            return true;
          }

          data_type get_delta_t(const typename segment_collection_type::iterator &it) const
          {
            assert (it != segments.end());

            typename segment_collection_type::iterator itnext = it;
            itnext++;

            data_type delta_t;

            if(itnext != segments.end())
              delta_t = itnext->first - it->first;
            else
              delta_t = tmax - it->first;

            return delta_t;
          }

          data_type get_delta_t(const typename segment_collection_type::const_iterator &it) const
          {
            assert (it != segments.end());

            typename segment_collection_type::const_iterator itnext = it;
            itnext++;

            data_type delta_t;

            if(itnext != segments.end())
              delta_t = itnext->first - it->first;
            else
              delta_t = tmax - it->first;

            return delta_t;
          }

          data_type get_delta_t(const typename segment_collection_type::reverse_iterator &it) const
          {
            assert (it != segments.rend());

            data_type delta_t;

            if(it != segments.rbegin())
            {
              typename segment_collection_type::reverse_iterator itprev = it;
              itprev--;
              delta_t = itprev->first - it->first;
            }
            else
            {
              delta_t = tmax - it->first;
            }

            return delta_t;
          }

          data_type get_delta_t(const typename segment_collection_type::const_reverse_iterator &it) const
          {
            assert (it != segments.rend());

            data_type delta_t;

            if(it != segments.rbegin())
            {
              typename segment_collection_type::reverse_iterator itprev = it;
              itprev--;
              delta_t = itprev->first - it->first;
            }
            else
            {
              delta_t = tmax - it->first;
            }

            return delta_t;
          }

          void find_segment(typename segment_collection_type::const_iterator &it, const index_type &index) const
          {
              // advance to desired index
              index_type i;
              for (i=0, it=segments.begin(); i<index; ++i, ++it) {}
          }

          void find_segment(typename segment_collection_type::iterator &it, const index_type &index)
          {
              // advance to desired index
              index_type i;
              for (i=0, it=segments.begin(); i<index; ++i, ++it) {}
          }

          void find_segment(typename segment_collection_type::const_iterator &it, data_type &tt, const data_type &t_in) const
          {
            tol__ tol;

            // handle the end of the piecewise curve specially
            if(tol.approximately_equal(t_in, tmax))
            {
              tt=static_cast<data_type>(1);
              it=segments.end();
              it--;
              return;
            }
            if(t_in>tmax)
            {
              tt=static_cast<data_type>(2);
              it=segments.end();
              return;
            }

            // catch cases that are before the beginning of the piecewise curve
            if(t_in<get_parameter_min())
            {
              tt=static_cast<data_type>(-1);
              it=segments.end();
              return;
            }

            // Use map::upper_bound for fast lookup of segment after t_in
            it=segments.upper_bound(t_in);

            // Decrement to segment containing t_in
            if(it != segments.begin())
            {
              it--;
            }

            // At start of segment
            if(tol.approximately_equal(t_in, it->first))
            {
              tt=static_cast<data_type>(0);
              return;
            }

            data_type delta_t = get_delta_t(it);

            // At end of segment
            if(tol.approximately_equal(t_in, it->first + delta_t))
            {
              tt=static_cast<data_type>(1);
              return;
            }

            // Typical case
            tt=(t_in-it->first)/delta_t;

            // Super careful checks
            if (tt>static_cast<data_type>(1))
            {
              tt=static_cast<data_type>(1);
            }
            if (tt<static_cast<data_type>(0))
            {
              tt=static_cast<data_type>(0);
            }
          }

          void find_segment(typename segment_collection_type::iterator &it, data_type &tt, const data_type &t_in)
          {
            tol__ tol;

            // handle the end of the piecewise curve specially
            if(tol.approximately_equal(t_in, tmax))
            {
              tt=static_cast<data_type>(1);
              it=segments.end();
              it--;
              return;
            }
            if(t_in>tmax)
            {
              tt=static_cast<data_type>(2);
              it=segments.end();
              return;
            }

            // catch cases that are before the beginning of the piecewise curve
            if(t_in<get_parameter_min())
            {
              tt=static_cast<data_type>(-1);
              it=segments.end();
              return;
            }

            // Use map::upper_bound for fast lookup of segment after t_in
            it=segments.upper_bound(t_in);

            // Decrement to segment containing t_in
            if(it != segments.begin())
              it--;

            // At start of segment
            if(tol.approximately_equal(t_in, it->first))
            {
              tt=static_cast<data_type>(0);
              return;
            }

            data_type delta_t = get_delta_t(it);

            // At end of segment
            if(tol.approximately_equal(t_in, it->first + delta_t))
            {
              tt=static_cast<data_type>(1);
              return;
            }

            // Typical case
            tt=(t_in-it->first)/delta_t;

            // Super careful checks
            if (tt>static_cast<data_type>(1))
            {
              tt=static_cast<data_type>(1);
            }
            if (tt<static_cast<data_type>(0))
            {
              tt=static_cast<data_type>(0);
            }

          }

          error_code replace_it(const piecewise<curve__, data_type, dim__> &p, typename segment_collection_type::iterator &scit)
          {
            typename segment_collection_type::iterator scito;

            // Find parameter span to insert
            data_type pt0(p.get_parameter_min()), ptmax(p.get_parameter_max()), ptspan(ptmax - pt0);

            // get the first and last curve to insert
            typename segment_collection_type::const_iterator itps, itpe;
            curve_type cs, ce;
            data_type dts, dte;

            itps = p.segments.begin();
            cs = itps->second;
            dts = p.get_delta_t(itps);

            itpe = p.segments.end(); --itpe;
            ce = itpe->second;
            dte = p.get_delta_t(itpe);

            // Find parameter span of segment to replace
            data_type dti, dto;
            dti = get_delta_t(scit);

            // Ratio of parameter lengths, replace/insert
            data_type pratio = dti/ptspan;

            // check the connectivity on adjacent nodes (if available)
            if (scit != segments.begin())
            {
              scito=scit;
              --scito;
              dto = get_delta_t(scito);
              if (!eli::geom::utility::check_point_continuity(scito->second, dto, cs, dts*pratio, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (scit != segments.end() )
            {
              scito=scit;
              ++scito;
              if (scito != segments.end() )
              {
                dto = get_delta_t(scito);
                if (!eli::geom::utility::check_point_continuity(ce, dte*pratio, scito->second, dto, eli::geom::general::C0, tol))
                {
                  return SEGMENT_NOT_CONNECTED;
                }
              }
            }

            typename segment_collection_type::const_iterator it=itps;
            typename segment_collection_type::iterator itguess = scit;

            data_type dtp = dts;

            data_type t = scit->first;

            // Substitute first curve
            scit->second = it->second;

            t += dtp*pratio;

            for(++it; it != p.segments.end(); ++it)
            {
              itguess = segments.insert(itguess, std::make_pair(t, it->second));
              dtp = p.get_delta_t(it);
              t += dtp*pratio;
            }

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERRORS;
          }

      };
    }
  }
}
#endif
