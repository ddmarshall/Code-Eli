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
              v1.normalize();
              v2.normalize();

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
              v1.normalize();
              v2.normalize();

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
              v1.normalize();
              v2.normalize();

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

          void get_bounding_box(bounding_box_type &bb) const
          {
            typename segment_collection_type::const_iterator it;
            bounding_box_type bb_local;

            bb.clear();

            // cycle through all segments to get each bounding box to add
            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              it->c.get_bounding_box(bb_local);
              bb.add(bb_local);
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
            typename segment_collection_type::const_iterator itlast, itfirst;

            itlast=segments.end(); --itlast;
            itfirst=segments.begin();
            return eli::geom::utility::check_point_continuity(itlast->c, itlast->delta_t, itfirst->c, itfirst->delta_t, eli::geom::general::C0, tol);
          }
          bool open() const
          {
            return !closed();
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
              if (!eli::geom::utility::check_point_continuity(curve, dt, segments.begin()->c, segments.begin()->delta_t, eli::geom::general::C0, tol))
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
              if (!eli::geom::utility::check_point_continuity(segments.rbegin()->c, segments.rbegin()->delta_t, curve, dt, eli::geom::general::C0, tol))
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
              if (!eli::geom::utility::check_point_continuity(scito->c, scito->delta_t, curve, scito->delta_t, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if ((index+1)<number_segments())
            {
              scito=scit;
              ++scito;
              if (!eli::geom::utility::check_point_continuity(curve, scito->delta_t, scito->c, scito->delta_t, eli::geom::general::C0, tol))
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
              if (!eli::geom::utility::check_point_continuity(scito->c, scito->delta_t, curve, dt, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if ((index+1)<number_segments())
            {
              scito=scit;
              ++scito;
              if (!eli::geom::utility::check_point_continuity(curve, dt, scito->c, scito->delta_t, eli::geom::general::C0, tol))
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

          /** index0 - start index of segments to replace
            * index1 - index of the first segment after the ones to be replaced
            */
          error_code replace(const curve_type &curve, const data_type &dt, const index_type &index0, const index_type &index1)
          {
            if (index0>=number_segments())
              return INVALID_INDEX;
            if (index1>number_segments())
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
              if (!eli::geom::utility::check_point_continuity(scito->c, scito->delta_t, curve, dt, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (index1<number_segments())
            {
              scito=scit1;
              if (!eli::geom::utility::check_point_continuity(curve, dt, scito->c, scito->delta_t, eli::geom::general::C0, tol))
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
            data_type dts, dte;
            p.get(cs, dts, 0);
            p.get(ce, dte, p.number_segments()-1);

            // check the connectivity on adjacent nodes (if available)
            if (index>0)
            {
              scito=scit;
              --scito;
              if (!eli::geom::utility::check_point_continuity(scito->c, scito->delta_t, cs, dts, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if ((index+1)<number_segments())
            {
              scito=scit;
              ++scito;
              if (!eli::geom::utility::check_point_continuity(ce, dte, scito->c, scito->delta_t, eli::geom::general::C0, tol))
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

            // advance to desired index
            index_type i;
            typename segment_collection_type::iterator scit0, scit1, scito;
            for (i=0, scit0=segments.begin(); i<index0; ++i, ++scit0) {}
            for (scit1=scit0; i<index1; ++i, ++scit1) {}

            // get the first and last curve
            curve_type cs, ce;
            data_type dts, dte;
            p.get(cs, dts, 0);
            p.get(ce, dte, p.number_segments()-1);

            // check the connectivity on adjacent nodes (if available)
            if (index0>0)
            {
              scito=scit0;
              --scito;
              if (!eli::geom::utility::check_point_continuity(scito->c, scito->delta_t, cs, dts, eli::geom::general::C0, tol))
              {
                return SEGMENT_NOT_CONNECTED;
              }
            }
            if (index1<number_segments())
            {
              scito=scit1;
              if (!eli::geom::utility::check_point_continuity(ce, dte, scito->c, scito->delta_t, eli::geom::general::C0, tol))
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
            data_type dtim1, dti, r, lenim1, leni, tim1_split(-1), ti_split(-1);
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
              ti_split=0;
              r=leni;
            }

            // find coordinate that corresponds to location of radius on each curve
            if (tim1_split<0)
            {
              // FIX: this is exact for straight lines and approximate for other curves
              tim1_split=1-r/lenim1;
            }
            if (ti_split<0)
            {
              // FIX: this is exact for straight lines and approximate for other curves
              ti_split=r/leni;
            }

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
            cim1.split(c, ctrim, tim1_split); cim1=c;
            ci.split(ctrim, c, ti_split); ci=c;

            // replace/add curves
            error_code ec;
            if (rounding_end)
            {
              curve_type arc1, arc2;
              data_type orig_t0(get_t0());

              // split the arc
              arc.split(arc1, arc2, static_cast<data_type>(0.5));

              // put the ith segment and arc2 onto the end of the list of curves
              ec=replace(cim1, dtim1*tim1_split, number_segments()-1);
              if (ec!=NO_ERROR)
              {
                assert(false);
                return false;
              }
              ec=push_back(arc1, dtim1*(static_cast<data_type>(1)-tim1_split));
              if (ec!=NO_ERROR)
              {
                assert(false);
                return false;
              }

              // put arc1 and the (i-1)st segment onto the front of the list of curves
              ec=replace(ci, dti*(static_cast<data_type>(1)-ti_split), 0);
              if (ec!=NO_ERROR)
              {
                assert(false);
                return false;
              }
              ec=push_front(arc2, dti*ti_split);
              if (ec!=NO_ERROR)
              {
                assert(false);
                return false;
              }
              set_t0(orig_t0);
            }
            else
            {
              piecewise<curve__, data_type, dim__> pct;

              ec=pct.push_back(cim1, dtim1*tim1_split);
              if (ec!=NO_ERROR)
              {
                assert(false);
                return false;
              }
              ec=pct.push_back(arc, dtim1*(1-tim1_split)+dti*ti_split);
              if (ec!=NO_ERROR)
              {
                assert(false);
                return false;
              }
              ec=pct.push_back(ci, dti*(1-ti_split));
              if (ec!=NO_ERROR)
              {
                assert(false);
                return false;
              }

              // replace the two segments with the piecewise curve
              ec=replace(pct, i-1, i+1);
              if (ec!=NO_ERROR)
              {
                assert(false);
                return false;
              }
            }

            return true;
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
              return eli::geom::general::NOT_CONNECTED;
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
                itfirst=it;
                itsecond=it; ++itsecond;
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

            // check the continuity of the two sections
            return eli::geom::utility::check_point_continuity(itfirst->c, itfirst->delta_t, itsecond->c, itsecond->delta_t, cont, tol);
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
                itfirst=it;
                itsecond=it; ++itsecond;
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

            // check the continuity of the two sections
            return eli::geom::utility::report_point_continuity(itfirst->c, itfirst->delta_t, itsecond->c, itsecond->delta_t, tol);
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
              if (!eli::geom::utility::check_point_continuity(itp->c, itp->delta_t, it->c, it->delta_t, cont, tol))
              {
                return false;
              }
            }

            return true;
          }

          void find_segment(typename segment_collection_type::const_iterator &it, data_type &tt, const data_type &t_in) const
          {
            tol__ tol;
            data_type t(t0);

            // check to see if have invalid t_in
            if (t_in<t0)
            {
              if (tol.approximately_equal(t_in, t0))
              {
                it=segments.begin();
                tt=0;
                return;
              }
              it=segments.end();
              return;
            }

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              if (tol.approximately_equal(t_in, t+it->delta_t))
              {
                tt=1;
                return;
              }

              if (t_in<=t+it->delta_t)
              {
                tt=(t_in-t)/it->delta_t;
                if (tt>static_cast<data_type>(1))
                  tt=static_cast<data_type>(1);
                if (tt<static_cast<data_type>(0))
                  tt=static_cast<data_type>(0);
                return;
              }
              t+=it->delta_t;
            }
          }

          void find_segment(typename segment_collection_type::iterator &it, data_type &tt, const data_type &t_in)
          {
            tol__ tol;
            data_type t(t0);

            // check to see if have invalid t_in
            if (t_in<t0)
            {
              if (tol.approximately_equal(t_in, t0))
              {
                it=segments.begin();
                tt=0;
                return;
              }
              it=segments.end();
              return;
            }

            for (it=segments.begin(); it!=segments.end(); ++it)
            {
              if (tol.approximately_equal(t_in, t+it->delta_t))
              {
                tt=1;
                return;
              }

              if (t_in<=t+it->delta_t)
              {
                tt=(t_in-t)/it->delta_t;
                if (tt>static_cast<data_type>(1))
                  tt=static_cast<data_type>(1);
                if (tt<static_cast<data_type>(0))
                  tt=static_cast<data_type>(0);
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
