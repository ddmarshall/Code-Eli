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

#ifndef eli_geom_curve_piecewise_spline_creator_hpp
#define eli_geom_curve_piecewise_spline_creator_hpp

#include <iterator>
#include <vector>

#include "eli/mutil/fd/d1o2.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_cubic_spline_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef data__  data_type;
          typedef int index_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef tol__ tolerance_type;

          piecewise_cubic_spline_creator() : piecewise_creator_base<data_type, dim__, tolerance_type>(0, 0), control_point(0) {}
          piecewise_cubic_spline_creator(const index_type &ns)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(ns, 0), control_point(3*ns+1) {}
          piecewise_cubic_spline_creator(const piecewise_circle_creator<data_type, dim__, tolerance_type> &pcc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(pcc), control_point(pcc.control_point) {}

          void get_segment_control_points(point_type &cp0, point_type &cp1,
                                          point_type &cp2, point_type &cp3, const index_type &i) const
          {
            if ((3*i+1)<static_cast<index_type>(control_point.size()))
            {
              cp0=control_point[3*i];
              cp1=control_point[3*i+1];
              cp2=control_point[3*i+2];
              cp3=control_point[3*i+3];
            }
          }

          void set_segment_control_points(const point_type &cp0, const point_type &cp1,
                                          const point_type &cp2, const point_type &cp3, const index_type &i)
          {
            if ((3*i+1)<static_cast<index_type>(control_point.size()))
            {
              control_point[3*i]  =cp0;
              control_point[3*i+1]=cp1;
              control_point[3*i+2]=cp2;
              control_point[3*i+3]=cp3;
            }
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;

            curve_type c(3);
            error_code err;
            index_type nsegs(this->get_number_segments());

            // do sanity check
            if (control_point.size()!=(3*static_cast<size_t>(nsegs)+1))
            {
              assert(false);
              return false;
            }

            // set the start parameter
            pc.set_t0(this->get_t0());

            // set each segment
            for (index_type i=0; i<nsegs; ++i)
            {
              c.set_control_point(control_point[3*i  ], 0);
              c.set_control_point(control_point[3*i+1], 1);
              c.set_control_point(control_point[3*i+2], 2);
              c.set_control_point(control_point[3*i+3], 3);
              err=pc.push_back(c, this->get_segment_dt(i));
              if (err!=piecewise_curve_type::NO_ERROR)
              {
                pc.clear();
                pc.set_t0(0);
                assert(false);
                return false;
              }
            }

            return true;
          }

          /**
           * This creates a 3rd order piecewise Bezier curve that interpolates the given
           * points with Piecewise Cubit Hermite Interpolating Polynomials. The slopes
           * at the joints are approximated via 2nd order finite differences.
           * Interior slopes are calculated using central differences. The end slopes
           * are either one-sided differences unless the end condition is C1-continuous, in which
           * case central difference is used. The resulting piecewise curves are C1 continuous.
           */
          template<typename point_it__>
          void set_chip(point_it__ itb, const eli::geom::general::continuity &end_cont)
          {
            index_type i, j, nsegs=this->get_number_segments(), npts;

            npts=nsegs;
            if (end_cont==eli::geom::general::NOT_CONNECTED)
            {
              ++npts;
            }

            // can't work with less than three points
            if (npts<3)
            {
              assert(false);
              return;
            }

            point_it__ it, itm1, itp1, ite;
            data_type tmp[3], t[3], dt;
            point_type m[2];
            eli::mutil::fd::d1o2<data_type> d1approx;

            // set the end iterator if needed
            ite=itb;
            switch(end_cont)
            {
              case(eli::geom::general::NOT_CONNECTED):
              {
                break;
              }
              case(eli::geom::general::C0):
              case(eli::geom::general::C1):
              case(eli::geom::general::G1):
              {
                std::advance(ite, npts);
                break;
              }
              default:
              {
                assert(false);
                return;
                break;
              }
            }

            // need to do first segment separately
            // this is the mapping between iterators and indexes
            // itm1 -> 0
            // it ---> 1
            // itp1 -> 2
            itm1=itb;
            it=itm1; ++it;
            itp1=it; ++itp1;
            i=0;

            // calculate the slope at the start of curve
            switch (end_cont)
            {
              case(eli::geom::general::NOT_CONNECTED):
              case(eli::geom::general::C0):
              {
                // set the parameter values
                t[0]=this->get_t0();
                t[1]=t[0]+this->get_segment_dt(0);
                t[2]=t[1]+this->get_segment_dt(1);

                d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::RIGHT);
                for (j=0; j<dim__; ++j)
                {
                  tmp[0]=(*itm1)(j);
                  tmp[1]=(*it)(j);
                  tmp[2]=(*itp1)(j);
                  d1approx.evaluate(m[0](j), tmp, t);
                }

                break;
              }
              case(eli::geom::general::C1):
              case(eli::geom::general::G1):
              {
                point_it__ item1(ite);
                --item1;

                // set the parameter values
                t[0]=this->get_t0()-this->get_segment_dt(nsegs-1);
                t[1]=this->get_t0();
                t[2]=t[1]+this->get_segment_dt(0);

                d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::CENTER);
                for (j=0; j<dim__; ++j)
                {
                  tmp[0]=(*item1)(j);
                  tmp[1]=(*itm1)(j);
                  tmp[2]=(*it)(j);
                  d1approx.evaluate(m[0](j), tmp, t);
                }
                break;
              }
              default:
              {
                return;
                break;
              }
            }

            // calculate the slope at end of first segment
            t[0]=this->get_t0();
            t[1]=t[0]+this->get_segment_dt(0);
            t[2]=t[1]+this->get_segment_dt(1);
            d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::CENTER);
            for (j=0; j<dim__; ++j)
            {
              tmp[0]=(*itm1)(j);
              tmp[1]=(*it)(j);
              tmp[2]=(*itp1)(j);
              d1approx.evaluate(m[1](j), tmp, t);
            }

            // calculate the first set of control points
            dt=this->get_segment_dt(i);
            set_segment_control_points(*itm1, (*itm1)+dt*m[0]/3, (*it)-dt*m[1]/3, *it, i);
            m[0]=m[1];

            // do all interior segments
            // this is the mapping between iterators and indexes
            // itm1 -> i-1
            // it ---> i
            // itp1 -> i+1
            for (i=1; i<npts-2; ++i, ++itm1, ++it, ++itp1)
            {
              t[0]=t[1];
              t[1]=t[2];
              t[2]+=this->get_segment_dt(i+1);

              for (j=0; j<dim__; ++j)
              {
                tmp[0]=(*itm1)(j);
                tmp[1]=(*it)(j);
                tmp[2]=(*itp1)(j);
                d1approx.evaluate(m[1](j), tmp, t);
              }
              dt=this->get_segment_dt(i);
              set_segment_control_points(*it, (*it)+dt*m[0]/3, (*itp1)-dt*m[1]/3, *itp1, i);
              m[0]=m[1];
            }

            // need to do the remaining segments separately
            // this is the mapping between iterators and indexes
            // itm1 -> npts-3
            // it ---> npts-2
            // itp1 -> npts-1
            dt=this->get_segment_dt(npts-2);
            switch (end_cont)
            {
              case(eli::geom::general::NOT_CONNECTED):
              {
                // last regular segment
                d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::LEFT);
                for (j=0; j<dim__; ++j)
                {
                  tmp[0]=(*itm1)(j);
                  tmp[1]=(*it)(j);
                  tmp[2]=(*itp1)(j);
                  d1approx.evaluate(m[1](j), tmp, t);
                }

                set_segment_control_points(*it, (*it)+dt*m[0]/3, (*itp1)-dt*m[1]/3, *itp1, i);
                break;
              }
              case(eli::geom::general::C0):
              case(eli::geom::general::C1):
              case(eli::geom::general::G1):
              {
                t[0]=t[1];
                t[1]=t[2];
                t[2]+=this->get_segment_dt(i+1);

                // need to do last regular segment separately
                d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::CENTER);
                for (j=0; j<dim__; ++j)
                {
                  tmp[0]=(*it)(j);
                  tmp[1]=(*itp1)(j);
                  tmp[2]=(*itb)(j);
                  d1approx.evaluate(m[1](j), tmp, t);
                }

                set_segment_control_points(*it, (*it)+dt*m[0]/3, (*itp1)-dt*m[1]/3, *itp1, i);
                m[0]=m[1];

                // need to do closing segment separately
                dt=this->get_segment_dt(i+1);
                if (end_cont==eli::geom::general::C0)
                {
                  d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::LEFT);
                  for (j=0; j<dim__; ++j)
                  {
                    tmp[0]=(*it)(j);
                    tmp[1]=(*itp1)(j);
                    tmp[2]=(*itb)(j);
                    d1approx.evaluate(m[1](j), tmp, t);
                  }
                }
                else
                {
                  point_it__ it2(itb);
                  ++it2;

                  t[0]=t[1];
                  t[1]=t[2];
                  t[2]+=this->get_segment_dt(0);

                  d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::CENTER);
                  for (j=0; j<dim__; ++j)
                  {
                    tmp[0]=(*itp1)(j);
                    tmp[1]=(*itb)(j);
                    tmp[2]=(*it2)(j);
                    d1approx.evaluate(m[1](j), tmp, t);
                  }
                }

                set_segment_control_points(*itp1, (*itp1)+dt*m[0]/3, (*itb)-dt*m[1]/3, *itb, i+1);
                break;
              }
              default:
              {
                assert(false);
                return;
                break;
              }
            }
          }

          /**
           * This creates a 3rd order piecewise Bezier curve that interpolates the given
           * points using a cardinal spline. The cardinal spline requires a tension  parameter
           * which controls to strength of the slopes. The end slopes use the same tension
           * term, but are one-sided differences unless the end condition is C1-continuous,
           * in which the standard slope algorithm is used. The resulting piecewise
           * curves are C1 continuous.
           */
          template<typename point_it__>
          void set_cardinal(point_it__ itb, const data__ &c, const eli::geom::general::continuity &end_cont)
          {
            index_type i, nsegs=this->get_number_segments(), npts;

            npts=nsegs;
            if (end_cont==eli::geom::general::NOT_CONNECTED)
            {
              ++npts;
            }

            // can't work with less than three points
            if (npts<3)
            {
              assert(false);
              return;
            }

            // check the parameter
            if ((c<0) || (c>=1))
            {
              assert(false);
              return;
            }

            point_it__ it, itm1, itp1, ite;
            data_type dt;
            point_type m[2];

            // set the end iterator if needed
            ite=itb;
            switch(end_cont)
            {
              case(eli::geom::general::NOT_CONNECTED):
              {
                break;
              }
              case(eli::geom::general::C0):
              case(eli::geom::general::C1):
              case(eli::geom::general::G1):
              {
                std::advance(ite, npts);
                break;
              }
              default:
              {
                assert(false);
                return;
                break;
              }
            }

            // need to do first segment separately
            // this is the mapping between iterators and indexes
            // itm1 -> 0
            // it ---> 1
            // itp1 -> 2
            itm1=itb;
            it=itm1; ++it;
            itp1=it; ++itp1;
            i=0;

            // calculate the slope at the start of curve
            switch (end_cont)
            {
              case(eli::geom::general::NOT_CONNECTED):
              case(eli::geom::general::C0):
              {
                // set the parameter values
                m[0]=(1-c)*((*it)-(*itm1))/this->get_segment_dt(0);

                break;
              }
              case(eli::geom::general::C1):
              case(eli::geom::general::G1):
              {
                point_it__ item1(ite);
                --item1;

                // set the parameter values
                m[0]=(1-c)*((*it)-(*item1))/(this->get_segment_dt(nsegs-1)+this->get_segment_dt(0));

                break;
              }
              default:
              {
                return;
                break;
              }
            }

            // calculate the slope at end of first segment
            m[1]=(1-c)*((*itp1)-(*itm1))/(this->get_segment_dt(0)+this->get_segment_dt(1));

            // calculate the first set of control points
            dt=this->get_segment_dt(i);
            set_segment_control_points(*itm1, (*itm1)+dt*m[0]/3, (*it)-dt*m[1]/3, *it, i);
            m[0]=m[1];

            // do all interior segments
            // this is the mapping between iterators and indexes
            // itm1 -> i-1
            // it ---> i
            // itp1 -> i+1
            for (i=1; i<npts-2; ++i, ++itm1, ++it, ++itp1)
            {
              m[1]=(1-c)*((*itp1)-(*itm1))/(this->get_segment_dt(i)+this->get_segment_dt(i+1));

              dt=this->get_segment_dt(i);
              set_segment_control_points(*it, (*it)+dt*m[0]/3, (*itp1)-dt*m[1]/3, *itp1, i);
              m[0]=m[1];
            }

            // need to do the remaining segments separately
            // this is the mapping between iterators and indexes
            // itm1 -> npts-3
            // it ---> npts-2
            // itp1 -> npts-1
            dt=this->get_segment_dt(npts-2);
            switch (end_cont)
            {
              case(eli::geom::general::NOT_CONNECTED):
              {
                // last regular segment
                m[1]=(1-c)*((*itp1)-(*it))/this->get_segment_dt(nsegs-1);

                set_segment_control_points(*it, (*it)+dt*m[0]/3, (*itp1)-dt*m[1]/3, *itp1, i);
                break;
              }
              case(eli::geom::general::C0):
              case(eli::geom::general::C1):
              case(eli::geom::general::G1):
              {
                // need to do last regular segment separately
                m[1]=(1-c)*((*itb)-(*it))/(this->get_segment_dt(nsegs-2)+this->get_segment_dt(nsegs-1));

                set_segment_control_points(*it, (*it)+dt*m[0]/3, (*itp1)-dt*m[1]/3, *itp1, i);
                m[0]=m[1];

                // need to do closing segment separately
                dt=this->get_segment_dt(i+1);
                if (end_cont==eli::geom::general::C0)
                {
                  m[1]=(1-c)*((*itb)-(*itp1))/this->get_segment_dt(nsegs-1);
                }
                else
                {
                  point_it__ it2(itb);
                  ++it2;

                  m[1]=(1-c)*((*it2)-(*itp1))/(this->get_segment_dt(nsegs-1)+this->get_segment_dt(0));
                }

                set_segment_control_points(*itp1, (*itp1)+dt*m[0]/3, (*itb)-dt*m[1]/3, *itb, i+1);
                break;
              }
              default:
              {
                assert(false);
                return;
                break;
              }
            }
          }

          /**
           * This create a 3rd order piecewise Bezier curve that interpolates the given
           * points using a Catmull-Rom spline. This is the same as the cardinal spline
           * with the tension term set to zero.
           */
          template<typename point_it__>
          void set_catmull_rom(point_it__ itb, const eli::geom::general::continuity &end_cont)
          {
            set_cardinal(itb, 0, end_cont);
          }

          /**
           * This creates a 3rd order piecewise Bezier curve that interpolates the given
           * points using a Kochanek-Bartels spline. The Kochanek-Bartels spline requires
           * three terms (tension, bias and continuity) that control the shape of the
           * curve near the knots. The end slopes use the same tension, bias and continuity
           * terms, but only use half of the slope term unless the end condition is C1-continuous,
           * in which the standard slope algorithm is used. The resulting piecewise
           * curves are C1 continuous.
           */
          template<typename point_it__>
          void set_kochanek_bartels(point_it__ itb, const data__ &tension, const data__ &bias,
                                    const data__ &continuity, const eli::geom::general::continuity &end_cont)
          {
            index_type i, nsegs=this->get_number_segments(), npts;

            npts=nsegs;
            if (end_cont==eli::geom::general::NOT_CONNECTED)
            {
              ++npts;
            }

            // can't work with less than three points
            if (npts<4)
            {
              assert(false);
              return;
            }

            // check some parameters
            if ( (tension<-1) || (tension>1) )
            {
              assert(false);
              return;
            }
            if ( (bias<-1) || (bias>1) )
            {
              assert(false);
              return;
            }
            if ( (continuity<-1) || (continuity>1) )
            {
              assert(false);
              return;
            }

            point_it__ it, itm1, itp1, itp2, ite;
            data_type dt;
            point_type m[2];

            // set the end iterator if needed
            ite=itb;
            switch(end_cont)
            {
              case(eli::geom::general::NOT_CONNECTED):
              {
                break;
              }
              case(eli::geom::general::C0):
              case(eli::geom::general::C1):
              case(eli::geom::general::G1):
              {
                std::advance(ite, npts);
                break;
              }
              default:
              {
                assert(false);
                return;
                break;
              }
            }

            // need to do first segment separately
            // this is the mapping between iterators and indexes
            // itm1 -> 0
            // it ---> 1
            // itp1 -> 2
            // itp1 -> 3
            itm1=itb;
            it=itm1; ++it;
            itp1=it; ++itp1;
            itp2=itp1; ++itp2;
            i=0;

            // calculate the slope at the start of curve
            switch (end_cont)
            {
              case(eli::geom::general::NOT_CONNECTED):
              case(eli::geom::general::C0):
              {
                // set the parameter values
                m[0]=(1-tension)*(1-bias)*(1-continuity)*((*it)-(*itm1))/this->get_segment_dt(0);

                break;
              }
              case(eli::geom::general::C1):
              case(eli::geom::general::G1):
              {
                point_it__ item1(ite);
                --item1;

                // set the parameter values
                m[0]=0.5*(1-tension)*(1+bias)*(1+continuity)*((*itm1)-(*item1))/(-this->get_segment_dt(nsegs-1))
                    +0.5*(1-tension)*(1-bias)*(1-continuity)*((*it)-(*itm1))/this->get_segment_dt(0);

                break;
              }
              default:
              {
                return;
                break;
              }
            }

            // calculate the slope at end of first segment
            m[1]=0.5*(1-tension)*(1+bias)*(1-continuity)*((*it)-(*itm1))/this->get_segment_dt(i)
                +0.5*(1-tension)*(1-bias)*(1+continuity)*((*itp1)-(*it))/this->get_segment_dt(i+1);

            // calculate the first set of control points
            dt=this->get_segment_dt(i);
            set_segment_control_points(*itm1, (*itm1)+dt*m[0]/3, (*it)-dt*m[1]/3, *it, i);

            // do all interior segments
            // this is the mapping between iterators and indexes
            // itm1 -> i-1
            // it ---> i
            // itp1 -> i+1
            for (i=1; i<npts-2; ++i, ++itm1, ++it, ++itp1, ++itp2)
            {
              dt=this->get_segment_dt(i);
              m[0]=0.5*(1-tension)*(1+bias)*(1+continuity)*((*it)-(*itm1))/this->get_segment_dt(i-1)
                  +0.5*(1-tension)*(1-bias)*(1-continuity)*((*itp1)-(*it))/dt;
              m[1]=0.5*(1-tension)*(1+bias)*(1-continuity)*((*itp1)-(*it))/dt
                  +0.5*(1-tension)*(1-bias)*(1+continuity)*((*itp2)-(*itp1))/this->get_segment_dt(i+1);

              set_segment_control_points(*it, (*it)+dt*m[0]/3, (*itp1)-dt*m[1]/3, *itp1, i);
            }

            // need to do the remaining segments separately
            // this is the mapping between iterators and indexes
            // itm1 -> npts-3
            // it ---> npts-2
            // itp1 -> npts-1
            dt=this->get_segment_dt(npts-2);
            switch (end_cont)
            {
              case(eli::geom::general::NOT_CONNECTED):
              {
                // last regular segment
                m[0]=0.5*(1-tension)*(1+bias)*(1+continuity)*((*it)-(*itm1))/this->get_segment_dt(nsegs-2)
                    +0.5*(1-tension)*(1-bias)*(1-continuity)*((*itp1)-(*it))/this->get_segment_dt(nsegs-1);
                m[1]=(1-tension)*(1+bias)*(1-continuity)*((*itp1)-(*it))/this->get_segment_dt(nsegs-1);

                set_segment_control_points(*it, (*it)+dt*m[0]/3, (*itp1)-dt*m[1]/3, *itp1, i);
                break;
              }
              case(eli::geom::general::C0):
              case(eli::geom::general::C1):
              case(eli::geom::general::G1):
              {
                // need to do last regular segment separately
                m[0]=0.5*(1-tension)*(1+bias)*(1+continuity)*((*it)-(*itm1))/this->get_segment_dt(nsegs-3)
                    +0.5*(1-tension)*(1-bias)*(1-continuity)*((*itp1)-(*it))/this->get_segment_dt(nsegs-2);
                m[1]=0.5*(1-tension)*(1+bias)*(1-continuity)*((*itp1)-(*it))/this->get_segment_dt(nsegs-2)
                    +0.5*(1-tension)*(1-bias)*(1+continuity)*((*itb)-(*itp1))/this->get_segment_dt(nsegs-1);

                set_segment_control_points(*it, (*it)+dt*m[0]/3, (*itp1)-dt*m[1]/3, *itp1, i);

                // need to do closing segment separately
                dt=this->get_segment_dt(i+1);
                m[0]=0.5*(1-tension)*(1+bias)*(1+continuity)*((*itp1)-(*it))/this->get_segment_dt(nsegs-2)
                    +0.5*(1-tension)*(1-bias)*(1-continuity)*((*itb)-(*itp1))/this->get_segment_dt(nsegs-1);
                if (end_cont==eli::geom::general::C0)
                {
                  m[1]=(1-tension)*(1+bias)*(1-continuity)*((*itb)-(*itp1))/this->get_segment_dt(nsegs-1);
                }
                else
                {
                  point_it__ it2(itb);
                  ++it2;

                  m[1]=0.5*(1-tension)*(1+bias)*(1-continuity)*((*itb)-(*itp1))/this->get_segment_dt(nsegs-1)
                      +0.5*(1-tension)*(1-bias)*(1+continuity)*((*it2)-(*itb))/this->get_segment_dt(0);
                }

                set_segment_control_points(*itp1, (*itp1)+dt*m[0]/3, (*itb)-dt*m[1]/3, *itb, i+1);
                break;
              }
              default:
              {
                assert(false);
                return;
                break;
              }
            }
          }

          /**
           * This creates a 3rd order piecewise Bezier curve that interpolates the given
           * points enforcing C1 and C2 constraints at the knots. The not-a-knot condition
           * is used to complete the specification of the curve. The resulting piecewise
           * curves are C2 continuous.
           */
          template<typename point_it__>
          void set_cubic_spline(point_it__ /*itb*/)
          {
            // TODO: NEED TO IMPLEMENT
            assert(false);
          }

          /**
           * This creates a 3rd order piecewise Bezier curve that interpolates the given
           * points enforcing C1 and C2 constraints at the knots. The slopes are set at
           * the ends to complete the specification of the curve. The resulting piecewise
           * curves are C2 continuous.
           */
          template<typename point_it__>
          void set_clamped_cubic_spline(point_it__ /*itb*/, const point_type &/*start_slope*/, const point_type &/*end_slope*/)
          {
            // TODO: NEED TO IMPLEMENT
            assert(false);
          }

          /**
           * This creates a 3rd order piecewise Bezier curve that interpolates the given
           * points enforcing C1 and C2 constraints at the knots. The natural condition
           * is used to complete the specification of the curve. The resulting piecewise
           * curves are C2 continuous.
           */
          template<typename point_it__>
          void set_natural_cubic_spline(point_it__ /*itb*/)
          {
            // TODO: NEED TO IMPLEMENT
            assert(false);
          }

          /**
           * This creates a 3rd order piecewise Bezier curve that interpolates the given
           * points enforcing C1 and C2 constraints at the knots. The closed condition
           * (with specified smoothness) is used to complete the specification of the
           * curve. The resulting piecewise curves are C2 continuous.
           */
          template<typename point_it__>
          void set_closed_cubic_spline(point_it__ /*itb*/, const eli::geom::general::continuity &/*end_cont*/)
          {
            // TODO: NEED TO IMPLEMENT
            assert(false);
          }

          /**
           * This creates a 3rd order piecewise Bezier curve that interpolates the given
           * points enforcing C1 and C2 constraints at the knots. The periodic condition
           * (f' and f'' are the same at both ends) is used to complete the specification of the
           * curve. The resulting piecewise curves are C2 continuous.
           */
          template<typename point_it__>
          void set_periodic_cubic_spline(point_it__ /*itb*/)
          {
            // TODO: NEED TO IMPLEMENT
            assert(false);
          }

        private:
          typedef std::vector<point_type, Eigen::aligned_allocator<point_type>> point_collection_type;

          void number_segments_changed() {control_point.resize(3*this->get_number_segments()+1);}

        private:
          point_collection_type control_point;
      };
    }
  }
}
#endif
