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
      /**
       * This creates a 3rd order piecewise Bezier curve that interpolates the given
       * points with Piecewise Cubit Hermite Interpolating Polynomials that have a specified
       * parameter spacing. The slopes at the joints are approximated via 2nd order finite
       * differences. Interior slopes are calculated using central differences. The end slopes
       * are either one-sided differences unless the end condition is C1-continuous, in which
       * case central difference is used. The resulting piecewise curves are C1 continuous.
       */
      template<typename point_it__, typename t_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_chip(piecewise<bezier, data__, dim__, tol__> &pc, point_it__ itb, point_it__ ite, t_it__ titb, const eli::geom::general::continuity &cont, const data__ &tend=0)
      {
        typedef piecewise<bezier, data__, dim__, tol__> piecewise_curve_type;
        typedef typename piecewise<bezier, data__, dim__, tol__>::curve_type curve_type;
        typedef typename curve_type::data_type data_type;
        typedef typename curve_type::index_type index_type;

        index_type i, j, npts=static_cast<index_type>(std::distance(itb, ite));
        typename piecewise_curve_type::error_code err;

        if (npts<3)
          return false;

        // cycle through the points and create the segments
        point_it__ it, itm1, itp1;
        t_it__ tit, titm1, titp1;
        curve_type c(3);
        data_type dt, tmp[3];
        typename curve_type::control_point_type cp[4], m[2];
        eli::mutil::fd::d1o2<data_type> d1approx;

        // set the starting time
        pc.set_t0(*titb);

        // need to do first segment separately
        itm1=itb;
        it=itm1; ++it;
        itp1=it; ++itp1;
        titm1=titb;
        tit=titm1; ++tit;
        titp1=tit; ++titp1;
        dt=(*tit)-(*titm1);
        if (dt<=0)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }
        switch (cont)
        {
          case(eli::geom::general::NOT_CONNECTED):
          case(eli::geom::general::C0):
          {
            d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::RIGHT);
            for (j=0; j<dim__; ++j)
            {
              tmp[0]=(*itm1)(j);
              tmp[1]=(*it)(j);
              tmp[2]=(*itp1)(j);
              d1approx.evaluate(m[0](j), tmp, titm1);
            }
            break;
          }
          case(eli::geom::general::C1):
          case(eli::geom::general::G1):
          {
            point_it__ item1(ite);
            t_it__ tite(titb);
            data_type ttmp[3];
            --item1;
            std::advance(tite, npts-1);
            ttmp[0]=(*titm1)-(tend-(*tite));
            ttmp[1]=(*titm1);
            ttmp[2]=(*tit);
            d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::CENTER);
            for (j=0; j<dim__; ++j)
            {
              tmp[0]=(*item1)(j);
              tmp[1]=(*itm1)(j);
              tmp[2]=(*it)(j);
              d1approx.evaluate(m[0](j), tmp, ttmp);
            }
            break;
          }
          default:
          {
            pc.clear();
            pc.set_t0(0);
            return false;
            break;
          }
        }
        d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::CENTER);
        for (j=0; j<dim__; ++j)
        {
          tmp[0]=(*itm1)(j);
          tmp[1]=(*it)(j);
          tmp[2]=(*itp1)(j);
          d1approx.evaluate(m[1](j), tmp, titm1);
        }
        cp[0]=*itm1;
        cp[3]=*it;
        cp[1]=cp[0]+dt*m[0]/3;
        cp[2]=cp[3]-dt*m[1]/3;
        for (j=0; j<4; ++j)
        {
          c.set_control_point(cp[j], j);
        }
        err=pc.push_back(c, dt);
        if (err!=piecewise_curve_type::NO_ERROR)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }

        m[0]=m[1];

        // do all interior segments
        for (i=1; i<npts-2; ++i, ++itm1, ++it, ++itp1, ++titm1, ++tit, ++titp1)
        {
          dt=(*titp1)-(*tit);
          if (dt<=0)
          {
            pc.clear();
            pc.set_t0(0);
            return false;
          }
          for (j=0; j<dim__; ++j)
          {
            tmp[0]=(*itm1)(j);
            tmp[1]=(*it)(j);
            tmp[2]=(*itp1)(j);
            d1approx.evaluate(m[1](j), tmp, titm1);
          }
          cp[0]=*it;
          cp[3]=*itp1;
          cp[1]=cp[0]+dt*m[0]/3;
          cp[2]=cp[3]-dt*m[1]/3;
          for (j=0; j<4; ++j)
          {
            c.set_control_point(cp[j], j);
          }
          err=pc.push_back(c, dt);
          if (err!=piecewise_curve_type::NO_ERROR)
          {
            pc.clear();
            pc.set_t0(0);
            return false;
          }

          m[0]=m[1];
        }

        dt=(*titp1)-(*tit);
        if (dt<=0)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }
        switch (cont)
        {
          case(eli::geom::general::NOT_CONNECTED):
          {
            // need to do last regular segment separately
            d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::LEFT);
            for (j=0; j<dim__; ++j)
            {
              tmp[0]=(*itm1)(j);
              tmp[1]=(*it)(j);
              tmp[2]=(*itp1)(j);
              d1approx.evaluate(m[1](j), tmp, titm1);
            }
            cp[0]=*it;
            cp[3]=*itp1;

            break;
          }
          case(eli::geom::general::C0):
          case(eli::geom::general::C1):
          case(eli::geom::general::G1):
          {
            // need to do last regular segment separately
            d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::CENTER);
            for (j=0; j<dim__; ++j)
            {
              data_type ttmp[3];
              ttmp[0]=(*tit);
              ttmp[1]=(*titp1);
              ttmp[2]=tend;
              tmp[0]=(*it)(j);
              tmp[1]=(*itp1)(j);
              tmp[2]=(*itb)(j);
              d1approx.evaluate(m[1](j), tmp, ttmp);
            }
            cp[0]=*it;
            cp[3]=*itp1;
            cp[1]=cp[0]+dt*m[0]/3;
            cp[2]=cp[3]-dt*m[1]/3;
            for (j=0; j<4; ++j)
            {
              c.set_control_point(cp[j], j);
            }
            err=pc.push_back(c, dt);
            if (err!=piecewise_curve_type::NO_ERROR)
            {
              pc.clear();
              pc.set_t0(0);
              return false;
            }
            m[0]=m[1];

            // need to do closing segment separately
            dt=tend-(*titp1);
            if (dt<=0)
            {
              pc.clear();
              pc.set_t0(0);
              return false;
            }
            if (cont==eli::geom::general::C0)
            {
              d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::LEFT);
              for (j=0; j<dim__; ++j)
              {
                data_type ttmp[3];
                ttmp[0]=(*tit);
                ttmp[1]=(*titp1);
                ttmp[2]=tend;
                tmp[0]=(*it)(j);
                tmp[1]=(*itp1)(j);
                tmp[2]=(*itb)(j);
                d1approx.evaluate(m[1](j), tmp, ttmp);
              }
            }
            else
            {
              point_it__ it2(itb);
              t_it__ tit2(titb);
              data_type ttmp[3];
              ++it2;
              ++tit2;
              ttmp[0]=(*titp1);
              ttmp[1]=tend;
              ttmp[2]=tend+(*tit2)-(*titb);
              d1approx.set_stencil(eli::mutil::fd::d1o2<data_type>::CENTER);
              for (j=0; j<dim__; ++j)
              {
                tmp[0]=(*itp1)(j);
                tmp[1]=(*itb)(j);
                tmp[2]=(*it2)(j);
                d1approx.evaluate(m[1](j), tmp, ttmp);
              }
            }
            cp[0]=*itp1;
            cp[3]=*itb;

            break;
          }
          default:
          {
            pc.clear();
            pc.set_t0(0);
            return false;
            break;
          }
        }
        cp[1]=cp[0]+dt*m[0]/3;
        cp[2]=cp[3]-dt*m[1]/3;
        for (j=0; j<4; ++j)
        {
          c.set_control_point(cp[j], j);
        }
        err=pc.push_back(c, dt);
        if (err!=piecewise_curve_type::NO_ERROR)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }

        assert((cont==eli::geom::general::NOT_CONNECTED) || pc.closed());

        return true;
      }

      /**
       * This creates a 3rd order piecewise Bezier curve that interpolates the given
       * points with Piecewise Cubit Hermite Interpolating Polynomials that have a unit
       * parameter spacing. The slopes at the joints are approximated via 2nd order finite
       * differences. Interior slopes are calculated using central differences. The end slopes
       * are either one-sided differences unless the end condition is C1-continuous, in which
       * case central difference is used. The resulting piecewise curves are C1 continuous.
       */
      template<typename point_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_chip(piecewise<bezier, data__, dim__, tol__> &pc, point_it__ itb, point_it__ ite, const eli::geom::general::continuity &cont)
      {
        std::vector<data__> dt(std::distance(itb, ite));
        data__ tend;

        for (size_t i=0; i<dt.size(); ++i)
        {
          dt[i]=static_cast<data__>(i);
        }
        tend=dt[dt.size()-1]+1;

        return create_piecewise_chip(pc, itb, ite, dt.begin(), cont, tend);
      }

      /**
       * This creates a 3rd order piecewise Bezier curve that interpolates the given
       * points using a cardinal spline with a specified parameter spacing. The end slopes
       * use the same tension term, but are one-sided differences unless the end condition
       * is C1-continuous, in which the standard slope algorithm is used. The resulting
       * piecewise curves are C1 continuous.
       */
      template<typename point_it__, typename t_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_cardinal_spline(piecewise<bezier, data__, dim__, tol__> &pc, point_it__ itb, point_it__ ite, t_it__ titb, const data__ &c, const eli::geom::general::continuity &cont, const data__ &tend=0)
      {
        // check the parameter
        if ( (c<0) || (c>1) )
          return false;


        typedef piecewise<bezier, data__, dim__, tol__> piecewise_curve_type;
        typedef typename piecewise<bezier, data__, dim__, tol__>::curve_type curve_type;
        typedef typename curve_type::data_type data_type;
        typedef typename curve_type::index_type index_type;

        index_type i, j, npts=static_cast<index_type>(std::distance(itb, ite));
        typename piecewise_curve_type::error_code err;

        if (npts<3)
          return false;

        // cycle through the points and create the segments
        point_it__ it, itm1, itp1;
        t_it__ tit, titm1, titp1;
        curve_type cc(3);
        data_type dt;
        typename curve_type::control_point_type cp[4], m[2];

        // set the starting time
        pc.set_t0(*titb);

        // need to do first segment separately
        itm1=itb;
        it=itm1; ++it;
        itp1=it; ++itp1;
        titm1=titb;
        tit=titm1; ++tit;
        titp1=tit; ++titp1;
        dt=(*tit)-(*titm1);
        if (dt<=0)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }
        switch (cont)
        {
          case(eli::geom::general::NOT_CONNECTED):
          case(eli::geom::general::C0):
          {
            m[0]=(1-c)*((*it)-(*itm1));
            break;
          }
          case(eli::geom::general::C1):
          case(eli::geom::general::G1):
          {
            point_it__ item1(ite);
            --item1;
            m[0]=(1-c)*((*it)-(*item1));
            break;
          }
          default:
          {
            pc.clear();
            pc.set_t0(0);
            return false;
            break;
          }
        }
        m[1]=(1-c)*((*itp1)-(*itm1));
        cp[0]=*itm1;
        cp[3]=*it;
        cp[1]=cp[0]+dt*m[0]/3;
        cp[2]=cp[3]-dt*m[1]/3;
        for (j=0; j<4; ++j)
        {
          cc.set_control_point(cp[j], j);
        }
        err=pc.push_back(cc, dt);
        if (err!=piecewise_curve_type::NO_ERROR)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }
        m[0]=m[1];

        // do all interior segments
        for (i=1; i<npts-2; ++i, ++itm1, ++it, ++itp1, ++titm1, ++tit, ++titp1)
        {
          dt=(*titp1)-(*tit);
          if (dt<=0)
          {
            pc.clear();
            pc.set_t0(0);
            return false;
          }
          m[1]=(1-c)*((*itp1)-(*itm1));
          cp[0]=*it;
          cp[3]=*itp1;
          cp[1]=cp[0]+dt*m[0]/3;
          cp[2]=cp[3]-dt*m[1]/3;
          for (j=0; j<4; ++j)
          {
            cc.set_control_point(cp[j], j);
          }
          err=pc.push_back(cc, dt);
          if (err!=piecewise_curve_type::NO_ERROR)
          {
            pc.clear();
            pc.set_t0(0);
            return false;
          }
          m[0]=m[1];
        }

        dt=(*titp1)-(*tit);
        if (dt<=0)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }
        switch (cont)
        {
          case(eli::geom::general::NOT_CONNECTED):
          {
            // need to do last regular segment separately
            m[1]=(1-c)*((*itp1)-(*itm1));
            cp[0]=*it;
            cp[3]=*itp1;

            break;
          }
          case(eli::geom::general::C0):
          case(eli::geom::general::C1):
          case(eli::geom::general::G1):
          {
            // need to do last regular segment separately
            m[1]=(1-c)*((*itb)-(*it));
            cp[0]=*it;
            cp[3]=*itp1;
            cp[1]=cp[0]+dt*m[0]/3;
            cp[2]=cp[3]-dt*m[1]/3;
            for (j=0; j<4; ++j)
            {
              cc.set_control_point(cp[j], j);
            }
            err=pc.push_back(cc, dt);
            if (err!=piecewise_curve_type::NO_ERROR)
            {
              pc.clear();
              pc.set_t0(0);
              return false;
            }
            m[0]=m[1];

            // need to do closing segment separately
            dt=tend-(*titp1);
            if (dt<=0)
            {
              pc.clear();
              pc.set_t0(0);
              return false;
            }
            if (cont==eli::geom::general::C0)
            {
              m[1]=(1-c)*((*itb)-(*itp1));
            }
            else
            {
              point_it__ it2(itb);
              ++it2;
              m[1]=(1-c)*((*it2)-(*itp1));
            }
            cp[0]=*itp1;
            cp[3]=*itb;

            break;
          }
          default:
          {
            pc.clear();
            pc.set_t0(0);
            return false;
            break;
          }
        }
        cp[1]=cp[0]+dt*m[0]/3;
        cp[2]=cp[3]-dt*m[1]/3;
        for (j=0; j<4; ++j)
        {
          cc.set_control_point(cp[j], j);
        }
        err=pc.push_back(cc, dt);
        if (err!=piecewise_curve_type::NO_ERROR)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }

        assert((cont==eli::geom::general::NOT_CONNECTED) || pc.closed());

        return true;
      }

      /**
       * This creates a 3rd order piecewise Bezier curve that interpolates the given
       * points using a cardinal spline with a unit parameter spacing. The end slopes
       * use the same tension term, but are one-sided differences unless the end condition
       * is C1-continuous, in which the standard slope algorithm is used. The resulting
       * piecewise curves are C1 continuous.
       */
      template<typename point_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_cardinal_spline(piecewise<bezier, data__, dim__, tol__> &pc, point_it__ itb, point_it__ ite, const data__ &c, const eli::geom::general::continuity &cont)
      {
        std::vector<data__> dt(std::distance(itb, ite));
        data__ tend;

        for (size_t i=0; i<dt.size(); ++i)
        {
          dt[i]=static_cast<data__>(i);
        }
        tend=dt[dt.size()-1]+1;

        return create_piecewise_cardinal_spline(pc, itb, ite, dt.begin(), c, cont, tend);
      }

      /**
       * This create a 3rd order piecewise Bezier curve that interpolates the given
       * points using a Catmull-Rom spline with a specified parameter spacing. This is
       * the same as the cardinal spline with the tension term set to zero.
       */
      template<typename point_it__, typename t_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_catmull_rom_spline(piecewise<bezier, data__, dim__, tol__> &pc, point_it__ itb, point_it__ ite, t_it__ titb, const eli::geom::general::continuity &cont, const data__ &tend=0)
      {
        return create_piecewise_cardinal_spline(pc, itb, ite, titb, static_cast<data__>(0), cont, tend);
      }

      /**
       * This create a 3rd order piecewise Bezier curve that interpolates the given
       * points using a Catmull-Rom spline with a unit parameter spacing. This is
       * the same as the cardinal spline with the tension term set to zero.
       */
      template<typename point_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_catmull_rom_spline(piecewise<bezier, data__, dim__, tol__> &pc, point_it__ itb, point_it__ ite, const eli::geom::general::continuity &cont)
      {
        return create_piecewise_cardinal_spline(pc, itb, ite, static_cast<data__>(0), cont);
      }

      template<typename point_it__, typename t_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_kochanek_bartels_spline(piecewise<bezier, data__, dim__, tol__> &pc, point_it__ itb, point_it__ ite, t_it__ titb, const data__ &tension, const data__ &bias, const data__ &continuity, const eli::geom::general::continuity &cont, const data__ &tend=0)
      {
        // check some parameters
        if ( (tension<-1) || (tension>1) )
          return false;
        if ( (bias<-1) || (bias>1) )
          return false;
        if ( (continuity<-1) || (continuity>1) )
          return false;

        typedef piecewise<bezier, data__, dim__, tol__> piecewise_curve_type;
        typedef typename piecewise<bezier, data__, dim__, tol__>::curve_type curve_type;
        typedef typename curve_type::data_type data_type;
        typedef typename curve_type::index_type index_type;

        index_type i, j, npts=static_cast<index_type>(std::distance(itb, ite));
        typename piecewise_curve_type::error_code err;

        if (npts<3)
          return false;

        // cycle through the points and create the segments
        point_it__ it, itm1, itp1, itp2;
        t_it__ tit, titm1, titp1, titp2;
        curve_type c(3);
        data_type dt;
        typename curve_type::control_point_type cp[4], m[2];

        // set the starting time
        pc.set_t0(*titb);

        // need to do first segment separately
        itm1=itb;
        it=itm1; ++it;
        itp1=it; ++itp1;
        itp2=itp1; ++itp2;
        titm1=titb;
        tit=titm1; ++tit;
        titp1=tit; ++titp1;
        titp2=titp1; ++titp2;
        dt=(*tit)-(*titm1);
        if (dt<=0)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }
        switch (cont)
        {
          case(eli::geom::general::NOT_CONNECTED):
          case(eli::geom::general::C0):
          {
            m[0]=(1-tension)*(1-bias)*(1-continuity)*((*it)-(*itm1))/((*tit)-(*titm1));
            break;
          }
          case(eli::geom::general::C1):
          case(eli::geom::general::G1):
          {
            point_it__ item1(ite);
            t_it__ tite(titb);

            --item1;
            std::advance(tite, npts-1);
            m[0]=0.5*(1-tension)*(1+bias)*(1+continuity)*((*itm1)-(*item1))/(tend-(*tite))
                +0.5*(1-tension)*(1-bias)*(1-continuity)*((*it)-(*itm1))/((*tit)-(*titm1));
            break;
          }
          default:
          {
            pc.clear();
            pc.set_t0(0);
            return false;
            break;
          }
        }
        m[1]=0.5*(1-tension)*(1+bias)*(1-continuity)*((*it)-(*itm1))/((*tit)-(*titm1))
            +0.5*(1-tension)*(1-bias)*(1+continuity)*((*itp1)-(*it))/((*titp1)-(*tit));
        cp[0]=*itm1;
        cp[3]=*it;
        cp[1]=cp[0]+dt*m[0]/3;
        cp[2]=cp[3]-dt*m[1]/3;
        for (j=0; j<4; ++j)
        {
          c.set_control_point(cp[j], j);
        }
        err=pc.push_back(c, dt);
        if (err!=piecewise_curve_type::NO_ERROR)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }

        m[0]=m[1];

        // do all interior segments
        for (i=1; i<npts-2; ++i, ++itm1, ++it, ++itp1, ++itp2, ++titm1, ++tit, ++titp1, ++titp2)
        {
          dt=(*titp1)-(*tit);
          if (dt<=0)
          {
            pc.clear();
            pc.set_t0(0);
            return false;
          }
          m[0]=0.5*(1-tension)*(1+bias)*(1+continuity)*((*it)-(*itm1))/((*tit)-(*titm1))
              +0.5*(1-tension)*(1-bias)*(1-continuity)*((*itp1)-(*it))/((*titp1)-(*tit));
          m[1]=0.5*(1-tension)*(1+bias)*(1-continuity)*((*itp1)-(*it))/((*titp1)-(*tit))
              +0.5*(1-tension)*(1-bias)*(1+continuity)*((*itp2)-(*itp1))/((*titp2)-(*titp1));
          cp[0]=*it;
          cp[3]=*itp1;
          cp[1]=cp[0]+dt*m[0]/3;
          cp[2]=cp[3]-dt*m[1]/3;
          for (j=0; j<4; ++j)
          {
            c.set_control_point(cp[j], j);
          }
          err=pc.push_back(c, dt);
          if (err!=piecewise_curve_type::NO_ERROR)
          {
            pc.clear();
            pc.set_t0(0);
            return false;
          }
        }

        dt=(*titp1)-(*tit);
        if (dt<=0)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }
        switch (cont)
        {
          case(eli::geom::general::NOT_CONNECTED):
          {
            // need to do last regular segment separately
            m[0]=0.5*(1-tension)*(1+bias)*(1+continuity)*((*it)-(*itm1))/((*tit)-(*titm1))
                +0.5*(1-tension)*(1-bias)*(1-continuity)*((*itp1)-(*it))/((*titp1)-(*tit));
            m[1]=(1-tension)*(1+bias)*(1-continuity)*((*itp1)-(*it))/((*titp1)-(*tit));
            cp[0]=*it;
            cp[3]=*itp1;

            break;
          }
          case(eli::geom::general::C0):
          case(eli::geom::general::C1):
          case(eli::geom::general::G1):
          {
            // need to do last regular segment separately
            m[0]=0.5*(1-tension)*(1+bias)*(1+continuity)*((*it)-(*itm1))/((*tit)-(*titm1))
                +0.5*(1-tension)*(1-bias)*(1-continuity)*((*itp1)-(*it))/((*titp1)-(*tit));
            m[1]=0.5*(1-tension)*(1+bias)*(1-continuity)*((*itp1)-(*it))/((*titp1)-(*tit))
                +0.5*(1-tension)*(1-bias)*(1+continuity)*((*itb)-(*itp1))/(tend-(*titp1));
            cp[0]=*it;
            cp[3]=*itp1;
            cp[1]=cp[0]+dt*m[0]/3;
            cp[2]=cp[3]-dt*m[1]/3;
            for (j=0; j<4; ++j)
            {
              c.set_control_point(cp[j], j);
            }
            err=pc.push_back(c, dt);
            if (err!=piecewise_curve_type::NO_ERROR)
            {
              pc.clear();
              pc.set_t0(0);
              return false;
            }

            // need to do closing segment separately
            dt=tend-(*titp1);
            if (dt<=0)
            {
              pc.clear();
              pc.set_t0(0);
              return false;
            }
            m[0]=0.5*(1-tension)*(1+bias)*(1+continuity)*((*itp1)-(*it))/((*titp1)-(*tit))
                +0.5*(1-tension)*(1-bias)*(1-continuity)*((*itb)-(*itp1))/(tend-(*titp1));
            if (cont==eli::geom::general::C0)
            {
              m[1]=(1-tension)*(1+bias)*(1-continuity)*((*itb)-(*itp1))/(tend-(*titp1));
            }
            else
            {
              point_it__ it2(itb);
              t_it__ tit2(titb);

              ++it2;
              ++tit2;
              m[1]=0.5*(1-tension)*(1+bias)*(1-continuity)*((*itb)-(*itp1))/(tend-(*titp1))
                  +0.5*(1-tension)*(1-bias)*(1+continuity)*((*it2)-(*itb))/((*tit2)-(*titb));
            }
            cp[0]=*itp1;
            cp[3]=*itb;

            break;
          }
          default:
          {
            pc.clear();
            pc.set_t0(0);
            return false;
            break;
          }
        }
        cp[1]=cp[0]+dt*m[0]/3;
        cp[2]=cp[3]-dt*m[1]/3;
        for (j=0; j<4; ++j)
        {
          c.set_control_point(cp[j], j);
        }
        err=pc.push_back(c, dt);
        if (err!=piecewise_curve_type::NO_ERROR)
        {
          pc.clear();
          pc.set_t0(0);
          return false;
        }

        assert((cont==eli::geom::general::NOT_CONNECTED) || pc.closed());

        return true;
      }

      template<typename point_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_kochanek_bartels_spline(piecewise<bezier, data__, dim__, tol__> &pc, point_it__ itb, point_it__ ite, const data__ &tension, const data__ &bias, const data__ &continuity, const eli::geom::general::continuity &cont)
      {
        std::vector<data__> dt(std::distance(itb, ite));
        data__ tend;

        for (size_t i=0; i<dt.size(); ++i)
        {
          dt[i]=static_cast<data__>(i);
        }
        tend=dt[dt.size()-1]+1;

        return create_piecewise_kochanek_bartels_spline(pc, itb, ite, dt.begin(), tension, bias, continuity, cont, tend);
      }



      template<typename point_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_piecewise_cubic_spline(piecewise<bezier, data__, dim__, tol__> &/*pc*/, point_it__ /*itb*/, point_it__ /*ite*/)
      {
        // this implements the not-a-knot end condition
#if 0
        size_type npts=std::distance(itb, ite);
        std::vector<typename point_it__::T> dt(npts);
        return create_piecewise_cubic(itb, ite, dt.begin());
#else
        // NOT IMPLEMENTED
        assert(false);
        return false;
#endif
      }

      template<typename point_it__, typename data__, unsigned short dim__, typename tol__, typename point__>
      bool create_clamped_piecewise_cubic_spline(piecewise<bezier, data__, dim__, tol__> &/*pc*/, point_it__ /*itb*/, point_it__ /*ite*/, const point__ &/*start_slope*/, const point__ &/*end_slope*/)
      {
        // NOT IMPLEMENTED
        assert(false);
        return false;
      }

      template<typename point_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_natural_piecewise_cubic_spline(piecewise<bezier, data__, dim__, tol__> &/*pc*/, point_it__ /*itb*/, point_it__ /*ite*/)
      {
        // NOT IMPLEMENTED
        assert(false);
        return false;
      }

      template<typename point_it__, typename data__, unsigned short dim__, typename tol__>
      bool create_closed_natural_piecewise_cubic_spline(piecewise<bezier, data__, dim__, tol__> &/*pc*/, point_it__ /*itb*/, point_it__ /*ite*/)
      {
        // NOT IMPLEMENTED
        assert(false);
        return false;
      }
    }
  }
}
#endif
