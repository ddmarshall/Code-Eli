/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - implementation of binary cubic creator
********************************************************************************/

#ifndef eli_geom_curve_piecewise_binary_cubic_creator_hpp
#define eli_geom_curve_piecewise_binary_cubic_creator_hpp

#include <iterator>
#include <vector>

#include "eli/code_eli.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/point/distance.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_binary_cubic_creator
      {
        public:
          typedef data__  data_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef typename point_type::Index index_type;
          typedef tol__ tolerance_type;

          typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
          typedef typename piecewise_curve_type::curve_type curve_type;

          piecewise_binary_cubic_creator() {}

          void setup(const piecewise_curve_type &pc, const data_type &t, const index_type &dmin, const index_type &dmax)
          {
            parent_curve = pc;
            parent_curve.get_pmap( parent_pmap );

            ttol = t;
            atol = -1.0;

            min_depth = dmin;
            max_depth = dmax;
          }

          void setup(const piecewise_curve_type &pc, const data_type &t, const data_type &a, const index_type &dmin, const index_type &dmax)
          {
            setup( pc, t, dmin, dmax );
            atol = a;
          }

          virtual bool corner_create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            std::vector<data_type> tdisc;
            point_type p0, m01, m02, p1, m11, m12;

            parent_curve.find_discontinuities( atol, tdisc );

            data_type t0, tmax, t1;
            t0 = parent_curve.get_t0();
            tmax = parent_curve.get_tmax();

            tdisc.push_back( tmax );

            pc.clear();

            // set the start parameter
            pc.set_t0( t0 );

            p0 = parent_curve.f(t0);
            parent_curve.fps(t0, m01, m02);

            for ( typename std::vector< data_type>::size_type i = 0; i < tdisc.size(); i++ )
            {
              t1 = tdisc[i];
              p1 = parent_curve.f(t1);
              parent_curve.fps(t1, m11, m12);

              // Build approximate curve.
              curve_type c;
              c = make_curve_point_slope(p0, m02, p1, m11, t1-t0);

              pc.push_back(c, t1-t0);

              adapt_pc( pc, t0, p0, m02, t1, p1, m11);

              t0 = t1;
              p0 = p1;
              m01 = m11;
              m02 = m12;
            }
            return true;
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef typename piecewise_curve_type::error_code error_code;

            data_type t0, t1;
            t0 = parent_curve.get_t0();
            t1 = parent_curve.get_tmax();

            point_type p0, m0, p1, m1;
            p0 = parent_curve.f(t0);
            m0 = parent_curve.fp(t0);
            p1 = parent_curve.f(t1);
            m1 = parent_curve.fp(t1);

            pc.clear();

            // set the start parameter
            pc.set_t0( t0 );

            // Build approximate curve.
            curve_type c;
            c = make_curve_point_slope(p0, m0, p1, m1, t1-t0);

            pc.push_back(c, t1-t0);

            adapt_pc( pc, t0, p0, m0, t1, p1, m1);

            return true;
          }

        protected:

          void adapt_pc( piecewise<bezier, data_type, dim__, tolerance_type> &pc, const data_type &t0, const point_type &p0, const point_type &m0,
                         const data_type &t1, const point_type &p1, const point_type &m1, index_type depth = 0 ) const
          {
            data_type tmid = ( t0 + t1 ) / 2.0;
            point_type pmid = parent_curve.f( tmid );

            bool adapt = false;

            if ( depth < min_depth )
            {
              adapt = true;
            }
            else if ( eli::geom::point::distance( pmid, pc.f( tmid ) ) > ttol )
            {
              adapt = true;
            }
            else if ( !test_match( pc, t0, t1 ) )
            {
              adapt = true;
            }

            if ( adapt && depth < max_depth )
            {
              point_type mmid = parent_curve.fp( tmid );

              piecewise<bezier, data_type, dim__, tolerance_type> insert;

              insert.set_t0( t0 );

              curve_type c1;
              c1 = make_curve_point_slope(p0, m0, pmid, mmid, tmid-t0);
              insert.push_back(c1, tmid-t0);

              curve_type c2;
              c2 = make_curve_point_slope(pmid, mmid, p1, m1, t1-tmid);
              insert.push_back(c2, t1-tmid);

              pc.replace_t( insert, t0 );

              adapt_pc( pc, t0, p0, m0, tmid, pmid, mmid, depth + 1 );
              adapt_pc( pc, tmid, pmid, mmid, t1, p1, m1, depth + 1 );
            }
          }

          bool test_match( const piecewise<bezier, data_type, dim__, tolerance_type> &pc, const data_type &tstart, const data_type &tend ) const
          {
            data_type tcheck;

            // Check start/end and midpoints of piecewise parent curve that overlap range.
            for ( typename std::vector< data_type>::size_type i = 0; i < parent_pmap.size() - 1; i++ )
            {
              if ( parent_pmap[ i ] <= tend && parent_pmap[ i + 1 ] >= tstart )
              {
                tcheck = parent_pmap[ i ];
                if ( tcheck >= tstart && tcheck <= tend )
                {
                  if ( eli::geom::point::distance( parent_curve.f( tcheck ), pc.f( tcheck ) ) > ttol )
                  {
                    return false;
                  }
                }

                tcheck = ( parent_pmap[ i ] + parent_pmap[ i + 1 ] ) / 2.0;
                if ( tcheck >= tstart && tcheck <= tend )
                {
                  if ( eli::geom::point::distance( parent_curve.f( tcheck ), pc.f( tcheck ) ) > ttol )
                  {
                    return false;
                  }
                }
              }
            }

            return true;
          }

          static curve_type make_curve_point_slope(const point_type &p0, const point_type &m0,
                                                   const point_type &p1, const point_type &m1, const data_type &dt)
          {
            curve_type c(3);
            point_type cp[4];

            cp[0]=p0;
            cp[1]=p0+(dt*m0/3.0);
            cp[2]=p1-(dt*m1/3.0);
            cp[3]=p1;

            for (index_type i=0; i<4; ++i)
            {
              c.set_control_point(cp[i], i);
            }

            return c;
          }

        private:
          piecewise_curve_type parent_curve;

          std::vector < data_type > parent_pmap;

          data_type ttol;
          data_type atol;

          index_type min_depth;
          index_type max_depth;

      };
    }
  }
}
#endif
