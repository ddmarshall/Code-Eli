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

#ifndef eli_geom_curve_length_hpp
#define eli_geom_curve_length_hpp

#include "eli/mutil/quad.hpp"

#include "eli/geom/curve/piecewise.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace internal
      {
        template <typename curve__>
        struct length_functor
        {
          const curve__ *pcurve;

          typename curve__::data_type operator()(const typename curve__::data_type &t)
          {
            return pcurve->fp(t).norm();
          }
        };
      }

      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__>
      void length(typename piecewise<curve__, data__, dim__, tol__>::data_type &len,
                  const piecewise<curve__, data__, dim__, tol__> &pc,
                  const typename piecewise<curve__, data__, dim__, tol__>::data_type &tol)
      {
        typedef piecewise<curve__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::data_type data_type;

        typename piecewise_type::segment_collection_type::const_iterator it;
        data_type seg_len;
        for (len=0, it=pc.segments.begin(); it!=pc.segments.end(); ++it)
        {
          length(seg_len, it->c, tol);
          len+=seg_len;
        }
      }
      template<template<typename, unsigned short, typename> class curve__, typename data__, unsigned short dim__, typename tol__>
      void length(typename piecewise<curve__, data__, dim__, tol__>::data_type &len,
                  const piecewise<curve__, data__, dim__, tol__> &pc,
                  const typename piecewise<curve__, data__, dim__, tol__>::data_type &t0,
                  const typename piecewise<curve__, data__, dim__, tol__>::data_type &t1,
                  const typename piecewise<curve__, data__, dim__, tol__>::data_type &tol)
      {
        typedef piecewise<curve__, data__, dim__, tol__> piecewise_type;
        typedef typename piecewise_type::data_type data_type;

        // short circuit for invalid parameters
        if (t0>=t1)
        {
          len=0;
          return;
        }

        typename piecewise_type::segment_collection_type::const_iterator it, it0, it1;
        data_type tt0, tt1, seg_len;

        // calculate the length of the start and end pieces
        pc.find_segment(it0, tt0, t0);
        pc.find_segment(it1, tt1, t1);
        if (it0==it1)
        {
          length(len, it0->c, tt0, tt1, tol);
          return;
        }
        length(seg_len, it0->c, tt0, 1, tol);
        len=seg_len;
        length(seg_len, it1->c, 0, tt1, tol);
        len+=seg_len;

        // add the length of all of the complete segments between the start and end pieces
        it=it0;
        for (++it; it!=it1; ++it)
        {
          length(seg_len, it->c, tol);
          len+=seg_len;
        }
      }

      template<typename curve__>
      void length(typename curve__::data_type &len, const curve__ &c, const typename curve__::data_type &tol)
      {
        length(len, c, static_cast<typename curve__::data_type>(0.0), static_cast<typename curve__::data_type>(1.0), tol);
      }
      template<typename curve__>
      void length(typename curve__::data_type &len, const curve__ &c, const typename curve__::data_type &t0, const typename curve__::data_type &t1, const typename curve__::data_type &tol)
      {
        eli::mutil::quad::simpson<typename curve__::data_type> quad;
        typename eli::mutil::quad::simpson<typename curve__::data_type>::adaptive_params ap;
        internal::length_functor<curve__> f;

        // short circuit for invalid parameters
        if (t0>=t1)
        {
          len=0;
          return;
        }

        f.pcurve=&c;

        // set the specified tolerance
        ap.tolerance=tol;
        len=quad(f, t0, t1, ap);
      }

// NOTE: These are here as a reference implementation of an algorithm for bezier curves.
//       It was written as methods to class, so will need to be rewritten to be an external
//       function. It didn't seem to be any better than the general algorithm.
#if 0
          data_type length_adaptive_subdivision() const
          {
            // NOTE: Implements Adaptive Subdivision method: Jens Gravesen. Adaptive Subdivision and
            //       the length and energy of Bezier curves. Computational Geometry, v8, pp. 13-31, 1997.
            return internal_length_adaptive_subdivision(SMALL_POS_FLOAT, 0);
          }

          data_type internal_length_adaptive_subdivision(const data_type &tol, const index_type &depth) const
          {
            data_type L1, Lc, Lp, err;
            index_type i, n(this->degree()), depth_max(20);

            // find the length of the control polygon and the chord
            Lc=dist(b[0], b[n]);
            Lp=0;
            for (i=0; i<n; ++i)
              Lp+=dist(b[i], b[i+1]);
            L1=(2*Lc+(n-1)*Lp)/(n+1);

            // calculate error
            // NOTE: Lp-Lc is order(h^2) approximation to error, so just dividing by h^2 to get
            //       approximation to tighter error bounds
            err=std::abs(Lp-Lc)/(1<<(2*depth));

            // if error small enough, converged
            if ( (err>tol) && (depth<depth_max) )
            {
              beziern_curve bc_l, bc_r;

              // split the curve in half
              split(bc_l, bc_r, 0.5);

              // calculate the length using the two segments and halving the tolerances
              L1=bc_l.internal_length_adaptive_subdivision(tol/2, depth+1)+bc_r.internal_length_adaptive_subdivision(tol/2, depth+1);
            }

            return L1;
          }
#endif
    }
  }
}
#endif
