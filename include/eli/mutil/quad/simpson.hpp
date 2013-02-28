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

#ifndef eli_mutil_quad_simpson_hpp
#define eli_mutil_quad_simpson_hpp

#include <cassert>

#include <limits>

namespace eli
{
  namespace mutil
  {
    namespace quad
    {
      template<typename data__>
      class simpson
      {
        public:
          struct adaptive_params
          {
            size_t function_count, recursion_depth, max_depth;
            data__ coarse_value, fine_value, tolerance, error_factor, tol_factor, approximate_error;

            adaptive_params()
            {
              tolerance=std::sqrt(std::numeric_limits<data__>::epsilon());
              recursion_depth=0;
              function_count=0;
              max_depth=30;
              tol_factor=1.25;
              error_factor=100;
              coarse_value=std::numeric_limits<data__>::quiet_NaN();
              fine_value=std::numeric_limits<data__>::quiet_NaN();
              approximate_error=std::numeric_limits<data__>::quiet_NaN();
            }
          };

        public:
          simpson() {};
          simpson(const simpson<data__> &) {};
          ~simpson() {};

          simpson<data__> & operator=(const simpson<data__> &) {return *this;};

          size_t order() const {return 4;}

          /** Performs quadrature for a uniform spacing of y points
           *
           *  @param yit__ - iterator type for y points
           *
           *  @param[in] dx - uniform grid spacing
           *  @param[in] yb - start iterator of y points
           *  @param[in] ye - end iterator of y points
           *
           *  @return result from quadrature
           */
          template<typename yit__>
          data__ operator()(const data__ &dx, yit__ yb, yit__ ye) const
          {
            yit__ ym1, y(yb);

            // short circuit special cases
            if (yb==ye)
              return static_cast<data__>(0);
            ++y;
            if (y==ye)
              return (*yb)/dx;
            ym1=y;
            ++y;
            if (y==ye)
              return ((*ym1)-(*yb))/dx;

            // cycle through the points until last odd point
            data__ rtnval((*yb));
            bool even_pts(true);

            while (y!=ye)
            {
              rtnval+=4*(*ym1)+2*(*y);
              ym1=y;
              ++y;
              if (y==ye)
              {
                // remove extra last odd point
                rtnval-=(*ym1);
                even_pts=false;
                break;
              }
              ym1=y;
              ++y;
            }

            if (even_pts) // have even number of points
            {
              yit__ ym2;

              // (re)set the iterators
              --y;
              --ym1;
              ym2=ym1; --ym2;

              // remove extra last odd point
              rtnval-=(*ym1);

              // add the last even segment
              rtnval+=5*(*y)/4+2*(*ym1)-(*ym2)/4;
            }

            // multiply by factor
            rtnval*=dx/static_cast<data__>(3);

            return rtnval;
          }

          /** Performs quadrature for a nonuniform spacing of y points
           *
           *  @param xit__ - iterator type for grid points
           *  @param yit__ - iterator type for y points
           *
           *  @param[in] x - start iterator of grid points
           *  @param[in] yb - start iterator of y points
           *  @param[in] ye - end iterator of y points
           *
           *  @return result from quadrature
           */
          template<typename xit__, typename yit__>
          data__ operator()(xit__ x, yit__ yb, yit__ ye) const
          {
            data__ rtnval(0), delta;
            xit__ xm1(x), xp1;
            yit__ ym1(yb), y(yb), yp1(yb);

            // set the initial iterators
            ++x; xp1=x;
            ++y; yp1=y;

            // short circuit special cases
            if (yb==ye)
              return static_cast<data__>(0);
            if (y==ye)
              return static_cast<data__>(0);
            ++xp1;
            ++yp1;
            if (yp1==ye)
              return ((*y)-(*yb))/((*x)-(*xm1));

            // cycle through the points
            bool even_pts(true);
            for (; (yp1!=ye) && (y!=ye); ++xm1, ++x, ++xp1, ++ym1, ++y, ++yp1)
            {
              delta=((*x)-(*xm1))/((*xp1)-(*x));
              rtnval+=(((*xp1)-(*x))/6)*(1+1/delta)*((delta*(2-delta))*(*yp1)+(1+delta)*(1+delta)*(*y)+(2*delta-1)*(*ym1));

              ++xm1; ++x; ++xp1;
              ++ym1; ++y; ++yp1;
              if (yp1==ye)
              {
                even_pts = false;
                break;
              }
            }

            // test if even number of point
            if (even_pts)
            {
              // reset the iterators for last segment
              --xm1; --x; --xp1;
              --ym1; --y; --yp1;

              // add last segment
              delta=((*x)-(*xm1))/((*xp1)-(*x));
              rtnval+=(((*xp1)-(*x))/6)*(((2+3*delta)/(1+delta))*(*yp1)+(3+1/delta)*(*y)-(1/(delta*(1+delta)))*(*ym1));
            }

            return rtnval;
          }

          /** Performs adaptive quadrature for a given functor
           *
           *  @param f__ - functor to be integrated
           *
           *  @param[in] x - start iterator of grid points
           *  @param[in] yb - start iterator of y points
           *  @param[in] ye - end iterator of y points
           *
           *  @return result from quadrature
           */
          template<typename f__>
          data__ operator()(f__ fun, const data__ &x0, const data__ &x1) const
          {
            typename simpson<data__>::adaptive_params ap;
            return operator()(fun, x0, x1, ap);
          }

          /** Performs adaptive quadrature for a given functor
           *
           *  @param f__ - functor to be integrated
           *
           *  @param[in] fun - instance of functor to be integrated
           *  @param[in] x0 - lower bounds of integration
           *  @param[in] x1 - upper bounds of integration
           *  @param[in] tol - error tolerance to use for convergence test
           *  @param[in] max_depth - maximum number of recursive depths to go before stopping
           *
           *  @return result from quadrature
           */
          template<typename f__>
          data__ operator()(f__ fun, const data__ &x0, const data__ &x1, const data__ &tol, const size_t &max_depth) const
          {
            typename simpson<data__>::adaptive_params ap;

            // set adaption parameters
            ap.tolerance=tol;
            ap.max_depth=max_depth;

            // do the real work
            return operator()(fun, x0, x1, ap);
          }

          /** Performs adaptive quadrature for a given functor
           *
           *  @param f__ - functor to be integrated
           *
           *  @param[in] fun - instance of functor to be integrated
           *  @param[in] x0 - lower bounds of integration
           *  @param[in] x1 - upper bounds of integration
           *  @param[in] ap - adaptive parameters
           *
           *  @return result from quadrature
           */
          template<typename f__>
          data__ operator()(f__ fun, const data__ &x0, const data__ &x1, typename simpson<data__>::adaptive_params &ap) const
          {
            data__ xc[3], fc[3];

            // create the first level evaluation of integral
            xc[0]=x0;
            fc[0]=fun(xc[0]);
            xc[1]=(x0+x1)/2;
            fc[1]=fun(xc[1]);
            xc[2]=x1;
            fc[2]=fun(xc[2]);
            ap.coarse_value=simple(x0, x1, fc);
            ap.recursion_depth=0;
            ap.function_count=3;

            // recursively subdivide the interval to get refined integral value
            if (ap.max_depth>0)
            {
              internal_recurse(fun, xc, fc, ap);
            }
            else
            {
              ap.fine_value=ap.coarse_value;
              return ap.fine_value;
            }

            // extrapolate final value of integral
            return ((1<<this->order())*ap.fine_value-ap.coarse_value)/((1<<this->order())-1);
          }

        private:
          data__ simple(const data__ &a, const data__ &b, const data__ f[]) const
          {
            return((b-a)/6)*(f[0]+4*f[1]+f[2]);
          }

          template<typename f__>
          void internal_recurse(f__ &fun, const data__ xc[], const data__ fc[], typename simpson<data__>::adaptive_params &ap) const
          {
            data__ xf[5], ff[5], sf1, sf2;

            ++ap.recursion_depth;
            // copy over and set the fine values
            // TODO: Could generalize the recursive integration if could abstract away
            //       this portion of code. Could do this with a base class that had a
            //       template parameter that is the number of points in the stencil.
            //       That way would not have to dynamically allocate any memory.
            xf[0]=xc[0];
            xf[1]=(xc[0]+xc[1])/2;
            xf[2]=xc[1];
            xf[3]=(xc[1]+xc[2])/2;
            xf[4]=xc[2];
            ff[0]=fc[0];
            ff[1]=fun(xf[1]);
            ff[2]=fc[1];
            ff[3]=fun(xf[3]);
            ff[4]=fc[2];
            ap.function_count+=2;

            // calculate the simpson rules
            sf1=simple(xf[0], xf[2], &(ff[0]));
            sf2=simple(xf[2], xf[4], &(ff[2]));
            ap.fine_value=sf1+sf2;

            // check to see if error tolerance met
            ap.approximate_error=std::abs(ap.coarse_value-ap.fine_value)/ap.error_factor;
            if ( (ap.approximate_error>ap.tolerance) && (ap.recursion_depth<ap.max_depth) )
            {
              // set up the adaption parameters based on the ones passed in
              typename simpson<data__>::adaptive_params ap1(ap), ap2(ap);
              ap1.function_count=0;
              ap1.tolerance/=ap1.tol_factor;
              ap1.coarse_value=sf1;
              ap2.function_count=0;
              ap2.tolerance/=ap2.tol_factor;
              ap2.coarse_value=sf2;

              // evaluate both sides of the interval
              internal_recurse(fun, &(xf[0]), &(ff[0]), ap1);
              internal_recurse(fun, &(xf[2]), &(ff[2]), ap2);

              // recover the adaption information
              ap.function_count+=ap1.function_count+ap2.function_count;
              ap.fine_value=ap1.fine_value+ap2.fine_value;
              ap.coarse_value=ap1.coarse_value+ap2.coarse_value;
              ap.recursion_depth=std::max(ap1.recursion_depth, ap2.recursion_depth);
              ap.approximate_error=ap1.approximate_error+ap2.approximate_error;
            }
          }
      };
    }
  }
}
#endif
