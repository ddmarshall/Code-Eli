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

#ifndef eli_quad_trapezoid_hpp
#define eli_quad_trapezoid_hpp

namespace eli
{
  namespace quad
  {
    template<typename data__>
    class trapezoid
    {
      public:
        trapezoid() {};
        trapezoid(const trapezoid<data__> &) {};
        ~trapezoid() {};

        trapezoid<data__> & operator=(const trapezoid<data__> &) {return *this;};

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
          data__ rtnval(-0.5*(*yb));
          yit__ y(yb);

          for (; y!=ye; ++y)
            rtnval+=(*y);

          --y;
          rtnval-=0.5*(*y);
          rtnval*=dx;

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
          xit__ xp(x);
          yit__ y(yb), yp(yb);
          data__ rtnval(0);

          for (++x, ++y; y!=ye; ++x, ++xp, ++y, ++yp)
            rtnval+=0.5*((*x)-(*xp))*((*y)+(*yp));

          return rtnval;
        }
    };
  }
}
#endif