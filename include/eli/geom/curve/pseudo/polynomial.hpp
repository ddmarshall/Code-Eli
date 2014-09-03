/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_curve_pseudo_polynomial_hpp
#define eli_geom_curve_pseudo_polynomial_hpp

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/mutil/poly.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace pseudo
      {
        template<typename data__, unsigned short dim__, typename tol__=eli::util::tolerance<data__> >
        class polynomial
        {
          public:
            typedef unsigned short dimension_type;
            typedef eli::mutil::poly::polynomial<data__> polynomial_type;
            typedef typename polynomial_type::data_type data_type;
            typedef typename polynomial_type::index_type index_type;
            typedef Eigen::Matrix<data_type, 1, dim__> point_type;
            typedef tol__ tolerance_type;
            typedef typename polynomial_type::coefficient_type coefficient_type;

          public:
            polynomial()
            {
              for (dimension_type i=0; i<dim__; ++i)
              {
                poly[i].set_roots(0);
              }
            }
            polynomial(const polynomial<data_type, dim__, tolerance_type> &pol)
            {
              coefficient_type a;
              for (dimension_type i=0; i<dim__; ++i)
              {
                pol.poly[i].get_coefficients(a);
                poly[i].set_coefficients(a);
              }
            }
            ~polynomial() {}

            bool operator==(const polynomial<data_type, dim__, tolerance_type> &pol) const
            {
              if (this==&pol)
                return true;
              for (dimension_type i=0; i<dim__; ++i)
              {
                if (pol.poly[i]!=poly[i])
                {
                  return false;
                }
              }
              return true;
            }

            bool operator!=(const polynomial<data_type, dim__, tolerance_type> &pol) const
            {
              return !operator==(pol);
            }

            static dimension_type dimension() {return dim__;}

            index_type degree(dimension_type i) const
            {
              if (i>=i)
                return -1;

              return poly[i].degree();
            }

            void set_coefficients(const coefficient_type &a, const dimension_type &i)
            {
              if (i<dim__)
              {
                poly[i].set_coefficients(a);
              }
            }

            void get_coefficients(coefficient_type &a, const dimension_type &i) const
            {
              if (i<dim__)
              {
                poly[i].get_coefficients(a);
              }
            }

            point_type f(const data_type &t) const
            {
              point_type rtn;

              for (dimension_type i=0; i<dim__; ++i)
              {
                rtn(0, i)=poly[i].f(t);
              }

              return rtn;
            }

            point_type fp(const data_type &t) const
            {
              point_type rtn;

              for (dimension_type i=0; i<dim__; ++i)
              {
                rtn(0, i)=poly[i].fp(t);
              }

              return rtn;
            }

            point_type fpp(const data_type &t) const
            {
              point_type rtn;

              for (dimension_type i=0; i<dim__; ++i)
              {
                rtn(0, i)=poly[i].fpp(t);
              }

              return rtn;
            }

            point_type fppp(const data_type &t) const
            {
              point_type rtn;

              for (dimension_type i=0; i<dim__; ++i)
              {
                rtn(0, i)=poly[i].fppp(t);
              }

              return rtn;
            }

            point_type tangent(const data_type &t) const
            {
              point_type tgt(fp(t));

              tgt.normalize();
              return tgt;
            }

          private:
            polynomial_type poly[dim__];
        };
      }
    }
  }
}
#endif
