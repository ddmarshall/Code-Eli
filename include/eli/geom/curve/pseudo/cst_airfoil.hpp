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

#ifndef eli_geom_curve_pseudo_cst_airfoil_hpp
#define eli_geom_curve_pseudo_cst_airfoil_hpp

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/geom/curve/pseudo/cst_base.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace pseudo
      {
        template<typename data__, typename tol__=eli::util::tolerance<data__> >
        class cst_airfoil_curve : public cst_base<data__, tol__>
        {
          private:
            typedef cst_base<data__, tol__> base_class_type;

          public:
            typedef typename base_class_type::data_type data_type;
            typedef typename base_class_type::point_type point_type;
            typedef typename base_class_type::dimension_type dimension_type;
            typedef typename base_class_type::control_point_type control_point_type;
            typedef typename base_class_type::index_type index_type;
            typedef typename base_class_type::tolerance_type tolerance_type;

          public:
            cst_airfoil_curve() : cst_base<data_type, tolerance_type> (0.5, 1, 1, true) {}
            cst_airfoil_curve(index_type dim, bool u) : cst_base<data_type, tolerance_type>(0.5, 1, dim, u) {}
            cst_airfoil_curve(const cst_airfoil_curve<data_type, tolerance_type> &cst) : cst_base<data_type, tolerance_type>(cst) {}
            ~cst_airfoil_curve() {}

            bool operator==(const cst_airfoil_curve<data_type, tolerance_type> &cst) const
            {
              return base_class_type::operator==(cst);
            }

            bool operator!=(const cst_airfoil_curve<data_type, tolerance_type> &cst) const
            {
              return base_class_type::operator!=(cst);
            }

            cst_airfoil_curve & operator=(const cst_airfoil_curve<data_type, tolerance_type> &cst)
            {
              base_class_type::operator=(cst);
              return (*this);
            }
        };

        /** This encapsulates the CST airfoil algorithms into the common airfoil interface. The
         *  common parameterization for airfoils is [-1, 1] from lower surface trailing edge to
         *  upper surface trailing edge. Thus, the upper surface is the traditional [0, 1]
         *  parameterization. For the lower surface, the user specifies the CST control points
         *  in the traditional [0, 1] leading edge to trailing edge order, however the function
         *  evalations of the lower surface will occur for parameterizations of [-1, 0] from
         *  trailing edge to leading edge. One artifact of this reparameterization is that the 
         *  first derivative will be pointing in the opposite direction.
         */
        template<typename data__, typename tol__=eli::util::tolerance<data__> >
        class cst_airfoil
        {
          public:
            typedef cst_airfoil_curve<data__, tol__> airfoil_curve_type;
            typedef typename airfoil_curve_type::data_type data_type;
            typedef typename airfoil_curve_type::point_type point_type;
            typedef typename airfoil_curve_type::dimension_type dimension_type;
            typedef typename airfoil_curve_type::index_type index_type;
            typedef typename airfoil_curve_type::control_point_type control_point_type;
            typedef typename airfoil_curve_type::tolerance_type tolerance_type;
            typedef typename airfoil_curve_type::monomial_coefficient_type monomial_coefficient_type;

          public:
            cst_airfoil() {}
            cst_airfoil(index_type deg) : upper(deg, true), lower(deg, false) {}
            cst_airfoil(index_type degu, index_type degl) : upper(degu, true), lower(degl, false) {}
            cst_airfoil(const cst_airfoil<data_type, tolerance_type> &a) : upper(a.upper), lower(a.lower) {}
            ~cst_airfoil(){}

            bool operator==(const cst_airfoil<data_type, tolerance_type> &a) const
            {
              if (upper!=a.upper)
                return false;
              return (lower==a.lower);
            }

            bool operator!=(const cst_airfoil<data_type, tolerance_type> &a) const
            {
              return !operator==(a);
            }

            cst_airfoil<data_type, tolerance_type> & operator=(const cst_airfoil<data_type, tolerance_type> &a)
            {
              if (this!=&a)
              {
                upper=a.upper;
                lower=a.lower;
              }

              return (*this);
            }

            static dimension_type dimension() {return 2;}

            void resize_upper(index_type deg)
            {
              upper.resize(deg);
            }
            index_type upper_degree()
            {
              return upper.degree();
            }

            void resize_lower(index_type deg)
            {
              lower.resize(deg);
            }
            index_type lower_degree()
            {
              return lower.degree();
            }

            void upper_degree_promote()
            {
              upper.degree_promote();
            }
            void lower_degree_promote()
            {
              lower.degree_promote();
            }

            bool upper_degree_demote(const geom::general::continuity &continuity_degree=geom::general::C0)
            {
              if (continuity_degree<eli::geom::general::C0)
                return false;

              geom::general::continuity explicit_cont(eli::geom::general::C0);
              switch (continuity_degree)
              {
                case (eli::geom::general::C0):
                {
                  explicit_cont=eli::geom::general::NOT_CONNECTED;
                  break;
                }
                case (eli::geom::general::G1):
                case (eli::geom::general::C1):
                {
                  explicit_cont=eli::geom::general::C0;
                  break;
                }
                case (eli::geom::general::G2):
                case (eli::geom::general::C2):
                {
                  explicit_cont=eli::geom::general::C1;
                  break;
                }
                case (eli::geom::general::NOT_CONNECTED):
                default:
                {
                  return false;
                  break;
                }
              }
              return upper.degree_demote(explicit_cont);
            }
            bool lower_degree_demote(const geom::general::continuity &continuity_degree=geom::general::C0)
            {
              if (continuity_degree<eli::geom::general::C0)
                return false;

              geom::general::continuity explicit_cont(eli::geom::general::C0);
              switch (continuity_degree)
              {
                case (eli::geom::general::C0):
                {
                  explicit_cont=eli::geom::general::NOT_CONNECTED;
                  break;
                }
                case (eli::geom::general::G1):
                case (eli::geom::general::C1):
                {
                  explicit_cont=eli::geom::general::C0;
                  break;
                }
                case (eli::geom::general::G2):
                case (eli::geom::general::C2):
                {
                  explicit_cont=eli::geom::general::C1;
                  break;
                }
                case (eli::geom::general::NOT_CONNECTED):
                default:
                {
                  return false;
                  break;
                }
              }
              return lower.degree_demote(explicit_cont);
            }

            data_type get_t0() const {return static_cast<data_type>(-1);}
            data_type get_tmax() const {return static_cast<data_type>(1);}

            data_type get_trailing_edge_thickness() const
            {
              return upper.get_trailing_edge_thickness()+lower.get_trailing_edge_thickness();
            }
            void set_trailing_edge_thickness(const data_type &dte)
            {
              upper.set_trailing_edge_thickness(0.5*dte);
              lower.set_trailing_edge_thickness(0.5*dte);
            }

            void set_upper_control_point(const control_point_type &cp_in, const index_type &i)
            {
              upper.set_control_point(cp_in, i);
            }
            control_point_type get_upper_control_point(const index_type &i) const
            {
              return upper.get_control_point(i);
            }

            void set_lower_control_point(const control_point_type &cp_in, const index_type &i)
            {
              lower.set_control_point(cp_in, i);
            }
            control_point_type get_lower_control_point(const index_type &i) const
            {
              return lower.get_control_point(i);
            }

            void get_upper_monomial_coefficients(monomial_coefficient_type &a) const
            {
              upper.get_monomial_coefficients(a);
            }

            void get_lower_monomial_coefficients(monomial_coefficient_type &a) const
            {
              lower.get_monomial_coefficients(a);
            }

            point_type f(const data_type &t) const
            {
              assert((t>=-1) && (t<=1));

              if (t<0)
              {
                return lower.f(-t);
              }
              else
              {
                return upper.f(t);
              }
            }

            point_type fp(const data_type &t) const
            {
              assert((t>=-1) && (t<=1));

              if (t<0)
              {
                return -lower.fp(-t);
              }
              else
              {
                return upper.fp(t);
              }
            }

            point_type fpp(const data_type &t) const
            {
              assert((t>=-1) && (t<=1));

              if (t<0)
              {
                return lower.fpp(-t);
              }
              else
              {
                return upper.fpp(t);
              }
            }
          private:
            cst_airfoil_curve<data_type> upper, lower;
        };
      }
    }
  }
}
#endif
