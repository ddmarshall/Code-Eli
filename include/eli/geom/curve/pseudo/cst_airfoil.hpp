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
      template<typename data__, typename tol__=eli::util::tolerance<data__> >
      class cst_airfoil : public cst_base<data__, tol__>
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
          cst_airfoil() : cst_base<data_type, tolerance_type> (0.5, 1, 1) {}
          cst_airfoil(index_type dim) : cst_base<data_type, tolerance_type>(0.5, 1, dim) {}
          cst_airfoil(const cst_airfoil<data_type, tolerance_type> &cst) : cst_base<data_type, tolerance_type>(cst) {}
          ~cst_airfoil() {}

          bool operator==(const cst_airfoil<data_type, tolerance_type> &cst) const
          {
            return base_class_type::operator==(cst);
          }

          bool operator!=(const cst_airfoil<data_type, tolerance_type> &cst) const
          {
            return base_class_type::operator!=(cst);
          }

          cst_airfoil & operator=(const cst_airfoil<data_type, tolerance_type> &cst)
          {
            base_class_type::operator=(cst);
            return (*this);
          }
      };
    }
  }
}
#endif
