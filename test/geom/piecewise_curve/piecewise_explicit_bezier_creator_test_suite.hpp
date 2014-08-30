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

#ifndef piecewise_explicit_bezier_creator_test_suite_hpp
#define piecewise_explicit_bezier_creator_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"
#include "eli/mutil/fd/d1o2.hpp"
#include "eli/mutil/fd/d2o2.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_explicit_bezier_creator.hpp"

#include "octave_helpers.hpp"

template<typename data__>
class piecewise_explicit_bezier_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
#if 0
    typedef eli::geom::curve::piecewise_explicit_bezier_creator<data__, 3, tolerance_type> explicit_bezier_type;
#endif

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_explicit_bezier_creator_test_suite<float>::create_curve_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_explicit_bezier_creator_test_suite<double>::create_curve_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_explicit_bezier_creator_test_suite<long double>::create_curve_test);
    }

  public:
    piecewise_explicit_bezier_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_explicit_bezier_creator_test_suite()
    {
    }

  private:

    void create_curve_test()
    {
#if 0
      point_creator_type af;

      data_type th, cam, cam_loc;
      bool rtn;

      // set airfoil thickness
      th=24;
      rtn=af.set_thickness(th);
      TEST_ASSERT(rtn);
      TEST_ASSERT(af.get_thickness()==th);

      // set airfoil camber
      cam=2;
      cam_loc=3;
      rtn=af.set_camber(cam, cam_loc);
      TEST_ASSERT(rtn);
      TEST_ASSERT(af.get_maximum_camber()==cam);
      TEST_ASSERT(af.get_maximum_camber_location()==cam_loc);

      // test the name
      std::string name, name_ref;

      af.set_sharp_trailing_edge(true);

      name_ref="NACA "+std::to_string(static_cast<int>(std::round(cam)))
                      +std::to_string(static_cast<int>(std::round(cam_loc)))
                      +std::to_string(static_cast<int>(std::round(th)));
      name=af.get_name();
      TEST_ASSERT(name==name_ref);


      piecewise_curve_type af_pwc;

      af.create(af_pwc);

//      octave_print(1, af_pwc);
#endif
    }
};

#endif
