/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald - initial code and implementation
********************************************************************************/

#ifndef piecewise_binary_cubic_creator_test_suite_hpp
#define piecewise_binary_cubic_creator_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_binary_cubic_creator.hpp"
#include "eli/geom/intersect/minimum_distance_curve.hpp"

#include "eli/geom/curve/piecewise_four_digit_creator.hpp"
#include "eli/geom/curve/pseudo/four_digit.hpp"

#include "octave_helpers.hpp"

template<typename data__>
class piecewise_binary_cubic_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;

    typedef eli::geom::curve::piecewise_binary_cubic_creator<data__, 3, tolerance_type> binary_creator_type;
    typedef eli::geom::curve::piecewise_four_digit_creator<data__, 3, tolerance_type> four_digit_type;



    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_binary_cubic_creator_test_suite<float>::test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_binary_cubic_creator_test_suite<double>::test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_binary_cubic_creator_test_suite<long double>::test);
    }

  public:
    piecewise_binary_cubic_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_binary_cubic_creator_test_suite()
    {
    }

  private:

    void test()
    {
      binary_creator_type bc;
      four_digit_type af;

      data_type th, cam, cam_loc;
      bool rtn;

      // set airfoil thickness
      th=24;
      rtn=af.set_thickness(th);
      TEST_ASSERT(rtn);

      // set airfoil camber
      cam=4;
      cam_loc=2;
      rtn=af.set_camber(cam, cam_loc);
      TEST_ASSERT(rtn);

      af.set_sharp_trailing_edge(true);

      piecewise_curve_type af_pc;
      bool fit_success;

      fit_success = af.create(af_pc);
      TEST_ASSERT(fit_success);

      piecewise_curve_type bin_pc;

      bc.setup( af_pc, 0.0001, 3, 10 );

      bc.create( bin_pc );

//      if (typeid(data_type)==typeid(double))
//      {
//        af_pc.parameter_report();
//        bin_pc.parameter_report();
//
//        eli::test::octave_start(1);
//        eli::test::octave_print( 1, af_pc, "four_series");
//        eli::test::octave_print( 1, bin_pc, "binary_cubic");
//        eli::test::octave_finish(1);
//      }

    }

};

#endif

