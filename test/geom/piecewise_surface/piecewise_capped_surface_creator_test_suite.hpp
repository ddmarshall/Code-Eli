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

#ifndef piecewise_capped_surface_creator_test_suite_hpp
#define piecewise_capped_surface_creator_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"

#include "eli/geom/curve/piecewise.hpp"

#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/piecewise_capped_surface_creator.hpp"

template<typename data__>
class piecewise_capped_surface_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;
    typedef typename piecewise_surface_type::tolerance_type tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_capped_surface_creator_test_suite<float>::create_flat_capped_surface_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_capped_surface_creator_test_suite<double>::create_flat_capped_surface_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_capped_surface_creator_test_suite<long double>::create_flat_capped_surface_test);
    }

  public:
    piecewise_capped_surface_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_capped_surface_creator_test_suite()
    {
    }

  private:
    void create_flat_capped_surface_test()
    {
#if 0
      typedef typename piecewise_curve_type::curve_type curve_type;
      typedef typename curve_type::control_point_type curve_control_point_type;

      // create geometry with default parameterizations
      {
        piecewise_curve_type pc;
        curve_type c(3);
        curve_control_point_type cp[4];
        data_type k=4*(eli::constants::math<data_type>::sqrt_two()-1)/3;
        index_type i;

        // create curve
        cp[0] << 1, 0, 0;
        cp[1] << 1, k, 0;
        cp[2] << k, 1, 0;
        cp[3] << 0, 1, 0;
        for (i=0; i<4; ++i)
        {
          c.set_control_point(cp[i], i);
        }
        TEST_ASSERT(pc.push_back(c, 0.25)==piecewise_curve_type::NO_ERRORS);

        // set 2nd quadrant curve
        cp[0] <<  0, 1, 0;
        cp[1] << -k, 1, 0;
        cp[2] << -1, k, 0;
        cp[3] << -1, 0, 0;
        for (i=0; i<4; ++i)
        {
          c.set_control_point(cp[i], i);
        }
        TEST_ASSERT(pc.push_back(c, 0.25)==piecewise_curve_type::NO_ERRORS);

        piecewise_surface_type ps;

        TEST_ASSERT(eli::geom::surface::create_body_of_revolution(ps, pc, 0, true));

//         if (typeid(data_type)==typeid(double))
//           octave_print(2, ps);

        TEST_ASSERT(ps.open_u());
        TEST_ASSERT(ps.closed_v());
      }
#endif
    }
};

#endif

