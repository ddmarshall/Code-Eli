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

#ifndef piecewise_four_digit_creator_test_suite_hpp
#define piecewise_four_digit_creator_test_suite_hpp

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
#include "eli/geom/curve/piecewise_four_digit_creator.hpp"
#include "eli/geom/curve/pseudo/four_digit.hpp"

#include "octave_helpers.hpp"

template<typename data__>
class piecewise_four_digit_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::piecewise_four_digit_creator<data__, 3, tolerance_type> four_digit_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_four_digit_creator_test_suite<float>::create_symmetric_airfoil_test);
      TEST_ADD(piecewise_four_digit_creator_test_suite<float>::create_general_airfoil_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_four_digit_creator_test_suite<double>::create_symmetric_airfoil_test);
      TEST_ADD(piecewise_four_digit_creator_test_suite<double>::create_general_airfoil_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_four_digit_creator_test_suite<long double>::create_symmetric_airfoil_test);
      TEST_ADD(piecewise_four_digit_creator_test_suite<long double>::create_general_airfoil_test);
    }

  public:
    piecewise_four_digit_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_four_digit_creator_test_suite()
    {
    }

  private:
    void create_symmetric_airfoil_test()
    {
      four_digit_type af;
      eli::geom::curve::pseudo::four_digit<data_type> four_series;
      typename eli::geom::curve::pseudo::four_digit<data_type>::point_type pt_ref;

      data_type th, cam, cam_loc, t, eps(std::numeric_limits<data_type>::epsilon());
      point_type pt;
      bool rtn;

      // set airfoil thickness
      th=24;
      rtn=af.set_thickness(th);
      TEST_ASSERT(rtn);
      TEST_ASSERT(af.get_thickness()==th);

      // set airfoil camber
      cam=0;
      cam_loc=0;
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

      piecewise_curve_type af_pc;
      bool fit_success;

      fit_success = af.create(af_pc);
      TEST_ASSERT(fit_success);

      four_series.set_sharp_trailing_edge(true);
      four_series.set_thickness(th);
      four_series.set_camber(cam, cam_loc);

//      if (typeid(data_type)==typeid(float))
//      {
//        std::cout.flush();
//        eli::test::octave_start(1);
//        eli::test::octave_print(1, four_series, "poly");
//        eli::test::octave_print(1, af_pc, "piecewise");
//        eli::test::octave_finish(1);
//      }

      // test lower surface trailing edge
      t=0;
      pt=af_pc.f(t);
      pt_ref=four_series.f(t-1);
      TEST_ASSERT((pt.x()==pt_ref.x()) && (pt.y()==pt_ref.y()));

      // test leading edge
      t=1;
      pt=af_pc.f(t);
      pt_ref=four_series.f(t-1);
      TEST_ASSERT((pt.x()==pt_ref.x()) && (pt.y()==pt_ref.y()));

      // test upper surface trailing edge
      t=2;
      pt=af_pc.f(t);
      pt_ref=four_series.f(t-1);
      TEST_ASSERT((pt.x()==pt_ref.x()) && (pt.y()==pt_ref.y()));

      // test point on aft lower surface
      t=1-0.16;
      pt=af_pc.f(t);
      pt_ref=four_series.f(-(1-t)*(1-t));
      if (typeid(data_type)==typeid(long double))
      {
        TEST_ASSERT(std::sqrt((pt.x()-pt_ref.x())*(pt.x()-pt_ref.x())+(pt.y()-pt_ref.y())*(pt.y()-pt_ref.y()))<30*eps);
      }
      else
      {
        TEST_ASSERT(std::sqrt((pt.x()-pt_ref.x())*(pt.x()-pt_ref.x())+(pt.y()-pt_ref.y())*(pt.y()-pt_ref.y()))<eps);
      }

      // test point on fore lower surface
      t=1-0.64;
      pt=af_pc.f(t);
      pt_ref=four_series.f(-(1-t)*(1-t));
      if (typeid(data_type)==typeid(long double))
      {
        TEST_ASSERT(std::sqrt((pt.x()-pt_ref.x())*(pt.x()-pt_ref.x())+(pt.y()-pt_ref.y())*(pt.y()-pt_ref.y()))<60*eps);
      }
      else
      {
        TEST_ASSERT(std::sqrt((pt.x()-pt_ref.x())*(pt.x()-pt_ref.x())+(pt.y()-pt_ref.y())*(pt.y()-pt_ref.y()))<eps);
      }

      // test point on fore upper surface
      t=1+0.16;
      pt=af_pc.f(t);
      pt_ref=four_series.f((t-1)*(t-1));
      if (typeid(data_type)==typeid(long double))
      {
        TEST_ASSERT(std::sqrt((pt.x()-pt_ref.x())*(pt.x()-pt_ref.x())+(pt.y()-pt_ref.y())*(pt.y()-pt_ref.y()))<30*eps);
      }
      else
      {
        TEST_ASSERT(std::sqrt((pt.x()-pt_ref.x())*(pt.x()-pt_ref.x())+(pt.y()-pt_ref.y())*(pt.y()-pt_ref.y()))<eps);
      }

      // test point on aft upper surface
      t=1+0.64;
      pt=af_pc.f(t);
      pt_ref=four_series.f((t-1)*(t-1));
      if (typeid(data_type)==typeid(long double))
      {
        TEST_ASSERT(std::sqrt((pt.x()-pt_ref.x())*(pt.x()-pt_ref.x())+(pt.y()-pt_ref.y())*(pt.y()-pt_ref.y()))<60*eps);
      }
      else
      {
        TEST_ASSERT(std::sqrt((pt.x()-pt_ref.x())*(pt.x()-pt_ref.x())+(pt.y()-pt_ref.y())*(pt.y()-pt_ref.y()))<eps);
      }
    }

    void create_general_airfoil_test()
    {
      four_digit_type af;

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

    }
};

#endif

