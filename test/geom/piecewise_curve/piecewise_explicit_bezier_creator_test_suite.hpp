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
    typedef eli::geom::curve::pseudo::explicit_bezier<data__> explicit_bezier_type;
    typedef typename explicit_bezier_type::control_point_type control_point_type;

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
      typedef eli::geom::curve::bezier<data__, 3> bezier_curve_type;

      eli::geom::curve::piecewise_explicit_bezier_creator<data_type, 3, tolerance_type> pebc;
      piecewise_curve_type pc1, pc_ref;
      control_point_type cntrl_in[5];
      typename bezier_curve_type::control_point_type bez_cntrl[5];
      explicit_bezier_type ebc;
      bezier_curve_type bc, bc_out;
      typename explicit_bezier_type::data_type t0, t1;
      bool rtn_flag;

      // set control points and create curves
      cntrl_in[0] << 2.0;
      cntrl_in[1] << 1.5;
      cntrl_in[2] << 0.0;
      cntrl_in[3] << 1.0;
      cntrl_in[4] << 0.5;
      bez_cntrl[0] << 0,    2,   0;
      bez_cntrl[1] << 0.25, 1.5, 0;
      bez_cntrl[2] << 0.5,  0,   0;
      bez_cntrl[3] << 0.75, 1,   0;
      bez_cntrl[4] << 1,    0.5, 0;

      ebc.resize(4);
      bc.resize(4);
      for (index_type i=0; i<5; ++i)
      {
        ebc.set_control_point(cntrl_in[i], i);
        bc.set_control_point(bez_cntrl[i], i);
      }

      // set the parameterization
      t0=-1;
      t1=1;

      // build reference piecewise curve
      pc_ref.clear();
      pc_ref.set_t0(t0);
      pc_ref.push_back(bc, t1-t0);

      // create curve
      rtn_flag=pebc.set_conditions(ebc, false);
      TEST_ASSERT(rtn_flag);
      pebc.set_t0(t0);
      pebc.set_segment_dt(t1-t0, 0);
      rtn_flag=pebc.create(pc1);
      TEST_ASSERT(rtn_flag);

      // extract control points and test
      TEST_ASSERT(pc1.number_segments()==1);
      TEST_ASSERT(pc1.get_t0()==t0);
      TEST_ASSERT(pc1.get_tmax()==t1);
      pc1.get(bc_out, 0);
      for (index_type i=0; i<5; ++i)
      {
        TEST_ASSERT(tol.approximately_equal(bc_out.get_control_point(i), bez_cntrl[i]));
      }

      if (typeid(data_type)==typeid(float))
      {
        std::cout.flush();
        eli::test::octave_start(1);
//        eli::test::octave_print(1, ebc, "explicit_bezier");
//        eli::test::octave_print(1, bc, "bezier");
//        eli::test::octave_print(1, pc_ref, "ref_piecewise");
        eli::test::octave_print(1, pc1, "exp_piecewise");
        eli::test::octave_finish(1);
      }
    }
};

#endif
