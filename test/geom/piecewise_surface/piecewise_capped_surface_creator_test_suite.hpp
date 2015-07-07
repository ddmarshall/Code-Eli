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
#include "eli/geom/curve/piecewise_circle_creator.hpp"

#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/piecewise_general_skinning_surface_creator.hpp"
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
      typedef eli::geom::surface::piecewise_capped_surface_creator<data_type, 3, tolerance_type> capped_creator_type;

      piecewise_surface_type s_orig;
      bool rtn_flag;

      // create cylinder with both ends open
      {
        typedef eli::geom::surface::piecewise_general_skinning_surface_creator<data_type, 3, tolerance_type> skinning_creator_type;
        typedef typename eli::geom::surface::connection_data<data__, 3, tolerance_type> rib_data_type;
        typedef typename rib_data_type::curve_type rib_curve_type;

        index_type nsegs(2);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename skinning_creator_type::index_type> max_degree(nsegs);
        std::vector<data_type> u(nsegs+1);
        rib_curve_type rc1, rc2, rc3;
        skinning_creator_type gc;

        // create the 3 ribs
        {
          eli::geom::curve::piecewise_circle_creator<data_type, 3, tolerance_type> circle_creator;
          point_type start, origin;

          // set the parameters for first circle
          u[0]=-1;
          start  << 1, 0, 0;
          origin << 0, 0, 0;

          // create the first circle
          circle_creator.set(start, origin);
          rtn_flag = circle_creator.create(rc1);
          TEST_ASSERT(rtn_flag);

          // set the parameters for second circle
          u[1]=1;
          start  << 1, 0, 1;
          origin << 0, 0, 1;

          // create the second circle
          circle_creator.set(start, origin);
          rtn_flag = circle_creator.create(rc2);
          TEST_ASSERT(rtn_flag);

          // set the parameters for third circle
          u[2]=2;
          start  << 1, 0, 3;
          origin << 0, 0, 3;

          // create the third circle
          circle_creator.set(start, origin);
          rtn_flag = circle_creator.create(rc3);
          TEST_ASSERT(rtn_flag);
        }

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[1].set_f(rc2);
        ribs[2].set_f(rc3);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create the cylinder
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u[0]);
        gc.set_segment_du(u[1]-u[0], 0);
        gc.set_segment_du(u[2]-u[1], 1);
        rtn_flag = gc.create(s_orig);
        TEST_ASSERT(rtn_flag);

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::test::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::test::octave_print(1, ribs[2].get_f(), "rib2", true);
//          eli::test::octave_print(1, s_orig, "surf", true);
//          eli::test::octave_finish(1, false);
//        }
      }

      // cap umin edge
      {
        piecewise_surface_type s_umin_cap(s_orig);
        capped_creator_type cc;

        rtn_flag = cc.set_conditions(s_umin_cap, 0.5, capped_creator_type::CAP_UMIN);
        TEST_ASSERT(rtn_flag);

        rtn_flag = cc.create(s_umin_cap);
        TEST_ASSERT(rtn_flag);

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, s_umin_cap, "surf", true);
//          eli::test::octave_finish(1, false);
//        }
      }

      // cap umax edge
      {
        piecewise_surface_type s_umax_cap(s_orig);
        capped_creator_type cc;

        rtn_flag = cc.set_conditions(s_umax_cap, 0.5, capped_creator_type::CAP_UMAX);
        TEST_ASSERT(rtn_flag);

        rtn_flag = cc.create(s_umax_cap);
        TEST_ASSERT(rtn_flag);

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, s_umax_cap, "surf", true);
//          eli::test::octave_finish(1, false);
//        }
      }

      // cap umin & umax edge
      {
        piecewise_surface_type s_uminmax_cap(s_orig);
        capped_creator_type cc;

        rtn_flag = cc.set_conditions(s_uminmax_cap, 0.5, capped_creator_type::CAP_UMIN);
        TEST_ASSERT(rtn_flag);

        rtn_flag = cc.create(s_uminmax_cap);
        TEST_ASSERT(rtn_flag);

        rtn_flag = cc.set_conditions(s_uminmax_cap, 0.5, capped_creator_type::CAP_UMAX);
        TEST_ASSERT(rtn_flag);

        rtn_flag = cc.create(s_uminmax_cap);
        TEST_ASSERT(rtn_flag);

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, s_uminmax_cap, "surf", true);
//          eli::test::octave_finish(1, false);
//        }
      }
    }
};

#endif

