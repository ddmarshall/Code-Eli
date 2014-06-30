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

#ifndef piecewise_general_skinning_surface_creator_test_suite_suite_hpp
#define piecewise_general_skinning_surface_creator_test_suite_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_linear_creator.hpp"

#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/piecewise_general_skinning_surface_creator.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits
#include <iterator> // std::insert_iterator

template<typename data__>
class piecewise_general_skinning_surface_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;
    typedef typename piecewise_surface_type::tolerance_type tolerance_type;
    typedef typename eli::geom::surface::connection_data<data__, 3, tolerance_type> rib_data_type;
    typedef typename rib_data_type::curve_type rib_curve_type;
    typedef typename eli::geom::curve::piecewise_linear_creator<data__, 3, tolerance_type> piecewise_line_creator_type;
    typedef typename eli::geom::surface::general_skinning_surface_creator<data__, 3, tolerance_type> general_creator_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<float>::create_rib_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<float>::create_single_surface_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<float>::create_multijoint_rib_single_surface_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<float>::create_multirib_surface_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<double>::create_rib_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<double>::create_single_surface_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<double>::create_multijoint_rib_single_surface_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<double>::create_multirib_surface_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<long double>::create_rib_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<long double>::create_single_surface_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<long double>::create_multijoint_rib_single_surface_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<long double>::create_multirib_surface_test);
    }

  public:
    piecewise_general_skinning_surface_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_general_skinning_surface_creator_test_suite()
    {
    }

  private:
    void create_rib_test()
    {
      rib_curve_type rc1, rc2, rc3;
      point_type p1, p2, p3;
      bool rtn_flag;

      p1 << 1, 0, 0;
      p2 << 2, 2, 2;
      p3 << 3, 2, 1;

      // create three rib curves
      piecewise_line_creator_type plc(1);

      plc.set_corner(p1, 0);
      plc.set_corner(p2, 1);
      plc.create(rc1);

      plc.set_corner(p1, 0);
      plc.set_corner(p3, 1);
      plc.create(rc2);

      plc.set_corner(p3, 0);
      plc.set_corner(p2, 1);
      plc.create(rc3);

      // set rib
      {
        rib_data_type rib;

        // should start in bad state
        rtn_flag=rib.check_state();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(!rib.use_f());

        // add rib
        rtn_flag=rib.set_f(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());

        // unset rib
        rtn_flag=rib.unset_f();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(!rib.use_f());
      }

      // set rib and fp
      {
        rib_data_type rib;

        // should start in bad state
        rtn_flag=rib.check_state();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(!rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());

        // add left fp
        rtn_flag=rib.set_left_fp(rc1);
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_left_fp()==rc1);
        TEST_ASSERT(!rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());

        // add rib
        rtn_flag=rib.set_f(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());

        // add right fp
        rtn_flag=rib.set_right_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_right_fp()==rc2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());

        // remove left fp
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());

        // remove left fp again
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());

        // remove right fp
        rtn_flag=rib.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());

        // add fp
        rtn_flag=rib.set_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_left_fp()==rc2);
        TEST_ASSERT(rib.get_right_fp()==rc2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());

        // remove right fp
        rtn_flag=rib.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
      }

      // set rib, fp and fpp
      {
        rib_data_type rib;

        // should start in bad state
        rtn_flag=rib.check_state();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(!rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add left fpp
        rtn_flag=rib.set_left_fpp(rc1);
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_left_fpp()==rc1);
        TEST_ASSERT(!rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add right fpp
        rtn_flag=rib.set_right_fpp(rc2);
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_right_fpp()==rc2);
        TEST_ASSERT(!rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // add rib
        rtn_flag=rib.set_f(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // remove right fpp
        rtn_flag=rib.unset_right_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fp());

        // remove right fpp again
        rtn_flag=rib.unset_right_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fp());

        // remove left fpp
        rtn_flag=rib.unset_left_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add fp
        rtn_flag=rib.set_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_left_fp()==rc2);
        TEST_ASSERT(rib.get_right_fp()==rc2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add fpp
        rtn_flag=rib.set_fpp(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_left_fpp()==rc1);
        TEST_ASSERT(rib.get_right_fpp()==rc1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // remove right fp
        rtn_flag=rib.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // remove right fpp
        rtn_flag=rib.unset_right_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // remove everything
        rtn_flag=rib.unset_fpp();
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.unset_fp();
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.unset_f();
        TEST_ASSERT(!rtn_flag);

        // add everything to right only
        rtn_flag=rib.set_f(rc2);
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.set_right_fp(rc1);
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.set_right_fpp(rc3);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.get_f()==rc2);
        TEST_ASSERT(rib.get_right_fp()==rc1);
        TEST_ASSERT(rib.get_right_fpp()==rc3);
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());
      }

      // set conditions with C1
      {
        rib_data_type rib;

        // start with point
        rtn_flag=rib.set_f(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // set continuity to C1
        rtn_flag=rib.set_continuity(rib_data_type::C1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add left fp
        rtn_flag=rib.set_left_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // unset left fp
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        rtn_flag=rib.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add right fp
        rtn_flag=rib.set_right_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // change continuity
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        rtn_flag=rib.set_continuity(rib_data_type::C0);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add fp
        rtn_flag=rib.set_fp(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // change continuity
        rtn_flag=rib.set_continuity(rib_data_type::C1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // unset left fp
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // set left fp to change left and right value
        rtn_flag=rib.set_left_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_left_fp()==rc2);
        TEST_ASSERT(rib.get_right_fp()==rc2);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // set left fp
        rtn_flag=rib.set_left_fp(rib.get_right_fp());
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());
      }

      // set point, fp, fpp and C2
      {
        rib_data_type rib;

        // start with point and fp
        rtn_flag=rib.set_f(rc1);
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.set_fp(rc2);
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.set_continuity(rib_data_type::C1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // set continuity to C2
        rtn_flag=rib.set_continuity(rib_data_type::C2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // unset left fp
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // unset right fp
        rtn_flag=rib.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // set fpp
        rtn_flag=rib.set_fpp(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // unset left fpp
        rtn_flag=rib.unset_left_fpp();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // set left fpp to set both left and right values
        rtn_flag=rib.set_left_fpp(rc3);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // set left fpp
        rtn_flag=rib.set_left_fpp(rib.get_right_fpp());
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // set fp
        rtn_flag=rib.set_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());
      }

      // copy constructor, equivalence and assignment operators
      {
        rib_data_type r1, r2, r3;

        // defaults should be the same
        TEST_ASSERT(r1==r2);

        // manually build identical ones
        r1.set_f(rc2);
        r1.set_left_fp(rc1);
        r1.set_right_fpp(rc3);
        r2.set_f(rc2);
        r2.set_left_fp(rc1);
        r2.set_right_fpp(rc3);
        TEST_ASSERT(r1==r2);

        // assignment operator
        r3=r2;
        TEST_ASSERT(r3==r1);
        TEST_ASSERT(r3==r2);

        // set value for r2's fpp and then don't use it
        r1.unset_fpp();
        r2.unset_fpp();
        r3.unset_fpp();
        r2.set_fpp(rc1);
        r2.unset_fpp();
        TEST_ASSERT(r3==r1);
        TEST_ASSERT(r3==r2);

        // copy constructor
        rib_data_type r4(r3);
        TEST_ASSERT(r4==r3);
      }

      // test the joint interface
      {
        rib_data_type r1;
        data_type s1, s2, s3, s4, s5;

        s1=static_cast<data_type>(0.25);
        s2=static_cast<data_type>(0.35);
        s3=static_cast<data_type>(0.45);
        s4=static_cast<data_type>(0.55);
        s5=static_cast<data_type>(0.70);
        rc1.split(s1);
        rc1.split(s4);
        rc2.split(s2);
        rc2.split(s3);
        rc2.split(s4);
        rc3.split(s3);
        rc3.split(s5);

        r1.set_f(rc1);
        r1.set_fp(rc2);
        r1.set_fpp(rc3);

        // get the joints from the rib
        std::vector<data_type> rjoints;
        r1.get_joints(std::back_inserter(rjoints));

        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }

        // split the ribs at the specified locations
        std::vector<index_type> r1degs;
        r1.split(rjoints.begin(), rjoints.end(), std::back_inserter(r1degs));

        // check the degree
        TEST_ASSERT(r1degs.size()==6);
        if (r1degs.size()==6)
        {
          TEST_ASSERT(r1degs[0]==1);
          TEST_ASSERT(r1degs[1]==1);
          TEST_ASSERT(r1degs[2]==1);
          TEST_ASSERT(r1degs[3]==1);
          TEST_ASSERT(r1degs[4]==1);
          TEST_ASSERT(r1degs[5]==1);
        }

        // check the joint locations
        rjoints.clear();
        r1.get_f().get_parameters(std::back_inserter(rjoints));
        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }
        rjoints.clear();
        r1.get_left_fp().get_parameters(std::back_inserter(rjoints));
        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }
        rjoints.clear();
        r1.get_right_fp().get_parameters(std::back_inserter(rjoints));
        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }
        rjoints.clear();
        r1.get_left_fpp().get_parameters(std::back_inserter(rjoints));
        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }
        rjoints.clear();
        r1.get_right_fpp().get_parameters(std::back_inserter(rjoints));
        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }
      }
    }

    void create_single_surface_test()
    {
      // simple surface connecting 2 lines
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;

        // create two rib curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(p00, 0);
        plc.set_corner(p01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc1);

        plc.set_corner(p10, 0);
        plc.set_corner(p11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[1].set_f(rc2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=(v0+v1)/2;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_finish(1);
//        }
      }

      // surface connecting 2 lines with specified slopes
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rs1, rs2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, s00, s01, s10, s11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;

        // create two rib curves && two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(p00, 0);
        plc.set_corner(p01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc1);

        plc.set_corner(s00, 0);
        plc.set_corner(s01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(p10, 0);
        plc.set_corner(p11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc2);

        plc.set_corner(s10, 0);
        plc.set_corner(s11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[1].set_f(rc2);
        ribs[1].set_left_fp(rs2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=(v0+v1)/2;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", false);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_finish(1);
//        }
      }

      // surface connecting 2 lines with specified 1st and 2nd derivatives
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rs1, rs2, ra1, ra2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, s00, s01, s10, s11, a00, a01, a10, a11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;
        a00 << 0, 1, 0;
        a01 << 1, 0, 1;
        a10 << 1, 0, 1;
        a11 << 0, 1, 0;

        // create two rib curves && two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(p00, 0);
        plc.set_corner(p01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc1);

        plc.set_corner(s00, 0);
        plc.set_corner(s01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(a00, 0);
        plc.set_corner(a01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra1);

        plc.set_corner(p10, 0);
        plc.set_corner(p11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc2);

        plc.set_corner(s10, 0);
        plc.set_corner(s11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        plc.set_corner(a10, 0);
        plc.set_corner(a11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[0].set_right_fpp(ra1);
        ribs[1].set_f(rc2);
        ribs[1].set_left_fp(rs2);
        ribs[1].set_left_fpp(ra2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=(v0+v1)/2;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", false);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_finish(1);
//        }
      }

      // surface connecting 2 cubic curves with specified 1st and 2nd derivatives
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rs1, rs2, ra1, ra2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, s00, s01, s10, s11, a00, a01, a10, a11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;
        a00 << 0, 1, 0;
        a01 << 1, 0, 1;
        a10 << 1, 0, 1;
        a11 << 0, 1, 0;

        // create two rib curves && two slope curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);

        // create 1st and 2nd derivative curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(s00, 0);
        plc.set_corner(s01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(a00, 0);
        plc.set_corner(a01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra1);

        plc.set_corner(s10, 0);
        plc.set_corner(s11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        plc.set_corner(a10, 0);
        plc.set_corner(a11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[0].set_right_fpp(ra1);
        ribs[1].set_f(rc2);
        ribs[1].set_left_fp(rs2);
        ribs[1].set_left_fpp(ra2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=(v0+v1)/2;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", false);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_finish(1);
//        }
      }

      // surface connecting 1 line and 1 cubic curves with specified 1st and 2nd derivatives
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rs1, rs2, ra1, ra2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, s00, s01, s10, s11, a00, a01, a10, a11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;
        a00 << 0, 1, 0;
        a01 << 1, 0, 1;
        a10 << 1, 0, 1;
        a11 << 0, 1, 0;

        // create two rib curves && two slope curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curve for rib
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);

        // create 1st and 2nd derivative curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(s00, 0);
        plc.set_corner(s01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(a00, 0);
        plc.set_corner(a01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra1);

        plc.set_corner(p10, 0);
        plc.set_corner(p11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc2);

        plc.set_corner(s10, 0);
        plc.set_corner(s11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        plc.set_corner(a10, 0);
        plc.set_corner(a11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[0].set_right_fpp(ra1);
        ribs[1].set_f(rc2);
        ribs[1].set_left_fp(rs2);
        ribs[1].set_left_fpp(ra2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=(v0+v1)/2;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_finish(1);
//        }
      }
    }

    void create_multijoint_rib_single_surface_test()
    {
      // surface connecting ribs made of 2 line segments each
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;

        // create two rib curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(p00, 0);
        plc.set_corner(p01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc1);
        rc1.split((v0+v1)/2);

        plc.set_corner(p10, 0);
        plc.set_corner(p11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc2);
        rc2.split((v0+v1)/2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[1].set_f(rc2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_finish(1);
//        }
      }

      // surface connecting ribs made of 2 cubic segments each
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;

        // create two rib curves && two slope curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);
        rc1.split((v0+v1)/2);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);
        rc2.split((v0+v1)/2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[1].set_f(rc2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_finish(1);
//        }
      }

      // surface connecting ribs made of 2 cubic segments each, specified 1st
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rs1, rs2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, s00, s01, s10, s11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;

        // create two rib curves && two slope curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);
        rc1.split((v0+v1)/2);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);
        rc2.split((v0+v1)/2);

        // create two two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(s00/2, 0);
        plc.set_corner(s01/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(s10/2, 0);
        plc.set_corner(s11/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[1].set_f(rc2);
        ribs[1].set_left_fp(rs2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_finish(1);
//        }
      }

      // surface connecting ribs made of 2 cubic segments each, specified 1st and 2nd derivatives
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rs1, rs2, ra1, ra2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, s00, s01, s10, s11, a00, a01, a10, a11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;
        a00 << 0, 1, 0;
        a01 << 1, 0, 1;
        a10 << 1, 0, 1;
        a11 << 0, 1, 0;

        // create two rib curves && two slope curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);
        rc1.split((v0+v1)/2);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);
        rc2.split((v0+v1)/2);

        // create two two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(s00/2, 0);
        plc.set_corner(s01/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(a00, 0);
        plc.set_corner(a01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra1);

        plc.set_corner(s10/2, 0);
        plc.set_corner(s11/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        plc.set_corner(a10, 0);
        plc.set_corner(a11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[0].set_right_fpp(ra1);
        ribs[1].set_f(rc2);
        ribs[1].set_left_fp(rs2);
        ribs[1].set_left_fpp(ra2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_finish(1);
//        }
      }

      // surface connecting 1 rib with 1 cubic segment and 1 rib made of 2 cubic segments, specified 1st and 2nd derivatives
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rs1, rs2, ra1, ra2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, s00, s01, s10, s11, a00, a01, a10, a11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;
        a00 << 0, 1, 0;
        a01 << 1, 0, 1;
        a10 << 1, 0, 1;
        a11 << 0, 1, 0;

        // create two rib curves && two slope curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);
        rc1.split((v0+v1)/2);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);

        // create two two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(s00/2, 0);
        plc.set_corner(s01/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(a00, 0);
        plc.set_corner(a01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra1);

        plc.set_corner(s10/2, 0);
        plc.set_corner(s11/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        plc.set_corner(a10, 0);
        plc.set_corner(a11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[0].set_right_fpp(ra1);
        ribs[1].set_f(rc2);
        ribs[1].set_left_fp(rs2);
        ribs[1].set_left_fpp(ra2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_finish(1);
//        }
      }

      // surface connecting ribs made of 2 cubic segments each with different split locations, specified 1st and 2nd derivatives
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rs1, rs2, ra1, ra2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, s00, s01, s10, s11, a00, a01, a10, a11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;
        a00 << 0, 1, 0;
        a01 << 1, 0, 1;
        a10 << 1, 0, 1;
        a11 << 0, 1, 0;

        // create two rib curves && two slope curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);
        rc1.split(v0+(v1-v0)/4);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);
        rc2.split(v0+3*(v1-v0)/4);

        // create two two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(s00/2, 0);
        plc.set_corner(s01/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(a00, 0);
        plc.set_corner(a01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra1);

        plc.set_corner(s10/2, 0);
        plc.set_corner(s11/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        plc.set_corner(a10, 0);
        plc.set_corner(a11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[0].set_right_fpp(ra1);
        ribs[1].set_f(rc2);
        ribs[1].set_left_fp(rs2);
        ribs[1].set_left_fpp(ra2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=v0+(v1-v0)/8;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_finish(1);
//        }
      }
    }

    void create_multirib_surface_test()
    {
      // surface connecting 3 ribs made of 1 line segment each
      {
        index_type nsegs(2);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rc3;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, p20, p21;
        data_type u0(2), v0(1), u1(4), v1(5), u2(7), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        p20 << 5, 2, 5;
        p21 << 5, 5, 6;

        // create two rib curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(p00, 0);
        plc.set_corner(p01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc1);

        plc.set_corner(p10, 0);
        plc.set_corner(p11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc2);

        plc.set_corner(p20, 0);
        plc.set_corner(p21, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc3);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[1].set_f(rc2);
        ribs[2].set_f(rc3);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        gc.set_segment_du(u2-u1, 1);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u2, v0);
        TEST_ASSERT(tol.approximately_equal(p20, p_test));
//        std::cout << "p=" << p20
//                  << "\tpc=" << ribs[2].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p20-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        p_test=s.f(u2, v1);
        TEST_ASSERT(tol.approximately_equal(p21, p_test));
//        std::cout << "p=" << p21
//                  << "\tpc=" << ribs[2].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p21-p_test).norm() << std::endl;
        v=v0+(v1-v0)/2;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[2].get_f(), "rib2", true);
//          eli::octave_finish(1);
//        }
      }

      // surface connecting 3 ribs made of 2 cubic segments each
      {
        index_type nsegs(2);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rc3;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, p20, p21;
        data_type u0(2), v0(1), u1(4), v1(5), u2(7), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        p20 << 5, 2, 5;
        p21 << 5, 5, 6;

        // create rib curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);
        rc1.split((v0+v1)/2);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);
        rc2.split((v0+v1)/2);

        cp[0]=p20;
        cp[1] << 4, 3, 5;
        cp[2] << 4, 4, 6;
        cp[3]=p21;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc3.push_front(curve, v1-v0);
        rc3.set_t0(v0);
        rc3.split((v0+v1)/2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[1].set_f(rc2);
        ribs[2].set_f(rc3);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        gc.set_segment_du(u2-u1, 1);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u2, v0);
        TEST_ASSERT(tol.approximately_equal(p20, p_test));
//        std::cout << "p=" << p20
//                  << "\tpc=" << ribs[2].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p20-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        p_test=s.f(u2, v1);
        TEST_ASSERT(tol.approximately_equal(p21, p_test));
//        std::cout << "p=" << p21
//                  << "\tpc=" << ribs[2].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p21-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[2].get_f(), "rib2", true);
//          eli::octave_finish(1);
//        }
      }

      // surface connecting ribs made of 3 cubic segments each, specified continuous 1st derivative
      {
        index_type nsegs(2);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rc3, rs1, rs2, rs3;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, p20, p21, s00, s01, s10, s11, s20, s21;
        data_type u0(2), v0(1), u1(4), v1(5), u2(7), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        p20 << 5, 2, 5;
        p21 << 5, 5, 6;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;
        s20 << 1, 0, 1;
        s21 << 1, 0, 1;

        // create rib curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);
        rc1.split((v0+v1)/2);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);
        rc2.split((v0+v1)/2);

        cp[0]=p20;
        cp[1] << 4, 3, 5;
        cp[2] << 4, 4, 6;
        cp[3]=p21;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc3.push_front(curve, v1-v0);
        rc3.set_t0(v0);
        rc3.split((v0+v1)/2);

        // create two two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(s00/2, 0);
        plc.set_corner(s01/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(s10/2, 0);
        plc.set_corner(s11/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        plc.set_corner(s20/2, 0);
        plc.set_corner(s21/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs3);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[1].set_f(rc2);
        ribs[1].set_fp(rs2);
        ribs[2].set_f(rc3);
        ribs[2].set_left_fp(rs3);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        gc.set_segment_du(u2-u1, 1);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u2, v0);
        TEST_ASSERT(tol.approximately_equal(p20, p_test));
//        std::cout << "p=" << p20
//                  << "\tpc=" << ribs[2].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p20-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        p_test=s.f(u2, v1);
        TEST_ASSERT(tol.approximately_equal(p21, p_test));
//        std::cout << "p=" << p21
//                  << "\tpc=" << ribs[2].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p21-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u2, v);
        p_ref=ribs[2].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u2, v);
        p_ref=ribs[2].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[2].get_f(), "rib2", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_right_fp(), "rve1");
//          eli::octave_print(1, ribs[2].get_f(), ribs[2].get_left_fp(), "lve2");
//          eli::octave_finish(1);
//        }
      }

      // surface connecting ribs made of 3 cubic segments each, specified discontinuous 1st derivative
      {
        index_type nsegs(2);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rc3, rs1, rs2l, rs2r, rs3;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, p20, p21, s00, s01, sl10, sl11, sr10, sr11, s20, s21;
        data_type u0(2), v0(1), u1(4), v1(5), u2(7), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        p20 << 5, 2, 5;
        p21 << 5, 5, 6;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        sl10 << 0, 1, 1;
        sl11 << 0, 0, 2;
        sr10 << 2, 1, 0;
        sr11 << 2, 0, 1;
        s20 << 1, 0, 1;
        s21 << 1, 0, 1;

        // create rib curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);
        rc1.split((v0+v1)/2);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);
        rc2.split((v0+v1)/2);

        cp[0]=p20;
        cp[1] << 4, 3, 5;
        cp[2] << 4, 4, 6;
        cp[3]=p21;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc3.push_front(curve, v1-v0);
        rc3.set_t0(v0);
        rc3.split((v0+v1)/2);

        // create two two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(s00/2, 0);
        plc.set_corner(s01/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(sl10/2, 0);
        plc.set_corner(sl11/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2l);

        plc.set_corner(sr10/2, 0);
        plc.set_corner(sr11/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2r);

        plc.set_corner(s20/2, 0);
        plc.set_corner(s21/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs3);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[1].set_f(rc2);
        ribs[1].set_left_fp(rs2l);
        ribs[1].set_right_fp(rs2r);
        ribs[2].set_f(rc3);
        ribs[2].set_left_fp(rs3);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        gc.set_segment_du(u2-u1, 1);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u2, v0);
        TEST_ASSERT(tol.approximately_equal(p20, p_test));
//        std::cout << "p=" << p20
//                  << "\tpc=" << ribs[2].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p20-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        p_test=s.f(u2, v1);
        TEST_ASSERT(tol.approximately_equal(p21, p_test));
//        std::cout << "p=" << p21
//                  << "\tpc=" << ribs[2].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p21-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1-small_num, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1+small_num, v);
        p_ref=ribs[1].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u2, v);
        p_ref=ribs[2].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1-small_num, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1+small_num, v);
        p_ref=ribs[1].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u2, v);
        p_ref=ribs[2].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[2].get_f(), "rib2", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_right_fp(), "rve1");
//          eli::octave_print(1, ribs[2].get_f(), ribs[2].get_left_fp(), "lve2");
//          eli::octave_finish(1);
//        }
      }

      // surface connecting ribs made of 3 cubic segments each, specified continuous 1st & 2nd derivatives
      {
        index_type nsegs(2);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rc3, rs1, rs2, rs3, ra1, ra2, ra3;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, p20, p21, s00, s01, s10, s11, s20, s21, a00, a01, a10, a11, a20, a21;
        data_type u0(2), v0(1), u1(4), v1(5), u2(7), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        p20 << 5, 2, 5;
        p21 << 5, 5, 6;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;
        s20 << 1, 0, 1;
        s21 << 1, 0, 1;
        a00 << 0, 1, 0;
        a01 << 1, 0, 1;
        a10 << 1, 0, 1;
        a11 << 0, 1, 0;
        a20 << 0, 1, 0;
        a21 << 1, 0, 1;

        // create rib curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);
        rc1.split((v0+v1)/2);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);
        rc2.split((v0+v1)/2);

        cp[0]=p20;
        cp[1] << 4, 3, 5;
        cp[2] << 4, 4, 6;
        cp[3]=p21;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc3.push_front(curve, v1-v0);
        rc3.set_t0(v0);
        rc3.split((v0+v1)/2);

        // create two two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(s00/2, 0);
        plc.set_corner(s01/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(a00, 0);
        plc.set_corner(a01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra1);

        plc.set_corner(s10/2, 0);
        plc.set_corner(s11/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        plc.set_corner(a10, 0);
        plc.set_corner(a11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra2);

        plc.set_corner(s20/2, 0);
        plc.set_corner(s21/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs3);

        plc.set_corner(a20, 0);
        plc.set_corner(a21, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra3);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[0].set_right_fpp(ra1);
        ribs[1].set_f(rc2);
        ribs[1].set_fp(rs2);
        ribs[1].set_fpp(ra2);
        ribs[2].set_f(rc3);
        ribs[2].set_left_fp(rs3);
        ribs[2].set_left_fpp(ra3);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        gc.set_segment_du(u2-u1, 1);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u2, v0);
        TEST_ASSERT(tol.approximately_equal(p20, p_test));
//        std::cout << "p=" << p20
//                  << "\tpc=" << ribs[2].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p20-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        p_test=s.f(u2, v1);
        TEST_ASSERT(tol.approximately_equal(p21, p_test));
//        std::cout << "p=" << p21
//                  << "\tpc=" << ribs[2].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p21-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u2, v);
        p_ref=ribs[2].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u2, v);
        p_ref=ribs[2].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u2, v);
        p_ref=ribs[2].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u2, v);
        p_ref=ribs[2].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[2].get_f(), "rib2", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_right_fp(), "rve1");
//          eli::octave_print(1, ribs[2].get_f(), ribs[2].get_left_fp(), "lve2");
//          eli::octave_finish(1);
//        }
      }

      // surface connecting ribs made of 3 cubic segments each, specified continuous 1st derivative & discontinuous 2nd derivative
      {
        index_type nsegs(2);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rc3, rs1, rs2, rs3, ra1, ra2l, ra2r, ra3;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, p20, p21, s00, s01, s10, s11, s20, s21, a00, a01, al10, al11, ar10, ar11, a20, a21;
        data_type u0(2), v0(1), u1(4), v1(5), u2(7), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        p20 << 5, 2, 5;
        p21 << 5, 5, 6;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;
        s20 << 1, 0, 1;
        s21 << 1, 0, 1;
        a00 << 0, 1, 0;
        a01 << 1, 0, 1;
        al10 << 1, 0, 1;
        al11 << 0, 1, 0;
        ar10 << -1, 0, -1;
        ar11 << 0, -1, 0;
        a20 << 0, 1, 0;
        a21 << 1, 0, 1;

        // create rib curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);
        rc1.split((v0+v1)/2);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);
        rc2.split((v0+v1)/2);

        cp[0]=p20;
        cp[1] << 4, 3, 5;
        cp[2] << 4, 4, 6;
        cp[3]=p21;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc3.push_front(curve, v1-v0);
        rc3.set_t0(v0);
        rc3.split((v0+v1)/2);

        // create two two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(s00/2, 0);
        plc.set_corner(s01/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_corner(a00, 0);
        plc.set_corner(a01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra1);

        plc.set_corner(s10/2, 0);
        plc.set_corner(s11/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        plc.set_corner(al10, 0);
        plc.set_corner(al11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra2l);

        plc.set_corner(ar10, 0);
        plc.set_corner(ar11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra2r);

        plc.set_corner(s20/2, 0);
        plc.set_corner(s21/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs3);

        plc.set_corner(a20, 0);
        plc.set_corner(a21, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra3);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[0].set_right_fpp(ra1);
        ribs[1].set_f(rc2);
        ribs[1].set_fp(rs2);
        ribs[1].set_left_fpp(ra2l);
        ribs[1].set_right_fpp(ra2r);
        ribs[2].set_f(rc3);
        ribs[2].set_left_fp(rs3);
        ribs[2].set_left_fpp(ra3);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        gc.set_segment_du(u2-u1, 1);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u2, v0);
        TEST_ASSERT(tol.approximately_equal(p20, p_test));
//        std::cout << "p=" << p20
//                  << "\tpc=" << ribs[2].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p20-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        p_test=s.f(u2, v1);
        TEST_ASSERT(tol.approximately_equal(p21, p_test));
//        std::cout << "p=" << p21
//                  << "\tpc=" << ribs[2].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p21-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u2, v);
        p_ref=ribs[2].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1-small_num, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1+small_num, v);
        p_ref=ribs[1].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u2, v);
        p_ref=ribs[2].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u2, v);
        p_ref=ribs[2].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1-small_num, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1+small_num, v);
        p_ref=ribs[1].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u2, v);
        p_ref=ribs[2].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[2].get_f(), "rib2", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_right_fp(), "rve1");
//          eli::octave_print(1, ribs[2].get_f(), ribs[2].get_left_fp(), "lve2");
//          eli::octave_finish(1);
//        }
      }

      // surface connecting ribs of variying order and segment count, specified continuous 1st & 2nd derivatives
      {
        index_type nsegs(2);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rc3, rs1, rs2, rs3, ra1, ra2, ra3;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, p20, p21, s00, s01, s10, s11, s20, s21, a00, a01, a10, a11, a20, a21;
        data_type u0(2), v0(1), u1(4), v1(5), u2(7), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        p20 << 5, 2, 5;
        p21 << 5, 5, 6;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;
        s20 << 1, 0, 1;
        s21 << 1, 0, 1;
        a00 << 0, 1, 0;
        a01 << 1, 0, 1;
        a10 << 1, 0, 1;
        a11 << 0, 1, 0;
        a20 << 0, 1, 0;
        a21 << 1, 0, 1;

        // create rib curves
        typename rib_curve_type::curve_type curve(3);
        typename rib_curve_type::control_point_type cp[4];

        // manually create cubic curves for ribs
        cp[0]=p00;
        cp[1] << 0, 2, 1;
        cp[2] << 1, 3, 1;
        cp[3]=p01;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc1.push_front(curve, v1-v0);
        rc1.set_t0(v0);

        cp[0]=p10;
        cp[1] << 1, 2, 3;
        cp[2] << 2, 3, 4;
        cp[3]=p11;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc2.push_front(curve, v1-v0);
        rc2.set_t0(v0);
        rc2.split(v0+(v1-v0)/4);

        cp[0]=p20;
        cp[1] << 4, 3, 5;
        cp[2] << 4, 4, 6;
        cp[3]=p21;
        for (index_type i=0; i<4; ++i)
        {
          curve.set_control_point(cp[i], i);
        }
        rc3.push_front(curve, v1-v0);
        rc3.set_t0(v0);
        rc3.split((v0+v1)/2);

        // create two two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_corner(p20, 0);
        plc.set_corner(p21, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc3);

        plc.set_corner(s00/2, 0);
        plc.set_corner(s01/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);
        rs1.split((v0+v1)/2);

        plc.set_corner(a00, 0);
        plc.set_corner(a01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra1);

        plc.set_corner(s10/2, 0);
        plc.set_corner(s11/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        plc.set_corner(a10, 0);
        plc.set_corner(a11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra2);
        ra2.split(v0+4*(v1-v0)/5);

        plc.set_corner(s20/2, 0);
        plc.set_corner(s21/2, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs3);

        plc.set_corner(a20, 0);
        plc.set_corner(a21, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(ra3);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[0].set_right_fpp(ra1);
        ribs[1].set_f(rc2);
        ribs[1].set_fp(rs2);
        ribs[1].set_fpp(ra2);
        ribs[2].set_f(rc3);
        ribs[2].set_left_fp(rs3);
        ribs[2].set_left_fpp(ra3);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_segment_du(u1-u0, 0);
        gc.set_segment_du(u2-u1, 1);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u2, v0);
        TEST_ASSERT(tol.approximately_equal(p20, p_test));
//        std::cout << "p=" << p20
//                  << "\tpc=" << ribs[2].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p20-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        p_test=s.f(u2, v1);
        TEST_ASSERT(tol.approximately_equal(p21, p_test));
//        std::cout << "p=" << p21
//                  << "\tpc=" << ribs[2].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p21-p_test).norm() << std::endl;
        v=v0+(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u2, v);
        p_ref=ribs[2].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u2, v);
        p_ref=ribs[2].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        v=v0+3*(v1-v0)/4;
        p_test=s.f(u0, v);
        p_ref=ribs[0].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u1, v);
        p_ref=ribs[1].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f(u2, v);
        p_ref=ribs[2].get_f().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u2, v);
        p_ref=ribs[2].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u0, v);
        p_ref=ribs[0].get_right_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u1, v);
        p_ref=ribs[1].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_uu(u2, v);
        p_ref=ribs[2].get_left_fpp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_print(1, ribs[2].get_f(), "rib2", true);
//          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
//          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_right_fp(), "rve1");
//          eli::octave_print(1, ribs[2].get_f(), ribs[2].get_left_fp(), "lve2");
//          eli::octave_finish(1);
//        }
      }
    }
};

#endif

