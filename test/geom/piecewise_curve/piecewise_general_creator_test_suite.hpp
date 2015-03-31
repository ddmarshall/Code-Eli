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

#ifndef piecewise_general_creator_test_suite_hpp
#define piecewise_general_creator_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_general_creator.hpp"
#include "eli/geom/intersect/minimum_distance_curve.hpp"

#include "octave_helpers.hpp"

template<typename data__>
class piecewise_general_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::piecewise_general_creator<data__, 3, tolerance_type> general_creator_type;
    typedef typename general_creator_type::joint_data joint_data_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_general_creator_test_suite<float>::create_joint_test);
      TEST_ADD(piecewise_general_creator_test_suite<float>::create_single_curve_test);
      TEST_ADD(piecewise_general_creator_test_suite<float>::create_multi_curve_no_dep_test);
      TEST_ADD(piecewise_general_creator_test_suite<float>::create_multi_curve_1st_deriv_dep_test);
      TEST_ADD(piecewise_general_creator_test_suite<float>::create_multi_curve_2nd_deriv_dep_test);
      TEST_ADD(piecewise_general_creator_test_suite<float>::create_single_curve_least_sq_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_general_creator_test_suite<double>::create_joint_test);
      TEST_ADD(piecewise_general_creator_test_suite<double>::create_single_curve_test);
      TEST_ADD(piecewise_general_creator_test_suite<double>::create_multi_curve_no_dep_test);
      TEST_ADD(piecewise_general_creator_test_suite<double>::create_multi_curve_1st_deriv_dep_test);
      TEST_ADD(piecewise_general_creator_test_suite<double>::create_multi_curve_2nd_deriv_dep_test);
      TEST_ADD(piecewise_general_creator_test_suite<double>::create_single_curve_least_sq_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_general_creator_test_suite<long double>::create_joint_test);
      TEST_ADD(piecewise_general_creator_test_suite<long double>::create_single_curve_test);
      TEST_ADD(piecewise_general_creator_test_suite<long double>::create_multi_curve_no_dep_test);
      TEST_ADD(piecewise_general_creator_test_suite<long double>::create_multi_curve_1st_deriv_dep_test);
      TEST_ADD(piecewise_general_creator_test_suite<long double>::create_multi_curve_2nd_deriv_dep_test);
      TEST_ADD(piecewise_general_creator_test_suite<long double>::create_single_curve_least_sq_test);
    }

  public:
    piecewise_general_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_general_creator_test_suite()
    {
    }

  private:
    void create_joint_test()
    {
      point_type p1, p2;
      bool rtn_flag;

      p1 << 1, 0, 0;
      p2 << 2, 2, 2;

      // set point
      {
        joint_data_type joint;

        // should start in bad state
        rtn_flag=joint.check_state();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C0);
        TEST_ASSERT(!joint.use_f());

        // add point
        rtn_flag=joint.set_f(p1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C0);
        TEST_ASSERT(joint.use_f());

        // unset point
        rtn_flag=joint.unset_f();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(!joint.use_f());
      }

      // set point and fp
      {
        joint_data_type joint;

        // should start in bad state
        rtn_flag=joint.check_state();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C0);
        TEST_ASSERT(!joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());

        // add left fp
        rtn_flag=joint.set_left_fp(p1);
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_left_fp()==p1);
        TEST_ASSERT(!joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());

        // add point
        rtn_flag=joint.set_f(p1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C0);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());

        // add right fp
        rtn_flag=joint.set_right_fp(p2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_right_fp()==p2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());

        // remove left fp
        rtn_flag=joint.unset_left_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());

        // remove left fp again
        rtn_flag=joint.unset_left_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());

        // remove right fp
        rtn_flag=joint.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());

        // add fp
        rtn_flag=joint.set_fp(p2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_left_fp()==p2);
        TEST_ASSERT(joint.get_right_fp()==p2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());

        // remove right fp
        rtn_flag=joint.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
      }

      // set point, fp and fpp
      {
        joint_data_type joint;

        // should start in bad state
        rtn_flag=joint.check_state();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C0);
        TEST_ASSERT(!joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // add left fpp
        rtn_flag=joint.set_left_fpp(p1);
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_left_fpp()==p1);
        TEST_ASSERT(!joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // add right fpp
        rtn_flag=joint.set_right_fpp(p2);
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_right_fpp()==p2);
        TEST_ASSERT(!joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(joint.use_right_fpp());

        // add point
        rtn_flag=joint.set_f(p1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C0);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(joint.use_right_fpp());

        // remove right fpp
        rtn_flag=joint.unset_right_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fp());

        // remove right fpp again
        rtn_flag=joint.unset_right_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fp());

        // remove left fpp
        rtn_flag=joint.unset_left_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // add fp
        rtn_flag=joint.set_fp(p2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_left_fp()==p2);
        TEST_ASSERT(joint.get_right_fp()==p2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // add fpp
        rtn_flag=joint.set_fpp(p1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_left_fpp()==p1);
        TEST_ASSERT(joint.get_right_fpp()==p1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(joint.use_right_fpp());

        // remove right fp
        rtn_flag=joint.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(joint.use_right_fpp());

        // remove right fpp
        rtn_flag=joint.unset_right_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // remove everything
        rtn_flag=joint.unset_fpp();
        TEST_ASSERT(rtn_flag);
        rtn_flag=joint.unset_fp();
        TEST_ASSERT(rtn_flag);
        rtn_flag=joint.unset_f();
        TEST_ASSERT(!rtn_flag);

        // add everything to right only
        rtn_flag=joint.set_f(p2);
        TEST_ASSERT(rtn_flag);
        rtn_flag=joint.set_right_fp(p1);
        TEST_ASSERT(rtn_flag);
        rtn_flag=joint.set_right_fpp(p2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.get_f()==p2);
        TEST_ASSERT(joint.get_right_fp()==p1);
        TEST_ASSERT(joint.get_right_fpp()==p2);
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(joint.use_right_fpp());
      }

      // set conditions with C1
      {
        joint_data_type joint;

        // start with point
        rtn_flag=joint.set_f(p1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C0);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // set continuity to C1
        rtn_flag=joint.set_continuity(general_creator_type::C1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // add left fp
        rtn_flag=joint.set_left_fp(p2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // unset left fp
        rtn_flag=joint.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        rtn_flag=joint.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // add right fp
        rtn_flag=joint.set_right_fp(p2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // change continuity
        rtn_flag=joint.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        rtn_flag=joint.set_continuity(general_creator_type::C0);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C0);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // add fp
        rtn_flag=joint.set_fp(p1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C0);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // change continuity
        rtn_flag=joint.set_continuity(general_creator_type::C1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // unset left fp
        rtn_flag=joint.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // set left fp to change left and right value
        rtn_flag=joint.set_left_fp(p2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_left_fp()==p2);
        TEST_ASSERT(joint.get_right_fp()==p2);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // set left fp
        rtn_flag=joint.set_left_fp(joint.get_right_fp());
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());
      }

      // set point, fp, fpp and C2
      {
        joint_data_type joint;

        // start with point and fp
        rtn_flag=joint.set_f(p1);
        TEST_ASSERT(rtn_flag);
        rtn_flag=joint.set_fp(p2);
        TEST_ASSERT(rtn_flag);
        rtn_flag=joint.set_continuity(general_creator_type::C1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // set continuity to C2
        rtn_flag=joint.set_continuity(general_creator_type::C2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // unset left fp
        rtn_flag=joint.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // unset right fp
        rtn_flag=joint.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // set fpp
        rtn_flag=joint.set_fpp(p1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(joint.use_right_fpp());

        // unset left fpp
        rtn_flag=joint.unset_left_fpp();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(joint.use_right_fpp());

        // set left fpp to set both left and right values
        rtn_flag=joint.set_left_fpp(p2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(joint.use_right_fpp());

        // set left fpp
        rtn_flag=joint.set_left_fpp(joint.get_right_fpp());
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(joint.use_right_fpp());

        // set fp
        rtn_flag=joint.set_fp(p2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C2);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(joint.use_left_fpp());
        TEST_ASSERT(joint.use_right_fpp());
      }

      // copy constructor, equivalence and assignment operators
      {
        joint_data_type j1, j2, j3;

        // defaults should be the same
        TEST_ASSERT(j1==j2);

        // manually build identical ones
        j1.set_f(p2);
        j1.set_left_fp(p1);
        j1.set_right_fpp(p2);
        j2.set_f(p2);
        j2.set_left_fp(p1);
        j2.set_right_fpp(p2);
        TEST_ASSERT(j1==j2);

        // assignment operator
        j3=j2;
        TEST_ASSERT(j3==j1);
        TEST_ASSERT(j3==j2);

        // set value for j2's fpp and then don't use it
        j1.unset_fpp();
        j2.unset_fpp();
        j3.unset_fpp();
        j2.set_fpp(p1);
        j2.unset_fpp();
        TEST_ASSERT(j3==j1);
        TEST_ASSERT(j3==j2);

        // copy constructor
        joint_data_type j4(j3);
        TEST_ASSERT(j4==j3);
      }
    }

    void create_single_curve_test()
    {
      // simple line connecting 2 points
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set the maximum degrees of each segment
        max_degree[0]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }

      // simple 2nd degree curve with 1st derivative specified
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 1st derivative
        p << 1, -1, 1;
        joints[0].set_right_fp(p);

        // set the maximum degrees of each segment
        max_degree[0]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }

      // simple 2nd degree curve with 2nd derivative specified
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 2nd derivative
        p << 1, -1, 1;
        joints[0].set_right_fpp(p);

        // set the maximum degrees of each segment
        max_degree[0]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fpp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fpp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fpp()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }

      // simple 3rd degree curve with both 1st derivatives specified
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 1st derivatives
        p << 1, -1, 1;
        joints[0].set_right_fp(p);
        p << 1, 0, 0;
        joints[1].set_left_fp(p);

        // set the maximum degrees of each segment to be too small
        max_degree[0]=2;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(!rtn_flag);

        // set the maximum degrees of each segment
        max_degree[0]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }

      // simple 5th degree curve with both 1st and 2nd derivatives specified
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 1st derivatives
        p << 1, -1, 1;
        joints[0].set_right_fp(p);
        p << 1, 0, 0;
        joints[1].set_left_fp(p);

        // set joint 2nd derivatives
        p << 0, -1, 0;
        joints[0].set_right_fpp(p);
        p << 0, 1, 0;
        joints[1].set_left_fpp(p);

        // set the maximum degrees of each segment
        max_degree[0]=6;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fpp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fpp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fpp()-p_test).norm() << std::endl;
        p_test=c.fpp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fpp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fpp()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }
    }

    void create_multi_curve_no_dep_test()
    {
      // create two simple segments
      {
        index_type nsegs(2);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 1, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);

        // set joint 1st derivatives
        p << 0, 1, 1;
        joints[0].set_right_fp(p);
        p << 0, 0, 1;
        joints[1].set_left_fp(p);
        p << 1, 0, 1;
        joints[1].set_right_fp(p);
        p << 0, 1, 1;
        joints[2].set_left_fp(p);

        // set the maximum degrees of each segment to be too small
        max_degree[0]=4;
        max_degree[1]=2;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        gc.set_segment_dt(t[1]-t[0], 0);
        gc.set_segment_dt(t[2]-t[1], 1);
        rtn_flag=gc.create(c);
        TEST_ASSERT(!rtn_flag);

        // set the maximum degrees of each segment
        max_degree[0]=4;
        max_degree[1]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        gc.set_segment_dt(t[1]-t[0], 0);
        gc.set_segment_dt(t[2]-t[1], 1);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_right_fp(), p_test));
//        std::cout << "p=" << joints[1].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }

      // create five simple segments
      {
        index_type i, nsegs(5);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;
        t[3]=8;
        t[4]=10;
        t[5]=12;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 1, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);
        p << 1, 1, 2;
        joints[3].set_f(p);
        p << 0, 2, 1;
        joints[4].set_f(p);
        p << 3, 3, 3;
        joints[5].set_f(p);

        // set joint 1st derivatives
        p << 0, 1, 1;
        joints[0].set_right_fp(p);
        p << 0, 0, 1;
        joints[1].set_left_fp(p);
        p << 1, 0, 1;
        joints[1].set_right_fp(p);
        p << 0, 1, 1;
        joints[2].set_left_fp(p);
        p << 1, 0, 1;
        joints[2].set_right_fp(p);
        p << 0, 1, 1;
        joints[3].set_left_fp(p);
        p << 0, 1, 1;
        joints[3].set_right_fp(p);
        p << 0, 0, 1;
        joints[4].set_left_fp(p);
        p << 1, 0, 1;
        joints[4].set_right_fp(p);
        p << 0, 0, 1;
        joints[5].set_left_fp(p);

        // set the maximum degrees of each segment
        max_degree[0]=4;
        max_degree[1]=4;
        max_degree[2]=4;
        max_degree[3]=4;
        max_degree[4]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        for (i=0; i<nsegs; ++i)
        {
          gc.set_segment_dt(t[i+1]-t[i], i);
        }
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[3]);
        TEST_ASSERT(tol.approximately_equal(joints[3].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[4]);
        TEST_ASSERT(tol.approximately_equal(joints[4].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[5]);
        TEST_ASSERT(tol.approximately_equal(joints[5].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_right_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[2]-small_num);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[2]+small_num);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_right_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[3]-small_num);
        TEST_ASSERT(tol.approximately_equal(joints[3].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[3]+small_num);
        TEST_ASSERT(tol.approximately_equal(joints[3].get_right_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[4]-small_num);
        TEST_ASSERT(tol.approximately_equal(joints[4].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[4]+small_num);
        TEST_ASSERT(tol.approximately_equal(joints[4].get_right_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[5]);
        TEST_ASSERT(tol.approximately_equal(joints[5].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }
    }

    void create_multi_curve_1st_deriv_dep_test()
    {
      // create two simple segments with 1st derivative continuity
      {
        index_type nsegs(2);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 1, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);

        // set joint 1st derivative smoothness
        joints[1].set_continuity(general_creator_type::C1);

        // set the maximum degrees of each segment
        max_degree[0]=4;
        max_degree[1]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        gc.set_segment_dt(t[1]-t[0], 0);
        gc.set_segment_dt(t[2]-t[1], 1);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_test1;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[2].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[2].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        p_test1=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }

      // create three segments where 1st interior joint has discontinuous, specified 1st derivative
      // and 2nd has continuous 1st derivative
      {
        index_type nsegs(3);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;
        t[3]=9;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 1, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);
        p << 1, 1, 2;
        joints[3].set_f(p);

        // set joint 1st derivative smoothness

        // set joint 1st derivatives
        p << 0, 1, 1;
        joints[0].set_right_fp(p);
        p << 0, 0, 1;
        joints[1].set_left_fp(p);
        p << -1, 0, 1;
        joints[1].set_right_fp(p);
        joints[2].set_continuity(general_creator_type::C1);
        p << 0, 1, 1;
        joints[3].set_left_fp(p);

        // set the maximum degrees of each segment
        max_degree[0]=4;
        max_degree[1]=4;
        max_degree[2]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        gc.set_segment_dt(t[1]-t[0], 0);
        gc.set_segment_dt(t[2]-t[1], 1);
        gc.set_segment_dt(t[3]-t[2], 2);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_test1;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[2].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[2].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[3]);
        TEST_ASSERT(tol.approximately_equal(joints[3].get_f(), p_test));
//        std::cout << "p=" << joints[3].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[3].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_right_fp(), p_test));
//        std::cout << "p=" << joints[1].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_right_fp()-p_test).norm() << std::endl;

        p_test=c.fp(t[2]-small_num);
        p_test1=c.fp(t[2]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//          octave_print(1, c);
      }

      // create three segments where 1st interior joint has discontinuous, specified 1st derivative
      // and 2nd has one side slope specified and other is continuous 1st derivative
      {
        index_type nsegs(3);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;
        t[3]=9;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 1, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);
        p << 1, 1, 2;
        joints[3].set_f(p);

        // set joint 1st derivative smoothness

        // set joint 1st derivatives
        p << 0, 1, 1;
        joints[0].set_right_fp(p);
        p << 0, 0, 1;
        joints[1].set_left_fp(p);
        p << -1, 0, 1;
        joints[1].set_right_fp(p);
        joints[2].set_continuity(general_creator_type::C1);
        p << 0, 1, 0;
        joints[2].set_right_fp(p);
        p << 0, 1, 1;
        joints[3].set_left_fp(p);

        // set the maximum degrees of each segment
        max_degree[0]=4;
        max_degree[1]=4;
        max_degree[2]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        gc.set_segment_dt(t[1]-t[0], 0);
        gc.set_segment_dt(t[2]-t[1], 1);
        gc.set_segment_dt(t[3]-t[2], 2);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_test1;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[2].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[2].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[3]);
        TEST_ASSERT(tol.approximately_equal(joints[3].get_f(), p_test));
//        std::cout << "p=" << joints[3].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[3].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_right_fp(), p_test));
//        std::cout << "p=" << joints[1].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_right_fp()-p_test).norm() << std::endl;

        p_test=c.fp(t[2]-small_num);
        p_test1=c.fp(t[2]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
        TEST_ASSERT(tol.approximately_equal(joints[2].get_right_fp(), p_test));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//          octave_print(1, c);
      }

      // create five segments with all interior joints continuous 1st derivatives
      {
        index_type i, nsegs(5);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;
        t[3]=8;
        t[4]=10;
        t[5]=12;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 1, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);
        p << 1, 1, 2;
        joints[3].set_f(p);
        p << 0, 2, 1;
        joints[4].set_f(p);
        p << 3, 3, 3;
        joints[5].set_f(p);

        // set joint 1st derivatives
        joints[1].set_continuity(general_creator_type::C1);
        joints[2].set_continuity(general_creator_type::C1);
        joints[3].set_continuity(general_creator_type::C1);
        joints[4].set_continuity(general_creator_type::C1);

        // set the maximum degrees of each segment
        max_degree[0]=4;
        max_degree[1]=4;
        max_degree[2]=4;
        max_degree[3]=4;
        max_degree[4]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        for (i=0; i<nsegs; ++i)
        {
          gc.set_segment_dt(t[i+1]-t[i], i);
        }
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_test1;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[3]);
        TEST_ASSERT(tol.approximately_equal(joints[3].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[4]);
        TEST_ASSERT(tol.approximately_equal(joints[4].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[5]);
        TEST_ASSERT(tol.approximately_equal(joints[5].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        p_test1=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[2]-small_num);
        p_test1=c.fp(t[2]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[3]-small_num);
        p_test1=c.fp(t[3]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[4]-small_num);
        p_test1=c.fp(t[4]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//          octave_print(1, c);
      }
    }
    void create_multi_curve_2nd_deriv_dep_test()
    {
      // create two simple segments with 2nd derivative continuity
      {
        index_type nsegs(2);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 1, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);

        // set joint 2nd derivative smoothness
        joints[1].set_continuity(general_creator_type::C2);

        // set the maximum degrees of each segment
        max_degree[0]=4;
        max_degree[1]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        gc.set_segment_dt(t[1]-t[0], 0);
        gc.set_segment_dt(t[2]-t[1], 1);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_test1;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[2].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[2].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        p_test1=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[1]-small_num);
        p_test1=c.fpp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }

      // create one line segment and another segment with 2nd derivative continuity
      {
        index_type nsegs(2);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 1, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);

        // set joint 2nd derivative smoothness
        joints[1].set_continuity(general_creator_type::C2);

        // set the maximum degrees of each segment
        max_degree[0]=1;
        max_degree[1]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        gc.set_segment_dt(t[1]-t[0], 0);
        gc.set_segment_dt(t[2]-t[1], 1);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_test1;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[2].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[2].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        p_test1=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[1]-small_num);
        p_test1=c.fpp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//          octave_print(1, c);
      }

      // create one line segment and another segment with 2nd derivative continuity
      {
        index_type nsegs(2);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 1, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);

        // set joint 2nd derivative smoothness
        joints[1].set_continuity(general_creator_type::C2);

        // set the maximum degrees of each segment
        max_degree[0]=4;
        max_degree[1]=1;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        gc.set_segment_dt(t[1]-t[0], 0);
        gc.set_segment_dt(t[2]-t[1], 1);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_test1;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[2].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[2].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        p_test1=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[1]-small_num);
        p_test1=c.fpp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//          octave_print(1, c);
      }

      // create five segments with all interior joints continuous 2nd derivatives
      {
        index_type i, nsegs(5);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;
        t[3]=8;
        t[4]=10;
        t[5]=12;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 1, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);
        p << 1, 1, 2;
        joints[3].set_f(p);
        p << 0, 2, 1;
        joints[4].set_f(p);
        p << 3, 3, 3;
        joints[5].set_f(p);

        // set joint 1st derivatives
        joints[1].set_continuity(general_creator_type::C2);
        joints[2].set_continuity(general_creator_type::C2);
        joints[3].set_continuity(general_creator_type::C2);
        joints[4].set_continuity(general_creator_type::C2);

        // set the maximum degrees of each segment
        max_degree[0]=4;
        max_degree[1]=4;
        max_degree[2]=4;
        max_degree[3]=4;
        max_degree[4]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        for (i=0; i<nsegs; ++i)
        {
          gc.set_segment_dt(t[i+1]-t[i], i);
        }
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_test1;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[3]);
        TEST_ASSERT(tol.approximately_equal(joints[3].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[4]);
        TEST_ASSERT(tol.approximately_equal(joints[4].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[5]);
        TEST_ASSERT(tol.approximately_equal(joints[5].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        p_test1=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[2]-small_num);
        p_test1=c.fp(t[2]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[3]-small_num);
        p_test1=c.fp(t[3]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[4]-small_num);
        p_test1=c.fp(t[4]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[1]-small_num);
        p_test1=c.fpp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[2]-small_num);
        p_test1=c.fpp(t[2]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[3]-small_num);
        p_test1=c.fpp(t[3]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[4]-small_num);
        p_test1=c.fpp(t[4]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//          octave_print(1, c);
      }

      // create five segments that are colinear and enforce 2nd derivative smoothness
      {
        index_type i, nsegs(5);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;
        t[3]=8;
        t[4]=10;
        t[5]=12;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 1, 1, 1;
        joints[1].set_f(p);
        p << 2, 2, 2;
        joints[2].set_f(p);
        p << 4, 4, 4;
        joints[3].set_f(p);
        p << 5, 5, 5;
        joints[4].set_f(p);
        p << 7, 7, 7;
        joints[5].set_f(p);

        // set joint 1st derivatives
        joints[1].set_continuity(general_creator_type::C2);
        joints[2].set_continuity(general_creator_type::C2);
        joints[3].set_continuity(general_creator_type::C2);
        joints[4].set_continuity(general_creator_type::C2);

        // set the maximum degrees of each segment
        max_degree[0]=4;
        max_degree[1]=4;
        max_degree[2]=4;
        max_degree[3]=4;
        max_degree[4]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        for (i=0; i<nsegs; ++i)
        {
          gc.set_segment_dt(t[i+1]-t[i], i);
        }
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_test1;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[3]);
        TEST_ASSERT(tol.approximately_equal(joints[3].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[4]);
        TEST_ASSERT(tol.approximately_equal(joints[4].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[5]);
        TEST_ASSERT(tol.approximately_equal(joints[5].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        p_test1=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[2]-small_num);
        p_test1=c.fp(t[2]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[3]-small_num);
        p_test1=c.fp(t[3]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[4]-small_num);
        p_test1=c.fp(t[4]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
       p_test=c.fpp(t[1]-small_num);
        p_test1=c.fpp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[2]-small_num);
        p_test1=c.fpp(t[2]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[3]-small_num);
        p_test1=c.fpp(t[3]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[4]-small_num);
        p_test1=c.fpp(t[4]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//          octave_print(1, c);
      }

      // create five segments with a variety of continuities
      {
        index_type i, nsegs(5);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        std::vector<data_type> t(nsegs+1);
        bool rtn_flag;

        // set the times
        t[0]=2;
        t[1]=4;
        t[2]=7;
        t[3]=8;
        t[4]=10;
        t[5]=12;

        // set the joints
        p << 0, 0, 0;
        joints[0].set_f(p);
        p << 0, 2, 0;
        joints[1].set_f(p);
        p << 1, 1, 1;
        joints[2].set_f(p);
        p << 1, 1, 2;
        joints[3].set_f(p);
        p << 0, 2, 1;
        joints[4].set_f(p);
        p << 3, 3, 3;
        joints[5].set_f(p);

        // set joint 1st derivatives
        joints[1].set_continuity(general_creator_type::C1);
        joints[2].set_continuity(general_creator_type::C0);
        joints[3].set_continuity(general_creator_type::C2);
        joints[4].set_continuity(general_creator_type::C2);

        // set the maximum degrees of each segment
        max_degree[0]=0;
        max_degree[1]=1;
        max_degree[2]=0;
        max_degree[3]=0;
        max_degree[4]=0;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t[0]);
        for (i=0; i<nsegs; ++i)
        {
          gc.set_segment_dt(t[i+1]-t[i], i);
        }
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_test1;
        data_type small_num(100*std::numeric_limits<data_type>::epsilon());

        p_test=c.f(t[0]);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[1]);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[2]);
        TEST_ASSERT(tol.approximately_equal(joints[2].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[3]);
        TEST_ASSERT(tol.approximately_equal(joints[3].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[4]);
        TEST_ASSERT(tol.approximately_equal(joints[4].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t[5]);
        TEST_ASSERT(tol.approximately_equal(joints[5].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t[1]-small_num);
        p_test1=c.fp(t[1]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[2]-small_num);
        p_test1=c.fp(t[2]+small_num);
        TEST_ASSERT(!tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[3]-small_num);
        p_test1=c.fp(t[3]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fp(t[4]-small_num);
        p_test1=c.fp(t[4]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
       p_test=c.fpp(t[1]-small_num);
        p_test1=c.fpp(t[1]+small_num);
        TEST_ASSERT(!tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[2]-small_num);
        p_test1=c.fpp(t[2]+small_num);
        TEST_ASSERT(!tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[3]-small_num);
        p_test1=c.fpp(t[3]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;
        p_test=c.fpp(t[4]-small_num);
        p_test1=c.fpp(t[4]+small_num);
        TEST_ASSERT(tol.approximately_equal(p_test, p_test1));
//        std::cout << "p_test=" << p_test
//                  << "\tp_test1=" << p_test1 << "\tdiff="
//                  << (p_test-p_test1).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//          octave_print(1, c);
      }
    }

    void create_single_curve_least_sq_test()
    {
      // simple min 2nd degree curve connecting 2 points and fitting 3 points with no free degrees to fit
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        std::vector<typename general_creator_type::fit_data> fit_points(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set the fit points
        p << 0.5, 0.5, 0.5;
        fit_points[0].add_point(p);
        p << 0.25, 0.25, 0.25;
        fit_points[0].add_point(p);
        p << -0.25, -0.5, 0;
        fit_points[0].add_point(p);

        // set the maximum degrees of each segment to be the degree needed to satisfy the joints (1)
        max_degree[0]=1;

        // create curve
        rtn_flag=gc.set_conditions(joints, fit_points, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;
        std::vector<index_type> c_degree(nsegs);

        c.degrees(c_degree.begin());
        TEST_ASSERT(c_degree[0]==max_degree[0]);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, fit_points[0].get_point(0), "fp0");
//          eli::test::octave_print(1, fit_points[0].get_point(1), "fp1");
//          eli::test::octave_print(1, fit_points[0].get_point(2), "fp2");
//          eli::test::octave_print(1, c, "piecewise");
//          eli::test::octave_finish(1);
//        }

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
      }

      // simple min 2nd degree curve connecting 2 points and fitting 3 points with max_degree too low
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        std::vector<typename general_creator_type::fit_data> fit_points(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set the fit points
        p << 0.5, 0.5, 0.5;
        fit_points[0].add_point(p);
        p << 0.25, 0.25, 0.25;
        fit_points[0].add_point(p);
        p << -0.25, -0.5, 0;
        fit_points[0].add_point(p);

        // set the maximum degrees of each segment to be higher than degree needed to
        // go through all fit points (4)
        max_degree[0]=3;

        // create curve
        rtn_flag=gc.set_conditions(joints, fit_points, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;
        std::vector<index_type> c_degree(nsegs);

        c.degrees(c_degree.begin());
        TEST_ASSERT(c_degree[0]==max_degree[0]);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, fit_points[0].get_point(0), "fp0");
//          eli::test::octave_print(1, fit_points[0].get_point(1), "fp1");
//          eli::test::octave_print(1, fit_points[0].get_point(2), "fp2");
//          eli::test::octave_print(1, c, "piecewise");
//          eli::test::octave_finish(1);
//        }

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
      }

      // simple min 2nd degree curve connecting 2 points and fitting 3 points with max_degree just right
      {
        index_type i, nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        std::vector<typename general_creator_type::fit_data> fit_points(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set the fit points
        p << 0.5, 0.5, 0.5;
        fit_points[0].add_point(p);
        p << 0.25, 0.25, 0.25;
        fit_points[0].add_point(p);
        p << -0.25, -0.5, 0;
        fit_points[0].add_point(p);

        // set the maximum degrees of each segment to be the same degree needed to
        // go through all fit points (4)
        max_degree[0]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, fit_points, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;
        data_type dist, dist_t;
        std::vector<index_type> c_degree(nsegs);

        c.degrees(c_degree.begin());
        TEST_ASSERT(c_degree[0]==max_degree[0]);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, fit_points[0].get_point(0), "fp0");
//          eli::test::octave_print(1, fit_points[0].get_point(1), "fp1");
//          eli::test::octave_print(1, fit_points[0].get_point(2), "fp2");
//          eli::test::octave_print(1, c, "piecewise");
//          eli::test::octave_finish(1);
//        }

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;

        i=0;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
        i=1;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
        i=2;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
      }

      // simple min 2nd degree curve connecting 2 points and fitting 3 points with max_degree too high
      {
        index_type i, nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        std::vector<typename general_creator_type::fit_data> fit_points(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set the fit points
        p << 0.5, 0.5, 0.5;
        fit_points[0].add_point(p);
        p << 0.25, 0.25, 0.25;
        fit_points[0].add_point(p);
        p << -0.25, -0.5, 0;
        fit_points[0].add_point(p);

        // set the maximum degrees of each segment to be higher than degree needed to
        // go through all fit points (4)
        max_degree[0]=5;

        // create curve
        rtn_flag=gc.set_conditions(joints, fit_points, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;
        data_type dist, dist_t;
        std::vector<index_type> c_degree(nsegs);

        c.degrees(c_degree.begin());
        TEST_ASSERT(c_degree[0]<max_degree[0]);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, fit_points[0].get_point(0), "fp0");
//          eli::test::octave_print(1, fit_points[0].get_point(1), "fp1");
//          eli::test::octave_print(1, fit_points[0].get_point(2), "fp2");
//          eli::test::octave_print(1, c, "piecewise");
//          eli::test::octave_finish(1);
//        }

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;

        i=0;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
        i=1;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
        i=2;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
      }

      // simple min 5th degree curve with both 1st and 2nd derivatives specified with no free degrees to fit
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        std::vector<typename general_creator_type::fit_data> fit_points(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 1st derivatives
        p << 1, -1, 1;
        joints[0].set_right_fp(p);
        p << 1, 0, 0;
        joints[1].set_left_fp(p);

        // set joint 2nd derivatives
        p << 0, -1, 0;
        joints[0].set_right_fpp(p);
        p << 0, 1, 0;
        joints[1].set_left_fpp(p);

        // set the fit points
        p << 0.5, 0.5, 0.5;
        fit_points[0].add_point(p);
        p << 0.25, 0.25, 0.25;
        fit_points[0].add_point(p);
        p << -0.25, -0.5, 0;
        fit_points[0].add_point(p);

        // set the maximum degrees of each segment to be the degree needed to satisfy the joints (5)
        max_degree[0]=5;

        // create curve
        rtn_flag=gc.set_conditions(joints, fit_points, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;
        std::vector<index_type> c_degree(nsegs);

        c.degrees(c_degree.begin());
        TEST_ASSERT(c_degree[0]==max_degree[0]);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, fit_points[0].get_point(0), "fp0");
//          eli::test::octave_print(1, fit_points[0].get_point(1), "fp1");
//          eli::test::octave_print(1, fit_points[0].get_point(2), "fp2");
//          eli::test::octave_print(1, c, "piecewise");
//          eli::test::octave_finish(1);
//        }

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fpp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fpp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fpp()-p_test).norm() << std::endl;
        p_test=c.fpp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fpp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fpp()-p_test).norm() << std::endl;
      }

      // simple min 5th degree curve with both 1st and 2nd derivatives specified with max_degree too low
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        std::vector<typename general_creator_type::fit_data> fit_points(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 1st derivatives
        p << 1, -1, 1;
        joints[0].set_right_fp(p);
        p << 1, 0, 0;
        joints[1].set_left_fp(p);

        // set joint 2nd derivatives
        p << 0, -1, 0;
        joints[0].set_right_fpp(p);
        p << 0, 1, 0;
        joints[1].set_left_fpp(p);

        // set the fit points
        p << 0.5, 0.5, 0.5;
        fit_points[0].add_point(p);
        p << 0.25, 0.25, 0.25;
        fit_points[0].add_point(p);
        p << -0.25, -0.5, 0;
        fit_points[0].add_point(p);

        // set the maximum degrees of each segment to be lower than degree needed to
        // go through all fit points (8)
        max_degree[0]=7;

        // create curve
        rtn_flag=gc.set_conditions(joints, fit_points, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, fit_points[0].get_point(0), "fp0");
//          eli::test::octave_print(1, fit_points[0].get_point(1), "fp1");
//          eli::test::octave_print(1, fit_points[0].get_point(2), "fp2");
//          eli::test::octave_print(1, c, "piecewise");
//          eli::test::octave_finish(1);
//        }

        // test the resulting curve
        point_type p_test;
        std::vector<index_type> c_degree(nsegs);

        c.degrees(c_degree.begin());
        TEST_ASSERT(c_degree[0]==max_degree[0]);

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fpp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fpp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fpp()-p_test).norm() << std::endl;
        p_test=c.fpp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fpp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fpp()-p_test).norm() << std::endl;
      }

      // simple min 5th degree curve with both 1st and 2nd derivatives specified with max_degree just right
      {
        index_type i, nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        std::vector<typename general_creator_type::fit_data> fit_points(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 1st derivatives
        p << 1, -1, 1;
        joints[0].set_right_fp(p);
        p << 1, 0, 0;
        joints[1].set_left_fp(p);

        // set joint 2nd derivatives
        p << 0, -1, 0;
        joints[0].set_right_fpp(p);
        p << 0, 1, 0;
        joints[1].set_left_fpp(p);

        // set the fit points
        p << 0.5, 0.5, 0.5;
        fit_points[0].add_point(p);
        p << 0.25, 0.25, 0.25;
        fit_points[0].add_point(p);
        p << -0.25, -0.5, 0;
        fit_points[0].add_point(p);

        // set the maximum degrees of each segment to be the same degree needed to
        // go through all fit points (8)
        max_degree[0]=8;

        // create curve
        rtn_flag=gc.set_conditions(joints, fit_points, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;
        data_type dist, dist_t;
        std::vector<index_type> c_degree(nsegs);

        c.degrees(c_degree.begin());
        TEST_ASSERT(c_degree[0]==max_degree[0]);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, fit_points[0].get_point(0), "fp0");
//          eli::test::octave_print(1, fit_points[0].get_point(1), "fp1");
//          eli::test::octave_print(1, fit_points[0].get_point(2), "fp2");
//          eli::test::octave_print(1, c, "piecewise");
//          eli::test::octave_finish(1);
//        }

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fpp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fpp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fpp()-p_test).norm() << std::endl;
        p_test=c.fpp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fpp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fpp()-p_test).norm() << std::endl;

        i=0;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
        i=1;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
        i=2;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
      }

      // simple min 5th degree curve with both 1st and 2nd derivatives specified with max_degree too high
      {
        index_type i, nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        std::vector<typename general_creator_type::fit_data> fit_points(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 1st derivatives
        p << 1, -1, 1;
        joints[0].set_right_fp(p);
        p << 1, 0, 0;
        joints[1].set_left_fp(p);

        // set joint 2nd derivatives
        p << 0, -1, 0;
        joints[0].set_right_fpp(p);
        p << 0, 1, 0;
        joints[1].set_left_fpp(p);

        // set the fit points
        p << 0.5, 0.5, 0.5;
        fit_points[0].add_point(p);
        p << 0.25, 0.25, 0.25;
        fit_points[0].add_point(p);
        p << -0.25, -0.5, 0;
        fit_points[0].add_point(p);

        // set the maximum degrees of each segment to be higher than degree needed to
        // go through all fit points (8)
        max_degree[0]=10;

        // create curve
        rtn_flag=gc.set_conditions(joints, fit_points, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;
        data_type dist, dist_t;
        std::vector<index_type> c_degree(nsegs);

        c.degrees(c_degree.begin());
        TEST_ASSERT(c_degree[0]<max_degree[0]);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          eli::test::octave_print(1, fit_points[0].get_point(0), "fp0");
//          eli::test::octave_print(1, fit_points[0].get_point(1), "fp1");
//          eli::test::octave_print(1, fit_points[0].get_point(2), "fp2");
//          eli::test::octave_print(1, c, "piecewise");
//          eli::test::octave_finish(1);
//        }

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fpp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fpp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fpp()-p_test).norm() << std::endl;
        p_test=c.fpp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fpp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fpp()-p_test).norm() << std::endl;

        i=0;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
        i=1;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
        i=2;
        dist=eli::geom::intersect::minimum_distance(dist_t, c, fit_points[0].get_point(i));
        TEST_ASSERT(tol.approximately_equal(dist, 0));
//        std::cout << "dist=" << dist << std::endl;
      }
    }
};

#endif

