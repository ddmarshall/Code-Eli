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

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_general_creator.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

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
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_general_creator_test_suite<double>::create_joint_test);
      TEST_ADD(piecewise_general_creator_test_suite<double>::create_single_curve_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_general_creator_test_suite<long double>::create_joint_test);
      TEST_ADD(piecewise_general_creator_test_suite<long double>::create_single_curve_test);
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
    void octave_print(int figno, const piecewise_curve_type &pc) const
    {
      index_type i, pp, ns;
      data_type tmin, tmax;

      ns=pc.number_segments();
      pc.get_parameter_min(tmin);
      pc.get_parameter_max(tmax);

      std::cout << "figure(" << figno << ");" << std::endl;

      // get control points and print
      std::cout << "cp_x=[";
      for (pp=0; pp<ns; ++pp)
      {
        curve_type bez;
        pc.get(bez, pp);
        for (i=0; i<=bez.degree(); ++i)
        {
          std::cout << bez.get_control_point(i).x();
          if (i<bez.degree())
            std::cout << ", ";
          else if (pp<ns-1)
            std::cout << "; ";
        }
        std::cout << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "cp_y=[";
      for (pp=0; pp<ns; ++pp)
      {
        curve_type bez;
        pc.get(bez, pp);
        for (i=0; i<=bez.degree(); ++i)
        {
          std::cout << bez.get_control_point(i).y();
          if (i<bez.degree())
            std::cout << ", ";
          else if (pp<ns-1)
            std::cout << "; ";
        }
        std::cout << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "cp_z=[";
      for (pp=0; pp<ns; ++pp)
      {
        curve_type bez;
        pc.get(bez, pp);
        for (i=0; i<=bez.degree(); ++i)
        {
          std::cout << bez.get_control_point(i).z();
          if (i<bez.degree())
            std::cout << ", ";
          else if (pp<ns-1)
            std::cout << "; ";
        }
        std::cout << std::endl;
      }
      std::cout << "];" << std::endl;

      // initialize the t parameters
      std::vector<data__> t(129);
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(t.size()-1);
      }

      // set the surface points
      std::cout << "surf_x=[";
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        std::cout << pc.f(t[i]).x();
        if (i<static_cast<index_type>(t.size()-1))
          std::cout << ", ";
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_y=[";
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        std::cout << pc.f(t[i]).y();
        if (i<static_cast<index_type>(t.size()-1))
          std::cout << ", ";
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_z=[";
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        std::cout << pc.f(t[i]).z();
        if (i<static_cast<index_type>(t.size()-1))
          std::cout << ", ";
      }
      std::cout << "];" << std::endl;

      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot3(surf_x, surf_y, surf_z, '-k');" << std::endl;
      std::cout << "hold on;" << std::endl;
      std::cout << "plot3(cp_x', cp_y', cp_z', '-ok', 'MarkerFaceColor', [0 0 0]);" << std::endl;
      std::cout << "hold off;" << std::endl;
    }

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
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // unset left fp
        rtn_flag=joint.unset_left_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(!joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // add right fp
        rtn_flag=joint.set_right_fp(p2);
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(joint.get_continuity()==general_creator_type::C1);
        TEST_ASSERT(joint.use_f());
        TEST_ASSERT(!joint.use_left_fp());
        TEST_ASSERT(joint.use_right_fp());
        TEST_ASSERT(!joint.use_left_fpp());
        TEST_ASSERT(!joint.use_right_fpp());

        // change continuity
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

        // set left fp to discontinuous value
        rtn_flag=joint.set_left_fp(p2);
        TEST_ASSERT(!rtn_flag);
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

        // set left fpp to discontinuous value
        rtn_flag=joint.set_left_fpp(p2);
        TEST_ASSERT(!rtn_flag);
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
};

#endif

