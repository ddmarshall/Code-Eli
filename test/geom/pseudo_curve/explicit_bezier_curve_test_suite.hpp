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

#ifndef explicit_bezier_curve_test_suite_hpp
#define explicit_bezier_curve_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/curve/length.hpp"
#include "eli/geom/curve/curvature.hpp"
#include "eli/geom/curve/pseudo/explicit_bezier.hpp"

template<typename data__>
class explicit_bezier_curve_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::pseudo::explicit_bezier<data__> curve_type;
    typedef typename curve_type::control_point_type control_point_type;
    typedef typename curve_type::point_type point_type;
    typedef typename curve_type::data_type data_type;
    typedef typename curve_type::index_type index_type;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(explicit_bezier_curve_test_suite<float>::assignment_test);
      TEST_ADD(explicit_bezier_curve_test_suite<float>::reverse_test);
      TEST_ADD(explicit_bezier_curve_test_suite<float>::evaluation_test);
      TEST_ADD(explicit_bezier_curve_test_suite<float>::derivative_test);
      TEST_ADD(explicit_bezier_curve_test_suite<float>::promotion_test);
      TEST_ADD(explicit_bezier_curve_test_suite<float>::demotion_test);
      TEST_ADD(explicit_bezier_curve_test_suite<float>::length_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(explicit_bezier_curve_test_suite<double>::assignment_test);
      TEST_ADD(explicit_bezier_curve_test_suite<double>::reverse_test);
      TEST_ADD(explicit_bezier_curve_test_suite<double>::evaluation_test);
      TEST_ADD(explicit_bezier_curve_test_suite<double>::derivative_test);
      TEST_ADD(explicit_bezier_curve_test_suite<double>::promotion_test);
      TEST_ADD(explicit_bezier_curve_test_suite<double>::demotion_test);
      TEST_ADD(explicit_bezier_curve_test_suite<double>::length_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(explicit_bezier_curve_test_suite<long double>::assignment_test);
      TEST_ADD(explicit_bezier_curve_test_suite<long double>::reverse_test);
      TEST_ADD(explicit_bezier_curve_test_suite<long double>::evaluation_test);
      TEST_ADD(explicit_bezier_curve_test_suite<long double>::derivative_test);
      TEST_ADD(explicit_bezier_curve_test_suite<long double>::promotion_test);
      TEST_ADD(explicit_bezier_curve_test_suite<long double>::demotion_test);
      TEST_ADD(explicit_bezier_curve_test_suite<long double>::length_test);
    }

  public:
    explicit_bezier_curve_test_suite()
    {
      AddTests(data__());
    }
    ~explicit_bezier_curve_test_suite()
    {
    }

  private:
    void octave_print(int figno, const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts, const curve_type &bez) const
    {
      size_t i;

      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << "xpts=[" << pts[0].x();
      for (i=1; i<pts.size(); ++i)
        std::cout << ", " << pts[i].x();
      std::cout << "];" << std::endl;
      std::cout << "ypts=[" << pts[0].y();
      for (i=1; i<pts.size(); ++i)
        std::cout << ", " << pts[i].y();
      std::cout << "];" << std::endl;

      std::vector<data_type> t(101);
      for (i=0; i<t.size(); ++i)
        t[i]=static_cast<data_type>(i)/(t.size()-1);

      std::cout << "xint=[" << bez.f(t[0])(0);
      for (i=1; i<t.size(); ++i)
        std::cout << ", " << bez.f(t[i])(0);
      std::cout << "];" << std::endl;
      std::cout << "yint=[" << bez.f(t[0])(1);
      for (i=1; i<t.size(); ++i)
        std::cout << ", " << bez.f(t[i])(1);
      std::cout << "];" << std::endl;

      std::cout << "plot(xpts, ypts, 'bo', xint, yint, 'k-');" << std::endl;
    }

    void assignment_test()
    {
      curve_type ebc1(3), ebc2;
      control_point_type cntrl_in[4];

      // test default constructor then set control points
      cntrl_in[0] << 2.0;
      cntrl_in[1] << 1.0;
      cntrl_in[2] << 3.5;
      cntrl_in[3] << 4.0;

      for (index_type i=0; i<4; ++i)
      {
        ebc1.set_control_point(cntrl_in[i], i);
      }

      for (index_type i=0; i<4; ++i)
      {
        TEST_ASSERT(ebc1.get_control_point(i)==cntrl_in[i]);
      }

      // test copy ctr
      curve_type ebcc1(ebc1);
      for (index_type i=0; i<4; ++i)
      {
        TEST_ASSERT(ebcc1.get_control_point(i)==cntrl_in[i]);
      }

      // test assignment operator
      ebc2=ebc1;
      for (index_type i=0; i<4; ++i)
      {
        TEST_ASSERT(ebc2.get_control_point(i)==cntrl_in[i]);
      }

      // test order
      TEST_ASSERT(ebc2.degree()==3);
    }

    void reverse_test()
    {
      curve_type ebc1(3), ebc2;
      control_point_type cntrl_in[4];

      // test default constructor then set control points
      cntrl_in[0] << 2.0;
      cntrl_in[1] << 1.0;
      cntrl_in[2] << 3.5;
      cntrl_in[3] << 4.0;

      for (index_type i=0; i<4; ++i)
      {
        ebc1.set_control_point(cntrl_in[i], i);
      }

      ebc2=ebc1;
      ebc1.reverse();
      for (index_type i=0; i<4; ++i)
      {
        TEST_ASSERT(ebc1.get_control_point(i)==ebc2.get_control_point(3-i));
      }
    }

    void evaluation_test()
    {
      typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
      data_type eps(std::numeric_limits<data__>::epsilon());

      control_point_type cntrl_in[5];
      typename bezier_curve_type::control_point_type bez_cntrl[5];
      curve_type ebc;
      bezier_curve_type bc;
      typename curve_type::point_type eval_out, eval_ref;
      typename curve_type::data_type t;

      // set control points and create curves
      cntrl_in[0] << 2.0;
      cntrl_in[1] << 1.5;
      cntrl_in[2] << 0.0;
      cntrl_in[3] << 1.0;
      cntrl_in[4] << 0.5;
      bez_cntrl[0] << 0,    2;
      bez_cntrl[1] << 0.25, 1.5;
      bez_cntrl[2] << 0.5,  0;
      bez_cntrl[3] << 0.75, 1;
      bez_cntrl[4] << 1,    0.5;

      ebc.resize(4);
      bc.resize(4);
      for (index_type i=0; i<5; ++i)
      {
        ebc.set_control_point(cntrl_in[i], i);
        bc.set_control_point(bez_cntrl[i], i);
      }

      // test evaluation at end points
      t=0;
      eval_out=ebc.f(t);
      eval_ref=bc.f(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.f(t);
      eval_ref=bc.f(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test evaluation at interior point
      t=static_cast<data__>(0.45);
      eval_out=ebc.f(t);
      eval_ref=bc.f(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<1.1*eps);
    }

    void derivative_test()
    {
      typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
      data_type eps(std::numeric_limits<data__>::epsilon());
      control_point_type cntrl_in[5];
      typename bezier_curve_type::control_point_type bez_cntrl[5];
      curve_type ebc;
      bezier_curve_type bc;
      typename curve_type::point_type eval_out, eval_ref;
      typename curve_type::data_type t;

      // set control points and create curves
      cntrl_in[0] << 2.0;
      cntrl_in[1] << 1.5;
      cntrl_in[2] << 0.0;
      cntrl_in[3] << 1.0;
      cntrl_in[4] << 0.5;
      bez_cntrl[0] << 0,    2;
      bez_cntrl[1] << 0.25, 1.5;
      bez_cntrl[2] << 0.5,  0;
      bez_cntrl[3] << 0.75, 1;
      bez_cntrl[4] << 1,    0.5;

      ebc.resize(4);
      bc.resize(4);
      for (index_type i=0; i<5; ++i)
      {
        ebc.set_control_point(cntrl_in[i], i);
        bc.set_control_point(bez_cntrl[i], i);
      }

      // test 1st derivative at end points
      t=0;
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 1st derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      if (typeid(data__)==typeid(float))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
      }
      else
      {
        TEST_ASSERT(eval_out==eval_ref);
      }

      // test 2nd derivative at end points
      t=0;
      eval_out=ebc.fpp(t);
      eval_ref=bc.fpp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.fpp(t);
      eval_ref=bc.fpp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 2nd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=ebc.fpp(t);
      eval_ref=bc.fpp(t);
      if (typeid(data__)==typeid(float))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<17*eps);
      }
      else
      {
        TEST_ASSERT(eval_out==eval_ref);
      }

      // test 3rd derivative at end points
      t=0;
      eval_out=ebc.fppp(t);
      eval_ref=bc.fppp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.fppp(t);
      eval_ref=bc.fppp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 3rd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=ebc.fppp(t);
      eval_ref=bc.fppp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test curvature at end points
      data_type curv_out, curv_ref;
      t=0;
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      TEST_ASSERT(curv_out==curv_ref);
      t=1;
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      TEST_ASSERT(curv_out==curv_ref);

      // test curvature at interior point
      t=static_cast<data__>(0.45);
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      if (typeid(data__)==typeid(float))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
      }
      else
      {
        TEST_ASSERT(curv_out==curv_ref);
      }
    }

    void promotion_test()
    {
      typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
      data_type eps(std::numeric_limits<data__>::epsilon());
      control_point_type cntrl_in[5];
      typename bezier_curve_type::control_point_type bez_cntrl[5];
      curve_type ebc;
      bezier_curve_type bc;
      typename curve_type::point_type eval_out, eval_ref;
      typename curve_type::data_type t, curv_out, curv_ref;

      // set control points and create curves
      cntrl_in[0] << 2.0;
      cntrl_in[1] << 1.5;
      cntrl_in[2] << 0.0;
      cntrl_in[3] << 1.0;
      cntrl_in[4] << 0.5;
      bez_cntrl[0] << 0,    2;
      bez_cntrl[1] << 0.25, 1.5;
      bez_cntrl[2] << 0.5,  0;
      bez_cntrl[3] << 0.75, 1;
      bez_cntrl[4] << 1,    0.5;

      ebc.resize(4);
      bc.resize(4);
      for (index_type i=0; i<5; ++i)
      {
        ebc.set_control_point(cntrl_in[i], i);
        bc.set_control_point(bez_cntrl[i], i);
      }

      ebc.degree_promote();
      bc.degree_promote();

      // test to see if degree has increased
      TEST_ASSERT(ebc.degree()==bc.degree());

      // test evaluation at end points
      t=0;
      eval_out=ebc.f(t);
      eval_ref=bc.f(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.f(t);
      eval_ref=bc.f(t);
      if (typeid(data_type)==typeid(double))
      {
        TEST_ASSERT((eval_out-eval_ref).norm()<5*eps);
      }
      else
      {
        TEST_ASSERT(eval_out==eval_ref);
      }

      // test evaluation at interior point
      t=static_cast<data__>(0.45);
      eval_out=ebc.f(t);
      eval_ref=bc.f(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<2.1*eps);

      // test 1st derivative at end points
      t=0;
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<25*eps);

      // test 1st derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=ebc.fp(t);
      eval_ref=bc.fp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<4.1*eps);

      // test 2nd derivative at end points
      t=0;
      eval_out=ebc.fpp(t);
      eval_ref=bc.fpp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=ebc.fpp(t);
      eval_ref=bc.fpp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<25*eps);

      // test 2nd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=ebc.fpp(t);
      eval_ref=bc.fpp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<17*eps);

      // test 3rd derivative at end points
      t=0;
      eval_out=ebc.fppp(t);
      eval_ref=bc.fppp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<132*eps);
      t=1;
      eval_out=ebc.fppp(t);
      eval_ref=bc.fppp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<151*eps);

      // test 3rd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=ebc.fppp(t);
      eval_ref=bc.fppp(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<65*eps);

      // test curvature at end points
      t=0;
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      TEST_ASSERT(curv_out==curv_ref);
      t=1;
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      TEST_ASSERT(std::abs(curv_out-curv_ref)<48*eps);

      // test curvature at interior point
      t=static_cast<data__>(0.45);
      eli::geom::curve::curvature(curv_out, ebc, t);
      eli::geom::curve::curvature(curv_ref, bc, t);
      TEST_ASSERT(std::abs(curv_out-curv_ref)<10*eps);
    }

    void demotion_test()
    {
      // no constraint
      {
        typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
        control_point_type cntrl_in[9];
        typename bezier_curve_type::control_point_type bez_cntrl[9];
        curve_type ebc;
        bezier_curve_type bc;

        // set control points and create curves
        cntrl_in[0] <<  2.0;
        cntrl_in[1] <<  1.5;
        cntrl_in[2] <<  1.0;
        cntrl_in[3] <<  0.5;
        cntrl_in[4] <<  0.0;
        cntrl_in[5] << -0.5;
        cntrl_in[6] << -1.0;
        cntrl_in[7] <<  1.0;
        cntrl_in[8] <<  0.5;
        bez_cntrl[0] << static_cast<data__>(0),     2.0;
        bez_cntrl[1] << static_cast<data__>(0.125), 1.5;
        bez_cntrl[2] << static_cast<data__>(0.25),  1.0;
        bez_cntrl[3] << static_cast<data__>(0.375), 0.5;
        bez_cntrl[4] << static_cast<data__>(0.5),   0.0;
        bez_cntrl[5] << static_cast<data__>(0.625),-0.5;
        bez_cntrl[6] << static_cast<data__>(0.75), -1.0;
        bez_cntrl[7] << static_cast<data__>(0.875), 1.0;
        bez_cntrl[8] << static_cast<data__>(1),     0.5;

        ebc.resize(8);
        bc.resize(8);
        for (index_type i=0; i<9; ++i)
        {
          ebc.set_control_point(cntrl_in[i], i);
          bc.set_control_point(bez_cntrl[i], i);
        }

        ebc.degree_demote(eli::geom::general::NOT_CONNECTED);
        bc.degree_demote(eli::geom::general::NOT_CONNECTED);
        for(index_type i=0; i<8; ++i)
        {
          if (typeid(data__)==typeid(float))
          {
            TEST_ASSERT((ebc.get_control_point(i).col(0)-bc.get_control_point(i).col(1)).norm()<36*std::numeric_limits<data__>::epsilon());
          }
          else
          {
            TEST_ASSERT(ebc.get_control_point(i).col(0)==bc.get_control_point(i).col(1));
          }
        }
      }

      // C0 constraint
      {
        typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
        control_point_type cntrl_in[9];
        typename bezier_curve_type::control_point_type bez_cntrl[9];
        curve_type ebc;
        bezier_curve_type bc;

        // set control points and create curves
        cntrl_in[0] <<  2.0;
        cntrl_in[1] <<  1.5;
        cntrl_in[2] <<  1.0;
        cntrl_in[3] <<  0.5;
        cntrl_in[4] <<  0.0;
        cntrl_in[5] << -0.5;
        cntrl_in[6] << -1.0;
        cntrl_in[7] <<  1.0;
        cntrl_in[8] <<  0.5;
        bez_cntrl[0] << static_cast<data__>(0),     2.0;
        bez_cntrl[1] << static_cast<data__>(0.125), 1.5;
        bez_cntrl[2] << static_cast<data__>(0.25),  1.0;
        bez_cntrl[3] << static_cast<data__>(0.375), 0.5;
        bez_cntrl[4] << static_cast<data__>(0.5),   0.0;
        bez_cntrl[5] << static_cast<data__>(0.625),-0.5;
        bez_cntrl[6] << static_cast<data__>(0.75), -1.0;
        bez_cntrl[7] << static_cast<data__>(0.875), 1.0;
        bez_cntrl[8] << static_cast<data__>(1),     0.5;

        ebc.resize(8);
        bc.resize(8);
        for (index_type i=0; i<9; ++i)
        {
          ebc.set_control_point(cntrl_in[i], i);
          bc.set_control_point(bez_cntrl[i], i);
        }

        ebc.degree_demote(eli::geom::general::C0);
        bc.degree_demote(eli::geom::general::C0);
        for(index_type i=0; i<8; ++i)
        {
          TEST_ASSERT(ebc.get_control_point(i).col(0)==bc.get_control_point(i).col(1));
        }
      }

      // C1 constraint
      {
        typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
        control_point_type cntrl_in[9];
        typename bezier_curve_type::control_point_type bez_cntrl[9];
        curve_type ebc;
        bezier_curve_type bc;

        // set control points and create curves
        cntrl_in[0] <<  2.0;
        cntrl_in[1] <<  1.5;
        cntrl_in[2] <<  1.0;
        cntrl_in[3] <<  0.5;
        cntrl_in[4] <<  0.0;
        cntrl_in[5] << -0.5;
        cntrl_in[6] << -1.0;
        cntrl_in[7] <<  1.0;
        cntrl_in[8] <<  0.5;
        bez_cntrl[0] << static_cast<data__>(0),     2.0;
        bez_cntrl[1] << static_cast<data__>(0.125), 1.5;
        bez_cntrl[2] << static_cast<data__>(0.25),  1.0;
        bez_cntrl[3] << static_cast<data__>(0.375), 0.5;
        bez_cntrl[4] << static_cast<data__>(0.5),   0.0;
        bez_cntrl[5] << static_cast<data__>(0.625),-0.5;
        bez_cntrl[6] << static_cast<data__>(0.75), -1.0;
        bez_cntrl[7] << static_cast<data__>(0.875), 1.0;
        bez_cntrl[8] << static_cast<data__>(1),     0.5;

        ebc.resize(8);
        bc.resize(8);
        for (index_type i=0; i<9; ++i)
        {
          ebc.set_control_point(cntrl_in[i], i);
          bc.set_control_point(bez_cntrl[i], i);
        }

        ebc.degree_demote(eli::geom::general::C1);
        bc.degree_demote(eli::geom::general::C1);
        for(index_type i=0; i<8; ++i)
        {
          TEST_ASSERT(ebc.get_control_point(i).col(0)==bc.get_control_point(i).col(1));
        }
      }

      // C2 constraint
      {
        typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
        control_point_type cntrl_in[9];
        typename bezier_curve_type::control_point_type bez_cntrl[9];
        curve_type ebc;
        bezier_curve_type bc;

        // set control points and create curves
        cntrl_in[0] <<  2.0;
        cntrl_in[1] <<  1.5;
        cntrl_in[2] <<  1.0;
        cntrl_in[3] <<  0.5;
        cntrl_in[4] <<  0.0;
        cntrl_in[5] << -0.5;
        cntrl_in[6] << -1.0;
        cntrl_in[7] <<  1.0;
        cntrl_in[8] <<  0.5;
        bez_cntrl[0] << static_cast<data__>(0),     2.0;
        bez_cntrl[1] << static_cast<data__>(0.125), 1.5;
        bez_cntrl[2] << static_cast<data__>(0.25),  1.0;
        bez_cntrl[3] << static_cast<data__>(0.375), 0.5;
        bez_cntrl[4] << static_cast<data__>(0.5),   0.0;
        bez_cntrl[5] << static_cast<data__>(0.625),-0.5;
        bez_cntrl[6] << static_cast<data__>(0.75), -1.0;
        bez_cntrl[7] << static_cast<data__>(0.875), 1.0;
        bez_cntrl[8] << static_cast<data__>(1),     0.5;

        ebc.resize(8);
        bc.resize(8);
        for (index_type i=0; i<9; ++i)
        {
          ebc.set_control_point(cntrl_in[i], i);
          bc.set_control_point(bez_cntrl[i], i);
        }

        ebc.degree_demote(eli::geom::general::C2);
        bc.degree_demote(eli::geom::general::C2);
        for(index_type i=0; i<8; ++i)
        {
          TEST_ASSERT(ebc.get_control_point(i).col(0)==bc.get_control_point(i).col(1));
        }
      }
    }

    void length_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
      typedef eli::geom::curve::bezier<data__, 2> bezier_curve_type;
      control_point_type cntrl_in[5];
      typename bezier_curve_type::control_point_type bez_cntrl[5];
      curve_type ebc;
      bezier_curve_type bc;
      typename curve_type::point_type eval_out, eval_ref;
      typename curve_type::data_type length_cal, length_ref;

      // set control points and create curves
      cntrl_in[0] << 2.0;
      cntrl_in[1] << 1.5;
      cntrl_in[2] << 0.0;
      cntrl_in[3] << 1.0;
      cntrl_in[4] << 0.5;
      bez_cntrl[0] << 0,    2;
      bez_cntrl[1] << 0.25, 1.5;
      bez_cntrl[2] << 0.5,  0;
      bez_cntrl[3] << 0.75, 1;
      bez_cntrl[4] << 1,    0.5;

      ebc.resize(4);
      bc.resize(4);
      for (index_type i=0; i<5; ++i)
      {
        ebc.set_control_point(cntrl_in[i], i);
        bc.set_control_point(bez_cntrl[i], i);
      }

      // calculate the length of curve
      data_type tol(std::sqrt(eps));
      eli::geom::curve::length(length_cal, ebc, tol);
      eli::geom::curve::length(length_ref, bc, tol);
      TEST_ASSERT(length_cal==length_ref);

      // test computing some segment length
      typename curve_type::data_type t0, t1;
      t0 = static_cast<data__>(0.2);
      t1 = static_cast<data__>(0.7);

      eli::geom::curve::length(length_cal, ebc, t0, t1, tol);
      eli::geom::curve::length(length_ref, bc, t0, t1, tol);
      TEST_ASSERT(length_cal==length_ref);
    }
};
#endif

