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

#ifndef bezier_curve_test_suite_hpp
#define bezier_curve_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/length.hpp"
#include "eli/geom/curve/curvature.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class bezier_curve_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::bezier<data__, 3> bezier_type;
    typedef typename bezier_type::fit_container_type fit_container_type;
    typedef typename fit_container_type::constraint_point_type constraint_point_type;
    typedef typename bezier_type::control_point_type control_point_type;
    typedef typename bezier_type::point_type point_type;
    typedef typename bezier_type::data_type data_type;


  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(bezier_curve_test_suite<float>::assignment_test);
      TEST_ADD(bezier_curve_test_suite<float>::evaluation_test);
      TEST_ADD(bezier_curve_test_suite<float>::derivative_test);
      TEST_ADD(bezier_curve_test_suite<float>::promotion_test);
      TEST_ADD(bezier_curve_test_suite<float>::demotion_test);
      TEST_ADD(bezier_curve_test_suite<float>::split_test);
      TEST_ADD(bezier_curve_test_suite<float>::length_test);
      TEST_ADD(bezier_curve_test_suite<float>::fit_free_ends_test);
      TEST_ADD(bezier_curve_test_suite<float>::fit_C0_ends_test);
      TEST_ADD(bezier_curve_test_suite<float>::fit_C1_ends_test);
      TEST_ADD(bezier_curve_test_suite<float>::fit_C2_ends_test);
      TEST_ADD(bezier_curve_test_suite<float>::fit_closed_test);
      TEST_ADD(bezier_curve_test_suite<float>::interpolate_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(bezier_curve_test_suite<double>::assignment_test);
      TEST_ADD(bezier_curve_test_suite<double>::evaluation_test);
      TEST_ADD(bezier_curve_test_suite<double>::derivative_test);
      TEST_ADD(bezier_curve_test_suite<double>::promotion_test);
      TEST_ADD(bezier_curve_test_suite<double>::demotion_test);
      TEST_ADD(bezier_curve_test_suite<double>::split_test);
      TEST_ADD(bezier_curve_test_suite<double>::length_test);
      TEST_ADD(bezier_curve_test_suite<double>::fit_free_ends_test);
      TEST_ADD(bezier_curve_test_suite<double>::fit_C0_ends_test);
      TEST_ADD(bezier_curve_test_suite<double>::fit_C1_ends_test);
      TEST_ADD(bezier_curve_test_suite<double>::fit_C2_ends_test);
      TEST_ADD(bezier_curve_test_suite<double>::fit_closed_test);
      TEST_ADD(bezier_curve_test_suite<double>::interpolate_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(bezier_curve_test_suite<long double>::assignment_test);
      TEST_ADD(bezier_curve_test_suite<long double>::evaluation_test);
      TEST_ADD(bezier_curve_test_suite<long double>::derivative_test);
      TEST_ADD(bezier_curve_test_suite<long double>::promotion_test);
      TEST_ADD(bezier_curve_test_suite<long double>::demotion_test);
      TEST_ADD(bezier_curve_test_suite<long double>::split_test);
      TEST_ADD(bezier_curve_test_suite<long double>::length_test);
      TEST_ADD(bezier_curve_test_suite<long double>::fit_free_ends_test);
      TEST_ADD(bezier_curve_test_suite<long double>::fit_C0_ends_test);
      TEST_ADD(bezier_curve_test_suite<long double>::fit_C1_ends_test);
      TEST_ADD(bezier_curve_test_suite<long double>::fit_C2_ends_test);
      TEST_ADD(bezier_curve_test_suite<long double>::fit_closed_test);
      TEST_ADD(bezier_curve_test_suite<long double>::interpolate_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(bezier_curve_test_suite<dd_real>::assignment_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::evaluation_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::derivative_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::promotion_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::demotion_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::split_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::length_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::fit_free_ends_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::fit_C0_ends_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::fit_C1_ends_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::fit_C2_ends_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::fit_closed_test);
      TEST_ADD(bezier_curve_test_suite<dd_real>::interpolate_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(bezier_curve_test_suite<qd_real>::assignment_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::evaluation_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::derivative_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::promotion_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::demotion_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::split_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::length_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::fit_free_ends_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::fit_C0_ends_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::fit_C1_ends_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::fit_C2_ends_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::fit_closed_test);
      TEST_ADD(bezier_curve_test_suite<qd_real>::interpolate_test);
    }
#endif

  public:
    bezier_curve_test_suite()
    {
      AddTests(data__());
    }
    ~bezier_curve_test_suite()
    {
    }

  private:
    void octave_print(int figno, const std::vector<point_type> &pts, const bezier_type &bez) const
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

      std::vector<data__> t(101);
      for (i=0; i<t.size(); ++i)
        t[i]=static_cast<data__>(i)/(t.size()-1);

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

    void create_circle(std::vector<point_type> &pts)
    {
      // NOTE: This will not create a closed circle, the last point will be
      //       the point just prior to the 2*pi point
      size_t n=pts.size();
      for (size_t i=0; i<n; ++i)
      {
        data__ theta(eli::constants::math<data__>::two_pi()*static_cast<data__>(i)/n);
        pts[i](0)=std::cos(theta);
        pts[i](1)=std::sin(theta);
        pts[i](2)=0;
      }
    }

    void assignment_test()
    {
      bezier_type bc1, bc2;
      control_point_type cntrl_in(4,3), cntrl_out;

      // test default constructor then set control points
      cntrl_in << 2.0, 2.0, 0.0,
                  1.0, 1.5, 0.0,
                  3.5, 0.0, 0.0,
                  4.0, 1.0, 0.0;
      bc1.set_control_points(cntrl_in);
      bc1.get_control_points(cntrl_out);
      TEST_ASSERT(cntrl_in==cntrl_out);
      cntrl_out.setZero();

      // test constructor with vector of control points
      bezier_type bcc1(cntrl_in);
      bcc1.get_control_points(cntrl_out);
      TEST_ASSERT(cntrl_in==cntrl_out);
      cntrl_out.setZero();

      // test copy ctr
      bezier_type bcc2(bc1);
      bcc2.get_control_points(cntrl_out);
      TEST_ASSERT(cntrl_in==cntrl_out);
      cntrl_out.setZero();

      // test assignment operator
      bc2=bc1;
      bc2.get_control_points(cntrl_out);
      TEST_ASSERT(cntrl_in==cntrl_out);
      cntrl_out.setZero();

      // test order
      TEST_ASSERT(bc2.degree()==cntrl_in.rows()-1);
    }

    void evaluation_test()
    {
      control_point_type cntrl_in(4,3);
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // set control points
      cntrl_in << 2.0, 2.0, 0.0,
                  1.0, 1.5, 0.0,
                  3.5, 0.0, 0.0,
                  4.0, 1.0, 0.0;

      bezier_type bc1(cntrl_in), bc2;
      point_type  eval_out, eval_ref;
      data_type  t;

      // test evaluation at end points
      t=0;
      eval_out=bc1.f(t);
      TEST_ASSERT(eval_out==cntrl_in.row(0));
      t=1;
      eval_out=bc1.f(t);
      TEST_ASSERT(eval_out==cntrl_in.row(3));

      // test evaluation at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc1.f(t);
      eval_ref << static_cast<data__>(2.2750625), static_cast<data__>(1.0364375), static_cast<data__>(0);
      TEST_ASSERT((eval_out-eval_ref).norm()<5e3*eps);
//       if ( (typeid(data__)==typeid(dd_real)) || (typeid(data__)==typeid(qd_real)) )
//       {
//         std::cout << "257 rat=" << (eval_out-eval_ref).norm()/std::numeric_limits<double>::epsilon() << std::endl;
//       }

    }

    void derivative_test()
    {
      control_point_type cntrl_in(4,3);
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
      // set control points
      cntrl_in << 2.0, 2.0, 0.0,
                  1.0, 1.5, 0.0,
                  3.5, 0.0, 0.0,
                  4.0, 1.0, 0.0;

      bezier_type bc1(cntrl_in), bc2;
      point_type eval_out, eval_ref;
      data_type t;

      // test 1st derivative at end points
      t=0;
      eval_out=bc1.fp(t);
      eval_ref << -3, -1.5, 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc1.fp(t);
      eval_ref << 1.5, 3, 0;
      TEST_ASSERT(eval_out==eval_ref);

      // test 1st derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc1.fp(t);
      eval_ref << static_cast<data__>(3.10875), static_cast<data__>(-2.07375), static_cast<data__>(0);
      TEST_ASSERT((eval_out-eval_ref).norm() < 5e3*eps);
//       if ( (typeid(data__)==typeid(dd_real)) || (typeid(data__)==typeid(qd_real)) )
//       {
//         std::cout << "293 rat=" << (eval_out-eval_ref).norm()/std::numeric_limits<double>::epsilon() << std::endl;
//       }

      // test 2nd derivative at end points
      t=0;
      eval_out=bc1.fpp(t);
      eval_ref << 21, -6, 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc1.fpp(t);
      eval_ref << -12, 15, 0;
      TEST_ASSERT(eval_out==eval_ref);

      // test 2nd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc1.fpp(t);
      eval_ref << static_cast<data__>(6.15), static_cast<data__>(3.45), static_cast<data__>(0);
      TEST_ASSERT_DELTA((eval_out-eval_ref).norm(), 0, 1e4*eps);
//       if ( (typeid(data__)==typeid(dd_real)) || (typeid(data__)==typeid(qd_real)) )
//       {
//         std::cout << "313 rat=" << (eval_out-eval_ref).norm()/std::numeric_limits<double>::epsilon() << std::endl;
//       }

      // test 3rd derivative at end points
      t=0;
      eval_out=bc1.fppp(t);
      eval_ref << -33, 21, 0;
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc1.fppp(t);
      eval_ref << -33, 21, 0;
      TEST_ASSERT(eval_out==eval_ref);

      // test 3rd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc1.fppp(t);
      eval_ref << -33, 21, 0;
      TEST_ASSERT(eval_out==eval_ref);

      // test curvature at end points
      data_type curv_out, curv_ref;
      t=0;
      eli::geom::curve::curvature(curv_out, bc1, t);
      curv_ref=static_cast<data__>(1.31182654679988);
      TEST_ASSERT_DELTA(curv_out, curv_ref, std::max(static_cast<data__>(1e-13), static_cast<data__>(1e2)*eps));
      t=1;
      eli::geom::curve::curvature(curv_out, bc1, t);
      curv_ref=static_cast<data__>(1.55034046439985);
      TEST_ASSERT_DELTA(curv_out, curv_ref, std::max(static_cast<data__>(1e-13), static_cast<data__>(1e2)*eps));

      // test curvature at interior point
      t=static_cast<data__>(0.45);
      eli::geom::curve::curvature(curv_out, bc1, t);
      curv_ref=static_cast<data__>(0.449908807121445);
      TEST_ASSERT_DELTA(curv_out, curv_ref, std::max(static_cast<data__>(1e-13), static_cast<data__>(1e2)*eps));
    }

    void promotion_test()
    {
      control_point_type cntrl_in(5,3), cntrl_out;
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // set control points
      cntrl_in << 0.0, 0.0, 0.0,
                  0.0, 4.0, 0.0,
                  2.0, 4.0, 0.0,
                  2.0, 3.0, 0.0,
                  1.5, 3.0, 0.0;

      bezier_type bc1(cntrl_in), bc2(bc1);
      point_type eval_out, eval_ref;
      data_type t, curv_out, curv_ref;

      bc2.degree_promote();

      // test to see if degree has increased
      TEST_ASSERT(bc1.degree()+1==bc2.degree());

      // test to see if get expected polygon
      bc2.get_control_points(cntrl_out);
      if (cntrl_out.rows()!=cntrl_in.rows()+1)
      {
        TEST_ASSERT(false);
      }
      else
      {
        control_point_type cntrl_ref(6,3);

        cntrl_ref << 0.0,                      0.0,                      0.0,
                     0.0,                      static_cast<data__>(3.2), 0.0,
                     static_cast<data__>(1.2), 4.0,                      0.0,
                     2.0,                      static_cast<data__>(3.6), 0.0,
                     static_cast<data__>(1.9), 3.0,                      0.0,
                     static_cast<data__>(1.5), 3.0,                      0.0;
        if (typeid(data__)==typeid(long double))
        {
          TEST_ASSERT((cntrl_out-cntrl_ref).norm()<std::sqrt(eps));
        }
#ifdef ELI_QD_FOUND
        else if ( (typeid(data__)==typeid(dd_real)) || (typeid(data__)==typeid(qd_real)) )
        {
          TEST_ASSERT((cntrl_out-cntrl_ref).norm()<2.0*eps);
//         std::cout << "407 rat=" << (cntrl_out-cntrl_ref).norm()/std::numeric_limits<double>::epsilon() << std::endl;
        }
#endif
        else
        {
          TEST_ASSERT(cntrl_out==cntrl_ref);
        }
      }

      // test evaluation at end points
      t=0;
      eval_out=bc2.f(t);
      eval_ref=bc1.f(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc2.f(t);
      eval_ref=bc1.f(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test evaluation at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc2.f(t);
      eval_ref=bc1.f(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 1st derivative at end points
      t=0;
      eval_out=bc2.fp(t);
      eval_ref=bc1.fp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc2.fp(t);
      eval_ref=bc1.fp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 1st derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc2.fp(t);
      eval_ref=bc1.fp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 2nd derivative at end points
      t=0;
      eval_out=bc2.fpp(t);
      eval_ref=bc1.fpp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc2.fpp(t);
      eval_ref=bc1.fpp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 2nd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc2.fpp(t);
      eval_ref=bc1.fpp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 3rd derivative at end points
      t=0;
      eval_out=bc2.fppp(t);
      eval_ref=bc1.fppp(t);
      TEST_ASSERT(eval_out==eval_ref);
      t=1;
      eval_out=bc2.fppp(t);
      eval_ref=bc1.fppp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test 3rd derivative at interior point
      t=static_cast<data__>(0.45);
      eval_out=bc2.fppp(t);
      eval_ref=bc1.fppp(t);
      TEST_ASSERT(eval_out==eval_ref);

      // test curvature at end points
      t=0;
      eli::geom::curve::curvature(curv_out, bc2, t);
      eli::geom::curve::curvature(curv_ref, bc1, t);
      TEST_ASSERT(curv_out==curv_ref);
      t=1;
      eli::geom::curve::curvature(curv_out, bc2, t);
      eli::geom::curve::curvature(curv_ref, bc1, t);
      TEST_ASSERT(curv_out==curv_ref);

      // test curvature at interior point
      t=static_cast<data__>(0.45);
      eli::geom::curve::curvature(curv_out, bc2, t);
      eli::geom::curve::curvature(curv_ref, bc1, t);
      TEST_ASSERT(curv_out==curv_ref);
    }

    void demotion_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      {
        // hand checked with "Degree Reduction of Bezier Curves" by Dave Morgan
        control_point_type cntrl_in(7,3), cntrl_out, cntrl_ref(6,3);

        cntrl_in << 0.0, 0.0, 0.0,
                    2.0, 6.0, 0.0,
                    3.0, 0.0, 0.0,
                    5.0, 4.0, 0.0,
                    7.0, 1.0, 0.0,
                    5.0, 5.0, 0.0,
                   10.0, 6.0, 0.0;

        bezier_type bc1(cntrl_in), bc2(bc1);
        bc1.degree_demote(eli::geom::general::NOT_CONNECTED);

        // test if the degree is correct
        TEST_ASSERT(bc1.degree()+1==bc2.degree());

        // test of get the correct control points
        cntrl_ref << static_cast<data__>(-0.00878906), static_cast<data__>(0.0610352),  0.0,
                     static_cast<data__>( 2.5177734),  static_cast<data__>(6.3821289),  0.0,
                     static_cast<data__>( 2.8060547), static_cast<data__>(-0.16982422), 0.0,
                     static_cast<data__>( 8.0060547), static_cast<data__>( 2.5301758),  0.0,
                     static_cast<data__>( 4.1177734), static_cast<data__>( 3.9821289),  0.0,
                     static_cast<data__>( 9.9912109), static_cast<data__>( 6.0610352),  0.0;
        bc1.get_control_points(cntrl_out);
        if (cntrl_out.rows()!=cntrl_ref.rows())
        {
          TEST_ASSERT(false);
        }
        else
        {
          TEST_ASSERT_DELTA((cntrl_out-cntrl_ref).norm(), 0, 1e-6);
        }
      }

      // test whether can promote and demote and get the same control points
      {
        control_point_type cntrl_in(5,3), cntrl_out, cntrl_ref;

        // set control points
        cntrl_in << 0.0, 0.0, 0.0,
                    0.0, 4.0, 0.0,
                    2.0, 4.0, 0.0,
                    2.0, 3.0, 0.0,
                    1.5, 3.0, 0.0;

        bezier_type bc1(cntrl_in), bc2(bc1);
        bc2.degree_promote();
        bc2.degree_demote(eli::geom::general::NOT_CONNECTED);

        // test to see if degree has changed
        TEST_ASSERT(bc1.degree()==bc2.degree());

        // test to see if get expected polygon
        bc2.get_control_points(cntrl_out);
        cntrl_ref=cntrl_in;
        if (cntrl_out.rows()!=cntrl_ref.rows())
        {
          TEST_ASSERT(false);
        }
        else
        {
          if (typeid(data__)==typeid(long double))
          {
            TEST_ASSERT((cntrl_out-cntrl_ref).norm()<std::sqrt(eps));
          }
          else
          {
            TEST_ASSERT(cntrl_out==cntrl_ref);
          }
        }
      }

      // test whether can demote and maintain continuity at end points going all the way
      // to the n/2 limit and checking how close an interior point is
      {
        control_point_type cntrl_in(7,3);

        cntrl_in << 0.0, 0.0, 0.0,
                    2.0, 6.0, 0.0,
                    3.0, 0.0, 0.0,
                    5.0, 4.0, 0.0,
                    7.0, 1.0, 0.0,
                    5.0, 5.0, 0.0,
                   10.0, 6.0, 0.0;

        bezier_type bc1(cntrl_in);

        // No End Constraint
        {
          bezier_type bc2(bc1);
          point_type eval_out[3], eval_ref[3];
          data_type d[3], tv[3]={0, 1, static_cast<data__>(0.45)};
          bool success;

          success=bc2.degree_demote(eli::geom::general::NOT_CONNECTED);
          TEST_ASSERT(success);

          for (size_t k=0; k<3; ++k)
          {
            data_type t(tv[k]);
            eval_out[0]=bc2.f(t);
            eval_out[1]=bc2.fp(t);
            eval_out[2]=bc2.fpp(t);
            eval_ref[0]=bc1.f(t);
            eval_ref[1]=bc1.fp(t);
            eval_ref[2]=bc1.fpp(t);
            eli::geom::point::distance(d[0], eval_out[0], eval_ref[0]);
            eli::geom::point::distance(d[1], eval_out[1], eval_ref[1]);
            eli::geom::point::distance(d[2], eval_out[2], eval_ref[2]);

            if (k<2)
            {
              TEST_ASSERT_DELTA(d[0], 0.061665, 5e-6);
              TEST_ASSERT_DELTA(d[1], 4.43986, 5e-5);
              TEST_ASSERT_DELTA(d[2], 103.597, 5e-3);
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 0.05086, 1e-5);
              TEST_ASSERT_DELTA(d[1], 0.4205, 5e-4);
              TEST_ASSERT_DELTA(d[2], 7.4826, 5e-4);
            }
          }
        }

        // Point End Constraint
        {
          bezier_type bc2(bc1);
          point_type eval_out[3], eval_ref[3];
          data_type d[3], tv[3]={0, 1, static_cast<data__>(0.45)};
          bool success;

          success=bc2.degree_demote(eli::geom::general::C0);
          TEST_ASSERT(success);

          for (size_t k=0; k<3; ++k)
          {
            data_type t(tv[k]);
            eval_out[0]=bc2.f(t);
            eval_out[1]=bc2.fp(t);
            eval_out[2]=bc2.fpp(t);
            eval_ref[0]=bc1.f(t);
            eval_ref[1]=bc1.fp(t);
            eval_ref[2]=bc1.fpp(t);
            eli::geom::point::distance(d[0], eval_out[0], eval_ref[0]);
            eli::geom::point::distance(d[1], eval_out[1], eval_ref[1]);
            eli::geom::point::distance(d[2], eval_out[2], eval_ref[2]);

            if (k<2)
            {
              TEST_ASSERT(d[0] < 1e2*eps);
              TEST_ASSERT_DELTA(d[1], 3.82695, 5e-5);
              TEST_ASSERT_DELTA(d[2], 99.5007, 5e-4);
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 0.048737, 5e-6);
              TEST_ASSERT_DELTA(d[1], 0.43030, 5e-5);
              TEST_ASSERT_DELTA(d[2], 7.64895, 1e-4);
            }
          }
        }

        // Slope End Constraint
        {
          bezier_type bc2(bc1);
          point_type eval_out[3], eval_ref[3];
          data_type d[3], tv[3]={0, 1, static_cast<data__>(0.45)};
          bool success;

          success=bc2.degree_demote(eli::geom::general::C1);
          TEST_ASSERT(success);

          for (size_t k=0; k<3; ++k)
          {
            data_type t(tv[k]);
            eval_out[0]=bc2.f(t);
            eval_out[1]=bc2.fp(t);
            eval_out[2]=bc2.fpp(t);
            eval_ref[0]=bc1.f(t);
            eval_ref[1]=bc1.fp(t);
            eval_ref[2]=bc1.fpp(t);
            eli::geom::point::distance(d[0], eval_out[0], eval_ref[0]);
            eli::geom::point::distance(d[1], eval_out[1], eval_ref[1]);
            eli::geom::point::distance(d[2], eval_out[2], eval_ref[2]);

            if (k<2)
            {
              TEST_ASSERT(d[0] < 5e2*eps);
              TEST_ASSERT(d[1] < 1e3*eps);
              TEST_ASSERT_DELTA(d[2], 57.4042, 5e-4);
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 0.156476, 5e-6);
              TEST_ASSERT_DELTA(d[1], 0.900073, 5e-5);
              TEST_ASSERT_DELTA(d[2], 16.6997, 5e-4);
            }
          }
        }

        // Second Derivative Constraint
        {
          bezier_type bc2(bc1);
          point_type eval_out[3], eval_ref[3];
          data_type d[3], tv[3]={0, 1, static_cast<data__>(0.45)};
          bool success;

          success=bc2.degree_demote(eli::geom::general::C2);
          TEST_ASSERT(success);

          for (size_t k=0; k<3; ++k)
          {
            data_type t(tv[k]);
            eval_out[0]=bc2.f(t);
            eval_out[1]=bc2.fp(t);
            eval_out[2]=bc2.fpp(t);
            eval_ref[0]=bc1.f(t);
            eval_ref[1]=bc1.fp(t);
            eval_ref[2]=bc1.fpp(t);
            eli::geom::point::distance(d[0], eval_out[0], eval_ref[0]);
            eli::geom::point::distance(d[1], eval_out[1], eval_ref[1]);
            eli::geom::point::distance(d[2], eval_out[2], eval_ref[2]);

            if (k<2)
            {
              TEST_ASSERT(d[0] < 5e2*eps);
              TEST_ASSERT(d[1] < 5e3*eps);
              TEST_ASSERT(d[2] < 5e3*eps);
            }
            else
            {
              TEST_ASSERT_DELTA(d[0], 1.91466, 5e-5);
              TEST_ASSERT_DELTA(d[1], 2.32081, 5e-5);
              TEST_ASSERT_DELTA(d[2], 44.5408, 5e-4);
            }
          }
        }

      }
    }

    void split_test()
    {
      control_point_type cntrl_in(4,3), cntrl_out, cntrl_ref(4,3);
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // set control points
      cntrl_in << 0.0, 0.0, 0.0,
                  0.0, 2.0, 0.0,
                  8.0, 2.0, 0.0,
                  4.0, 0.0, 0.0;

      bezier_type bc1(cntrl_in), bc1l, bc1r;
      point_type eval_out, eval_ref;
      data_type t;

      // test split with known control points
      t=0.5;
      bc1.split(bc1l, bc1r, t);
      cntrl_ref << 0.0, 0.0, 0.0,
                   0.0, 1.0, 0.0,
                   2.0, 1.5, 0.0,
                   3.5, 1.5, 0.0;
      bc1l.get_control_points(cntrl_out);
      TEST_ASSERT(cntrl_out==cntrl_ref);
      cntrl_ref << 3.5, 1.5, 0.0,
                   5.0, 1.5, 0.0,
                   6.0, 1.0, 0.0,
                   4.0, 0.0, 0.0;
      bc1r.get_control_points(cntrl_out);
      TEST_ASSERT(cntrl_out==cntrl_ref);

      // split the curve and check the evaluations
      data_type tl, tr, ts;
      tl=static_cast<data__>(0.3);
      tr=static_cast<data__>(0.87);
      ts=static_cast<data__>(0.586);

      bc1.split(bc1l, bc1r, ts);

      // check the left curve
      t=tl*ts;
      eval_out=bc1l.f(tl);
      eval_ref=bc1.f(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1l.fp(tl);
      eval_ref=bc1.fp(t)*ts;
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1l.fpp(tl);
      eval_ref=bc1.fpp(t)*ts*ts;
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1l.fppp(tl);
      eval_ref=bc1.fppp(t)*ts*ts*ts;
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);

      // check the right curve
      t=ts+tr*(1-ts);
      eval_out=bc1r.f(tr);
      eval_ref=bc1.f(t);
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1r.fp(tr);
      eval_ref=bc1.fp(t)*(1-ts);
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1r.fpp(tr);
      eval_ref=bc1.fpp(t)*(1-ts)*(1-ts);
      TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      eval_out=bc1r.fppp(tr);
      eval_ref=bc1.fppp(t)*(1-ts)*(1-ts)*(1-ts);
      TEST_ASSERT((eval_out-eval_ref).norm()<1.21e2*eps);
    }

    void length_test()
    {
      control_point_type  cntrl_in(4,3);
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // set control points
      cntrl_in << 0.0, 0.0, 0.0,
                  0.0, 2.0, 0.0,
                  8.0, 2.0, 0.0,
                  4.0, 0.0, 0.0;

      bezier_type bc1(cntrl_in);
      data_type length_cal, length_ref;

      // calculate the reference length using simpson's rule
      size_t i, npts(1001), n(npts-1);
      std::vector<data__> speed(npts);
      for (i=0; i<npts; ++i)
      {
        data_type t(i/static_cast<data__>(n));
        speed[i]=bc1.fp(t).norm();
      }
      length_ref=speed[0];
      for (i=1; i<=n-1; i+=2)
      {
        length_ref+=4*speed[i];
        length_ref+=2*speed[i+1];
      }
      length_ref-=speed[n];
      length_ref*=static_cast<data__>(1-0)/n/3;

      // compute length and compare
      data_type tol(std::sqrt(eps));
      length(length_cal, bc1, tol);
      TEST_ASSERT_DELTA(1, length_cal/length_ref, tol);

      // test computing some segment length
      data_type length01_cal, length12_cal, t0, t1, t2;
      t0=0;
      t1=static_cast<data__>(0.3);
      t2=1;

      length(length01_cal, bc1, t0, t1, tol);
      length(length12_cal, bc1, t1, t2, tol);
      TEST_ASSERT_DELTA(1, (length01_cal+length12_cal)/length_cal, tol);
    }

    void fit_free_ends_test()
    {
      // open curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::NOT_CONNECTED);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.110);

        // check if went through points
        TEST_ASSERT(pts[0]!=bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT(bez.f(0)!=bez.f(1));

//         octave_print(1, pts, bez);
      }

      // closed curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.765);

        // check if went through points
        TEST_ASSERT(pts[0]!=bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT(bez.f(0)!=bez.f(1));

//         octave_print(2, pts, bez);
      }
    }

    void fit_C0_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // open curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::NOT_CONNECTED);
        fcon.add_start_C0_constraint();
        fcon.add_end_C0_constraint();

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.105);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT((pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()<70*eps);
        TEST_ASSERT(bez.f(0)!=bez.f(1));
//         if ( (typeid(data__)==typeid(float)) || (typeid(data__)==typeid(long double)) )
//         {
//           std::cout << "934 rat=" << (pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()/eps << std::endl;
//         }

//         octave_print(3, pts, bez);
      }

      // closed curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);
        fcon.add_start_C0_constraint();

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.207);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<172*eps);

//         octave_print(4, pts, bez);
      }

      // closed curve with additional point constraint
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);
        fcon.add_start_C0_constraint();
        fcon.add_C0_constraint(5);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.233);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<197*eps);
        TEST_ASSERT((bez.f(t[5])-pts[5]).norm()<14*eps);

//         octave_print(5, pts, bez);
      }
    }

    void fit_C1_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // open curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        constraint_point_type fp1, fp2;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fp2(0)=-eli::constants::math<data_type>::two_pi()*pts[pts.size()-1](1);
        fp2(1)=eli::constants::math<data_type>::two_pi()*pts[pts.size()-1](0);
        fp2(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::NOT_CONNECTED);
        fcon.add_start_C1_constraint(fp1);
        fcon.add_end_C1_constraint(fp2);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.366);

        // check if went through points
#ifdef ELI_QD_FOUND
        if (typeid(data__)==typeid(qd_real))
        {
          TEST_ASSERT((pts[0]-bez.f(t[0])).norm()<2.0*eps);
//         std::cout << "1132 rat=" << (pts[0]-bez.f(t[0])).norm()/eps << std::endl;
        }
        else
#endif
        {
          TEST_ASSERT(pts[0]==bez.f(t[0]));
        }
        TEST_ASSERT(fp1==bez.fp(t[0]));
        TEST_ASSERT((pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()<88*eps);
        TEST_ASSERT((fp2-bez.fp(t[t.size()-1])).norm()<351*eps);
        TEST_ASSERT(bez.f(0)!=bez.f(1));
//         if (typeid(data__)==typeid(float))
//         {
//           std::cout << "1071 rat=" << (pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()/eps << std::endl;
//           std::cout << "1072 rat=" << (fp2-bez.fp(t[t.size()-1])).norm()/eps << std::endl;
//         }

//         octave_print(6, pts, bez);
      }

      // closed curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        constraint_point_type fp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);
        fcon.add_start_C1_constraint(fp1);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.331);

        // check if went through points
#ifdef ELI_QD_FOUND
        if (typeid(data__)==typeid(qd_real))
        {
          TEST_ASSERT((pts[0]-bez.f(t[0])).norm()<2.0*eps);
//         std::cout << "1182 rat=" << (pts[0]-bez.f(t[0])).norm()/eps << std::endl;
        }
        else
#endif
        {
          TEST_ASSERT(pts[0]==bez.f(t[0]));
        }
        TEST_ASSERT(fp1==bez.fp(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<66*eps);

//         octave_print(7, pts, bez);
      }

      // closed curve with additional point constraint
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        constraint_point_type fp1, fp2;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fp2(0)=-eli::constants::math<data_type>::two_pi()*pts[5](1);
        fp2(1)=eli::constants::math<data_type>::two_pi()*pts[5](0);
        fp2(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);
        fcon.add_start_C1_constraint(fp1);
        fcon.add_C1_constraint(5, fp2);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.418);

        // check if went through points
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<5*eps);
        TEST_ASSERT(fp1==bez.fp(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<118*eps);
        TEST_ASSERT((bez.f(t[5])-pts[5]).norm()<19*eps);
        TEST_ASSERT((bez.fp(t[5])-fp2).norm()<73*eps);
//         if (typeid(data__)==typeid(long double))
//         {
//           std::cout << "1260 rat=" << (bez.f(t[5])-pts[5]).norm()/eps << std::endl;
//         }

//         octave_print(8, pts, bez);
      }
    }

    void fit_C2_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // open curve
      {
        fit_container_type fcon;
        size_t i, deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        constraint_point_type fp1, fpp1, fp2, fpp2;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fpp1(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](0);
        fpp1(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](1);
        fpp1(2)=0;
        fp2(0)=-eli::constants::math<data_type>::two_pi()*pts[pts.size()-1](1);
        fp2(1)=eli::constants::math<data_type>::two_pi()*pts[pts.size()-1](0);
        fp2(2)=0;
        fpp2(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[pts.size()-1](0);
        fpp2(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[pts.size()-1](1);
        fpp2(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::NOT_CONNECTED);
        fcon.add_start_C2_constraint(fp1, fpp1);
        fcon.add_end_C2_constraint(fp2, fpp2);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.372);

        // check if went through points
        TEST_ASSERT((pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()<422*eps);
        TEST_ASSERT((fp2-bez.fp(t[t.size()-1])).norm()<1.88e3*eps);
        TEST_ASSERT((fpp2-bez.fpp(t[t.size()-1])).norm()<7.31e3*eps);
        TEST_ASSERT(bez.f(0)!=bez.f(1));
        TEST_ASSERT((pts[0]-bez.f(t[0])).norm()<4*eps);
        TEST_ASSERT((fp1-bez.fp(t[0])).norm()<9*eps);
        TEST_ASSERT(fpp1==bez.fpp(t[0]));
//         if (typeid(data__)==typeid(long double))
//         {
//           std::cout << "1223 rat=" << (pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()/eps << std::endl;
//           std::cout << "1224 rat=" << (fp2-bez.fp(t[t.size()-1])).norm()/eps << std::endl;
//           std::cout << "1225 rat=" << (fpp2-bez.fpp(t[t.size()-1])).norm()/eps << std::endl;
//         }

//         octave_print(9, pts, bez);
      }

      // closed curve
      {
        fit_container_type fcon;
        size_t i, deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        constraint_point_type fp1, fpp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fpp1(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](0);
        fpp1(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](1);
        fpp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);
        fcon.add_start_C2_constraint(fp1, fpp1);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0479);

        // check if went through points
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<387*eps);
        TEST_ASSERT((pts[0]-bez.f(t[0])).norm()<3*eps);
        TEST_ASSERT((fp1-bez.fp(t[0])).norm()<17*eps);
        TEST_ASSERT((fpp1-bez.fpp(t[0])).norm()<65*eps);
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));

//         octave_print(10, pts, bez);
      }

      // closed curve with additional point constraint
      {
        fit_container_type fcon;
        size_t i, deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        constraint_point_type fp1, fpp1, fp2, fpp2;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fpp1(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](0);
        fpp1(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](1);
        fpp1(2)=0;
        fp2(0)=-eli::constants::math<data_type>::two_pi()*pts[10](1);
        fp2(1)=eli::constants::math<data_type>::two_pi()*pts[10](0);
        fp2(2)=0;
        fpp2(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[10](0);
        fpp2(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[10](1);
        fpp2(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);
        fcon.add_start_C2_constraint(fp1, fpp1);
        fcon.add_C2_constraint(10, fp2, fpp2);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0785);

        // check if went through points
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<420*eps);
        TEST_ASSERT((bez.f(t[10])-pts[10]).norm()<39*eps);
        TEST_ASSERT((bez.fp(t[10])-fp2).norm()<165*eps);
        TEST_ASSERT((bez.fpp(t[10])-fpp2).norm()<1.44e3*eps);
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<8*eps);
        TEST_ASSERT((bez.fp(t[0])-fp1).norm()<9*eps);
        TEST_ASSERT((bez.fpp(t[0])-fpp1).norm()<33*eps);
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
//         if ((typeid(data_type)==typeid(long double))||(typeid(data_type)==typeid(float)))
//         {
//           std::cout << "1441 rat=" << (bez.fp(t[10])-fp2).norm()/eps << std::endl;
//         }

//         octave_print(11, pts, bez);
      }
    }

    void fit_closed_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // fit no constraints C0 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.393);

        // check if went through points
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<380*eps);
        TEST_ASSERT(bez.fp(0)!=bez.fp(1));
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));
//         if (typeid(data__)==typeid(double))
//         {
//           std::cout << "1484 rat=" << (bez.f(1)-bez.f(0)).norm()/eps << std::endl;
//         }

//         octave_print(12, pts, bez);
      }

      // fit no constraints C1 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C1);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0312);

        // check if went through points
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<531*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<2.61e3*eps);
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));
//         if (typeid(data__)==typeid(double))
//         {
//           std::cout << "1518 rat=" << (bez.f(1)-bez.f(0)).norm()/eps << std::endl;
//         }

//         octave_print(13, pts, bez);
      }

      // fit no constraints C2 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C2);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0343);

        // check if went through points
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<398*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<1.87e3*eps);
        TEST_ASSERT((bez.fpp(1)-bez.fpp(0)).norm()<1.06e4*eps);
//         if (typeid(data__)==typeid(long double))
//         {
//           std::cout << "1443 rat=" << (bez.fp(1)-bez.fp(0)).norm()/eps << std::endl;
//         }
//         if (typeid(data__)==typeid(float))
//         {
//           std::cout << "err=" << err << std::endl;
//           std::cout << "1451 rat=" << (bez.fpp(1)-bez.fpp(0)).norm()/eps << std::endl;
//         }

//         octave_print(14, pts, bez);
      }

      // fit start constraints C0 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);
        fcon.add_start_C0_constraint();

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0185);

        // check if went through points
        TEST_ASSERT(bez.f(0)==pts[0]);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<431*eps);
        TEST_ASSERT(bez.fp(0)!=bez.fp(1));
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

//         octave_print(15, pts, bez);
      }

      // fit start constraints C1 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        constraint_point_type fp0;
        fp0(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp0(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp0(2)=eli::constants::math<data_type>::two_pi()*pts[0](2);
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C1);
        fcon.add_start_C1_constraint(fp0);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0395);

        // check if went through points
        TEST_ASSERT((bez.f(0)-pts[0]).norm()<2*eps);
        TEST_ASSERT(bez.fp(0)==fp0);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<483*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<2.77e3*eps);
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

//         octave_print(16, pts, bez);
      }

      // fit start constraints C2 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        constraint_point_type fp0, fpp0;
        fp0(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp0(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp0(2)=eli::constants::math<data_type>::two_pi()*pts[0](2);
        fpp0(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](0);
        fpp0(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](1);
        fpp0(2)=4*eli::constants::math<data_type>::pi_squared()*pts[0](2);
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C2);
        fcon.add_start_C2_constraint(fp0, fpp0);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.141);

        // check if went through points
        TEST_ASSERT((bez.f(0)-pts[0]).norm()<2*eps);
        TEST_ASSERT((bez.fp(0)-fp0).norm()<9*eps);
        TEST_ASSERT(bez.fpp(0)==fpp0);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<470*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<4.70e3*eps);
        TEST_ASSERT((bez.fpp(1)-bez.fpp(0)).norm()<1.00e4*eps);
//           std::cout << "1669 rat=" << (bez.f(1)-bez.f(0)).norm()/eps << std::endl;
//           std::cout << "1670 rat=" << (bez.fpp(1)-bez.fpp(0)).norm()/eps << std::endl;

//         octave_print(17, pts, bez);
      }

      // fit middle constraints C0 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        int ci=5;
        data_type t5;
        point_type f5;
        f5=pts[ci];
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);
        fcon.add_C0_constraint(ci);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0180);

        // check if went through points
        t5=t[ci];
        TEST_ASSERT((bez.f(t5)-f5).norm()<6*eps);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<594*eps);
        TEST_ASSERT(bez.fp(0)!=bez.fp(1));
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

//         octave_print(18, pts, bez);
      }

      // fit middle constraints C1 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        int ci=5;
        data_type t5;
        point_type f5;
        f5=pts[ci];
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C1);
        fcon.add_C0_constraint(ci);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0316);

        // check if went through points
        t5=t[ci];
        TEST_ASSERT((bez.f(t5)-f5).norm()<4*eps);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<391*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<2.48e3*eps);
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));
//         if (typeid(data__)==typeid(float))
//         {
//           std::cout << "1747 rat=" << (bez.f(1)-bez.f(0)).norm()/eps << std::endl;
//         }

//         octave_print(19, pts, bez);
      }

      // fit middle constraints C2 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        data_type err;
        create_circle(pts);

        // configure fit container
        int ci=5;
        data_type t5;
        point_type f5;
        f5=pts[ci];
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C2);
        fcon.add_C0_constraint(ci);

        // fit points
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0353);

        // check if went through points
        t5=t[ci];
        TEST_ASSERT((bez.f(t5)-f5).norm()<4*eps);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<422*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<2.37e3*eps);
        TEST_ASSERT((bez.fpp(1)-bez.fpp(0)).norm()<1.17e4*eps);
//         if ((typeid(data__)==typeid(float))||(typeid(data__)==typeid(double)))
//         {
//           std::cout << "1788 rat=" << (bez.f(1)-bez.f(0)).norm()/eps << std::endl;
//           std::cout << "1789 rat=" << (bez.fp(1)-bez.fp(0)).norm()/eps << std::endl;
//           std::cout << "1790 rat=" << (bez.fpp(1)-bez.fpp(0)).norm()/eps << std::endl;
//         }
//         if (typeid(data__)==typeid(long double))
//         {
//           std::cout << "1794 rat=" << (bez.fp(1)-bez.fp(0)).norm()/eps << std::endl;
//         }

//         octave_print(20, pts, bez);
      }
    }

    void interpolate_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // interpolate through all points open
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size());
        TEST_ASSERT(bez.open());

        // check if went through points
        TEST_ASSERT(bez.f(t[0])==pts[0]);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<4*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<7*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<13*eps);
//           std::cout << "1823 rat=" << (bez.f(t[1])-pts[1]).norm()/eps << std::endl;

//         octave_print(18, pts, bez);
      }

      // interpolate through all points closed
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(eli::geom::general::C0);

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size()+1);
        TEST_ASSERT(bez.closed());

        // check if went through points
        TEST_ASSERT(bez.f(t[0])==pts[0]);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<3*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<9*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<37*eps);
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<49*eps);

//         octave_print(19, pts, bez);
      }

      // interpolate through all points and specify C1 open
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[1](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[1](0);
        fp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.add_C1_constraint(1, fp1);

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size()+1);
        TEST_ASSERT(bez.open());

        // check if went through points
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<4*eps);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<10*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<28*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<63*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<17*eps);
//         if (typeid(data__)==typeid(float))
//         {
//           std::cout << "1890 rat=" << (bez.f(t[0])-pts[0]).norm()/eps << std::endl;
//           std::cout << "1891 rat=" << (bez.f(t[1])-pts[1]).norm()/eps << std::endl;
//           std::cout << "1892 rat=" << (bez.f(t[2])-pts[2]).norm()/eps << std::endl;
//         }
//         if (typeid(data__)==typeid(long double))
//         {
//           std::cout << "1896 rat=" << (bez.f(t[0])-pts[0]).norm()/eps << std::endl;
//         }

//         octave_print(20, pts, bez);
      }

      // interpolate through all points and specify C1 closed
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[1](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[1](0);
        fp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.add_C1_constraint(1, fp1);
        fcon.set_end_flag(eli::geom::general::C1);

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size()+3);
        TEST_ASSERT(bez.closed());

        // check if went through points
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<5*eps);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<3*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<12*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<89*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<22*eps);
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<202*eps);
        TEST_ASSERT((bez.fp(0)-bez.fp(1)).norm()<995*eps);
//           std::cout << "1926 rat=" << (bez.f(t[3])-pts[3]).norm()/eps << std::endl;
//           std::cout << "1927 rat=" << (bez.fp(0)-bez.fp(1)).norm()/eps << std::endl;

//         if (typeid(data__)==typeid(long double))
//         {
//           std::cout << "1797 rat=" << (bez.fp(0)-bez.fp(1)).norm()/eps << std::endl;
//         }

//         octave_print(21, pts, bez);
      }

      // interpolate through all points and specify C2 open
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1, fpp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[1](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[1](0);
        fp1(2)=0;
        fpp1(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[1](0);
        fpp1(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[1](1);
        fpp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.add_C2_constraint(1, fp1, fpp1);

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size()+2);
        TEST_ASSERT(bez.open());

        // check if went through points
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<9*eps);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<9*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<84*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<162*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<68*eps);
        TEST_ASSERT((bez.fpp(t[1])-fpp1).norm()<353*eps);
//           std::cout << "1969 rat=" << (bez.fpp(t[1])-fpp1).norm()/eps << std::endl;

//         octave_print(22, pts, bez);
      }

      // interpolate through all points and specify C2 closed
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1, fpp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[1](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[1](0);
        fp1(2)=0;
        fpp1(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[1](0);
        fpp1(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[1](1);
        fpp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.add_C2_constraint(1, fp1, fpp1);
        fcon.set_end_flag(eli::geom::general::C2);

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size()+5);
        TEST_ASSERT(bez.closed());

        // check if went through points
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<4*eps);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<4*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<32*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<143*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<39*eps);
        TEST_ASSERT((bez.fpp(t[1])-fpp1).norm()<391*eps);
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<605*eps);
        TEST_ASSERT((bez.fp(0)-bez.fp(1)).norm()<3.76e3*eps);
        TEST_ASSERT((bez.fpp(0)-bez.fpp(1)).norm()<16.7e3*eps);
//           std::cout << "2016 rat=" << (bez.f(t[3])-pts[3]).norm()/eps << std::endl;
//           std::cout << "2017 rat=" << (bez.f(0)-bez.f(1)).norm()/eps << std::endl;
//           std::cout << "2018 rat=" << (bez.fpp(0)-bez.fpp(1)).norm()/eps << std::endl;

//         if (typeid(data__)==typeid(long double))
//         {
//           std::cout << "1877 rat=" << (bez.f(t[2])-pts[2]).norm()/eps << std::endl;
//           std::cout << "1878 rat=" << (bez.f(t[3])-pts[3]).norm()/eps << std::endl;
//           std::cout << "1879 rat=" << (bez.fpp(t[1])-fpp1).norm()/eps << std::endl;
//           std::cout << "1880 rat=" << (bez.f(0)-bez.f(1)).norm()/eps << std::endl;
//           std::cout << "1881 rat=" << (bez.fp(0)-bez.fp(1)).norm()/eps << std::endl;
//           std::cout << "1882 rat=" << (bez.fpp(0)-bez.fpp(1)).norm()/eps << std::endl;
//         }

//         octave_print(23, pts, bez);
      }
    }
};

#endif

