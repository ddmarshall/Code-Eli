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

#ifndef bezier_curve_fit_test_suite_hpp
#define bezier_curve_fit_test_suite_hpp

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
class bezier_curve_fit_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::bezier<data__, 3> bezier_type;
    typedef typename bezier_type::fit_container_type fit_container_type;
    typedef typename fit_container_type::constraint_point_type constraint_point_type;
    typedef typename bezier_type::index_type index_type;
    typedef typename bezier_type::point_type point_type;
    typedef typename bezier_type::data_type data_type;
    typedef typename bezier_type::tolerance_type tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(bezier_curve_fit_test_suite<float>::fit_free_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<float>::fit_C0_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<float>::fit_C1_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<float>::fit_C2_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<float>::fit_closed_test);
      TEST_ADD(bezier_curve_fit_test_suite<float>::interpolate_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(bezier_curve_fit_test_suite<double>::fit_free_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<double>::fit_C0_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<double>::fit_C1_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<double>::fit_C2_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<double>::fit_closed_test);
      TEST_ADD(bezier_curve_fit_test_suite<double>::interpolate_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(bezier_curve_fit_test_suite<long double>::fit_free_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<long double>::fit_C0_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<long double>::fit_C1_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<long double>::fit_C2_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<long double>::fit_closed_test);
      TEST_ADD(bezier_curve_fit_test_suite<long double>::interpolate_test);
    }
#ifdef ELI_USING_QD
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(bezier_curve_fit_test_suite<dd_real>::fit_free_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<dd_real>::fit_C0_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<dd_real>::fit_C1_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<dd_real>::fit_C2_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<dd_real>::fit_closed_test);
      TEST_ADD(bezier_curve_fit_test_suite<dd_real>::interpolate_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(bezier_curve_fit_test_suite<qd_real>::fit_free_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<qd_real>::fit_C0_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<qd_real>::fit_C1_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<qd_real>::fit_C2_ends_test);
      TEST_ADD(bezier_curve_fit_test_suite<qd_real>::fit_closed_test);
      TEST_ADD(bezier_curve_fit_test_suite<qd_real>::interpolate_test);
    }
#endif

  public:
    bezier_curve_fit_test_suite()
    {
      AddTests(data__());
    }
    ~bezier_curve_fit_test_suite()
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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.765);

        // check if went through points
        TEST_ASSERT(pts[0]!=bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT(bez.f(0)==bez.f(1));

//         octave_print(2, pts, bez);
      }
    }

    void fit_C0_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_USING_QD
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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.105);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT((pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()<81*eps);
        TEST_ASSERT(bez.f(0)!=bez.f(1));
//         if (typeid(data__)==typeid(double))
//         {
//           std::cout << "998 rat=" << (pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()/eps << std::endl;
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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.207);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT(bez.f(0)==bez.f(1));

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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.233);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT(bez.f(0)==bez.f(1));
        TEST_ASSERT((bez.f(t[5])-pts[5]).norm()<15*eps);

//         octave_print(5, pts, bez);
      }
    }

    void fit_C1_ends_test()
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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.366);

        // check if went through points
        TEST_ASSERT(tol.approximately_equal(pts[0], bez.f(t[0])));
        TEST_ASSERT(tol.approximately_equal(fp1, bez.fp(t[0])));
        TEST_ASSERT(tol.approximately_equal(pts[pts.size()-1], bez.f(t[t.size()-1])));
        TEST_ASSERT(tol.approximately_equal(fp2, bez.fp(t[t.size()-1])));
        TEST_ASSERT(bez.f(0)!=bez.f(1));

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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.331);

        // check if went through points
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT(tol.approximately_equal(pts[0], bez.f(t[0])));
        TEST_ASSERT(tol.approximately_equal(fp1, bez.fp(t[0])));
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));

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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.418);

        // check if went through points
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT(tol.approximately_equal(pts[0], bez.f(t[0])));
        // precision in float calculation is totally lost
        if (typeid(data_type)!=typeid(float))
        {
          TEST_ASSERT(tol.approximately_equal(fp1, bez.fp(t[0])));
        }
        TEST_ASSERT(tol.approximately_equal(pts[5], bez.f(t[5])));
        TEST_ASSERT(tol.approximately_equal(fp2, bez.fp(t[5])));
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));

//         octave_print(8, pts, bez);
      }
    }

    void fit_C2_ends_test()
    {
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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.372);

        // check if went through points
        TEST_ASSERT(tol.approximately_equal(pts[0], bez.f(t[0])));
        TEST_ASSERT(tol.approximately_equal(fp1, bez.fp(t[0])));
        TEST_ASSERT(tol.approximately_equal(fpp1, bez.fpp(t[0])));
        TEST_ASSERT(tol.approximately_equal(pts[pts.size()-1], bez.f(t[t.size()-1])));
        TEST_ASSERT(tol.approximately_equal(fp2, bez.fp(t[t.size()-1])));
        TEST_ASSERT(tol.approximately_equal(fpp2, bez.fpp(t[t.size()-1])));
        TEST_ASSERT(bez.f(0)!=bez.f(1));

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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0479);

        // check if went through points
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
        TEST_ASSERT(tol.approximately_equal(pts[0], bez.f(t[0])));
        TEST_ASSERT(tol.approximately_equal(fp1, bez.fp(t[0])));
        TEST_ASSERT(tol.approximately_equal(fpp1, bez.fpp(t[0])));
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
          err+=eli::geom::point::distance(pts[i], bez.f(t[i]));
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0785);

        // check if went through points
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
        TEST_ASSERT(tol.approximately_equal(pts[0], bez.f(t[0])));
        TEST_ASSERT(tol.approximately_equal(fp1, bez.fp(t[0])));
        TEST_ASSERT(tol.approximately_equal(fpp1, bez.fpp(t[0])));
        TEST_ASSERT(tol.approximately_equal(pts[10], bez.f(t[10])));
        TEST_ASSERT(tol.approximately_equal(fp2, bez.fp(t[10])));
        TEST_ASSERT(tol.approximately_equal(fpp2, bez.fpp(t[10])));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));

//         octave_print(11, pts, bez);
      }
    }

    void fit_closed_test()
    {
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
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
        TEST_ASSERT(bez.fp(0)!=bez.fp(1));
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

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
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
        TEST_ASSERT(tol.approximately_equal(bez.fp(0), bez.fp(1)));
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

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
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
        TEST_ASSERT(tol.approximately_equal(bez.fp(0), bez.fp(1)));
        TEST_ASSERT(tol.approximately_equal(bez.fpp(0), bez.fpp(1)));

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
        TEST_ASSERT(tol.approximately_equal(bez.f(0), pts[0]));
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
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
        TEST_ASSERT(tol.approximately_equal(bez.f(0), pts[0]));
        TEST_ASSERT(tol.approximately_equal(bez.fp(0), fp0));
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
        TEST_ASSERT(tol.approximately_equal(bez.fp(0), bez.fp(1)));
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
        TEST_ASSERT(tol.approximately_equal(bez.f(0), pts[0]));
        TEST_ASSERT(tol.approximately_equal(bez.fp(0), fp0));
        TEST_ASSERT(tol.approximately_equal(bez.fpp(0), fpp0));
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
        TEST_ASSERT(tol.approximately_equal(bez.fp(0), bez.fp(1)));
        TEST_ASSERT(tol.approximately_equal(bez.fpp(0), bez.fpp(1)));

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
        TEST_ASSERT(tol.approximately_equal(bez.f(t5), f5));
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
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
        TEST_ASSERT(tol.approximately_equal(bez.f(t5), f5));
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
        TEST_ASSERT(tol.approximately_equal(bez.fp(0), bez.fp(1)));
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

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
        TEST_ASSERT(tol.approximately_equal(bez.f(t5), f5));
        TEST_ASSERT(tol.approximately_equal(bez.f(0), bez.f(1)));
        TEST_ASSERT(tol.approximately_equal(bez.fp(0), bez.fp(1)));
        TEST_ASSERT(tol.approximately_equal(bez.fpp(0), bez.fpp(1)));

//         TEST_ASSERT((bez.f(t5)-f5).norm()<4*eps);
//         TEST_ASSERT(bez.f(0)==bez.f(1));
//         TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<2.72e3*eps);
//         TEST_ASSERT((bez.fpp(1)-bez.fpp(0)).norm()<1.19e4*eps);
//         if (typeid(data__)==typeid(long double))
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
#ifdef ELI_USING_QD
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
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<13*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<37*eps);
        TEST_ASSERT(bez.f(0)==bez.f(1));

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
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<21*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<89*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<22*eps);
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<202*eps);
        TEST_ASSERT(bez.f(0)==bez.f(1));

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
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<2*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<77*eps);
        TEST_ASSERT((bez.fpp(t[1])-fpp1).norm()<353*eps);

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
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<93*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<45*eps);
        TEST_ASSERT((bez.fpp(t[1])-fpp1).norm()<398*eps);
        TEST_ASSERT(bez.f(0)==bez.f(1));
        TEST_ASSERT((bez.fp(0)-bez.fp(1)).norm()<5.30e3*eps);
        TEST_ASSERT((bez.fpp(0)-bez.fpp(1)).norm()<24.1e3*eps);

//         octave_print(23, pts, bez);
      }
    }
};

#endif

