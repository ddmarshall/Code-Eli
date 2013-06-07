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

#ifndef explicit_bezier_curve_fit_test_suite_hpp
#define explicit_bezier_curve_fit_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/curve/explicit_bezier.hpp"
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
class explicit_bezier_curve_fit_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::explicit_bezier<data__> curve_type;
    typedef typename curve_type::control_point_type control_point_type;
    typedef typename curve_type::point_type point_type;
    typedef typename curve_type::data_type data_type;
    typedef typename curve_type::index_type index_type;
    typedef typename curve_type::fit_container_type fit_container_type;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(explicit_bezier_curve_fit_test_suite<float>::fit_free_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<float>::fit_C0_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<float>::fit_C1_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<float>::fit_C2_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<float>::interpolate_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(explicit_bezier_curve_fit_test_suite<double>::fit_free_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<double>::fit_C0_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<double>::fit_C1_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<double>::fit_C2_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<double>::interpolate_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(explicit_bezier_curve_fit_test_suite<long double>::fit_free_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<long double>::fit_C0_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<long double>::fit_C1_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<long double>::fit_C2_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<long double>::interpolate_test);
    }
#ifdef ELI_USING_QD
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(explicit_bezier_curve_fit_test_suite<dd_real>::fit_free_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<dd_real>::fit_C0_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<dd_real>::fit_C1_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<dd_real>::fit_C2_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<dd_real>::interpolate_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(explicit_bezier_curve_fit_test_suite<qd_real>::fit_free_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<qd_real>::fit_C0_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<qd_real>::fit_C1_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<qd_real>::fit_C2_ends_test);
      TEST_ADD(explicit_bezier_curve_fit_test_suite<qd_real>::interpolate_test);
    }
#endif

  public:
    explicit_bezier_curve_fit_test_suite()
    {
      AddTests(data__());
    }
    ~explicit_bezier_curve_fit_test_suite()
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

    void create_circle(std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts)
    {
      // NOTE: This will create a semi-circle
      size_t n=pts.size();
      for (size_t i=0; i<n; ++i)
      {
        data__ theta(eli::constants::math<data__>::pi()*static_cast<data__>(i)/(n-1));
        pts[i](0)=(1-std::cos(theta))/2;
        pts[i](1)=std::sin(theta);
      }
    }

    void fit_free_ends_test()
    {
      fit_container_type fcon;
      size_t i, deg(5);
      std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(10);
      curve_type ebez;
      std::vector<data_type> t;
      data_type err_out, err_ref;

      // get points
      create_circle(pts);

      // configure fit container
      fcon.set_points(pts.begin(), pts.end());

      // fit points
      err_out=ebez.fit(t, fcon, deg);

      // calculate the error at the point
      for (err_ref=0, i=0; i<pts.size(); ++i)
      {
        err_ref+=eli::geom::point::distance(pts[i], ebez.f(t[i]));
      }

      TEST_ASSERT(err_ref==err_out);
      TEST_ASSERT(err_out < 0.430);

      // check if went through points
      TEST_ASSERT(pts[0]!=ebez.f(t[0]));
      TEST_ASSERT(pts[pts.size()-1]!=ebez.f(t[t.size()-1]));

//       octave_print(1, pts, ebez);
    }

    void fit_C0_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_USING_QD
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
      fit_container_type fcon;
      size_t i, deg(5);
      std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(10);
      curve_type ebez;
      std::vector<data_type> t;
      data_type err_out, err_ref;

      // get points
      create_circle(pts);

      // configure fit container
      fcon.set_points(pts.begin(), pts.end());
      fcon.add_start_C0_constraint();
      fcon.add_end_C0_constraint();

      // fit points
      err_out=ebez.fit(t, fcon, deg);

      // calculate the error at the point
      for (err_ref=0, i=0; i<pts.size(); ++i)
      {
        err_ref+=eli::geom::point::distance(pts[i], ebez.f(t[i]));
      }

      TEST_ASSERT(err_out==err_ref);
      TEST_ASSERT(err_out < 0.431);

      // check if went through points
      TEST_ASSERT(pts[0]==ebez.f(t[0]));
      TEST_ASSERT((pts[pts.size()-1]-ebez.f(t[t.size()-1])).norm()<65*eps);

//       octave_print(2, pts, ebez);
    }

    void fit_C1_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_USING_QD
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
      fit_container_type fcon;
      size_t i, deg(5);
      std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(10);
      curve_type ebez;
      std::vector<data_type> t;
      data_type err_out, err_ref;

      // get points
      create_circle(pts);

      // configure fit container
      point_type fp;
      typename fit_container_type::constraint_info ci;
      fcon.set_points(pts.begin(), pts.end());
      fcon.add_start_C0_constraint();
      fcon.add_end_C0_constraint();
      fcon.add_C1_constraint(4);
      fcon.get_constraint(4, ci);
      fp << 1, ci.get_fp()(0);

      // fit points
      err_out=ebez.fit(t, fcon, deg);

      // calculate the error at the point
      for (err_ref=0, i=0; i<pts.size(); ++i)
      {
        err_ref+=eli::geom::point::distance(pts[i], ebez.f(t[i]));
      }

      TEST_ASSERT(err_out==err_ref);
      TEST_ASSERT(err_out < 0.682);

      // check if went through points
      TEST_ASSERT(pts[0]==ebez.f(t[0]));
      TEST_ASSERT((pts[pts.size()-1]-ebez.f(t[t.size()-1])).norm()<40*eps);
      TEST_ASSERT((pts[4]-ebez.f(t[4])).norm()<35*eps);
      TEST_ASSERT((fp-ebez.fp(t[4])).norm()<48*eps);

//       octave_print(3, pts, ebez);
    }

    void fit_C2_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_USING_QD
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
      fit_container_type fcon;
      size_t i, deg(7);
      std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(10);
      curve_type ebez;
      std::vector<data_type> t;
      data_type err_out, err_ref;

      // get points
      create_circle(pts);

      // configure fit container
      point_type fp, fpp;
      typename fit_container_type::constraint_info ci;
      fcon.set_points(pts.begin(), pts.end());
      fcon.add_start_C0_constraint();
      fcon.add_end_C0_constraint();
      fcon.add_C2_constraint(4);
      fcon.get_constraint(4, ci);
      fp << 1, ci.get_fp()(0);
      fpp << 0, ci.get_fpp()(0);

      // fit points
      err_out=ebez.fit(t, fcon, deg);

      // calculate the error at the point
      for (err_ref=0, i=0; i<pts.size(); ++i)
      {
        err_ref+=eli::geom::point::distance(pts[i], ebez.f(t[i]));
      }

      TEST_ASSERT(err_out==err_ref);
      TEST_ASSERT(err_out < 0.310);

      // check if went through points
      TEST_ASSERT((pts[0]-ebez.f(t[0])).norm()<6*eps);
      TEST_ASSERT((pts[pts.size()-1]-ebez.f(t[t.size()-1])).norm()<512*eps);
      TEST_ASSERT((pts[4]-ebez.f(t[4])).norm()<35*eps);
      TEST_ASSERT((fp-ebez.fp(t[4])).norm()<240*eps);
      TEST_ASSERT((fpp-ebez.fpp(t[4])).norm()<1.54e3*eps);

//       octave_print(4, pts, ebez);
    }

    void interpolate_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_USING_QD
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
      // interpolate through all points
      {
        fit_container_type fcon;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        curve_type ebez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());

        // fit points
        ebez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(ebez.degree()+1)==pts.size());

        // check if went through points
        TEST_ASSERT(ebez.f(t[0])==pts[0]);
        TEST_ASSERT((ebez.f(t[1])-pts[1]).norm()<3*eps);
        TEST_ASSERT((ebez.f(t[2])-pts[2]).norm()<7*eps);
        TEST_ASSERT((ebez.f(t[3])-pts[3]).norm()<23*eps);

//         octave_print(5, pts, ebez);
      }

      // interpolate through all points and specify C1
      {
        fit_container_type fcon;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        curve_type ebez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp;
        typename fit_container_type::constraint_info ci;
        fcon.set_points(pts.begin(), pts.end());
        fcon.add_C1_constraint(1);
        fcon.get_constraint(1, ci);
        fp << 1, ci.get_fp()(0);

        // fit points
        ebez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(ebez.degree()+1)==(pts.size()+1));

        // check if went through points
        TEST_ASSERT((ebez.f(t[0])-pts[0]).norm()<2*eps);
        TEST_ASSERT((ebez.f(t[1])-pts[1]).norm()<3*eps);
        TEST_ASSERT((ebez.f(t[2])-pts[2]).norm()<9*eps);
        TEST_ASSERT((ebez.f(t[3])-pts[3]).norm()<44*eps);
        TEST_ASSERT((fp-ebez.fp(t[1])).norm()<48*eps);

//         octave_print(6, pts, ebez);
      }

      // interpolate through all points and specify C2
      {
        fit_container_type fcon;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        curve_type ebez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp, fpp;
        typename fit_container_type::constraint_info ci;
        fcon.set_points(pts.begin(), pts.end());
        fcon.add_C2_constraint(1);
        fcon.get_constraint(1, ci);
        fp << 1, ci.get_fp()(0);
        fpp << 0, ci.get_fpp()(0);

        // fit points
        ebez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(ebez.degree()+1)==(pts.size()+2));

        // check if went through points
        TEST_ASSERT((ebez.f(t[0])-pts[0]).norm()<3*eps);
        TEST_ASSERT((ebez.f(t[1])-pts[1]).norm()<3*eps);
        TEST_ASSERT((ebez.f(t[2])-pts[2]).norm()<9*eps);
        TEST_ASSERT((ebez.f(t[3])-pts[3]).norm()<70*eps);
        TEST_ASSERT((fp-ebez.fp(t[1])).norm()<48*eps);
        TEST_ASSERT((fpp-ebez.fpp(t[1])).norm()<81*eps);

//         octave_print(7, pts, ebez);
      }
    }
};
#endif

