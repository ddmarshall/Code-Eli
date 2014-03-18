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

#ifndef piecewise_four_series_creator_test_suite_hpp
#define piecewise_four_series_creator_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/mutil/fd/d1o2.hpp"
#include "eli/mutil/fd/d2o2.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_four_series_creator.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_four_series_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::piecewise_four_series_creator<data__, 3, tolerance_type> point_creator_type;

    typedef eli::geom::airfoil::four_series<data_type> airfoil_type;
    typedef typename airfoil_type::coefficient_type airfoil_thickness_coefficient_type;
    typedef typename airfoil_type::point_type airfoil_point_type;
    typedef typename airfoil_type::param_type aifoil_parameter_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_four_series_creator_test_suite<float>::airfoil_coefficients_test);
      TEST_ADD(piecewise_four_series_creator_test_suite<float>::create_airfoil_test);
      TEST_ADD(piecewise_four_series_creator_test_suite<float>::thickness_test);
      TEST_ADD(piecewise_four_series_creator_test_suite<float>::thickness_derivatives_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_four_series_creator_test_suite<double>::airfoil_coefficients_test);
      TEST_ADD(piecewise_four_series_creator_test_suite<double>::create_airfoil_test);
      TEST_ADD(piecewise_four_series_creator_test_suite<double>::thickness_test);
      TEST_ADD(piecewise_four_series_creator_test_suite<double>::thickness_derivatives_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_four_series_creator_test_suite<long double>::airfoil_coefficients_test);
      TEST_ADD(piecewise_four_series_creator_test_suite<long double>::create_airfoil_test);
      TEST_ADD(piecewise_four_series_creator_test_suite<long double>::thickness_test);
      TEST_ADD(piecewise_four_series_creator_test_suite<long double>::thickness_derivatives_test);
    }

  public:
    piecewise_four_series_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_four_series_creator_test_suite()
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

    void airfoil_coefficients_test()
    {
      airfoil_type af;
      airfoil_thickness_coefficient_type a, a_ref;

      af.set_sharp_trailing_edge(false);
      a_ref << 0.2969, -0.1260, -0.3516, 0.2843, -0.1015;
      a=af.get_thickness_coefficients();
      TEST_ASSERT((a-a_ref).norm()<1e-4);

      af.set_sharp_trailing_edge(true);
      a_ref << 0.2983, -0.1325, -0.3286, 0.2442, -0.0815;
      a=af.get_thickness_coefficients();
      TEST_ASSERT((a-a_ref).norm()<1e-4);
    }

    void create_airfoil_test()
    {
      airfoil_type af;
      aifoil_parameter_type th, cam, cam_loc;
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

      name_ref="NACA "+std::to_string(cam)+std::to_string(cam_loc)+std::to_string(th);
      name=af.get_name();
      TEST_ASSERT(name==name_ref);
    }

    void thickness_test()
    {
      // non-sharp trailing edge
      {
        airfoil_type af;
        std::vector<data_type> x(18), y(18);
        airfoil_point_type xout;
        const data_type tol(1e-5);

        // create 0021 airfoil to test against Abbott & von Doenhoff data NACA Report 824 p. 330
        af.set_thickness(21);
        x[0]  = static_cast<data_type>(0);       y[0]  = static_cast<data_type>(0);
        x[1]  = static_cast<data_type>(1.25e-2); y[1]  = static_cast<data_type>(3.315e-2);
        x[2]  = static_cast<data_type>(2.5e-2);  y[2]  = static_cast<data_type>(4.576e-2);
        x[3]  = static_cast<data_type>(5e-2);    y[3]  = static_cast<data_type>(6.221e-2);
        x[4]  = static_cast<data_type>(7.5e-2);  y[4]  = static_cast<data_type>(7.350e-2);
        x[5]  = static_cast<data_type>(10e-2);   y[5]  = static_cast<data_type>(8.195e-2);
        x[6]  = static_cast<data_type>(15e-2);   y[6]  = static_cast<data_type>(9.354e-2);
        x[7]  = static_cast<data_type>(20e-2);   y[7]  = static_cast<data_type>(10.040e-2);
        x[8]  = static_cast<data_type>(25e-2);   y[8]  = static_cast<data_type>(10.397e-2);
        x[9]  = static_cast<data_type>(30e-2);   y[9]  = static_cast<data_type>(10.504e-2);
        x[10] = static_cast<data_type>(40e-2);   y[10] = static_cast<data_type>(10.156e-2);
        x[11] = static_cast<data_type>(50e-2);   y[11] = static_cast<data_type>(9.265e-2);
        x[12] = static_cast<data_type>(60e-2);   y[12] = static_cast<data_type>(7.986e-2);
        x[13] = static_cast<data_type>(70e-2);   y[13] = static_cast<data_type>(6.412e-2);
        x[14] = static_cast<data_type>(80e-2);   y[14] = static_cast<data_type>(4.591e-2);
        x[15] = static_cast<data_type>(90e-2);   y[15] = static_cast<data_type>(2.534e-2);
        x[16] = static_cast<data_type>(95e-2);   y[16] = static_cast<data_type>(1.412e-2);
        x[17] = static_cast<data_type>(1);       y[17] = static_cast<data_type>(0.221e-2);

        // test the airfoil coordinates
        for (size_t i=0; i<x.size(); ++i)
        {
          xout = af.f(x[i]);
          TEST_ASSERT_DELTA(x[i], xout(0), tol);
          TEST_ASSERT_DELTA(y[i], xout(1), tol);
        }
      }

      // sharp trailing edge
      {
        airfoil_type af;
        std::vector<data_type> x(18), y(18);
        airfoil_point_type xout;
        const data_type tol(1e-5);

        // create 0021 airfoil to test against Abbott & von Doenhoff data NACA Report 824 p. 330
        af.set_thickness(21);
        af.set_sharp_trailing_edge(true);
        x[0]  = static_cast<data_type>(0);       y[0]  = static_cast<data_type>(0);
        x[1]  = static_cast<data_type>(1.25e-2); y[1]  = static_cast<data_type>(3.323e-2);
        x[2]  = static_cast<data_type>(2.5e-2);  y[2]  = static_cast<data_type>(4.584e-2);
        x[3]  = static_cast<data_type>(5e-2);    y[3]  = static_cast<data_type>(6.226e-2);
        x[4]  = static_cast<data_type>(7.5e-2);  y[4]  = static_cast<data_type>(7.352e-2);
        x[5]  = static_cast<data_type>(10e-2);   y[5]  = static_cast<data_type>(8.195e-2);
        x[6]  = static_cast<data_type>(15e-2);   y[6]  = static_cast<data_type>(9.352e-2);
        x[7]  = static_cast<data_type>(20e-2);   y[7]  = static_cast<data_type>(10.039e-2);
        x[8]  = static_cast<data_type>(25e-2);   y[8]  = static_cast<data_type>(10.397e-2);
        x[9]  = static_cast<data_type>(30e-2);   y[9]  = static_cast<data_type>(10.503e-2);
        x[10] = static_cast<data_type>(40e-2);   y[10] = static_cast<data_type>(10.150e-2);
        x[11] = static_cast<data_type>(50e-2);   y[11] = static_cast<data_type>(9.241e-2);
        x[12] = static_cast<data_type>(60e-2);   y[12] = static_cast<data_type>(7.929e-2);
        x[13] = static_cast<data_type>(70e-2);   y[13] = static_cast<data_type>(6.308e-2);
        x[14] = static_cast<data_type>(80e-2);   y[14] = static_cast<data_type>(4.434e-2);
        x[15] = static_cast<data_type>(90e-2);   y[15] = static_cast<data_type>(2.333e-2);
        x[16] = static_cast<data_type>(95e-2);   y[16] = static_cast<data_type>(1.196e-2);
        x[17] = static_cast<data_type>(1);       y[17] = static_cast<data_type>(0);

        // test the airfoil coordinates
        for (size_t i=0; i<x.size(); ++i)
        {
          xout = af.f(x[i]);
          TEST_ASSERT_DELTA(x[i], xout(0), tol);
          TEST_ASSERT_DELTA(y[i], xout(1), tol);
//          std::cout << "i=" << i << "\tycal=" << xout(1) << "\tyref=" << y[i] << "\tdy=" << xout(1)-y[i] << std::endl;
        }
      }
    }

    void thickness_derivatives_test()
    {
      // test the known derivatives of thick trailing edge geometry
      {
        airfoil_type af;
        airfoil_point_type x, xp, x_ref, xp_ref;
        data_type xi;

        af.set_thickness(20);

        // lower surface trailing edge
        xi=static_cast<data_type>(-1);
        x_ref << -xi, static_cast<data_type>(-0.002);
        xp_ref << static_cast<data_type>(-1), static_cast<data_type>(-0.234);
        x=af.f(xi);
        xp=af.fp(xi);
        TEST_ASSERT((x-x_ref).norm()<2e-4);
        TEST_ASSERT((xp-xp_ref).norm()<2e-4);

        // lower surface max thickness
        xi=static_cast<data_type>(-0.3);
        x_ref << -xi, static_cast<data_type>(-0.1);
        xp_ref << static_cast<data_type>(-1), static_cast<data_type>(0);
        x=af.f(xi);
        xp=af.fp(xi);
        TEST_ASSERT((x-x_ref).norm()<3e-5);
        TEST_ASSERT((xp-xp_ref).norm()<2e-4);

        // leading edge
        xi=static_cast<data_type>(0);
        x_ref << xi, static_cast<data_type>(0);
        xp_ref << static_cast<data_type>(0), static_cast<data_type>(1)/std::numeric_limits<data_type>::epsilon();
        x=af.f(xi);
        xp=af.fp(xi);
        TEST_ASSERT(x==x_ref);
        TEST_ASSERT(xp==xp_ref);

        // upper surface max thickness
        xi=static_cast<data_type>(0.3);
        x_ref << xi, static_cast<data_type>(0.1);
        xp_ref << static_cast<data_type>(1), static_cast<data_type>(0);
        x=af.f(xi);
        xp=af.fp(xi);
        TEST_ASSERT((x-x_ref).norm()<3e-5);
        TEST_ASSERT((xp-xp_ref).norm()<2e-4);

        // upper surface trailing edge
        xi=static_cast<data_type>(1);
        x_ref << xi, static_cast<data_type>(0.002);
        xp_ref << static_cast<data_type>(1), static_cast<data_type>(-0.234);
        x=af.f(xi);
        xp=af.fp(xi);
        TEST_ASSERT((x-x_ref).norm()<2e-4);
        TEST_ASSERT((xp-xp_ref).norm()<2e-4);
//        std::cout << "xi=" << xi;
//        std::cout << "\tx=" << x << "\tx_ref=" << x_ref;
//        std::cout << "\tdx=" << x-x_ref << "\tdx.norm()=" << (x-x_ref).norm();
//        std::cout << "\txp=" << xp << "\txp_ref=" << xp_ref;
//        std::cout << "\tdxp=" << xp-xp_ref << "\tdxp.norm()=" << (xp-xp_ref).norm();
//        std::cout << std::endl;
      }

      // use finite differences to calculate the derivatives
      {
        airfoil_type af;
        airfoil_point_type xr, xp, xpp, xp_ref, xpp_ref;
        data_type xi, x[3], y[3];
        const data_type dxi(1e2*std::sqrt(std::numeric_limits<data_type>::epsilon()));
        eli::mutil::fd::d1o2<data_type> d1_calc;
        eli::mutil::fd::d2o2<data_type> d2_calc;

        af.set_thickness(12);

        // lower surface
        xi=static_cast<data_type>(-0.8);
        xr=af.f(xi-dxi); x[0]=xr(0); y[0]=xr(1);
        xr=af.f(xi);     x[1]=xr(0); y[1]=xr(1);
        xr=af.f(xi+dxi); x[2]=xr(0); y[2]=xr(1);
        d1_calc.evaluate(xp_ref(0), x, dxi);
        d1_calc.evaluate(xp_ref(1), y, dxi);
        d2_calc.evaluate(xpp_ref(0), x, dxi);
        d2_calc.evaluate(xpp_ref(1), y, dxi);
        xp=af.fp(xi);
        xpp=af.fpp(xi);
        if (typeid(data_type)==typeid(float))
        {
          TEST_ASSERT((xp-xp_ref).norm()<1e-3);
          TEST_ASSERT((xpp-xpp_ref).norm()<3e-3);
        }
        else
        {
          TEST_ASSERT((xp-xp_ref).norm()<2e-4);
          TEST_ASSERT((xpp-xpp_ref).norm()<2e-4);
        }

        xi=static_cast<data_type>(-0.25);
        xr=af.f(xi-dxi); x[0]=xr(0); y[0]=xr(1);
        xr=af.f(xi);     x[1]=xr(0); y[1]=xr(1);
        xr=af.f(xi+dxi); x[2]=xr(0); y[2]=xr(1);
        d1_calc.evaluate(xp_ref(0), x, dxi);
        d1_calc.evaluate(xp_ref(1), y, dxi);
        d2_calc.evaluate(xpp_ref(0), x, dxi);
        d2_calc.evaluate(xpp_ref(1), y, dxi);
        xp=af.fp(xi);
        xpp=af.fpp(xi);
        if (typeid(data_type)==typeid(float))
        {
          TEST_ASSERT((xp-xp_ref).norm()<1e-3);
          TEST_ASSERT((xpp-xpp_ref).norm()<3e-3);
        }
        else
        {
          TEST_ASSERT((xp-xp_ref).norm()<2e-4);
          TEST_ASSERT((xpp-xpp_ref).norm()<2e-4);
        }

        // upper surface
        xi=static_cast<data_type>(0.15);
        xr=af.f(xi-dxi); x[0]=xr(0); y[0]=xr(1);
        xr=af.f(xi);     x[1]=xr(0); y[1]=xr(1);
        xr=af.f(xi+dxi); x[2]=xr(0); y[2]=xr(1);
        d1_calc.evaluate(xp_ref(0), x, dxi);
        d1_calc.evaluate(xp_ref(1), y, dxi);
        d2_calc.evaluate(xpp_ref(0), x, dxi);
        d2_calc.evaluate(xpp_ref(1), y, dxi);
        xp=af.fp(xi);
        xpp=af.fpp(xi);
        if (typeid(data_type)==typeid(float))
        {
          TEST_ASSERT((xp-xp_ref).norm()<2e-3);
          TEST_ASSERT((xpp-xpp_ref).norm()<2e-2);
        }
        else
        {
          TEST_ASSERT((xp-xp_ref).norm()<2e-4);
          TEST_ASSERT((xpp-xpp_ref).norm()<2e-4);
        }

        xi=static_cast<data_type>(0.7);
        xr=af.f(xi-dxi); x[0]=xr(0); y[0]=xr(1);
        xr=af.f(xi);     x[1]=xr(0); y[1]=xr(1);
        xr=af.f(xi+dxi); x[2]=xr(0); y[2]=xr(1);
        d1_calc.evaluate(xp_ref(0), x, dxi);
        d1_calc.evaluate(xp_ref(1), y, dxi);
        d2_calc.evaluate(xpp_ref(0), x, dxi);
        d2_calc.evaluate(xpp_ref(1), y, dxi);
        xp=af.fp(xi);
        xpp=af.fpp(xi);
        TEST_ASSERT((xp-xp_ref).norm()<2e-4);
        TEST_ASSERT((xpp-xpp_ref).norm()<3e-4);
//        std::cout << "xi=" << xi;
//        std::cout << "\tdx=" << (x[2]-x[0])/(2*dxi) << "\tdy=" << (y[2]-y[0])/(2*dxi);
//        std::cout << std::endl << "     ";
//        std::cout << "ddx=" << (x[0]-2*x[1]+x[2])/(dxi*dxi) << "\tddy=" << (y[0]-2*y[1]+y[2])/(dxi*dxi);
//        std::cout << std::endl << "     ";
//        std::cout << "xp=" << xp;
//        std::cout << std::endl << "     ";
//        std::cout << "xp_ref=" << xp_ref;
//        std::cout << std::endl << "     ";
//        std::cout << "dxp=" << xp-xp_ref << "dxp.norm()=" << (xp-xp_ref).norm();
//        std::cout << std::endl << "     ";
//        std::cout << "xpp=    " << xpp;
//        std::cout << std::endl << "     ";
//        std::cout << "xpp_ref=" << xpp_ref;
//        std::cout << std::endl << "     ";
//        std::cout << "dxpp=" << xpp-xpp_ref << "\tdxpp.norm()=" << (xpp-xpp_ref).norm();
//        std::cout << std::endl;
      }
    }
};

#endif

