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

#ifndef four_digit_test_suite_hpp
#define four_digit_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/mutil/fd/d1o2.hpp"
#include "eli/mutil/fd/d2o2.hpp"
#include "eli/geom/curve/pseudo/four_digit.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class four_digit_test_suite : public Test::Suite
{
  private:
    typedef data__ data_type;
    typedef eli::geom::curve::pseudo::four_digit<data_type> airfoil_type;
    typedef typename airfoil_type::coefficient_type airfoil_thickness_coefficient_type;
    typedef typename airfoil_type::point_type airfoil_point_type;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(four_digit_test_suite<float>::airfoil_coefficients_test);
      TEST_ADD(four_digit_test_suite<float>::create_airfoil_test);
      TEST_ADD(four_digit_test_suite<float>::thickness_test);
      TEST_ADD(four_digit_test_suite<float>::thickness_derivatives_test);
      TEST_ADD(four_digit_test_suite<float>::camber_test);
      TEST_ADD(four_digit_test_suite<float>::camber_derivatives_test);
      TEST_ADD(four_digit_test_suite<float>::airfoil_test);
      TEST_ADD(four_digit_test_suite<float>::airfoil_derivatives_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(four_digit_test_suite<double>::airfoil_coefficients_test);
      TEST_ADD(four_digit_test_suite<double>::create_airfoil_test);
      TEST_ADD(four_digit_test_suite<double>::thickness_test);
      TEST_ADD(four_digit_test_suite<double>::thickness_derivatives_test);
      TEST_ADD(four_digit_test_suite<double>::camber_test);
      TEST_ADD(four_digit_test_suite<double>::camber_derivatives_test);
      TEST_ADD(four_digit_test_suite<double>::airfoil_test);
      TEST_ADD(four_digit_test_suite<double>::airfoil_derivatives_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(four_digit_test_suite<long double>::airfoil_coefficients_test);
      TEST_ADD(four_digit_test_suite<long double>::create_airfoil_test);
      TEST_ADD(four_digit_test_suite<long double>::thickness_test);
      TEST_ADD(four_digit_test_suite<long double>::thickness_derivatives_test);
      TEST_ADD(four_digit_test_suite<long double>::camber_test);
      TEST_ADD(four_digit_test_suite<long double>::camber_derivatives_test);
      TEST_ADD(four_digit_test_suite<long double>::airfoil_test);
      TEST_ADD(four_digit_test_suite<long double>::airfoil_derivatives_test);
    }

  public:
    four_digit_test_suite()
    {
      AddTests(data__());
    }
    ~four_digit_test_suite()
    {
    }

  private:
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
      data_type th, cam, cam_loc;
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

      name_ref="NACA "+std::to_string((int)round(cam))+std::to_string((int)round(cam_loc))+std::to_string((int)round(th));
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
          TEST_ASSERT((xp-xp_ref).norm()<1e-5);
          TEST_ASSERT((xpp-xpp_ref).norm()<3e-3);
        }
        else
        {
          TEST_ASSERT((xp-xp_ref).norm()<1e-10);
          TEST_ASSERT((xpp-xpp_ref).norm()<5e-5);
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
          TEST_ASSERT((xp-xp_ref).norm()<6e-4);
          TEST_ASSERT((xpp-xpp_ref).norm()<3e-3);
        }
        else
        {
          TEST_ASSERT((xp-xp_ref).norm()<1e-10);
          TEST_ASSERT((xpp-xpp_ref).norm()<8e-6);
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
          TEST_ASSERT((xp-xp_ref).norm()<1e-10);
          TEST_ASSERT((xpp-xpp_ref).norm()<5e-6);
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
        if (typeid(data_type)==typeid(float))
        {
          TEST_ASSERT((xp-xp_ref).norm()<5e-5);
          TEST_ASSERT((xpp-xpp_ref).norm()<5e-4);
        }
        else
        {
          TEST_ASSERT((xp-xp_ref).norm()<1e-10);
          TEST_ASSERT((xpp-xpp_ref).norm()<8e-5);
        }
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

    void camber_test()
    {
      airfoil_type af;
      std::vector<data_type> x(18), y(18);
      airfoil_point_type xout;
      const data_type tol(1e-5);

      // create 2300 airfoil to test against externally calculated data
      af.set_camber(2, 3);
      af.set_thickness(00);
      x[ 0]=static_cast<data_type>(0);        y[ 0]=static_cast<data_type>(0);
      x[ 1]=static_cast<data_type>(1.25e-2);  y[ 1]=static_cast<data_type>(0.1632e-2);
      x[ 2]=static_cast<data_type>(2.50e-2);  y[ 2]=static_cast<data_type>(0.3194e-2);
      x[ 3]=static_cast<data_type>(5.00e-2);  y[ 3]=static_cast<data_type>(0.6111e-2);
      x[ 4]=static_cast<data_type>(7.50e-2);  y[ 4]=static_cast<data_type>(0.8750e-2);
      x[ 5]=static_cast<data_type>(10.00e-2); y[ 5]=static_cast<data_type>(1.1111e-2);
      x[ 6]=static_cast<data_type>(15.00e-2); y[ 6]=static_cast<data_type>(1.5000e-2);
      x[ 7]=static_cast<data_type>(20.00e-2); y[ 7]=static_cast<data_type>(1.7778e-2);
      x[ 8]=static_cast<data_type>(25.00e-2); y[ 8]=static_cast<data_type>(1.9444e-2);
      x[ 9]=static_cast<data_type>(30.00e-2); y[ 9]=static_cast<data_type>(2.0000e-2);
      x[10]=static_cast<data_type>(40.00e-2); y[10]=static_cast<data_type>(1.9592e-2);
      x[11]=static_cast<data_type>(50.00e-2); y[11]=static_cast<data_type>(1.8367e-2);
      x[12]=static_cast<data_type>(60.00e-2); y[12]=static_cast<data_type>(1.6327e-2);
      x[13]=static_cast<data_type>(70.00e-2); y[13]=static_cast<data_type>(1.3469e-2);
      x[14]=static_cast<data_type>(80.00e-2); y[14]=static_cast<data_type>(0.9796e-2);
      x[15]=static_cast<data_type>(90.00e-2); y[15]=static_cast<data_type>(0.5306e-2);
      x[16]=static_cast<data_type>(95.00e-2); y[16]=static_cast<data_type>(0.2755e-2);
      x[17]=static_cast<data_type>(1);        y[17]=static_cast<data_type>(0);

      // test the airfoil coordinates
      for (size_t i=0; i<x.size(); ++i)
      {
        xout = af.f(x[i]);
        TEST_ASSERT_DELTA(x[i], xout(0), tol);
        TEST_ASSERT_DELTA(y[i], xout(1), tol);
//        std::cout << "i=" << i << "\tycal=" << xout(1) << "\tyref=" << y[i] << "\tdy=" << xout(1)-y[i] << std::endl;
      }
    }

    void camber_derivatives_test()
    {
      // test the known derivatives
      {
        airfoil_type af;
        airfoil_point_type x, xp, x_ref, xp_ref;
        data_type xi;

        // configure camber line
        af.set_camber(3, 2);
        af.set_thickness(00);

        // lower surface max camber
        xi=static_cast<data_type>(-0.2);
        x_ref << -xi, static_cast<data_type>(0.03);
        xp_ref << static_cast<data_type>(-1), static_cast<data_type>(0);
        x=af.f(xi);
        xp=af.fp(xi);
        TEST_ASSERT((x-x_ref).norm()<1e-8);
        TEST_ASSERT((xp-xp_ref).norm()<1e-8);

        // upper surface max thickness
        xi=static_cast<data_type>(0.2);
        x_ref << xi, static_cast<data_type>(0.03);
        xp_ref << static_cast<data_type>(1), static_cast<data_type>(0);
        x=af.f(xi);
        xp=af.fp(xi);
        TEST_ASSERT((x-x_ref).norm()<1e-8);
        TEST_ASSERT((xp-xp_ref).norm()<1e-8);
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

        // configure camber line
        af.set_camber(3, 2);
        af.set_thickness(00);

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
          TEST_ASSERT((xp-xp_ref).norm()<1e-6);
          TEST_ASSERT((xpp-xpp_ref).norm()<1e-6);
        }
        else
        {
          TEST_ASSERT((xp-xp_ref).norm()<1e-10);
          TEST_ASSERT((xpp-xpp_ref).norm()<1e-10);
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
          TEST_ASSERT((xp-xp_ref).norm()<1e-6);
          TEST_ASSERT((xpp-xpp_ref).norm()<5e-6);
        }
        else
        {
          TEST_ASSERT((xp-xp_ref).norm()<1e-10);
          TEST_ASSERT((xpp-xpp_ref).norm()<5e-6);
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
          TEST_ASSERT((xp-xp_ref).norm()<1e-6);
          TEST_ASSERT((xpp-xpp_ref).norm()<5e-6);
        }
        else
        {
          TEST_ASSERT((xp-xp_ref).norm()<1e-10);
          TEST_ASSERT((xpp-xpp_ref).norm()<2e-6);
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
        if (typeid(data_type)==typeid(float))
        {
          TEST_ASSERT((xp-xp_ref).norm()<2e-3);
          TEST_ASSERT((xpp-xpp_ref).norm()<2e-2);
        }
        else
        {
          TEST_ASSERT((xp-xp_ref).norm()<1e-10);
          TEST_ASSERT((xpp-xpp_ref).norm()<7e-5);
        }
//        std::cout << "xi=" << xi;
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

    void airfoil_test()
    {
      airfoil_type af;
      std::vector<data_type> xi(36), x(36), y(36);
      airfoil_point_type xout;
      const data_type tol(3e-3);

      // create 2412 airfoil to test against Abbott & von Doenhoff data NACA Report 824 p. 358
      af.set_camber(2, 4);
      af.set_thickness(12);
      xi[ 0] = static_cast<data_type>(       0.0e-2); x[ 0] = static_cast<data_type>(   0.0e-2); y[ 0] = static_cast<data_type>(  0.0e-2);
      xi[ 1] = static_cast<data_type>( -1.078935e-2); x[ 1] = static_cast<data_type>(  1.25e-2); y[ 1] = static_cast<data_type>(-1.65e-2); 
      xi[ 2] = static_cast<data_type>( -2.265269e-2); x[ 2] = static_cast<data_type>(  2.50e-2); y[ 2] = static_cast<data_type>(-2.27e-2); 
      xi[ 3] = static_cast<data_type>( -4.695760e-2); x[ 3] = static_cast<data_type>(  5.00e-2); y[ 3] = static_cast<data_type>(-3.01e-2); 
      xi[ 4] = static_cast<data_type>( -7.162585e-2); x[ 4] = static_cast<data_type>(  7.50e-2); y[ 4] = static_cast<data_type>(-3.46e-2); 
      xi[ 5] = static_cast<data_type>( -9.650263e-2); x[ 5] = static_cast<data_type>( 10.00e-2); y[ 5] = static_cast<data_type>(-3.75e-2);
      xi[ 6] = static_cast<data_type>(-14.664316e-2); x[ 6] = static_cast<data_type>( 15.00e-2); y[ 6] = static_cast<data_type>(-4.10e-2); 
      xi[ 7] = static_cast<data_type>(-19.710203e-2); x[ 7] = static_cast<data_type>( 20.00e-2); y[ 7] = static_cast<data_type>(-4.23e-2); 
      xi[ 8] = static_cast<data_type>(-24.774236e-2); x[ 8] = static_cast<data_type>( 25.00e-2); y[ 8] = static_cast<data_type>(-4.22e-2); 
      xi[ 9] = static_cast<data_type>(-29.847722e-2); x[ 9] = static_cast<data_type>( 30.00e-2); y[ 9] = static_cast<data_type>(-4.12e-2);
      xi[10] = static_cast<data_type>(-40.000000e-2); x[10] = static_cast<data_type>( 40.00e-2); y[10] = static_cast<data_type>(-3.80e-2); 
      xi[11] = static_cast<data_type>(-50.059125e-2); x[11] = static_cast<data_type>( 50.00e-2); y[11] = static_cast<data_type>(-3.34e-2);
      xi[12] = static_cast<data_type>(-60.101712e-2); x[12] = static_cast<data_type>( 60.00e-2); y[12] = static_cast<data_type>(-2.76e-2);
      xi[13] = static_cast<data_type>(-70.122161e-2); x[13] = static_cast<data_type>( 70.00e-2); y[13] = static_cast<data_type>(-2.14e-2);
      xi[14] = static_cast<data_type>(-80.116232e-2); x[14] = static_cast<data_type>( 80.00e-2); y[14] = static_cast<data_type>(-1.50e-2); 
      xi[15] = static_cast<data_type>(-90.079880e-2); x[15] = static_cast<data_type>( 90.00e-2); y[15] = static_cast<data_type>(-0.82e-2);
      xi[16] = static_cast<data_type>(-95.048848e-2); x[16] = static_cast<data_type>( 95.00e-2); y[16] = static_cast<data_type>(-0.48e-2); 
      xi[17] = static_cast<data_type>(      -1.0e+0); x[17] = static_cast<data_type>( 99.99e-2); y[17] = static_cast<data_type>(-0.13e-2); 
      xi[18] = static_cast<data_type>(       0.0e-2); x[18] = static_cast<data_type>(   0.0e-2); y[18] = static_cast<data_type>(  0.0e-2);
      xi[19] = static_cast<data_type>(  1.444524e-2); x[19] = static_cast<data_type>(  1.25e-2); y[19] = static_cast<data_type>( 2.15e-2);
      xi[20] = static_cast<data_type>(  2.753309e-2); x[20] = static_cast<data_type>(  2.50e-2); y[20] = static_cast<data_type>( 2.99e-2); 
      xi[21] = static_cast<data_type>(  5.315146e-2); x[21] = static_cast<data_type>(  5.00e-2); y[21] = static_cast<data_type>( 4.13e-2); 
      xi[22] = static_cast<data_type>(  7.842503e-2); x[22] = static_cast<data_type>(  7.50e-2); y[22] = static_cast<data_type>( 4.96e-2);
      xi[23] = static_cast<data_type>( 10.350449e-2); x[23] = static_cast<data_type>( 10.00e-2); y[23] = static_cast<data_type>( 5.63e-2); 
      xi[24] = static_cast<data_type>( 15.331062e-2); x[24] = static_cast<data_type>( 15.00e-2); y[24] = static_cast<data_type>( 6.61e-2); 
      xi[25] = static_cast<data_type>( 20.283261e-2); x[25] = static_cast<data_type>( 20.00e-2); y[25] = static_cast<data_type>( 7.26e-2); 
      xi[26] = static_cast<data_type>( 25.219585e-2); x[26] = static_cast<data_type>( 25.00e-2); y[26] = static_cast<data_type>( 7.67e-2); 
      xi[27] = static_cast<data_type>( 30.147780e-2); x[27] = static_cast<data_type>( 30.00e-2); y[27] = static_cast<data_type>( 7.88e-2); 
      xi[28] = static_cast<data_type>( 40.000000e-2); x[28] = static_cast<data_type>( 40.00e-2); y[28] = static_cast<data_type>( 7.80e-2); 
      xi[29] = static_cast<data_type>( 49.941485e-2); x[29] = static_cast<data_type>( 50.00e-2); y[29] = static_cast<data_type>( 7.24e-2); 
      xi[30] = static_cast<data_type>( 59.898946e-2); x[30] = static_cast<data_type>( 60.00e-2); y[30] = static_cast<data_type>( 6.36e-2);
      xi[31] = static_cast<data_type>( 69.878040e-2); x[31] = static_cast<data_type>( 70.00e-2); y[31] = static_cast<data_type>( 5.18e-2); 
      xi[32] = static_cast<data_type>( 79.883299e-2); x[32] = static_cast<data_type>( 80.00e-2); y[32] = static_cast<data_type>( 3.75e-2); 
      xi[33] = static_cast<data_type>( 89.919268e-2); x[33] = static_cast<data_type>( 90.00e-2); y[33] = static_cast<data_type>( 2.08e-2);
      xi[34] = static_cast<data_type>( 94.950448e-2); x[34] = static_cast<data_type>( 95.00e-2); y[34] = static_cast<data_type>( 1.14e-2);
      xi[35] = static_cast<data_type>(       1.0e+0); x[35] = static_cast<data_type>(100.01e-2); y[35] = static_cast<data_type>( 0.13e-2); 

      // test the airfoil coordinates
      for (size_t i=0; i<x.size(); ++i)
      {
        xout = af.f(xi[i]);
        TEST_ASSERT_DELTA(x[i], xout(0), tol);
        TEST_ASSERT_DELTA(y[i], xout(1), tol);
//        std::cout << "i=" << i << std::endl;
//        std::cout << "      " << "\txcal=" << xout(0)
//                  << "\txref=" << x[i]
//                  << "\tdx=" << xout(0)-x[i] << std::endl;
//        std::cout << "      " << "\tycal=" << xout(1)
//                  << "\tyref=" << y[i]
//                  << "\tdy=" << xout(1)-y[i] << std::endl;
      }
    }

    void airfoil_derivatives_test()
    {
      // use finite differences to calculate the derivatives
      airfoil_type af;
      airfoil_point_type xr, xp, xpp, xp_ref, xpp_ref;
      data_type xi, x[3], y[3];
      const data_type dxi(1e2*std::sqrt(std::numeric_limits<data_type>::epsilon()));
      eli::mutil::fd::d1o2<data_type> d1_calc;
      eli::mutil::fd::d2o2<data_type> d2_calc;

      // configure camber line
      af.set_camber(3, 2);
      af.set_thickness(24);

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
        TEST_ASSERT((xp-xp_ref).norm()<2e-5);
        TEST_ASSERT((xpp-xpp_ref).norm()<5e-4);
      }else
      {
        TEST_ASSERT((xp-xp_ref).norm()<1e-10);
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
        TEST_ASSERT((xp-xp_ref).norm()<5e-3);
        TEST_ASSERT((xpp-xpp_ref).norm()<5e-3);
      }
      else
      {
        TEST_ASSERT((xp-xp_ref).norm()<1e-10);
        TEST_ASSERT((xpp-xpp_ref).norm()<5e-5);
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
        TEST_ASSERT((xp-xp_ref).norm()<5e-3);
        TEST_ASSERT((xpp-xpp_ref).norm()<3e-2);
      }
      else
      {
        TEST_ASSERT((xp-xp_ref).norm()<1e-10);
        TEST_ASSERT((xpp-xpp_ref).norm()<6e-3);
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
      if (typeid(data_type)==typeid(float))
      {
        TEST_ASSERT((xp-xp_ref).norm()<1e-4);
        TEST_ASSERT((xpp-xpp_ref).norm()<5e-4);
      }
      else
      {
        TEST_ASSERT((xp-xp_ref).norm()<1e-10);
        TEST_ASSERT((xpp-xpp_ref).norm()<1e-4);
      }
//      std::cout << "xi=" << xi;
//      std::cout << "\tdx=" << (x[2]-x[0])/(2*dxi) << "\tdy=" << (y[2]-y[0])/(2*dxi);
//      std::cout << std::endl << "     ";
//      std::cout << "ddx=" << (x[0]-2*x[1]+x[2])/(dxi*dxi) << "\tddy=" << (y[0]-2*y[1]+y[2])/(dxi*dxi);
//      std::cout << std::endl << "     ";
//      std::cout << "xp=" << xp;
//      std::cout << std::endl << "     ";
//      std::cout << "xp_ref=" << xp_ref;
//      std::cout << std::endl << "     ";
//      std::cout << "dxp=" << xp-xp_ref << "\tdxp.norm()=" << (xp-xp_ref).norm();
//      std::cout << std::endl << "     ";
//      std::cout << "xpp=    " << xpp;
//      std::cout << std::endl << "     ";
//      std::cout << "xpp_ref=" << xpp_ref;
//      std::cout << std::endl << "     ";
//      std::cout << "dxpp=" << xpp-xpp_ref << "\tdxpp.norm()=" << (xpp-xpp_ref).norm();
//      std::cout << std::endl;
    }

};

#endif

