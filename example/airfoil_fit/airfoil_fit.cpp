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

#include <cstdlib>  // EXIT_SUCCESS
#include <cmath>    // std::sin, std::cos
#include <iostream> // std::iostream
#include <iomanip>
#include <string>
#include <vector>

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/bezier.hpp"

typedef eli::geom::curve::bezier2d bezier_type;
typedef bezier_type::data_type data_type;
typedef bezier_type::point_type point_type;
typedef bezier_type::fit_container_type fit_container_type;
typedef fit_container_type::constraint_info constraint_info_type;

void create_airfoil(std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts, int &le)
{
  pts.resize(40);
  // lower surface
  pts[ 0] << 1.0000,-0.0020;
  pts[ 1] << 0.9012,-0.0087;
  pts[ 2] << 0.8078,-0.0247;
  pts[ 3] << 0.7017,-0.0453;
  pts[ 4] << 0.6178,-0.0572;
  pts[ 5] << 0.4123,-0.0684;
  pts[ 6] << 0.3545,-0.0678;
  pts[ 7] << 0.2986,-0.0657;
  pts[ 8] << 0.2453,-0.0621;
  pts[ 9] << 0.1499,-0.0511;
  pts[10] << 0.1029,-0.0441;
  pts[11] << 0.0741,-0.0365;
  pts[12] << 0.0451,-0.0284;
  pts[13] << 0.0143,-0.0112;
  pts[14] << 0.0076,-0.0116;
  pts[15] << 0.0029,-0.0077;
  pts[16] << 0.0000, 0.0000;
  le=16;
  // upper surface
  pts[17] << 0.0013, 0.0030;
  pts[18] << 0.0045, 0.0094;
  pts[19] << 0.0096, 0.0138;
  pts[20] << 0.0165, 0.0183;
  pts[21] << 0.0252, 0.0228;
  pts[22] << 0.0477, 0.0315;
  pts[23] << 0.0767, 0.0397;
  pts[24] << 0.1118, 0.0473;
  pts[25] << 0.1524, 0.0539;
  pts[26] << 0.1980, 0.0593;
  pts[27] << 0.2979, 0.0636;
  pts[28] << 0.3015, 0.0665;
  pts[29] << 0.3578, 0.0680;
  pts[30] << 0.4160, 0.0680;
  pts[31] << 0.4455, 0.0675;
  pts[32] << 0.5049, 0.0652;
  pts[33] << 0.5930, 0.0585;
  pts[34] << 0.6501, 0.0514;
  pts[35] << 0.7050, 0.0416;
  pts[36] << 0.7623, 0.0297;
  pts[37] << 0.8168, 0.0221;
  pts[38] << 0.9074, 0.0108;
  pts[39] << 1.0000, 0.0050;
}

void create_circle(std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts)
{
  // NOTE: This will not create a closed circle, the last point will be
  //       the point just prior to the 2*pi point
  pts.resize(20);

  for (size_t i=0; i<pts.size(); ++i)
  {
    data_type theta(eli::constants::math<data_type>::two_pi()*static_cast<data_type>(i)/pts.size());
    pts[i](0)=std::cos(theta);
    pts[i](1)=std::sin(theta);
  }
}

void octave_print(std::ostream &ostr, const std::string &varname, const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts)
{
  size_t i;

  ostr << varname << "_x=[" << pts[0](0);
  for (i=1; i<pts.size(); ++i)
    ostr << ", " << pts[i](0);
  ostr << "];" << std::endl;
  ostr << varname << "_y=[" << pts[0](1);
  for (i=1; i<pts.size(); ++i)
    ostr << ", " << pts[i](1);
  ostr << "];" << std::endl;
}

void octave_print(std::ostream &ostr, const std::string &varname, const bezier_type &bez)
{
  size_t i;

  std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(101);
  for (i=0; i<pts.size(); ++i)
    pts[i]=bez.f(static_cast<double>(i)/(pts.size()-1));

  octave_print(ostr, varname, pts);
}

data_type example_fit_1(bezier_type &bez, const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts, const size_t dim)
{
  // set up fit container
  fit_container_type fcon;
  fcon.set_points(pts.begin(), pts.end());

  // do fit
  return bez.fit(fcon, dim);
}

data_type example_fit_2(bezier_type &bez, const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts, const size_t dim)
{
  // set up fit container
  fit_container_type fcon;
  fcon.set_points(pts.begin(), pts.end());
  fcon.add_start_C0_constraint();
  fcon.add_end_C0_constraint();

  // do fit
  return bez.fit(fcon, dim);
}

data_type example_fit_3(bezier_type &bez, constraint_info_type ci[2], const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts, const size_t dim)
{
  // set up fit container
  fit_container_type fcon;
  fcon.set_points(pts.begin(), pts.end());
  fcon.add_start_C1_constraint();
  fcon.add_end_C1_constraint();
  fcon.get_start_constraint(ci[0]);
  fcon.get_end_constraint(ci[1]);

  // do fit
  return bez.fit(fcon, dim);
}

data_type example_fit_4(bezier_type &bez, constraint_info_type ci[2], const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts, const size_t dim, const int &index)
{
  // set up fit container
  fit_container_type fcon;
  fcon.set_points(pts.begin(), pts.end());
  fcon.add_start_C1_constraint();
  fcon.add_end_C1_constraint();
  fcon.add_C0_constraint(index);
  fcon.get_start_constraint(ci[0]);
  fcon.get_end_constraint(ci[1]);

  // do fit
  return bez.fit(fcon, dim);
}

data_type example_fit_5(bezier_type &bez, constraint_info_type ci[2], const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts, const size_t dim, const int index[3])
{
  // set up fit container
  fit_container_type fcon;
  fcon.set_points(pts.begin(), pts.end());
  fcon.add_start_C1_constraint();
  fcon.add_end_C1_constraint();
  fcon.add_C0_constraint(index[0]);
  fcon.add_C0_constraint(index[1]);
  fcon.add_C0_constraint(index[2]);
  fcon.get_start_constraint(ci[0]);
  fcon.get_end_constraint(ci[1]);

  // do fit
  return bez.fit(fcon, dim);
}

void example_fit_6(bezier_type bez[4], data_type err[4], bezier_type::fit_container_type::constraint_info ci[6], const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts, const bezier_type::dimension_type dim[4])
{
  size_t i;

  // default fit
  i=0;
  {
    // set up fit container
    fit_container_type fcon;
    fcon.set_points(pts.begin(), pts.end());

    // do fit
    err[i]=bez[i].fit(fcon, dim[i]);
  }

  // default fit
  i=1;
  {
    // set up fit container
    fit_container_type fcon;
    fcon.set_points(pts.begin(), pts.end());
    fcon.set_end_flag(eli::geom::general::C2);

    // do fit
    err[i]=bez[i].fit(fcon, dim[i]);
  }

  // C2-closed fit
  i=2;
  {
    // set up fit container
    fit_container_type fcon;
    fcon.set_points(pts.begin(), pts.end());
    fcon.set_end_flag(eli::geom::general::C2);

    // do fit
    err[i]=bez[i].fit(fcon, dim[i]);
  }

  // C2-closed fit with C1 at start specified
  i=3;
  {
    // set the first and second derivative values
    point_type fp0;
    fp0(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
    fp0(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);

    // set up fit container
    fit_container_type fcon;
    fcon.set_points(pts.begin(), pts.end());
    fcon.set_end_flag(eli::geom::general::C2);
    fcon.add_start_C1_constraint(fp0);
    fcon.get_constraint(0, ci[i]); // constraint on point 0

    // do fit
    err[i]=bez[i].fit(fcon, dim[i]);
  }

  // C2-closed fit with C2 at start specified
  i=4;
  {
    // set the first and second derivative values
    point_type fp0, fpp0;
    fp0(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
    fp0(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
    fpp0(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](0);
    fpp0(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](1);

    // set up fit container
    fit_container_type fcon;
    fcon.set_points(pts.begin(), pts.end());
    fcon.set_end_flag(eli::geom::general::C2);
    fcon.add_start_C2_constraint(fp0, fpp0);
    fcon.get_constraint(0, ci[i]); // constraint on point 0

    // do fit
    err[i]=bez[i].fit(fcon, dim[i]);
  }

  // C2-closed fit with C2 at start specified
  i=5;
  {
    // set the first and second derivative values
    point_type fp0, fpp0;
    fp0(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
    fp0(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
    fpp0(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](0);
    fpp0(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](1);

    // set up fit container
    fit_container_type fcon;
    fcon.set_points(pts.begin(), pts.end());
    fcon.set_end_flag(eli::geom::general::C2);
    fcon.add_start_C2_constraint(fp0, fpp0);
    fcon.get_constraint(0, ci[i]); // constraint on point 0

    // do fit
    err[i]=bez[i].fit(fcon, dim[i]);
  }

}

int main(int /*argc*/, char * /*argv*/[])
{
  std::cout << "figure(1);" << std::endl;
  // Use Case 1: Simple fit with no constraints
  //
  //             This would correspond to the most basic fitting that a user would perform.
  //             They would input a set of points and then fit a bezier curve to the points.
  //             The use would have to specify the degree for the bezier curve.
  {
    int le;
    bezier_type::dimension_type d(8);
    bezier_type bez;
    std::vector<point_type, Eigen::aligned_allocator<point_type> > pts;
    data_type err;

    // get points to print out
    create_airfoil(pts, le);

    // get the bezier curve from the fit
    err=example_fit_1(bez, pts, d);

    std::cout << "subplot(5, 1, 1);" << std::endl;
    octave_print(std::cout, "af", pts);
    octave_print(std::cout, "free", bez);
    std::cout << "plot(af_x, af_y, 'ok', free_x, free_y, '-k');" << std::endl;
    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
    std::cout << "title('" << d << "th Order Bezier - Free Fit');" << std::endl;
    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;
  }

  // Use Case 2: Fit with fixed end points
  //
  //             This corresponds to a case were the user wants to fit the points and ensure
  //             that the curve goes through the first and last points.
  {
    int le;
    bezier_type::dimension_type d(8);
    bezier_type bez;
    std::vector<point_type, Eigen::aligned_allocator<point_type> > pts;
    data_type err;

    // get points to print out
    create_airfoil(pts, le);

    // get the bezier curve from the fit
    err=example_fit_2(bez, pts, d);

    std::cout << "subplot(5, 1, 2);" << std::endl;
    octave_print(std::cout, "af", pts);
    octave_print(std::cout, "free", bez);
    std::cout << "plot(af_x, af_y, 'ok', free_x, free_y, '-k',";
    std::cout << " [" << pts[0](0) << "," << pts[pts.size()-1](0) << "],";
    std::cout << " [" << pts[0](1) << "," << pts[pts.size()-1](1) << "], 'or');" << std::endl;
    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
    std::cout << "title('" << d << "th Order Bezier - Fit C0 Ends');" << std::endl;
    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;
  }

  // Use Case 3: Fixed end points and slope
  //
  //             This corresponds to the use case where the user wants to ensure that the
  //             curve goes through the first and last points and that the slope of the curve
  //             at those two points without knowing the slope. A 2nd order finite difference
  //             is used of the 1st derivative vector is not specified.
  {
    int le;
    bezier_type::dimension_type d(8);
    bezier_type::fit_container_type::constraint_info ci[2];
    bezier_type bez;
    std::vector<point_type, Eigen::aligned_allocator<point_type> > pts;
    data_type err;

    // get points to print out
    create_airfoil(pts, le);

    // get the bezier curve from the fit
    err=example_fit_3(bez, ci, pts, d);

    std::cout << "subplot(5, 1, 3);" << std::endl;
    octave_print(std::cout, "af", pts);
    octave_print(std::cout, "free", bez);
    std::cout << "plot(af_x, af_y, 'ok', free_x, free_y, '-k',";
    std::cout << " [" << pts[0](0) << "," << pts[pts.size()-1](0) << "],";
    std::cout << " [" << pts[0](1) << "," << pts[pts.size()-1](1) << "], 'or');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "quiver([" << pts[0](0) << ", " << pts[pts.size()-1](0) << "], [" << pts[0](1) << ", " << pts[pts.size()-1](1) << "], 0.05*["
                            << ci[0].get_fp()(0) << ", " << ci[1].get_fp()(0) << "], 0.05*["<< ci[0].get_fp()(1) << ", " << ci[1].get_fp()(1) << "], 0, 'r');" << std::endl;
    std::cout << "hold off;" << std::endl;
    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
    std::cout << "title('" << d << "th Order Bezier - Fit C1 Ends');" << std::endl;
    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;
    std::cout << "text(1.03, 0.015, '\\color{red}dr/dt');" << std::endl;
    std::cout << "text(0.95, -0.015, '\\color{red}dr/dt');" << std::endl;
  }

  // Use Case 4: Fixed end points and fixed leading edge
  //
  //             This corresponds to the use case where the user picks a point and wants to
  //             improve the quality of the fit by ensuring that the curve goes through that
  //             point. For this airfoil the leading edge goes through point 16.
  {
    bezier_type::dimension_type d(8);
    bezier_type::fit_container_type::constraint_info ci[2];
    int le;
    bezier_type bez;
    std::vector<point_type, Eigen::aligned_allocator<point_type> > pts;
    data_type err;

    // get points to print out
    create_airfoil(pts, le);

    // get the bezier curve from the fit
    err=example_fit_4(bez, ci, pts, d, le);

    std::cout << "subplot(5, 1, 4);" << std::endl;
    octave_print(std::cout, "af", pts);
    octave_print(std::cout, "free", bez);
    std::cout << "plot(af_x, af_y, 'ok', free_x, free_y, '-k',";
    std::cout << " [" << pts[0](0) << "," << pts[pts.size()-1](0) << "," << pts[le](0) << "],";
    std::cout << " [" << pts[0](1) << "," << pts[pts.size()-1](1) << "," << pts[le](1)  << "], 'or');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "quiver([" << pts[0](0) << ", " << pts[pts.size()-1](0) << "], [" << pts[0](1) << ", " << pts[pts.size()-1](1) << "], 0.05*["
                            << ci[0].get_fp()(0) << ", " << ci[1].get_fp()(0) << "], 0.05*["<< ci[0].get_fp()(1) << ", " << ci[1].get_fp()(1) << "], 0, 'r');" << std::endl;
    std::cout << "hold off;" << std::endl;
    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
    std::cout << "title('" << d << "th Order Bezier - Fit C1 Ends, C0 L.E.');" << std::endl;
    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;
    std::cout << "text(1.03, 0.015, '\\color{red}dr/dt');" << std::endl;
    std::cout << "text(0.95, -0.015, '\\color{red}dr/dt');" << std::endl;
  }

  // Use Case 5: Fixed end points, fixed leading edge and two other points on airfoil
  //
  //             This corresponds to the use case where the user picks a point and wants to
  //             further improve the quality of the fit by ensuring that the curve goes through
  //             several points. For this airfoil the leading edge goes through point 16.
  //             The location of maximum thickness is 5 on lower surface and 30 on upper
  //             surface.
  {
    bezier_type::dimension_type d(8);
    constraint_info_type ci[3];
    int le, index[3];
    bezier_type bez;
    std::vector<point_type, Eigen::aligned_allocator<point_type> > pts;
    data_type err;

    // get points to print out
    create_airfoil(pts, le);

    // get the bezier curve from the fit
    index[0]=le;
    index[1]=5;
    index[2]=30;
    err=example_fit_5(bez, ci, pts, d, index);

    std::cout << "subplot(5, 1, 5);" << std::endl;
    octave_print(std::cout, "af", pts);
    octave_print(std::cout, "free", bez);
    std::cout << "plot(af_x, af_y, 'ok', free_x, free_y, '-k',";
    std::cout << " [" << pts[0](0) << "," << pts[pts.size()-1](0) << "," << pts[index[0]](0) << "," << pts[index[1]](0) << "," << pts[index[2]](0) << "],";
    std::cout << " [" << pts[0](1) << "," << pts[pts.size()-1](1) << "," << pts[index[0]](1) << "," << pts[index[1]](1) << "," << pts[index[2]](1) << "], 'or');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "quiver([" << pts[0](0) << ", " << pts[pts.size()-1](0) << "], [" << pts[0](1) << ", " << pts[pts.size()-1](1) << "], 0.05*["
                            << ci[0].get_fp()(0) << ", " << ci[1].get_fp()(0) << "], 0.05*["<< ci[0].get_fp()(1) << ", " << ci[1].get_fp()(1) << "], 0, 'r');" << std::endl;
    std::cout << "hold off;" << std::endl;
    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
    std::cout << "title('" << d << "th Order Bezier - Fit C1 Ends, C0 L.E. & 2 Other Points');" << std::endl;
    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;
    std::cout << "text(1.03, 0.015, '\\color{red}dr/dt');" << std::endl;
    std::cout << "text(0.95, -0.015, '\\color{red}dr/dt');" << std::endl;
  }

  std::cout << "figure(2);" << std::endl;
  // Use Case 6: Fit curve using known derivative information and a closed curve
  //
  //             This case creates several fits to a circle which cannot be exactly represented
  //             by a bezier curve. Several different constraints are used to improve the
  //             quality of the fit.
  {
    bezier_type::dimension_type d[6]={5, 5, 7, 7, 9, 12};
    bezier_type::fit_container_type::constraint_info ci[6];
    bezier_type bez[6];
    std::vector<point_type, Eigen::aligned_allocator<point_type> > pts;
    data_type err[6];

    // get points to print out
    create_circle(pts);

    // get the bezier curve from the fit
    example_fit_6(bez, err, ci, pts, d);

    std::cout << "subplot(3, 2, 1);" << std::endl;
    octave_print(std::cout, "cir", pts);
    octave_print(std::cout, "bez", bez[0]);
    std::cout << "plot(cir_x, cir_y, 'ok', bez_x, bez_y, '-k');" << std::endl;
    std::cout << "axis([-1.2,1.2,-1.2,1.2]);" << std::endl;
    std::cout << "title('" << d[0] << "th Order Bezier - Free Fit');" << std::endl;
    std::cout << "text(-0.3, 0.2, '\\Sigma d^2=" << err[0] << "');" << std::endl;

    std::cout << "subplot(3, 2, 2);" << std::endl;
    octave_print(std::cout, "cir", pts);
    octave_print(std::cout, "bez", bez[1]);
    std::cout << "plot(cir_x, cir_y, 'ok', bez_x, bez_y, '-k');" << std::endl;
    std::cout << "axis([-1.2,1.2,-1.2,1.2]);" << std::endl;
    std::cout << "title('" << d[1] << "th Order Bezier - C2 Closed Fit');" << std::endl;
    std::cout << "text(-0.3, 0.2, '\\Sigma d^2=" << err[1] << "');" << std::endl;

    std::cout << "subplot(3, 2, 3);" << std::endl;
    octave_print(std::cout, "cir", pts);
    octave_print(std::cout, "bez", bez[2]);
    std::cout << "plot(cir_x, cir_y, 'ok', bez_x, bez_y, '-k');" << std::endl;
    std::cout << "axis([-1.2,1.2,-1.2,1.2]);" << std::endl;
    std::cout << "title('" << d[2] << "th Order Bezier - C2 Closed Fit');" << std::endl;
    std::cout << "text(-0.3, 0.2, '\\Sigma d^2=" << err[2] << "');" << std::endl;

    std::cout << "subplot(3, 2, 4);" << std::endl;
    octave_print(std::cout, "cir", pts);
    octave_print(std::cout, "bez", bez[3]);
    std::cout << "plot(cir_x, cir_y, 'ok', bez_x, bez_y, '-k',";
    std::cout << " [" << pts[0](0) << "],";
    std::cout << " [" << pts[0](1) << "], 'or');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "quiver([" << pts[0](0) << "], ["<< pts[0](1) << "], [" << ci[3].get_fp()(0)/10 << "], ["<< ci[3].get_fp()(1)/10 << "], 'r');" << std::endl;
    std::cout << "hold off;" << std::endl;
    std::cout << "axis([-1.2,1.2,-1.2,1.2]);" << std::endl;
    std::cout << "title('" << d[3] << "th Order Bezier - C2 Closed Fit with C1 Set at (1,0)');" << std::endl;
    std::cout << "text(-0.3, 0.2, '\\Sigma d^2=" << err[3] << "');" << std::endl;
    std::cout << "text(0.9, 0.7, '\\color{red}dr/dt');" << std::endl;

    std::cout << "subplot(3, 2, 5);" << std::endl;
    octave_print(std::cout, "cir", pts);
    octave_print(std::cout, "bez", bez[4]);
    std::cout << "plot(cir_x, cir_y, 'ok', bez_x, bez_y, '-k',";
    std::cout << " [" << pts[0](0) << "],";
    std::cout << " [" << pts[0](1) << "], 'or');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "quiver([" << pts[0](0) << ", " << pts[0](0) << "], [" << pts[0](1) << ", " << pts[0](1) << "], ["
                            << ci[4].get_fp()(0)/10 << ", " << ci[4].get_fpp()(0)/80 << "], ["<< ci[4].get_fp()(1)/10 << ", " << ci[4].get_fpp()(1)/80 << "], 'r');" << std::endl;
    std::cout << "hold off;" << std::endl;
    std::cout << "axis([-1.2,1.2,-1.2,1.2]);" << std::endl;
    std::cout << "title('" << d[4] << "th Order Bezier - C2 Closed Fit with C2 Set at (1,0)');" << std::endl;
    std::cout << "text(-0.3, 0.2, '\\Sigma d^2=" << err[4] << "');" << std::endl;
    std::cout << "text(0.9, 0.7, '\\color{red}dr/dt');" << std::endl;
    std::cout << "text(0.4, -0.1, '\\color{red}d^2r/dt^2');" << std::endl;

    std::cout << "subplot(3, 2, 6);" << std::endl;
    octave_print(std::cout, "cir", pts);
    octave_print(std::cout, "bez", bez[5]);
    std::cout << "plot(cir_x, cir_y, 'ok', bez_x, bez_y, '-k',";
    std::cout << " [" << pts[0](0) << "],";
    std::cout << " [" << pts[0](1) << "], 'or');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "quiver([" << pts[0](0) << ", " << pts[0](0) << "], [" << pts[0](1) << ", " << pts[0](1) << "], ["
                            << ci[5].get_fp()(0)/10 << ", " << ci[5].get_fpp()(0)/80 << "], ["<< ci[5].get_fp()(1)/10 << ", " << ci[5].get_fpp()(1)/80 << "], 'r');" << std::endl;
    std::cout << "hold off;" << std::endl;
    std::cout << "axis([-1.2,1.2,-1.2,1.2]);" << std::endl;
    std::cout << "title('" << d[5] << "th Order Bezier - C2 Closed Fit with C2 Set at (1,0)');" << std::endl;
    std::cout << "text(-0.3, 0.2, '\\Sigma d^2=" << err[5] << "');" << std::endl;
    std::cout << "text(0.9, 0.7, '\\color{red}dr/dt');" << std::endl;
    std::cout << "text(0.4, -0.1, '\\color{red}d^2r/dt^2');" << std::endl;
  }

  std::cout << "figure(3);" << std::endl;
  {
    int le;
    bezier_type::dimension_type d(8);
    bezier_type bez;
    std::vector<point_type, Eigen::aligned_allocator<point_type> > pts;
    std::vector<data_type> t;
    data_type err;
    std::vector<point_type, Eigen::aligned_allocator<point_type> > af_pts;
    bezier_type::fit_container_type::constraint_info ci[2];
    // get points to print out
    create_airfoil(pts, le);

    // set up fit container
    {
      fit_container_type fcon;
      fcon.set_points(pts.begin(), pts.end());

      // do fit
      err=bez.fit(t, fcon, d);
    }

    af_pts.resize(t.size());

    for (size_t i=0; i<t.size(); ++i)
      af_pts[i]=bez.f(t[i]);

    std::cout << "subplot(5, 1, 1);" << std::endl;
    octave_print(std::cout, "af", pts);
    octave_print(std::cout, "free", bez);
    octave_print(std::cout, "af_pts", af_pts);
    std::cout << "plot(af_x, af_y, 'ob', free_x, free_y, '-k', af_pts_x, af_pts_y, 'ok');" << std::endl;
    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
    std::cout << "title('" << d << "th Order Bezier - Free Fit');" << std::endl;
    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;

    // set up fit container
    {
      fit_container_type fcon;
      fcon.set_points(pts.begin(), pts.end());
      fcon.add_start_C0_constraint();
      fcon.add_end_C0_constraint();

      // do fit
      err=bez.fit(fcon, d);
    }

    af_pts.resize(t.size());

    for (size_t i=0; i<t.size(); ++i)
      af_pts[i]=bez.f(t[i]);

    std::cout << "subplot(5, 1, 2);" << std::endl;
    octave_print(std::cout, "af", pts);
    octave_print(std::cout, "free", bez);
    octave_print(std::cout, "af_pts", af_pts);
    std::cout << "plot(af_x, af_y, 'ob', free_x, free_y, '-k', af_pts_x, af_pts_y, 'ok',";
    std::cout << " [" << pts[0](0) << "," << pts[pts.size()-1](0) << "],";
    std::cout << " [" << pts[0](1) << "," << pts[pts.size()-1](1) << "], 'or');" << std::endl;
    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
    std::cout << "title('" << d << "th Order Bezier - Fit C0 Ends');" << std::endl;
    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;

    // set up fit container
    {
      fit_container_type fcon;
      fcon.set_points(pts.begin(), pts.end());
      fcon.add_start_C1_constraint();
      fcon.add_end_C1_constraint();
      fcon.get_start_constraint(ci[0]);
      fcon.get_end_constraint(ci[1]);

      // do fit
      err=bez.fit(fcon, d);
    }

    af_pts.resize(t.size());

    for (size_t i=0; i<t.size(); ++i)
      af_pts[i]=bez.f(t[i]);

    std::cout << "subplot(5, 1, 3);" << std::endl;
    octave_print(std::cout, "af", pts);
    octave_print(std::cout, "free", bez);
    octave_print(std::cout, "af_pts", af_pts);
    std::cout << "plot(af_x, af_y, 'ob', free_x, free_y, '-k', af_pts_x, af_pts_y, 'ok',";
    std::cout << " [" << pts[0](0) << "," << pts[pts.size()-1](0) << "],";
    std::cout << " [" << pts[0](1) << "," << pts[pts.size()-1](1) << "], 'or');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "quiver([" << pts[0](0) << ", " << pts[pts.size()-1](0) << "], [" << pts[0](1) << ", " << pts[pts.size()-1](1) << "], 0.05*["
                            << ci[0].get_fp()(0) << ", " << ci[1].get_fp()(0) << "], 0.05*["<< ci[0].get_fp()(1) << ", " << ci[1].get_fp()(1) << "], 0, 'r');" << std::endl;
    std::cout << "hold off;" << std::endl;
    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
    std::cout << "title('" << d << "th Order Bezier - Fit C1 Ends');" << std::endl;
    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;
    std::cout << "text(1.03, 0.015, '\\color{red}dr/dt');" << std::endl;
    std::cout << "text(0.95, -0.015, '\\color{red}dr/dt');" << std::endl;

    // set up fit container
    {
      fit_container_type fcon;
      fcon.set_points(pts.begin(), pts.end());
      fcon.add_start_C1_constraint();
      fcon.add_end_C1_constraint();
      fcon.add_C0_constraint(le);
      fcon.get_start_constraint(ci[0]);
      fcon.get_end_constraint(ci[1]);

      // do fit
      err=bez.fit(fcon, d);
    }

    af_pts.resize(t.size());

    for (size_t i=0; i<t.size(); ++i)
      af_pts[i]=bez.f(t[i]);

    std::cout << "subplot(5, 1, 4);" << std::endl;
    octave_print(std::cout, "af", pts);
    octave_print(std::cout, "free", bez);
    octave_print(std::cout, "af_pts", af_pts);
    std::cout << "plot(af_x, af_y, 'ob', free_x, free_y, '-k', af_pts_x, af_pts_y, 'ok',";
    std::cout << " [" << pts[0](0) << "," << pts[pts.size()-1](0) << "," << pts[le](0) << "],";
    std::cout << " [" << pts[0](1) << "," << pts[pts.size()-1](1) << "," << pts[le](1)  << "], 'or');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "quiver([" << pts[0](0) << ", " << pts[pts.size()-1](0) << "], [" << pts[0](1) << ", " << pts[pts.size()-1](1) << "], 0.05*["
                            << ci[0].get_fp()(0) << ", " << ci[1].get_fp()(0) << "], 0.05*["<< ci[0].get_fp()(1) << ", " << ci[1].get_fp()(1) << "], 0, 'r');" << std::endl;
    std::cout << "hold off;" << std::endl;
    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
    std::cout << "title('" << d << "th Order Bezier - Fit C1 Ends, C0 L.E.');" << std::endl;
    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;
    std::cout << "text(1.03, 0.015, '\\color{red}dr/dt');" << std::endl;
    std::cout << "text(0.95, -0.015, '\\color{red}dr/dt');" << std::endl;

    // set up fit container
    {
      fit_container_type fcon;
      fcon.set_points(pts.begin(), pts.end());
      fcon.add_start_C1_constraint();
      fcon.add_end_C1_constraint();
      fcon.add_C0_constraint(le);
      fcon.add_C0_constraint(5);
      fcon.add_C0_constraint(30);
      fcon.get_start_constraint(ci[0]);
      fcon.get_end_constraint(ci[1]);

      // do fit
      err=bez.fit(fcon, d);
    }

    af_pts.resize(t.size());

    for (size_t i=0; i<t.size(); ++i)
      af_pts[i]=bez.f(t[i]);

    std::cout << "subplot(5, 1, 5);" << std::endl;
    octave_print(std::cout, "af", pts);
    octave_print(std::cout, "free", bez);
    octave_print(std::cout, "af_pts", af_pts);
    std::cout << "plot(af_x, af_y, 'ob', free_x, free_y, '-k', af_pts_x, af_pts_y, 'ok',";
    std::cout << " [" << pts[0](0) << "," << pts[pts.size()-1](0) << "," << pts[le](0) << "," << pts[5](0) << "," << pts[30](0) << "],";
    std::cout << " [" << pts[0](1) << "," << pts[pts.size()-1](1) << "," << pts[le](1) << "," << pts[5](1) << "," << pts[30](1) << "], 'or');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "quiver([" << pts[0](0) << ", " << pts[pts.size()-1](0) << "], [" << pts[0](1) << ", " << pts[pts.size()-1](1) << "], 0.05*["
                            << ci[0].get_fp()(0) << ", " << ci[1].get_fp()(0) << "], 0.05*["<< ci[0].get_fp()(1) << ", " << ci[1].get_fp()(1) << "], 0, 'r');" << std::endl;
    std::cout << "hold off;" << std::endl;
    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
    std::cout << "title('" << d << "th Order Bezier - Fit C1 Ends, C0 L.E. & 2 Other Points');" << std::endl;
    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;
    std::cout << "text(1.03, 0.015, '\\color{red}dr/dt');" << std::endl;
    std::cout << "text(0.95, -0.015, '\\color{red}dr/dt');" << std::endl;
  }

  return EXIT_SUCCESS;
}

