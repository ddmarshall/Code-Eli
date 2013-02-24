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
#include "eli/geom/surface/bezier.hpp"

typedef eli::geom::curve::bezier3d bezier_curve_type;
typedef bezier_curve_type::control_point_type curve_control_point_type;
typedef eli::geom::surface::bezier3d bezier_surface_type;
typedef bezier_curve_type::data_type data_type;
typedef bezier_curve_type::point_type point_type;
typedef bezier_curve_type::index_type index_type;

void octave_print(std::ostream &ostr, const std::string &varname, const std::vector<point_type, Eigen::aligned_allocator<point_type> > &pts)
{
  index_type i;

  ostr << varname << "_x=[" << pts[0](0);
  for (i=1; i<static_cast<index_type>(pts.size()); ++i)
    ostr << ", " << pts[i](0);
  ostr << "];" << std::endl;
  ostr << varname << "_y=[" << pts[0](1);
  for (i=1; i<static_cast<index_type>(pts.size()); ++i)
    ostr << ", " << pts[i](1);
  ostr << "];" << std::endl;
  ostr << varname << "_z=[" << pts[0](2);
  for (i=1; i<static_cast<index_type>(pts.size()); ++i)
    ostr << ", " << pts[i](2);
  ostr << "];" << std::endl;
}

void octave_print(std::ostream &ostr, const std::string &varname, const bezier_curve_type &bez)
{
  index_type i;

  std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(101);
  for (i=0; i<static_cast<index_type>(pts.size()); ++i)
    pts[i]=bez.f(static_cast<data_type>(i)/(pts.size()-1));

  octave_print(ostr, varname, pts);
}

void octave_print(std::ostream &ostr, const std::string &varname, const bezier_surface_type &bez)
{
  index_type i, j;

    // set the u and v sample values
    std::vector<data_type> u(51), v(51);
    for (i=0; i<static_cast<index_type>(u.size()); ++i)
    {
      u[i]=static_cast<data_type>(i)/(u.size()-1);
    }
    for (j=0; j<static_cast<index_type>(v.size()); ++j)
    {
      v[j]=static_cast<data_type>(j)/(v.size()-1);
    }

    // set the surface points
    ostr << varname << "_x=[";
    for (i=0; i<static_cast<index_type>(u.size()); ++i)
    {
      ostr << bez.f(u[i], v[0]).x();
      for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
      {
        ostr << ", " << bez.f(u[i], v[j]).x();
      }
      j=static_cast<index_type>(v.size()-1);
      ostr << ", " << bez.f(u[i], v[j]).x();
      if (i<static_cast<index_type>(u.size()-1))
        ostr << "; " << std::endl;
    }
    ostr << "];" << std::endl;

    ostr << varname << "_y=[";
    for (i=0; i<static_cast<index_type>(u.size()); ++i)
    {
      ostr << bez.f(u[i], v[0]).y();
      for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
      {
        ostr << ", " << bez.f(u[i], v[j]).y();
      }
      j=static_cast<index_type>(v.size()-1);
      ostr << ", " << bez.f(u[i], v[j]).y();
      if (i<static_cast<index_type>(u.size()-1))
        ostr << "; " << std::endl;
    }
    ostr << "];" << std::endl;

    ostr << varname << "_z=[";
    for (i=0; i<static_cast<index_type>(u.size()); ++i)
    {
      ostr << bez.f(u[i], v[0]).z();
      for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
      {
        ostr << ", " << bez.f(u[i], v[j]).z();
      }
      j=static_cast<index_type>(v.size()-1);
      ostr << ", " << bez.f(u[i], v[j]).z();
      if (i<static_cast<index_type>(u.size()-1))
        ostr << "; " << std::endl;
    }
    ostr << "];" << std::endl;
}

int main(int /*argc*/, char * /*argv*/[])
{
  data_type fine_ratio;
  data_type len;

  //////////////
  // these are the user inputs
  //////////////
  fine_ratio=15;
  len=10;

  //////////////
  // build the curve through the points
  //////////////
  double dia;
  curve_control_point_type cp(7,3);
  bezier_curve_type bez_c, bez_c2;

  dia = len/fine_ratio;
  // FIX: THIS SHOULD BE TWO 3rd ORDER BEZIER CURVES
  cp << 0.00*len, 0, 0.00*dia,
        0.05*len, 0, 0.95*dia,
        0.20*len, 0, 1.00*dia,
        0.50*len, 0, 1.00*dia,
        0.60*len, 0, 1.00*dia,
        0.95*len, 0, 0.30*dia,
        1.00*len, 0, 0.00*dia;

  // create top reference curve
  bez_c.set_control_points(cp);

  // create bottom reference curve
  cp.col(2)*=-1;
  bez_c2.set_control_points(cp);
  cp.col(2)*=-1;

  // create bezier surface by rotating the control points of curve
  // FIX: THIS NEEDS TO BE A PIECEWISE
  index_type i, j, n(6), m(9);
  bezier_surface_type bez_s(n, m);
  point_type tmp;

  for (j=0; j<=m; ++j)
  {
    data_type alpha=(static_cast<data_type>(j)/m)*eli::constants::math<data_type>::two_pi();
    for (i=0; i<=n; ++i)
    {
      tmp(0)=cp(i, 0);
      tmp(1)= std::cos(alpha)*cp(i,1)+std::sin(alpha)*cp(i,2);
      tmp(2)=-std::sin(alpha)*cp(i,1)+std::cos(alpha)*cp(i,2);
      bez_s.set_control_point(tmp, i, j);
    }
  }

  // print out octave code
  octave_print(std::cout, "bez_c", bez_c);
  octave_print(std::cout, "bez_c2", bez_c2);
  octave_print(std::cout, "bez_s", bez_s);
  std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
  std::cout << "mesh(bez_s_x, bez_s_y, bez_s_z, zeros(size(bez_s_z)), 'EdgeColor', [0 0 0]);" << std::endl;
  std::cout << "hold on;" << std::endl;
  std::cout << "plot3(bez_c_x, bez_c_y, bez_c_z, '-k');" << std::endl;
  std::cout << "plot3(bez_c2_x, bez_c2_y, bez_c2_z, '-k');" << std::endl;
  std::cout << "hold off;" << std::endl;
#if 0
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
#endif

  return EXIT_SUCCESS;
}

