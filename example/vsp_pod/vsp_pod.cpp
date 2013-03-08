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
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_creator.hpp"

#include "eli/geom/surface/bezier.hpp"
#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/piecewise_creator.hpp"

typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, double, 3> piecewise_curve_type;
typedef piecewise_curve_type::curve_type curve_type;
typedef curve_type::control_point_type curve_control_point_type;

typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, double, 3> piecewise_surface_type;
typedef piecewise_surface_type::surface_type surface_type;
typedef surface_type::control_point_type surface_control_point_type;

typedef curve_type::point_type point_type;
typedef curve_type::data_type data_type;
typedef curve_type::index_type index_type;

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

void octave_print(std::ostream &ostr, const std::string &varname, const piecewise_curve_type &pc)
{
  index_type i;
  data_type tmin, tmax;

  pc.get_parameter_min(tmin);
  pc.get_parameter_max(tmax);

  // initialize the t parameters
  std::vector<data_type> t(101);
  for (i=0; i<static_cast<index_type>(t.size()); ++i)
  {
    t[i]=tmin+(tmax-tmin)*static_cast<data_type>(i)/(t.size()-1);
  }

  std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(t.size());
  for (i=0; i<static_cast<index_type>(pts.size()); ++i)
    pts[i]=pc.f(t[i]);

  octave_print(ostr, varname, pts);
}

void octave_print(std::ostream &ostr, const std::string &varname, const curve_type &bez)
{
  index_type i;

  std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(101);
  for (i=0; i<static_cast<index_type>(pts.size()); ++i)
    pts[i]=bez.f(static_cast<data_type>(i)/(pts.size()-1));

  octave_print(ostr, varname, pts);
}

void octave_print(std::ostream &ostr, const std::string &varname, const piecewise_surface_type &ps)
{
  index_type i, j;
  data_type umin, vmin, umax, vmax;

  ps.get_parameter_min(umin, vmin);
  ps.get_parameter_max(umax, vmax);

  // initialize the u & v parameters
  std::vector<data_type> u(21), v(21);
  for (i=0; i<static_cast<index_type>(u.size()); ++i)
  {
    u[i]=umin+(umax-umin)*static_cast<data_type>(i)/(u.size()-1);
  }
  for (j=0; j<static_cast<index_type>(v.size()); ++j)
  {
    v[j]=vmin+(vmax-vmin)*static_cast<data_type>(j)/(v.size()-1);
  }

  // set the surface points
  ostr << varname << "_x=[";
  for (i=0; i<static_cast<index_type>(u.size()); ++i)
  {
    ostr << ps.f(u[i], v[0]).x();
    for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
    {
      ostr << ", " << ps.f(u[i], v[j]).x();
    }
    j=static_cast<index_type>(v.size()-1);
    ostr << ", " << ps.f(u[i], v[j]).x();
    if (i<static_cast<index_type>(u.size()-1))
      std::cout << "; " << std::endl;
  }
  ostr << "];" << std::endl;

  ostr << varname << "_y=[";
  for (i=0; i<static_cast<index_type>(u.size()); ++i)
  {
    ostr << ps.f(u[i], v[0]).y();
    for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
    {
      ostr << ", " << ps.f(u[i], v[j]).y();
    }
    j=static_cast<index_type>(v.size()-1);
    ostr << ", " << ps.f(u[i], v[j]).y();
    if (i<static_cast<index_type>(u.size()-1))
      ostr << "; " << std::endl;
  }
  ostr << "];" << std::endl;

  ostr << varname << "_z=[";
  for (i=0; i<static_cast<index_type>(u.size()); ++i)
  {
    ostr << ps.f(u[i], v[0]).z();
    for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
    {
      ostr << ", " << ps.f(u[i], v[j]).z();
    }
    j=static_cast<index_type>(v.size()-1);
    ostr << ", " << ps.f(u[i], v[j]).z();
    if (i<static_cast<index_type>(u.size()-1))
      ostr << "; " << std::endl;
  }
  ostr << "];" << std::endl;
}

void octave_print(std::ostream &ostr, const std::string &varname, const surface_type &bez)
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
  curve_control_point_type cp[7];
  piecewise_curve_type pc1, pc2;
  curve_type c;

  dia = len/fine_ratio;
  cp[0] << 0.00*len, 0, 0.00*dia;
  cp[1] << 0.05*len, 0, 0.95*dia;
  cp[2] << 0.20*len, 0, 1.00*dia;
  cp[3] << 0.50*len, 0, 1.00*dia;
  cp[4] << 0.60*len, 0, 1.00*dia;
  cp[5] << 0.95*len, 0, 0.30*dia;
  cp[6] << 1.00*len, 0, 0.00*dia;

  // create top reference curve
  c.resize(3);
  for (index_type i=0; i<4; ++i)
  {
    c.set_control_point(cp[i], i);
  }
  pc1.push_back(c);
  for (index_type i=3; i<7; ++i)
  {
    c.set_control_point(cp[i], i-3);
  }
  pc1.push_back(c);

  // create bottom reference curve
  for (index_type i=0; i<7; ++i)
  {
    cp[i].z()=-cp[i].z();
  }
  for (index_type i=0; i<4; ++i)
  {
    c.set_control_point(cp[i], i);
  }
  pc2.push_back(c);
  for (index_type i=3; i<7; ++i)
  {
    c.set_control_point(cp[i], i-3);
  }
  pc2.push_back(c);

  piecewise_surface_type ps;

  // create the pod surface from the piecewise curve
  eli::geom::surface::create_body_of_revolution(ps, pc1, 0);

  // print out octave code
  octave_print(std::cout, "bez_c", pc1);
  octave_print(std::cout, "bez_c2", pc2);
  octave_print(std::cout, "bez_s", ps);
  std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
  std::cout << "mesh(bez_s_x, bez_s_y, bez_s_z, zeros(size(bez_s_z)), 'EdgeColor', [0 0 0]);" << std::endl;
  std::cout << "hold on;" << std::endl;
  std::cout << "plot3(bez_c_x, bez_c_y, bez_c_z, '-r');" << std::endl;
  std::cout << "plot3(bez_c2_x, bez_c2_y, bez_c2_z, '-r');" << std::endl;
  std::cout << "hold off;" << std::endl;
  std::cout << "axis equal;" << std::endl;

  return EXIT_SUCCESS;
}

