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

#ifndef piecewise_surface_creator_test_suite_suite_hpp
#define piecewise_surface_creator_test_suite_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"

#include "eli/geom/curve/piecewise.hpp"

#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/piecewise_creator.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_surface_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;
    typedef typename piecewise_surface_type::tolerance_type tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_creator_test_suite<float>::create_body_of_revolution_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_creator_test_suite<double>::create_body_of_revolution_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_surface_creator_test_suite<long double>::create_body_of_revolution_test);
    }

  public:
    piecewise_surface_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_surface_creator_test_suite()
    {
    }

  private:
    void octave_print(int figno, const piecewise_surface_type &ps) const
    {
      index_type i, j, pp, qq, nup, nvp;
      data_type umin, vmin, umax, vmax;

      nup=ps.number_u_patches();
      nvp=ps.number_v_patches();
      ps.get_parameter_min(umin, vmin);
      ps.get_parameter_max(umax, vmax);

      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << "cp_x=[" << std::endl;
      for (pp=0; pp<nup; ++pp)
      {
        for (qq=0; qq<nvp; ++qq)
        {
          surface_type bez;
          ps.get(bez, pp, qq);
          for (i=0; i<=bez.degree_u(); ++i)
          {
            std::cout << bez.get_control_point(i, 0).x();
            for (j=1; j<bez.degree_v(); ++j)
            {
              std::cout << ", " << bez.get_control_point(i, j).x();
            }
            j=bez.degree_v();
            std::cout << ", " << bez.get_control_point(i, j).x();
            if (i<bez.degree_u())
              std::cout << "; ";
            else if (pp<nup-1)
              std::cout << "; ";
            else if (qq<nvp-1)
              std::cout << "; ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << "];" << std::endl;

      std::cout << "cp_y=[";
      for (pp=0; pp<nup; ++pp)
      {
        for (qq=0; qq<nvp; ++qq)
        {
          surface_type bez;
          ps.get(bez, pp, qq);
          for (i=0; i<=bez.degree_u(); ++i)
          {
            std::cout << bez.get_control_point(i, 0).y();
            for (j=1; j<bez.degree_v(); ++j)
            {
              std::cout << ", " << bez.get_control_point(i, j).y();
            }
            j=bez.degree_v();
            std::cout << ", " << bez.get_control_point(i, j).y();
            if (i<bez.degree_u())
              std::cout << "; ";
            else if (pp<nup-1)
              std::cout << "; ";
            else if (qq<nvp-1)
              std::cout << "; ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << "];" << std::endl;

      std::cout << "cp_z=[";
      for (pp=0; pp<nup; ++pp)
      {
        for (qq=0; qq<nvp; ++qq)
        {
          surface_type bez;
          ps.get(bez, pp, qq);
          for (i=0; i<=bez.degree_u(); ++i)
          {
            std::cout << bez.get_control_point(i, 0).z();
            for (j=1; j<bez.degree_v(); ++j)
            {
              std::cout << ", " << bez.get_control_point(i, j).z();
            }
            j=bez.degree_v();
            std::cout << ", " << bez.get_control_point(i, j).z();
            if (i<bez.degree_u())
              std::cout << "; ";
            else if (pp<nup-1)
              std::cout << "; ";
            else if (qq<nvp-1)
              std::cout << "; ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << "];" << std::endl;

      // initialize the u & v parameters
      std::vector<data__> u(11), v(11);
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        u[i]=umin+(umax-umin)*static_cast<data__>(i)/(u.size()-1);
      }
      for (j=0; j<static_cast<index_type>(v.size()); ++j)
      {
        v[j]=vmin+(vmax-vmin)*static_cast<data__>(j)/(v.size()-1);
      }

      // set the surface points
      std::cout << "surf_x=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << ps.f(u[i], v[0]).x();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << ps.f(u[i], v[j]).x();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << ps.f(u[i], v[j]).x();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_y=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << ps.f(u[i], v[0]).y();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << ps.f(u[i], v[j]).y();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << ps.f(u[i], v[j]).y();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "surf_z=[";
      for (i=0; i<static_cast<index_type>(u.size()); ++i)
      {
        std::cout << ps.f(u[i], v[0]).z();
        for (j=1; j<static_cast<index_type>(v.size()-1); ++j)
        {
          std::cout << ", " << ps.f(u[i], v[j]).z();
        }
        j=static_cast<index_type>(v.size()-1);
        std::cout << ", " << ps.f(u[i], v[j]).z();
        if (i<static_cast<index_type>(u.size()-1))
          std::cout << "; " << std::endl;
      }
      std::cout << "];" << std::endl;

      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "mesh(surf_x, surf_y, surf_z, zeros(size(surf_z)), 'EdgeColor', [0 0 0]);" << std::endl;
      std::cout << "hold on;" << std::endl;
      std::cout << "plot3(cp_x, cp_y, cp_z, 'ok', 'MarkerFaceColor', [0 0 0]);" << std::endl;
      std::cout << "hold off;" << std::endl;
    }

    void create_body_of_revolution_test()
    {
      typedef typename piecewise_curve_type::curve_type curve_type;
      typedef typename curve_type::control_point_type curve_control_point_type;

      // create geometry with default parameterizations
      {
        piecewise_curve_type pc;
        curve_type c(3);
        curve_control_point_type cp[4];
        data_type k=4*(eli::constants::math<data_type>::sqrt_two()-1)/3;
        index_type i;

        // create curve
        cp[0] << 1, 0, 0;
        cp[1] << 1, k, 0;
        cp[2] << k, 1, 0;
        cp[3] << 0, 1, 0;
        for (i=0; i<4; ++i)
        {
          c.set_control_point(cp[i], i);
        }
        TEST_ASSERT(pc.push_back(c, 0.25)==piecewise_curve_type::NO_ERROR);

        // set 2nd quadrant curve
        cp[0] <<  0, 1, 0;
        cp[1] << -k, 1, 0;
        cp[2] << -1, k, 0;
        cp[3] << -1, 0, 0;
        for (i=0; i<4; ++i)
        {
          c.set_control_point(cp[i], i);
        }
        TEST_ASSERT(pc.push_back(c, 0.25)==piecewise_curve_type::NO_ERROR);

        piecewise_surface_type ps;

        TEST_ASSERT(eli::geom::surface::create_body_of_revolution(ps, pc, 0, eli::geom::surface::OUTWARD_NORMAL));

//         if (typeid(data_type)==typeid(double))
//           octave_print(2, ps);

        TEST_ASSERT(ps.open_u());
        TEST_ASSERT(ps.closed_v());
      }
    }
};

#endif

