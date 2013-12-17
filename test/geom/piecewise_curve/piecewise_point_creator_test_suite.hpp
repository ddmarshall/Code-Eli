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

#ifndef piecewise_point_creator_test_suite_hpp
#define piecewise_point_creator_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_point_creator.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_point_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::piecewise_point_creator<data__, 3, tolerance_type> point_creator_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_point_creator_test_suite<float>::create_point_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_point_creator_test_suite<double>::create_point_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_point_creator_test_suite<long double>::create_point_test);
    }

  public:
    piecewise_point_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_point_creator_test_suite()
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

    void create_point_test()
    {
      // create point with specified parameterization
      {
        piecewise_curve_type pc;
        point_creator_type pt_creator(4);
        data_type dt0(3), dt1(2), dt2(3), dt3(2), t0(-1), dt;
        point_type p0, ptemp;

        // set the point
        p0 << 1, 2, 3;
        pt_creator.set_point(p0);

        // set the times
        pt_creator.set_t0(t0);
        pt_creator.set_segment_dt(dt0, 0);
        pt_creator.set_segment_dt(dt1, 1);
        pt_creator.set_segment_dt(dt2, 2);
        pt_creator.set_segment_dt(dt3, 3);

        // test point setting
        ptemp=pt_creator.get_point();
        TEST_ASSERT(p0==ptemp);

        // test time step settings
        TEST_ASSERT(pt_creator.get_t0()==t0);
        dt=pt_creator.get_segment_dt(0);
        TEST_ASSERT(dt==dt0);
        dt=pt_creator.get_segment_dt(1);
        TEST_ASSERT(dt==dt1);
        dt=pt_creator.get_segment_dt(2);
        TEST_ASSERT(dt==dt2);
        dt=pt_creator.get_segment_dt(3);
        TEST_ASSERT(dt==dt3);

        // create the polygon
        TEST_ASSERT(pt_creator.create(pc));
      }

      // create point with default parameterization
      {
        piecewise_curve_type pc;
        point_creator_type pt_creator(4);
        data_type dt;
        point_type p0, ptemp;

        // set the corners
        p0 << 1, 0, 0;
        pt_creator.set_point(p0);

        // test corner settings
        ptemp=pt_creator.get_point();
        TEST_ASSERT(p0==ptemp);

        // test time step settings
        TEST_ASSERT(pt_creator.get_t0()==0);
        dt=pt_creator.get_segment_dt(0);
        TEST_ASSERT(dt==1);
        dt=pt_creator.get_segment_dt(1);
        TEST_ASSERT(dt==1);
        dt=pt_creator.get_segment_dt(2);
        TEST_ASSERT(dt==1);
        dt=pt_creator.get_segment_dt(3);
        TEST_ASSERT(dt==1);

        // create the polygon
        TEST_ASSERT(pt_creator.create(pc));
      }
    }
};

#endif

