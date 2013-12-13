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

#ifndef piecewise_superellipse_creator_test_suite_hpp
#define piecewise_superellipse_creator_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_superellipse_creator.hpp"

#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_superellipse_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::piecewise_superellipse_creator<data__, 3, tolerance_type> superellipse_creator_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_superellipse_creator_test_suite<float>::create_degenerate_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<float>::create_4seg3deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<float>::create_4seg6deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<float>::create_6seg3deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<float>::create_6seg6deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<float>::create_8seg3deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<float>::create_8seg6deg_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_superellipse_creator_test_suite<double>::create_degenerate_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<double>::create_4seg3deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<double>::create_4seg6deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<double>::create_6seg3deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<double>::create_6seg6deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<double>::create_8seg3deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<double>::create_8seg6deg_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_superellipse_creator_test_suite<long double>::create_degenerate_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<long double>::create_4seg3deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<long double>::create_4seg6deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<long double>::create_6seg3deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<long double>::create_6seg6deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<long double>::create_8seg3deg_test);
      TEST_ADD(piecewise_superellipse_creator_test_suite<long double>::create_8seg6deg_test);
    }

  public:
    piecewise_superellipse_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_superellipse_creator_test_suite()
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

      // initialize the t parameters
      std::vector<data__> t(129);
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        t[i]=tmin+(tmax-tmin)*static_cast<data__>(i)/(t.size()-1);
      }

      // set the curve points
      std::cout << "curve_x=[";
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        std::cout << pc.f(t[i]).x();
        if (i<static_cast<index_type>(t.size()-1))
          std::cout << ", ";
      }
      std::cout << "];" << std::endl;

      std::cout << "curve_y=[";
      for (i=0; i<static_cast<index_type>(t.size()); ++i)
      {
        std::cout << pc.f(t[i]).y();
        if (i<static_cast<index_type>(t.size()-1))
          std::cout << ", ";
      }
      std::cout << "];" << std::endl;

      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "plot(curve_x, curve_y, '-k');" << std::endl;
      std::cout << "hold on;" << std::endl;
      std::cout << "plot(cp_x', cp_y', '-or', 'MarkerFaceColor', [0 0 0]);" << std::endl;
      std::cout << "hold off;" << std::endl;
    }

    void create_degenerate_test()
    {
      superellipse_creator_type he_creator(4);
      data_type dt0(3), dt1(2), dt2(3), dt3(2), t0(-1);
      point_type f, fref;

      // set the times
      he_creator.set_t0(t0);
      he_creator.set_segment_dt(dt0, 0);
      he_creator.set_segment_dt(dt1, 1);
      he_creator.set_segment_dt(dt2, 2);
      he_creator.set_segment_dt(dt3, 3);

      he_creator.set_exponents(2., 2.);
      he_creator.set_max_degree(3);

      // create an x-degenerate ellipse
      {
        piecewise_curve_type pc;

        he_creator.set_axis(0, 3);
        TEST_ASSERT(he_creator.create(pc));

        fref << 0, 1.5, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
      }

      // create an y-degenerate ellipse
      {
        piecewise_curve_type pc;

        he_creator.set_axis(2, 0);
        TEST_ASSERT(he_creator.create(pc));

        fref << 1, 0, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
      }

      // create an x- and y-degenerate ellipse
      {
        piecewise_curve_type pc;

        he_creator.set_axis(0, 0);
        TEST_ASSERT(he_creator.create(pc));

        fref << 0, 0, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
      }
    }

    void create_4seg3deg_test()
    {
      superellipse_creator_type he_creator(4);
      data_type dt0(3), dt1(2), dt2(3), dt3(2), t0(-1);
      point_type f, fref;

      // set the times
      he_creator.set_t0(t0);
      he_creator.set_segment_dt(dt0, 0);
      he_creator.set_segment_dt(dt1, 1);
      he_creator.set_segment_dt(dt2, 2);
      he_creator.set_segment_dt(dt3, 3);

      he_creator.set_axis(2, 3);
      he_creator.set_max_degree(3);

      // create an ellipse
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(2., 2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.547913, 1.903984, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a moderately concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./2., 1./2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 0.397104, 1.058210, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a severely concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./5., 1./5.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 0.134030, 1.149813, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a moderately convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1.2, 1.2);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.173693, 1.617232, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a severely convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(5., 5.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.76445, 1.897658, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./3., 3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.1199889, 1.604568, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(3., 1./3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.149061, 1.645629, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
    }

    void create_4seg6deg_test()
    {
      superellipse_creator_type he_creator(4);
      data_type dt0(3), dt1(2), dt2(3), dt3(2), t0(-1);
      point_type f, fref;

      // set the times
      he_creator.set_t0(t0);
      he_creator.set_segment_dt(dt0, 0);
      he_creator.set_segment_dt(dt1, 1);
      he_creator.set_segment_dt(dt2, 2);
      he_creator.set_segment_dt(dt3, 3);

      he_creator.set_axis(2, 3);
      he_creator.set_max_degree(6);

      // create an ellipse
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(2., 2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.544666, 1.904765, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a moderately concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./2., 1./2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 0.360349, 1.00602, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a severely concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./5., 1./5.);

        TEST_ASSERT(he_creator.create(pc));

        if (typeid(data_type)==typeid(float))
        {
          fref << 47.022301, 3.482945, 0;
        }
        else
        {
          fref << 46.87239, 3.468185, 0;
        }
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a moderately convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1.2, 1.2);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.168557, 1.614039, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a severely convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(5., 5.);

        TEST_ASSERT(he_creator.create(pc));

        if (typeid(data_type)==typeid(float))
        {
          fref << 1.351337, 2.2768044, 0;
        }
        else
        {
          fref << 1.349747, 2.277204, 0;
        }
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 6e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./3., 3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.203145, 1.619022, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(3., 1./3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.135717, 1.623672, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;

//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << std::setprecision(6);
//           octave_print(1, pc);
//         }
      }
    }

    void create_6seg3deg_test()
    {
      superellipse_creator_type he_creator(6);
      data_type dt0(3), dt1(2), dt2(3), dt3(2), dt4(3), dt5(2), t0(-1);
      point_type f, fref;

      // set the times
      he_creator.set_t0(t0);
      he_creator.set_segment_dt(dt0, 0);
      he_creator.set_segment_dt(dt1, 1);
      he_creator.set_segment_dt(dt2, 2);
      he_creator.set_segment_dt(dt3, 3);
      he_creator.set_segment_dt(dt4, 4);
      he_creator.set_segment_dt(dt5, 5);

      he_creator.set_axis(2, 3);
      he_creator.set_max_degree(3);

      // create an ellipse
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(2., 2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.819996, 1.243964, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
      }

      // create a moderately concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./2., 1./2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.019236, 0.267005, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
      }

      // create a severely concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./5., 1./5.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.108429, 0.525021, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
      }

      // create a moderately convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1.2, 1.2);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.530603, 1.024579, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
      }

      // create a severely convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(5., 5.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.926631, 1.313695, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./3., 3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.664021, 1.092592, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(3., 1./3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.333809, 0.985460, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
      }
    }

    void create_6seg6deg_test()
    {
      superellipse_creator_type he_creator(6);
      data_type dt0(3), dt1(2), dt2(3), dt3(2), dt4(3), dt5(2), t0(-1);
      point_type f, fref;

      // set the times
      he_creator.set_t0(t0);
      he_creator.set_segment_dt(dt0, 0);
      he_creator.set_segment_dt(dt1, 1);
      he_creator.set_segment_dt(dt2, 2);
      he_creator.set_segment_dt(dt3, 3);
      he_creator.set_segment_dt(dt4, 4);
      he_creator.set_segment_dt(dt5, 5);

      he_creator.set_axis(2, 3);
      he_creator.set_max_degree(6);

      // create an ellipse
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(2., 2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.791588, 1.333412, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(float))
//         {
//           std::cout << std::setprecision(6);
//           octave_print(1, pc);
//         }
      }

      // create a moderately concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./2., 1./2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 0.931743, 0.302345, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << std::setprecision(6);
//           octave_print(1, pc);
//         }
      }

      // create a severely concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./5., 1./5.);

        TEST_ASSERT(he_creator.create(pc));

        if (typeid(data_type)==typeid(float))
        {
          fref << -19.968788, 9.096475, 0;
        }
        else
        {
          fref << -14.07074, 6.209714, 0;
        }
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << std::setprecision(6);
//           octave_print(1, pc);
//         }
      }

      // create a moderately convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1.2, 1.2);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.471385, 1.124692, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << std::setprecision(6);
//           octave_print(1, pc);
//         }
      }

      // create a severely convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(5., 5.);

        TEST_ASSERT(he_creator.create(pc));

        if (typeid(data_type)==typeid(float))
        {
          fref << 1.928891, 1.436470, 0;
        }
        else
        {
          fref << 1.935626, 1.432817, 0;
        }
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 6e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << std::setprecision(6);
//           octave_print(1, pc);
//         }
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./3., 3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.647633, 1.190923, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << std::setprecision(6);
//           octave_print(1, pc);
//         }
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(3., 1./3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.314240, 1.100464, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << std::setprecision(6);
//           octave_print(1, pc);
//         }
      }
    }

    void create_8seg3deg_test()
    {
      superellipse_creator_type he_creator(8);
      data_type dt0(3), dt1(2), dt2(3), dt3(2), dt4(2), dt5(1), dt6(2), dt7(1), t0(-1);
      point_type f, fref;

      // set the times
      he_creator.set_t0(t0);
      he_creator.set_segment_dt(dt0, 0);
      he_creator.set_segment_dt(dt1, 1);
      he_creator.set_segment_dt(dt2, 2);
      he_creator.set_segment_dt(dt3, 3);
      he_creator.set_segment_dt(dt4, 4);
      he_creator.set_segment_dt(dt5, 5);
      he_creator.set_segment_dt(dt6, 6);
      he_creator.set_segment_dt(dt7, 7);

      he_creator.set_axis(2, 3);
      he_creator.set_max_degree(3);

      // create an ellipse
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(2., 2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.859473, 1.104687, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a moderately concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./2., 1./2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.147263, 0.188791, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a severely concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./5., 1./5.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.209256, 0.430468, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a moderately convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1.2, 1.2);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.616116, 0.869995, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a severely convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(5., 5.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.937615, 1.281272, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./3., 3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.774508, 0.938826, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(3., 1./3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.407137, 0.785128, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

//         if (typeid(data_type)==typeid(double))
//           octave_print(1, pc);
    }

    void create_8seg6deg_test()
    {
      superellipse_creator_type he_creator(8);
      data_type dt0(3), dt1(2), dt2(3), dt3(2), dt4(2), dt5(1), dt6(2), dt7(1), t0(-1);
      point_type f, fref;

      // set the times
      he_creator.set_t0(t0);
      he_creator.set_segment_dt(dt0, 0);
      he_creator.set_segment_dt(dt1, 1);
      he_creator.set_segment_dt(dt2, 2);
      he_creator.set_segment_dt(dt3, 3);
      he_creator.set_segment_dt(dt4, 4);
      he_creator.set_segment_dt(dt5, 5);
      he_creator.set_segment_dt(dt6, 6);
      he_creator.set_segment_dt(dt7, 7);

      he_creator.set_axis(2, 3);
      he_creator.set_max_degree(6);

      // create an ellipse
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(2., 2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.859616, 1.104135, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a moderately concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./2., 1./2.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.148073, 0.176278, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a severely concave
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./5., 1./5.);

        TEST_ASSERT(he_creator.create(pc));

        if (typeid(data_type)==typeid(float))
        {
          fref << -1.324945, 3.099777, 0;
        }
        else
        {
          fref << -1.799115, 3.268398, 0;
        }
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a moderately convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1.2, 1.2);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.615210, 0.869275, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a severely convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(5., 5.);

        TEST_ASSERT(he_creator.create(pc));

        if (typeid(data_type)==typeid(float))
        {
          fref << 1.975022, 1.332086, 0;
        }
        else
        {
          fref << 1.971965, 1.328557, 0;
        }
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 6e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(1./3., 3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.822104, 0.938021, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(3., 1./3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 1.423473, 0.786496, 0;
        f=pc.f(t0+dt0/2);
        TEST_ASSERT((f-fref).norm() < 5e-6);
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;

//         if (typeid(data_type)==typeid(double))
//         {
//           std::cout << std::setprecision(6);
//           octave_print(1, pc);
//         }
      }
    }
};

#endif

