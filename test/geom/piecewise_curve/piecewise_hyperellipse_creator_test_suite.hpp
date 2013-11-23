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

#ifndef piecewise_hyperellipse_creator_test_suite_hpp
#define piecewise_hyperellipse_creator_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_hyperellipse_creator.hpp"

#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_hyperellipse_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::piecewise_hyperellipse_creator<data__, 3, tolerance_type> hyperellipse_creator_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_hyperellipse_creator_test_suite<float>::create_4seg3deg_test);
      TEST_ADD(piecewise_hyperellipse_creator_test_suite<float>::create_4seg6deg_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_hyperellipse_creator_test_suite<double>::create_4seg3deg_test);
      TEST_ADD(piecewise_hyperellipse_creator_test_suite<double>::create_4seg6deg_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_hyperellipse_creator_test_suite<long double>::create_4seg3deg_test);
      TEST_ADD(piecewise_hyperellipse_creator_test_suite<long double>::create_4seg6deg_test);
    }
#ifdef ELI_USING_QD
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_hyperellipse_creator_test_suite<dd_real>::create_4seg3deg_test);
      TEST_ADD(piecewise_hyperellipse_creator_test_suite<dd_real>::create_4seg6deg_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_hyperellipse_creator_test_suite<qd_real>::create_4seg3deg_test);
      TEST_ADD(piecewise_hyperellipse_creator_test_suite<qd_real>::create_4seg6deg_test);
    }
#endif
  public:
    piecewise_hyperellipse_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_hyperellipse_creator_test_suite()
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

    void create_4seg3deg_test()
    {
      hyperellipse_creator_type he_creator(4);
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

        fref << 1.54543, 1.91076, 0;
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

        fref << 0.390768, 1.00505, 0;
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

        fref << 0.020186, 0.66746, 0;
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

        fref << 1.1795355, 1.6218755, 0;
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

        fref << 1.905055, 2.234862, 0;
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

        fref << 1.219502, 1.657993, 0;
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

        fref << 1.21118, 1.647418, 0;
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
      hyperellipse_creator_type he_creator(4);
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

        fref << 1.5440357, 1.906499, 0;
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

        fref << 0.350574, 1.00775, 0;
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

        fref << -0.04397, 0.52341, 0;
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

        fref << 1.168494, 1.614415, 0;
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

        fref << 1.89498, 2.2451515, 0;
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

        fref << 1.083538, 1.6076799, 0;
        f=pc.f(t0+dt0/2);
        if(typeid(data_type)==typeid(float))
        {
          TEST_ASSERT((f-fref).norm() < 5e-4);
        }
        else
        {
          TEST_ASSERT((f-fref).norm() < 5e-6);
        }
//         std::cout << "f=" << std::setprecision(12) << f << std::endl;
//         std::cout << "diff=" << std::setprecision(12) << (f-fref).norm() << std::endl;
      }

      // create a mixed concave-convex
      {
        piecewise_curve_type pc;
        he_creator.set_exponents(3., 1./3.);

        TEST_ASSERT(he_creator.create(pc));

        fref << 0.928627, 1.545853, 0;
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

