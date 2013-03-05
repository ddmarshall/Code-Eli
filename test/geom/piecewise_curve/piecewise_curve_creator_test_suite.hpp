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

#ifndef piecewise_curve_creator_test_suite_hpp
#define piecewise_curve_creator_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_creator.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_curve_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_creator_test_suite<float>::create_cubic_spline_fd_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_creator_test_suite<double>::create_cubic_spline_fd_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_creator_test_suite<long double>::create_cubic_spline_fd_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_creator_test_suite<dd_real>::create_cubic_spline_fd_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_creator_test_suite<qd_real>::create_cubic_spline_fd_test);
    }
#endif
  public:
    piecewise_curve_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_curve_creator_test_suite()
    {
    }

  private:
    void create_cubic_spline_fd_test()
    {
      // create with specified times
      {
        piecewise_curve_type c1;
        std::vector<point_type, Eigen::aligned_allocator<point_type> > pts(5);
        std::vector<data_type> t(5);
        index_type i, nseg;

        // set the points and times
        pts[0] <<  1, 3,  2;
        pts[1] <<  0, 2,  3;
        pts[2] << -1, 4,  1;
        pts[3] <<  0, 5,  0;
        pts[4] <<  1, 6, -1;
        t[0]=1;
        t[1]=3;
        t[2]=4;
        t[3]=7;
        t[4]=9;

        // create the piecewise curve
        TEST_ASSERT(create_piecewise_cubic(c1, pts.begin(), pts.end(), t.begin()));

        // check the number segments
        nseg=c1.number_segments();
        TEST_ASSERT(nseg==static_cast<index_type>(pts.size()-1));

        // check the starting time
        TEST_ASSERT(tol.approximately_equal(t[0], c1.get_t0()));

        // check the continuity at each point
        curve_type seg1, seg2;
        data_type t_out;
        c1.get(seg1, t_out, 0);
        TEST_ASSERT(tol.approximately_equal(t_out, t[1]-t[0]));
        for (i=1; i<nseg; ++i)
        {
          c1.get(seg2, t_out, i);
          TEST_ASSERT(tol.approximately_equal(t_out, t[i+1]-t[i]));
          TEST_ASSERT(tol.approximately_equal(seg1.fp(1), seg2.fp(0)));
          TEST_ASSERT(!tol.approximately_equal(seg1.fpp(1), seg2.fpp(0)));
          seg1=seg2;
        }
      }

      // create with default times
      {
      }
    }
};

#endif

