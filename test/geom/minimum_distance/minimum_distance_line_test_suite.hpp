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

#ifndef minimum_distance_line_test_suite_hpp
#define minimum_distance_line_test_suite_hpp

#include <cassert>  // assert()
#include <cmath>    // cos(), sin()

#include <typeinfo> // typeid

#include "Eigen/Eigen"

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/geom/intersect/minimum_distance_line.hpp"

template<typename data__>
class minimum_distance_line_test_suite : public Test::Suite
{
  private:
    typedef data__ data_type;
    typedef Eigen::Matrix<data_type, 1, 2> point_type2;
    typedef Eigen::Matrix<data_type, 1, 3> point_type3;
    typedef eli::util::tolerance<data_type> tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(minimum_distance_line_test_suite<float>::point_2d_test);
      TEST_ADD(minimum_distance_line_test_suite<float>::point_3d_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_line_test_suite<double>::point_2d_test);
      TEST_ADD(minimum_distance_line_test_suite<double>::point_3d_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_line_test_suite<long double>::point_2d_test);
      TEST_ADD(minimum_distance_line_test_suite<long double>::point_3d_test);
    }

  public:
    minimum_distance_line_test_suite()
    {
      AddTests(data__());
    }
    ~minimum_distance_line_test_suite()
    {
    }

  private:
    void point_2d_test()
    {
      point_type2 pt, a0, a1;
      data_type dist, t, dist_ref, t_ref;

      // x-axis as line
      a0 << 0, 0;
      a1 << 1, 0;
      pt << 1, 2;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=1;
      dist_ref=2;
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // y-axis as line
      a0 << 0, 0;
      a1 << 0, 2;
      pt << 3, 4;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=2;
      dist_ref=3;
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // line as a point
      a0 << 0, 0;
      a1 << 0, 0;
      pt << 3, 4;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=0;
      dist_ref=pt.norm();
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general line from origin intersecting point
      a0 << 0, 0;
      a1 << 1, 3;
      pt << 2, 6;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=2;
      dist_ref=0;
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general line from origin
      a0 << 0, 0;
      a1 << 3, 4;
      pt << -1, 32; pt/=10;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=0.5;
      dist_ref=2;
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general line intersecting point
      a0 << 2,-1;
      a1 << 4, 2;
      pt << 3, -0.5;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=0.25;
      dist_ref=0;
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general line
      a0 << 2, 1;
      a1 << 3,-4;
      pt << 57, -56; pt/=10;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=1.5;
      dist_ref=1;
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
    }

    void point_3d_test()
    {
      point_type3 a0, a1, pt;
      data_type dist, t, dist_ref, t_ref;

      // x-axis as line
      a0 << 0, 0, 0;
      a1 << 1, 0, 0;
      pt << 1, 2, 2;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=1;
      dist_ref=std::sqrt(8);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // y-axis as line
      a0 << 0, 0, 0;
      a1 << 0, 2, 0;
      pt << 3, 4, -1;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=2;
      dist_ref=std::sqrt(9+1);
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // line as a point
      a0 << 0, 0, 0;
      a1 << 0, 0, 0;
      pt << 3, 4, 1;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=0;
      dist_ref=pt.norm();
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general line from origin intersecting point
      a0 << 0, 0, 0;
      a1 << 1, 3, 2;
      pt << 2, 6, 4;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=2;
      dist_ref=0;
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general line from origin
      a0 << 0, 0, 0;
      a1 << 1, 2, 2;
      pt << 3, 3, 0;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=1;
      dist_ref=3;
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general line intersecting point
      a0 << 2,-1, 2;
      a1 << 4, 2, 6;
      pt << 3, -0.5, 3.5;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=0.25;
      dist_ref=0;
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general line
      a0 << 2, 1, 1;
      a1 << 1, 2, 2;
      pt <<17,13,10; pt/=6;
      dist=eli::geom::intersect::minimum_distance(t, a0, a1, pt);
      t_ref=0.5;
      dist_ref=0.5;
      TEST_ASSERT(tol.approximately_equal(t, t_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
    }
};

#endif

