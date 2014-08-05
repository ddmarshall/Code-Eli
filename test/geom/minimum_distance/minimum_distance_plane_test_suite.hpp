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

#ifndef minimum_distance_plane_test_suite_hpp
#define minimum_distance_plane_test_suite_hpp

#include <cmath>    // cos(), sin()

#include <typeinfo> // typeid

#include "eli/util/tolerance.hpp"

#include "eli/geom/intersect/minimum_distance_plane.hpp"

template<typename data__>
class minimum_distance_plane_test_suite : public Test::Suite
{
  private:
    typedef data__ data_type;
    typedef Eigen::Matrix<data_type, 1, 3> point_type;
    typedef eli::util::tolerance<data_type> tolerance_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(minimum_distance_plane_test_suite<float>::point_3d_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_plane_test_suite<double>::point_3d_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_plane_test_suite<long double>::point_3d_test);
    }

  public:
    minimum_distance_plane_test_suite()
    {
      AddTests(data__());
    }
    ~minimum_distance_plane_test_suite()
    {
    }

  private:
    void point_3d_test()
    {
      point_type a, b, c, pt;
      data_type dist, u, v, dist_ref, u_ref, v_ref;

      // xy-plane as plane
      a << 1, 0, 0;
      b << 0, 1, 0;
      c << 0, 0, 0;
      pt << 2, 4, 2;
      dist=eli::geom::intersect::minimum_distance(u, v, a, b, c, pt);
      u_ref=2;
      v_ref=4;
      dist_ref=2;
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // xz-plane as plane
      a << 1, 0, 0;
      b << 0, 0, 1;
      c << 0, 0, 0;
      pt << 2, 4, 2;
      dist=eli::geom::intersect::minimum_distance(u, v, a, b, c, pt);
      u_ref=2;
      v_ref=2;
      dist_ref=4;
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // yz-plane as plane
      a << 0, 1, 0;
      b << 0, 0, 1;
      c << 0, 0, 0;
      pt << 2, 4, 2;
      dist=eli::geom::intersect::minimum_distance(u, v, a, b, c, pt);
      u_ref=4;
      v_ref=2;
      dist_ref=2;
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // plane degenerates as point
      a << 0, 0, 0;
      b << 0, 0, 0;
      c << 0, 0, 0;
      pt << 2, 4, 2;
      dist=eli::geom::intersect::minimum_distance(u, v, a, b, c, pt);
      u_ref=0;
      v_ref=0;
      dist_ref=pt.norm();
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // plane degenerates as line
      a << 1, 2, 1;
      b << 2, 4, 2;
      c << 1, 1, 1;
      pt << 2, 4, 2;
      dist=eli::geom::intersect::minimum_distance(u, v, a, b, c, pt);
      v_ref=0;
      dist_ref=eli::geom::intersect::minimum_distance(u_ref, c, a, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general plane at origin intersecting point
      dist_ref=0;
      u_ref=0.5;
      v_ref=1;
      a << 1, 2,-1;
      b <<-2, 1, 2;
      c << 0, 0, 0;
      pt=u_ref*a+v_ref*b+c;
      dist=eli::geom::intersect::minimum_distance(u, v, a, b, c, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general plane at origin
      dist_ref=1.5;
      u_ref=0.5;
      v_ref=1;
      a << 1, 2,-1;
      b <<-2, 1, 2;
      c << 0, 0, 0;
      pt=(u_ref*a+v_ref*b+c)+dist_ref*a.cross(b)/a.cross(b).norm();
      dist=eli::geom::intersect::minimum_distance(u, v, a, b, c, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general plane intersecting point
      dist_ref=0;
      u_ref=0.5;
      v_ref=1;
      a << 1, 2,-1;
      b <<-2, 1, 2;
      c << 1,-1, 3;
      pt=(u_ref*a+v_ref*b+c)+dist_ref*a.cross(b)/a.cross(b).norm();
      dist=eli::geom::intersect::minimum_distance(u, v, a, b, c, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // general plane
      dist_ref=1;
      u_ref=2.5;
      v_ref=0.5;
      a << 1, 2,-1;
      b <<-2, 1, 2;
      c << 1,-1, 3;
      pt=(u_ref*a+v_ref*b+c)+dist_ref*a.cross(b)/a.cross(b).norm();
      dist=eli::geom::intersect::minimum_distance(u, v, a, b, c, pt);
      TEST_ASSERT(tol.approximately_equal(u, u_ref));
      TEST_ASSERT(tol.approximately_equal(v, v_ref));
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
    }
};

#endif

