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

#ifndef minimum_distance_point_test_suite_hpp
#define minimum_distance_point_test_suite_hpp

#include <cassert>  // assert()
#include <cmath>    // cos(), sin()

#include <typeinfo> // typeid

#include "Eigen/Eigen"

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/geom/intersect/minimum_distance_point.hpp"

template<typename data__>
class minimum_distance_point_test_suite : public Test::Suite
{
  private:
    typedef data__ data_type;
    typedef Eigen::Matrix<data_type, 1, 1> point_type1;
    typedef Eigen::Matrix<data_type, 1, 2> point_type2;
    typedef Eigen::Matrix<data_type, 1, 3> point_type3;

    typedef typename point_type1::Index index_type;

    eli::util::tolerance<data_type> tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(minimum_distance_point_test_suite<float>::simple_1d_test);
      TEST_ADD(minimum_distance_point_test_suite<float>::simple_2d_test);
      TEST_ADD(minimum_distance_point_test_suite<float>::simple_3d_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_point_test_suite<double>::simple_1d_test);
      TEST_ADD(minimum_distance_point_test_suite<double>::simple_2d_test);
      TEST_ADD(minimum_distance_point_test_suite<double>::simple_3d_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_point_test_suite<long double>::simple_1d_test);
      TEST_ADD(minimum_distance_point_test_suite<long double>::simple_2d_test);
      TEST_ADD(minimum_distance_point_test_suite<long double>::simple_3d_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(minimum_distance_point_test_suite<dd_real>::simple_1d_test);
      TEST_ADD(minimum_distance_point_test_suite<dd_real>::simple_2d_test);
      TEST_ADD(minimum_distance_point_test_suite<dd_real>::simple_3d_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(minimum_distance_point_test_suite<qd_real>::simple_1d_test);
      TEST_ADD(minimum_distance_point_test_suite<qd_real>::simple_2d_test);
      TEST_ADD(minimum_distance_point_test_suite<qd_real>::simple_3d_test);
    }
#endif
  public:
    minimum_distance_point_test_suite()
    {
      AddTests(data__());
    }
    ~minimum_distance_point_test_suite()
    {
    }

  private:
    void simple_1d_test()
    {
      point_type1 pt1, pt2;
      data_type dist;

      pt1 << 1;
      pt2 << 2;

      dist=(pt1-pt2).norm();
      TEST_ASSERT(tol.approximately_equal(dist, eli::geom::intersect::minimum_distance(pt1, pt2)));
    }

    void simple_2d_test()
    {
      point_type2 pt1, pt2;
      data_type dist;

      pt1 << 1, 1;
      pt2 << 2, 2;

      dist=(pt1-pt2).norm();
      TEST_ASSERT(tol.approximately_equal(dist, eli::geom::intersect::minimum_distance(pt1, pt2)));
    }

    void simple_3d_test()
    {
      point_type3 pt1, pt2;
      data_type dist;

      pt1 << 1, 1, 2;
      pt2 << 2, 2, 3;

      dist=(pt1-pt2).norm();
      TEST_ASSERT(tol.approximately_equal(dist, eli::geom::intersect::minimum_distance(pt1, pt2)));
    }
};

#endif

