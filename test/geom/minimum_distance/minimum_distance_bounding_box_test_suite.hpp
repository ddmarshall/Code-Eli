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

#ifndef minimum_distance_bounding_box_test_suite_hpp
#define minimum_distance_bounding_box_test_suite_hpp

#include <cassert>  // assert()
#include <cmath>    // cos(), sin()

#include <typeinfo> // typeid

#include "Eigen/Eigen"

#include "eli/code_eli.hpp"

#include "eli/geom/intersect/minimum_distance_bounding_box.hpp"
#include "eli/geom/general/bounding_box.hpp"

template<typename data__>
class minimum_distance_bounding_box_test_suite : public Test::Suite
{
  private:
    typedef data__ data_type;
    typedef eli::geom::general::bounding_box<data__, 1> bounding_box_type1;
    typedef typename bounding_box_type1::point_type point_type1;
    typedef eli::geom::general::bounding_box<data__, 2> bounding_box_type2;
    typedef typename bounding_box_type2::point_type point_type2;
    typedef eli::geom::general::bounding_box<data__, 3> bounding_box_type3;
    typedef typename bounding_box_type3::point_type point_type3;

    typedef typename bounding_box_type1::tolerance_type tolerance_type;
    typedef typename bounding_box_type1::index_type index_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(minimum_distance_bounding_box_test_suite<float>::point_1d_test);
      TEST_ADD(minimum_distance_bounding_box_test_suite<float>::point_2d_test);
      TEST_ADD(minimum_distance_bounding_box_test_suite<float>::point_3d_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_bounding_box_test_suite<double>::point_1d_test);
      TEST_ADD(minimum_distance_bounding_box_test_suite<double>::point_2d_test);
      TEST_ADD(minimum_distance_bounding_box_test_suite<double>::point_3d_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(minimum_distance_bounding_box_test_suite<long double>::point_1d_test);
      TEST_ADD(minimum_distance_bounding_box_test_suite<long double>::point_2d_test);
      TEST_ADD(minimum_distance_bounding_box_test_suite<long double>::point_3d_test);
    }
#ifdef ELI_USING_QD
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(minimum_distance_bounding_box_test_suite<dd_real>::point_1d_test);
      TEST_ADD(minimum_distance_bounding_box_test_suite<dd_real>::point_2d_test);
      TEST_ADD(minimum_distance_bounding_box_test_suite<dd_real>::point_3d_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(minimum_distance_bounding_box_test_suite<qd_real>::point_1d_test);
      TEST_ADD(minimum_distance_bounding_box_test_suite<qd_real>::point_2d_test);
      TEST_ADD(minimum_distance_bounding_box_test_suite<qd_real>::point_3d_test);
    }
#endif
  public:
    minimum_distance_bounding_box_test_suite()
    {
      AddTests(data__());
    }
    ~minimum_distance_bounding_box_test_suite()
    {
    }

  private:
    void point_1d_test()
    {
      bounding_box_type1 bb;
      point_type1 pt;
      data_type dist, dist_ref;

      pt << -1;
      bb.add(pt);
      pt << 1;
      bb.add(pt);

      // test lower than min
      pt << -3;
      dist_ref=std::abs(pt.x()-bb.get_min().x());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test higher than max
      pt << 3;
      dist_ref=std::abs(pt.x()-bb.get_max().x());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside
      pt << 0;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
      pt << 0.5;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test on point
      pt << 1;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
    }

    void point_2d_test()
    {
      bounding_box_type2 bb;
      point_type2 pt, corner;
      data_type dist, dist_ref;

      pt << -1, -1;
      bb.add(pt);
      pt << 1, 1;
      bb.add(pt);

      // test lower than x-min, lower than y-min
      pt << -3, -3;
      corner << bb.get_min().x(), bb.get_min().y();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-min, higher than y-max
      pt << -3, 3;
      corner << bb.get_min().x(), bb.get_max().y();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test higher than x-max, lower than y-min
      pt << 3, -3;
      corner << bb.get_max().x(), bb.get_min().y();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test higher than x-max, higher than y-max
      pt << 3, 3;
      corner << bb.get_max().x(), bb.get_max().y();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-min, inside y
      pt << -3, 0;
      dist_ref=std::abs(pt.x()-bb.get_min().x());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test higher than x-max, inside y
      pt << 3, 0;
      dist_ref=std::abs(pt.x()-bb.get_max().x());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside x, lower than y-min
      pt << 0, -4;
      dist_ref=std::abs(pt.y()-bb.get_min().y());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside x, higher than y-max
      pt << 0, 4;
      dist_ref=std::abs(pt.y()-bb.get_max().y());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside
      pt << 0, 0;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
      pt << 0.5, 0.5;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test on point
      pt << 1, -1;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
      pt << -1, -1;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
    }

    void point_3d_test()
    {
      bounding_box_type3 bb;
      point_type3 pt, corner;
      data_type dist, dist_ref;

      pt << -1, -1, -1;
      bb.add(pt);
      pt << 1, 1, 1;
      bb.add(pt);

      // test lower than x-min, lower than y-min, lower than z-min
      pt << -3, -3, -3;
      corner << bb.get_min().x(), bb.get_min().y(), bb.get_min().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-min, higher than y-max, lower than z-min
      pt << -3, 3, -3;
      corner << bb.get_min().x(), bb.get_max().y(), bb.get_min().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test higher than x-max, lower than y-min, lower than z-min
      pt << 3, -3, -3;
      corner << bb.get_max().x(), bb.get_min().y(), bb.get_min().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test higher than x-max, lower than y-max, lower than z-min
      pt << 3, 3, -3;
      corner << bb.get_max().x(), bb.get_max().y(), bb.get_min().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-min, lower than y-min, higher than z-max
      pt << -3, -3, 3;
      corner << bb.get_min().x(), bb.get_min().y(), bb.get_max().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-min, higher than y-max, higher than z-max
      pt << -3, 3, 3;
      corner << bb.get_min().x(), bb.get_max().y(), bb.get_max().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test higher than x-max, lower than y-min, higher than z-max
      pt << 3, -3, 3;
      corner << bb.get_max().x(), bb.get_min().y(), bb.get_max().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test higher than x-max, lower than y-max, lower than z-max
      pt << 3, 3, 3;
      corner << bb.get_max().x(), bb.get_max().y(), bb.get_max().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside x, lower than y-min, lower than z-min
      pt << 0, -3, -3;
      corner << 0, bb.get_min().y(), bb.get_min().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside x, lower than y-max, lower than z-min
      pt << 0, 3, -3;
      corner << 0, bb.get_max().y(), bb.get_min().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside x, lower than y-min, lower than z-max
      pt << 0, -3, 3;
      corner << 0, bb.get_min().y(), bb.get_max().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside x, lower than y-max, lower than z-max
      pt << 0, 3, 3;
      corner << 0, bb.get_max().y(), bb.get_max().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-min, inside y, lower than z-min
      pt << -3, 0, -3;
      corner << bb.get_min().x(), 0, bb.get_min().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-min, inside y, lower than z-max
      pt << -3, 0, 3;
      corner << bb.get_min().x(), 0, bb.get_max().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-max, inside y, lower than z-min
      pt << 3, 0, -3;
      corner << bb.get_max().x(), 0, bb.get_min().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-max, inside y, lower than z-max
      pt << 3, 0, 3;
      corner << bb.get_max().x(), 0, bb.get_max().z();
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-min, lower than y-min, inside z
      pt <<-3, -3, 0;
      corner <<bb.get_min().x(), bb.get_min().y(), 0;
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-max, lower than y-min, inside z
      pt <<3, -3, 0;
      corner <<bb.get_max().x(), bb.get_min().y(), 0;
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-min, lower than y-max, inside z
      pt <<-3, 3, 0;
      corner <<bb.get_min().x(), bb.get_max().y(), 0;
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-max, lower than y-max, inside z
      pt << 3, 3, 0;
      corner <<bb.get_max().x(), bb.get_max().y(), 0;
      dist_ref=eli::geom::point::distance(pt, corner);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test lower than x-min, inside y, inside z
      pt << -3, 0, 0;
      dist_ref=std::abs(pt.x()-bb.get_min().x());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test higher than x-max, inside y, inside z
      pt << 3, 0, 0;
      dist_ref=std::abs(pt.x()-bb.get_max().x());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside x, lower than y-min, inside z
      pt << 0, -3, 0;
      dist_ref=std::abs(pt.y()-bb.get_min().y());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside x, higher than y-max, inside z
      pt << 0, 3, 0;
      dist_ref=std::abs(pt.y()-bb.get_max().y());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside x, inside y, lower than z-min
      pt << 0, 0, -3;
      dist_ref=std::abs(pt.z()-bb.get_min().z());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside x, inside y, higher than z-max
      pt << 0, 0, 3;
      dist_ref=std::abs(pt.z()-bb.get_max().z());
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test inside
      pt << 0, 0, 0;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
      pt << 0.5, 0.5, 0.5;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
      pt <<-0.5, 0.25, 0;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));

      // test on point
      pt << 1, -1, 1;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
      pt << -1, -1, -1;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
      pt << 1, 0.5, -1;
      dist_ref=static_cast<data_type>(0);
      dist=eli::geom::intersect::minimum_distance(bb, pt);
      TEST_ASSERT(tol.approximately_equal(dist, dist_ref));
    }
};

#endif

