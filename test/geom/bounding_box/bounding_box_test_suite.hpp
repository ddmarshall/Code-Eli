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

#ifndef bounding_box_test_suite_hpp
#define bounding_box_test_suite_hpp

#include <cmath>    // cos(), sin()

#include <typeinfo> // typeid

#include "eli/geom/general/bounding_box.hpp"

template<typename data__>
class bounding_box_test_suite : public Test::Suite
{
  private:
    typedef data__ data_type;
    typedef eli::geom::general::bounding_box<data__, 1> bounding_box_type1;
    typedef typename bounding_box_type1::point_type point_type1;
    typedef eli::geom::general::bounding_box<data__, 2> bounding_box_type2;
    typedef typename bounding_box_type2::point_type point_type2;
    typedef eli::geom::general::bounding_box<data__, 3> bounding_box_type3;
    typedef typename bounding_box_type3::point_type point_type3;

    typedef typename bounding_box_type1::index_type index_type;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(bounding_box_test_suite<float>::construction_test);
      TEST_ADD(bounding_box_test_suite<float>::add_test);
      TEST_ADD(bounding_box_test_suite<float>::add_bb_test);
      TEST_ADD(bounding_box_test_suite<float>::inside_test);
      TEST_ADD(bounding_box_test_suite<float>::intersect_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(bounding_box_test_suite<double>::construction_test);
      TEST_ADD(bounding_box_test_suite<double>::add_test);
      TEST_ADD(bounding_box_test_suite<double>::add_bb_test);
      TEST_ADD(bounding_box_test_suite<double>::intersect_test);
      TEST_ADD(bounding_box_test_suite<double>::inside_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(bounding_box_test_suite<long double>::construction_test);
      TEST_ADD(bounding_box_test_suite<long double>::add_test);
      TEST_ADD(bounding_box_test_suite<long double>::add_bb_test);
      TEST_ADD(bounding_box_test_suite<long double>::inside_test);
      TEST_ADD(bounding_box_test_suite<long double>::intersect_test);
    }

  public:
    bounding_box_test_suite()
    {
      AddTests(data__());
    }
    ~bounding_box_test_suite()
    {
    }

  private:
    void construction_test()
    {
      // test 1D
      {
        bounding_box_type1 bb1, bb2;

        TEST_ASSERT(bb1.empty_set());

        bb2=bb1;
        TEST_ASSERT(bb2==bb1);

        point_type1 pt1;
        pt1 << 1;

        bounding_box_type1 bb3(pt1);

        TEST_ASSERT(!bb3.empty_set());
        TEST_ASSERT(bb3.get_min()==pt1);
        TEST_ASSERT(bb3.get_max()==pt1);

        bounding_box_type1 bb4(bb1);

        TEST_ASSERT(bb1==bb4);
      }

      // test 2D
      {
        bounding_box_type2 bb1, bb2;

        TEST_ASSERT(bb1.empty_set());

        bb2=bb1;
        TEST_ASSERT(bb2==bb1);

        point_type2 pt1;
        pt1 << 1, 2;

        bounding_box_type2 bb3(pt1);

        TEST_ASSERT(!bb3.empty_set());
        TEST_ASSERT(bb3.get_min()==pt1);
        TEST_ASSERT(bb3.get_max()==pt1);

        bounding_box_type2 bb4(bb1);

        TEST_ASSERT(bb1==bb4);
      }

      // test 3D
      {
        bounding_box_type3 bb1, bb2;

        TEST_ASSERT(bb1.empty_set());

        bb2=bb1;
        TEST_ASSERT(bb2==bb1);

        point_type3 pt1;
        pt1 << 1, 2,  3;

        bounding_box_type3 bb3(pt1);

        TEST_ASSERT(!bb3.empty_set());
        TEST_ASSERT(bb3.get_min()==pt1);
        TEST_ASSERT(bb3.get_max()==pt1);

        bounding_box_type3 bb4(bb1);

        TEST_ASSERT(bb1==bb4);
      }
    }

    void add_test()
    {
      // test 1D
      {
        bounding_box_type1 bb;
        point_type1 pt1, pt2, pt3, pt4, pt5;
        bool changed;

        pt1 << 1;
        pt2 << 3;
        pt3 << 4;
        pt4 << -1;
        pt5 << 0;

        TEST_ASSERT(bb.empty_set());

        changed=bb.add(pt1);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min()==pt1);
        TEST_ASSERT(bb.get_max()==pt1);

        changed=bb.add(pt2);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min()==pt1);
        TEST_ASSERT(bb.get_max()==pt2);

        changed=bb.add(pt3);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min()==pt1);
        TEST_ASSERT(bb.get_max()==pt3);

        changed=bb.add(pt4);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min()==pt4);
        TEST_ASSERT(bb.get_max()==pt3);

        changed=bb.add(pt5);
        TEST_ASSERT(!changed);
        TEST_ASSERT(bb.get_min()==pt4);
        TEST_ASSERT(bb.get_max()==pt3);
      }

      // test 2D
      {
        bounding_box_type2 bb;
        point_type2 pt1, pt2, pt3, pt4, pt5;
        bool changed;

        pt1 << 1, 2;
        pt2 << 3, 4;
        pt3 << 4, 1;
        pt4 << -1, 5;
        pt5 << 0, 3;

        TEST_ASSERT(bb.empty_set());

        changed=bb.add(pt1);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min()==pt1);
        TEST_ASSERT(bb.get_max()==pt1);

        changed=bb.add(pt2);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min()==pt1);
        TEST_ASSERT(bb.get_max()==pt2);

        changed=bb.add(pt3);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min().x()==pt1.x());
        TEST_ASSERT(bb.get_min().y()==pt3.y());
        TEST_ASSERT(bb.get_max().x()==pt3.x());
        TEST_ASSERT(bb.get_max().y()==pt2.y());

        changed=bb.add(pt4);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min().x()==pt4.x());
        TEST_ASSERT(bb.get_min().y()==pt3.y());
        TEST_ASSERT(bb.get_max().x()==pt3.x());
        TEST_ASSERT(bb.get_max().y()==pt4.y());

        changed=bb.add(pt5);
        TEST_ASSERT(!changed);
        TEST_ASSERT(bb.get_min().x()==pt4.x());
        TEST_ASSERT(bb.get_min().y()==pt3.y());
        TEST_ASSERT(bb.get_max().x()==pt3.x());
        TEST_ASSERT(bb.get_max().y()==pt4.y());
      }

      // test 3D
      {
        bounding_box_type3 bb;
        point_type3 pt1, pt2, pt3, pt4, pt5;
        bool changed;

        pt1 << 1, 2, 1;
        pt2 << 3, 4, 2;
        pt3 << 4, 1, -1;
        pt4 << -1, 5, 4;
        pt5 << 0, 3, 0;

        TEST_ASSERT(bb.empty_set());

        changed=bb.add(pt1);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min()==pt1);
        TEST_ASSERT(bb.get_max()==pt1);

        changed=bb.add(pt2);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min()==pt1);
        TEST_ASSERT(bb.get_max()==pt2);

        changed=bb.add(pt3);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min().x()==pt1.x());
        TEST_ASSERT(bb.get_min().y()==pt3.y());
        TEST_ASSERT(bb.get_min().z()==pt3.z());
        TEST_ASSERT(bb.get_max().x()==pt3.x());
        TEST_ASSERT(bb.get_max().y()==pt2.y());
        TEST_ASSERT(bb.get_max().z()==pt2.z());

        changed=bb.add(pt4);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb.get_min().x()==pt4.x());
        TEST_ASSERT(bb.get_min().y()==pt3.y());
        TEST_ASSERT(bb.get_min().z()==pt3.z());
        TEST_ASSERT(bb.get_max().x()==pt3.x());
        TEST_ASSERT(bb.get_max().y()==pt4.y());
        TEST_ASSERT(bb.get_max().z()==pt4.z());

        changed=bb.add(pt5);
        TEST_ASSERT(!changed);
        TEST_ASSERT(bb.get_min().x()==pt4.x());
        TEST_ASSERT(bb.get_min().y()==pt3.y());
        TEST_ASSERT(bb.get_min().z()==pt3.z());
        TEST_ASSERT(bb.get_max().x()==pt3.x());
        TEST_ASSERT(bb.get_max().y()==pt4.y());
        TEST_ASSERT(bb.get_max().z()==pt4.z());
      }
    }

    void add_bb_test()
    {
      // test 1D
      {
        bounding_box_type1 bb1, bb2, bb3;
        point_type1 pt1, pt2, pt3, pt4, pt5, pt6;
        bool changed;

        pt1 << 1;
        pt2 << 3;
        pt3 << 4;
        pt4 << -1;
        pt5 << 0;
        pt6 << 1;

        bb1.add(pt1);
        bb1.add(pt2);

        bb2.add(pt3);
        bb2.add(pt4);

        bb3.add(pt5);
        bb3.add(pt6);

        // test add that changes bbox
        changed=bb1.add(bb2);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb1.get_min()==pt4);
        TEST_ASSERT(bb1.get_max()==pt3);

        // test add that does change bbox
        changed=bb1.add(bb3);
        TEST_ASSERT(!changed);
        TEST_ASSERT(bb1.get_min()==pt4);
        TEST_ASSERT(bb1.get_max()==pt3);
      }

      // test 2D
      {
        bounding_box_type2 bb1, bb2;
        point_type2 pt1, pt2, pt3, pt4, pt5;
        bool changed;

        pt1 << 1, 2;
        pt2 << 3, 4;
        pt3 << 4, 1;
        pt4 << -1, 5;
        pt5 << 0, 3;

        bb1.add(pt1);
        bb1.add(pt2);

        bb2.add(pt3);
        bb2.add(pt4);

        changed=bb1.add(bb2);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb1.get_min().x()==pt4.x());
        TEST_ASSERT(bb1.get_min().y()==pt3.y());
        TEST_ASSERT(bb1.get_max().x()==pt3.x());
        TEST_ASSERT(bb1.get_max().y()==pt4.y());
      }

      // test 3D
      {
        bounding_box_type3 bb1, bb2;
        point_type3 pt1, pt2, pt3, pt4, pt5;
        bool changed;

        pt1 << 1, 2, 1;
        pt2 << 3, 4, 2;
        pt3 << 4, 1, -1;
        pt4 << -1, 5, 4;
        pt5 << 0, 3, 0;

        bb1.add(pt1);
        bb1.add(pt2);

        bb2.add(pt3);
        bb2.add(pt4);

        changed=bb1.add(bb2);
        TEST_ASSERT(changed);
        TEST_ASSERT(bb1.get_min().x()==pt4.x());
        TEST_ASSERT(bb1.get_min().y()==pt3.y());
        TEST_ASSERT(bb1.get_min().z()==pt3.z());
        TEST_ASSERT(bb1.get_max().x()==pt3.x());
        TEST_ASSERT(bb1.get_max().y()==pt4.y());
        TEST_ASSERT(bb1.get_max().z()==pt4.z());
      }
    }

    void inside_test()
    {
      // test 1D
      {
        bounding_box_type1 bb;
        point_type1 pt;

        pt << -1;
        bb.set_min(pt);
        TEST_ASSERT(bb.get_min()==pt);
        TEST_ASSERT(bb.get_max()==pt);

        pt << 1;
        bb.set_max(pt);
        TEST_ASSERT(bb.get_max()==pt);

        // test case where there is intersection
        pt << static_cast<data_type>(0.9);
        TEST_ASSERT(bb.inside(pt));

        // test case where there is no intersection
        pt << static_cast<data_type>(1.1);
        TEST_ASSERT(!bb.inside(pt));

        // test case where there is contact intersection
        pt << static_cast<data_type>(1);
        TEST_ASSERT(bb.inside(pt));
      }

      // test 2D
      {
        bounding_box_type2 bb;
        point_type2 pt;

        pt << -1, -1;
        bb.set_min(pt);
        TEST_ASSERT(bb.get_min()==pt);
        TEST_ASSERT(bb.get_max()==pt);

        pt << 1, 1;
        bb.set_max(pt);
        TEST_ASSERT(bb.get_max()==pt);

        // test case where there is intersection
        pt << static_cast<data_type>(0.9), static_cast<data_type>(0.9);
        TEST_ASSERT(bb.inside(pt));

        // test case where there is no intersection
        pt << static_cast<data_type>(1.1), static_cast<data_type>(1.1);
        TEST_ASSERT(!bb.inside(pt));

        // test case where there is contact intersection
        pt << static_cast<data_type>(1), static_cast<data_type>(0);
        TEST_ASSERT(bb.inside(pt));
      }

      // test 3D
      {
        bounding_box_type3 bb;
        point_type3 pt;

        pt << -1, -1, -1;
        bb.set_min(pt);
        TEST_ASSERT(bb.get_min()==pt);
        TEST_ASSERT(bb.get_max()==pt);

        pt << 1, 1, 1;
        bb.set_max(pt);
        TEST_ASSERT(bb.get_max()==pt);

        // test case where there is intersection
        pt << static_cast<data_type>(0.9), static_cast<data_type>(0.9), static_cast<data_type>(0.9);
        TEST_ASSERT(bb.inside(pt));

        // test case where there is no intersection
        pt << static_cast<data_type>(1.1), static_cast<data_type>(1.1), static_cast<data_type>(1.1);
        TEST_ASSERT(!bb.inside(pt));

        // test case where there is contact intersection
        pt << static_cast<data_type>(1), static_cast<data_type>(0), static_cast<data_type>(0);
        TEST_ASSERT(bb.inside(pt));
      }
    }

    void intersect_test()
    {
      // test 1D
      {
        bounding_box_type1 bb1, bb2;
        point_type1 pt;

        pt << -1;
        bb1.set_min(pt);
        pt << 1;
        bb1.set_max(pt);

        // test case where simple intersection
        pt << static_cast<data_type>(0);
        bb2.set_min(pt);
        pt << static_cast<data_type>(2);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));

        // test case where no intersection
        pt << static_cast<data_type>(-3);
        bb2.set_min(pt);
        pt << static_cast<data_type>(-2);
        bb2.set_max(pt);
        TEST_ASSERT(!bb1.intersect(bb2));

        // test case where contact intersection
        pt << static_cast<data_type>(1);
        bb2.set_min(pt);
        pt << static_cast<data_type>(2);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));

        // test case where one contains other
        pt << static_cast<data_type>(0);
        bb2.set_min(pt);
        pt << static_cast<data_type>(0.5);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));

        // test case where one is contained by other
        pt << static_cast<data_type>(-2);
        bb2.set_min(pt);
        pt << static_cast<data_type>(2);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));
      }

      // test 2D
      {
        bounding_box_type2 bb1, bb2;
        point_type2 pt;

        pt << -1, -1;
        bb1.set_min(pt);
        pt << 1, 1;
        bb1.set_max(pt);

        // test case where simple intersection
        pt << static_cast<data_type>(0), static_cast<data_type>(0);
        bb2.set_min(pt);
        pt << static_cast<data_type>(2), static_cast<data_type>(2);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));

        // test case where no intersection
        pt << static_cast<data_type>(-3), static_cast<data_type>(-4);
        bb2.set_min(pt);
        pt << static_cast<data_type>(-2), static_cast<data_type>(0);
        bb2.set_max(pt);
        TEST_ASSERT(!bb1.intersect(bb2));

        // test case where contact intersection
        pt << static_cast<data_type>(1), static_cast<data_type>(0);
        bb2.set_min(pt);
        pt << static_cast<data_type>(2), static_cast<data_type>(2);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));

        // test case where one contains other
        pt << static_cast<data_type>(0), static_cast<data_type>(0);
        bb2.set_min(pt);
        pt << static_cast<data_type>(0.5), static_cast<data_type>(0.5);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));

        // test case where one is contained by other
        pt << static_cast<data_type>(-2), static_cast<data_type>(-2);
        bb2.set_min(pt);
        pt << static_cast<data_type>(2), static_cast<data_type>(2);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));
      }

      // test 3D
      {
        bounding_box_type3 bb1, bb2;
        point_type3 pt;

        pt << -1, -1, -1;
        bb1.set_min(pt);
        pt << 1, 1, 1;
        bb1.set_max(pt);

        // test case where simple intersection
        pt << static_cast<data_type>(0), static_cast<data_type>(0), static_cast<data_type>(0);
        bb2.set_min(pt);
        pt << static_cast<data_type>(2), static_cast<data_type>(2), static_cast<data_type>(2);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));

        // test case where no intersection
        pt << static_cast<data_type>(-3), static_cast<data_type>(-4), static_cast<data_type>(-3);
        bb2.set_min(pt);
        pt << static_cast<data_type>(-2), static_cast<data_type>(0), static_cast<data_type>(0);
        bb2.set_max(pt);
        TEST_ASSERT(!bb1.intersect(bb2));

        // test case where contact intersection
        pt << static_cast<data_type>(1), static_cast<data_type>(0), static_cast<data_type>(0);
        bb2.set_min(pt);
        pt << static_cast<data_type>(2), static_cast<data_type>(2), static_cast<data_type>(2);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));

        // test case where one contains other
        pt << static_cast<data_type>(0), static_cast<data_type>(0), static_cast<data_type>(0);
        bb2.set_min(pt);
        pt << static_cast<data_type>(0.5), static_cast<data_type>(0.5), static_cast<data_type>(0.5);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));

        // test case where one is contained by other
        pt << static_cast<data_type>(-2), static_cast<data_type>(-2), static_cast<data_type>(-2);
        bb2.set_min(pt);
        pt << static_cast<data_type>(2), static_cast<data_type>(2), static_cast<data_type>(2);
        bb2.set_max(pt);
        TEST_ASSERT(bb1.intersect(bb2));
      }
    }
};

#endif

