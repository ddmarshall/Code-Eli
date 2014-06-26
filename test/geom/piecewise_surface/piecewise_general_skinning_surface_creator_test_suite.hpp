/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_curve_piecewise_line_creator_hpp
#define eli_geom_curve_piecewise_line_creator_hpp

#include <vector>

#include "Eigen/Eigen"

#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_line_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;

          piecewise_line_creator() : piecewise_creator_base<data_type, dim__, tolerance_type>(2, 0), pt(2) {}
          piecewise_line_creator(const index_type &ns) : piecewise_creator_base<data_type, dim__, tolerance_type>(ns, 0), pt(ns+1) {}
          piecewise_line_creator(const piecewise_polygon_creator<data_type, dim__, tolerance_type> &ppc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(ppc), pt(ppc.corner) {}
          ~piecewise_line_creator() {}

          void set_point(const point_type &c, const index_type &i)
          {
            if ((i>=0) && (i<static_cast<index_type>(pt.size())))
              pt[i]=c;
            else
              assert(false);
          }
          point_type get_point(const index_type &i) const
          {
            if ((i<0) || (i>=static_cast<index_type>(pt.size())))
            {
              return pt[0];
              assert(false);
            }
            return pt[i];
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;

            pc.clear();

            curve_type c(1);
            error_code err;
            index_type nsegs(this->get_number_segments());

            // do sanity check
            if (pt.size()!=static_cast<size_t>(nsegs+1))
            {
              assert(false);
              return false;
            }

            // set the start parameter
            pc.set_t0(this->get_t0());

            // set the lines
            for (index_type i=0; i<nsegs; ++i)
            {
              c.set_control_point(pt[i], 0);
              c.set_control_point(pt[i+1], 1);
              err=pc.push_back(c, this->get_segment_dt(i));
              if (err!=piecewise_curve_type::NO_ERRORS)
              {
                pc.clear();
                pc.set_t0(0);
                return false;
              }
            }

            return true;
          }

        private:
          void number_segments_changed() {pt.resize(this->get_number_segments()+1);}

        private:
          std::vector<point_type, Eigen::aligned_allocator<point_type>> pt;
      };
    }
  }
}
#endif

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

#ifndef piecewise_general_skinning_surface_creator_test_suite_suite_hpp
#define piecewise_general_skinning_surface_creator_test_suite_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"

#include "eli/geom/curve/piecewise.hpp"
// TODO: Add this to type of piecewise curve creators
//#include "eli/geom/curve/piecewise_line_creator.hpp"

#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/piecewise_general_skinning_surface_creator.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits
#include <iterator> // std::insert_iterator

template<typename data__>
class piecewise_general_skinning_surface_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::surface::piecewise<eli::geom::surface::bezier, data__, 3> piecewise_surface_type;
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_surface_type::surface_type surface_type;
    typedef typename piecewise_surface_type::point_type point_type;
    typedef typename piecewise_surface_type::data_type data_type;
    typedef typename piecewise_surface_type::index_type index_type;
    typedef typename piecewise_surface_type::tolerance_type tolerance_type;
    typedef typename eli::geom::surface::connection_data<data__, 3, tolerance_type> rib_data_type;
    typedef typename rib_data_type::curve_type rib_curve_type;
    typedef typename eli::geom::curve::piecewise_line_creator<data__, 3, tolerance_type> piecewise_line_creator_type;
    typedef typename eli::geom::surface::general_skinning_surface_creator<data__, 3, tolerance_type> general_creator_type;

    tolerance_type tol;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<float>::create_rib_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<float>::create_single_surface_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<double>::create_rib_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<double>::create_single_surface_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<long double>::create_rib_test);
      TEST_ADD(piecewise_general_skinning_surface_creator_test_suite<long double>::create_single_surface_test);
    }

  public:
    piecewise_general_skinning_surface_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_general_skinning_surface_creator_test_suite()
    {
    }

  private:
    void create_rib_test()
    {
      rib_curve_type rc1, rc2, rc3;
      point_type p1, p2, p3;
      bool rtn_flag;

      p1 << 1, 0, 0;
      p2 << 2, 2, 2;
      p3 << 3, 2, 1;

      // create three rib curves
      piecewise_line_creator_type plc(1);

      plc.set_point(p1, 0);
      plc.set_point(p2, 1);
      plc.create(rc1);

      plc.set_point(p1, 0);
      plc.set_point(p3, 1);
      plc.create(rc2);

      plc.set_point(p3, 0);
      plc.set_point(p2, 1);
      plc.create(rc3);

      // set rib
      {
        rib_data_type rib;

        // should start in bad state
        rtn_flag=rib.check_state();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(!rib.use_f());

        // add rib
        rtn_flag=rib.set_f(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());

        // unset rib
        rtn_flag=rib.unset_f();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(!rib.use_f());
      }

      // set rib and fp
      {
        rib_data_type rib;

        // should start in bad state
        rtn_flag=rib.check_state();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(!rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());

        // add left fp
        rtn_flag=rib.set_left_fp(rc1);
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_left_fp()==rc1);
        TEST_ASSERT(!rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());

        // add rib
        rtn_flag=rib.set_f(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());

        // add right fp
        rtn_flag=rib.set_right_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_right_fp()==rc2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());

        // remove left fp
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());

        // remove left fp again
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());

        // remove right fp
        rtn_flag=rib.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());

        // add fp
        rtn_flag=rib.set_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_left_fp()==rc2);
        TEST_ASSERT(rib.get_right_fp()==rc2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());

        // remove right fp
        rtn_flag=rib.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
      }

      // set rib, fp and fpp
      {
        rib_data_type rib;

        // should start in bad state
        rtn_flag=rib.check_state();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(!rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add left fpp
        rtn_flag=rib.set_left_fpp(rc1);
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_left_fpp()==rc1);
        TEST_ASSERT(!rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add right fpp
        rtn_flag=rib.set_right_fpp(rc2);
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_right_fpp()==rc2);
        TEST_ASSERT(!rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // add rib
        rtn_flag=rib.set_f(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // remove right fpp
        rtn_flag=rib.unset_right_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fp());

        // remove right fpp again
        rtn_flag=rib.unset_right_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fp());

        // remove left fpp
        rtn_flag=rib.unset_left_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add fp
        rtn_flag=rib.set_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_left_fp()==rc2);
        TEST_ASSERT(rib.get_right_fp()==rc2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add fpp
        rtn_flag=rib.set_fpp(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_left_fpp()==rc1);
        TEST_ASSERT(rib.get_right_fpp()==rc1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // remove right fp
        rtn_flag=rib.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // remove right fpp
        rtn_flag=rib.unset_right_fpp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // remove everything
        rtn_flag=rib.unset_fpp();
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.unset_fp();
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.unset_f();
        TEST_ASSERT(!rtn_flag);

        // add everything to right only
        rtn_flag=rib.set_f(rc2);
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.set_right_fp(rc1);
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.set_right_fpp(rc3);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.get_f()==rc2);
        TEST_ASSERT(rib.get_right_fp()==rc1);
        TEST_ASSERT(rib.get_right_fpp()==rc3);
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());
      }

      // set conditions with C1
      {
        rib_data_type rib;

        // start with point
        rtn_flag=rib.set_f(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // set continuity to C1
        rtn_flag=rib.set_continuity(rib_data_type::C1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add left fp
        rtn_flag=rib.set_left_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // unset left fp
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        rtn_flag=rib.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add right fp
        rtn_flag=rib.set_right_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // change continuity
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        rtn_flag=rib.set_continuity(rib_data_type::C0);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // add fp
        rtn_flag=rib.set_fp(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C0);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // change continuity
        rtn_flag=rib.set_continuity(rib_data_type::C1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // unset left fp
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // set left fp to change left and right value
        rtn_flag=rib.set_left_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_left_fp()==rc2);
        TEST_ASSERT(rib.get_right_fp()==rc2);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // set left fp
        rtn_flag=rib.set_left_fp(rib.get_right_fp());
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());
      }

      // set point, fp, fpp and C2
      {
        rib_data_type rib;

        // start with point and fp
        rtn_flag=rib.set_f(rc1);
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.set_fp(rc2);
        TEST_ASSERT(rtn_flag);
        rtn_flag=rib.set_continuity(rib_data_type::C1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C1);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // set continuity to C2
        rtn_flag=rib.set_continuity(rib_data_type::C2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // unset left fp
        rtn_flag=rib.unset_left_fp();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // unset right fp
        rtn_flag=rib.unset_right_fp();
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(!rib.use_right_fpp());

        // set fpp
        rtn_flag=rib.set_fpp(rc1);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // unset left fpp
        rtn_flag=rib.unset_left_fpp();
        TEST_ASSERT(!rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(!rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // set left fpp to set both left and right values
        rtn_flag=rib.set_left_fpp(rc3);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // set left fpp
        rtn_flag=rib.set_left_fpp(rib.get_right_fpp());
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(!rib.use_left_fp());
        TEST_ASSERT(!rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());

        // set fp
        rtn_flag=rib.set_fp(rc2);
        TEST_ASSERT(rtn_flag);
        TEST_ASSERT(rib.get_continuity()==rib_data_type::C2);
        TEST_ASSERT(rib.use_f());
        TEST_ASSERT(rib.use_left_fp());
        TEST_ASSERT(rib.use_right_fp());
        TEST_ASSERT(rib.use_left_fpp());
        TEST_ASSERT(rib.use_right_fpp());
      }

      // copy constructor, equivalence and assignment operators
      {
        rib_data_type r1, r2, r3;

        // defaults should be the same
        TEST_ASSERT(r1==r2);

        // manually build identical ones
        r1.set_f(rc2);
        r1.set_left_fp(rc1);
        r1.set_right_fpp(rc3);
        r2.set_f(rc2);
        r2.set_left_fp(rc1);
        r2.set_right_fpp(rc3);
        TEST_ASSERT(r1==r2);

        // assignment operator
        r3=r2;
        TEST_ASSERT(r3==r1);
        TEST_ASSERT(r3==r2);

        // set value for r2's fpp and then don't use it
        r1.unset_fpp();
        r2.unset_fpp();
        r3.unset_fpp();
        r2.set_fpp(rc1);
        r2.unset_fpp();
        TEST_ASSERT(r3==r1);
        TEST_ASSERT(r3==r2);

        // copy constructor
        rib_data_type r4(r3);
        TEST_ASSERT(r4==r3);
      }

      // test the joint interface
      {
        rib_data_type r1;
        data_type s1, s2, s3, s4, s5;

        s1=static_cast<data_type>(0.25);
        s2=static_cast<data_type>(0.35);
        s3=static_cast<data_type>(0.45);
        s4=static_cast<data_type>(0.55);
        s5=static_cast<data_type>(0.70);
        rc1.split(s1);
        rc1.split(s4);
        rc2.split(s2);
        rc2.split(s3);
        rc2.split(s4);
        rc3.split(s3);
        rc3.split(s5);

        r1.set_f(rc1);
        r1.set_fp(rc2);
        r1.set_fpp(rc3);

        // get the joints from the rib
        std::vector<data_type> rjoints;
        r1.get_joints(std::back_inserter(rjoints));

        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }

        // split the ribs at the specified locations
        std::vector<index_type> r1degs;
        r1.split(rjoints.begin(), rjoints.end(), std::back_inserter(r1degs));

        // check the degree
        TEST_ASSERT(r1degs.size()==6);
        if (r1degs.size()==6)
        {
          TEST_ASSERT(r1degs[0]==1);
          TEST_ASSERT(r1degs[1]==1);
          TEST_ASSERT(r1degs[2]==1);
          TEST_ASSERT(r1degs[3]==1);
          TEST_ASSERT(r1degs[4]==1);
          TEST_ASSERT(r1degs[5]==1);
        }

        // check the joint locations
        rjoints.clear();
        r1.get_f().get_parameters(std::back_inserter(rjoints));
        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }
        rjoints.clear();
        r1.get_left_fp().get_parameters(std::back_inserter(rjoints));
        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }
        rjoints.clear();
        r1.get_right_fp().get_parameters(std::back_inserter(rjoints));
        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }
        rjoints.clear();
        r1.get_left_fpp().get_parameters(std::back_inserter(rjoints));
        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }
        rjoints.clear();
        r1.get_right_fpp().get_parameters(std::back_inserter(rjoints));
        TEST_ASSERT(rjoints.size()==7);
        if (rjoints.size()==7)
        {
          TEST_ASSERT(tol.approximately_equal(rc1.get_t0(),   rjoints[0]));
          TEST_ASSERT(tol.approximately_equal(s1,             rjoints[1]));
          TEST_ASSERT(tol.approximately_equal(s2,             rjoints[2]));
          TEST_ASSERT(tol.approximately_equal(s3,             rjoints[3]));
          TEST_ASSERT(tol.approximately_equal(s4,             rjoints[4]));
          TEST_ASSERT(tol.approximately_equal(s5,             rjoints[5]));
          TEST_ASSERT(tol.approximately_equal(rc1.get_tmax(), rjoints[6]));
        }
      }
    }

    void create_single_surface_test()
    {
      // simple surface connecting 2 lines
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11;
        data_type u0(2), v0(1), u1(4), v1(5);
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;

        // create two rib curves
        piecewise_line_creator_type plc(1);

        plc.set_point(p00, 0);
        plc.set_point(p01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc1);

        plc.set_point(p10, 0);
        plc.set_point(p11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[1].set_f(rc2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_v0(v0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(float)))
//        {
//          std::cout.flush();
//          eli::octave_start(1);
//          eli::octave_print(1, s, "surf", true);
//          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
//          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
//          eli::octave_finish(1);
//        }
      }

      // simple surface connecting 2 lines
      {
        index_type nsegs(1);
        std::vector<rib_data_type> ribs(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        rib_curve_type rc1, rc2, rs1, rs2;
        general_creator_type gc;
        piecewise_surface_type s;
        point_type p00, p01, p10, p11, s00, s01, s10, s11;
        data_type u0(2), v0(1), u1(4), v1(5), v;
        bool rtn_flag;

        // four corners of surface
        p00 << 1, 1, 1;
        p01 << 2, 3, 1;
        p10 << 3, 2, 3;
        p11 << 4, 4, 4;
        s00 << 2, 2, 0;
        s01 << 1, 3, 2;
        s10 << 2, 1, 0;
        s11 << 2, 0, 1;

        // create two rib curves && two slope curves
        piecewise_line_creator_type plc(1);

        plc.set_point(p00, 0);
        plc.set_point(p01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc1);

        plc.set_point(s00, 0);
        plc.set_point(s01, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs1);

        plc.set_point(p10, 0);
        plc.set_point(p11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rc2);

        plc.set_point(s10, 0);
        plc.set_point(s11, 1);
        plc.set_t0(v0);
        plc.set_segment_dt(v1-v0, 0);
        plc.create(rs2);

        // set the rib data
        ribs[0].set_f(rc1);
        ribs[0].set_right_fp(rs1);
        ribs[1].set_f(rc2);
        ribs[1].set_left_fp(rs2);

        // set the maximum degrees of each segment
        max_degree[0]=0;

        // create surface
        rtn_flag=gc.set_conditions(ribs, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_u0(u0);
        gc.set_v0(v0);
        gc.set_segment_du(u1-u0, 0);
        rtn_flag=gc.create(s);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test, p_ref;

        p_test=s.f(u0, v0);
        TEST_ASSERT(tol.approximately_equal(p00, p_test));
//        std::cout << "p=" << p00
//                  << "\tpc=" << ribs[0].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (ribs[0].get_f().f(v0)-p_test).norm() << std::endl;
        p_test=s.f(u1, v0);
        TEST_ASSERT(tol.approximately_equal(p10, p_test));
//        std::cout << "p=" << p10
//                  << "\tpc=" << ribs[1].get_f().f(v0)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p10-p_test).norm() << std::endl;
        p_test=s.f(u0, v1);
        TEST_ASSERT(tol.approximately_equal(p01, p_test));
//        std::cout << "p=" << p01
//                  << "\tpc=" << ribs[0].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p01-p_test).norm() << std::endl;
        p_test=s.f(u1, v1);
        TEST_ASSERT(tol.approximately_equal(p11, p_test));
//        std::cout << "p=" << p11
//                  << "\tpc=" << ribs[1].get_f().f(v1)
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p11-p_test).norm() << std::endl;
        v=(v0+v1)/2;
        p_test=s.f_u(u0, v);
        p_ref=ribs[0].get_right_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;
        p_test=s.f_u(u1, v);
        p_ref=ribs[1].get_left_fp().f(v);
        TEST_ASSERT(tol.approximately_equal(p_test, p_ref));
//        std::cout << "p=" << p_ref
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (p_ref-p_test).norm() << std::endl;

        if (rtn_flag && (typeid(data_type)==typeid(float)))
        {
          std::cout.flush();
          eli::octave_start(1);
          eli::octave_print(1, s, "surf", false);
          eli::octave_print(1, ribs[0].get_f(), "rib0", true);
          eli::octave_print(1, ribs[1].get_f(), "rib1", true);
          eli::octave_print(1, ribs[0].get_f(), ribs[0].get_right_fp(), "rve0");
          eli::octave_print(1, ribs[1].get_f(), ribs[1].get_left_fp(), "lve1");
          eli::octave_finish(1);
        }
      }

#if 0
      // simple 2nd degree curve with 1st derivative specified
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 1st derivative
        p << 1, -1, 1;
        joints[0].set_right_fp(p);

        // set the maximum degrees of each segment
        max_degree[0]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }

      // simple 2nd degree curve with 2nd derivative specified
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 2nd derivative
        p << 1, -1, 1;
        joints[0].set_right_fpp(p);

        // set the maximum degrees of each segment
        max_degree[0]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fpp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fpp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fpp()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }

      // simple 3rd degree curve with both 1st derivatives specified
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 1st derivatives
        p << 1, -1, 1;
        joints[0].set_right_fp(p);
        p << 1, 0, 0;
        joints[1].set_left_fp(p);

        // set the maximum degrees of each segment to be too small
        max_degree[0]=2;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(!rtn_flag);

        // set the maximum degrees of each segment
        max_degree[0]=4;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }

      // simple 5th degree curve with both 1st and 2nd derivatives specified
      {
        index_type nsegs(1);
        std::vector<typename general_creator_type::joint_data> joints(nsegs+1);
        std::vector<typename general_creator_type::index_type> max_degree(nsegs);
        point_type p;
        general_creator_type gc;
        piecewise_curve_type c;
        data_type t0(2), t1(4);
        bool rtn_flag;

        // set the joints
        p << 1, 1, 0;
        joints[0].set_f(p);
        p << 0, 0, -1;
        joints[1].set_f(p);

        // set joint 1st derivatives
        p << 1, -1, 1;
        joints[0].set_right_fp(p);
        p << 1, 0, 0;
        joints[1].set_left_fp(p);

        // set joint 2nd derivatives
        p << 0, -1, 0;
        joints[0].set_right_fpp(p);
        p << 0, 1, 0;
        joints[1].set_left_fpp(p);

        // set the maximum degrees of each segment
        max_degree[0]=6;

        // create curve
        rtn_flag=gc.set_conditions(joints, max_degree, false);
        TEST_ASSERT(rtn_flag);
        gc.set_t0(t0);
        gc.set_segment_dt(t1-t0, 0);
        rtn_flag=gc.create(c);
        TEST_ASSERT(rtn_flag);

        // test the resulting curve
        point_type p_test;

        p_test=c.f(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_f(), p_test));
//        std::cout << "p=" << joints[0].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_f()-p_test).norm() << std::endl;
        p_test=c.f(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_f(), p_test));
//        std::cout << "p=" << joints[1].get_f()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_f()-p_test).norm() << std::endl;
        p_test=c.fp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fp()-p_test).norm() << std::endl;
        p_test=c.fp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fp()-p_test).norm() << std::endl;
        p_test=c.fpp(t0);
        TEST_ASSERT(tol.approximately_equal(joints[0].get_right_fpp(), p_test));
//        std::cout << "p=" << joints[0].get_right_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[0].get_right_fpp()-p_test).norm() << std::endl;
        p_test=c.fpp(t1);
        TEST_ASSERT(tol.approximately_equal(joints[1].get_left_fpp(), p_test));
//        std::cout << "p=" << joints[1].get_left_fpp()
//                  << "\tp_test=" << p_test << "\tdiff="
//                  << (joints[1].get_left_fpp()-p_test).norm() << std::endl;

//        if (rtn_flag && (typeid(data_type)==typeid(double)))
//          octave_print(1, c);
      }
#endif
    }

};

#endif

