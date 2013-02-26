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

#ifndef piecewise_curve_test_suite_hpp
#define piecewise_curve_test_suite_hpp

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/length.hpp"
#include "eli/geom/curve/curvature.hpp"
#include "eli/geom/curve/piecewise.hpp"

#include <cmath>    // std::pow, std::exp
#include <cassert>  // assert()

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

template<typename data__>
class piecewise_curve_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_test_suite<float>::creation_test);
      TEST_ADD(piecewise_curve_test_suite<float>::reverse_test);
      TEST_ADD(piecewise_curve_test_suite<float>::replace_test);
      TEST_ADD(piecewise_curve_test_suite<float>::evaluation_test);
      TEST_ADD(piecewise_curve_test_suite<float>::split_test);
      TEST_ADD(piecewise_curve_test_suite<float>::length_test);
//       TEST_ADD(piecewise_curve_test_suite<float>::fit_free_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<float>::fit_C0_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<float>::fit_C1_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<float>::fit_C2_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<float>::fit_closed_test);
//       TEST_ADD(piecewise_curve_test_suite<float>::interpolate_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_test_suite<double>::creation_test);
      TEST_ADD(piecewise_curve_test_suite<double>::replace_test);
      TEST_ADD(piecewise_curve_test_suite<double>::evaluation_test);
      TEST_ADD(piecewise_curve_test_suite<double>::split_test);
      TEST_ADD(piecewise_curve_test_suite<double>::length_test);
//       TEST_ADD(piecewise_curve_test_suite<double>::fit_free_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<double>::fit_C0_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<double>::fit_C1_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<double>::fit_C2_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<double>::fit_closed_test);
//       TEST_ADD(piecewise_curve_test_suite<double>::interpolate_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_test_suite<long double>::creation_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::replace_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::evaluation_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::split_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::length_test);
//       TEST_ADD(piecewise_curve_test_suite<long double>::fit_free_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<long double>::fit_C0_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<long double>::fit_C1_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<long double>::fit_C2_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<long double>::fit_closed_test);
//       TEST_ADD(piecewise_curve_test_suite<long double>::interpolate_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_test_suite<dd_real>::creation_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::replace_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::evaluation_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::split_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::length_test);
//       TEST_ADD(piecewise_curve_test_suite<dd_real>::fit_free_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<dd_real>::fit_C0_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<dd_real>::fit_C1_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<dd_real>::fit_C2_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<dd_real>::fit_closed_test);
//       TEST_ADD(piecewise_curve_test_suite<dd_real>::interpolate_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_test_suite<qd_real>::creation_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::replace_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::evaluation_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::split_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::length_test);
//       TEST_ADD(piecewise_curve_test_suite<qd_real>::fit_free_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<qd_real>::fit_C0_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<qd_real>::fit_C1_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<qd_real>::fit_C2_ends_test);
//       TEST_ADD(piecewise_curve_test_suite<qd_real>::fit_closed_test);
//       TEST_ADD(piecewise_curve_test_suite<qd_real>::interpolate_test);
    }
#endif
  public:
    piecewise_curve_test_suite()
    {
      AddTests(data__());
    }
    ~piecewise_curve_test_suite()
    {
    }

  private:
    void creation_test()
    {
      piecewise_curve_type c1, c2;
      curve_type bc[3], bc2, bc_out;
      data_type dt[3], dt_out;
      index_type i;
      typename curve_type::control_point_type cntrl1_in[4], cntrl2_in[5], cntrl3_in[3], cntrl2a_in[2];
      typename piecewise_curve_type::error_code err;

      // create bezier curves
      cntrl1_in[0] << 2.0, 2.0, 0.0;
      cntrl1_in[1] << 1.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 0.0, 0.0;
      cntrl1_in[3] << 4.0, 1.0, 0.0;
      dt[0]=0.5;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 4.0, 1.0, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 1.0, 0.0;
      cntrl2_in[3] << 6.0, 0.0, 0.0;
      cntrl2_in[4] << 6.5,-0.5, 0.0;
      dt[1]=2.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 6.5,-0.5, 0.0;
      cntrl3_in[1] << 6.0,-1.0, 0.0;
      cntrl3_in[2] << 5.5,-2.0, 0.0;
      dt[2]=1.5;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }

      // initialize by passing iterators to curve collection
      err=c1.set(bc, bc+2);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==2);
      err=c1.get(bc_out, dt_out, 0);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc[0]==bc_out);
      TEST_ASSERT(dt_out==1);
      err=c1.get(bc_out, dt_out, 1);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc[1]==bc_out);
      TEST_ASSERT(dt_out==1);

      // initialize by setting one at a time
      c1.clear();
      err=c1.push_back(bc[0], dt[0]);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==1);
      err=c1.push_back(bc[1], dt[1]);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==2);
      err=c1.get(bc_out, dt_out, 0);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc[0]==bc_out);
      TEST_ASSERT(dt[0]==dt_out);
      err=c1.get(bc_out, dt_out, 1);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc[1]==bc_out);
      TEST_ASSERT(dt[1]==dt_out);

      // try to get segment that doesn't exist
      err=c1.get(bc_out, dt_out, 3);
      TEST_ASSERT(err==piecewise_curve_type::INVALID_INDEX);

      // try to set segment that doesn't connect to neighbor(s)
      err=c1.push_back(bc[0]);
      TEST_ASSERT(err==piecewise_curve_type::SEGMENT_NOT_CONNECTED);

      // append segment to start
      cntrl2a_in[0] << 1.0, 1.0, 0.0;
      cntrl2a_in[1] << 2.0, 2.0, 0.0;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.push_front(bc2);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      err=c1.get(bc_out, dt_out, 0);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc2==bc_out);
      TEST_ASSERT(dt_out==1);

      // test changing t0
      c1.set_t0(-1.5);
      TEST_ASSERT(c1.get_t0()==-1.5);

      // test copy ctr
      piecewise_curve_type c1copy(c1);
      TEST_ASSERT(c1==c1copy);

      // test assignment operator
      c2=c1;
      TEST_ASSERT(c2==c1);
    }

    void replace_test()
    {
      piecewise_curve_type c1, c1c, c2, c3;
      curve_type bc[3], bc2, bc_out;
      data_type dt[3], dt2, dt_out;
      index_type i;
      typename curve_type::control_point_type cntrl1_in[4], cntrl2_in[5], cntrl3_in[3], cntrl2a_in[2];
      typename piecewise_curve_type::error_code err;

      // create bezier curves
      cntrl1_in[0] << 2.0, 2.0, 0.0;
      cntrl1_in[1] << 1.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 0.0, 0.0;
      cntrl1_in[3] << 4.0, 1.0, 0.0;
      dt[0]=0.5;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 4.0, 1.0, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 1.0, 0.0;
      cntrl2_in[3] << 6.0, 0.0, 0.0;
      cntrl2_in[4] << 6.5,-0.5, 0.0;
      dt[1]=2.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 6.5,-0.5, 0.0;
      cntrl3_in[1] << 6.0,-1.0, 0.0;
      cntrl3_in[2] << 5.5,-2.0, 0.0;
      dt[2]=1.5;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }

      // initialize by passing iterators to curve collection
      err=c1.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);

      // replace segment correctly
      cntrl2a_in[0] << 4.0, 1.0, 0.0;
      cntrl2a_in[1] << 6.5,-0.5, 0.0;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, 1);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      err=c1.get(bc_out, dt_out, 1);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc2==bc_out);
      TEST_ASSERT(dt[1]==dt_out);

      // replace segment with segment that doesn't connect to neighbor(s)
      cntrl2a_in[0] << 4.0, 1.0, 0.0;
      cntrl2a_in[1] << 5.5,-0.5, 0.0;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, 1);
      TEST_ASSERT(err==piecewise_curve_type::SEGMENT_NOT_CONNECTED);
      cntrl2a_in[0] << 4.0, 2.0, 0.0;
      cntrl2a_in[1] << 6.5,-0.5, 0.0;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, 1);
      TEST_ASSERT(err==piecewise_curve_type::SEGMENT_NOT_CONNECTED);

      // replace one segment with several segments
      cntrl1_in[0] << 4.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=0.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 2.0, 0.0;
      cntrl2_in[3] << 5.0, 2.0, 0.0;
      cntrl2_in[4] << 5.0, 0.5, 0.0;
      dt[1]=1.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 5.0, 0.5, 0.0;
      cntrl3_in[1] << 5.5, 0.0, 0.0;
      cntrl3_in[2] << 6.5,-0.5, 0.0;
      dt[2]=3.0;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      c1c=c1;
      err=c1.replace(c2, 1);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==5);
      err=c1.get(bc_out, dt_out, 1);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc[0]==bc_out);
      TEST_ASSERT(dt[0]==dt_out);
      err=c1.get(bc_out, dt_out, 2);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc[1]==bc_out);
      TEST_ASSERT(dt[1]==dt_out);
      err=c1.get(bc_out, dt_out, 3);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc[2]==bc_out);
      TEST_ASSERT(dt[2]==dt_out);
      c3=c1;

      // replace segment with piecewise that doesn't connect to neighbor(s)
      c1=c1c;
      cntrl1_in[0] << 5.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=0.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 2.0, 0.0;
      cntrl2_in[3] << 5.0, 2.0, 0.0;
      cntrl2_in[4] << 5.0, 0.5, 0.0;
      dt[1]=1.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 5.0, 0.5, 0.0;
      cntrl3_in[1] << 5.5, 0.0, 0.0;
      cntrl3_in[2] << 6.5,-0.5, 0.0;
      dt[2]=3.0;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      err=c1.replace(c2, 1);
      TEST_ASSERT(err==piecewise_curve_type::SEGMENT_NOT_CONNECTED);
      cntrl1_in[0] << 4.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=0.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 2.0, 0.0;
      cntrl2_in[3] << 5.0, 2.0, 0.0;
      cntrl2_in[4] << 5.0, 0.5, 0.0;
      dt[1]=1.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 5.0, 0.5, 0.0;
      cntrl3_in[1] << 5.5, 0.0, 0.0;
      cntrl3_in[2] << 7.5,-0.5, 0.0;
      dt[2]=3.0;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      err=c1.replace(c2, 1);
      TEST_ASSERT(err==piecewise_curve_type::SEGMENT_NOT_CONNECTED);

      // replace several segments with one segment
      c1=c3;
      c1c=c1;
      cntrl2a_in[0] << 4.0, 1.0, 0.0;
      cntrl2a_in[1] << 5.0, 0.5, 0.0;
      dt2=1.25;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, dt2, 1, 3);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==4);
      err=c1.get(bc_out, dt_out, 1);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc2==bc_out);
      TEST_ASSERT(dt2==dt_out);
      c3=c1;

      // replace segments with segment that doesn't connect to neighbor(s)
      c1=c1c;
      cntrl2a_in[0] << 4.5, 1.0, 0.0;
      cntrl2a_in[1] << 5.0, 0.5, 0.0;
      dt2=1.25;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, dt2, 1, 3);
      TEST_ASSERT(err==piecewise_curve_type::SEGMENT_NOT_CONNECTED);
      cntrl2a_in[0] << 4.0, 1.0, 0.0;
      cntrl2a_in[1] << 5.5, 0.5, 0.0;
      dt2=1.25;
      bc2.resize(1);
      for (i=0; i<2; ++i)
      {
        bc2.set_control_point(cntrl2a_in[i], i);
      }
      err=c1.replace(bc2, dt2, 1, 3);
      TEST_ASSERT(err==piecewise_curve_type::SEGMENT_NOT_CONNECTED);

      // replace several segments with several segments
      c1=c3;
      c1c=c1;
      cntrl1_in[0] << 4.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=0.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 4.5, 2.0, 0.0;
      cntrl2_in[3] << 4.0, 2.0, 0.0;
      cntrl2_in[4] << 3.5, 0.5, 0.0;
      dt[1]=1.0;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 3.5, 0.5, 0.0;
      cntrl3_in[1] << 4.5, 0.0, 0.0;
      cntrl3_in[2] << 6.5,-0.5, 0.0;
      dt[2]=3.0;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      err=c1.replace(c2, 1, 3);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==5);
      err=c1.get(bc_out, dt_out, 1);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc[0]==bc_out);
      TEST_ASSERT(dt[0]==dt_out);
      err=c1.get(bc_out, dt_out, 2);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc[1]==bc_out);
      TEST_ASSERT(dt[1]==dt_out);
      err=c1.get(bc_out, dt_out, 3);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(bc[2]==bc_out);
      TEST_ASSERT(dt[2]==dt_out);

      // replace segments with piecewise that doesn't connect to neighbor(s)
      c1=c1c;
      cntrl1_in[0] << 3.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=1.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 4.5, 2.0, 0.0;
      cntrl2_in[3] << 4.0, 2.0, 0.0;
      cntrl2_in[4] << 3.5, 0.5, 0.0;
      dt[1]=1.5;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 3.5, 0.5, 0.0;
      cntrl3_in[1] << 4.5, 0.0, 0.0;
      cntrl3_in[2] << 6.5,-0.5, 0.0;
      dt[2]=3.5;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      err=c1.replace(c2, 1, 3);
      TEST_ASSERT(err==piecewise_curve_type::SEGMENT_NOT_CONNECTED);
      cntrl1_in[0] << 4.0, 1.0, 0.0;
      cntrl1_in[1] << 3.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 2.0, 0.0;
      cntrl1_in[3] << 5.5, 2.5, 0.0;
      dt[0]=1.25;
      bc[0].resize(3);
      for (i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      cntrl2_in[0] << 5.5, 2.5, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 4.5, 2.0, 0.0;
      cntrl2_in[3] << 4.0, 2.0, 0.0;
      cntrl2_in[4] << 3.5, 0.5, 0.0;
      dt[1]=1.5;
      bc[1].resize(4);
      for (i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      cntrl3_in[0] << 3.5, 0.5, 0.0;
      cntrl3_in[1] << 4.5, 0.0, 0.0;
      cntrl3_in[2] << 5.5,-0.5, 0.0;
      dt[2]=3.5;
      bc[2].resize(2);
      for (i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      err=c2.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      err=c1.replace(c2, 1, 3);
      TEST_ASSERT(err==piecewise_curve_type::SEGMENT_NOT_CONNECTED);
    }

    void evaluation_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // test two curves with delta t=1
      {
        piecewise_curve_type pwc;
        typename curve_type::control_point_type cntrl_in[4];
        curve_type bc1, bc1l, bc1r;
        point_type eval_out, eval_ref;
        data_type t, tl, tr, ts;
        tl = static_cast<data__>(0.3);
        tr = static_cast<data__>(0.87);
        ts = static_cast<data__>(0.586);

        // set control points
        cntrl_in[0] << 0, 0, 0;
        cntrl_in[1] << 0, 2, 0;
        cntrl_in[2] << 8, 2, 0;
        cntrl_in[3] << 4, 0, 0;
        bc1.resize(3);
        for (index_type i=0; i<4; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        // split curve and create piecewise
        bc1.split(bc1l, bc1r, ts);
        pwc.push_back(bc1l);
        pwc.push_back(bc1r);
        TEST_ASSERT(pwc.number_segments()==2);

        // check the left curve
        t=tl*ts;
        eval_out=pwc.f(tl);
        eval_ref=bc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fp(tl);
        eval_ref=bc1.fp(t)*ts;
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fpp(tl);
        eval_ref=bc1.fpp(t)*ts*ts;
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fppp(tl);
        eval_ref=bc1.fppp(t)*ts*ts*ts;
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);

        // check the right curve
        t=ts+tr*(1-ts);
        eval_out=pwc.f(1+tr);
        eval_ref=bc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fp(1+tr);
        eval_ref=bc1.fp(t)*(1-ts);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fpp(1+tr);
        eval_ref=bc1.fpp(t)*(1-ts)*(1-ts);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fppp(1+tr);
        eval_ref=bc1.fppp(t)*(1-ts)*(1-ts)*(1-ts);
        TEST_ASSERT((eval_out-eval_ref).norm()<1.21e2*eps);
      }

      // test two curves with delta t!=1
      {
        piecewise_curve_type pwc;
        typename curve_type::control_point_type cntrl_in[4];
        curve_type bc1, bc1l, bc1r;
        point_type eval_out, eval_ref;
        data_type tl, tr, ts;
        tl = static_cast<data__>(0.3);
        tr = static_cast<data__>(0.87);
        ts = static_cast<data__>(0.586);

        // set control points
        cntrl_in[0] << 0, 0, 0;
        cntrl_in[1] << 0, 2, 0;
        cntrl_in[2] << 8, 2, 0;
        cntrl_in[3] << 4, 0, 0;
        bc1.resize(3);
        for (index_type i=0; i<4; ++i)
        {
          bc1.set_control_point(cntrl_in[i], i);
        }

        // split curve and create piecewise
        bc1.split(bc1l, bc1r, ts);
        pwc.push_back(bc1l, ts);
        pwc.push_back(bc1r, 1-ts);
        TEST_ASSERT(pwc.number_segments()==2);

        // check the left curve
        eval_out=pwc.f(tl);
        eval_ref=bc1.f(tl);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fp(tl);
        eval_ref=bc1.fp(tl);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fpp(tl);
        eval_ref=bc1.fpp(tl);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fppp(tl);
        eval_ref=bc1.fppp(tl);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);

        // check the right curve
        eval_out=pwc.f(tr);
        eval_ref=bc1.f(tr);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fp(tr);
        eval_ref=bc1.fp(tr);
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
        eval_out=pwc.fpp(tr);
        eval_ref=bc1.fpp(tr);
        TEST_ASSERT((eval_out-eval_ref).norm()<2.6e2*eps);
        eval_out=pwc.fppp(tr);
        eval_ref=bc1.fppp(tr);
        TEST_ASSERT((eval_out-eval_ref).norm()<1.8e3*eps);
      }
    }

    void split_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
        piecewise_curve_type pwc0, pwc1;
        typename curve_type::control_point_type cntrl_in[4];
        typename piecewise_curve_type::error_code err;
        curve_type bc;
        point_type eval_out, eval_ref;
        data_type ts, t;
        ts=static_cast<data__>(1.56);

        // build piecewise curve
        cntrl_in[0] << 0, 0, 0;
        cntrl_in[1] << 0, 2, 0;
        cntrl_in[2] << 8, 2, 0;
        cntrl_in[3] << 4, 0, 0;
        bc.resize(3);
        for (index_type i=0; i<4; ++i)
        {
          bc.set_control_point(cntrl_in[i], i);
        }
        err=pwc0.push_back(bc);
        TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
        cntrl_in[0] << 4,  0,   0;
        cntrl_in[1] << 3, -0.5, 0;
        cntrl_in[2] << 2, -1,   0;
        cntrl_in[3] << 1, -1,   0;
        bc.resize(3);
        for (index_type i=0; i<4; ++i)
        {
          bc.set_control_point(cntrl_in[i], i);
        }
        err=pwc0.push_back(bc);
        TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
        TEST_ASSERT(pwc0.number_segments()==2);

        // split curve and create piecewise
        pwc1=pwc0;
        err=pwc0.split(ts);
        TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
        TEST_ASSERT(pwc0.number_segments()==pwc1.number_segments()+1);

        t=0;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT(eval_out==eval_ref);

        t=0.5;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT(eval_out==eval_ref);

        t=1;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT(eval_out==eval_ref);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT(eval_out==eval_ref);

        t=1.25;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<18*eps);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<138*eps);

        t=1.5;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<7*eps);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<35*eps);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<138*eps);

        t=1.75;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<2*eps);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<3*eps);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<28*eps);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<142*eps);

        t=2;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<1*eps);
        eval_out=pwc0.fp(t);
        eval_ref=pwc1.fp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<9*eps);
        eval_out=pwc0.fpp(t);
        eval_ref=pwc1.fpp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<63*eps);
        eval_out=pwc0.fppp(t);
        eval_ref=pwc1.fppp(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<142*eps);
    }

    void length_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
      piecewise_curve_type c1;
      curve_type bc[3];
      data_type dt[3], len, bc_len[3], ref_len, t0, t1, temp0, temp1;
      typename curve_type::control_point_type cntrl1_in[4], cntrl2_in[5], cntrl3_in[3];
      typename piecewise_curve_type::error_code err;
      data_type tol(std::sqrt(eps));

      // create piecewise curve
      cntrl1_in[0] << 2.0, 2.0, 0.0;
      cntrl1_in[1] << 1.0, 1.5, 0.0;
      cntrl1_in[2] << 3.5, 0.0, 0.0;
      cntrl1_in[3] << 4.0, 1.0, 0.0;
      dt[0]=0.5;
      bc[0].resize(3);
      for (index_type i=0; i<4; ++i)
      {
        bc[0].set_control_point(cntrl1_in[i], i);
      }
      eli::geom::curve::length(bc_len[0], bc[0], tol);
      cntrl2_in[0] << 4.0, 1.0, 0.0;
      cntrl2_in[1] << 5.0, 2.5, 0.0;
      cntrl2_in[2] << 5.5, 1.0, 0.0;
      cntrl2_in[3] << 6.0, 0.0, 0.0;
      cntrl2_in[4] << 6.5,-0.5, 0.0;
      dt[1]=2.0;
      bc[1].resize(4);
      for (index_type i=0; i<5; ++i)
      {
        bc[1].set_control_point(cntrl2_in[i], i);
      }
      eli::geom::curve::length(bc_len[1], bc[1], tol);
      cntrl3_in[0] << 6.5,-0.5, 0.0;
      cntrl3_in[1] << 6.0,-1.0, 0.0;
      cntrl3_in[2] << 5.5,-2.0, 0.0;
      dt[2]=1.5;
      bc[2].resize(2);
      for (index_type i=0; i<3; ++i)
      {
        bc[2].set_control_point(cntrl3_in[i], i);
      }
      eli::geom::curve::length(bc_len[2], bc[2], tol);
      err=c1.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(c1.number_segments()==3);

      // create two segment curve calc length of each segment to compare
      eli::geom::curve::length(len, c1, tol);
      ref_len=bc_len[0]+bc_len[1]+bc_len[2];
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
      {
        TEST_ASSERT(std::abs(len-ref_len)<3*eps);
      }
      else
#endif
      {
        TEST_ASSERT(len==ref_len);
      }

      // choose part of first segment to calc length and compare
      t0=0.125;
      t1=0.375;
      eli::geom::curve::length(len, c1, t0, t1, tol);
      eli::geom::curve::length(ref_len, bc[0], 0.25, 0.75, tol);
      TEST_ASSERT(len==ref_len);

      // choose part of second segment to calc length and compare
      t0=0.25;
      t1=1.5;
      eli::geom::curve::length(len, c1, t0, t1, tol);
      eli::geom::curve::length(temp0, bc[0], 0.5, 1, tol);
      eli::geom::curve::length(temp1, bc[1], 0, 0.5, tol);
      ref_len=temp0+temp1;
      TEST_ASSERT(std::abs(len-ref_len)<2*tol);

      // choose part of third segment to calc length and compare
      t0=0.25;
      t1=3.25;
      eli::geom::curve::length(len, c1, t0, t1, tol);
      eli::geom::curve::length(temp0, bc[0], 0.5, 1, tol);
      eli::geom::curve::length(temp1, bc[2], 0, 0.5, tol);
      ref_len=temp0+bc_len[1]+temp1;
      TEST_ASSERT(std::abs(len-ref_len)<166*tol);
    }

#if 0
    void fit_free_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
      // open curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::NOT_CONNECTED);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.110);

        // check if went through points
        TEST_ASSERT(pts[0]!=bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT(bez.f(0)!=bez.f(1));

//         octave_print(1, pts, bez);
      }

      // closed curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.765);

        // check if went through points
        TEST_ASSERT(pts[0]!=bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT(bez.f(0)!=bez.f(1));

//         octave_print(2, pts, bez);
      }
    }

    void fit_C0_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
      // open curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::NOT_CONNECTED);
        fcon.add_start_C0_constraint();
        fcon.add_end_C0_constraint();

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          eli::geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.105);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT((pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()<35*eps);
        TEST_ASSERT(bez.f(0)!=bez.f(1));

//         octave_print(3, pts, bez);
      }

      // closed curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);
        fcon.add_start_C0_constraint();

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.207);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<172*eps);

//         octave_print(4, pts, bez);
      }

      // closed curve with additional point constraint
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);
        fcon.add_start_C0_constraint();
        fcon.add_C0_constraint(5);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.233);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<197*eps);
        TEST_ASSERT((bez.f(t[5])-pts[5]).norm()<14*eps);

//         octave_print(5, pts, bez);
      }
    }

    void fit_C1_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // open curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1, fp2;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fp2(0)=-eli::constants::math<data_type>::two_pi()*pts[pts.size()-1](1);
        fp2(1)=eli::constants::math<data_type>::two_pi()*pts[pts.size()-1](0);
        fp2(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::NOT_CONNECTED);
        fcon.add_start_C1_constraint(fp1);
        fcon.add_end_C1_constraint(fp2);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.366);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT(fp1==bez.fp(t[0]));
        TEST_ASSERT((pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()<69*eps);
        TEST_ASSERT((fp2-bez.fp(t[t.size()-1])).norm()<289*eps);
        TEST_ASSERT(bez.f(0)!=bez.f(1));

//         octave_print(6, pts, bez);
      }

      // closed curve
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);
        fcon.add_start_C1_constraint(fp1);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.331);

        // check if went through points
        TEST_ASSERT(pts[0]==bez.f(t[0]));
        TEST_ASSERT(fp1==bez.fp(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<66*eps);

//         octave_print(7, pts, bez);
      }

      // closed curve with additional point constraint
      {
        fit_container_type fcon;
        size_t i, deg(5);
        std::vector<point_type> pts(10);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1, fp2;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fp2(0)=-eli::constants::math<data_type>::two_pi()*pts[5](1);
        fp2(1)=eli::constants::math<data_type>::two_pi()*pts[5](0);
        fp2(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);
        fcon.add_start_C1_constraint(fp1);
        fcon.add_C1_constraint(5, fp2);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.418);

        // check if went through points
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<5*eps);
        TEST_ASSERT(fp1==bez.fp(t[0]));
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<97*eps);
        TEST_ASSERT((bez.f(t[5])-pts[5]).norm()<6*eps);
        TEST_ASSERT((bez.fp(t[5])-fp2).norm()<73*eps);

//         octave_print(8, pts, bez);
      }
    }

    void fit_C2_ends_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // open curve
      {
        fit_container_type fcon;
        size_t i, deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1, fpp1, fp2, fpp2;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fpp1(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](0);
        fpp1(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](1);
        fpp1(2)=0;
        fp2(0)=-eli::constants::math<data_type>::two_pi()*pts[pts.size()-1](1);
        fp2(1)=eli::constants::math<data_type>::two_pi()*pts[pts.size()-1](0);
        fp2(2)=0;
        fpp2(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[pts.size()-1](0);
        fpp2(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[pts.size()-1](1);
        fpp2(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::NOT_CONNECTED);
        fcon.add_start_C2_constraint(fp1, fpp1);
        fcon.add_end_C2_constraint(fp2, fpp2);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.open());
        TEST_ASSERT(err < 0.372);

        // check if went through points
        TEST_ASSERT((pts[pts.size()-1]-bez.f(t[t.size()-1])).norm()<192*eps);
        TEST_ASSERT((fp2-bez.fp(t[t.size()-1])).norm()<956*eps);
        TEST_ASSERT((fpp2-bez.fpp(t[t.size()-1])).norm()<3.94e3*eps);
        TEST_ASSERT(bez.f(0)!=bez.f(1));
        TEST_ASSERT((pts[0]-bez.f(t[0])).norm()<4*eps);
        TEST_ASSERT((fp1-bez.fp(t[0])).norm()<9*eps);
        TEST_ASSERT(fpp1==bez.fpp(t[0]));

//         octave_print(9, pts, bez);
      }

      // closed curve
      {
        fit_container_type fcon;
        size_t i, deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1, fpp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fpp1(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](0);
        fpp1(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](1);
        fpp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);
        fcon.add_start_C2_constraint(fp1, fpp1);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0479);

        // check if went through points
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<387*eps);
        TEST_ASSERT((pts[0]-bez.f(t[0])).norm()<3*eps);
        TEST_ASSERT((fp1-bez.fp(t[0])).norm()<17*eps);
        TEST_ASSERT((fpp1-bez.fpp(t[0])).norm()<65*eps);
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));

//         octave_print(10, pts, bez);
      }

      // closed curve with additional point constraint
      {
        fit_container_type fcon;
        size_t i, deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1, fpp1, fp2, fpp2;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp1(2)=0;
        fpp1(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](0);
        fpp1(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](1);
        fpp1(2)=0;
        fp2(0)=-eli::constants::math<data_type>::two_pi()*pts[10](1);
        fp2(1)=eli::constants::math<data_type>::two_pi()*pts[10](0);
        fp2(2)=0;
        fpp2(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[10](0);
        fpp2(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[10](1);
        fpp2(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);
        fcon.add_start_C2_constraint(fp1, fpp1);
        fcon.add_C2_constraint(10, fp2, fpp2);

        // fit points
        bez.fit(t, fcon, deg);

        // calculate the error at the point
        data_type err(0);
        for (i=0; i<pts.size(); ++i)
        {
          data_type d;
          geom::point::distance(d, pts[i], bez.f(t[i]));
          err+=d;
        }

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0785);

        // check if went through points
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<420*eps);
        TEST_ASSERT((bez.f(t[10])-pts[10]).norm()<39*eps);
        TEST_ASSERT((bez.fp(t[10])-fp2).norm()<142*eps);
        TEST_ASSERT((bez.fpp(t[10])-fpp2).norm()<1.26e3*eps);
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<5*eps);
        TEST_ASSERT((bez.fp(t[0])-fp1).norm()<5*eps);
        TEST_ASSERT((bez.fpp(t[0])-fpp1).norm()<33*eps);
        TEST_ASSERT(pts[pts.size()-1]!=bez.f(t[t.size()-1]));

//         octave_print(11, pts, bez);
      }
    }

    void fit_closed_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // fit no constraints C0 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.393);

        // check if went through points
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<245*eps);
        TEST_ASSERT(bez.fp(0)!=bez.fp(1));
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

//         octave_print(12, pts, bez);
      }

      // fit no constraints C1 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C1);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0312);

        // check if went through points
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<492*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<2.56e3*eps);
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

//         octave_print(13, pts, bez);
      }

      // fit no constraints C2 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C2);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0342);

        // check if went through points
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<367*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<1.41e3*eps);
        TEST_ASSERT((bez.fpp(1)-bez.fpp(0)).norm()<7.97e3*eps);

//         octave_print(14, pts, bez);
      }

      // fit start constraints C0 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);
        fcon.add_start_C0_constraint();

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0185);

        // check if went through points
        TEST_ASSERT(bez.f(0)==pts[0]);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<431*eps);
        TEST_ASSERT(bez.fp(0)!=bez.fp(1));
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

//         octave_print(15, pts, bez);
      }

      // fit start constraints C1 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp0;
        fp0(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp0(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp0(2)=eli::constants::math<data_type>::two_pi()*pts[0](2);
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C1);
        fcon.add_start_C1_constraint(fp0);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0395);

        // check if went through points
        TEST_ASSERT((bez.f(0)-pts[0]).norm()<2*eps);
        TEST_ASSERT(bez.fp(0)==fp0);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<483*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<2.77e3*eps);
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

//         octave_print(16, pts, bez);
      }

      // fit start constraints C2 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp0, fpp0;
        fp0(0)=-eli::constants::math<data_type>::two_pi()*pts[0](1);
        fp0(1)=eli::constants::math<data_type>::two_pi()*pts[0](0);
        fp0(2)=eli::constants::math<data_type>::two_pi()*pts[0](2);
        fpp0(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](0);
        fpp0(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[0](1);
        fpp0(2)=4*eli::constants::math<data_type>::pi_squared()*pts[0](2);
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C2);
        fcon.add_start_C2_constraint(fp0, fpp0);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.141);

        // check if went through points
        TEST_ASSERT((bez.f(0)-pts[0]).norm()<2*eps);
        TEST_ASSERT((bez.fp(0)-fp0).norm()<9*eps);
        TEST_ASSERT(bez.fpp(0)==fpp0);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<409*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<2.36e3*eps);
        TEST_ASSERT((bez.fpp(1)-bez.fpp(0)).norm()<6.80e3*eps);

//         octave_print(17, pts, bez);
      }

      // fit middle constraints C0 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        int ci=5;
        data_type t5;
        point_type f5;
        f5=pts[ci];
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);
        fcon.add_C0_constraint(ci);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0180);

        // check if went through points
        t5=t[ci];
        TEST_ASSERT((bez.f(t5)-f5).norm()<6*eps);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<594*eps);
        TEST_ASSERT(bez.fp(0)!=bez.fp(1));
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

//         octave_print(18, pts, bez);
      }

      // fit middle constraints C1 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        int ci=5;
        data_type t5;
        point_type f5;
        f5=pts[ci];
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C1);
        fcon.add_C0_constraint(ci);

        // fit points
        data_type err;
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0316);

        // check if went through points
        t5=t[ci];
        TEST_ASSERT((bez.f(t5)-f5).norm()<4*eps);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<357*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<1.69e3*eps);
        TEST_ASSERT(bez.fpp(0)!=bez.fpp(1));

//         octave_print(19, pts, bez);
      }

      // fit middle constraints C2 closed
      {
        fit_container_type fcon;
        size_t deg(7);
        std::vector<point_type> pts(20);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        data_type err;
        create_circle(pts);

        // configure fit container
        int ci=5;
        data_type t5;
        point_type f5;
        f5=pts[ci];
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C2);
        fcon.add_C0_constraint(ci);

        // fit points
        err=bez.fit(t, fcon, deg);

        TEST_ASSERT(bez.closed());
        TEST_ASSERT(err < 0.0353);

        // check if went through points
        t5=t[ci];
        TEST_ASSERT((bez.f(t5)-f5).norm()<3*eps);
        TEST_ASSERT((bez.f(1)-bez.f(0)).norm()<322*eps);
        TEST_ASSERT((bez.fp(1)-bez.fp(0)).norm()<1.40e3*eps);
        TEST_ASSERT((bez.fpp(1)-bez.fpp(0)).norm()<9.10e3*eps);

//         octave_print(20, pts, bez);
      }
    }

    void interpolate_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_QD_FOUND
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      // interpolate through all points open
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size());
        TEST_ASSERT(bez.open());

        // check if went through points
        TEST_ASSERT(bez.f(t[0])==pts[0]);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<3*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<7*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<13*eps);

//         octave_print(18, pts, bez);
      }

      // interpolate through all points closed
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        fcon.set_points(pts.begin(), pts.end());
        fcon.set_end_flag(fit_container_type::end_continuity_type::C0);

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size()+1);
        TEST_ASSERT(bez.closed());

        // check if went through points
        TEST_ASSERT(bez.f(t[0])==pts[0]);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<3*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<9*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<37*eps);
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<49*eps);

//         octave_print(19, pts, bez);
      }

      // interpolate through all points and specify C1 open
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[1](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[1](0);
        fp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.add_C1_constraint(1, fp1);

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size()+1);
        TEST_ASSERT(bez.open());

        // check if went through points
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<2*eps);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<5*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<21*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<63*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<17*eps);

//         octave_print(20, pts, bez);
      }

      // interpolate through all points and specify C1 closed
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[1](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[1](0);
        fp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.add_C1_constraint(1, fp1);
        fcon.set_end_flag(fit_container_type::end_continuity_type::C1);

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size()+3);
        TEST_ASSERT(bez.closed());

        // check if went through points
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<5*eps);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<3*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<12*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<50*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<22*eps);
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<202*eps);
        TEST_ASSERT((bez.fp(0)-bez.fp(1)).norm()<761*eps);

//         octave_print(21, pts, bez);
      }

      // interpolate through all points and specify C2 open
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1, fpp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[1](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[1](0);
        fp1(2)=0;
        fpp1(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[1](0);
        fpp1(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[1](1);
        fpp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.add_C2_constraint(1, fp1, fpp1);

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size()+2);
        TEST_ASSERT(bez.open());

        // check if went through points
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<9*eps);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<9*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<84*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<162*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<68*eps);
        TEST_ASSERT((bez.fpp(t[1])-fpp1).norm()<255*eps);

//         octave_print(22, pts, bez);
      }

      // interpolate through all points and specify C2 closed
      {
        fit_container_type fcon;
        std::vector<point_type> pts(4);
        bezier_type bez;
        std::vector<data_type> t;

        // get points
        create_circle(pts);

        // configure fit container
        point_type fp1, fpp1;
        fp1(0)=-eli::constants::math<data_type>::two_pi()*pts[1](1);
        fp1(1)=eli::constants::math<data_type>::two_pi()*pts[1](0);
        fp1(2)=0;
        fpp1(0)=-4*eli::constants::math<data_type>::pi_squared()*pts[1](0);
        fpp1(1)=-4*eli::constants::math<data_type>::pi_squared()*pts[1](1);
        fpp1(2)=0;
        fcon.set_points(pts.begin(), pts.end());
        fcon.add_C2_constraint(1, fp1, fpp1);
        fcon.set_end_flag(fit_container_type::end_continuity_type::C2);

        // fit points
        bez.interpolate(t, fcon);

        TEST_ASSERT(static_cast<size_t>(bez.degree()+1)==pts.size()+5);
        TEST_ASSERT(bez.closed());

        // check if went through points
        TEST_ASSERT((bez.f(t[0])-pts[0]).norm()<4*eps);
        TEST_ASSERT((bez.f(t[1])-pts[1]).norm()<4*eps);
        TEST_ASSERT((bez.f(t[2])-pts[2]).norm()<20*eps);
        TEST_ASSERT((bez.f(t[3])-pts[3]).norm()<100*eps);
        TEST_ASSERT((bez.fp(t[1])-fp1).norm()<39*eps);
        TEST_ASSERT((bez.fpp(t[1])-fpp1).norm()<269*eps);
        TEST_ASSERT((bez.f(0)-bez.f(1)).norm()<489*eps);
        TEST_ASSERT((bez.fp(0)-bez.fp(1)).norm()<3160*eps);
        TEST_ASSERT((bez.fpp(0)-bez.fpp(1)).norm()<14613*eps);

//         octave_print(23, pts, bez);
      }
    }
#endif
};

#endif

