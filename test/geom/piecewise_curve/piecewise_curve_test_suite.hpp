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
      TEST_ADD(piecewise_curve_test_suite<float>::bounding_box_test);
      TEST_ADD(piecewise_curve_test_suite<float>::reverse_test);
      TEST_ADD(piecewise_curve_test_suite<float>::replace_test);
      TEST_ADD(piecewise_curve_test_suite<float>::transformation_test);
      TEST_ADD(piecewise_curve_test_suite<float>::evaluation_test);
      TEST_ADD(piecewise_curve_test_suite<float>::split_test);
      TEST_ADD(piecewise_curve_test_suite<float>::to_cubic_test);
      TEST_ADD(piecewise_curve_test_suite<float>::length_test);
      TEST_ADD(piecewise_curve_test_suite<float>::round_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_test_suite<double>::creation_test);
      TEST_ADD(piecewise_curve_test_suite<double>::bounding_box_test);
      TEST_ADD(piecewise_curve_test_suite<double>::reverse_test);
      TEST_ADD(piecewise_curve_test_suite<double>::replace_test);
      TEST_ADD(piecewise_curve_test_suite<double>::transformation_test);
      TEST_ADD(piecewise_curve_test_suite<double>::evaluation_test);
      TEST_ADD(piecewise_curve_test_suite<double>::split_test);
      TEST_ADD(piecewise_curve_test_suite<double>::to_cubic_test);
      TEST_ADD(piecewise_curve_test_suite<double>::length_test);
      TEST_ADD(piecewise_curve_test_suite<double>::round_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_test_suite<long double>::creation_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::bounding_box_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::reverse_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::replace_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::transformation_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::evaluation_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::split_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::to_cubic_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::length_test);
      TEST_ADD(piecewise_curve_test_suite<long double>::round_test);
    }
#ifdef ELI_USING_QD
    void AddTests(const dd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_test_suite<dd_real>::creation_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::bounding_box_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::reverse_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::replace_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::transformation_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::evaluation_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::split_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::to_cubic_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::length_test);
      TEST_ADD(piecewise_curve_test_suite<dd_real>::round_test);
    }

    void AddTests(const qd_real &)
    {
      // add the tests
      TEST_ADD(piecewise_curve_test_suite<qd_real>::creation_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::bounding_box_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::reverse_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::replace_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::transformation_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::evaluation_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::split_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::to_cubic_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::length_test);
      TEST_ADD(piecewise_curve_test_suite<qd_real>::round_test);
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
  void octave_print(int figno, const piecewise_curve_type &pc)
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

    std::cout << "cp_z=[";
    for (pp=0; pp<ns; ++pp)
    {
      curve_type bez;
      pc.get(bez, pp);
      for (i=0; i<=bez.degree(); ++i)
      {
        std::cout << bez.get_control_point(i).z();
        if (i<bez.degree())
          std::cout << ", ";
        else if (pp<ns-1)
          std::cout << "; ";
      }
      std::cout << std::endl;
    }
    std::cout << "];" << std::endl;

    // initialize the t parameters
    std::vector<data_type> t(129);
    for (i=0; i<static_cast<index_type>(t.size()); ++i)
    {
      t[i]=tmin+(tmax-tmin)*static_cast<data_type>(i)/(t.size()-1);
    }

    // set the surface points
    std::cout << "surf_x=[";
    for (i=0; i<static_cast<index_type>(t.size()); ++i)
    {
      std::cout << pc.f(t[i]).x();
      if (i<static_cast<index_type>(t.size()-1))
        std::cout << ", ";
    }
    std::cout << "];" << std::endl;

    std::cout << "surf_y=[";
    for (i=0; i<static_cast<index_type>(t.size()); ++i)
    {
      std::cout << pc.f(t[i]).y();
      if (i<static_cast<index_type>(t.size()-1))
        std::cout << ", ";
    }
    std::cout << "];" << std::endl;

    std::cout << "surf_z=[";
    for (i=0; i<static_cast<index_type>(t.size()); ++i)
    {
      std::cout << pc.f(t[i]).z();
      if (i<static_cast<index_type>(t.size()-1))
        std::cout << ", ";
    }
    std::cout << "];" << std::endl;

    std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
    std::cout << "plot3(surf_x, surf_y, surf_z, '-k');" << std::endl;
    std::cout << "hold on;" << std::endl;
    std::cout << "plot3(cp_x', cp_y', cp_z', '-ok', 'MarkerFaceColor', [0 0 0]);" << std::endl;
    std::cout << "hold off;" << std::endl;
  }

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

    void bounding_box_test()
    {
      piecewise_curve_type pc;
      curve_type bc[3];
      data_type dt[3];
      index_type i;
      typename curve_type::control_point_type cntrl1_in[4], cntrl2_in[5], cntrl3_in[3];
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
      err=pc.set(bc, bc+3, dt);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);

      // test the bounding box
      typename piecewise_curve_type::bounding_box_type bb;
      point_type pmin_ref, pmax_ref;

      pc.get_bounding_box(bb);
      pmin_ref << 1, -2, 0;
      pmax_ref << 6.5, 2.5, 0;
      TEST_ASSERT(bb.get_min()==pmin_ref);
      TEST_ASSERT(bb.get_max()==pmax_ref);
    }

    void reverse_test()
    {
      piecewise_curve_type c1, c2;
      curve_type bc[3], bc1_out, bc2_out;
      data_type dt[3], dt1_out, dt2_out;
      index_type i;
      typename curve_type::control_point_type cntrl1_in[4], cntrl2_in[5], cntrl3_in[3];
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

      c1=c2;
      c2.reverse();
      for (i=0; i<c1.number_segments(); ++i)
      {
        c1.get(bc1_out, dt1_out, i);
        c2.get(bc2_out, dt2_out, c1.number_segments()-i-1);
        for (index_type ii=0; ii<=bc1_out.degree(); ++ii)
        {
          TEST_ASSERT(bc1_out.get_control_point(ii)==bc2_out.get_control_point(bc1_out.degree()-ii));
        }
      }
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

    void transformation_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_USING_QD
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif

      piecewise_curve_type pwc, pwc2;
      typename curve_type::control_point_type cntrl_in[4];
      curve_type bc1, bc1l, bc1r;
      point_type eval_out, eval_ref;
      data_type tl, ts;
      tl = static_cast<data__>(0.3);
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

      // test translation
      {
        point_type trans;

        // set up translation vector and apply
        pwc2=pwc;
        trans << 2, 1, 3;
        pwc2.translate(trans);

        // check the left curve
        eval_out=pwc2.f(tl);
        eval_ref=pwc.f(tl)+trans;
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      }

      // test rotation about origin
      {
        typename piecewise_curve_type::rotation_matrix_type rmat;

        // set up rotation and apply
        pwc2=pwc;
        rmat << cos(1), 0, -sin(1),
                     0, 1,       0,
                sin(1), 0,  cos(1);
        pwc2.rotate(rmat);

        // check the left curve
        eval_out=pwc2.f(tl);
        eval_ref=pwc.f(tl)*rmat.transpose();
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      }

      // test rotation about point
      {
        typename piecewise_curve_type::rotation_matrix_type rmat;
        point_type rorig;

        // set up rotation and apply
        pwc2=pwc;
        rorig << 2, 1, 3;
        rmat << cos(1), 0, -sin(1),
                     0, 1,       0,
                sin(1), 0,  cos(1);
        pwc2.rotate(rmat, rorig);

        // check the left curve
        eval_out=pwc2.f(tl);
        eval_ref=rorig+(pwc.f(tl)-rorig)*rmat.transpose();
        TEST_ASSERT((eval_out-eval_ref).norm()<1e2*eps);
      }
    }

    void evaluation_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_USING_QD
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
#ifdef ELI_USING_QD
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

    void to_cubic_test()
    {
      // Two-segment piecewise 4th order.
      {
        piecewise_curve_type pwc0, pwc1;
        typename curve_type::control_point_type cntrl_in[5];
        typename piecewise_curve_type::error_code err;
        curve_type bc;
        point_type eval_out, eval_ref;
        data_type t;

        // build piecewise curve
        cntrl_in[0] << 0, 0, 0;
        cntrl_in[1] << 0, 2, 0;
        cntrl_in[2] << 8, 2, 0;
        cntrl_in[3] << 4, 0, 0;
        cntrl_in[4] << 6, 3, 0;
        bc.resize(4);
        for (index_type i=0; i<5; ++i)
          {
            bc.set_control_point(cntrl_in[i], i);
          }
        err=pwc0.push_back(bc);
        TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);

        cntrl_in[0] << 6,  3,   0;
        cntrl_in[1] << 3, -0.5, 0;
        cntrl_in[2] << 2, -1,   0;
        cntrl_in[3] << 1, -1,   0;
        cntrl_in[4] << 2,  1,   0;
        bc.resize(4);
        for (index_type i=0; i<5; ++i)
          {
            bc.set_control_point(cntrl_in[i], i);
          }
        err=pwc0.push_back(bc);
        TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
        TEST_ASSERT(pwc0.number_segments()==2);

        index_type mind, maxd;
        pwc0.degree(mind, maxd);

        TEST_ASSERT(mind==4);
        TEST_ASSERT(maxd==4);

        TEST_ASSERT(pwc0.number_segments()==2);

        pwc1=pwc0;

        data_type ttol = 0.0001;
        pwc1.to_cubic(ttol);

//        std::cout << pwc1.number_segments() << std::endl;
//        octave_print(1, pwc0);
//        octave_print(2, pwc1);

        pwc1.degree(mind, maxd);
        TEST_ASSERT(mind==3);
        TEST_ASSERT(maxd==3);

        t=0;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=0.5;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=1;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=1.25;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=1.5;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=1.75;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=2;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);
      }

      // One segment cubic, promoted.
      {
        piecewise_curve_type pwc0, pwc1;
        typename curve_type::control_point_type cntrl_in[4];
        typename piecewise_curve_type::error_code err;
        curve_type bc;
        point_type eval_out, eval_ref;
        point_type p, offset;
        data_type t;

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
        bc.degree_promote_to(bc.degree()+4);

        err=pwc0.push_back(bc);
        TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);

        TEST_ASSERT(pwc0.number_segments()==1);

        index_type mind, maxd;
        pwc0.degree(mind, maxd);

        TEST_ASSERT(mind==7);
        TEST_ASSERT(maxd==7);

        pwc1=pwc0;

        data_type ttol = 0.0001;
        pwc1.to_cubic(ttol);

//        std::cout << pwc1.number_segments() << std::endl;
//        octave_print(1, pwc0);
//        octave_print(2, pwc1);

        TEST_ASSERT(pwc1.number_segments()==1);  // Original was cubic, then promoted.

        pwc1.degree(mind, maxd);
        TEST_ASSERT(mind==3);
        TEST_ASSERT(maxd==3);

        t=0;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=0.235;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=1;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);
      }

      // One segment cubic, promoted and perturbed.
      {
        piecewise_curve_type pwc0, pwc1, pwc2;
        typename curve_type::control_point_type cntrl_in[4];
        typename piecewise_curve_type::error_code err;
        curve_type bc;
        point_type eval_out, eval_ref;
        point_type p, offset;
        data_type t;

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

        err=pwc2.push_back(bc);  // Keep unperturbed low-order curve for reference.

        bc.degree_promote_to(bc.degree()+4);

        p=bc.get_control_point(5); // Pick control point near middle.
        offset << -1, -1, 0;
        p=p+offset; // Perturb point.
        bc.set_control_point(p, 5);  // Set perturbed value.

        err=pwc0.push_back(bc);
        TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);

        TEST_ASSERT(pwc0.number_segments()==1);

        index_type mind, maxd;
        pwc0.degree(mind, maxd);

        TEST_ASSERT(mind==7);
        TEST_ASSERT(maxd==7);

        pwc1=pwc0;

        data_type ttol = 0.0001;
        pwc1.to_cubic(ttol);

//        std::cout << pwc1.number_segments() << std::endl;
//        octave_print(1, pwc2);
//        octave_print(2, pwc0);
//        octave_print(3, pwc1);

        pwc1.degree(mind, maxd);
        TEST_ASSERT(mind==3);
        TEST_ASSERT(maxd==3);

        t=0;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=0.235;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=1;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);
      }

      // One segment cubic
      {
        piecewise_curve_type pwc0, pwc1;
        typename curve_type::control_point_type cntrl_in[4];
        typename piecewise_curve_type::error_code err;
        curve_type bc;
        point_type eval_out, eval_ref;
        data_type t;

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

        TEST_ASSERT(pwc0.number_segments()==1);

        index_type mind, maxd;
        pwc0.degree(mind, maxd);

        TEST_ASSERT(mind==3);
        TEST_ASSERT(maxd==3);

        pwc1=pwc0;

        data_type ttol = 0.0001;
        pwc1.to_cubic(ttol);

//        std::cout << pwc1.number_segments() << std::endl;
//        octave_print(1, pwc0);
//        octave_print(2, pwc1);

        TEST_ASSERT(pwc1.number_segments()==1);  // Original was cubic.

        pwc1.degree(mind, maxd);
        TEST_ASSERT(mind==3);
        TEST_ASSERT(maxd==3);

        t=0;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=0.235;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=1;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);
      }

      // One segment quadratic
      {
        piecewise_curve_type pwc0, pwc1;
        typename curve_type::control_point_type cntrl_in[3];
        typename piecewise_curve_type::error_code err;
        curve_type bc;
        point_type eval_out, eval_ref;
        data_type t;

        // build piecewise curve
        cntrl_in[0] << 0, 0, 0;
        cntrl_in[1] << 0, 2, 0;
        cntrl_in[2] << 8, 2, 0;
        bc.resize(2);
        for (index_type i=0; i<3; ++i)
          {
            bc.set_control_point(cntrl_in[i], i);
          }
        err=pwc0.push_back(bc);
        TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);

        TEST_ASSERT(pwc0.number_segments()==1);

        index_type mind, maxd;
        pwc0.degree(mind, maxd);

        TEST_ASSERT(mind==2);
        TEST_ASSERT(maxd==2);

        pwc1=pwc0;

        data_type ttol = 0.0001;
        pwc1.to_cubic(ttol);

//        std::cout << pwc1.number_segments() << std::endl;
//        octave_print(1, pwc0);
//        octave_print(2, pwc1);

        TEST_ASSERT(pwc1.number_segments()==1);  // Original was quadratic.

        pwc1.degree(mind, maxd);
        TEST_ASSERT(mind==3);
        TEST_ASSERT(maxd==3);

        t=0;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=0.235;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=1;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);
      }

      // One segment linear
      {
        piecewise_curve_type pwc0, pwc1;
        typename curve_type::control_point_type cntrl_in[2];
        typename piecewise_curve_type::error_code err;
        curve_type bc;
        point_type eval_out, eval_ref;
        data_type t;

        // build piecewise curve
        cntrl_in[0] << 0, 0, 0;
        cntrl_in[1] << 1.4, 2, 0;
        bc.resize(1);
        for (index_type i=0; i<2; ++i)
          {
            bc.set_control_point(cntrl_in[i], i);
          }
        err=pwc0.push_back(bc);
        TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);

        TEST_ASSERT(pwc0.number_segments()==1);

        index_type mind, maxd;
        pwc0.degree(mind, maxd);

        TEST_ASSERT(mind==1);
        TEST_ASSERT(maxd==1);

        pwc1=pwc0;

        data_type ttol = 0.0001;
        pwc1.to_cubic(ttol);

//        std::cout << pwc1.number_segments() << std::endl;
//        octave_print(1, pwc0);
//        octave_print(2, pwc1);

        TEST_ASSERT(pwc1.number_segments()==1);  // Original was linear.

        pwc1.degree(mind, maxd);
        TEST_ASSERT(mind==3);
        TEST_ASSERT(maxd==3);

        t=0;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=0.235;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);

        t=1;
        eval_out=pwc0.f(t);
        eval_ref=pwc1.f(t);
        TEST_ASSERT((eval_out-eval_ref).norm()<ttol);
      }
    }

    void length_test()
    {
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_USING_QD
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
#ifdef ELI_USING_QD
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

    void round_test()
    {
      // create closed piecewise curve
      data_type eps(std::numeric_limits<data__>::epsilon());
#ifdef ELI_USING_QD
      if ( (typeid(data_type)==typeid(dd_real)) || (typeid(data_type)==typeid(qd_real)) )
        eps=std::numeric_limits<double>::epsilon();
#endif
      piecewise_curve_type pwc0, pwc1;
      typename curve_type::control_point_type cntrl_in[3];
      typename piecewise_curve_type::error_code err;
      curve_type bc;
      point_type pt0, pt1, ptref;
      data_type t[4], rad(0.2);
      index_type i, ip1;

      // build piecewise curve
      cntrl_in[0] << 0, 0, 0;
      cntrl_in[1] << 0, 2, 0;
      cntrl_in[2] << 1, 1, 0;
      t[0]=-0.5;
      t[1]=0.5;
      t[2]=1.5;
      t[3]=2.5;
      bc.resize(1);
      pwc0.set_t0(t[0]);
      i=0;
      ip1=1;
      bc.set_control_point(cntrl_in[i], 0);
      bc.set_control_point(cntrl_in[ip1], 1);
      err=pwc0.push_back(bc, t[i+1]-t[i]);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      i=1;
      ip1=2;
      bc.set_control_point(cntrl_in[i], 0);
      bc.set_control_point(cntrl_in[ip1], 1);
      err=pwc0.push_back(bc, t[i+1]-t[i]);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      i=2;
      ip1=0;
      bc.set_control_point(cntrl_in[i], 0);
      bc.set_control_point(cntrl_in[ip1], 1);
      err=pwc0.push_back(bc, t[i+1]-t[i]);
      TEST_ASSERT(err==piecewise_curve_type::NO_ERROR);
      TEST_ASSERT(pwc0.number_segments()==3);
      pwc1=pwc0;

      // round internal joint
      TEST_ASSERT(pwc0.round(rad, 1));
      TEST_ASSERT(pwc0.number_segments()==4);
      pt0=pwc0.f(t[1]);
      pt1=pwc1.f(t[1]);
      ptref << 0.0291414, 1.89247, 0;
      TEST_ASSERT(pt0!=pt1);
      TEST_ASSERT((pt0-ptref).norm()<1e-4);

      // round joint 0
      pwc0=pwc1;
      TEST_ASSERT(pwc0.round(rad, 0));
      pt0=pwc0.f(t[0]);
      pt1=pwc1.f(t[0]);
      ptref << 0.0414214, 0.1, 0;
      data_type v; pwc0.get_parameter_max(v);
      TEST_ASSERT(pt0!=pt1);
      TEST_ASSERT((pt0-ptref).norm()<1e-4);

      // round last joint
      pwc0=pwc1;
      TEST_ASSERT(pwc0.round(rad, 3));
      pt0=pwc0.f(t[3]);
      pt1=pwc1.f(t[3]);
      ptref << 0.0414214, 0.1, 0;
      TEST_ASSERT(pt0!=pt1);
      TEST_ASSERT((pt0-ptref).norm()<1e-4);

      // round all joints
      pwc0=pwc1;
      pwc0.round(rad);
      pt0=pwc0.f(t[0]);
      pt1=pwc1.f(t[0]);
      ptref << 0.0414214, 0.1, 0;
      TEST_ASSERT(pt0!=pt1);
      TEST_ASSERT((pt0-ptref).norm()<1e-4);
      pt0=pwc0.f(t[1]);
      pt1=pwc1.f(t[1]);
      ptref << 0.0291414, 1.89247, 0;
      TEST_ASSERT(pt0!=pt1);
      TEST_ASSERT((pt0-ptref).norm()<1e-4);
      pt0=pwc0.f(t[2]);
      pt1=pwc1.f(t[2]);
      ptref << 0.9171573, 1, 0;
      TEST_ASSERT(pt0!=pt1);
      TEST_ASSERT((pt0-ptref).norm()<1e-4);
      pt0=pwc0.f(t[3]);
      pt1=pwc1.f(t[3]);
      ptref << 0.0414214, 0.1, 0;
      TEST_ASSERT(pt0!=pt1);
      TEST_ASSERT((pt0-ptref).norm()<1e-4);
    }
};
#endif

