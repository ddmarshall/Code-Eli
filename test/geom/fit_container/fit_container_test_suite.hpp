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

#ifndef fit_container_test_suite_hpp
#define fit_container_test_suite_hpp

#include <cassert>  // assert()
#include <cmath>    // cos(), sin()

#include <typeinfo> // typeid

#include "Eigen/Eigen"

#include "eli/code_eli.hpp"

#include "eli/constants/math.hpp"
#include "eli/geom/curve/fit_container.hpp"

template<typename data__>
class fit_container_test_suite : public Test::Suite
{
  private:
    typedef Eigen::Matrix<data__, 1, 3> vector_type;
    typedef eli::geom::curve::fit_container<data__, int, 3, 3> fit_container_type;
    typedef typename fit_container_type::constraint_info constraint_info;
    typedef typename fit_container_type::point_type point_type;
    typedef typename fit_container_type::error_code error_code;

  protected:
    void AddTests(const float &)
    {
      // add the tests
      TEST_ADD(fit_container_test_suite<float>::construction_test);
      TEST_ADD(fit_container_test_suite<float>::add_C0_constraints_test);
      TEST_ADD(fit_container_test_suite<float>::add_C1_constraints_test);
      TEST_ADD(fit_container_test_suite<float>::add_C2_constraints_test);
      TEST_ADD(fit_container_test_suite<float>::add_end_constraints_test);
      TEST_ADD(fit_container_test_suite<float>::remove_constraints_test);
      TEST_ADD(fit_container_test_suite<float>::list_constraints_test);
      TEST_ADD(fit_container_test_suite<float>::closed_constraints_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(fit_container_test_suite<double>::construction_test);
      TEST_ADD(fit_container_test_suite<double>::add_C0_constraints_test);
      TEST_ADD(fit_container_test_suite<double>::add_C1_constraints_test);
      TEST_ADD(fit_container_test_suite<double>::add_C2_constraints_test);
      TEST_ADD(fit_container_test_suite<double>::add_end_constraints_test);
      TEST_ADD(fit_container_test_suite<double>::remove_constraints_test);
      TEST_ADD(fit_container_test_suite<double>::list_constraints_test);
      TEST_ADD(fit_container_test_suite<double>::closed_constraints_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(fit_container_test_suite<long double>::construction_test);
      TEST_ADD(fit_container_test_suite<long double>::add_C0_constraints_test);
      TEST_ADD(fit_container_test_suite<long double>::add_C1_constraints_test);
      TEST_ADD(fit_container_test_suite<long double>::add_C2_constraints_test);
      TEST_ADD(fit_container_test_suite<long double>::add_end_constraints_test);
      TEST_ADD(fit_container_test_suite<long double>::remove_constraints_test);
      TEST_ADD(fit_container_test_suite<long double>::list_constraints_test);
      TEST_ADD(fit_container_test_suite<long double>::closed_constraints_test);
    }

  public:
    fit_container_test_suite()
    {
      AddTests(data__());
    }
    ~fit_container_test_suite()
    {
    }

  private:
    template<typename it__>
    void create_points(it__ it, size_t npts)
    {
      for (size_t i=0; i<npts; ++i, ++it)
      {
        typename it__::value_type p(3);
        p(0)=std::cos(2*eli::constants::math<data__>::pi()*static_cast<data__>(i)/npts);
        p(1)=std::sin(2*eli::constants::math<data__>::pi()*static_cast<data__>(i)/npts);
        p(2)=0;
        if (i+6>npts)
          p(0)+=0.5;
        (*it)=p;
      }
    }

    void construction_test()
    {
      fit_container_type ccon;
      std::vector<point_type> points(20);

      // set points
      create_points(points.begin(), points.size());
      ccon.set_points(points.begin(), points.end());

      TEST_ASSERT(ccon.number_points()==20);
      TEST_ASSERT(ccon.number_constraint_points()==0);
      TEST_ASSERT(ccon.number_constraints()==0);
    }

    void add_C0_constraints_test()
    {
      int index1(2), index2(5);
      fit_container_type ccon;
      error_code ec;
      constraint_info ciout;
      std::vector<point_type> points(20);

      // set points
      create_points(points.begin(), points.size());
      ccon.set_points(points.begin(), points.end());

      // add constraint
      ec=ccon.add_C0_constraint(index1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==1);
      TEST_ASSERT(ccon.number_constraints()==1);

      // check for the constraint
      ec=ccon.get_constraint(index1, ciout);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::NOT_USED);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // check for constraint that does not exist
      ec=ccon.get_constraint(index2, ciout);
      TEST_ASSERT(ec==fit_container_type::INDEX_NOT_FOUND);

      // add another constraint
      ec=ccon.add_C0_constraint(index2);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==2);

      // check for the constraint
      ec=ccon.get_constraint(index2, ciout);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::NOT_USED);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // replace constraint
      ec=ccon.add_C0_constraint(index1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==2);

      // check for the constraint
      ec=ccon.get_constraint(index1, ciout);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::NOT_USED);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);
    }

    void add_C1_constraints_test()
    {
      int index1(2), index2(5), index3(4);
      vector_type v(3), vout(3);
      fit_container_type ccon;
      error_code ec;
      constraint_info ciout;
      std::vector<point_type> points(20);

      // set points
      create_points(points.begin(), points.size());
      ccon.set_points(points.begin(), points.end());

      // add constraint
      v[0]=1; v[1]=2; v[2]=3;
      ec=ccon.add_C1_constraint(index1, v);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==1);
      TEST_ASSERT(ccon.number_constraints()==2);

      // check for the constraint
      ec=ccon.get_constraint(index1, ciout);
      vout=ciout.get_fp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((vout[0]==v[0]) && (vout[1]==v[1]) && (vout[2]==v[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // check for constraint that does not exist
      ec=ccon.get_constraint(index2, ciout);
      TEST_ASSERT(ec==fit_container_type::INDEX_NOT_FOUND);

      // add another constraint
      ec=ccon.add_C1_constraint(index2, v);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==4);

      // check for the constraint
      ec=ccon.get_constraint(index2, ciout);
      vout=ciout.get_fp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((vout[0]==v[0]) && (vout[1]==v[1]) && (vout[2]==v[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // replace constraint
      ec=ccon.add_C1_constraint(index1, v);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==4);

      // check for the constraint
      ec=ccon.get_constraint(index1, ciout);
      vout=ciout.get_fp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((vout[0]==v[0]) && (vout[1]==v[1]) && (vout[2]==v[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // add finite difference constraint
      ec=ccon.add_C1_constraint(index3);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==3);
      TEST_ASSERT(ccon.number_constraints()==6);

      // check for the constraint
      data__ t=2*eli::constants::math<data__>::pi()*static_cast<data__>(index3)/points.size();
      v[0]=-std::sin(t);
      v[1]=std::cos(t);
      v[2]=0;
      ec=ccon.get_constraint(index3, ciout);
      vout=ciout.get_fp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
      TEST_ASSERT_DELTA(vout[0], v[0], 2e-2);
      TEST_ASSERT_DELTA(vout[1], v[1], 5e-3);
      TEST_ASSERT(vout[2]==v[2]);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);
    }

    void add_C2_constraints_test()
    {
      int index1(2), index2(5), index3(4), index4(6);
      vector_type v1(3), v2(3), v1out(3), v2out(3);
      fit_container_type ccon;
      error_code ec;
      constraint_info ciout;
      std::vector<point_type> points(20);

      // set points
      create_points(points.begin(), points.size());
      ccon.set_points(points.begin(), points.end());

      // add constraint
      v1[0]=1; v1[1]=2; v1[2]=3;
      v2[0]=0; v2[1]=-1; v2[2]=-2;
      ec=ccon.add_C2_constraint(index1, v1, v2);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==1);
      TEST_ASSERT(ccon.number_constraints()==3);

      // check for the constraint
      ec=ccon.get_constraint(index1, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);
      TEST_ASSERT((v2out[0]==v2[0]) && (v2out[1]==v2[1]) && (v2out[2]==v2[2]));

      // check for constraint that does not exist
      ec=ccon.get_constraint(index2, ciout);
      TEST_ASSERT(ec==fit_container_type::INDEX_NOT_FOUND);

      // add another constraint
      ec=ccon.add_C2_constraint(index2, v1, v2);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==6);

      // check for the constraint
      ec=ccon.get_constraint(index2, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);
      TEST_ASSERT((v2out[0]==v2[0]) && (v2out[1]==v2[1]) && (v2out[2]==v2[2]));

      // replace constraint
      ec=ccon.add_C2_constraint(index1, v1, v2);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==6);

      // check for the constraint
      ec=ccon.get_constraint(index1, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);
      TEST_ASSERT((v2out[0]==v2[0]) && (v2out[1]==v2[1]) && (v2out[2]==v2[2]));

      // add finite difference 2nd derivative constraint
      ec=ccon.add_C2_constraint(index3, v1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==3);
      TEST_ASSERT(ccon.number_constraints()==9);

      // check for the constraint
      data__ t=2*eli::constants::math<data__>::pi()*static_cast<data__>(index3)/points.size();
      v2[0]=-std::cos(t);
      v2[1]=-std::sin(t);
      v2[2]=0;
      ec=ccon.get_constraint(index3, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v2out[0], v2[0], 4*std::numeric_limits<data__>::epsilon());
      TEST_ASSERT_DELTA(v2out[1], v2[1], 11*std::numeric_limits<data__>::epsilon());
      TEST_ASSERT(v2out[2]==v2[2]);

      // add finite difference 1st and 2nd derivative constraint
      ec=ccon.add_C2_constraint(index4);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==4);
      TEST_ASSERT(ccon.number_constraints()==12);

      // check for the constraint
      t=2*eli::constants::math<data__>::pi()*static_cast<data__>(index4)/points.size();
      v1[0]=-std::sin(t);
      v1[1]=std::cos(t);
      v1[2]=0;
      v2[0]=-std::cos(t);
      v2[1]=-std::sin(t);
      v2[2]=0;
      ec=ccon.get_constraint(index4, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v1out[0], v1[0], 2e-2);
      TEST_ASSERT_DELTA(v1out[1], v1[1], 5e-3);
      TEST_ASSERT(v1out[2]==v1[2]);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v2out[0], v2[0], 11*std::numeric_limits<data__>::epsilon());
      TEST_ASSERT_DELTA(v2out[1], v2[1], 15*std::numeric_limits<data__>::epsilon());
      TEST_ASSERT(v2out[2]==v2[2]);
    }

    void add_end_constraints_test()
    {
      int index;
      vector_type v1(3), v2(3), v1out(3), v2out(3);
      fit_container_type ccon;
      error_code ec;
      constraint_info ciout;
      std::vector<point_type> points(20);
      data__ ts(0), te(2*eli::constants::math<data__>::pi()*static_cast<data__>(points.size()-1)/points.size());

      v1[0]=1; v1[1]=2; v1[2]=3;
      v2[0]=0; v2[1]=-1; v2[2]=-2;

      // set points
      create_points(points.begin(), points.size());
      ccon.set_points(points.begin(), points.end());

      // add C0 constraint at start
      ec=ccon.add_start_C0_constraint();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==1);
      TEST_ASSERT(ccon.number_constraints()==1);

      // check for the constraint
      index=0;
      ec=ccon.get_constraint(index, ciout);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::NOT_USED);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // add C0 constraint at end
      ec=ccon.add_end_C0_constraint();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==2);

      // check for the constraint
      index=static_cast<int>(ccon.number_points())-1;
      ec=ccon.get_constraint(index, ciout);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::NOT_USED);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // add C1 constraint at start
      ec=ccon.add_start_C1_constraint(v1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==3);

      // check for the constraint
      index=0;
      ec=ccon.get_constraint(index, ciout);
      v1out=ciout.get_fp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // add C1 constraint at end
      ec=ccon.add_end_C1_constraint(v1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==4);

      // check for the constraint
      index=static_cast<int>(ccon.number_points())-1;
      ec=ccon.get_constraint(index, ciout);
      v1out=ciout.get_fp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // add C2 constraint at start
      ec=ccon.add_start_C2_constraint(v1, v2);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==5);

      // check for the constraint
      index=0;
      ec=ccon.get_constraint(index, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);
      TEST_ASSERT((v2out[0]==v2[0]) && (v2out[1]==v2[1]) && (v2out[2]==v2[2]));

      // add C2 constraint at end
      ec=ccon.add_end_C2_constraint(v1, v2);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==6);

      // check for the constraint
      index=static_cast<int>(ccon.number_points())-1;
      ec=ccon.get_constraint(index, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);
      TEST_ASSERT((v2out[0]==v2[0]) && (v2out[1]==v2[1]) && (v2out[2]==v2[2]));

      // add C1 start finite difference constraint
      ec=ccon.add_start_C1_constraint();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==5);

      // check for the constraint
      index=0;
      v1[0]=-std::sin(ts);
      v1[1]=std::cos(ts);
      v1[2]=0;
      ec=ccon.get_constraint(index, ciout);
      v1out=ciout.get_fp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v1out[0], v1[0], 8e-3);
      TEST_ASSERT_DELTA(v1out[1], v1[1], 4e-2);
      TEST_ASSERT(v1out[2]==v1[2]);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // add C1 end finite difference constraint
      ec=ccon.add_end_C1_constraint();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==4);

      // check for the constraint
      index=static_cast<int>(ccon.number_points())-1;
      v1[0]=-std::sin(te);
      v1[1]=std::cos(te);
      v1[2]=0;
      ec=ccon.get_constraint(index, ciout);
      v1out=ciout.get_fp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v1out[0], v1[0], 2e-2);
      TEST_ASSERT_DELTA(v1out[1], v1[1], 4e-2);
      TEST_ASSERT(v1out[2]==v1[2]);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

      // add C2 start finite difference 2nd derivative constraint
      v1[0]=1; v1[1]=2; v1[2]=3;
      ec=ccon.add_start_C2_constraint(v1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==5);

      // check for the constraint
      index=0;
      v2[0]=-std::cos(ts);
      v2[1]=-std::sin(ts);
      v2[2]=0;
      ec=ccon.get_constraint(index, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v2out[0], v2[0], 1e-1);
      TEST_ASSERT_DELTA(v2out[1], v2[1], 4e-2);
      TEST_ASSERT(v2out[2]==v2[2]);

      // add C2 end finite difference 2nd derivative constraint
      v1[0]=1; v1[1]=2; v1[2]=3;
      ec=ccon.add_end_C2_constraint(v1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==6);

      // check for the constraint
      index=static_cast<int>(ccon.number_points())-1;
      v2[0]=-std::cos(te);
      v2[1]=-std::sin(te);
      v2[2]=0;
      ec=ccon.get_constraint(index, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
      TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
      TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v2out[0], v2[0], 8e-2);
      TEST_ASSERT_DELTA(v2out[1], v2[1], 6e-2);
      TEST_ASSERT(v2out[2]==v2[2]);

      // add C2 start finite difference 1st and 2nd derivative constraint
      ec=ccon.add_start_C2_constraint();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==6);

      // check for the constraint
      index=0;
      v1[0]=-std::sin(ts);
      v1[1]=std::cos(ts);
      v1[2]=0;
      v2[0]=-std::cos(ts);
      v2[1]=-std::sin(ts);
      v2[2]=0;
      ec=ccon.get_constraint(index, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v1out[0], v1[0], 8e-3);
      TEST_ASSERT_DELTA(v1out[1], v1[1], 4e-2);
      TEST_ASSERT(v1out[2]==v1[2]);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v2out[0], v2[0], 1e-1);
      TEST_ASSERT_DELTA(v2out[1], v2[1], 4e-2);
      TEST_ASSERT(v2out[2]==v2[2]);

      // add C2 end finite difference 1st and 2nd derivative constraint
      ec=ccon.add_end_C2_constraint();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==2);
      TEST_ASSERT(ccon.number_constraints()==6);

      // check for the constraint
      index=static_cast<int>(ccon.number_points())-1;
      v1[0]=-std::sin(te);
      v1[1]=std::cos(te);
      v1[2]=0;
      v2[0]=-std::cos(te);
      v2[1]=-std::sin(te);
      v2[2]=0;
      ec=ccon.get_constraint(index, ciout);
      v1out=ciout.get_fp();
      v2out=ciout.get_fpp();
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v1out[0], v1[0], 2e-2);
      TEST_ASSERT_DELTA(v1out[1], v1[1], 4e-2);
      TEST_ASSERT(v1out[2]==v1[2]);
      TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
      TEST_ASSERT_DELTA(v2out[0], v2[0], 8e-2);
      TEST_ASSERT_DELTA(v2out[1], v2[1], 6e-2);
      TEST_ASSERT(v2out[2]==v2[2]);
    }

    void closed_constraints_test()
    {
      // test no FD constraints
      {
        int index1(2), index;
        vector_type v1(3), v2(3), v1out(3), v2out(3);
        fit_container_type ccon;
        error_code ec;
        constraint_info ciout;
        std::vector<point_type> points(20);

        // set points
        create_points(points.begin(), points.size());
        ccon.set_points(points.begin(), points.end());

        // add constraint
        v1[0]=1; v1[1]=2; v1[2]=3;
        v2[0]=0; v2[1]=-1; v2[2]=-2;
        ec=ccon.add_C2_constraint(index1, v1, v2);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==1);
        TEST_ASSERT(ccon.number_constraints()==3);

        // add C2 constraint at start
        ec=ccon.add_start_C2_constraint(v1, v2);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==2);
        TEST_ASSERT(ccon.number_constraints()==6);

        // add C2 constraint at end
        ec=ccon.add_end_C2_constraint(v1, v2);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==3);
        TEST_ASSERT(ccon.number_constraints()==9);

        // check for the constraint
        ec=ccon.get_constraint(index1, ciout);
        v1out=ciout.get_fp();
        v2out=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);
        TEST_ASSERT((v2out[0]==v2[0]) && (v2out[1]==v2[1]) && (v2out[2]==v2[2]));

        // check for the constraint
        index=0;
        ec=ccon.get_constraint(index, ciout);
        v1out=ciout.get_fp();
        v2out=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);
        TEST_ASSERT((v2out[0]==v2[0]) && (v2out[1]==v2[1]) && (v2out[2]==v2[2]));

        // check for the constraint
        index=static_cast<int>(ccon.number_points())-1;
        ec=ccon.get_constraint(index, ciout);
        v1out=ciout.get_fp();
        v2out=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);
        TEST_ASSERT((v2out[0]==v2[0]) && (v2out[1]==v2[1]) && (v2out[2]==v2[2]));

        // change to closed curve
        ccon.set_end_flag(eli::geom::general::C0);
        TEST_ASSERT(ccon.closed());

        // check for the constraint
        ec=ccon.get_constraint(index1, ciout);
        v1out=ciout.get_fp();
        v2out=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);

        // check for the constraint
        index=0;
        ec=ccon.get_constraint(index, ciout);
        v1out=ciout.get_fp();
        v2out=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);
        TEST_ASSERT((v2out[0]==v2[0]) && (v2out[1]==v2[1]) && (v2out[2]==v2[2]));

        // check for the constraint
        index=static_cast<int>(ccon.number_points())-1;
        ec=ccon.get_constraint(index, ciout);
        v1out=ciout.get_fp();
        v2out=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((v1out[0]==v1[0]) && (v1out[1]==v1[1]) && (v1out[2]==v1[2]));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::SET);
        TEST_ASSERT((v2out[0]==v2[0]) && (v2out[1]==v2[1]) && (v2out[2]==v2[2]));
      }

      // test only FD C1 constraints
      {
        int index0, index1, indexnm2, indexnm1;
        vector_type vp(3), vp0(3), vp1(3), vpnm2(3), vpnm1(3);
        fit_container_type ccon;
        error_code ec;
        constraint_info ciout;
        std::vector<point_type> points(20);
        data__ t0(0), t1(2*eli::constants::math<data__>::pi()/points.size()), tnm2(2*eli::constants::math<data__>::pi()*static_cast<data__>(points.size()-1)/points.size()), tnm1(2*eli::constants::math<data__>::pi()*static_cast<data__>(points.size()-1)/points.size());

        // set points
        create_points(points.begin(), points.size());
        ccon.set_points(points.begin(), points.end());

        // set the reference values
        index0=0;
        index1=1;
        indexnm2=static_cast<int>(points.size())-2;
        indexnm1=static_cast<int>(points.size())-1;
        t0  =2*eli::constants::math<data__>::pi()*static_cast<data__>(index0)/points.size();
        t1  =2*eli::constants::math<data__>::pi()*static_cast<data__>(index1)/points.size();
        tnm2=2*eli::constants::math<data__>::pi()*static_cast<data__>(indexnm2)/points.size();
        tnm1=2*eli::constants::math<data__>::pi()*static_cast<data__>(indexnm1)/points.size();
        vp0(0)=-std::sin(t0);
        vp0(1)=std::cos(t0);
        vp0(2)=0;
        vp1(0)=-std::sin(t1);
        vp1(1)=std::cos(t1);
        vp1(2)=0;
        vpnm2(0)=-std::sin(tnm2);
        vpnm2(1)=std::cos(tnm2);
        vpnm2(2)=0;
        vpnm1(0)=-std::sin(tnm1);
        vpnm1(1)=std::cos(tnm1);
        vpnm1(2)=0;

        // add C1 constraint at start
        ec=ccon.add_start_C1_constraint();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==1);
        TEST_ASSERT(ccon.number_constraints()==2);

        // add C1 constraint at start+1
        ec=ccon.add_C1_constraint(1);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==2);
        TEST_ASSERT(ccon.number_constraints()==4);

        // add C1 constraint at end-1
        ec=ccon.add_C1_constraint(static_cast<int>(points.size())-2);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==3);
        TEST_ASSERT(ccon.number_constraints()==6);

        // add C1 constraint at end
        ec=ccon.add_end_C1_constraint();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==4);
        TEST_ASSERT(ccon.number_constraints()==8);

        // check all derivatives
        ec=ccon.get_constraint(index0, ciout);
        vp=ciout.get_fp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vp0(0), 8e-3);
        TEST_ASSERT_DELTA(vp(1), vp0(1), 4e-2);
        TEST_ASSERT(vp(2)==vp0(2));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);
        ec=ccon.get_constraint(index1, ciout);
        vp=ciout.get_fp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vp1(0), 4e-3);
        TEST_ASSERT_DELTA(vp(1), vp1(1), 2e-2);
        TEST_ASSERT(vp(2)==vp1(2));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);
        ec=ccon.get_constraint(indexnm2, ciout);
        vp=ciout.get_fp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vpnm2(0), 8e-3);
        TEST_ASSERT_DELTA(vp(1), vpnm2(1), 1e-2);
        TEST_ASSERT(vp(2)==vpnm2(2));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);
        ec=ccon.get_constraint(indexnm1, ciout);
        vp=ciout.get_fp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vpnm1(0), 2e-2);
        TEST_ASSERT_DELTA(vp(1), vpnm1(1), 4e-2);
        TEST_ASSERT(vp(2)==vpnm1(2));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);

        // change to closed curve
        ccon.set_end_flag(eli::geom::general::C0);
        TEST_ASSERT(ccon.closed());

        // check all derivatives
        ec=ccon.get_constraint(index0, ciout);
        vp=ciout.get_fp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vp0(0), 4e-1);
        TEST_ASSERT_DELTA(vp(1), vp0(1), 2e-1);
        TEST_ASSERT(vp(2)==vp0(2));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);
        ec=ccon.get_constraint(index1, ciout);
        vp=ciout.get_fp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vp1(0), 4e-3);
        TEST_ASSERT_DELTA(vp(1), vp1(1), 2e-2);
        TEST_ASSERT(vp(2)==vp1(2));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);
        ec=ccon.get_constraint(indexnm2, ciout);
        vp=ciout.get_fp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vpnm2(0), 8e-3);
        TEST_ASSERT_DELTA(vp(1), vpnm2(1), 1e-2);
        TEST_ASSERT(vp(2)==vpnm2(2));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);
        ec=ccon.get_constraint(indexnm1, ciout);
        vp=ciout.get_fp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vpnm1(0), 4e-1);
        TEST_ASSERT_DELTA(vp(1), vpnm1(1), 2e-1);
        TEST_ASSERT(vp(2)==vpnm1(2));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::NOT_USED);
      }

      // test FD C2 & set C1 constraints
      {
        int index0, index1, indexnm2, indexnm1;
        vector_type v1(3), vp(3), vpp(3), vpp0(3), vpp1(3), vppnm2(3), vppnm1(3);
        fit_container_type ccon;
        error_code ec;
        constraint_info ciout;
        std::vector<point_type> points(20);
        data__ t0(0), t1(2*eli::constants::math<data__>::pi()/points.size()), tnm2(2*eli::constants::math<data__>::pi()*static_cast<data__>(points.size()-1)/points.size()), tnm1(2*eli::constants::math<data__>::pi()*static_cast<data__>(points.size()-1)/points.size());

        // set points
        create_points(points.begin(), points.size());
        ccon.set_points(points.begin(), points.end());

        // set the reference values
        v1(0)=1; v1(1)=2; v1(2)=3;
        index0=0;
        index1=1;
        indexnm2=static_cast<int>(points.size())-2;
        indexnm1=static_cast<int>(points.size())-1;
        t0  =2*eli::constants::math<data__>::pi()*static_cast<data__>(index0)/points.size();
        t1  =2*eli::constants::math<data__>::pi()*static_cast<data__>(index1)/points.size();
        tnm2=2*eli::constants::math<data__>::pi()*static_cast<data__>(indexnm2)/points.size();
        tnm1=2*eli::constants::math<data__>::pi()*static_cast<data__>(indexnm1)/points.size();
        vpp0(0)=-std::cos(t0);
        vpp0(1)=-std::sin(t0);
        vpp0(2)=0;
        vpp1(0)=-std::cos(t1);
        vpp1(1)=-std::sin(t1);
        vpp1(2)=0;
        vppnm2(0)=-std::cos(tnm2);
        vppnm2(1)=-std::sin(tnm2);
        vppnm2(2)=0;
        vppnm1(0)=-std::cos(tnm1);
        vppnm1(1)=-std::sin(tnm1);
        vppnm1(2)=0;

        // add C2 constraint at start
        ec=ccon.add_start_C2_constraint(v1);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==1);
        TEST_ASSERT(ccon.number_constraints()==3);

        // add C2 constraint at start+1
        ec=ccon.add_C2_constraint(1, v1);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==2);
        TEST_ASSERT(ccon.number_constraints()==6);

        // add C2 constraint at end-1
        ec=ccon.add_C2_constraint(static_cast<int>(points.size())-2, v1);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==3);
        TEST_ASSERT(ccon.number_constraints()==9);

        // add C2 constraint at end
        ec=ccon.add_end_C2_constraint(v1);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==4);
        TEST_ASSERT(ccon.number_constraints()==12);

        // check all derivatives
        ec=ccon.get_constraint(index0, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((vp(0)==v1(0)) && (vp(1)==v1(1)) && (vp(2)==v1(2)));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vpp0(0), 1e-1);
        TEST_ASSERT_DELTA(vpp(1), vpp0(1), 4e-2);
        TEST_ASSERT(vpp(2)==vpp0(2));
        ec=ccon.get_constraint(index1, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((vp(0)==v1(0)) && (vp(1)==v1(1)) && (vp(2)==v1(2)));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vpp1(0), 16*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(vpp(1), vpp1(1), 11*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT(vpp(2)==vpp1(2));
        ec=ccon.get_constraint(indexnm2, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((vp(0)==v1(0)) && (vp(1)==v1(1)) && (vp(2)==v1(2)));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vppnm2(0), 33*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(vpp(1), vppnm2(1), 33*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT(vpp(2)==vppnm2(2));
        ec=ccon.get_constraint(indexnm1, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((vp(0)==v1(0)) && (vp(1)==v1(1)) && (vp(2)==v1(2)));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vppnm1(0), 8e-2);
        TEST_ASSERT_DELTA(vpp(1), vppnm1(1), 6e-2);
        TEST_ASSERT(vpp(2)==vppnm1(2));

        // change to closed curve
        ccon.set_end_flag(eli::geom::general::C0);
        TEST_ASSERT(ccon.closed());

        // check all derivatives
        ec=ccon.get_constraint(index0, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((vp(0)==v1(0)) && (vp(1)==v1(1)) && (vp(2)==v1(2)));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vpp0(0), 3.5);
        TEST_ASSERT_DELTA(vpp(1), vpp0(1), 1.4);
        TEST_ASSERT(vpp(2)==vpp0(2));
        ec=ccon.get_constraint(index1, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((vp(0)==v1(0)) && (vp(1)==v1(1)) && (vp(2)==v1(2)));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vpp1(0), 16*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(vpp(1), vpp1(1), 16*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT(vpp(2)==vpp1(2));
        ec=ccon.get_constraint(indexnm2, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((vp(0)==v1(0)) && (vp(1)==v1(1)) && (vp(2)==v1(2)));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vppnm2(0), 33*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(vpp(1), vppnm2(1), 33*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT(vpp(2)==vppnm2(2));
        ec=ccon.get_constraint(indexnm1, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::SET);
        TEST_ASSERT((vp(0)==v1(0)) && (vp(1)==v1(1)) && (vp(2)==v1(2)));
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vppnm1(0), 3.0);
        TEST_ASSERT_DELTA(vpp(1), vppnm1(1), 1.5);
        TEST_ASSERT(vpp(2)==vppnm1(2));
      }

      // test FD C1 and C2 constraints
      {
        int index0, index1, indexnm2, indexnm1;
        vector_type vp(3), vp0(3), vp1(3), vpnm2(3), vpnm1(3), vpp(3), vpp0(3), vpp1(3), vppnm2(3), vppnm1(3);
        fit_container_type ccon;
        error_code ec;
        constraint_info ciout;
        std::vector<point_type> points(20);
        data__ t0(0), t1(2*eli::constants::math<data__>::pi()/points.size()), tnm2(2*eli::constants::math<data__>::pi()*static_cast<data__>(points.size()-1)/points.size()), tnm1(2*eli::constants::math<data__>::pi()*static_cast<data__>(points.size()-1)/points.size());

        // set points
        create_points(points.begin(), points.size());
        ccon.set_points(points.begin(), points.end());

        // set the reference values
        index0=0;
        index1=1;
        indexnm2=static_cast<int>(points.size())-2;
        indexnm1=static_cast<int>(points.size())-1;
        t0  =2*eli::constants::math<data__>::pi()*static_cast<data__>(index0)/points.size();
        t1  =2*eli::constants::math<data__>::pi()*static_cast<data__>(index1)/points.size();
        tnm2=2*eli::constants::math<data__>::pi()*static_cast<data__>(indexnm2)/points.size();
        tnm1=2*eli::constants::math<data__>::pi()*static_cast<data__>(indexnm1)/points.size();
        vp0(0)=-std::sin(t0);
        vp0(1)=std::cos(t0);
        vp0(2)=0;
        vp1(0)=-std::sin(t1);
        vp1(1)=std::cos(t1);
        vp1(2)=0;
        vpnm2(0)=-std::sin(tnm2);
        vpnm2(1)=std::cos(tnm2);
        vpnm2(2)=0;
        vpnm1(0)=-std::sin(tnm1);
        vpnm1(1)=std::cos(tnm1);
        vpnm1(2)=0;
        vpp0(0)=-std::cos(t0);
        vpp0(1)=-std::sin(t0);
        vpp0(2)=0;
        vpp1(0)=-std::cos(t1);
        vpp1(1)=-std::sin(t1);
        vpp1(2)=0;
        vppnm2(0)=-std::cos(tnm2);
        vppnm2(1)=-std::sin(tnm2);
        vppnm2(2)=0;
        vppnm1(0)=-std::cos(tnm1);
        vppnm1(1)=-std::sin(tnm1);
        vppnm1(2)=0;

        // add C2 constraint at start
        ec=ccon.add_start_C2_constraint();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==1);
        TEST_ASSERT(ccon.number_constraints()==3);

        // add C2 constraint at start+1
        ec=ccon.add_C2_constraint(1);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==2);
        TEST_ASSERT(ccon.number_constraints()==6);

        // add C2 constraint at end-1
        ec=ccon.add_C2_constraint(static_cast<int>(points.size())-2);
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==3);
        TEST_ASSERT(ccon.number_constraints()==9);

        // add C2 constraint at end
        ec=ccon.add_end_C2_constraint();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ccon.number_constraint_points()==4);
        TEST_ASSERT(ccon.number_constraints()==12);

        // check all derivatives
        ec=ccon.get_constraint(index0, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vp0(0), 8e-3);
        TEST_ASSERT_DELTA(vp(1), vp0(1), 4e-2);
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vpp0(0), 1e-1);
        TEST_ASSERT_DELTA(vpp(1), vpp0(1), 4e-2);
        TEST_ASSERT(vpp(2)==vpp0(2));
        ec=ccon.get_constraint(index1, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vp1(0), 4e-3);
        TEST_ASSERT_DELTA(vp(1), vp1(1), 2e-2);
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vpp1(0), 16*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(vpp(1), vpp1(1), 11*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT(vpp(2)==vpp1(2));
        ec=ccon.get_constraint(indexnm2, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vpnm2(0), 8e-3);
        TEST_ASSERT_DELTA(vp(1), vpnm2(1), 1e-2);
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vppnm2(0), 33*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(vpp(1), vppnm2(1), 33*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT(vpp(2)==vppnm2(2));
        ec=ccon.get_constraint(indexnm1, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vpnm1(0), 2e-2);
        TEST_ASSERT_DELTA(vp(1), vpnm1(1), 4e-2);
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vppnm1(0), 8e-2);
        TEST_ASSERT_DELTA(vpp(1), vppnm1(1), 6e-2);
        TEST_ASSERT(vpp(2)==vppnm1(2));

        // change to closed curve
        ccon.set_end_flag(eli::geom::general::C0);
        TEST_ASSERT(ccon.closed());

        // check all derivatives
        ec=ccon.get_constraint(index0, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vp0(0), 4e-1);
        TEST_ASSERT_DELTA(vp(1), vp0(1), 2e-1);
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vpp0(0), 3.5);
        TEST_ASSERT_DELTA(vpp(1), vpp0(1), 1.4);
        TEST_ASSERT(vpp(2)==vpp0(2));
        ec=ccon.get_constraint(index1, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vp1(0), 4e-3);
        TEST_ASSERT_DELTA(vp(1), vp1(1), 2e-2);
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vpp1(0), 16*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(vpp(1), vpp1(1), 16*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT(vpp(2)==vpp1(2));
        ec=ccon.get_constraint(indexnm2, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vpnm2(0), 8e-3);
        TEST_ASSERT_DELTA(vp(1), vpnm2(1), 1e-2);
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vppnm2(0), 33*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(vpp(1), vppnm2(1), 33*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT(vpp(2)==vppnm2(2));
        ec=ccon.get_constraint(indexnm1, ciout);
        vp=ciout.get_fp();
        vpp=ciout.get_fpp();
        TEST_ASSERT(ec==fit_container_type::NO_ERROR);
        TEST_ASSERT(ciout.using_fp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vp(0), vpnm1(0), 4e-1);
        TEST_ASSERT_DELTA(vp(1), vpnm1(1), 2e-1);
        TEST_ASSERT(ciout.using_fpp()==constraint_info::FD);
        TEST_ASSERT_DELTA(vpp(0), vppnm1(0), 3.0);
        TEST_ASSERT_DELTA(vpp(1), vppnm1(1), 1.5);
        TEST_ASSERT(vpp(2)==vppnm1(2));
      }
    }

    void remove_constraints_test()
    {
      int index1(2), index2(5);
      vector_type v1(3), v2(3);
      fit_container_type ccon;
      error_code ec;
      constraint_info ciout;
      std::vector<point_type> points(20);

      // set points
      create_points(points.begin(), points.size());
      ccon.set_points(points.begin(), points.end());

      // add some constraints
      ec=ccon.add_C0_constraint(index1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      ec=ccon.add_C0_constraint(index2);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      ec=ccon.add_C1_constraint(index2-index1, v1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      ec=ccon.add_C2_constraint(index2+index1, v1, v2);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==4);
      TEST_ASSERT(ccon.number_constraints()==7);

      // remove constraint
      ec=ccon.remove_constraint(index2-index1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==3);
      TEST_ASSERT(ccon.number_constraints()==5);
      ec=ccon.get_constraint(index2-index1, ciout);
      TEST_ASSERT(ec==fit_container_type::INDEX_NOT_FOUND);

      // try to remove constraint that is not there
      ec=ccon.remove_constraint(index2-1);
      TEST_ASSERT(ec==fit_container_type::INDEX_NOT_FOUND);
    }

    void list_constraints_test()
    {
      int index1(2), index2(5);
      vector_type v1(3), v2(3);
      fit_container_type ccon;
      error_code ec;
      constraint_info ciout;
      std::vector<point_type> points(20);

      v1[0]=0; v2[1]=2; v2[2]=4;
      v2[0]=1; v2[1]=3; v2[2]=5;

      // set points
      create_points(points.begin(), points.size());
      ccon.set_points(points.begin(), points.end());

      // add some constraints
      ec=ccon.add_C0_constraint(index1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      ec=ccon.add_C0_constraint(index2);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      ec=ccon.add_C1_constraint(index2-index1, v1);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      ec=ccon.add_C2_constraint(index2+index1, v1, v2);
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(ccon.number_constraint_points()==4);
      TEST_ASSERT(ccon.number_constraints()==7);

      // get list of constraint indexes
      std::vector<int> indexes(ccon.number_constraint_points());
      ec=ccon.get_constraint_indexes(indexes.begin());
      TEST_ASSERT(ec==fit_container_type::NO_ERROR);
      TEST_ASSERT(indexes[0]==index1);
      TEST_ASSERT(indexes[1]==(index2-index1));
      TEST_ASSERT(indexes[2]==index2);
      TEST_ASSERT(indexes[3]==index2+index1);
    }
};

#endif

