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

#ifndef least_squares_test_suite_hpp
#define least_squares_test_suite_hpp

#include "Eigen/Eigen"

#include "eli/code_eli.hpp"

#include "eli/opt/least_squares.hpp"

#include <cmath>      // std::sqrt()
#include <limits>     // std::numeric_limits<>

template<typename data__>
class least_squares_test_suite : public Test::Suite
{
  protected:
    void AddTests(const float &)
    {
      TEST_ADD(least_squares_test_suite<float>::least_squares_unc_test);
      TEST_ADD(least_squares_test_suite<float>::least_squares_eqcon_test);
      TEST_ADD(least_squares_test_suite<float>::least_squares_ineqcon_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(least_squares_test_suite<double>::least_squares_unc_test);
      TEST_ADD(least_squares_test_suite<double>::least_squares_eqcon_test);
      TEST_ADD(least_squares_test_suite<double>::least_squares_ineqcon_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(least_squares_test_suite<long double>::least_squares_unc_test);
      TEST_ADD(least_squares_test_suite<long double>::least_squares_eqcon_test);
      TEST_ADD(least_squares_test_suite<long double>::least_squares_ineqcon_test);
    }

#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      TEST_ADD(least_squares_test_suite<dd_real>::least_squares_unc_test);
      TEST_ADD(least_squares_test_suite<dd_real>::least_squares_eqcon_test);
      TEST_ADD(least_squares_test_suite<dd_real>::least_squares_ineqcon_test);
    }

    void AddTests(const qd_real &)
    {
      TEST_ADD(least_squares_test_suite<qd_real>::least_squares_unc_test);
      TEST_ADD(least_squares_test_suite<qd_real>::least_squares_eqcon_test);
      TEST_ADD(least_squares_test_suite<qd_real>::least_squares_ineqcon_test);
    }
#endif

  public:
    least_squares_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~least_squares_test_suite()
    {
    }

  private:
    void least_squares_unc_test()
    {
      Eigen::Matrix<data__, 4, 3> A;
      Eigen::Matrix<data__, 4, 1> r;
      Eigen::Matrix<data__, 3, 1> x, x_ans;
      data__ delta(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // set the coefficient matrix
      A(0,0)=2.0;
      A(0,1)=-1.0;
      A(0,2)=1.0;
      A(1,0)=1.0;
      A(1,1)=-5.0;
      A(1,2)=2.0;
      A(2,0)=-3.0;
      A(2,1)=1.0;
      A(2,2)=-4.0;
      A(3,0)=1.0;
      A(3,1)=-1.0;
      A(3,2)=1.0;

      // set the rhs vector and answer vector
      r(0)=-4.0;
      r(1)=2.0;
      r(2)=5.0;
      r(3)=-1.0;
      x_ans(0)=-static_cast<data__>(18.0)/7.0;
      x_ans(1)=-static_cast<data__>(151.0)/210.0;
      x_ans(2)=static_cast<data__>(107.0)/210.0;

      // calculate the answer and compare
      eli::opt::least_squares_uncon(x, A, r);
      TEST_ASSERT((x_ans-x).norm()<=delta);

      // set the rhs matrix and answer matrix
      Eigen::Matrix<data__, 4, 3> rmat;
      Eigen::Matrix<data__, 3, 3> xmat, xmat_ans;
      rmat.col(0)=r;
      rmat.col(1)=r;
      rmat.col(2)=r;
      xmat_ans.col(0)=x_ans;
      xmat_ans.col(1)=x_ans;
      xmat_ans.col(2)=x_ans;

      // calculate the answer and compare
      eli::opt::least_squares_uncon(xmat, A, rmat);
      TEST_ASSERT((xmat_ans-xmat).norm()<=delta);
    }

    void least_squares_eqcon_test()
    {
      Eigen::Matrix<data__, 5, 4> A;
      Eigen::Matrix<data__, 3, 4> B;
      Eigen::Matrix<data__, 5, 1> b;
      Eigen::Matrix<data__, 3, 1> d;
      Eigen::Matrix<data__, 4, 1> x, x_ans;
      data__ delta(std::sqrt(std::numeric_limits<data__>::epsilon()));

      // set coefficient matrices
      A(0,0)=1.0;
      A(0,1)=1.0;
      A(0,2)=1.0;
      A(0,3)=1.0;
      A(1,0)=1.0;
      A(1,1)=3.0;
      A(1,2)=1.0;
      A(1,3)=1.0;
      A(2,0)=1.0;
      A(2,1)=-1.0;
      A(2,2)=3.0;
      A(2,3)=1.0;
      A(3,0)=1.0;
      A(3,1)=1.0;
      A(3,2)=1.0;
      A(3,3)=3.0;
      A(4,0)=1.0;
      A(4,1)=1.0;
      A(4,2)=1.0;
      A(4,3)=-1.0;
      B(0,0)=1.0;
      B(0,1)=1.0;
      B(0,2)=1.0;
      B(0,3)=-1.0;
      B(1,0)=1.0;
      B(1,1)=-1.0;
      B(1,2)=1.0;
      B(1,3)=1.0;
      B(2,0)=1.0;
      B(2,1)=1.0;
      B(2,2)=-1.0;
      B(2,3)=1.0;

      // set right hand side vectors and solution vector
      b(0,0)=2.0;
      b(1,0)=1.0;
      b(2,0)=6.0;
      b(3,0)=3.0;
      b(4,0)=1.0;
      d(0,0)=1.0;
      d(1,0)=3.0;
      d(2,0)=-1.0;
      x_ans(0,0)=0.5;
      x_ans(1,0)=-0.5;
      x_ans(2,0)=1.5;
      x_ans(3,0)=0.5;

      // calculate the answer and compare
      eli::opt::least_squares_eqcon(x, A, b, B, d);
      TEST_ASSERT((x_ans-x).norm()<=delta);

      // set the rhs matrix and answer matrix
      Eigen::Matrix<data__, 5, 3> bmat;
      Eigen::Matrix<data__, 3, 3> dmat;
      Eigen::Matrix<data__, 4, 3> xmat, xmat_ans;
      bmat.col(0)=b;
      bmat.col(1)=b;
      bmat.col(2)=b;
      dmat.col(0)=d;
      dmat.col(1)=d;
      dmat.col(2)=d;
      xmat_ans.col(0)=x_ans;
      xmat_ans.col(1)=x_ans;
      xmat_ans.col(2)=x_ans;

      // calculate the answer and compare
      eli::opt::least_squares_eqcon(xmat, A, bmat, B, dmat);
      TEST_ASSERT((xmat_ans-xmat).norm()<=delta);
    }

    void least_squares_ineqcon_test()
    {
    }
};

#endif
