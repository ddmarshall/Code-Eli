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

#ifndef eli_mutil_opt_least_squares_hpp
#define eli_mutil_opt_least_squares_hpp

#ifdef Success  // X11 #define collides with Eigen
#undef Success
#endif

#include "Eigen/Eigen"

namespace eli
{
  namespace mutil
  {
    namespace opt
    {
      template<typename data1__, typename data2__, typename data3__>
      void least_squares_uncon(Eigen::MatrixBase<data1__> &x, const Eigen::MatrixBase<data2__> &A, const Eigen::MatrixBase<data3__> &r)
      {
        unsigned int flags(0);
        if (Eigen::MatrixBase<data2__>::ColsAtCompileTime==Eigen::Dynamic)
          flags=Eigen::ComputeThinU | Eigen::ComputeThinV;
        else
          flags=Eigen::ComputeFullU | Eigen::ComputeFullV;

        Eigen::JacobiSVD<typename Eigen::MatrixBase<data2__>::PlainObject > svd(A, flags);
        x=svd.solve(r);
      }

      template<typename data1__, typename data2__, typename data3__, typename data4__, typename data5__>
      void least_squares_eqcon(Eigen::MatrixBase<data1__>  &x, const Eigen::MatrixBase<data2__> &A, const Eigen::MatrixBase<data3__> &b, const Eigen::MatrixBase<data4__> &B, const Eigen::MatrixBase<data5__> &d)
      {
        typedef Eigen::Matrix<typename Eigen::MatrixBase<data4__>::Scalar, Eigen::Dynamic, Eigen::Dynamic> B_matrixD;
        B_matrixD Q, R, temp_B, temp_R;
        typedef Eigen::Matrix<typename Eigen::MatrixBase<data2__>::Scalar, Eigen::Dynamic, Eigen::Dynamic> A_matrixD;
        A_matrixD A1, A2, temp_A;
        typedef Eigen::Matrix<typename Eigen::MatrixBase<data5__>::Scalar, Eigen::Dynamic, Eigen::Dynamic> d_matrixD;
        d_matrixD y;
        typedef Eigen::Matrix<typename Eigen::MatrixBase<data3__>::Scalar, Eigen::Dynamic, Eigen::Dynamic> b_matrixD;
        b_matrixD z, rhs;
        typedef Eigen::Matrix<typename Eigen::MatrixBase<data1__>::Scalar, Eigen::Dynamic, Eigen::Dynamic> x_matrixD;
        x_matrixD x_temp(A.cols(), b.cols());
        typename A_matrixD::Index p(B.rows()), n(A.cols());
  #ifdef DEBUG
        typename A_matrixD::Index m(b.rows());
  #endif

        // build Q and R
        Eigen::HouseholderQR<B_matrixD> qr(B.transpose());
        Q=qr.householderQ();
        temp_B=qr.matrixQR();
        temp_R=Eigen::TriangularView<B_matrixD, Eigen::Upper>(temp_B);
        R=temp_R.topRows(p);
        assert((R.rows()==p) && (R.cols()==p));
        assert((Q.rows()==n) && (Q.cols()==n));

        // build A1 and A2
        temp_A=A*Q;
        A1=temp_A.leftCols(p);
        A2=temp_A.rightCols(n-p);
#ifdef DEBUG
        assert((A1.rows()==m) && (A1.cols()==p));
        assert((A2.rows()==m) && (A2.cols()==n-p));
#endif
        assert(A1.cols()==p);
        assert(A2.cols()==n-p);

        // solve for y
        y=R.transpose().lu().solve(d);

        // setup the unconstrained optimization
        rhs=b-A1*y;
        least_squares_uncon(z, A2, rhs);

        // build the solution
        x_temp.topRows(p)=y;
        x_temp.bottomRows(n-p)=z;
        x=Q*x_temp;
      }
    }
  }
}

#endif
