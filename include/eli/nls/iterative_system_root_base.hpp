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

#ifndef eli_nls_iterative_root_system_method_hpp
#define eli_nls_iterative_root_system_method_hpp

#include "Eigen/Eigen"

#include "eli/nls/iterative_root_base.hpp"

namespace eli
{
  namespace nls
  {
    template<typename data__, size_t N__, size_t NSOL__>
    class iterative_system_root_base : public iterative_root_base<data__>
    {
      public:
        typedef Eigen::Matrix<data__, N__, NSOL__> solution_matrix;
        typedef Eigen::Matrix<data__, N__, N__>    jacobian_matrix;

        enum system_norm_type
        {
          L1=100,            // for vectors this is sum of maxes, for matrices this is maximum absolute column sum
          L2=200,            // for vectors this is sqrt of sum of squares, for matrices this will be equivalent to the Frobenius norm
          Linf=300,          // for vectors this is the largest (in magnitude) terms, for matrices this is maximum abolute row sum
          max_norm=400,      // for vectors this is equivalent to Linf, for matrices it is the largest element in matrix
          Frobenius_norm=500 // for vectors this is L2, for matrices this is the sqrt of sum of squares of each term
        };

      private:
        system_norm_type norm_type;

      protected:
        data__ calculate_norm(const solution_matrix &mat) const
        {
          // handle the cases that do not need to distinguish between vectors and matrices
          switch (get_norm_type())
          {
            case(L1):            // for vectors this is sum of maxes, for matrices this is maximum absolute column sum
            {
              data__ rtn_val(-1), tmp;

              for (typename solution_matrix::Index nc=0; nc<mat.cols(); ++nc)
              {
                tmp=mat.col(nc).cwiseAbs().sum();
                if (tmp>rtn_val)
                  rtn_val=tmp;
              }

              return rtn_val;
              break;
            }
            case(Linf):          // for vectors this is the largest (in magnitude) terms, for matrices this is maximum absolute row sum
            {
              data__ rtn_val(-1), tmp;

              for (typename solution_matrix::Index nr=0; nr<mat.rows(); ++nr)
              {
                tmp=mat.row(nr).cwiseAbs().sum();
                if (tmp>rtn_val)
                  rtn_val=tmp;
              }

              return rtn_val;
              break;
            }
            case(max_norm):      // for vectors this is equivalent to Linf, for matrices it is the largest element in matrix
            {
              return mat.maxCoeff();
              break;
            }
            case(L2):             // for vectors this is sqrt of sum of squares, for matrices this will be equivalent to the Frobenius norm
            case(Frobenius_norm): // for vectors this is L2, for matrices this is the sqrt of sum of squares of each term
            {
              data__ rtn_val(0);

              for (typename solution_matrix::Index nc=0; nc<mat.cols(); ++nc)
                rtn_val+=mat.col(nc).squaredNorm();

              return std::sqrt(rtn_val);
              break;
            }
            default:
            {
              // unknown norm type
              assert(false);
              return static_cast<data__>(-1);
            }
          }

          return static_cast<data__>(-1);
        }

      public:
        iterative_system_root_base()
        : iterative_root_base<data__>(), norm_type(iterative_system_root_base<data__, N__, NSOL__>::max_norm)
        {
        }

        iterative_system_root_base(const iterative_system_root_base<data__, N__, NSOL__>& isrb)
        : iterative_root_base<data__>(isrb), norm_type(isrb.norm_type)
        {
        }

        ~iterative_system_root_base()
        {
        }

        system_norm_type get_norm_type() const
        {
          return norm_type;
        }

        void set_norm_type(system_norm_type snt)
        {
          norm_type = snt;
        }
    };
  }
}
#endif
