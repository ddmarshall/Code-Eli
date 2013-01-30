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

#ifndef nls_test_suite_hpp
#define nls_test_suite_hpp

#include "Eigen/Eigen"

#include "eli/code_eli.hpp"

#include "eli/nls.hpp"
#include "eli/constants/math.hpp"

#include <cmath>      // cos()
#include <functional> // std::ptr_fun()
#include <limits>     // numeric_limits

template<typename data__>
data__ my_function(const data__ &x)
{
  return std::cos(x);
}

template<typename data__>
data__ my_function_derivative(const data__ &x)
{
  return -std::sin(x);
}

template<typename data__>
struct my_functor
{
    data__ operator()(const data__ &x) const
    {
      return my_function(x);
    }
};

template<typename data__>
struct my_functor_derivative
{
    data__ operator()(const data__ &x) const
    {
      return my_function_derivative(x);
    }
};

template<typename data__>
Eigen::Matrix<data__, 3, 1> my_decoupled_system_function(const Eigen::Matrix<data__, 3, 1> &x)
{
  Eigen::Matrix<data__, 3, 1> rtn_vec;

  rtn_vec(0,0)=std::exp(x(0));
  rtn_vec(1,0)=std::exp(x(1));
  rtn_vec(2,0)=std::exp(x(2));

  return rtn_vec;
}

template<typename data__>
Eigen::Matrix<data__, 3, 3> my_decoupled_system_function_derivative(const Eigen::Matrix<data__, 3, 1> &x)
{
  Eigen::Matrix<data__, 3, 3> rtn_mat;

  rtn_mat(0,0)=std::exp(x(0));
  rtn_mat(0,1)=0.0;
  rtn_mat(0,2)=0.0;
  rtn_mat(1,0)=0.0;
  rtn_mat(1,1)=std::exp(x(1));
  rtn_mat(1,2)=0.0;
  rtn_mat(2,0)=0.0;
  rtn_mat(2,1)=0.0;
  rtn_mat(2,2)=std::exp(x(2));

  return rtn_mat;
}

template<typename data__>
struct my_decoupled_system_functor
{
    Eigen::Matrix<data__, 3, 1> operator()(const Eigen::Matrix<data__, 3, 1> &x) const
    {
      return my_decoupled_system_function(x);
    }
};

template<typename data__>
struct my_decoupled_system_functor_derivative
{
    Eigen::Matrix<data__, 3, 3> operator()(const Eigen::Matrix<data__, 3, 1> &x) const
    {
      return my_decoupled_system_function_derivative(x);
    }
};

template<typename data__>
Eigen::Matrix<data__, 3, 1> my_coupled_linear_system_function(const Eigen::Matrix<data__, 3, 1> &x)
{
  Eigen::Matrix<data__, 3, 1> rtn_vec;

  rtn_vec(0,0)= 3*x(0)+2*x(1)+4*x(2);
  rtn_vec(1,0)= 1*x(0)-4*x(1)+2*x(2);
  rtn_vec(2,0)=-2*x(0)+3*x(1)-1*x(2);

  return rtn_vec;
}

template<typename data__>
Eigen::Matrix<data__, 3, 3> my_coupled_linear_system_function_derivative(const Eigen::Matrix<data__, 3, 1> &)
{
  Eigen::Matrix<data__, 3, 3> rtn_mat;

  rtn_mat(0,0)= 3.0;
  rtn_mat(0,1)= 2.0;
  rtn_mat(0,2)= 4.0;
  rtn_mat(1,0)= 1.0;
  rtn_mat(1,1)=-4.0;
  rtn_mat(1,2)= 2.0;
  rtn_mat(2,0)=-2.0;
  rtn_mat(2,1)= 3.0;
  rtn_mat(2,2)=-1.0;

  return rtn_mat;
}

template<typename data__>
struct my_coupled_linear_system_functor
{
    Eigen::Matrix<data__, 3, 1> operator()(const Eigen::Matrix<data__, 3, 1> &x) const
    {
      return my_coupled_linear_system_function(x);
    }
};

template<typename data__>
struct my_coupled_linear_system_functor_derivative
{
    Eigen::Matrix<data__, 3, 3> operator()(const Eigen::Matrix<data__, 3, 1> &x) const
    {
      return my_coupled_linear_system_function_derivative(x);
    }
};

template<typename data__>
Eigen::Matrix<data__, 3, 1> my_coupled_nonlinear_system_function(const Eigen::Matrix<data__, 3, 1> &x)
{
  Eigen::Matrix<data__, 3, 1> rtn_vec;

  rtn_vec(0,0)=std::exp(x(0))*cos(x(1));
  rtn_vec(1,0)=x(0)*x(0)*std::log(x(1));
  rtn_vec(2,0)=x(2);

  return rtn_vec;
}

template<typename data__>
Eigen::Matrix<data__, 3, 3> my_coupled_nonlinear_system_function_derivative(const Eigen::Matrix<data__, 3, 1> &x)
{
  Eigen::Matrix<data__, 3, 3> rtn_mat;

  rtn_mat(0,0)=std::exp(x(0))*std::cos(x(1));
  rtn_mat(0,1)=-std::exp(x(0))*std::sin(x(1));
  rtn_mat(0,2)=0.0;
  rtn_mat(1,0)=2*x(0)*std::log(x(1));
  rtn_mat(1,1)=x(0)*x(0)/x(1);
  rtn_mat(1,2)=0.0;
  rtn_mat(2,0)=0.0;
  rtn_mat(2,1)=0.0;
  rtn_mat(2,2)=1.0;

  return rtn_mat;
}

template<typename data__>
struct my_coupled_nonlinear_system_functor
{
    Eigen::Matrix<data__, 3, 1> operator()(const Eigen::Matrix<data__, 3, 1> &x) const
    {
      return my_coupled_nonlinear_system_function(x);
    }
};

template<typename data__>
struct my_coupled_nonlinear_system_functor_derivative
{
    Eigen::Matrix<data__, 3, 3> operator()(const Eigen::Matrix<data__, 3, 1> &x) const
    {
      return my_coupled_nonlinear_system_function_derivative(x);
    }
};

template<typename data__>
class nls_test_suite : public Test::Suite
{
  protected:
    void AddTests(const float &)
    {
      TEST_ADD(nls_test_suite<float>::bisection_method_test);
      TEST_ADD(nls_test_suite<float>::newton_raphson_method_test);
      TEST_ADD(nls_test_suite<float>::secant_method_test);
      TEST_ADD(nls_test_suite<float>::newton_raphson_system_method_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(nls_test_suite<double>::bisection_method_test);
      TEST_ADD(nls_test_suite<double>::newton_raphson_method_test);
      TEST_ADD(nls_test_suite<double>::secant_method_test);
      TEST_ADD(nls_test_suite<double>::newton_raphson_system_method_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(nls_test_suite<long double>::bisection_method_test);
      TEST_ADD(nls_test_suite<long double>::newton_raphson_method_test);
      TEST_ADD(nls_test_suite<long double>::secant_method_test);
      TEST_ADD(nls_test_suite<long double>::newton_raphson_system_method_test);
    }
#ifdef ELI_QD_FOUND
    void AddTests(const dd_real &)
    {
      TEST_ADD(nls_test_suite<dd_real>::bisection_method_test);
      TEST_ADD(nls_test_suite<dd_real>::newton_raphson_method_test);
      TEST_ADD(nls_test_suite<dd_real>::secant_method_test);
      TEST_ADD(nls_test_suite<dd_real>::newton_raphson_system_method_test);
    }

    void AddTests(const qd_real &)
    {
      TEST_ADD(nls_test_suite<qd_real>::bisection_method_test);
      TEST_ADD(nls_test_suite<qd_real>::newton_raphson_method_test);
      TEST_ADD(nls_test_suite<qd_real>::secant_method_test);
      TEST_ADD(nls_test_suite<qd_real>::newton_raphson_system_method_test);
    }
#endif

  public:
    nls_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~nls_test_suite()
    {
    }

  private:
    void bisection_method_test()
    {
      data__ delta(std::sqrt(std::numeric_limits<data__>::epsilon())), rhs(0.5), root;
      eli::nls::bisection_method<data__> bm;
      typename eli::nls::bisection_method<data__>::status stat;

      bm.set_absolute_tolerance(delta);
      bm.set_max_iteration(200);
      bm.set_bounds(0, eli::constants::math<data__>::pi());

      // test using user defined function
      bm.set_absolute_tolerance(delta);
      bm.set_max_iteration(200);
      bm.set_bounds(0, eli::constants::math<data__>::pi());

      stat = bm.find_root(root, std::ptr_fun(my_function<data__>), rhs);
      TEST_ASSERT(stat==eli::nls::bisection_method<data__>::converged);
      TEST_ASSERT_DELTA(root, std::acos(rhs), 2*delta);

      // test using functor
      bm.set_absolute_tolerance(delta);
      bm.set_max_iteration(200);
      bm.set_bounds(0, eli::constants::math<data__>::pi());

      stat = bm.find_root(root, my_functor<data__>(), rhs);
      TEST_ASSERT(stat==eli::nls::bisection_method<data__>::converged);
      TEST_ASSERT_DELTA(root, std::acos(rhs), 2*delta);
    }

    void newton_raphson_method_test()
    {
      data__ delta(std::sqrt(std::numeric_limits<data__>::epsilon())), rhs(0.5), root;
      eli::nls::newton_raphson_method<data__> nrm;
      typename eli::nls::newton_raphson_method<data__>::status stat;

      nrm.set_absolute_tolerance(delta);
      nrm.set_max_iteration(200);
      nrm.set_initial_guess(0.3*eli::constants::math<data__>::pi_by_four());

      // test using user defined functions
      stat = nrm.find_root(root, std::ptr_fun(my_function<data__>), std::ptr_fun(my_function_derivative<data__>), rhs);
      TEST_ASSERT(stat==eli::nls::newton_raphson_method<data__>::converged);
      TEST_ASSERT_DELTA(root, std::acos(rhs), 2*delta);

      nrm.set_max_iteration(2);
      stat = nrm.find_root(root, std::ptr_fun(my_function<data__>), std::ptr_fun(my_function_derivative<data__>), rhs);
      TEST_ASSERT(stat==eli::nls::newton_raphson_method<data__>::max_iteration);

      nrm.set_max_iteration(200);
      nrm.set_initial_guess(0);
      stat = nrm.find_root(root, std::ptr_fun(my_function<data__>), std::ptr_fun(my_function_derivative<data__>), rhs);
      TEST_ASSERT(stat==eli::nls::newton_raphson_method<data__>::no_root_found);

      // test using functor
      nrm.set_absolute_tolerance(delta);
      nrm.set_max_iteration(200);
      nrm.set_initial_guess(eli::constants::math<data__>::pi_by_four());

      stat = nrm.find_root(root, my_functor<data__>(), my_functor_derivative<data__>(), rhs);
      TEST_ASSERT(stat==eli::nls::newton_raphson_method<data__>::converged);
      TEST_ASSERT_DELTA(root, std::acos(rhs), 2*delta);
    }

    void secant_method_test()
    {
      data__ delta(std::sqrt(std::numeric_limits<data__>::epsilon())), rhs(0.5), root(0);
      eli::nls::secant_method<data__> sm;
      typename eli::nls::secant_method<data__>::status stat;

      sm.set_absolute_tolerance(delta);
      sm.set_max_iteration(200);
      sm.set_initial_guesses(eli::constants::math<data__>::pi_by_four(), eli::constants::math<data__>::pi_by_two());

      // test using user defined function
      sm.set_absolute_tolerance(delta);
      sm.set_max_iteration(200);
      sm.set_initial_guesses(eli::constants::math<data__>::pi_by_four(), eli::constants::math<data__>::pi_by_two());

      stat = sm.find_root(root, std::ptr_fun(my_function<data__>), rhs);
      TEST_ASSERT(stat==eli::nls::secant_method<data__>::converged);
      TEST_ASSERT_DELTA(root, std::acos(rhs), 2*delta);

      // test using functor
      sm.set_absolute_tolerance(delta);
      sm.set_max_iteration(200);
      sm.set_initial_guesses(eli::constants::math<data__>::pi_by_four(), eli::constants::math<data__>::pi_by_two());

      stat = sm.find_root(root, my_functor<data__>(), rhs);
      TEST_ASSERT(stat==eli::nls::secant_method<data__>::converged);
      TEST_ASSERT_DELTA(root, std::acos(rhs), 2*delta);
    }

    void newton_raphson_system_method_test()
    {
      typedef eli::nls::newton_raphson_system_method<data__, 3, 1> nr_system;

      data__ delta(std::sqrt(std::numeric_limits<data__>::epsilon()));
      nr_system nrm;
      typename nr_system::solution_matrix rhs, root, x0, x_exact;
      typename nr_system::status stat;

      nrm.set_absolute_tolerance(delta);
      nrm.set_max_iteration(200);
      nrm.set_norm_type(nr_system::max_norm);

      // decoupled system
      // set right hand side, initial guess & exact answer
      x_exact << 2, 3, 1;
      rhs = my_decoupled_system_function<data__>(x_exact);
      x0 << 1.5, 2.5, 1.5;
      nrm.set_initial_guess(x0);

      // test using user defined functions
      stat = nrm.find_root(root, std::ptr_fun(my_decoupled_system_function<data__>), std::ptr_fun(my_decoupled_system_function_derivative<data__>), rhs);
      TEST_ASSERT(stat==nr_system::converged);
      TEST_ASSERT(nrm.get_iteration_count()<nrm.get_max_iteration());
      TEST_ASSERT((x_exact-root).norm()<=2*delta);

      nrm.set_max_iteration(2);
      stat = nrm.find_root(root, std::ptr_fun(my_decoupled_system_function<data__>), std::ptr_fun(my_decoupled_system_function_derivative<data__>), rhs);
      TEST_ASSERT(stat==nr_system::max_iteration);

      // test using functor
      nrm.set_absolute_tolerance(delta);
      nrm.set_max_iteration(200);
      nrm.set_initial_guess(x0);

      stat = nrm.find_root(root, my_decoupled_system_functor<data__>(), my_decoupled_system_functor_derivative<data__>(), rhs);
      TEST_ASSERT(stat==nr_system::converged);
      TEST_ASSERT(nrm.get_iteration_count()<nrm.get_max_iteration());
      TEST_ASSERT((x_exact-root).norm()<=2*delta);

      // linear coupled system
      // set right hand side, initial guess & exact answer
      x_exact << 2, 3, 1;
      rhs = my_coupled_linear_system_function<data__>(x_exact);
      x0 << 1.5, 2.5, 1.5;
      nrm.set_initial_guess(x0);

      // test using user defined functions
      stat = nrm.find_root(root, std::ptr_fun(my_coupled_linear_system_function<data__>), std::ptr_fun(my_coupled_linear_system_function_derivative<data__>), rhs);
      TEST_ASSERT(stat==nr_system::converged);
      TEST_ASSERT(nrm.get_iteration_count()<2);
      TEST_ASSERT((x_exact-root).norm()<delta);

      // test using functor
      nrm.set_absolute_tolerance(delta);
      nrm.set_max_iteration(200);
      nrm.set_initial_guess(x0);

      stat = nrm.find_root(root, my_coupled_linear_system_functor<data__>(), my_coupled_linear_system_functor_derivative<data__>(), rhs);
      TEST_ASSERT(stat==nr_system::converged);
      TEST_ASSERT(nrm.get_iteration_count()<2);
      TEST_ASSERT((x_exact-root).norm()<delta);

      // nonlinear coupled system
      // set right hand side, initial guess & exact answer
      x_exact << 2, 3, 1;
      rhs = my_coupled_nonlinear_system_function<data__>(x_exact);
      x0 << 1.5, 2.5, 1.5;
      nrm.set_initial_guess(x0);

      // test using user defined functions
      stat = nrm.find_root(root, std::ptr_fun(my_coupled_nonlinear_system_function<data__>), std::ptr_fun(my_coupled_nonlinear_system_function_derivative<data__>), rhs);
      TEST_ASSERT(stat==nr_system::converged);
      TEST_ASSERT(nrm.get_iteration_count()<nrm.get_max_iteration());
      TEST_ASSERT((x_exact-root).norm()<delta);

      // test using functor
      nrm.set_absolute_tolerance(delta);
      nrm.set_max_iteration(200);
      nrm.set_initial_guess(x0);

      stat = nrm.find_root(root, my_coupled_nonlinear_system_functor<data__>(), my_coupled_nonlinear_system_functor_derivative<data__>(), rhs);
      TEST_ASSERT(stat==nr_system::converged);
      TEST_ASSERT(nrm.get_iteration_count()<nrm.get_max_iteration());
      TEST_ASSERT((x_exact-root).norm()<delta);
    }
};

#endif
