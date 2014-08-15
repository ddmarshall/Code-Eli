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

#ifndef poly_root_test_suite_hpp
#define poly_root_test_suite_hpp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/mutil/poly/polynomial.hpp"
#include "eli/mutil/poly/root/descartes_rule.hpp"
#include "eli/mutil/poly/root/sturm_count.hpp"
#include "eli/mutil/poly/root/closed_form.hpp"

template<typename data__>
class poly_root_test_suite : public Test::Suite
{
  protected:
    void AddTests(const float &)
    {
      TEST_ADD(poly_root_test_suite<float>::descartes_test);
      TEST_ADD(poly_root_test_suite<float>::sturm_test);
      TEST_ADD(poly_root_test_suite<float>::closed_form_root_test_0);
      TEST_ADD(poly_root_test_suite<float>::closed_form_root_test_1);
      TEST_ADD(poly_root_test_suite<float>::closed_form_root_test_2);
      TEST_ADD(poly_root_test_suite<float>::closed_form_root_test_3);
      TEST_ADD(poly_root_test_suite<float>::closed_form_root_test_4);
      TEST_ADD(poly_root_test_suite<float>::root_finder_test);
    }

    void AddTests(const double &)
    {
      TEST_ADD(poly_root_test_suite<double>::descartes_test);
      TEST_ADD(poly_root_test_suite<double>::sturm_test);
      TEST_ADD(poly_root_test_suite<double>::closed_form_root_test_0);
      TEST_ADD(poly_root_test_suite<double>::closed_form_root_test_1);
      TEST_ADD(poly_root_test_suite<double>::closed_form_root_test_2);
      TEST_ADD(poly_root_test_suite<double>::closed_form_root_test_3);
      TEST_ADD(poly_root_test_suite<double>::closed_form_root_test_4);
      TEST_ADD(poly_root_test_suite<double>::root_finder_test);
    }

    void AddTests(const long double &)
    {
      TEST_ADD(poly_root_test_suite<long double>::descartes_test);
      TEST_ADD(poly_root_test_suite<long double>::sturm_test);
      TEST_ADD(poly_root_test_suite<long double>::closed_form_root_test_0);
      TEST_ADD(poly_root_test_suite<long double>::closed_form_root_test_1);
      TEST_ADD(poly_root_test_suite<long double>::closed_form_root_test_2);
      TEST_ADD(poly_root_test_suite<long double>::closed_form_root_test_3);
      TEST_ADD(poly_root_test_suite<long double>::closed_form_root_test_4);
      TEST_ADD(poly_root_test_suite<long double>::root_finder_test);
    }

  public:
    poly_root_test_suite()
    {
      // add the tests
      AddTests(data__());
    }
    ~poly_root_test_suite()
    {
    }

  private:
    void descartes_test()
    {
      // this is case where know the Sturm functions and the root count in several intervals
      // the roots are between (-4, -3), (-3, -2), (-1, 0), (0, 1), (1, 2)
      {
        typename eli::mutil::poly::polynomial<data__>::coefficient_type coef(6), coef_out, coef_ref;

        // set coefficients
        coef << 2, -10, -20, 0, 5, 1;

        eli::mutil::poly::polynomial<data__> f(coef);

        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f)==5);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, true)==2);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, false)==3);
      }

      // this is case where know the Sturm functions and the real roots, -1 and 1. The
      // other two roots are complex conjugates.
      {
        typename eli::mutil::poly::polynomial<data__>::coefficient_type coef(5), coef_out, coef_ref;

        // set coefficients
        coef << -1, -1, 0, 1, 1;

        eli::mutil::poly::polynomial<data__> f(coef);

        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f)==4);        // should be 2 but extra roots reported on neg. side
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, true)==1);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, false)==3); // should be 1 but extra roots reported
      }

      // this is case with a manually built polynomial with no repeated roots
      {
        data__ roots[7] = {-6, -4, -2, 2, 4, 6, 8};
        eli::mutil::poly::polynomial<data__> f(roots, roots+7);

        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f)==7);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, true)==4);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, false)==3);
      }

      // this is case with a manually built polynomial with repeated root -2 (twice) with
      // -4, 2 and 5 as roots as well
      {
        data__ roots[5] = {-2, -2, -4, 2, 5};
        eli::mutil::poly::polynomial<data__> f(roots, roots+5);

        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f)==5);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, true)==2);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, false)==3);
      }

      // this is case with a manually built polynomial with repeated root -2 (thrice) with
      // -4, 2 and 5 as roots as well
      {
        data__ roots[6] = {-2, -2, -2, -4, 2, 5};
        eli::mutil::poly::polynomial<data__> f(roots, roots+6);

        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f)==6);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, true)==2);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, false)==4);
      }

      // this is case with a manually built polynomial with repeated roots -2 (thrice) and
      // 2 (twice) with -4 and 5 as roots as well
      {
        data__ roots[7] = {-2, -2, -2, 2, 2, -4, 5};
        eli::mutil::poly::polynomial<data__> f(roots, roots+7);

        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f)==7);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, true)==3);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, false)==4);
      }

      // this is case with a manually built polynomial with repeated roots -2 (twice) and
      // 5 (thrice) with -4 and 2 as roots as well
      {
        data__ roots[7] = {-2, -2, 2, 2, 5, -4, 2};
        eli::mutil::poly::polynomial<data__> f(roots, roots+7);

        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f)==7);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, true)==4);
        TEST_ASSERT(eli::mutil::poly::root::descartes_rule(f, false)==3);
      }
    }

    void sturm_test()
    {
      // this is case where know the Sturm functions and the root count in several intervals
      // the roots are between (-4, -3), (-3, -2), (-1, 0), (0, 1), (1, 2)
      {
        typename eli::mutil::poly::polynomial<data__>::coefficient_type coef(6), coef_out, coef_ref;

        // set coefficients
        coef << 2, -10, -20, 0, 5, 1;

        eli::mutil::poly::polynomial<data__> f(coef);
        std::vector< eli::mutil::poly::polynomial<data__> > p;

        // get the Sturm functions and manipulate to match reference values
        eli::mutil::poly::root::sturm_functions(p, f);
        p[3].multiply(3);
        p[4].multiply(17);

        // test p[0]
        coef_ref.resize(6);
        coef_ref << 2, -10, -20, 0, 5, 1;
        p[0].get_coefficients(coef_out);
        TEST_ASSERT(coef_out==coef_ref);

        // test p[1]
        coef_ref.resize(5);
        coef_ref << -2, -8, 0, 4, 1;
        p[1].get_coefficients(coef_out);
        TEST_ASSERT(coef_out==coef_ref);

        // test p[2]
        coef_ref.resize(4);
        coef_ref << -1, 0, 3, 1;
        p[2].get_coefficients(coef_out);
        TEST_ASSERT(coef_out==coef_ref);

        // test p[3]
        coef_ref.resize(3);
        coef_ref << 1, 7, 3;
        p[3].get_coefficients(coef_out);
        if ((typeid(data__)==typeid(float)) || (typeid(data__)==typeid(double)))
        {
          TEST_ASSERT((coef_out-coef_ref).norm()<5*std::numeric_limits<data__>::epsilon());
        }
        else
        {
          TEST_ASSERT(coef_out==coef_ref);
        }

        // test p[4]
        coef_ref.resize(2);
        coef_ref << 11, 17;
        p[4].get_coefficients(coef_out);
        TEST_ASSERT((coef_out-coef_ref).norm()<4986*std::numeric_limits<data__>::epsilon());

        // test p[5]
        TEST_ASSERT(p[5].degree()==0);
        TEST_ASSERT(p[5].f(0)>0);

        // calculate the root count and get back the Sturm functions
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -5,  5)==5); // total real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -5,  0)==3); // negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  0,  5)==2); // positive real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -4, -3)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -3, -2)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -1,  0)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  0,  1)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  1,  2)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -2, -1)==0); // no root

        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p)==5);         // total number of real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, false)==3);  // total number of negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, true)==2);   // total number of positive real roots
      }

      // this is case where know the Sturm functions and the real roots, -1 and 1. The
      // other two roots are complex conjugates.
      {
        typename eli::mutil::poly::polynomial<data__>::coefficient_type coef(5), coef_out, coef_ref;

        // set coefficients
        coef << -1, -1, 0, 1, 1;

        eli::mutil::poly::polynomial<data__> f(coef);
        std::vector< eli::mutil::poly::polynomial<data__> > p;

        // get the Sturm functions and manipulate to match reference values
        eli::mutil::poly::root::sturm_functions(p, f);
        p[1].multiply(4);

        // test p[0]
        coef_ref.resize(5);
        coef_ref << -1, -1, 0, 1, 1;
        p[0].get_coefficients(coef_out);
        TEST_ASSERT(coef_out==coef_ref);

        // test p[1]
        coef_ref.resize(4);
        coef_ref << -1, 0, 3, 4;
        p[1].get_coefficients(coef_out);
        TEST_ASSERT(coef_out==coef_ref);

        // test p[2]
        coef_ref.resize(3);
        coef_ref << 5, 4, 1;
        p[2].get_coefficients(coef_out);
        TEST_ASSERT(coef_out==coef_ref);

        // test p[3]
        coef_ref.resize(2);
        coef_ref << -2, -1;
        p[3].get_coefficients(coef_out);
        TEST_ASSERT(coef_out==coef_ref);

        // test p[4]
        coef_ref.resize(1);
        coef_ref << -1;
        p[4].get_coefficients(coef_out);
        TEST_ASSERT(coef_out==coef_ref);

        // calculate the root count and get back the Sturm functions
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, -5,  5)==2); // total real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, -5,  0)==1); // negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f,  0,  5)==1); // positive real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, -2,  0)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f,  0,  2)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, -3, -2)==0); // no root

        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f)==2);         // total number of real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, false)==1);  // total number of negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, true)==1);   // total number of positive real roots
      }

      // this is case with a manually built polynomial with no repeated roots
      {
        data__ roots[7] = {-6, -4, -2, 2, 4, 6, 8};
        eli::mutil::poly::polynomial<data__> f(roots, roots+7);
        std::vector< eli::mutil::poly::polynomial<data__> > p;

        // get the Sturm functions and manipulate to match reference values
        eli::mutil::poly::root::sturm_functions(p, f);

        // calculate the root count and get back the Sturm functions
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -9,  9)==7); // total real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -9,  0)==3); // negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  0,  9)==4); // positive real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -7, -5)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -5, -3)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -3, -1)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  1,  3)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  3,  5)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  5,  7)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  7,  9)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -8, -7)==0); // no root

        // switched to passing the polynomial for no performance reason, just testing variety
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f)==7);         // total number of real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, false)==3);  // total number of negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, true)==4);   // total number of positive real roots
      }

      // this is case with a manually built polynomial with repeated root -2 (twice) with
      // -4, 2 and 5 as roots as well
      {
        data__ roots[5] = {-2, -2, -4, 2, 5};
        eli::mutil::poly::polynomial<data__> f(roots, roots+5);
        std::vector< eli::mutil::poly::polynomial<data__> > p;

        // get the Sturm functions and manipulate to match reference values
        eli::mutil::poly::root::sturm_functions(p, f);

        // calculate the root count and get back the Sturm functions
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -9,  9)==4); // total real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -9,  0)==2); // negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  0,  9)==2); // positive real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -5, -3)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -3, -1)==1); // double root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  1,  3)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  4,  6)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  3,  4)==0); // no root

        // switched to passing the polynomial for no performance reason, just testing variety
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f)==4);         // total number of real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, false)==2);  // total number of negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, true)==2);   // total number of positive real roots
      }

      // this is case with a manually built polynomial with repeated root -2 (thrice) with
      // -4, 2 and 5 as roots as well
      {
        data__ roots[6] = {-2, -2, -2, -4, 2, 5};
        eli::mutil::poly::polynomial<data__> f(roots, roots+6);
        std::vector< eli::mutil::poly::polynomial<data__> > p;

        // get the Sturm functions and manipulate to match reference values
        eli::mutil::poly::root::sturm_functions(p, f);

        // calculate the root count and get back the Sturm functions
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -9,  9)==4); // total real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -9,  0)==2); // negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  0,  9)==2); // positive real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -5, -3)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -3, -1)==1); // double root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  1,  3)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  4,  6)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  3,  4)==0); // no root

        // switched to passing the polynomial for no performance reason, just testing variety
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f)==4);         // total number of real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, false)==2);  // total number of negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, true)==2);   // total number of positive real roots
      }

      // this is case with a manually built polynomial with repeated roots -2 (thrice) and
      // 2 (twice) with -4 and 5 as roots as well
      {
        data__ roots[7] = {-2, -2, -2, 2, 2, -4, 5};
        eli::mutil::poly::polynomial<data__> f(roots, roots+7);
        std::vector< eli::mutil::poly::polynomial<data__> > p;

        // get the Sturm functions and manipulate to match reference values
        eli::mutil::poly::root::sturm_functions(p, f);

        // calculate the root count and get back the Sturm functions
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -9,  9)==4); // total real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -9,  0)==2); // negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  0,  9)==2); // positive real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -5, -3)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -3, -1)==1); // triple root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  1,  3)==1); // double root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  4,  6)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  3,  4)==0); // no root

        // switched to passing the polynomial for no performance reason, just testing variety
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f)==4);         // total number of real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, false)==2);  // total number of negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, true)==2);   // total number of positive real roots
      }

      // this is case with a manually built polynomial with repeated roots -2 (twice) and
      // 5 (thrice) with -4 and 2 as roots as well
      {
        data__ roots[7] = {-2, -2, 2, 2, 5, -4, 2};
        eli::mutil::poly::polynomial<data__> f(roots, roots+7);
        std::vector< eli::mutil::poly::polynomial<data__> > p;

        // get the Sturm functions and manipulate to match reference values
        eli::mutil::poly::root::sturm_functions(p, f);

        // calculate the root count and get back the Sturm functions
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -9,  9)==4); // total real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -9,  0)==2); // negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  0,  9)==2); // positive real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -5, -3)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p, -3, -1)==1); // double root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  1,  3)==1); // root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  4,  6)==1); // triple root bracket
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(p,  3,  4)==0); // no root

        // switched to passing the polynomial for no performance reason, just testing variety
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f)==4);         // total number of real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, false)==2);  // total number of negative real roots
        TEST_ASSERT(eli::mutil::poly::root::sturm_count(f, true)==2);   // total number of positive real roots
      }
    }

    // test 0th order polynomials (non-zero and zero)
    void closed_form_root_test_0()
    {
      eli::mutil::poly::polynomial<data__> f;

      typename eli::mutil::poly::polynomial<data__>::coefficient_type a;
      int nroot;
      std::vector<data__> root;

      // zero case
      f=0;
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), 0);
      TEST_ASSERT(nroot==std::numeric_limits<int>::max());
      TEST_ASSERT(root.size()==0);

      // non-zero case
      f=1;
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), 0);
      TEST_ASSERT(nroot==0);
      TEST_ASSERT(root.size()==0);
    }

    // test 1st order polynomials
    void closed_form_root_test_1()
    {
      eli::mutil::poly::polynomial<data__> f;

      typename eli::mutil::poly::polynomial<data__>::coefficient_type a;
      data__ root_val;
      int nroot, deg(1);
      std::vector<data__> root;

      // degenerate polynomial
      a.resize(deg+1);
      a.setZero();
      a[deg-1]=1;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==0);
      TEST_ASSERT(root.size()==0);

      // root at zero case
      root_val=static_cast<data__>(0);
      f.set_roots(root_val);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==1);
      TEST_ASSERT(root.size()==1);
      if (root.size()==1)
      {
        TEST_ASSERT(root[0]==root_val);
      }
      root.resize(0);

      // root away from zero case
      root_val=2;
      f.set_roots(root_val);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==1);
      TEST_ASSERT(root.size()==1);
      if (root.size()==1)
      {
        TEST_ASSERT(root[0]==root_val);
      }
    }

    // test 2nd order polynomials
    // (no real roots, one real double root, two real distinct roots)
    void closed_form_root_test_2()
    {
      eli::mutil::poly::polynomial<data__> f;

      typename eli::mutil::poly::polynomial<data__>::coefficient_type a;
      int nroot, deg(2);
      data__ root_val[2];
      std::vector<data__> root;

      // degenerate polynomial
      a.resize(deg+1);
      a.setZero();
      a[deg-1]=1;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==1);
      TEST_ASSERT(root.size()==1);
      if (root.size()==1)
      {
        TEST_ASSERT(root[0]==0);
      }
      root.resize(0);

      // two real roots case
      root_val[0]=-1;
      root_val[1]=2;
      f.set_roots(root_val, root_val+deg);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==2);
      TEST_ASSERT(root.size()==2);
      if (root.size()==2)
      {
        TEST_ASSERT(root[0]==root_val[0]);
        TEST_ASSERT(root[1]==root_val[1]);
      }
      root.resize(0);

      // double root case
      root_val[0]=3;
      root_val[1]=3;
      f.set_roots(root_val, root_val+deg);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==2);
      TEST_ASSERT(root.size()==2);
      if (root.size()==2)
      {
        TEST_ASSERT(root[0]==root_val[0]);
        TEST_ASSERT(root[1]==root_val[1]);
      }
      root.resize(0);

      // no root case
      a.resize(deg+1);
      a << 1, 1, 1;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==0);
      TEST_ASSERT(root.size()==0);
    }

    // test 3rd order polynomials
    // (one real triple root, one unique real & one real double root, three real unique roots, one real root)
    void closed_form_root_test_3()
    {
      eli::mutil::poly::polynomial<data__> f;

      typename eli::mutil::poly::polynomial<data__>::coefficient_type a;
      int nroot, deg(3);
      data__ root_val[3];
      std::vector<data__> root;

      // degenerate polynomial
      a.resize(deg+1);
      a.setZero();
      a[deg-1]=1;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==2);
      TEST_ASSERT(root.size()==2);
      if (root.size()==2)
      {
        TEST_ASSERT(root[0]==0);
        TEST_ASSERT(root[1]==0);
      }
      root.resize(0);

      // caused problems
      a.resize(deg+1);
      a << -54, -27, 0, 1;
      root_val[0]=-3;
      root_val[1]=-3;
      root_val[2]=6;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==3);
      TEST_ASSERT(root.size()==3);
      if (root.size()==3)
      {
        TEST_ASSERT(root[0]==root_val[0]);
        TEST_ASSERT(root[1]==root_val[1]);
        TEST_ASSERT(root[2]==root_val[2]);
      }
      root.resize(0);

      // triple root case
      root_val[0]=2;
      root_val[1]=2;
      root_val[2]=2;
      f.set_roots(root_val, root_val+deg);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==3);
      TEST_ASSERT(root.size()==3);
      if (root.size()==3)
      {
        TEST_ASSERT(root[0]==root_val[0]);
        TEST_ASSERT(root[1]==root_val[1]);
        TEST_ASSERT(root[2]==root_val[2]);
      }
      root.resize(0);

      // one unique & double root case
      root_val[0]=-2;
      root_val[1]=1;
      root_val[2]=1;
      f.set_roots(root_val, root_val+deg);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==3);
      TEST_ASSERT(root.size()==3);
      if (root.size()==3)
      {
        TEST_ASSERT(root[0]==root_val[0]);
        TEST_ASSERT(root[1]==root_val[1]);
        TEST_ASSERT(root[2]==root_val[2]);
      }
      root.resize(0);

      // three unique root case
      root_val[0]=-1;
      root_val[1]=2;
      root_val[2]=1;
      f.set_roots(root_val, root_val+deg);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==3);
      TEST_ASSERT(root.size()==3);
      if (root.size()==3)
      {
        TEST_ASSERT_DELTA(root[0], root_val[0], 2*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(root[1], root_val[2], 2*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(root[2], root_val[1], 2*std::numeric_limits<data__>::epsilon());
      }
      root.resize(0);

      // one real root case
      a.resize(deg+1);
      a << 1, 2, 2, 1;
      root_val[0]=-1;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==1);
      TEST_ASSERT(root.size()==1);
      if (root.size()==1)
      {
        if (typeid(data__)==typeid(double))
        {
          TEST_ASSERT_DELTA(root[0], root_val[0], 2*std::numeric_limits<data__>::epsilon());
        }
        else if (typeid(data__)==typeid(long double))
        {
          TEST_ASSERT_DELTA(root[0], root_val[0], 513*std::numeric_limits<data__>::epsilon());
        }
        else
        {
          TEST_ASSERT_DELTA(root[0], root_val[0], std::numeric_limits<data__>::epsilon());
        }
      }
    }

    // test 4th order polynomials
    // (no real roots, one real quadruple root, one unique real & one triple root, two real double roots, two unique real & one double root, one real double root, two unique real roots)
    void closed_form_root_test_4()
    {
      eli::mutil::poly::polynomial<data__> f;

      typename eli::mutil::poly::polynomial<data__>::coefficient_type a;
      int nroot, deg(4);
      data__ root_val[4];
      std::vector<data__> root;

      // degenerate polynomial
      a.resize(deg+1);
      a.setZero();
      a[deg-1]=1;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==3);
      TEST_ASSERT(root.size()==3);
      TEST_ASSERT(root[0]==0);
      TEST_ASSERT(root[1]==0);
      TEST_ASSERT(root[2]==0);
      root.resize(0);

      // 4 real root case sample calculation
      a.resize(deg+1);
      a << 4, 0, -5, 0, 1;
      root_val[0]=-2;
      root_val[1]=-1;
      root_val[2]=1;
      root_val[3]=2;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==4);
      TEST_ASSERT(root.size()==4);
      if (root.size()==4)
      {
        TEST_ASSERT_DELTA(root[0], root_val[0], 3*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(root[1], root_val[1], 3*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(root[2], root_val[2], 3*std::numeric_limits<data__>::epsilon());
        TEST_ASSERT_DELTA(root[3], root_val[3], 3*std::numeric_limits<data__>::epsilon());
      }
      root.resize(0);

      // 2 real root case sample calculation
      a.resize(deg+1);
      a << -1, 2, 0, 2, 1;
      {
        data__ one(1), two(2);
        root_val[0]=-one-std::sqrt(two);
        root_val[1]=-one+std::sqrt(two);
      }
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==2);
      TEST_ASSERT(root.size()==2);
      if (root.size()==2)
      {
        if (typeid(data__)==typeid(long double))
        {
          TEST_ASSERT_DELTA(root[0], root_val[0], 17*std::numeric_limits<data__>::epsilon());
          TEST_ASSERT_DELTA(root[1], root_val[1], 17*std::numeric_limits<data__>::epsilon());
        }
        else
        {
          TEST_ASSERT(root[0]==root_val[0]);
          TEST_ASSERT(root[1]==root_val[1]);
        }
      }
      root.resize(0);

      // no real root case sample calculation
      a.resize(deg+1);
      a << 4, 0, 5, 0, 1;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==0);
      TEST_ASSERT(root.size()==0);

      // quadruple root case
      root_val[0]=2;
      root_val[1]=2;
      root_val[2]=2;
      root_val[3]=2;
      f.set_roots(root_val, root_val+deg);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==4);
      TEST_ASSERT(root.size()==4);
      if (root.size()==4)
      {
        TEST_ASSERT(root[0]==root_val[0]);
        TEST_ASSERT(root[1]==root_val[1]);
        TEST_ASSERT(root[2]==root_val[2]);
        TEST_ASSERT(root[3]==root_val[3]);
      }
      root.resize(0);

      // one unique & triple root case
      root_val[0]=-2;
      root_val[1]=1;
      root_val[2]=1;
      root_val[3]=1;
      f.set_roots(root_val, root_val+deg);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==4);
      TEST_ASSERT(root.size()==4);
      if (root.size()==4)
      {
        TEST_ASSERT(root[0]==root_val[0]);
        TEST_ASSERT(root[1]==root_val[1]);
        TEST_ASSERT(root[2]==root_val[2]);
        TEST_ASSERT(root[3]==root_val[3]);
      }
      root.resize(0);

      // two real double root case
      root_val[0]=-2;
      root_val[1]=-2;
      root_val[2]=1;
      root_val[3]=1;
      f.set_roots(root_val, root_val+deg);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==4);
      TEST_ASSERT(root.size()==4);
      if (root.size()==4)
      {
        TEST_ASSERT(root[0]==root_val[0]);
        TEST_ASSERT(root[1]==root_val[1]);
        TEST_ASSERT(root[2]==root_val[2]);
        TEST_ASSERT(root[3]==root_val[3]);
      }
      root.resize(0);

      // two unique real & real double root case
      root_val[0]=-2;
      root_val[1]=-2;
      root_val[2]=1;
      root_val[3]=4;
      f.set_roots(root_val, root_val+deg);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==4);
      TEST_ASSERT(root.size()==4);
      if (root.size()==4)
      {
        if ((typeid(data__)==typeid(float)) || (typeid(data__)==typeid(double)) || (typeid(data__)==typeid(long double)))
        {
          TEST_ASSERT(root[0]==root_val[0]);
          TEST_ASSERT(root[1]==root_val[1]);
          TEST_ASSERT_DELTA(root[2], root_val[2], std::numeric_limits<data__>::epsilon());
          TEST_ASSERT(root[3]==root_val[3]);
        }
        else
        {
          TEST_ASSERT(root[0]==root_val[0]);
          TEST_ASSERT(root[1]==root_val[1]);
          TEST_ASSERT(root[2]==root_val[2]);
          TEST_ASSERT(root[3]==root_val[3]);
        }
      }
      root.resize(0);

      // four unique root case
      root_val[0]=-1;
      root_val[1]=2;
      root_val[2]=1;
      root_val[3]=4;
      f.set_roots(root_val, root_val+deg);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==4);
      TEST_ASSERT(root.size()==4);
      if (root.size()==4)
      {
        if ((typeid(data__)==typeid(float)) || (typeid(data__)==typeid(long double)))
        {
          TEST_ASSERT_DELTA(root[0], root_val[0], std::numeric_limits<data__>::epsilon());
          TEST_ASSERT_DELTA(root[1], root_val[2], std::numeric_limits<data__>::epsilon());
          TEST_ASSERT(root[2]==root_val[1]);
          TEST_ASSERT(root[3]==root_val[3]);
        }
        else
        {
          TEST_ASSERT(root[0]==root_val[0]);
          TEST_ASSERT(root[1]==root_val[2]);
          TEST_ASSERT(root[2]==root_val[1]);
          TEST_ASSERT(root[3]==root_val[3]);
        }
      }
      root.resize(0);

      // one real double root
      root_val[0]=1;
      root_val[1]=1;
      a.resize(deg+1);
      a << root_val[0]*root_val[1],
           root_val[0]*root_val[1]-root_val[0]-root_val[1],
           root_val[0]*root_val[1]-root_val[0]-root_val[1]+1,
           1-root_val[0]-root_val[1],
           1;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==2);
      TEST_ASSERT(root.size()==2);
      if (root.size()==2)
      {
        TEST_ASSERT(root[0]==root_val[0]);
        TEST_ASSERT(root[1]==root_val[1]);
      }
      root.resize(0);

      // two real root
      root_val[0]=2;
      root_val[1]=1;
      a.resize(deg+1);
      a << root_val[0]*root_val[1],
           root_val[0]*root_val[1]-root_val[0]-root_val[1],
           root_val[0]*root_val[1]-root_val[0]-root_val[1]+1,
           1-root_val[0]-root_val[1],
           1;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==2);
      TEST_ASSERT(root.size()==2);
      if (root.size()==2)
      {
        TEST_ASSERT(root[0]==root_val[1]);
        TEST_ASSERT(root[1]==root_val[0]);
      }
      root.resize(0);

      // no real root
      a.resize(deg+1);
      a << 2, 4, 5, 3, 1;
      f.set_coefficients(a);
      f.get_coefficients(a);
      nroot=eli::mutil::poly::root::closed_form(root, a.data(), deg);
      TEST_ASSERT(nroot==0);
      TEST_ASSERT(root.size()==0);
    }

    void root_finder_test()
    {
      // -1, 1, -1, 1, -1, 1 has one real root which is 1
    }
};

#endif
