/*********************************************************************************
* Copyright (c) 2013-2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_curve_piecewise_four_digit_creator_hpp
#define eli_geom_curve_piecewise_four_digit_creator_hpp

#include "eli/code_eli.hpp"

#include "eli/mutil/dm/binomial_coefficient.hpp"
#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/piecewise_polynomial_creator.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/pseudo/polynomial.hpp"
#include "eli/geom/curve/pseudo/four_digit.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_four_digit_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef eli::geom::curve::pseudo::four_digit<data_type> airfoil_type;

        private:
          typedef eli::mutil::poly::polynomial<data_type> polynomial_type;

        public:
          piecewise_four_digit_creator() : piecewise_creator_base<data_type, dim__, tolerance_type>(4, 0), trig_degree(6) {}
          piecewise_four_digit_creator(const piecewise_four_digit_creator<data_type, dim__, tolerance_type> &ppc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(ppc), af(ppc.af), trig_degree(ppc.trig_degree) {}
          ~piecewise_four_digit_creator() {}

          void set_sharp_trailing_edge(bool fl)
          {
            af.set_sharp_trailing_edge(fl);
          }
          bool sharp_trailing_edge() const
          {
            return af.sharp_trailing_edge();
          }

          bool set_thickness(const data_type &t)
          {
            return af.set_thickness(t);
          }

          data_type get_thickness() const
          {
            return af.get_thickness();
          }

          bool set_camber(const data_type &cam, const data_type &cam_loc)
          {
            return af.set_camber(cam, cam_loc);
          }

          data_type get_maximum_camber() const
          {
            return af.get_maximum_camber();
          }

          data_type get_maximum_camber_location() const
          {
            return af.get_maximum_camber_location();
          }

          bool set_name(const std::string &name)
          {
            return af.set_name(name);
          }

          std::string get_name() const
          {
            return af.get_name();
          }

          void set_trig_terms(index_type nt)
          {
            if ( (nt>2) && (nt<11) )
            {
              trig_degree=nt;
            }
          }

          index_type get_trig_terms() const
          {
            return trig_degree;
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;
            typedef typename airfoil_type::coefficient_type coefficient_type;

            // get the airfoil information
            data_type t, m, p;
            coefficient_type a(af.get_thickness_coefficients());
            bool symmetric;

            t=get_thickness()/100;
            m=get_maximum_camber()/100;
            p=get_maximum_camber_location()/10;
            symmetric = ((m<=0) || (p<=0));

            // build thickness polynomial coefficients
            polynomial_type thickness;
            {
              typename polynomial_type::coefficient_type delta(9, 1);

              delta << 0, a(0), a(1), 0, a(2), 0, a(3), 0, a(4);
              delta *= (t/0.20);
              thickness.set_coefficients(delta);
            }

            // build the camber terms
            polynomial_type camber_x, camber_front_y, camber_back_y;
            data_type t_split(std::sqrt(p));
            {
              typename polynomial_type::coefficient_type x(3, 1), y(5, 1);

              x << 0, 0, 1;
              camber_x.set_coefficients(x);

              if (symmetric)
              {
                y.setZero();
              }
              else
              {
                y << 0, 0, 2*m/p, 0, -m/p/p;
              }
              camber_front_y.set_coefficients(y);

              if (symmetric)
              {
                y.setZero();
              }
              else
              {
                y << m*(1-2*p), 0, 2*m*p, 0, -m;
                y /= (1-p)*(1-p);
              }
              camber_back_y.set_coefficients(y);
            }

            polynomial_type camber_front_cos, camber_front_sin, camber_back_cos, camber_back_sin;
            {
              if (symmetric)
              {
                typename polynomial_type::coefficient_type cterm(1,1), sterm(1,1);

                cterm(0)=1;
                sterm(0)=0;

                camber_front_cos.set_coefficients(cterm);
                camber_front_sin.set_coefficients(sterm);
                camber_back_cos.set_coefficients(cterm);
                camber_back_sin.set_coefficients(sterm);
              }
              else
              {
                typename polynomial_type::coefficient_type cterm, sterm;
                data_type a;

                // build the front terms
                a=2*m/p/p;
                build_cos_term(cterm, a, p, 0, trig_degree);
                build_sin_term(sterm, a, p, cterm);
                camber_front_cos.set_coefficients(cterm);
                camber_front_sin.set_coefficients(sterm);
//
//                std::cout << "camber_front_cos=" << cterm << std::endl;
//                std::cout << "camber_front_sin=" << sterm << std::endl;

                // build the back terms
                a=2*m/(1-p)/(1-p);
                build_cos_term(cterm, a, p, 1, trig_degree);
                build_sin_term(sterm, a, p, cterm);
                camber_back_cos.set_coefficients(cterm);
                camber_back_sin.set_coefficients(sterm);
//
//                std::cout << "camber_back_cos=" << cterm << std::endl;
//                std::cout << "camber_back_sin=" << sterm << std::endl;
              }
            }

            // combine the thickness and camber terms
            pseudo::polynomial<data_type, dim__> cu_front, cl_front, cu_back, cl_back;
            {
              polynomial_type front_upper_x, front_upper_y, front_lower_x, front_lower_y;
              polynomial_type back_upper_x, back_upper_y, back_lower_x, back_lower_y;
              polynomial_type delta_cos, delta_sin;

              // build front curves
              delta_cos.multiply(thickness, camber_front_cos);
              delta_sin.multiply(thickness, camber_front_sin);
              front_upper_x.subtract(camber_x, delta_sin);
              front_upper_y.add(camber_front_y, delta_cos);
              front_lower_x.add(camber_x, delta_sin);
              front_lower_y.subtract(camber_front_y, delta_cos);
//              if (!symmetric)
//              {
//                typename polynomial_type::coefficient_type co;
//                camber_front_y.get_coefficients(co);
//                std::cout << "camber_front_y=" << std::endl << co << std::endl;
//                delta_cos.get_coefficients(co);
//                std::cout << "delta_cos=" << std::endl << co << std::endl;
//                front_lower_x.get_coefficients(co);
//                std::cout << "front_lower_x=" << std::endl << co << std::endl;
//                front_lower_y.get_coefficients(co);
//                std::cout << "front_lower_y=" << std::endl << co << std::endl;
//              }

              // build back curves
              delta_cos.multiply(thickness, camber_back_cos);
              delta_sin.multiply(thickness, camber_back_sin);
              back_upper_x.subtract(camber_x, delta_sin);
              back_upper_y.add(camber_back_y, delta_cos);
              back_lower_x.add(camber_x, delta_sin);
              back_lower_y.subtract(camber_back_y, delta_cos);
//              if (!symmetric)
//              {
//                typename polynomial_type::coefficient_type co;
//                back_upper_x.get_coefficients(co);
//                std::cout << "back_lower_x=" << std::endl << co << std::endl;
//                back_upper_y.get_coefficients(co);
//                std::cout << "back_lower_y=" << std::endl << co << std::endl;
//              }

              // set the polynomial pseudo-curve coefficients
              typename polynomial_type::coefficient_type co;
              front_upper_x.get_coefficients(co); cu_front.set_coefficients(co, 0);
              front_upper_y.get_coefficients(co); cu_front.set_coefficients(co, 1);
              back_upper_x.get_coefficients(co);  cu_back.set_coefficients(co, 0);
              back_upper_y.get_coefficients(co);  cu_back.set_coefficients(co, 1);
              front_lower_x.get_coefficients(co); cl_front.set_coefficients(co, 0);
              front_lower_y.get_coefficients(co); cl_front.set_coefficients(co, 1);
              back_lower_x.get_coefficients(co);  cl_back.set_coefficients(co, 0);
              back_lower_y.get_coefficients(co);  cl_back.set_coefficients(co, 1);
//              if (!symmetric && (typeid(data_type)==typeid(float)))
//              {
//                pseudo::polynomial<data_type, dim__> temp_front, temp_back;
//
//                camber_x.get_coefficients(co);
//                temp_front.set_coefficients(co, 0);
//                camber_front_y.get_coefficients(co);
//                temp_front.set_coefficients(co, 1);
//                camber_x.get_coefficients(co);
//                temp_back.set_coefficients(co, 0);
//                camber_back_y.get_coefficients(co);
//                temp_back.set_coefficients(co, 1);
//                std::cout.flush();
//                eli::test::octave_start(1);
//                eli::test::octave_print(1, temp_front, "c_front");
//                eli::test::octave_print(1, temp_back, "c_back");
//                eli::test::octave_print(1, cu_front, "u_front");
//                eli::test::octave_print(1, cu_back, "u_back");
//                eli::test::octave_print(1, cl_front, "l_front");
//                eli::test::octave_print(1, cl_back, "l_back");
//                eli::test::octave_finish(1);
//              }
            }

            // create airfoil
            piecewise_polynomial_creator<data_type, dim__, tolerance_type> poly_creator;
            piecewise_curve_type pc_temp;
            typename curve_type::control_point_type cp, cp_split;
            curve_type c;
            bool rtn_flag;
            error_code er;

            //
            // lower surface
            //

            // build lower aft curve
            rtn_flag=poly_creator.set_conditions(cl_back);
            if (!rtn_flag)
            {
              assert(false);
              return false;
            }
            poly_creator.set_t0(0);
            poly_creator.set_segment_dt(1, 0);
            rtn_flag=poly_creator.create(pc_temp);
            if (!rtn_flag)
            {
              assert(false);
              return false;
            }
            pc_temp.reverse();
            // if symmetric airfoil, then only one curve for lower surface
            if (symmetric)
            {
              pc_temp.get(c, 0);
            }
            // else need to split curve and extract the aft portion for airfoil
            else
            {
              er=pc_temp.split(1-t_split);
              if (er!=piecewise_curve_type::NO_ERRORS)
              {
                pc.clear();
                return false;
              }
              pc_temp.get(c, 0);
              cp_split=c.get_control_point(c.degree());
            }
            if (sharp_trailing_edge())
            {
              cp.setZero();
              cp(0)=1;
              c.set_control_point(cp, 0);
            }
            er=pc.push_back(c, 1-t_split);
            if (er!=piecewise_curve_type::NO_ERRORS)
            {
              pc.clear();
              return false;
            }
//            if (!symmetric && (typeid(data_type)==typeid(float)))
//            {
//              pseudo::polynomial<data_type, dim__> temp_front, temp_back;
//
//              std::cout.flush();
//              eli::test::octave_start(1);
//              eli::test::octave_print(1, c, "lower");
//              eli::test::octave_print(1, pc, "lower");
//              eli::test::octave_print(1, cl_back, "l_back");
//              eli::test::octave_finish(1);
//            }

            // if not symmetric airfoil then need to get the front portion of curve
            if (!symmetric)
            {
              rtn_flag=poly_creator.set_conditions(cl_front);
              if (!rtn_flag)
              {
                assert(false);
                return false;
              }
              poly_creator.set_t0(0);
              poly_creator.set_segment_dt(1, 0);
              rtn_flag=poly_creator.create(pc_temp);
              if (!rtn_flag)
              {
                assert(false);
                return false;
              }
              pc_temp.reverse();
              er=pc_temp.split(1-t_split);
              if (er!=piecewise_curve_type::NO_ERRORS)
              {
                pc.clear();
                return false;
              }
              pc_temp.get(c, 1);
              c.set_control_point(cp_split, 0);
              er=pc.push_back(c, t_split);
              if (er!=piecewise_curve_type::NO_ERRORS)
              {
                pc.clear();
                return false;
              }
            }
//            if (!symmetric && (typeid(data_type)==typeid(float)))
//            {
//              pseudo::polynomial<data_type, dim__> temp_front, temp_back;
//
//              std::cout.flush();
//              eli::test::octave_start(1);
//              eli::test::octave_print(1, pc, "lower");
//              eli::test::octave_print(1, cl_front, "l_front");
//              eli::test::octave_print(1, cl_back, "l_back");
//              eli::test::octave_finish(1);
//            }

            //
            // upper surface
            //

            // build upper front curve
            rtn_flag=poly_creator.set_conditions(cu_front);
            if (!rtn_flag)
            {
              assert(false);
              return false;
            }
            poly_creator.set_t0(0);
            poly_creator.set_segment_dt(1, 0);
            rtn_flag=poly_creator.create(pc_temp);
            if (!rtn_flag)
            {
              assert(false);
              return false;
            }
            // if symmetric airfoil, then only one curve for upper surface
            if (symmetric)
            {
              pc_temp.get(c, 0);
              if (sharp_trailing_edge())
              {
                cp.setZero();
                cp(0)=1;
                c.set_control_point(cp, c.degree());
              }
              er=pc.push_back(c, 1);
            }
            // else need to split curve and extract the front portion for airfoil
            else
            {
              er=pc_temp.split(t_split);
              if (er!=piecewise_curve_type::NO_ERRORS)
              {
                pc.clear();
                return false;
              }
              pc_temp.get(c, 0);
              cp_split=c.get_control_point(c.degree());
              er=pc.push_back(c, t_split);
            }
            if (er!=piecewise_curve_type::NO_ERRORS)
            {
              pc.clear();
              return false;
            }

            // if not symmetric airfoil then need to get the aft portion of curve
            if (!symmetric)
            {
              rtn_flag=poly_creator.set_conditions(cu_back);
              if (!rtn_flag)
              {
                assert(false);
                return false;
              }
              poly_creator.set_t0(0);
              poly_creator.set_segment_dt(1, 0);
              rtn_flag=poly_creator.create(pc_temp);
              if (!rtn_flag)
              {
                assert(false);
                return false;
              }
              er=pc_temp.split(t_split);
              if (er!=piecewise_curve_type::NO_ERRORS)
              {
                pc.clear();
                return false;
              }
              pc_temp.get(c, 1);
              c.set_control_point(cp_split, 0);
              if (sharp_trailing_edge())
              {
                cp.setZero();
                cp(0)=1;
                c.set_control_point(cp, c.degree());
              }
              er=pc.push_back(c, 1-t_split);
              if (er!=piecewise_curve_type::NO_ERRORS)
              {
                pc.clear();
                return false;
              }
            }

            return true;
          }

        private:
          static void build_cos_term(typename polynomial_type::coefficient_type &c, const data_type &a, const data_type &p, const data_type &xi0, index_type deg)
          {
            typename polynomial_type::coefficient_type cc(deg+1, 1), ctemp;
            data_type p_rel(p-xi0), temp(a*p_rel), d(1+temp*temp), d12(std::sqrt(d)), a2(a*a);
            index_type n, i, k;

            // create Taylor series around xi0
            n=0; cc(n)=1/d12;
            n=1; cc(n)=a2*p_rel/(d*d12);
            for (n=2; n<=deg; ++n)
            {
              cc(n)=(a2/(n*d))*((2*n-1)*p_rel*cc(n-1)-(n-1)*cc(n-2));
            }

            // convert expansion to standard form around zero
            ctemp=cc;
            cc.setZero();
            for (i=0; i<=deg; ++i)
            {
              for (k=i; k<=deg; ++k)
              {
                // NOTE: could probably get rid of the std::pow call by multiplying by -xi0 in loop
                eli::mutil::dm::n_choose_k(temp, k, i);
                cc(i)+=temp*std::pow(-xi0, k-i)*ctemp(k);
              }
            }

            // reparameterize on t^2
            c.resize(2*deg+1, 1);
            c.setZero();
            for (i=0; i<=deg; ++i)
            {
              c(2*i)=cc(i);
            }
          }

          static void build_sin_term(typename polynomial_type::coefficient_type &sterm, const data_type &a, const data_type &p, const typename polynomial_type::coefficient_type &cterm)
          {
            index_type n, deg((cterm.rows()-1)/2);
            sterm.resize(2*deg+3, 1); // sin term has one extra term in taylor series expansion

            sterm.setZero();
            n=0; sterm(n)=a*p*cterm(n);
            for (n=1; n<=deg; ++n)
            {
              sterm(2*n)=a*(p*cterm(2*n)-cterm(2*n-2));
            }
            n=deg+1; sterm(2*n)=-a*cterm(2*n-2);
          }

        private:
          airfoil_type af;
          index_type trig_degree;
      };
    }
  }
}
#endif
