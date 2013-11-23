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

#ifndef eli_geom_curve_piecewise_hyperellipse_creator_hpp
#define eli_geom_curve_piecewise_hyperellipse_creator_hpp

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
      class piecewise_hyperellipse_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef data__  data_type;
          typedef int index_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef tol__ tolerance_type;

          piecewise_hyperellipse_creator() : piecewise_creator_base<data_type, dim__, tolerance_type>(4, 0), a(1), b(1), m(2), n(2), max_degree(3) {}
          piecewise_hyperellipse_creator(const index_type &ns) : piecewise_creator_base<data_type, dim__, tolerance_type>(ns, 0), a(1), b(1), m(2), n(2), max_degree(3) {}
          piecewise_hyperellipse_creator(const piecewise_hyperellipse_creator<data_type, dim__, tolerance_type> &ppc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(ppc), a(ppc.a), b(ppc.b), m(ppc.m), n(ppc.n), max_degree(ppc.max_degree) {}
          ~piecewise_hyperellipse_creator() {}

          void set_axis(const data_type &aa, const data_type &bb)
          {
            set_a_axis(aa);
            set_b_axis(bb);
          }
          void get_axis(data_type &aa, data_type &bb) const
          {
            aa = a;
            bb = b;
          }

          void set_a_axis(const data_type &aa)
          {
            if (aa>0)
            {
              a=aa;
            }
          }
          const data_type & get_a_axis() const {return a;}

          void set_b_axis(const data_type &bb)
          {
            if (bb>0)
            {
              b=bb;
            }
          }
          const data_type & get_b_axis() const {return b;}

          void set_exponents(const data_type &mm, const data_type &nn)
          {
            set_m_exponent(mm);
            set_n_exponent(nn);
          }
          void get_exponents(data_type &mm, data_type &nn) const
          {
            mm = m;
            nn = n;
          }

          void set_m_exponent(const data_type &mm)
          {
            if (mm>0)
            {
              m=mm;
            }
          }
          const data_type & get_m_exponent() const {return m;}

          void set_n_exponent(const data_type &nn)
          {
            if (nn>0)
            {
              n=nn;
            }
          }
          const data_type & get_n_exponent() const {return n;}

          void set_max_degree(const index_type &md)
          {
            if (md>2)
            {
              max_degree=md;
            }
          }
          index_type get_max_degree() const {return max_degree;}

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::control_point_type control_point_type;
            typedef typename piecewise_curve_type::error_code error_code;
            typedef typename curve_type::fit_container_type fit_container_type;

            pc.clear();

            // TODO: should this be an adaptive thing instead of just start with max?

            index_type i, nsegs(this->get_number_segments()), nsample_pts(max_degree+1);
            error_code err;

            // set the start parameter
            pc.set_t0(this->get_t0());

            // TODO: REMOVE THIS AFTER TESTING
            assert(nsegs==4);

            // only need to do one quarter and then mirror
            if (nsegs%4==0)
            {
              curve_type c0(max_degree), c1(max_degree), c2(max_degree), c3(max_degree);
              data_type t;
              control_point_type cp;
              fit_container_type fcon;
              std::vector<point_type> f(nsample_pts);

              for (i=0; i<nsample_pts; ++i)
              {
                t=(0.25*i)/(nsample_pts-1);
                fun(f[i], t);
              }
              f[0].x()=a;
              f[0].y()=0;
              f[nsample_pts-1].x()=0;
              f[nsample_pts-1].y()=b;

              fcon.set_points(f.begin(), f.end());

              c0.interpolate(fcon);
              err=pc.push_back(c0, this->get_segment_dt(0));
              if (err!=piecewise_curve_type::NO_ERROR)
              {
                assert(false);
                pc.clear();
                pc.set_t0(0);
                return false;
              }

              // mirror c0 about y-axis
              for (i=0; i<=max_degree; ++i)
              {
                cp=c0.get_control_point(i);
                cp.x()=-cp.x();
                c1.set_control_point(cp, i);
              }
              c1.reverse();
              err=pc.push_back(c1, this->get_segment_dt(1));
              if (err!=piecewise_curve_type::NO_ERROR)
              {
                std::cout << "error code=" << err << std::endl;
                assert(false);
                pc.clear();
                pc.set_t0(0);
                return false;
              }

              // mirror c1 about x-axis
              for (i=0; i<=max_degree; ++i)
              {
                cp=c1.get_control_point(i);
                cp.y()=-cp.y();
                c2.set_control_point(cp, i);
              }
              c2.reverse();
              err=pc.push_back(c2, this->get_segment_dt(2));
              if (err!=piecewise_curve_type::NO_ERROR)
              {
                std::cout << "error code=" << err << std::endl;
                assert(false);
                pc.clear();
                pc.set_t0(0);
                return false;
              }

              // mirror c2 about y-axis
              for (i=0; i<=max_degree; ++i)
              {
                cp=c2.get_control_point(i);
                cp.x()=-cp.x();
                c3.set_control_point(cp, i);
              }
              c3.reverse();
              err=pc.push_back(c3, this->get_segment_dt(3));
              if (err!=piecewise_curve_type::NO_ERROR)
              {
                std::cout << "error code=" << err << std::endl;
                assert(false);
                pc.clear();
                pc.set_t0(0);
                return false;
              }
            }
            // only need to do one half and then mirror
            else if (nsegs%2==2)
            {
              // FIX: NOT IMPLEMENTED
              assert(false);
#if 0
              std::cout << "fx=[";
              for (i=0; i<nsample_pts; ++i)
              {
                std::cout << f[i].x() << " ";
              }
              std::cout << "];" << std::endl;
              std::cout << "fy=[";
              for (i=0; i<nsample_pts; ++i)
              {
                std::cout << f[i].y() << " ";
              }
              std::cout << "];" << std::endl;
              std::cout << "fz=[";
              for (i=0; i<nsample_pts; ++i)
              {
                std::cout << f[i].z() << " ";
              }
              std::cout << "];" << std::endl;
              std::cout << "fpx=[";
              for (i=0; i<nsample_pts; ++i)
              {
                std::cout << fp[i].x() << " ";
              }
              std::cout << "];" << std::endl;
              std::cout << "fpy=[";
              for (i=0; i<nsample_pts; ++i)
              {
                std::cout << fp[i].y() << " ";
              }
              std::cout << "];" << std::endl;
              std::cout << "fpz=[";
              for (i=0; i<nsample_pts; ++i)
              {
                std::cout << fp[i].z() << " ";
              }
              std::cout << "];" << std::endl;
#endif
            }
            // need to do the entire
            else
            {
              // FIX: NOT IMPLEMENTED
              assert(false);
            }

            return true;
          }

        private:
          void fun(point_type &f, const data_type &t) const
          {
            // short circuit if given bad parameter
            if ((t<0) || (t>1))
              return;

            // zero out the points
            f.setZero();

            // determine the sign terms
            data_type sign_cos(1), sign_sin(1);

            if ((t>0.25) && (t<0.5))
            {
              sign_cos=-1;
            }
            else if ((t>0.5) && (t>0.75))
            {
              sign_cos=-1;
              sign_sin=-1;
            }
            else if (t>0.75)
            {
              sign_sin=-1;
            }

            data_type theta(eli::constants::math<data__>::two_pi()*t);
            data_type abs_cos(std::abs(std::cos(theta)));
            data_type abs_sin(std::abs(std::sin(theta)));

            // calculate the function
            f.x()=a*sign_cos*std::pow(abs_cos, 2/m);
            f.y()=b*sign_sin*std::pow(abs_sin, 2/n);

#if 0
            // calculate the 1st derivatives
            fp.x()=(-2*eli::constants::math<data__>::two_pi()*a/m)*std::sin(theta)*std::pow(abs_cos, 2/m-1);
            fp.y()=(2*eli::constants::math<data__>::two_pi()*b/n)*std::cos(theta)*std::pow(abs_sin, 2/n-1);
#endif
          }

        private:
          data_type a, b, m, n;
          index_type max_degree;
      };
    }
  }
}
#endif
