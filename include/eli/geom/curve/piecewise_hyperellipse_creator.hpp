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
#include <algorithm>

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

            curve_type crv(max_degree);
            index_type nsegs(this->get_number_segments());
            error_code err;

            // set the start parameter
            pc.set_t0(this->get_t0());

            // TODO: REMOVE THIS AFTER TESTING
            assert(nsegs%4==0);

            // if even number of segments then can create top half
            if (nsegs%2==0)
            {
              piecewise_curve_type pc_bottom;
              index_type iseg, nseg_top(nsegs/2);

              // if have multiple of 4 segments then can create first quarter
              if (nsegs%4==0)
              {
                piecewise_curve_type pc_left;
                index_type nseg_first(nseg_top/2), nsample_pts(nseg_first*(max_degree)+1);
                std::vector<point_type> f(nsample_pts);

                // ensure that the first and last points are set appropriately
                f[0].setZero();
                f[0].x()=a;
                f[nsample_pts-1].setZero();
                f[nsample_pts-1].y()=b;

                // create rest of sample points
                for (index_type i=1; i<(nsample_pts-1); ++i)
                {
                  data_type t, argx, argy;

                  // to get the rapid parameterization variation, find the parameters that
                  // produce uniform spacing in x and y and then average the two
                  argx=std::pow(1-static_cast<data_type>(i)/(nsample_pts-1), m/2);
                  argx=std::min(argx, static_cast<data_type>(1));
                  argx=std::max(argx, static_cast<data_type>(0));
                  argy=std::pow(static_cast<data_type>(i)/(nsample_pts-1), n/2);
                  argy=std::min(argy, static_cast<data_type>(1));
                  argy=std::max(argy, static_cast<data_type>(0));

                  // TODO: Fix this to get the parameter that better captures extreme n and m cases
                  t=(std::acos(argx)+std::asin(argy))/(2*eli::constants::math<data__>::two_pi());

                  fun(f[i], t);
                }
#if 0
                if (nsegs==8)
                {
                  int i;
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
//                   assert(false);
                }
#endif

                // create first quarter
                for (iseg=0; iseg<nseg_first; ++iseg)
                {
                  control_point_type cp;
                  fit_container_type fcon;
                  crv.clear();

                  fcon.set_points(f.begin()+(iseg*max_degree), f.begin()+((iseg+1)*max_degree+1));

                  crv.interpolate(fcon);

                  // check the first slope to make sure it is reasonable
                  if (iseg==0)
                  {
                    bool need_set(false);
                    control_point_type cp(crv.get_control_point(1));

                    if (cp.y()<0)
                    {
                      cp.y()=0;
                      need_set=true;
                    }
                    if (cp.x()>a)
                    {
                      cp.x()=a;
                      need_set=true;
                    }
                    if (need_set)
                    {
                      crv.set_control_point(cp, 1);
                    }
                  }

                  // check the first slope to make sure it is reasonable
                  if (iseg==(nseg_first-1))
                  {
                    bool need_set(false);
                    control_point_type cp(crv.get_control_point(max_degree-1));

                    if (cp.x()<0)
                    {
                      cp.x()=0;
                      need_set=true;
                    }
                    if (cp.y()>b)
                    {
                      cp.y()=b;
                      need_set=true;
                    }
                    if (need_set)
                    {
                      crv.set_control_point(cp, max_degree-1);
                    }
                  }

                  err=pc.push_back(crv, this->get_segment_dt(iseg));
                  if (err!=piecewise_curve_type::NO_ERROR)
                  {
                    std::cout << "error number: " << err << std::endl;
                    assert(false);
                    pc.clear();
                    pc.set_t0(0);
                    return false;
                  }
                }

                // mirror to create second quarter
                pc_left=pc;
                pc_left.reflect_yz();
                pc_left.reverse();

//                   for (int ii=0; ii<=max_degree; ++ii)
//                   {
//                     std::cout << "ii = " << ii << "\tf_r=" << pc.f(ii/(1.0*max_degree))
//                                                << "\tf_l=" << pc_left.f(1-ii/(1.0*max_degree)) << std::endl;
//                   }
//                   assert(false);

                // push back left side
                err=pc.push_back(pc_left);
                if (err!=piecewise_curve_type::NO_ERROR)
                {
                  std::cout << "error number: " << err << std::endl;
                  assert(false);
                  pc.clear();
                  pc.set_t0(0);
                  return false;
                }
              }
              // else create entire top half
              else
              {
                assert(false);
                return false;
              }

              // mirror for bottom half
              pc_bottom=pc;
              pc_bottom.reflect_xz();
              pc_bottom.reverse();

              // push back the bottom curve
              err=pc.push_back(pc_bottom);
              if (err!=piecewise_curve_type::NO_ERROR)
              {
                std::cout << "error number: " << err << std::endl;
                assert(false);
                pc.clear();
                pc.set_t0(0);
                return false;
              }
            }
            // else odd number of segments
            else
            {
                assert(false);
                return false;
            }

            return true;
#if 0
            // only need to do one quarter and then mirror
            if (nsegs%4==0)
            {
              curve_type c0(max_degree), c1(max_degree), c2(max_degree), c3(max_degree);
              data_type t, argx, argy;
              control_point_type cp;
              fit_container_type fcon;
              std::vector<point_type> f(nsample_pts);

              for (i=0; i<nsample_pts; ++i)
              {
                // to get the rapid parameterization variation, find the parameters that
                // produce uniform spacing in x and y and then average the two
                argx=std::pow(1-static_cast<data_type>(i)/(nsample_pts-1), m/2);
                argx=std::min(argx, static_cast<data_type>(1));
                argx=std::max(argx, static_cast<data_type>(0));
                argy=std::pow(static_cast<data_type>(i)/(nsample_pts-1), n/2);
                argy=std::min(argy, static_cast<data_type>(1));
                argy=std::max(argy, static_cast<data_type>(0));
                t=(std::acos(argx)+std::asin(argy))/(2*eli::constants::math<data__>::two_pi());

                fun(f[i], t);
              }
              f[0].x()=a;
              f[0].y()=0;
              f[nsample_pts-1].x()=0;
              f[nsample_pts-1].y()=b;

              fcon.set_points(f.begin(), f.end());

              c0.interpolate(fcon);

              // make sure the second and second from last control points don't cross
              // x-axis and y-axis, respectively
              bool need_set(false);
              cp=c0.get_control_point(1);
              if (cp.y()<0)
              {
                cp.y()=0;
                need_set=true;
              }
              if (cp.x()>a)
              {
                cp.x()=a;
                need_set=true;
              }
              if (need_set)
              {
                c0.set_control_point(cp, 1);
              }

              need_set=false;
              cp=c0.get_control_point(max_degree-1);
              if (cp.x()<0)
              {
                cp.x()=0;
                need_set=true;
              }
              if (cp.y()>b)
              {
                cp.y()=b;
                need_set=true;
              }
              if (need_set)
              {
                c0.set_control_point(cp, max_degree-1);
              }
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
            }
            // need to do the entire
            else
            {
              // FIX: NOT IMPLEMENTED
              assert(false);
            }
#endif
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
