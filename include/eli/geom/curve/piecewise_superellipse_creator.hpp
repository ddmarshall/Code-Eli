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

#ifndef eli_geom_curve_piecewise_superellipse_creator_hpp
#define eli_geom_curve_piecewise_superellipse_creator_hpp

#include "eli/code_eli.hpp"

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
      class piecewise_2d_curve_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;

          piecewise_2d_curve_creator() : piecewise_creator_base<data_type, dim__, tolerance_type>(4, 0)
          {
            origin.setZero();
          }
          piecewise_2d_curve_creator(const index_type &ns, const data_type &tt0) : piecewise_creator_base<data_type, dim__, tolerance_type>(ns, tt0)
          {
            origin.setZero();
          }
          piecewise_2d_curve_creator(const piecewise_2d_curve_creator<data_type, dim__, tolerance_type> &p2dc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(p2dc), origin(p2dc.origin) {}
          ~piecewise_2d_curve_creator() {}

          void set_origin(const point_type &orig) {origin=orig;}
          point_type get_origin() const {return origin;}

        protected:
          virtual void fun(point_type &f, const data_type &t) const = 0;
          virtual void fun(point_type &f, point_type &fp, const data_type &t) const = 0;

        private:
          point_type origin;
      };


      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_superellipse_creator : public piecewise_2d_curve_creator<data__, dim__, tol__>
      {
        public:
          typedef typename piecewise_2d_curve_creator<data__, dim__, tol__>::data_type data_type;
          typedef typename piecewise_2d_curve_creator<data__, dim__, tol__>::index_type index_type;
          typedef typename piecewise_2d_curve_creator<data__, dim__, tol__>::point_type point_type;
          typedef typename piecewise_2d_curve_creator<data__, dim__, tol__>::tolerance_type tolerance_type;

          piecewise_superellipse_creator() : piecewise_2d_curve_creator<data_type, dim__, tolerance_type>(4, 0), a(1), b(1), m(2), n(2), m_bot(2), n_bot(2), max_degree(3) {}
          piecewise_superellipse_creator(const index_type &ns) : piecewise_2d_curve_creator<data_type, dim__, tolerance_type>(ns, 0), a(1), b(1), m(2), n(2), m_bot(2), n_bot(2), max_degree(3) {}
          piecewise_superellipse_creator(const piecewise_superellipse_creator<data_type, dim__, tolerance_type> &ppc)
            : piecewise_2d_curve_creator<data_type, dim__, tolerance_type>(ppc), a(ppc.a), b(ppc.b), m(ppc.m), n(ppc.n), m_bot(ppc.m_bot), n_bot(ppc.n_bot), max_degree(ppc.max_degree) {}
          ~piecewise_superellipse_creator() {}

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
            if (aa>=0)
            {
              a=aa;
            }
          }
          const data_type & get_a_axis() const {return a;}

          void set_b_axis(const data_type &bb)
          {
            if (bb>=0)
            {
              b=bb;
            }
          }
          const data_type & get_b_axis() const {return b;}

          void set_top_bot_sym(const data_type &tbsym)
          {
            top_bot_sym = tbsym;
          }

          void get_top_bot_sym(data_type &tbsym)
          {
            tbsym = top_bot_sym;
          }

          void set_exponents(const data_type &mm, const data_type &nn)
          {
            set_m_exponent(mm);
            set_n_exponent(nn);
          }

          void set_exponents_bot(const data_type &mm, const data_type &nn)
          {
            set_m_exponent_bot(mm);
            set_n_exponent_bot(nn);
          }

          void get_exponents(data_type &mm, data_type &nn) const
          {
            mm = m;
            nn = n;
          }

          void get_exponents_bot(data_type &mm, data_type &nn) const
          {
            mm = m_bot;
            nn = n_bot;
          }

          void set_m_exponent(const data_type &mm)
          {
            if (mm>0)
            {
              m=mm;
            }
          }
          const data_type & get_m_exponent() const {return m;}

          void set_m_exponent_bot(const data_type &mm)
          {
            if (mm > 0)
            {
              m_bot = mm;
            }
          }

          const data_type & get_m_exponent_bot() const
          {
            return m_bot;
          }

          void set_n_exponent(const data_type &nn)
          {
            if (nn>0)
            {
              n=nn;
            }
          }
          const data_type & get_n_exponent() const {return n;}

          void set_n_exponent_bot(const data_type &nn)
          {
            if (nn > 0)
            {
              n_bot = nn;
            }
          }

          const data_type & get_n_exponent_bot() const
          {
            return n_bot;
          }

          void set_max_degree(const index_type &md)
          {
            if (md>2)
            {
              max_degree=md;
            }
          }
          index_type get_max_degree() const {return max_degree;}

          void set_max_width_loc(const data_type &mw)
          {
            max_width_loc = mw;
          }

          const data_type & get_max_width_loc() const
          {
            return max_width_loc;
          }

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
            assert(nsegs%2==0);

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
                f[0].y() = max_width_loc;
                f[nsample_pts-1].setZero();
                f[nsample_pts-1].y()=b;

                // create rest of sample points
                for (index_type i=1; i<(nsample_pts-1); ++i)
                {
                  data_type t, argx, argy, tmp;

                  // to get the rapid parameterization variation, find the parameters that
                  // produce uniform spacing in x and y and then average the two
                  tmp=static_cast<data_type>(i)/(nsample_pts-1);
                  argx=std::pow(1-tmp, m/2);
                  argx=std::min(argx, static_cast<data_type>(1));
                  argx=std::max(argx, static_cast<data_type>(0));
                  argy=std::pow(tmp, n/2);
                  argy=std::min(argy, static_cast<data_type>(1));
                  argy=std::max(argy, static_cast<data_type>(0));

                  // TODO: Fix this to get the parameter that better captures extreme n and m cases
                  t=(std::acos(argx)+std::asin(argy))/(2*eli::constants::math<data__>::two_pi());

                  fun(f[i], t);
                }

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

                    if (cp.y() < max_width_loc)
                    {
                      cp.y() = max_width_loc;
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

                  // check the last slope to make sure it is reasonable
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
                  if (err!=piecewise_curve_type::NO_ERRORS)
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

                // push back left side
                curve_type c;
                index_type rnseg(pc.number_segments());
                for (iseg=0; iseg<pc_left.number_segments(); ++iseg)
                {
                  pc_left.get(c, iseg);
                  err=pc.push_back(c, this->get_segment_dt(rnseg+iseg));
                  if (err!=piecewise_curve_type::NO_ERRORS)
                  {
                    std::cout << "error number: " << err << std::endl;
                    assert(false);
                    pc.clear();
                    pc.set_t0(0);
                    return false;
                  }
                }
              }
              // else create entire top half
              else
              {
                index_type nsample_pts(nseg_top*(max_degree)+1), nhalf(nsample_pts/2);
                std::vector<point_type> f(nsample_pts);

                // ensure that the first and last points are set appropriately
                f[0].setZero();
                f[0].x() = a;
                f[0].y() = max_width_loc;
                f[nsample_pts - 1].setZero();
                f[nsample_pts - 1].x() = -a;
                f[nsample_pts - 1].y() = max_width_loc;

                if (nsample_pts%2==1)
                {
                  f[nhalf].setZero();
                  f[nhalf].y()=b;
                }

                // create rest of sample points
                for (index_type i=1; i<nhalf; ++i)
                {
                  data_type t, argx, argy, tmp;

                  // to get the rapid parameterization variation, find the parameters that
                  // produce uniform spacing in x and y and then average the two
                  tmp=static_cast<data_type>(i)/(nhalf);
                  argx=std::pow(1-tmp, m/2);
                  argx=std::min(argx, static_cast<data_type>(1));
                  argx=std::max(argx, static_cast<data_type>(0));
                  argy=std::pow(tmp, n/2);
                  argy=std::min(argy, static_cast<data_type>(1));
                  argy=std::max(argy, static_cast<data_type>(0));

                  // TODO: Fix this to get the parameter that better captures extreme n and m cases
                  t=(std::acos(argx)+std::asin(argy))/(2*eli::constants::math<data__>::two_pi());
//                     std::cout << "argx=" << argx << "\targy=" << argy << "\tt=" << t << std::endl;

                  fun(f[i], t);
                  fun(f[nsample_pts-i-1], static_cast<data_type>(0.5)-t);
                }

                // create first quarter
                for (iseg=0; iseg<nseg_top; ++iseg)
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

                    if (cp.y()<max_width_loc)
                    {
                      cp.y() = max_width_loc;
                      need_set = true;
                    }
                    if (cp.x()>a)
                    {
                      cp.x()=a;
                      need_set = true;
                    }
                    if (need_set)
                    {
                      crv.set_control_point(cp, 1);
                    }
                  }

                  // check the last slope to make sure it is reasonable
                  if (iseg==(nseg_top-1))
                  {
                    bool need_set(false);
                    control_point_type cp(crv.get_control_point(max_degree-1));

                    if (cp.y() < max_width_loc)
                    {
                      cp.y() = max_width_loc;
                      need_set=true;
                    }
                    if (cp.x()<-a)
                    {
                      cp.x()=-a;
                      need_set=true;
                    }
                    if (need_set)
                    {
                      crv.set_control_point(cp, max_degree-1);
                    }
                  }

                  err=pc.push_back(crv, this->get_segment_dt(iseg));
                  if (err!=piecewise_curve_type::NO_ERRORS)
                  {
                    std::cout << "error number: " << err << std::endl;
                    assert(false);
                    pc.clear();
                    pc.set_t0(0);
                    return false;
                  }
                }
              }

              // if top_bot_sym is true, reflect top curve
              if (top_bot_sym)
              {
                // mirror for bottom half
                pc_bottom = pc;
                pc_bottom.reflect_xz();
                pc_bottom.reverse();

                // push back the bottom curve
                curve_type c;
                index_type rtseg(pc.number_segments());
                for (iseg = 0; iseg < pc_bottom.number_segments(); ++iseg)
                {
                  pc_bottom.get(c, iseg);
                  err = pc.push_back(c, this->get_segment_dt(rtseg + iseg));
                  if (err != piecewise_curve_type::NO_ERRORS)
                  {
                    std::cout << "error number: " << err << std::endl;
                    assert(false);
                    pc.clear();
                    pc.set_t0(0);
                    return false;
                  }
                }
              }
              else
              {
                // otherwise, generate bot curve separately
                index_type nseg_bot(nsegs - nseg_top), nsample_pts_bot(nseg_bot * (max_degree) + 1), nhalf_bot(nsample_pts_bot / 2);

                std::vector<point_type> f_bot(nsample_pts_bot);

                f_bot[0].setZero();
                f_bot[0].x() = -a;
                f_bot[0].y() = max_width_loc;

                //final point should be at origin
                f_bot[nsample_pts_bot - 1].setZero();
                f_bot[nsample_pts_bot - 1].x() = a;
                f_bot[nsample_pts_bot - 1].y() = max_width_loc;

                if (nsample_pts_bot % 2 == 1)
                {
                  f_bot[nhalf_bot].setZero();

                   //half way will be at -b
                  f_bot[nhalf_bot].y() = -b;
                }

                for (index_type i = 1; i < nhalf_bot; ++i)
                {
                  data_type t, argx, argy, tmp;

                  // to get the rapid parameterization variation, find the parameters that
                  // produce uniform spacing in x and y and then average the two
                  tmp = static_cast<data_type> (i) / (nhalf_bot);

                  argx = std::pow(1 - tmp, m_bot / 2);
                  argx = std::min(argx, static_cast<data_type> (1));
                  argx = std::max(argx, static_cast<data_type> (0));
                  argy = std::pow(tmp, n_bot / 2);
                  argy = std::min(argy, static_cast<data_type> (1));
                  argy = std::max(argy, static_cast<data_type> (0));

                  // TODO: Fix this to get the parameter that better captures extreme n and m cases
                  t = (std::acos(argx) + std::asin(argy)) / (2 * eli::constants::math<data__>::two_pi());
                  //std::cout << "argx=" << argx << "\targy=" << argy << "\tt=" << t << std::endl;

                  // parameterization will need to change as shown below
                  fun_bot(f_bot[i], static_cast<data_type> (0.5) + t);
                  fun_bot(f_bot[nsample_pts_bot - i - 1], static_cast<data_type> (1.0) - t);

                }
                //Create bottom half
                for (iseg = 0; iseg < nseg_bot; ++iseg)
                {
                  control_point_type cp;
                  fit_container_type fcon_bot;

                  crv.clear(); //crv is a max_degree bezier curve here

                  //sets the control points in segments of max_degree (so for a 3rd degree, 4 points per)
                  fcon_bot.set_points(f_bot.begin()+(iseg * max_degree), f_bot.begin()+((iseg + 1) * max_degree + 1));

                  //this applies the blending function to interpolate and produce a smooth curve
                  crv.interpolate(fcon_bot);

                  // check the first slope to make sure it is reasonable
                  if (iseg == 0)
                  {
                    bool need_set(false);
                    control_point_type cp_bot(crv.get_control_point(0));

                    if (cp_bot.y() < max_width_loc)
                    {
                      cp_bot.y() = max_width_loc;
                      need_set = true;
                    }
                    if (cp_bot.x() < -a)
                    {
                      cp_bot.x() = -a;
                      need_set = true;
                    }
                    if (need_set) {
                      crv.set_control_point(cp_bot, 1);
                    }
                  }

                  // check the last slope to make sure it is reasonable
                  if (iseg == (nseg_bot - 1))
                  {
                    bool need_set(false);
                    control_point_type cp_bot(crv.get_control_point(max_degree));

                    if (cp_bot.y() < max_width_loc)
                    {
                      cp_bot.y() = max_width_loc;
                      need_set = true;
                    }
                    if (cp_bot.x() > a)
                    {
                      cp_bot.x() = a;
                       need_set = true;
                    }
                    if (need_set)
                    {
                      crv.set_control_point(cp_bot, max_degree - 1);
                    }
                  }

                  err = pc.push_back(crv, this->get_segment_dt(iseg));
                  if (err != piecewise_curve_type::NO_ERRORS)
                  {
                    std::cout << "error number: " << err << std::endl;
                    assert(false);
                    pc.clear();
                    pc.set_t0(0);
                    return false;
                  }
                }
              }
            }
            // else odd number of segments
            else
            {
                assert(false);
                return false;
            }

            // translate to center on origin
            pc.translate(this->get_origin());

            return true;
          }

        protected:
          // TODO: These should be in the base class that is shared with circle, ellipse, and other 2D curves
          void set_x_radius(const data_type &xr)
          {
            if(xr>=0)
            {
              xradius=xr;
            }
            else
            {
              assert(false);
            }
          }
          const data_type & get_x_radius() const {return xradius;}
          void set_y_radius(const data_type &yr)
          {
            if (yr>=0)
            {
              yradius=yr;
            }
            else
            {
              assert(false);
            }
          }
          const data_type & get_y_radius() const {return yradius;}
          // TODO: end note

        protected:
          virtual void fun(point_type &f, const data_type &t) const
          {
            // short circuit if given bad parameter
            if ((t<0) || (t>1))
              return;

            // zero out point
            f.setZero();

            // determine the sign terms
            data_type sign_cos(1), sign_sin(1);

            if ((t > 0.25) && (t <= 0.5))
            {
              sign_cos=-1;
            }
            else if ((t > 0.5) && (t <= 0.75))
            {
              sign_cos=-1;
              sign_sin=-1;
            }
            else if (t > 0.75)
            {
              sign_sin=-1;
            }

            data_type theta(eli::constants::math<data__>::two_pi()*t);
            data_type abs_cos(std::abs(std::cos(theta)));
            data_type abs_sin(std::abs(std::sin(theta)));

            // calculate the function
            f.x() = a * sign_cos * std::pow(abs_cos, 2 / m);
            f.y() = max_width_loc + (b - max_width_loc) * sign_sin * std::pow(abs_sin, 2 / n);
          }

          virtual void fun(point_type &f, point_type &fp, const data_type &t) const
          {
            // short circuit if given bad parameter
            if ((t<0) || (t>1))
              return;

            // zero out point
            f.setZero();
            fp.setZero();

            // determine the sign terms
            data_type sign_cos(1), sign_sin(1);

            if ((t>0.25) && (t<=0.5))
            {
              sign_cos=-1;
            }
            else if ((t>0.5) && (t<=0.75))
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

            // calculate the 1st derivatives
            fp.x()=(-2*eli::constants::math<data__>::two_pi()*a/m)*std::sin(theta)*std::pow(abs_cos, 2/m-1);
            fp.y()=(2*eli::constants::math<data__>::two_pi()*b/n)*std::cos(theta)*std::pow(abs_sin, 2/n-1);
          }

          virtual void fun_bot(point_type &f, const data_type &t) const
          {
            // short circuit if given bad parameter
            if ((t < 0) || (t > 1))
            {
              return;
            }

            // zero out point
            f.setZero();

            // determine the sign terms
            data_type sign_cos(1), sign_sin(1);

            if ((t>0.25) && (t<=0.5))
            {
              sign_cos=-1;
            }
            else  if ((t>0.5) && (t<=0.75))
            {
              sign_cos=-1;
              sign_sin=-1;
            }
            else if (t>0.75)
            {
              sign_sin=-1;
            }

            data_type theta(eli::constants::math<data__>::two_pi() * t);
            data_type abs_cos(std::abs(std::cos(theta)));
            data_type abs_sin(std::abs(std::sin(theta)));

            // calculate the function
            f.x() = a * sign_cos * std::pow(abs_cos, 2 / m_bot);
            f.y() = max_width_loc + (b + max_width_loc) * sign_sin * std::pow(abs_sin, 2 / n_bot);
            //std::cout << "f.x=" << f.x() << "\tf.y=" << f.y() << std::endl;
          }

        private:

          data_type xradius, yradius;
          data_type a, b, m, n;
          data_type m_bot, n_bot;
          data_type max_width_loc;
          data_type top_bot_sym;
          index_type max_degree;
      };
    }
  }
}
#endif
