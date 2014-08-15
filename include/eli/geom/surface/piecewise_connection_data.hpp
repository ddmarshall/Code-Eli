/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_surface_piecewise_connection_data_hpp
#define eli_geom_surface_piecewise_connection_data_hpp

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise.hpp"

#include <iterator>

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class connection_data
      {
        public:
          enum
          {
            CONNECTION_SET=0x000001, // must be set for valid joint
            LEFT_FP_SET   =0x000010,
            RIGHT_FP_SET  =0x000100,
            LEFT_FPP_SET  =0x001000,
            RIGHT_FPP_SET =0x010000
          };

          enum connection_continuity
          {
            C0=general::C0,
            C1=general::C1,
            C2=general::C2
          };

          typedef data__ data_type;
          typedef tol__ tolerance_type;
          typedef typename eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, dim__, tol__> curve_type;
          typedef typename curve_type::index_type index_type;

        public:
          connection_data() : conditions(0), continuity(C0)
          {
          }

          connection_data(const connection_data &cd)
            : f(cd.f), fp_left(cd.fp_left), fp_right(cd.fp_right), fpp_left(cd.fpp_left),
              fpp_right(cd.fpp_right), conditions(cd.conditions), continuity(cd.continuity)
          {
          }

          ~connection_data() {}

          const connection_data & operator=(const connection_data &cd)
          {
            if (this!=&cd)
            {
              f=cd.f;
              fp_left=cd.fp_left;
              fp_right=cd.fp_right;
              fpp_left=cd.fpp_left;
              fpp_right=cd.fpp_right;
              conditions=cd.conditions;
              continuity=cd.continuity;
            }

            return (*this);
          }

          bool operator==(const connection_data &cd) const
          {
            tolerance_type tol;

            if (conditions!=cd.conditions)
              return false;
            if (continuity!=cd.continuity)
              return false;
            if (cd.use_f())
            {
              if (f!=cd.f)
                return false;
            }
            if (cd.use_left_fp())
            {
              if (fp_left!=cd.fp_left)
                return false;
            }
            if (cd.use_right_fp())
            {
              if (fp_right!=cd.fp_right)
                return false;
            }
            if (cd.use_left_fpp())
            {
              if (fpp_left!=cd.fpp_left)
                return false;
            }
            if (cd.use_right_fpp())
            {
              if (fpp_right!=cd.fpp_right)
                return false;
            }

            return true;
          }

          bool operator!=(const connection_data &cd) const {return !operator==(cd);}

          data_type get_t0() const
          {
            return f.get_t0();
          }
          data_type get_tmax() const
          {
            return f.get_tmax();
          }

          // get joints for all active curves
          template<typename output_it__>
          void get_joints(output_it__ it_in) const
          {
            tolerance_type tol;
            auto comp = [&tol](const data_type &x1, const data_type &x2)->bool
            {
              return tol.approximately_less_than(x1, x2);
            };

            std::vector<data_type> jts1, jts, jts_out;
            if (use_f())
            {
              f.get_parameters(std::back_inserter(jts1));
              if (jts.empty())
              {
                std::swap(jts, jts1);
              }
              else
              {
                std::set_union(jts.begin(), jts.end(), jts1.begin(), jts1.end(), std::back_inserter(jts_out), comp);
                std::swap(jts, jts_out);
                jts1.clear();
                jts_out.clear();
              }
            }
            if (use_left_fp())
            {
              fp_left.get_parameters(std::back_inserter(jts1));
              if (jts.empty())
              {
                std::swap(jts, jts1);
              }
              else
              {
                std::set_union(jts.begin(), jts.end(), jts1.begin(), jts1.end(), std::back_inserter(jts_out), comp);
                std::swap(jts, jts_out);
                jts1.clear();
                jts_out.clear();
              }
            }
            if (use_right_fp())
            {
              fp_right.get_parameters(std::back_inserter(jts1));
              if (jts.empty())
              {
                std::swap(jts, jts1);
              }
              else
              {
                std::set_union(jts.begin(), jts.end(), jts1.begin(), jts1.end(), std::back_inserter(jts_out), comp);
                std::swap(jts, jts_out);
                jts1.clear();
                jts_out.clear();
              }
            }
            if (use_left_fpp())
            {
              fpp_left.get_parameters(std::back_inserter(jts1));
              if (jts.empty())
              {
                std::swap(jts, jts1);
              }
              else
              {
                std::set_union(jts.begin(), jts.end(), jts1.begin(), jts1.end(), std::back_inserter(jts_out), comp);
                std::swap(jts, jts_out);
                jts1.clear();
                jts_out.clear();
              }
            }
            if (use_right_fpp())
            {
              fpp_right.get_parameters(std::back_inserter(jts1));
              if (jts.empty())
              {
                std::swap(jts, jts1);
              }
              else
              {
                std::set_union(jts.begin(), jts.end(), jts1.begin(), jts1.end(), std::back_inserter(jts_out), comp);
                std::swap(jts, jts_out);
                jts1.clear();
                jts_out.clear();
              }
            }

            std::copy(jts.begin(), jts.end(), it_in);
          }

          // split all curves
          template<typename it1__, typename it2__>
          void split(it1__ itb, it1__ ite, it2__ itd)
          {
            index_type i, njoints(0);
            for (it1__ it=itb; it!=ite; ++it, ++njoints)
            {
              if (use_f())
              {
                f.split(*it);
              }
              if (use_left_fp())
              {
                fp_left.split(*it);
              }
              if (use_right_fp())
              {
                fp_right.split(*it);
              }
              if (use_left_fpp())
              {
                fpp_left.split(*it);
              }
              if (use_right_fpp())
              {
                fpp_right.split(*it);
              }
            }

            index_type nsegs(njoints-1);
            std::vector<index_type> deg(nsegs), max_deg(nsegs, 0);
            if (use_f())
            {
              f.degrees(deg.begin());
              for (i=0; i<nsegs; ++i)
              {
                if (deg[i]>max_deg[i])
                {
                  max_deg[i]=deg[i];
                }
              }
            }
            if (use_left_fp())
            {
              fp_left.degrees(deg.begin());
              for (i=0; i<nsegs; ++i)
              {
                if (deg[i]>max_deg[i])
                {
                  max_deg[i]=deg[i];
                }
              }
            }
            if (use_right_fp())
            {
              fp_right.degrees(deg.begin());
              for (i=0; i<nsegs; ++i)
              {
                if (deg[i]>max_deg[i])
                {
                  max_deg[i]=deg[i];
                }
              }
            }
            if (use_left_fpp())
            {
              fpp_left.degrees(deg.begin());
              for (i=0; i<nsegs; ++i)
              {
                if (deg[i]>max_deg[i])
                {
                  max_deg[i]=deg[i];
                }
              }
            }
            if (use_right_fpp())
            {
              fpp_right.degrees(deg.begin());
              for (i=0; i<nsegs; ++i)
              {
                if (deg[i]>max_deg[i])
                {
                  max_deg[i]=deg[i];
                }
              }
            }

            for (i=0; i<nsegs; ++i, ++itd)
            {
              (*itd)=max_deg[i];
            }
          }

          template<typename it__>
          bool promote(it__ itb, it__ ite)
          {
            // NOTE: the piecewise curve class will return error if degree vector is
            //       larger than number of segments
            typename curve_type::error_code ec;
            it__ it;
            index_type i;
            for (i=0, it=itb; it!=ite; ++it, ++i)
            {
              if (use_f())
              {
                ec=f.degree_promote_to(i, (*it));
                if (ec!=curve_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              if (use_left_fp())
              {
                ec=fp_left.degree_promote_to(i, (*it));
                if (ec!=curve_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              if (use_right_fp())
              {
                ec=fp_right.degree_promote_to(i, (*it));
                if (ec!=curve_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              if (use_left_fpp())
              {
                ec=fpp_left.degree_promote_to(i, (*it));
                if (ec!=curve_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
              if (use_right_fpp())
              {
                ec=fpp_right.degree_promote_to(i, (*it));
                if (ec!=curve_type::NO_ERRORS)
                {
                  assert(false);
                  return false;
                }
              }
            }

            return true;
          }
          
          // connection interface
          bool set_f(const curve_type &ff)
          {
            f=ff;
            conditions|=CONNECTION_SET;
            return check_state();
          }
          const curve_type & get_f() const {return f;}
          bool unset_f()
          {
            conditions&=~CONNECTION_SET;
            return check_state();
          }
          bool use_f() const
          {
            return ((conditions & CONNECTION_SET) == CONNECTION_SET);
          }

          // first derivative interface
          bool set_left_fp(const curve_type &fpl)
          {
            fp_left=fpl;
            conditions|=LEFT_FP_SET;

            // if already wanting C1 continuous then set right as well
            if (get_continuity()>C0)
            {
              fp_right=fpl;
              conditions|=RIGHT_FP_SET;
            }

            return check_state();
          }
          bool set_right_fp(const curve_type &fpr)
          {
            fp_right=fpr;
            conditions|=RIGHT_FP_SET;

            // if already wanting C1 continuous then set left as well
            if (get_continuity()>C0)
            {
              fp_left=fpr;
              conditions|=LEFT_FP_SET;
            }

            return check_state();
          }
          bool set_fp(const curve_type &p)
          {
            conditions|=(LEFT_FP_SET | RIGHT_FP_SET);
            fp_left=p;
            fp_right=p;
            return check_state();
          }

          const curve_type & get_left_fp() const
          {
            return fp_left;
          }
          const curve_type & get_right_fp() const
          {
            return fp_right;
          }
          void get_fp(curve_type &fpl, curve_type &fpr) const
          {
            fpl=fp_left;
            fpr=fp_right;
          }

          bool unset_fp()
          {
            conditions&=~(LEFT_FP_SET | RIGHT_FP_SET);
            return check_state();
          }
          bool unset_left_fp()
          {
            conditions&=~LEFT_FP_SET;

            return check_state();
          }
          bool unset_right_fp()
          {
            conditions&=~RIGHT_FP_SET;

            return check_state();
          }

          bool use_left_fp() const
          {
            return ((conditions & LEFT_FP_SET) == LEFT_FP_SET);
          }
          bool use_right_fp() const
          {
            return ((conditions & RIGHT_FP_SET) == RIGHT_FP_SET);
          }

          // second derivative interface
          bool set_left_fpp(const curve_type &fppl)
          {
            fpp_left=fppl;
            conditions|=LEFT_FPP_SET;

            // if already wanting C2 continuous then set right as well
            if (get_continuity()>C1)
            {
              fpp_right=fppl;
              conditions|=RIGHT_FPP_SET;
            }

            return check_state();
          }
          bool set_right_fpp(const curve_type &fppr)
          {
            fpp_right=fppr;
            conditions|=RIGHT_FPP_SET;

            // if already wanting C2 continuous then set left as well
            if (get_continuity()>C1)
            {
              fpp_left=fppr;
              conditions|=LEFT_FPP_SET;
            }

            return check_state();
          }
          bool set_fpp(const curve_type &p)
          {
            conditions|=(LEFT_FPP_SET | RIGHT_FPP_SET);
            fpp_left=p;
            fpp_right=p;
            return check_state();
          }

          const curve_type & get_left_fpp() const
          {
            return fpp_left;
          }
          const curve_type & get_right_fpp() const
          {
            return fpp_right;
          }
          void get_fpp(curve_type &fppl, curve_type &fppr) const
          {
            fppl=fpp_left;
            fppr=fpp_right;
          }

          bool unset_fpp()
          {
            conditions&=~(LEFT_FPP_SET | RIGHT_FPP_SET);
            return check_state();
          }
          bool unset_left_fpp()
          {
            conditions&=~LEFT_FPP_SET;

            return check_state();
          }
          bool unset_right_fpp()
          {
            conditions&=~RIGHT_FPP_SET;

            return check_state();
          }

          bool use_left_fpp() const
          {
            return ((conditions & LEFT_FPP_SET) == LEFT_FPP_SET);
          }
          bool use_right_fpp() const
          {
            return ((conditions & RIGHT_FPP_SET) == RIGHT_FPP_SET);
          }

          bool set_continuity(connection_continuity jc)
          {
            continuity=jc;
            return check_state();
          }
          connection_continuity get_continuity() const
          {
            return continuity;
          }

          bool check_state() const
          {
            tolerance_type tol;

            // check if point set
            if ((conditions & CONNECTION_SET)==0)
              return false;

            // if highest continuity is C0 then done
            if (continuity==C0)
              return true;

            // check first derivatives match on both sides
            if ((conditions & (LEFT_FP_SET | RIGHT_FP_SET)) == (LEFT_FP_SET | RIGHT_FP_SET))
            {
              // TODO: make this comparison an apprimately_equal() call
//              if (!tol.approximately_equal(fp_left, fp_right))
              if (fp_left!=fp_right)
                return false;
            }
            // check neither side first derivatives set set
            else if ((conditions & (LEFT_FP_SET | RIGHT_FP_SET)) != 0)
            {
              return false;
            }

            // check that the parameterization of fp is same as f
            data_type t0(get_t0()), tmax(get_tmax());
            if (use_left_fp())
            {
              if (!tol.approximately_equal(t0, fp_left.get_t0()) || !tol.approximately_equal(tmax, fp_left.get_tmax()))
                return false;
            }
            if (use_right_fp())
            {
              if (!tol.approximately_equal(t0, fp_right.get_t0()) || !tol.approximately_equal(tmax, fp_right.get_tmax()))
                return false;
            }

            // if highest continuity is C1 then done
            if (continuity==C1)
              return true;

            // check second derivatives match on both sides
            if ((conditions & (LEFT_FPP_SET | RIGHT_FPP_SET)) == (LEFT_FPP_SET | RIGHT_FPP_SET))
            {
              // TODO: make this comparison an apprimately_equal() call
//              if (!tol.approximately_equal(fpp_left, fpp_right))
              if (fpp_left!=fpp_right)
                return false;
            }
            // check neither side second derivatives set set
            else if ((conditions & (LEFT_FPP_SET | RIGHT_FPP_SET)) != 0)
            {
              return false;
            }

            // check that the parameterization of fpp is same as f
            if (use_left_fpp())
            {
              if (!tol.approximately_equal(t0, fpp_left.get_t0()) || !tol.approximately_equal(tmax, fpp_left.get_tmax()))
                return false;
            }
            if (use_right_fpp())
            {
              if (!tol.approximately_equal(t0, fpp_right.get_t0()) || !tol.approximately_equal(tmax, fpp_right.get_tmax()))
                return false;
            }

            // since highest continuity can only be C2 then done
            return true;
          }

        private:
          curve_type f;
          curve_type fp_left, fp_right;
          curve_type fpp_left, fpp_right;
          unsigned int conditions;
          connection_continuity continuity;
      };
    }
  }
}
#endif
