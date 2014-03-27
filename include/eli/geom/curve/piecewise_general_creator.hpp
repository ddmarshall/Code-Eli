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

#ifndef eli_geom_curve_piecewise_general_creator_hpp
#define eli_geom_curve_piecewise_general_creator_hpp

#include <vector>

#include "eli/geom/general/continuity.hpp"

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
      class piecewise_general_creator : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;

          enum
          {
            POINT_SET    =0x000001, // must be set for valid joint
            LEFT_FP_SET  =0x000010,
            RIGHT_FP_SET =0x000100,
            LEFT_FPP_SET =0x001000,
            RIGHT_FPP_SET=0x010000
          };

          enum joint_continuity
          {
            C0=general::C0,
            C1=general::C1,
            C2=general::C2
          };

          class joint_data
          {
            public:
              joint_data() : conditions(0), continuity(C0)
              {
              }

              joint_data(const joint_data &jd)
                : f(jd.f), fp_left(jd.fp_left), fp_right(jd.fp_right), fpp_left(jd.fpp_left),
                  fpp_right(jd.fpp_right), conditions(jd.conditions), continuity(jd.continuity)
              {
              }

              ~joint_data() {}

              const joint_data & operator=(const joint_data &jd)
              {
                if (this!=&jd)
                {
                  f=jd.f;
                  fp_left=jd.fp_left;
                  fp_right=jd.fp_right;
                  fpp_left=jd.fpp_left;
                  fpp_right=jd.fpp_right;
                  conditions=jd.conditions;
                  continuity=jd.continuity;
                }

                return (*this);
              }

              bool operator==(const joint_data &jd) const
              {
                tolerance_type tol;

                if (conditions!=jd.conditions)
                  return false;
                if (continuity!=jd.continuity)
                  return false;
                if (jd.use_f())
                {
                  if (!tol.approximately_equal(f, jd.f))
                    return false;
                }
                if (jd.use_left_fp())
                {
                  if (!tol.approximately_equal(fp_left, jd.fp_left))
                    return false;
                }
                if (jd.use_right_fp())
                {
                  if (!tol.approximately_equal(fp_right, jd.fp_right))
                    return false;
                }
                if (jd.use_left_fpp())
                {
                  if (!tol.approximately_equal(fpp_left, jd.fpp_left))
                    return false;
                }
                if (jd.use_right_fpp())
                {
                  if (!tol.approximately_equal(fpp_right, jd.fpp_right))
                    return false;
                }

                return true;
              }

              bool operator!=(const joint_data &jd) const {return !operator==(jd);}

              // point interface
              bool set_f(const point_type &p)
              {
                f=p;
                conditions|=POINT_SET;
                return check_state();
              }
              point_type get_f() const {return f;}
              bool unset_f()
              {
                conditions&=~POINT_SET;
                return check_state();
              }
              bool use_f() const
              {
                return ((conditions & POINT_SET) == POINT_SET);
              }

              // first derivative interface
              bool set_left_fp(const point_type &fpl)
              {
                fp_left=fpl;
                conditions|=LEFT_FP_SET;
                return check_state();
              }
              bool set_right_fp(const point_type &fpr)
              {
                fp_right=fpr;
                conditions|=RIGHT_FP_SET;
                return check_state();
              }
              bool set_fp(const point_type &p)
              {
                conditions|=(LEFT_FP_SET | RIGHT_FP_SET);
                fp_left=p;
                fp_right=p;
                return check_state();
              }

              point_type get_left_fp() const
              {
                return fp_left;
              }
              point_type get_right_fp() const
              {
                return fp_right;
              }
              void get_fp(point_type &fpl, point_type &fpr) const
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
              bool set_left_fpp(const point_type &fppl)
              {
                fpp_left=fppl;
                conditions|=LEFT_FPP_SET;
                return check_state();
              }
              bool set_right_fpp(const point_type &fppr)
              {
                fpp_right=fppr;
                conditions|=RIGHT_FPP_SET;
                return check_state();
              }
              bool set_fpp(const point_type &p)
              {
                conditions|=(LEFT_FPP_SET | RIGHT_FPP_SET);
                fpp_left=p;
                fpp_right=p;
                return check_state();
              }

              point_type get_left_fpp() const
              {
                return fpp_left;
              }
              point_type get_right_fpp() const
              {
                return fpp_right;
              }
              void get_fpp(point_type &fppl, point_type &fppr) const
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

              bool set_continuity(joint_continuity jc)
              {
                continuity=jc;
                return check_state();
              }
              joint_continuity get_continuity() const
              {
                return continuity;
              }

              bool check_state() const
              {
                tolerance_type tol;

                // check if point set
                if ((conditions & POINT_SET)==0)
                  return false;

                // if highest continuity is C0 then done
                if (continuity==C0)
                  return true;

                // check first derivatives match on both sides
                if ((conditions & (LEFT_FP_SET | RIGHT_FP_SET)) == (LEFT_FP_SET | RIGHT_FP_SET))
                {
                  if (!tol.approximately_equal(fp_left, fp_right))
                    return false;
                }
                // check neither side first derivatives set set
                else if ((conditions & (LEFT_FP_SET | RIGHT_FP_SET)) != 0)
                {
                  return false;
                }

                // if highest continuity is C1 then done
                if (continuity==C1)
                  return true;

                // check second derivatives match on both sides
                if ((conditions & (LEFT_FPP_SET | RIGHT_FPP_SET)) == (LEFT_FPP_SET | RIGHT_FPP_SET))
                {
                  if (!tol.approximately_equal(fpp_left, fpp_right))
                    return false;
                }
                // check neither side second derivatives set set
                else if ((conditions & (LEFT_FPP_SET | RIGHT_FPP_SET)) != 0)
                {
                  return false;
                }

                // since highest continuity can only be C2 then done
                return true;
              }

            private:
              point_type f;
              point_type fp_left, fp_right;
              point_type fpp_left, fpp_right;
              unsigned int conditions;
              joint_continuity continuity;
          };

        public:
          piecewise_general_creator()
            : piecewise_creator_base<data_type, dim__, tolerance_type>(1, 0), joints(2),
              joint_cont(1), max_degree(1), closed(false)
          {
          }
          piecewise_general_creator(const index_type &ns)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(ns, 0), joints(ns+1),
              joint_cont(ns), max_degree(ns), closed_cont(C0), closed(false)
          {
          }
          piecewise_general_creator(const piecewise_general_creator<data_type, dim__, tolerance_type> &pcc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(pcc), joints(pcc.joints),
              joint_cont(pcc.joint_cont), max_degree(pcc.max_degree()), closed(pcc.closed)
          {
          }
          virtual ~piecewise_general_creator()
          {
          };

          void set_closed_continuity(const joint_continuity &cnt)
          {
            closed=true;
            closed_cont=cnt;
          }
          void set_open() {closed=false;}
          bool is_closed() const {return closed;}
          bool is_open() const {return !closed;}


          bool set_conditions(const std::vector<joint_data> &jnts, const std::vector<joint_continuity> &cnt,
                              const std::vector<index_type> &maxd, bool cl=false)
          {
            index_type nsegs(static_cast<index_type>(maxd.size()));

            // ensure input vectors are correct size
            if (jnts.size()!=(nsegs+1))
              return false;
            if (nsegs!=(cnt.size()+1))
              return false;

            // reset the number of segments
            set_number_segments(nsegs);

            joints=jnts;
            joint_cont=cnt;
            closed=cl;

            return true;
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;
            typedef typename curve_type::control_point_type control_point_type;

            index_type nsegs(this->get_number_segments()), i;
            std::vector<index_type> seg_degree(nsegs);

            pc.clear();
#if 0
            // cycle through each segment to find the minimum degree for each segment
            for (i=0; i<nsegs; ++i)
            {
              seg_degree[i]=1;
              switch (joints[i].fp_state())
              {
                case(NOT_SET):
                {
                  break;
                }
                case(VALUE_SET):
                {
                  seg_degree[i]+=1;
                  break;
                }
                default:
                {
                  // should not get here
                  assert(false);

                  return false;
                  break;
                }
              }
              switch (joints[i].fpp_state())
              {
                case(NOT_SET):
                {
                  break;
                }
                case(VALUE_SET):
                {
                  seg_degree[i]+=1;
                  break;
                }
                default:
                {
                  // should not get here
                  assert(false);

                  return false;
                  break;
                }
              }

              if (seg_degree[i]>max_degree[i])
              {
                std::cerr << "Required degree for segment " << i << " is " << seg_degree[i]
                          << " but maximum requested degree is only " << max_degree[i] << std::endl;
                return false;
              }
            }
#endif
            return false;
          }

        private:
          std::vector<joint_data> joints;
          std::vector<joint_continuity> joint_cont;
          std::vector<index_type> max_degree;
          joint_continuity closed_cont;
          bool closed;
      };
    }
  }
}
#endif
