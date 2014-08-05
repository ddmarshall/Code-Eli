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

#include "eli/code_eli.hpp"

#include "eli/geom/general/continuity.hpp"

#include "eli/geom/utility/bezier.hpp"

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

                // if already wanting C1 continuous then set right as well
                if (get_continuity()>C0)
                {
                  fp_right=fpl;
                  conditions|=RIGHT_FP_SET;
                }

                return check_state();
              }
              bool set_right_fp(const point_type &fpr)
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

                // if already wanting C2 continuous then set right as well
                if (get_continuity()>C1)
                {
                  fpp_right=fppl;
                  conditions|=RIGHT_FPP_SET;
                }

                return check_state();
              }
              bool set_right_fpp(const point_type &fppr)
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
              max_degree(1), closed(false)
          {
          }
          piecewise_general_creator(const index_type &ns)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(ns, 0), joints(ns+1),
              max_degree(ns), closed(false)
          {
          }
          piecewise_general_creator(const piecewise_general_creator<data_type, dim__, tolerance_type> &pcc)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(pcc), joints(pcc.joints),
              max_degree(pcc.max_degree()), closed(pcc.closed)
          {
          }
          virtual ~piecewise_general_creator()
          {
          }

          void set_closed() {closed=true;}
          void set_open() {closed=false;}
          bool is_closed() const {return closed;}
          bool is_open() const {return !closed;}

          bool set_conditions(const std::vector<joint_data> &jnts, const std::vector<index_type> &maxd, bool cl=false)
          {
            index_type i, nsegs(static_cast<index_type>(maxd.size())), njnts(jnts.size());

            // ensure input vectors are correct size
            if (!cl && (njnts!=(nsegs+1)))
              return false;
            if (cl && (njnts!=nsegs))
              return false;

            // check to make sure have valid end conditions
            if (!cl)
            {
              if (jnts[0].use_left_fp() || jnts[0].use_left_fpp() || jnts[0].get_continuity()!=C0)
              {
                return false;
              }
              if (jnts[nsegs].use_right_fp() || jnts[nsegs].use_right_fpp() || jnts[nsegs].get_continuity()!=C0)
              {
                return false;
              }
            }

            // make sure joints are in valid state
            for (i=0; i<njnts; ++i)
            {
              if (!jnts[i].check_state())
                return false;
            }

            // reset the number of segments
            this->set_number_segments(nsegs);

            joints=jnts;
            max_degree=maxd;
            closed=cl;

            return true;
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            index_type nsegs(this->get_number_segments()), i;
            std::vector<index_type> seg_degree(nsegs);
            std::vector<joint_data> joint_states(joints);

            pc.clear();

            // fix: cannot handle closed curves
            assert(!closed);

            // cycle through each segment to find the minimum degree for each segment
            for (i=0; i<nsegs; ++i)
            {
              seg_degree[i]=1;
              if (joint_states[i].use_right_fp())
                seg_degree[i]+=1;
              if (joint_states[i+1].use_left_fp())
                seg_degree[i]+=1;
              if (joint_states[i].use_right_fpp())
                seg_degree[i]+=1;
              if (joint_states[i+1].use_left_fpp())
                seg_degree[i]+=1;

              // check against maximum degree
              if (!valid_degree(seg_degree[i], max_degree[i]))
              {
                return false;
              }
            }

            // cycle through each joint to find final degrees
            for (i=0; i<=nsegs; ++i)
            {
              switch (joint_states[i].get_continuity())
              {
                // don't need to do anything for C0 case
                case (C0):
                {
                  break;
                }
                // increase degrees by at most two
                case (C2):
                // increase degrees by at most one
                case (C1):
                {
                  // handle C1 continuity
                  if (joint_states[i].use_left_fp())
                  {
                    // only need to increase degree if right is not in use
                    if (!joint_states[i].use_right_fp())
                    {
                      joint_states[i].set_right_fp(joint_states[i].get_left_fp());
                      // NOTE: Here is a place that needs to be modified if closed curve
                      seg_degree[i]+=1;
                    }
                  }
                  else
                  {
                    if (joint_states[i].use_right_fp())
                    {
                      joint_states[i].set_left_fp(joint_states[i].get_right_fp());
                      // NOTE: Here is a place that needs to be modified if closed curve
                      seg_degree[i-1]+=1;
                    }
                    else
                    {
                      bool hit_max_im1(max_degree[i-1]>0 && seg_degree[i-1]>=max_degree[i-1]);
                      bool hit_max_i(max_degree[i]>0 && seg_degree[i]>=max_degree[i]);

                      // raise the degree of the lower degreed segment
                      if ((i==0) & (!closed))
                      {
                        seg_degree[i]+=1;
                      }
                      else if ((i==nsegs) && (!closed))
                      {
                        seg_degree[i-1]+=1;
                      }
                      // NOTE: Here is a place that needs to be modified if closed curve
                      else if ( hit_max_i || (!hit_max_im1 && (seg_degree[i-1]<=seg_degree[i])) )
                      {
                        seg_degree[i-1]+=1;
                      }
                      else
                      {
                        seg_degree[i]+=1;
                      }
                    }
                  }

                  if (joint_states[i].get_continuity()==C1)
                    break;

                  // handle C2 continuity
                  if (joint_states[i].use_left_fpp())
                  {
                    // only need to increase degree if right is not in use
                    if (!joint_states[i].use_right_fpp())
                    {
                      joint_states[i].set_right_fpp(joint_states[i].get_left_fpp());
                      // NOTE: Here is a place that needs to be modified if closed curve
                      seg_degree[i]+=1;
                    }
                  }
                  else
                  {
                    if (joint_states[i].use_right_fpp())
                    {
                      joint_states[i].set_left_fpp(joint_states[i].get_right_fpp());
                      // NOTE: Here is a place that needs to be modified if closed curve
                      seg_degree[i-1]+=1;
                    }
                    else
                    {
                      bool hit_max_im1(max_degree[i-1]>0 && seg_degree[i-1]>=max_degree[i-1]);
                      bool hit_max_i(max_degree[i]>0 && seg_degree[i]>=max_degree[i]);

                      // raise the degree of the lower degreed segment
                      if ((i==0) && (!closed))
                      {
                        seg_degree[i]+=1;
                      }
                      else if ((i==nsegs) && (!closed))
                      {
                        seg_degree[i-1]+=1;
                      }
                      // NOTE: Here is a place that needs to be modified if closed curve
                      if ( hit_max_i || (!hit_max_im1 && (seg_degree[i-1]<=seg_degree[i])) )
                      {
                        seg_degree[i-1]+=1;
//                        assert(max_degree[i-1]<=0 || (seg_degree[i-1]<=max_degree[i-1]));
                      }
                      else
                      {
                        seg_degree[i]+=1;
//                        assert(max_degree[i]<=0 || (seg_degree[i]<=max_degree[i]));
                      }
                    }
                  }

                  break;
                }
                // all cases should be handled by above
                default:
                {
                  assert(false);
                  break;
                }
              }
            }

            // final check of maximum degree and determine number of unknowns
            std::vector<index_type> seg_ind(nsegs+1);

            seg_ind[0]=0;
            for (i=0; i<nsegs; ++i)
            {
              if (!valid_degree(seg_degree[i], max_degree[i]))
              {
                return false;
              }

              seg_ind[i+1]=seg_ind[i]+seg_degree[i]+1;
            }

            // build segments based on joint information
            Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> coef(seg_ind[nsegs]*dim__, seg_ind[nsegs]*dim__), rows(dim__, seg_ind[nsegs]*dim__);
            Eigen::Matrix<data_type, Eigen::Dynamic, 1> x(seg_ind[nsegs]*dim__, 1), rhs(seg_ind[nsegs]*dim__, 1), rhs_seg(dim__, 1);
            index_type cond_no(0);

            for (i=0; i<nsegs; ++i)
            {
              // set the end point conditions
              assert(cond_no<coef.rows());
              set_point_condition(rows, rhs_seg, seg_ind[i], seg_degree[i], joint_states[i].get_f(), true);
              coef.block(cond_no*dim__, 0, dim__, coef.cols())=rows;
              rhs.block(cond_no*dim__, 0, dim__, 1)=rhs_seg;
              ++cond_no;

              assert(cond_no<coef.rows());
              set_point_condition(rows, rhs_seg, seg_ind[i], seg_degree[i], joint_states[i+1].get_f(), false);
              coef.block(cond_no*dim__, 0, dim__, coef.cols())=rows;
              rhs.block(cond_no*dim__, 0, dim__, 1)=rhs_seg;
              ++cond_no;

              // set the end point 1st derivative conditions
              if (joint_states[i].use_right_fp())
              {
                assert(cond_no<coef.rows());
                set_fp_condition(rows, rhs_seg, seg_ind[i], seg_degree[i], joint_states[i].get_right_fp(), this->get_segment_dt(i), true);
                coef.block(cond_no*dim__, 0, dim__, coef.cols())=rows;
                rhs.block(cond_no*dim__, 0, dim__, 1)=rhs_seg;
                ++cond_no;
              }

              if (joint_states[i+1].use_left_fp())
              {
                assert(cond_no<coef.rows());
                set_fp_condition(rows, rhs_seg, seg_ind[i], seg_degree[i], joint_states[i+1].get_left_fp(), this->get_segment_dt(i), false);
                coef.block(cond_no*dim__, 0, dim__, coef.cols())=rows;
                rhs.block(cond_no*dim__, 0, dim__, 1)=rhs_seg;
                ++cond_no;
              }

              // set the end point 2nd derivative conditions
              if (joint_states[i].use_right_fpp())
              {
                assert(cond_no<coef.rows());
                set_fpp_condition(rows, rhs_seg, seg_ind[i], seg_degree[i], joint_states[i].get_right_fpp(), this->get_segment_dt(i), true);
                coef.block(cond_no*dim__, 0, dim__, coef.cols())=rows;
                rhs.block(cond_no*dim__, 0, dim__, 1)=rhs_seg;
                ++cond_no;
              }

              if (joint_states[i+1].use_left_fpp())
              {
                assert(cond_no<coef.rows());
                set_fpp_condition(rows, rhs_seg, seg_ind[i], seg_degree[i], joint_states[i+1].get_left_fpp(), this->get_segment_dt(i), false);
                coef.block(cond_no*dim__, 0, dim__, coef.cols())=rows;
                rhs.block(cond_no*dim__, 0, dim__, 1)=rhs_seg;
                ++cond_no;
              }
            }

            // cycle through each interior joint to set derivative continuity conditions
            // FIX: This needs to be changed for closed curves
            for (i=0; i<nsegs; ++i)
            {
              // set the 1st derivative continuous without specifying value
              if ( (joint_states[i].get_continuity()>C0) && !joint_states[i].use_left_fp() && !joint_states[i].use_right_fp() )
              {
                assert(cond_no<coef.rows());
                set_fp_continuous_condition(rows, rhs_seg, seg_ind[i-1], seg_degree[i-1], seg_degree[i], this->get_segment_dt(i-1), this->get_segment_dt(i));
                coef.block(cond_no*dim__, 0, dim__, coef.cols())=rows;
                rhs.block(cond_no*dim__, 0, dim__, 1)=rhs_seg;
                ++cond_no;
              }

              // set the 2nd derivative continuous without specifying value
              if ( (joint_states[i].get_continuity()>C1) && !joint_states[i].use_left_fpp() && !joint_states[i].use_right_fpp() )
              {
                assert(cond_no<coef.rows());
                set_fpp_continuous_condition(rows, rhs_seg, seg_ind[i-1], seg_degree[i-1], seg_degree[i], this->get_segment_dt(i-1), this->get_segment_dt(i));
                coef.block(cond_no*dim__, 0, dim__, coef.cols())=rows;
                rhs.block(cond_no*dim__, 0, dim__, 1)=rhs_seg;
                ++cond_no;
              }
            }

            assert(cond_no*dim__==coef.rows());

            // solve for the control points
            x=coef.lu().solve(rhs);

            // extract them into control points for each segment
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;
            typedef typename piecewise_curve_type::error_code error_code;
            typedef typename curve_type::control_point_type control_point_type;

            // set curve t0
            pc.set_t0(this->get_t0());

            for (i=0; i<nsegs; ++i)
            {
              curve_type c(seg_degree[i]);
              control_point_type cp;
              error_code err;

              for (index_type j=0; j<=seg_degree[i]; ++j)
              {
                cp=x.block((seg_ind[i]+j)*dim__, 0, dim__, 1).transpose();
                c.set_control_point(cp, j);
              }

              err=pc.push_back(c, this->get_segment_dt(i));
              if (err!=piecewise_curve_type::NO_ERRORS)
              {
                pc.clear();
                pc.set_t0(0);
                return false;
              }
            }

            return true;
          }

        protected:
          static bool valid_degree(const index_type &deg, const index_type &max_deg)
          {
            if (max_deg<=0)
              return true;
            if (deg<=max_deg)
              return true;

            return false;
          }

          template<typename Derived1, typename Derived2>
          void set_point_condition(Eigen::MatrixBase<Derived1> &rows, Eigen::MatrixBase<Derived2> &rhs,
                                   const index_type start_index, const index_type &seg_degree,
                                   const point_type &p, bool segment_start) const
          {
            // set terms
            index_type ind;

            if (segment_start)
            {
              ind=start_index*dim__;
            }
            else
            {
              ind=(start_index+seg_degree)*dim__;
            }

            rows.setConstant(0);
            rows.block(0, ind, dim__, dim__).setIdentity();
            rhs=p.transpose();
          }

          template<typename Derived1, typename Derived2>
          void set_fp_condition(Eigen::MatrixBase<Derived1> &rows, Eigen::MatrixBase<Derived2> &rhs,
                                const index_type start_index, const index_type &seg_degree,
                                const point_type &fp, const data_type &dt, bool segment_start) const
          {
            assert(seg_degree>1);

            // set terms
            index_type ind;
            Eigen::Matrix<data_type, dim__, dim__> coef;

            // initialize rows & rhs to zero
            rhs.setConstant(0);
            rows.setConstant(0);

            if (segment_start)
            {
              ind=start_index*dim__;
            }
            else
            {
              ind=(start_index+seg_degree-1)*dim__;
            }

            coef.setIdentity();
            coef*=seg_degree/dt;
            rows.setConstant(0);
            rows.block(0, ind, dim__, dim__)=-coef;
            rows.block(0, ind+dim__, dim__, dim__)=coef;
            rhs=fp.transpose();
          }

          template<typename Derived1, typename Derived2>
          void set_fp_continuous_condition(Eigen::MatrixBase<Derived1> &rows, Eigen::MatrixBase<Derived2> &rhs,
                                           const index_type start_index, const index_type &l_seg_degree,
                                           const index_type &r_seg_degree, const data_type &l_dt,
                                           const data_type &r_dt) const
          {
            assert( ((l_seg_degree>1) && (r_seg_degree>=1)) || ((r_seg_degree>1) && (l_seg_degree>=1)) );

            // set terms
            index_type l_ind, r_ind;
            Eigen::Matrix<data_type, dim__, dim__> coef;

            // initialize rows & rhs to zero
            rhs.setConstant(0);
            rows.setConstant(0);

            l_ind=(start_index+l_seg_degree-1)*dim__;
            r_ind=(start_index+l_seg_degree+1)*dim__;

            coef.setIdentity();
            coef*=l_seg_degree/l_dt;
            rows.setConstant(0);
            rows.block(0, l_ind, dim__, dim__)=-coef;
            rows.block(0, l_ind+dim__, dim__, dim__)=coef;
            coef.setIdentity();
            coef*=r_seg_degree/r_dt;
            rows.block(0, r_ind, dim__, dim__)=coef;
            rows.block(0, r_ind+dim__, dim__, dim__)=-coef;
          }

          template<typename Derived1, typename Derived2>
          void set_fpp_condition(Eigen::MatrixBase<Derived1> &rows, Eigen::MatrixBase<Derived2> &rhs,
                                 const index_type start_index, const index_type &seg_degree,
                                 const point_type &fpp, const data_type &dt, bool segment_start) const
          {
            assert(seg_degree>1);

            // set terms
            index_type ind;
            Eigen::Matrix<data_type, dim__, dim__> coef;

            // initialize rows & rhs to zero
            rhs.setConstant(0);
            rows.setConstant(0);

            if (segment_start)
            {
              ind=start_index*dim__;
            }
            else
            {
              ind=(start_index+seg_degree-2)*dim__;
            }

            coef.setIdentity();
            coef*=seg_degree*(seg_degree-1)/dt/dt;
            rows.setConstant(0);
            rows.block(0, ind, dim__, dim__)=coef;
            rows.block(0, ind+dim__, dim__, dim__)=-2*coef;
            rows.block(0, ind+2*dim__, dim__, dim__)=coef;
            rhs=fpp.transpose();
          }

          template<typename Derived1, typename Derived2>
          void set_fpp_continuous_condition(Eigen::MatrixBase<Derived1> &rows, Eigen::MatrixBase<Derived2> &rhs,
                                            const index_type &start_index, const index_type &l_seg_degree,
                                            const index_type &r_seg_degree, const data_type &l_dt,
                                            const data_type &r_dt) const
          {
            assert( ((l_seg_degree>1) && (r_seg_degree>=1)) || ((r_seg_degree>1) && (l_seg_degree>=1)) );

            // set terms
            index_type l_ind, r_ind;
            Eigen::Matrix<data_type, dim__, dim__> coef;

            // initialize rows & rhs to zero
            rhs.setConstant(0);
            rows.setConstant(0);

            l_ind=(start_index+l_seg_degree-2)*dim__;
            r_ind=(start_index+l_seg_degree+1)*dim__;

            // only set this if have enough control points to calculate the second derivative
            if (l_seg_degree>1)
            {
              coef.setIdentity();
              coef*=l_seg_degree*(l_seg_degree-1)/l_dt/l_dt;
              rows.block(0, l_ind, dim__, dim__)=coef;
              rows.block(0, l_ind+dim__, dim__, dim__)=-2*coef;
              rows.block(0, l_ind+2*dim__, dim__, dim__)=coef;
            }

            // only set this if have enough control points to calculate the second derivative
            if (r_seg_degree>1)
            {
              coef.setIdentity();
              coef*=r_seg_degree*(r_seg_degree-1)/r_dt/r_dt;
              rows.block(0, r_ind, dim__, dim__)=-coef;
              rows.block(0, r_ind+dim__, dim__, dim__)=2*coef;
              rows.block(0, r_ind+2*dim__, dim__, dim__)=-coef;
            }
          }

        private:
          std::vector<joint_data> joints;
          std::vector<index_type> max_degree;
          bool closed;
      };
    }
  }
}
#endif
