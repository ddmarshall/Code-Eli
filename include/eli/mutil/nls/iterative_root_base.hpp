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

#ifndef eli_mutil_nls_iterative_root_base_hpp
#define eli_mutil_nls_iterative_root_base_hpp

#include <utility> // std::pair

namespace eli
{
  namespace mutil
  {
    namespace nls
    {
      template<typename data__>
      class convergence_tester
      {
        public:
          typedef data__ data_type;

          enum activity_state
          {
            inactive       = 0,
            equal          = 1,
            less           = 2,
            less_equal     = 3,
            greater        = 4,
            greater_equal  = 5
          };

        public:
          convergence_tester() : satisfy_both_flag(false)
          {
            rel_tol_info=std::make_pair(static_cast<data__>(0), inactive);
            abs_tol_info=std::make_pair(static_cast<data__>(0), inactive);
          }

          convergence_tester(const convergence_tester<data__> &ct)
            : satisfy_both_flag(ct.satisfy_both_flag), rel_tol_info(ct.rel_tol_info), abs_tol_info(ct.abs_tol_info)
          {
          }

          virtual ~convergence_tester()
          {
          }

          void satisfy_both(bool sb)
          {
            satisfy_both_flag=sb;

            // need to make sure that both tolerances are activated
            if (abs_tol_info.second==inactive)
              set_absolute_tolerance_info(get_absolute_tolerance());
            if (rel_tol_info.second==inactive)
              set_relative_tolerance_info(get_relative_tolerance());
          }

          bool satisfy_both() const
          {
            return satisfy_both_flag;
          }

          void set_relative_tolerance_info(const data__ &rt, const activity_state &as = less)
          {
            rel_tol_info=std::make_pair(rt, as);
          }

          data__ get_relative_tolerance() const
          {
            return rel_tol_info.first;
          }

          const activity_state & get_relative_tolerance_state() const
          {
            return rel_tol_info.second;
          }

          void set_absolute_tolerance_info(const data__ &at, const activity_state &as = less)
          {
            abs_tol_info=std::make_pair(at, as);
          }

          data__ get_absolute_tolerance() const
          {
            return abs_tol_info.first;
          }

          const activity_state & get_absolute_tolerance_state() const
          {
            return abs_tol_info.second;
          }

          bool converged(const data__ &rel_val, const data__ &abs_val) const
          {
            bool rel_state(compare(rel_val, rel_tol_info)), abs_state(compare(abs_val, abs_tol_info));

            if (satisfy_both_flag)
              return rel_state && abs_state;
            else
              return rel_state || abs_state;
          }

        private:
          bool compare(const data__ &val, const std::pair<data__, activity_state> &ref) const
          {
            // test the relative tolerance terms
            switch(ref.second)
            {
              case(inactive):
                  return false;
              break;
              case(equal):
                  return (val==ref.first);
              break;
              case(less):
                  return (val<ref.first);
              break;
              case(less_equal):
                  return (val<=ref.first);
              break;
              case(greater):
                  return (val>ref.first);
              break;
              case(greater_equal):
                  return (val>=ref.first);
              break;
              default:
                assert(false);
                return false;
            }
          }

        private:
          bool satisfy_both_flag;
          std::pair<data__, activity_state> rel_tol_info;
          std::pair<data__, activity_state> abs_tol_info;
      };

      template<typename data__>
      class iterative_root_base
      {
        public:
          enum status
          {
            converged      = 0,
            max_iteration  = 1,
            no_root_found  = 2,
            hit_constraint = 3
          };
          typedef convergence_tester<data__> error_tolerance_type;
          typedef convergence_tester<size_t> max_iteration_type;
          typedef typename error_tolerance_type::data_type tolerance_type;
          typedef typename max_iteration_type::data_type iteration_type;

        public:
          iterative_root_base() : itcnt(0)
          {
            conv.set_absolute_tolerance_info(static_cast<tolerance_type>(1e-8), error_tolerance_type::less);
            conv.set_relative_tolerance_info(static_cast<tolerance_type>(0),    error_tolerance_type::inactive);
            itmax.set_absolute_tolerance_info(static_cast<iteration_type>(200), max_iteration_type::greater);
            itmax.set_relative_tolerance_info(static_cast<iteration_type>(0),   max_iteration_type::inactive);
          }

          iterative_root_base(const iterative_root_base<data__> &irb)
          : conv(irb.conv), itmax(irb.itmax), itcnt(irb.itcnt)
          {
          }

          virtual ~iterative_root_base()
          {
          }

          void set_relative_tolerance(const tolerance_type &rel_tol)
          {
            if (rel_tol<=0)
              conv.set_relative_tolerance_info(0, error_tolerance_type::inactive);
            else
              conv.set_relative_tolerance_info(rel_tol, error_tolerance_type::less);
          }

          tolerance_type get_relative_tolerance() const
          {
            conv.get_relative_tolerance();
          }

          void set_absolute_tolerance(const tolerance_type &abs_tol)
          {
            if (abs_tol<=0)
              conv.set_absolute_tolerance_info(0, error_tolerance_type::inactive);
            else
              conv.set_absolute_tolerance_info(abs_tol, error_tolerance_type::less);
          }

          tolerance_type get_absolute_tolerance() const
          {
            return conv.get_absolute_tolerance();
          }

          void set_max_iteration(const iteration_type &mi)
          {
            if (mi==0)
              itmax.set_absolute_tolerance_info(0, max_iteration_type::inactive);
            else
              itmax.set_absolute_tolerance_info(mi, max_iteration_type::greater);
          }

          iteration_type get_max_iteration() const
          {
            return itmax.get_absolute_tolerance();
          }

          void enforce_both_tolerance(bool ebt)
          {
            conv.satisfy_both(ebt);
          }

          bool enforce_both_tolerance() const
          {
            return conv.satisfy_both();
          }

          void disable_relative_tolerance()
          {
            set_relative_tolerance(-1);
          }

          void disable_absolute_tolerance()
          {
            set_absolute_tolerance(-1);
          }

          const error_tolerance_type & get_tolerance_tester() const
          {
            return conv;
          }

          const max_iteration_type & get_iteration_count_tester() const
          {
            return itmax;
          }

          const iteration_type & get_iteration_count() const
          {
            return itcnt;
          }

        protected:
          bool test_converged(const iteration_type &it, const tolerance_type &relv, const tolerance_type &absv) const
          {
            return conv.converged(relv, absv) || max_iteration_reached(it);
          }

          bool max_iteration_reached(const iteration_type &it) const
          {
            itcnt = it;
            return itmax.converged(0, it);
          }

        private:
          error_tolerance_type conv;
          max_iteration_type itmax;
          mutable iteration_type itcnt;
      };
    }
  }
}
#endif
