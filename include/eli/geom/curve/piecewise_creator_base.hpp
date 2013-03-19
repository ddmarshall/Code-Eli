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

#ifndef eli_geom_curve_piecewise_creator_base_hpp
#define eli_geom_curve_piecewise_creator_base_hpp

#include <vector>

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__>
      class piecewise_creator_base
      {
        private:
          typedef data__  data_type;
          typedef int index_type;

        public:
          piecewise_creator_base(int n, const data__ &tt0) : dt(n), t0(tt0)
          {
            for (index_type i=0; i<static_cast<index_type>(dt.size()); ++i)
              dt[i]=1;
          }
          piecewise_creator_base(const piecewise_creator_base<data__> &pcb) : dt(pcb.dt), t0(pcb.t0) {}

          index_type get_number_segments() const
          {
            return static_cast<index_type>(dt.size());
          }
          void set_number_segments(const index_type &ns)
          {
            dt.resize(ns);
            for (index_type i=0; i<ns; ++i)
              dt[i]=1;

            // tell child classes that number of segments changed
            number_segments_changed();
          }

          void set_t0(const data_type &tt0) {t0=tt0;}
          data_type get_t0() const {return t0;}

          void set_segment_dt(const data_type &dtt, const index_type &i)
          {
            if ((dtt>0) && (i>=0) && (i<static_cast<index_type>(dt.size())))
              dt[i]=dtt;
            else
              assert(false);
          }
          data_type get_segment_dt(const index_type &i) const
          {
            if ((i<0) || (i>=static_cast<index_type>(dt.size())))
            {
              assert(false);
              return static_cast<data_type>(-1);
            }

            return dt[i];
          }

        private:
          virtual void number_segments_changed() {};

        private:
          std::vector<data_type> dt;
          data_type t0;
      };
    }
  }
}
#endif
