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

#ifndef eli_geom_surface_piecewise_creator_base_hpp
#define eli_geom_surface_piecewise_creator_base_hpp

#include <list>
#include <iterator>

#include "eli/code_eli.hpp"

#include "eli/util/tolerance.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_creator.hpp"

#include "eli/geom/surface/piecewise.hpp"
#include "eli/geom/surface/bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_creator_base
      {
        public:
          typedef data__  data_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;
          typedef typename point_type::Index index_type;
          typedef tol__ tolerance_type;

        public:
          piecewise_creator_base(const data_type &uu0, const data_type &vv0)
            : du(0), dv(0), u0(uu0), v0(vv0)
          {
          }
          piecewise_creator_base(const piecewise_creator_base<data_type, dim__, tolerance_type> &pcb)
            : du(pcb.du), dv(pcb.dv), u0(pcb.u0), v0(pcb.v0) {}
          virtual ~piecewise_creator_base() {}

          index_type get_number_u_segments() const
          {
            return static_cast<index_type>(du.size());
          }

          index_type get_number_v_segments() const
          {
            return static_cast<index_type>(dv.size());
          }

          data_type get_u0() const {return u0;}
          data_type get_v0() const {return v0;}

          data_type get_segment_du(const index_type &i) const
          {
            if ((i<0) || (i>=static_cast<index_type>(du.size())))
            {
              assert(false);
              return static_cast<data_type>(-1);
            }

            return du[i];
          }

          data_type get_segment_dv(const index_type &i) const
          {
            if ((i<0) || (i>=get_number_v_segments()))
            {
              assert(false);
              return static_cast<data_type>(-1);
            }

            return dv[i];
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const = 0;

        protected:

          typedef std::vector<data_type> param_container_type;
          typedef typename param_container_type::const_iterator param_iterator;

          param_iterator du_begin() const {return du.begin();}
          param_iterator du_end() const {return du.end();}
          param_iterator dv_begin() const {return dv.begin();}
          param_iterator dv_end() const {return dv.end();}

          void set_number_u_segments(const index_type &ns)
          {
            size_t old_size(du.size());

            du.resize(ns);
            for (index_type i=old_size; i<ns; ++i)
              du[i]=1;
          }
          void set_number_v_segments(const index_type &ns)
          {
            size_t old_size(dv.size());

            dv.resize(ns);
            for (index_type i=old_size; i<ns; ++i)
              dv[i]=1;
          }

          /** insert du before element i
           */
          void insert_du(const data_type &duu, const index_type &i)
          {
            if ((duu>0) && (i>=0) && (i<static_cast<index_type>(du.size())))
            {
              du.insert(du.begin()+i, duu);
            }
            else
            {
              assert(false);
            }
          }
          /** insert dv before element j
           */
          void insert_dv(const data_type &dvv, const index_type &j)
          {
            if ((dvv>0) && (j>=0) && (j<static_cast<index_type>(dv.size())))
            {
              dv.insert(dv.begin()+j, dvv);
            }
            else
            {
              assert(false);
            }
          }

          void append_du(const data_type &duu)
          {
            if (duu>0)
            {
              du.append(duu);
            }
            else
            {
              assert(false);
            }
          }

          void append_dv(const data_type &dvv)
          {
            if (dvv>0)
            {
              dv.append(dvv);
            }
            else
            {
              assert(false);
            }
          }

          void set_initial_u(const data_type &uu0) {u0=uu0;}
          void set_initial_v(const data_type &vv0) {v0=vv0;}

          void set_du(const data_type &duu, const index_type &i)
          {
            if ((duu>0) && (i>=0) && (i<static_cast<index_type>(du.size())))
              du[i]=duu;
            else
              assert(false);
          }
          void set_dv(const data_type &dvv, const index_type &j)
          {
            if ((dvv>0) && (j>=0) && (j<static_cast<index_type>(dv.size())))
              dv[j]=dvv;
            else
              assert(false);
          }

        private:
          param_container_type du, dv;
          data_type u0, v0;
      };
    }
  }
}
#endif
