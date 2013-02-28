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

#ifndef geom_curve_fit_container_hpp
#define geom_curve_fit_container_hpp

#include <cassert>
#include <map>

#include "eli/mutil/fd/d1o2.hpp"
#include "eli/mutil/fd/d2o2.hpp"

#include "eli/geom/point/distance.hpp"
#include "eli/geom/general/continuity.hpp"

#include "Eigen/StdVector"

namespace eli
{
  namespace geom
  {
    namespace curve
    {

      template<typename data_type__, typename index_type__, size_t dim__, size_t con_dim__>
      class fit_container
      {
        public:
          typedef data_type__ data_type;
          typedef index_type__ index_type;
          typedef Eigen::Matrix<data_type, 1, dim__> point_type;

          enum error_code
          {
            NO_ERROR=0,
            INVALID_INDEX=1,
            INDEX_NOT_FOUND=2,
            TOO_FEW_POINTS=3,
            UNKNOWN_ERROR=999
          };

          class constraint_info
          {
            public:
              typedef Eigen::Matrix<data_type, 1, con_dim__> point_type;

              enum use_states
              {
                NOT_USED =-1,
                SET      = 0,
                FD       = 1
              };

            private:
              point_type fp, fpp;
              use_states use_fp, use_fpp;

            public:
              constraint_info() : use_fp(NOT_USED), use_fpp(NOT_USED) {}
              constraint_info(const constraint_info &ci) : fp(ci.fp), fpp(ci.fpp), use_fp(ci.use_fp), use_fpp(ci.use_fpp) {}
              ~constraint_info() {}

              constraint_info & operator=(const constraint_info &ci)
              {
                if (this!=&ci)
                {
                  fp=ci.fp;
                  fpp=ci.fpp;
                  use_fp=ci.use_fp;
                  use_fpp=ci.use_fpp;
                }

                return *this;
              }

              point_type get_fp() const {return fp;}
              point_type get_fpp() const {return fpp;}
              use_states using_fp() const {return use_fp;}
              use_states using_fpp() const {return use_fpp;}

              void set_fp(const point_type &p, bool fd)
              {
                if (fd)
                  use_fp=FD;
                else
                  use_fp=SET;

                fp=p;
              }

              void set_fpp(const point_type &p, bool pfd, const point_type &pp, bool ppfd)
              {
                if (pfd)
                  use_fp=FD;
                else
                  use_fp=SET;
                fp=p;

                if (ppfd)
                  use_fpp=FD;
                else
                  use_fpp=SET;
                fpp=pp;
              }
          };

          typedef typename constraint_info::point_type constraint_point_type;

        private:
          eli::geom::general::continuity end_continuity;
          std::vector<point_type, Eigen::aligned_allocator<point_type> > points;
          typedef std::map<index_type, constraint_info, std::less<index_type>, Eigen::aligned_allocator<std::pair<const index_type, point_type> > > constraint_collection;

          constraint_collection constraints;

        private:
          bool valid_index(const index_type &i) const
          {
            if ( (i>=0) && (i<static_cast<index_type>(number_points())) )
              return true;
            else
              return false;
          }

          void evaluate_fp(constraint_point_type &fp, const index_type &i) const
          {
            // need to do finite difference to calculate 1st derivative vector
            typedef  eli::mutil::fd::d1o2<data_type> fd1_type;
            fd1_type d;
            point_type vec[3];
            data_type t[3];

            if (i==0)
            {
              if (open())
              {
                d.set_stencil(fd1_type::RIGHT);

                // set points
                vec[0]=points[0];
                vec[1]=points[1];
                vec[2]=points[2];
              }
              else
              {
                d.set_stencil(fd1_type::CENTER);

                // set points
                vec[0]=points[number_points()-1];
                vec[1]=points[0];
                vec[2]=points[1];
              }
            }
            else if (i==static_cast<index_type>(number_points())-1)
            {
              if (open())
              {
                d.set_stencil(fd1_type::LEFT);

                // set points
                vec[0]=points[number_points()-3];
                vec[1]=points[number_points()-2];
                vec[2]=points[number_points()-1];
              }
              else
              {
                d.set_stencil(fd1_type::CENTER);

                // set points
                vec[0]=points[number_points()-2];
                vec[1]=points[number_points()-1];
                vec[2]=points[0];
              }
            }
            else
            {
              d.set_stencil(fd1_type::CENTER);

              // set points
              vec[0]=points[i-1];
              vec[1]=points[i];
              vec[2]=points[i+1];
            }

            // calculate the distances between points
            size_t joff;
            if (dim__==con_dim__)
            {
              joff=0;
              t[0]=0;
              eli::geom::point::distance(t[1], vec[1], vec[0]);
              eli::geom::point::distance(t[2], vec[2], vec[1]);
              t[2]+=t[1];
            }
            else
            {
              joff=1;
              t[0]=vec[0](0);
              t[1]=vec[1](0);
              t[2]=vec[2](0);
            }

            // calculate the approximate derivative
            data_type y[3], fpj(0);
            for (size_t j=0; j<con_dim__; ++j)
            {
              y[0]=vec[0](j+joff);
              y[1]=vec[1](j+joff);
              y[2]=vec[2](j+joff);
              d.evaluate(fpj, y, t);
              fp(j)=fpj;
            }
          }

          void evaluate_fpp(constraint_point_type &fpp, const index_type &i) const
          {
            // need to do finite difference to calculate 1st derivative vector
            typedef  eli::mutil::fd::d2o2<data_type> fd2_type;
            fd2_type d;
            point_type vec[4];
            data_type t[4];

            if (i==0)
            {
              if (open())
              {
                d.set_stencil(fd2_type::RIGHT);

                // set points
                vec[0]=points[0];
                vec[1]=points[1];
                vec[2]=points[2];
                vec[3]=points[3];
              }
              else
              {
                d.set_stencil(fd2_type::LEFT_BIASED);

                // set points
                vec[0]=points[number_points()-2];
                vec[1]=points[number_points()-1];
                vec[2]=points[0];
                vec[3]=points[1];
              }
            }
            else if (i==1)
            {
              if (open())
              {
                d.set_stencil(fd2_type::RIGHT_BIASED);

                // set points
                vec[0]=points[0];
                vec[1]=points[1];
                vec[2]=points[2];
                vec[3]=points[3];
              }
              else
              {
                d.set_stencil(fd2_type::LEFT_BIASED);

                // set points
                vec[0]=points[number_points()-1];
                vec[1]=points[0];
                vec[2]=points[1];
                vec[3]=points[2];
              }
            }
            else if (i==static_cast<index_type>(number_points())-2)
            {
              if (open())
              {
                d.set_stencil(fd2_type::LEFT_BIASED);

                // set points
                vec[0]=points[number_points()-4];
                vec[1]=points[number_points()-3];
                vec[2]=points[number_points()-2];
                vec[3]=points[number_points()-1];
              }
              else
              {
                d.set_stencil(fd2_type::RIGHT_BIASED);

                // set points
                vec[0]=points[number_points()-3];
                vec[1]=points[number_points()-2];
                vec[2]=points[number_points()-1];
                vec[3]=points[0];
              }
            }
            else if (i==static_cast<index_type>(number_points())-1)
            {
              if (open())
              {
                d.set_stencil(fd2_type::LEFT);

                // set points
                vec[0]=points[number_points()-4];
                vec[1]=points[number_points()-3];
                vec[2]=points[number_points()-2];
                vec[3]=points[number_points()-1];
              }
              else
              {
                d.set_stencil(fd2_type::RIGHT_BIASED);

                // set points
                vec[0]=points[number_points()-2];
                vec[1]=points[number_points()-1];
                vec[2]=points[0];
                vec[3]=points[1];
              }
            }
            else
            {
              d.set_stencil(fd2_type::RIGHT_BIASED);

              // set points
              vec[0]=points[i-1];
              vec[1]=points[i];
              vec[2]=points[i+1];
              vec[3]=points[i+2];
            }

            // calculate the distances between points
            size_t joff;
            if (dim__==con_dim__)
            {
              joff=0;
              t[0]=0;
              eli::geom::point::distance(t[1], vec[1], vec[0]);
              eli::geom::point::distance(t[2], vec[2], vec[1]);
              t[2]+=t[1];
              eli::geom::point::distance(t[3], vec[3], vec[2]);
              t[3]+=t[2];
            }
            else
            {
              joff=1;
              t[0]=vec[0](0);
              t[1]=vec[1](0);
              t[2]=vec[2](0);
              t[3]=vec[3](0);
            }

            // calculate the approximate derivative
            data_type y[4], fppj(0);
            for (size_t j=0; j<con_dim__; ++j)
            {
              y[0]=vec[0](j+joff);
              y[1]=vec[1](j+joff);
              y[2]=vec[2](j+joff);
              y[3]=vec[3](j+joff);
              d.evaluate(fppj, y, t);
              fpp(j)=fppj;
            }
          }

        public:
          fit_container() : end_continuity(eli::geom::general::NOT_CONNECTED) {}
          ~fit_container() {}

          size_t number_points() const {return points.size();}
          size_t number_constraint_points() const {return constraints.size();}

          size_t number_constraints(bool fit=true) const
          {
            typename constraint_collection::const_iterator it;
            size_t ncon(0), nfit=fit?1:0;

            for (it=constraints.begin(); it!=constraints.end(); ++it)
            {
              if (it->second.using_fp()==constraint_info::NOT_USED)
                ncon+=nfit;
              else if (it->second.using_fpp()==constraint_info::NOT_USED)
                ncon+=nfit+1;
              else
                ncon+=nfit+2;
            }

            return ncon;
          }

          bool open() const {return end_continuity==eli::geom::general::NOT_CONNECTED;}
          bool closed() const {return end_continuity!=eli::geom::general::NOT_CONNECTED;}
          const eli::geom::general::continuity & get_end_flag() const {return end_continuity;}

          void set_end_flag(const eli::geom::general::continuity &op)
          {
            // if changing flag, then need to check start and end constraints to see if the
            // finite differences need to be redone
            if (op!=get_end_flag())
            {
              index_type indx[4];
              typename constraint_collection::const_iterator it;

              // set the open flag
              end_continuity=op;

              // these are the four indexes that need to be checked for
              indx[0]=0;
              indx[1]=1;
              indx[2]=static_cast<index_type>(number_points())-2;
              indx[3]=static_cast<index_type>(number_points())-1;

              // check the start and end constraints need to be recalculated
              for (size_t i=0; i<4; ++i)
              {
                it=constraints.find(indx[i]);
                if (it!=constraints.end())
                {
                  // if first derivative is finite difference, then either both are, or the
                  // second derivative is not set
                  if (it->second.using_fp()==constraint_info::FD)
                  {
                    // check if second derivative is set
                    if (it->second.using_fpp()==constraint_info::FD)
                      add_C2_constraint(indx[i]);
                    else
                      add_C1_constraint(indx[i]);
                  }
                  // if second derivative is finite difference, then must have set first derivative
                  else if (it->second.using_fpp()==constraint_info::FD)
                  {
                    add_C2_constraint(indx[i], it->second.get_fp());
                  }
                }
              }
            }
          }

          template<typename it__>
          void set_points(it__ itb, it__ ite)
          {
            points.clear();
            points.insert(points.begin(), itb, ite);
          }

          template<typename it__>
          void get_points(it__ itb) const
          {
            std::copy(points.begin(), points.end(), itb);
          }

          error_code get_start_constraint(constraint_info &ci) const
          {
            return get_constraint(0, ci);
          }

          error_code get_start_constraint(constraint_info &ci, point_type &pt) const
          {
            return get_constraint(0, ci, pt);
          }

          error_code get_end_constraint(constraint_info &ci) const
          {
            return get_constraint(number_points()-1, ci);
          }

          error_code get_end_constraint(constraint_info &ci, point_type &pt) const
          {
            return get_constraint(number_points()-1, ci, pt);
          }

          error_code get_constraint(const index_type &i, constraint_info &ci) const
          {
            // check if have valid index
            if (!valid_index(i))
              return INVALID_INDEX;

            typename constraint_collection::const_iterator it=constraints.find(i);
            if (it==constraints.end())
              return INDEX_NOT_FOUND;

            ci=it->second;

            return NO_ERROR;
          }

          error_code get_constraint(const index_type &i, constraint_info &ci, point_type &pt) const
          {
            error_code ec;

            ec=get_constraint(i, ci);
            if (ec==NO_ERROR)
              pt=points[i];

            return NO_ERROR;
          }

          error_code add_start_C0_constraint()
          {
            return add_C0_constraint(0);
          }

          error_code add_start_C1_constraint()
          {
            return add_C1_constraint(0);
          }

          error_code add_start_C1_constraint(const constraint_point_type &fp)
          {
            return add_C1_constraint(0, fp);
          }

          error_code add_start_C2_constraint()
          {
            return add_C2_constraint(0);
          }

          error_code add_start_C2_constraint(const constraint_point_type &fp)
          {
            return add_C2_constraint(0, fp);
          }

          error_code add_start_C2_constraint(const constraint_point_type &fp, const constraint_point_type &fpp)
          {
            return add_C2_constraint(0, fp, fpp);
          }

          error_code add_end_C0_constraint()
          {
            return add_C0_constraint(static_cast<int>(number_points())-1);
          }

          error_code add_end_C1_constraint()
          {
            return add_C1_constraint(static_cast<int>(number_points())-1);
          }

          error_code add_end_C1_constraint(const constraint_point_type &fp)
          {
            return add_C1_constraint(static_cast<int>(number_points())-1, fp);
          }

          error_code add_end_C2_constraint()
          {
            return add_C2_constraint(static_cast<int>(number_points())-1);
          }

          error_code add_end_C2_constraint(const constraint_point_type &fp)
          {
            return add_C2_constraint(static_cast<int>(number_points())-1, fp);
          }

          error_code add_end_C2_constraint(const constraint_point_type &fp, const constraint_point_type &fpp)
          {
            return add_C2_constraint(static_cast<int>(number_points())-1, fp, fpp);
          }

          error_code add_C0_constraint(const index_type &i)
          {
            // check if have valid index
            if (!valid_index(i))
              return INVALID_INDEX;

            // this will either create or replace
            constraint_info ci;
            constraints[i]=ci;

            return NO_ERROR;
          }

          error_code add_C1_constraint(const index_type &i)
          {
            // check if have valid index
            if (!valid_index(i))
              return INVALID_INDEX;

            // check if have enough points (3 for finite difference of 1st derivative)
            if (number_points()<3)
              return TOO_FEW_POINTS;

            // need to do finite difference to calculate 1st derivative vector
            constraint_point_type fp;
            evaluate_fp(fp, i);

            // this will either create or replace
            constraint_info ci;
            ci.set_fp(fp, true);
            constraints[i]=ci;

            return NO_ERROR;
          }

          error_code add_C1_constraint(const index_type &i, const constraint_point_type &fp)
          {
            // check if have valid index
            if (!valid_index(i))
              return INVALID_INDEX;

            // this will either create or replace
            constraint_info ci;
            ci.set_fp(fp, false);
            constraints[i]=ci;

            return NO_ERROR;
          }

          error_code add_C2_constraint(const index_type &i)
          {
            // check if have valid index
            if (!valid_index(i))
              return INVALID_INDEX;

            // check if have enough points (4 for finite difference of 2nd derivative)
            if (number_points()<4)
              return TOO_FEW_POINTS;

            // need to do finite difference to calculate 1st derivative vector
            constraint_point_type fp;
            evaluate_fp(fp, i);

            // need to do finite difference to calculate 2nd derivative vector
            constraint_point_type fpp;
            evaluate_fpp(fpp, i);

            // this will either create or replace
            constraint_info ci;
            ci.set_fpp(fp, true, fpp, true);
            constraints[i]=ci;

            return NO_ERROR;
          }

          error_code add_C2_constraint(const index_type &i, const constraint_point_type &fp)
          {
            // check if have valid index
            if (!valid_index(i))
              return INVALID_INDEX;

            // check if have enough points (4 for finite difference of 2nd derivative)
            if (number_points()<4)
              return TOO_FEW_POINTS;

            // need to do finite difference to calculate 2nd derivative vector
            constraint_point_type fpp;
            evaluate_fpp(fpp, i);

            // this will either create or replace
            constraint_info ci;
            ci.set_fpp(fp, false, fpp, true);
            constraints[i]=ci;

            return NO_ERROR;
          }

          error_code add_C2_constraint(const index_type &i, const constraint_point_type &fp, const constraint_point_type &fpp)
          {
            // check if have valid index
            if (!valid_index(i))
              return INVALID_INDEX;

            // this will either create or replace
            constraint_info ci;
            ci.set_fpp(fp, false, fpp, false);
            constraints[i]=ci;

            return NO_ERROR;
          }

          error_code remove_constraint(const index_type &i)
          {
            // check if have valid index
            if (!valid_index(i))
              return INVALID_INDEX;

            typename constraint_collection::iterator it=constraints.find(i);
            if (it==constraints.end())
              return INDEX_NOT_FOUND;

            constraints.erase(it);
            return NO_ERROR;
          }

          template<typename it__>
          error_code get_constraint_indexes(it__ itout) const
          {
            typename constraint_collection::const_iterator it;
            for (it=constraints.begin(); it!=constraints.end(); ++it, ++itout)
              (*itout)=it->first;

            return NO_ERROR;
          }
      };
    }
  }
}

#endif
