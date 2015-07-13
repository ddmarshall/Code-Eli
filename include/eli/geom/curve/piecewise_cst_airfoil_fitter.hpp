/*********************************************************************************
* Copyright (c) 2015 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef piecewise_cst_airfoil_fitter_hpp
#define piecewise_cst_airfoil_fitter_hpp

#include <vector>
#include <iterator>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/opt/least_squares.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/piecewise_cst_airfoil_creator.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_cst_airfoil_fitter : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef unsigned short dimension_type;
          typedef eli::geom::curve::pseudo::cst_airfoil<data_type> cst_airfoil_type;

          piecewise_cst_airfoil_fitter()
            : piecewise_creator_base<data_type, dim__, tolerance_type>(2, 0),
              upper_pt(0), upper_degree(0), lower_pt(0), lower_degree(0), approx_te_slope(true)
          {
          }
          piecewise_cst_airfoil_fitter(const piecewise_cst_airfoil_fitter<data_type, dim__, tolerance_type> &pca)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(pca),
              upper_pt(pca.upper_pt), upper_degree(pca.upper_degree), lower_pt(pca.lower_pt), lower_degree(pca.lower_degree),
              approx_te_slope(pca.approx_te_slope)
          {
          }
          ~piecewise_cst_airfoil_fitter() {}

          template<typename it_type>
          void get_airfoil_points(it_type itu, it_type itl) const
          {
            std::copy(upper_pt.begin(), upper_pt.end(), itu);
            std::copy(lower_pt.begin(), lower_pt.end(), itl);
          }

          bool approximate_trailing_edge_slope() const
          {
            return approx_te_slope;
          }
          void set_approximate_trailing_edge_slope(bool ates)
          {
            approx_te_slope = ates;
          }

          // * This method assumes that the leading edge is located on the line y=0.
          // * This method assumes that the points are specified from lower trailing edge to upper.
          // * This method assumes that the airfoil points are in the x-y plane, all other
          //   dimensions will be set to zero.
          //   trailing edge with x monotonically decreasing to leading edge and montonically
          //   increasing to upper trailing edge
          // * This method works best when (1) the leading edge is at x=0 & x=1, respectively and
          //   (2) the leading edge is a point in the collection of airfoil points. If (1) is not
          //   the case, the points will be scales and translated to get the leading edge and
          //   trailing in the correct place. This will most likely result in the resulting airfoil
          //   having a different thickness; the airfoil will certainly have a different chord and
          //   x- and y-coordinates. If (2) is not the case, the leading edge point will be found
          //   by fitting a circle to the three airfoil points closest to the leading edge (with at
          //   least one point being on the upper and lower surface. This leading edge point will
          //   be inserted into the collection of fitting points as the actual leading edge
          //   location.
          // * The airfoil points are processed before being stored, so it is not possible to
          //   retrieve the raw points back from the class.
#if 0
          template<typename it_type>
          void set_airfoil_points(it_type itp, index_type npts)
          {
            assert(false);

            it_type it, le_lower(ite), le_upper(ite);
            data_type lower_xmax, upper_xmax, xmin;

            // find leading edge
            {
              it_type le_third(ite);

              // find leading edge
              it=itb;
              lower_xmax=it.x();
              it=ite; --it;
              upper_xmax=it.x();
              for (it=itb; it!=ite; ++it)
              {
                if (it->x()<=xmin)
                {
                  xmin=it->x();
                  le_lower=it;
                }
                else
                {
                  // if haven't set leading edge upper then can find it now
                  if (le_upper==ite)
                  {
                    // make sure that lower surface leading edge is valid
                    if (le_lower->y()==0)
                    {
                      le_upper=le_lower;
                      break;
                    }
                    else
                    {
                      if (le_lower->y()>0)
                      {
                        // set the upper and lower iterators
                        le_upper=le_lower;
                        le_lower=le_upper; --le_lower;
                        if ((le_lower->y()>=0) || (le_upper->y()<=0))
                        {
                          // invalid points: cannot find leading edge
                          assert(false);
                          return;
                        }

                        // find the third point to use in circle fitting to find leading edge
                        it_type ittemp(le_lower); --le_lower;

                        if (ittemp->x() <= le_upper->x())
                        {
                          le_third=ittemp;
                          break;
                        }
                        else
                        {
                          data_type xlow(ittemp->x());

                          ittemp=le_upper; ++le_upper;
                          if (le_upper->x() < xlow())
                          {
                            le_third=le_upper;
                          }
                          else
                          {
                            le_third=le_lower;
                          }
                        }
                      }
                    }
                  }
                }
              }

              // find leading edge, if needed
              if (le_lower!=le_upper)
              {
              }
            }

            // make sure found leading edge
            if ( (le_lower==ite) || (le_upper==ite) )
            {
              assert(false);
              return;
            }

            // resize point collections and fill them
            // if upper and lower leading edge points are different, then there is no point at the
            // leading edge and need to insert one.
            {
              size_t i, nupper, nlower;
              bool add_le_point(le_lower!=le_upper);

              nlower=std::distance(itb, le_lower)+1;
              nupper=std::distance(le_upper, ite);
              if (add_le_point)
              {
                ++nlower;
                ++nupper;
              }

              lower_pt.resize(nlower);
              for (it=itb, i=0; i<nlower; ++it, ++i)
              {
                // set point, scale & translate as needed
                lower_pt[nlower-i-1]=(*it);
              }
              lower_pt[0].SetZero();

              if (add_le_point)
              {
                i=1;
              }
              else
              {
                i=0;
              }
              upper_pt.resize(nupper);
              for (it=le_upper; i<nupper; ++it, ++i)
              {
                upper_pt[i]=(*it);
              }
              if (add_le_point)
              {
                upper_pt[0].SetZero();
              }
            }

            // determine the trailing edge spacings
            upper_dzi=upper_pt.back().y();
            lower_dzi=lower_pt.back().y();
          }
#endif

          template<typename it_type>
          void set_airfoil_points(it_type itu, index_type nupper, it_type itl, index_type nlower)
          {
            it_type ite;

            // fill upper points
            ite = itu; std::advance(ite, nupper);
            upper_pt.resize(static_cast<size_t>(nupper));
            std::copy(itu, ite, upper_pt.begin());

            // fill lower points
            ite = itl; std::advance(ite, nlower);
            lower_pt.resize(static_cast<size_t>(nlower));
            std::copy(itl, ite, lower_pt.begin());
          }

          index_type get_upper_shape_degree() const
          {
            return upper_degree;
          }
          void set_upper_shape_degree(index_type udeg)
          {
            if (udeg>0)
            {
              upper_degree=udeg;
            }
          }

          index_type get_lower_shape_degree() const
          {
            return lower_degree;
          }
          void set_lower_shape_degree(index_type ldeg)
          {
            if (ldeg>0)
            {
              lower_degree=ldeg;
            }
          }

          template<typename it_type>
          bool set_conditions(it_type itb, it_type ite, index_type udeg, index_type ldeg, bool ates)
          {
            set_airfoil_points(itb, ite);
            set_upper_shape_degree(udeg);
            set_lower_shape_degree(ldeg);
            set_approximate_trailing_edge_slope(ates);

            return true;
          }

          template<typename it_type>
          bool set_conditions(it_type itu, index_type nupper, index_type udeg,
                              it_type itl, index_type nlower, index_type ldeg, bool ates)
          {
            set_airfoil_points(itu, nupper, itl, nlower);
            set_upper_shape_degree(udeg);
            set_lower_shape_degree(ldeg);
            set_approximate_trailing_edge_slope(ates);

            return true;
          }

          // This will create a CST airfoil for the given points. Note that CST airfoils
          // are defined with leading edge (x,y)=(0,0) and trailing edge x=1, so if
          // airfoil points are not in that range, then the results CST airfoil will
          // represent the original airfoil but translated, rotated, & scaled to get to CST's
          // cannonical representation.
          bool create(cst_airfoil_type &cst, point_type &le_pt, data_type &theta, data_type &scale_factor) const
          {
            const int LOWER(0), UPPER(1);
            tolerance_type tol;

            // store the trailing edge points for later use
            point_type te_pt[2], te_pt_ave;

            te_pt[LOWER].resize(1,dim__);
            te_pt[UPPER].resize(1,dim__);
            te_pt[LOWER]=lower_pt.back();
            te_pt[UPPER]=upper_pt.back();
            te_pt_ave=0.5*(te_pt[LOWER]+te_pt[UPPER]);
            le_pt=lower_pt[LOWER];
            if (!tol.approximately_equal(le_pt, upper_pt[0]))
            {
              // upper and lower surface have to have same leading edge point
              assert(false);
              return false;
            }

            // fill sample point matrix
            index_type i, iupper(lower_pt.size()-2), npts(upper_pt.size()-2+lower_pt.size()-2);
            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> S(npts,dim__);

            for (i=1; i<static_cast<index_type>(lower_pt.size())-1; ++i)
            {
              for (index_type j=0; j<dim__; ++j)
              {
                S(i-1, j) = lower_pt[i](j);
              }
            }

            for (i=1; i<static_cast<index_type>(upper_pt.size())-1; ++i)
            {
              for (index_type j=0; j<dim__; ++j)
              {
                S(iupper+i-1, j) = upper_pt[i](j);
              }
            }

            // * Find chord line (leading edge and trailing edge)
            // * translate, rotate, & scale so that points go from (0,0) to x=1
            //   (don't forget te_pt[] which holds the trailing edge
            // Need to store the translation, rotation, and scaling in member variables
            theta=-atan2(te_pt_ave.y()-le_pt.y(), te_pt_ave.x()-le_pt.x());
            scale_factor=(te_pt_ave-le_pt).norm();

            // translate leading edge to origin
            data_type ct(std::cos(theta)), st(std::sin(theta));

            // any theta less than ~0.5730 deg. will be considered zero and scale factors
            // within between 0.999 and 1.001 will be considered one
            if ((std::abs(theta) < 0.01) && (std::abs(scale_factor-1) < 0.001))
            {
              theta=0;
              scale_factor=1;
            }
            else
            {
              point_type pt;

              for (i=0; i<npts; ++i)
              {
                // translate to get leading edge at origin
                S.row(i)-=le_pt;

                // rotate back to chord line along x-axis
                pt = S.row(i);
                S.row(i) << (pt.x()*ct-pt.y()*st), (pt.x()*st+pt.y()*ct), pt.z();

                // scale
                S.row(i) /= scale_factor;

//                std::cout << "pt[" << std::setw(2) << i << "]=" << S.row(i) << std::endl;
              }

              // need to transform trailing edge points
              te_pt[LOWER]-=le_pt;
              te_pt[UPPER]-=le_pt;
              pt=te_pt[LOWER]; te_pt[LOWER] << (pt.x()*ct-pt.y()*st), (pt.x()*st+pt.y()*ct), pt.z();
              pt=te_pt[UPPER]; te_pt[UPPER] << (pt.x()*ct-pt.y()*st), (pt.x()*st+pt.y()*ct), pt.z();
              te_pt[LOWER]/= scale_factor;
              te_pt[UPPER]/= scale_factor;
            }

            // need to return the angle needed to rotate a canonical airfoil to this orientation
            theta=-theta;

            // extract trailing edge thicknesses and slopes for fitting
            data_type dte[2], te_slope[2];

            dte[LOWER]=te_pt[LOWER].y();
            dte[UPPER]=te_pt[UPPER].y();
            te_slope[LOWER]=(te_pt[LOWER].y()-S(iupper-1,1))/(te_pt[LOWER].x()-S(iupper-1,0));
            te_slope[UPPER]=(te_pt[UPPER].y()-S(npts-1  ,1))/(te_pt[UPPER].x()-S(npts-1  ,0));

            // record the indices for the control point vector
            index_type n[2], ile[2], ite[2], ncp;

            n[LOWER] = get_lower_shape_degree();
            n[UPPER] = get_upper_shape_degree();
            ile[LOWER] = 0;
            ile[UPPER] = n[LOWER]+1;
            ite[LOWER] = n[LOWER];
            ite[UPPER] = n[LOWER]+n[UPPER]+1;
            ncp        = ite[UPPER]+1;

            // Transform airfoil y-coordinates to the Shape function form
            for (i=0; i<npts; ++i)
            {
              data_type x;
              index_type itemp((i<iupper)?(LOWER):(UPPER));

              x=S(i,0);
              S(i,1)=(S(i,1)-x*dte[itemp])/(sqrt(x)*(1-x));
            }

            // Fit points to 1-D control points
            // * Make sure that the two Shape functions satisfy the curvature continuity relation P0,u=-P0,l
            // * Make sure that the two Shape functions satisfy the trailing edge slope condition
            index_type neqcond;
            neqcond=(approx_te_slope) ? (3) : (1);

            // create new coefficient matrices and rhs vectors
            Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> A(npts,    ncp);
            Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic> B(neqcond, ncp);

            Eigen::Matrix<data_type, Eigen::Dynamic, 1> b_rhs(npts,    1), cp(ncp, 1);
            Eigen::Matrix<data_type, Eigen::Dynamic, 1> d_rhs(neqcond, 1);

            // fill the 3 equality conditions
            B.setZero();
            d_rhs.setZero();
            B(0,ile[LOWER])=1; B(0,ile[UPPER])=1; d_rhs(0)=0;
            if (approx_te_slope)
            {
              B(1,ite[LOWER])=1; d_rhs(1)=dte[LOWER]-te_slope[LOWER];
              B(2,ite[UPPER])=1; d_rhs(2)=dte[UPPER]-te_slope[UPPER];
            }

            // fill the least squares constraints
            {
              data_type coef, k, tau, t;

              A.setZero();
              b_rhs.setZero();
              for (i=0; i<iupper; ++i)
              {
                coef=1;
                k=1;
                t=S(i,0);
                tau=std::pow(1-t, n[LOWER]);
                A(i,0)=coef*k*tau;
                for (index_type j=1; j<=n[LOWER]; ++j)
                {
                  k*=static_cast<data_type>(n[LOWER]-j+1)/j;
                  tau*=t/(1-t);
                  A(i,j)=coef*k*tau;
                }
                b_rhs(i)=S(i,1);
              }
              for (i=iupper; i<npts; ++i)
              {
                coef=1;
                k=1;
                t=S(i,0);
                tau=std::pow(1-t, n[UPPER]);
                A(i,n[LOWER]+1)=coef*k*tau;
                for (index_type j=1; j<=n[UPPER]; ++j)
                {
                  k*=static_cast<data_type>(n[UPPER]-j+1)/j;
                  tau*=t/(1-t);
                  A(i,j+n[LOWER]+1)=coef*k*tau;
                }
                b_rhs(i)=S(i,1);
              }
            }

            eli::mutil::opt::least_squares_eqcon(cp, A, b_rhs, B, d_rhs);

//            std::cout << "cp=" << cp << std::endl;

            // Create CST airfoil from control points
            cst.resize(n[UPPER], n[LOWER]);
            {
              typedef typename cst_airfoil_type::control_point_type cst_control_point_type;

              cst_control_point_type cp_temp[1];

              for (i=0; i<=cst.lower_degree(); ++i)
              {
                cp_temp[0] << cp(i);
                cst.set_lower_control_point(cp_temp[0], i);
              }
              for (i=0; i<=cst.upper_degree(); ++i)
              {
                cp_temp[0] << cp(i+n[LOWER]+1);
                cst.set_upper_control_point(cp_temp[0], i);
              }

              // set the trailing edge thickness of CST airfoil
              cst.set_trailing_edge_thickness(dte[UPPER], -dte[LOWER]);
            }

            return true;
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            cst_airfoil_type cst;
            point_type le;
            data_type theta, scale_factor;

            // create the CST airfoil
            if (!create(cst, le, theta, scale_factor))
            {
              assert(false);
              return false;
            }

            // build piecewise bezier curve using CST Airfoil creator
            {
              typedef eli::geom::curve::piecewise_cst_airfoil_creator<data__, 3, tolerance_type> airfoil_creator_type;

              airfoil_creator_type pcst;
              bool rtn_flag;

              // create curve
              rtn_flag=pcst.set_conditions(cst);
              if (!rtn_flag)
              {
                return false;
              }
              pcst.set_t0(this->get_t0());
              pcst.set_segment_dt(this->get_segment_dt(0), 0);
              pcst.set_segment_dt(this->get_segment_dt(1), 1);
              rtn_flag=pcst.create(pc);
              if (!rtn_flag)
              {
                return false;
              }
            }

            // scale, rotate, & translate back to original
            // Need to use the translation, rotation, and scaling as member variables. Get
            // them from the get operators.
            typename piecewise<bezier, data_type, dim__, tolerance_type>::rotation_matrix_type rmat;

            rmat.setIdentity();
            rmat(0,0)=std::cos(theta); rmat(0,1)=-std::sin(theta);
            rmat(1,0)=-rmat(0,1);      rmat(1,1)=rmat(0,0);
            pc.scale(scale_factor);
            pc.rotate(rmat);
            pc.translate(le);

            return true;
          }

        private:
          typedef std::vector<point_type> airfoil_point_collection_type;

          airfoil_point_collection_type upper_pt;
          index_type upper_degree;
          airfoil_point_collection_type lower_pt;
          index_type lower_degree;
          bool approx_te_slope;
      };
    }
  }
}
#endif
