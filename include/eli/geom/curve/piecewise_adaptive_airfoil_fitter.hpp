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

#ifndef piecewise_adaptive_airfoil_fitter_hpp
#define piecewise_adaptive_airfoil_fitter_hpp

#include <vector>
#include <list>
#include <iterator>
#include <algorithm>

#include "eli/code_eli.hpp"

#include "eli/mutil/opt/least_squares.hpp"

#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_creator_base.hpp"
#include "eli/geom/curve/piecewise_adaptive_airfoil_fitter.hpp"
#include "eli/geom/curve/piecewise_general_creator.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      template<typename data__, unsigned short dim__, typename tol__>
      class piecewise_adaptive_airfoil_fitter : public piecewise_creator_base<data__, dim__, tol__>
      {
        public:
          enum subdivide_method
          {
            BISECTION=1,
            MAX_ERROR=2
          };

        public:
          typedef piecewise_creator_base<data__, dim__, tol__> base_class_type;
          typedef typename base_class_type::data_type data_type;
          typedef typename base_class_type::point_type point_type;
          typedef typename base_class_type::index_type index_type;
          typedef typename base_class_type::tolerance_type tolerance_type;
          typedef unsigned short dimension_type;

          piecewise_adaptive_airfoil_fitter()
            : piecewise_creator_base<data_type, dim__, tolerance_type>(2, 0),
              upper_pt(0), lower_pt(0), sub_method(MAX_ERROR), fit_tolerance(1e-6), max_degree(7), approx_te_slope(true)
          {
          }
          piecewise_adaptive_airfoil_fitter(const piecewise_adaptive_airfoil_fitter<data_type, dim__, tolerance_type> &pca)
            : piecewise_creator_base<data_type, dim__, tolerance_type>(pca),
              upper_pt(pca.upper_pt), lower_pt(pca.lower_pt), sub_method(pca.sub_method), fit_tolerance(pca.fit_tolerance),
              max_degree(pca.max_degree), approx_te_slope(pca.approx_te_slope)
          {
          }
          ~piecewise_adaptive_airfoil_fitter() {}

          template<typename it_type>
          void get_airfoil_points(it_type itu, it_type itl) const
          {
            std::copy(upper_pt.begin(), upper_pt.end(), itu);
            std::copy(lower_pt.begin(), lower_pt.end(), itl);
          }

          data_type get_fit_tolerance() const
          {
            return fit_tolerance;
          }
          void set_fit_tolerance(const data_type &ft)
          {
            if (ft>0)
            {
              fit_tolerance=ft;
            }
          }

          index_type get_max_degree() const
          {
            return max_degree;
          }
          void set_max_degree(index_type md)
          {
            // note: fixed maximum degree for fitting
            if (md<1)
            {
              max_degree=md;
            }
            else
            {
              assert(false);
            }
          }

          bool approximate_trailing_edge_slope() const
          {
            return approx_te_slope;
          }
          void set_approximate_trailing_edge_slope(bool ates)
          {
            approx_te_slope = ates;
          }

          subdivide_method get_subdivision_method() const
          {
            return sub_method;
          }
          void set_subdivision_method(subdivide_method sm)
          {
            sub_method=sm;
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

          template<typename it_type>
          bool set_conditions(it_type itb, it_type ite, subdivide_method sm, const data_type &ft, index_type maxd, bool ates)
          {
            set_airfoil_points(itb, ite);
            set_subdivision_method(sm);
            set_fit_tolerance(ft);
            set_max_degree(maxd);
            set_approximate_trailing_edge_slope(ates);

            return true;
          }

          template<typename it_type>
          bool set_conditions(it_type itu, index_type nupper, it_type itl, index_type nlower,
                              subdivide_method sm, const data_type &ft, bool ates)
          {
            set_airfoil_points(itu, nupper, itl, nlower);
            set_subdivision_method(sm);
            set_fit_tolerance(ft);
//            set_max_degree(5);
            set_approximate_trailing_edge_slope(ates);

            return true;
          }

#if 0
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
#endif

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            // typedefs for the piecewise curve and curve
            typedef piecewise<bezier, data_type, dim__, tolerance_type> piecewise_curve_type;
            typedef typename piecewise_curve_type::curve_type curve_type;

            // typedefs for airfoil point collection
            typedef typename airfoil_point_collection_type::const_iterator airfoil_collection_const_interator;
            typedef typename airfoil_point_collection_type::const_reverse_iterator airfoil_collection_const_reverse_interator;

            // typedefs for the general creator
            typedef piecewise_general_creator<data_type, dim__, tolerance_type> general_creator_type;
            typedef typename general_creator_type::fit_data general_fit_data_type;
            typedef typename general_creator_type::joint_data general_joint_data_type;
            typedef typename general_creator_type::joint_continuity general_continuity_type;

            // typedefs for the temporary collection of points needed for the fitting
            typedef std::list<point_type> point_collection_type;
            typedef typename point_collection_type::iterator point_collection_iterator;
            typedef typename point_collection_type::const_iterator point_collection_const_iterator;

            // typdefs for the info needed to build the joint information
            typedef std::list<std::pair<data_type, point_collection_const_iterator>> joint_info_collection_type;
            typedef typename joint_info_collection_type::iterator joint_info_collection_type_iterator;

            index_type i, nsegs(2);
            bool rtn_flag;
            std::vector<index_type> max_deg;                 // collection of maximum degrees for each segment
            std::vector<general_fit_data_type> fit_data;     // collection of fitting info passed into creator
            std::vector<general_joint_data_type> joint_data; // collection of joint data passed into creator
            joint_info_collection_type joint_info;           // collection of iterators to points in fit point
                                                             // collection that represent the joints
            general_creator_type gc;                         // creator used to create each attempted piecewise curve
            point_type te_slope_upper, te_slope_lower;       // trailing edge slopes for upper and lower surface
            index_type max_splitting_iteration_count(12);    // maximum number of splitting iteration are allowed
                                                             // worst case is have 2^max times more segments than
                                                             // was originally desired

            // calculate the trailing edge slope
            if (approx_te_slope)
            {
// FIX: Need better way of determining what the rate of parameter change is. Right now it doesn't
//      change as the sections change. Alternatively, could specify dy/dx, but that would alter the
//      fitting matrices.
              index_type il, ilm1;
              data_type dist;

              il=lower_pt.size()-1;
              ilm1=il-1;
              dist=eli::geom::point::distance(lower_pt[il], lower_pt[ilm1]);
              te_slope_lower=-(lower_pt[il]-lower_pt[ilm1])/(4*dist);

              il=upper_pt.size()-1;
              ilm1=il-1;
              dist=eli::geom::point::distance(upper_pt[il], upper_pt[ilm1]);
              te_slope_upper=(upper_pt[il]-upper_pt[ilm1])/(4*dist);
            }

            // combine points into one vector
            point_collection_type pt;

            // add lower surface points in reverse order
            airfoil_collection_const_reverse_interator afc_rit;
            for (afc_rit=lower_pt.crbegin(); afc_rit!=lower_pt.crend(); ++afc_rit)
            {
              pt.push_back(*afc_rit);
            }
            joint_info.push_back(std::make_pair(0, pt.cbegin()));
            joint_info.push_back(std::make_pair(2, std::prev(pt.cend())));

            // add upper surface points
            airfoil_collection_const_interator afc_it;
            afc_it=upper_pt.cbegin();
            for (++afc_it; afc_it!=upper_pt.cend(); ++afc_it)
            {
              pt.push_back(*afc_it);
            }
            joint_info.push_back(std::make_pair(4, std::prev(pt.cend())));

            assert(joint_info.size()==(nsegs+1));

            index_type iter_count(0);
            bool still_splitting(true);
            while (still_splitting & (iter_count<max_splitting_iteration_count))
            {
              // rebuild data structures for general fitting of new segments
              joint_info_collection_type_iterator jic_it;

              nsegs=joint_info.size()-1;

              joint_data.resize(nsegs+1);

              // set the joint info
              jic_it=joint_info.begin();
              joint_data[0].reset();
              joint_data[0].set_continuity(general_continuity_type::C0);
              joint_data[0].set_f(*(jic_it->second));
              if (approx_te_slope)
              {
                joint_data[0].set_right_fp(te_slope_lower);
              }
              for (++jic_it, i=1; i<nsegs; ++i, ++jic_it)
              {
                joint_data[i].reset();
                joint_data[i].set_continuity(general_continuity_type::C2);
                joint_data[i].set_f(*(jic_it->second));
              }
              joint_data[nsegs].reset();
              joint_data[nsegs].set_continuity(general_continuity_type::C0);
              joint_data[nsegs].set_f(*(jic_it->second));
              if (approx_te_slope)
              {
                joint_data[nsegs].set_left_fp(te_slope_upper);
              }

              // set the maximum degree of each segment
              max_deg.resize(nsegs);
              std::fill(max_deg.begin(), max_deg.end(), max_degree);

              // set the fit info
              point_collection_const_iterator pt_it;

              fit_data.resize(nsegs);

              for (jic_it=joint_info.begin(), i=0; jic_it!=std::prev(joint_info.end()); ++jic_it, ++i)
              {
                fit_data[i].reset();
                for (pt_it=std::next(jic_it->second); pt_it!=std::next(jic_it)->second; ++pt_it)
                {
                  fit_data[i].add_point(*pt_it);
                }
              }

              rtn_flag=gc.set_conditions(joint_data, fit_data, max_deg, false);
              if (!rtn_flag)
              {
                assert(false);
                return false;
              }

              // set the times for each segment
              data_type prev_t;

              jic_it=joint_info.begin();
              gc.set_t0(jic_it->first);
              prev_t=jic_it->first;
              for (++jic_it, i=0; i<nsegs; ++i, ++jic_it)
              {
                gc.set_segment_dt(jic_it->first-prev_t, i);
                prev_t=jic_it->first;
              }
              rtn_flag=gc.create(pc, fit_data);
              if (!rtn_flag)
              {
                assert(false);
                return false;
              }

              // check error of each segment and split those that are above max error
              std::vector<bool>  split(nsegs);
              std::vector<index_type> split_index(nsegs);
              std::vector<data_type> t(nsegs);
              for (jic_it=joint_info.begin(), i=0; i<fit_data.size(); ++jic_it, ++i)
              {
                data_type dist2, ft2(fit_tolerance*fit_tolerance), maxdist2(static_cast<data_type>(-1));

                split[i]=false;
                split_index[i]=-1;
                t[i]=-1;
                for (index_type j=0; j<fit_data[i].npoints(); ++j)
                {
                  dist2=(pc.f(fit_data[i].get_parameter(j))-fit_data[i].get_point(j)).squaredNorm();

//                  std::cout << "segment " << i << ": fit_point[" << j << "].t=" << fit_data[i].get_parameter(j)
//                            << "\t f=" << fit_data[i].get_point(j) << "\t dist=" << std::sqrt(dist2) << std::endl;

                  bool done(false);
                  if (dist2>ft2)
                  {
                    split[i]=true;
                    switch(sub_method)
                    {
                      case(MAX_ERROR):
                      {
                        if (dist2>maxdist2)
                        {
                          split_index[i]=j;
                          maxdist2=dist2;
                        }
                        done=false;
                        break;
                      }
                      case(BISECTION):
                      {
                        index_type jj;
                        data_type tn, tnp1, delt, tfit, tfitp1;


                        // find fit point that is first one past the split parameter
                        tn=jic_it->first;
                        tnp1=std::next(jic_it)->first;
                        delt=tnp1-tn;
                        t[i]=0.5*(tn+tnp1);
                        for (jj=0; jj<fit_data[i].npoints(); ++jj)
                        {
                            if (fit_data[i].get_parameter(jj)>t[i])
                            {
                            if (jj>0)
                            {
                              tfit=fit_data[i].get_parameter(jj-1);
                            }
                            else
                            {
                              tfit=tn;
                            }
                            if (jj==(fit_data[i].npoints()-1))
                            {
                              tfitp1=tnp1;
                            }
                            else
                            {
                              tfitp1=fit_data[i].get_parameter(jj);
                            }

                            break;
                          }
                        }
                        split_index[i]=jj;

                        // make previous fit point a joint
                        if (std::abs((tfit-t[i])/delt)<0.1)
                        {
                          split_index[i]-=1;
                        }
                        // make fit point a joint
                        else if (std::abs((tfitp1-t[i])/delt)<0.1)
                        {
                        }
                        // add point to list
                        else
                        {
                          pt_it=jic_it->second;
                          std::advance(pt_it, split_index[i]+1); // +1 because need to account for joint point that is not counted in fitting
                          pt.insert(pt_it, pc.f(t[i]));
                        }

                        done=true;
                        break;
                      }
                      default:
                      {
                        assert(false);
                        return false;
                      }
                    }
                  }
                  if (done)
                  {
                    // rebuild the joint information
                    break;
                  }
                }
                std::cout << "iteration " << std::setw(2) << iter_count << " section " << std::setw(2) << i;
                if (split[i])
                {
                  switch(sub_method)
                  {
                    case(MAX_ERROR):
                    {
                      std::cout << " index of max error is " << std::setw(3) << split_index[i] << " of " << std::sqrt(maxdist2);
                      break;
                    }
                    case(BISECTION):
                    {
                      std::cout << " index of new bisection point is " << std::setw(3) << split_index[i];
                      break;
                    }
                    default:
                    {
                      assert(false);
                      return false;
                    }
                  }
                  std::cout << std::endl;
                }
                else
                {
                  std::cout << " meets tolerances!" << std::endl;
                }
              }

              // if need to split then do it
              still_splitting=false;
              for (jic_it=joint_info.begin(), i=0; jic_it!=std::prev(joint_info.end()); ++jic_it, ++i)
              {
                if (split[i])
                {
//TODO: Can we consolidate this for all cases to use same code?
                  switch(sub_method)
                  {
                    case(MAX_ERROR):
                    {
                      pt_it=jic_it->second;
                      std::advance(pt_it, split_index[i]+1); // +1 because split_index counts from first fit data, but need to include joint here
                      joint_info.insert(std::next(jic_it), std::make_pair(fit_data[i].get_parameter(split_index[i]), pt_it));
                      ++jic_it;
                      still_splitting=true;
                      break;
                    }
                    case(BISECTION):
                    {
                      pt_it=jic_it->second;
                      std::advance(pt_it, split_index[i]+1); // +1 because split_index counts from first fit data, but need to include joint here
                      joint_info.insert(std::next(jic_it), std::make_pair(t[i], pt_it));
                      ++jic_it;
                      still_splitting=true;
                      break;
                    }
                    default:
                    {
                      assert(false);
                      return false;
                    }
                  }
                }
              }

              ++iter_count;
            }


            if (typeid(data_type)==typeid(float))
            {
              std::cout.flush();
              eli::test::octave_start(1);
              // print out the upper surface points
              for (size_t n=0; n<upper_pt.size(); ++n)
              {
                std::string name("upt"); name+=std::to_string(n);
                eli::test::octave_print(1, upper_pt[n], name);
              }
              // print out the lower surface points
              for (size_t n=0; n<lower_pt.size(); ++n)
              {
                std::string name("lpt"); name+=std::to_string(n);
                eli::test::octave_print(1, lower_pt[n], name);
              }
              eli::test::octave_print(1, pc, "af", false, true);
              eli::test::octave_finish(1, true);
            }

            return false;
          }

        private:
          typedef std::vector<point_type, Eigen::aligned_allocator<point_type>> airfoil_point_collection_type;

          airfoil_point_collection_type upper_pt;
          airfoil_point_collection_type lower_pt;
          subdivide_method sub_method;
          data_type fit_tolerance;
          index_type max_degree;
          bool approx_te_slope;
      };
    }
  }
}
#endif
