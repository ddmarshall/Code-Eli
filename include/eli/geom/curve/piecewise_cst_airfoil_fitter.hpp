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
          bool create(cst_airfoil_type &cst, Eigen::Matrix<data_type, dim__, dim__> &transform_out, point_type &translate_out, data_type &actual_leading_edge_t) const
          {
            const index_type LOWER(0), UPPER(1);
            tolerance_type tol;

            Eigen::Matrix<data_type, Eigen::Dynamic, dim__> fit_point;
            index_type i, n, n_te[2], n_le_des, n_le_act;

            // fill matrix with points
            // * store the index for lower trailing edge, upper trailing edge, and designated leading edge
            // * might need to specify actual leading edge later if points near leading edge are negative
            fit_point.resize(static_cast<index_type>(lower_pt.size()+upper_pt.size()-1), dim__);
            n_te[LOWER]=0;
            n_te[UPPER]=fit_point.rows()-1;
            n_le_des=static_cast<index_type>(lower_pt.size()-1);
            for (n=0, i=lower_pt.size()-1; i>0; --i, ++n)
            {
              fit_point.row(n)=lower_pt[i];
            }
            for (i=0; i<upper_pt.size(); ++i, ++n)
            {
              fit_point.row(n)=upper_pt[i];
            }

            assert(n==fit_point.rows());

            // find transformation parameters
            point_type te_pt_ave;

            te_pt_ave=0.5*(fit_point.row(n_te[LOWER])+fit_point.row(n_te[UPPER]));

            // build transformation matrix and vector and apply
            {
              Eigen::Matrix<data_type, Eigen::Dynamic, dim__> trans_vec(fit_point.rows(), dim__);
              point_type le_pt(fit_point.row(n_le_des));
              data_type scale_factor, theta, ct, st;

              scale_factor=(te_pt_ave-le_pt).norm();
              theta=atan2(te_pt_ave.y()-le_pt.y(), te_pt_ave.x()-le_pt.x());

              // any theta less than ~0.5730 deg. will be considered zero
              if (std::abs(theta) < 0.01)
              {
                theta=0;
              }

              // scale factors within between 0.999 and 1.001 will be considered one
              if (std::abs(scale_factor-1) < 0.001)
              {
                scale_factor=1;
              }
              ct=std::cos(-theta);
              st=std::sin(-theta);


              transform_out << ct/scale_factor, -st/scale_factor, 0,
                               st/scale_factor,  ct/scale_factor, 0,
                               0,                0,               1/scale_factor;
              translate_out=-le_pt*transform_out.transpose();
              trans_vec.col(0)=Eigen::Matrix<data_type, Eigen::Dynamic, 1>::Constant(fit_point.rows(), 1, translate_out.x());
              trans_vec.col(1)=Eigen::Matrix<data_type, Eigen::Dynamic, 1>::Constant(fit_point.rows(), 1, translate_out.y());
              trans_vec.col(2)=Eigen::Matrix<data_type, Eigen::Dynamic, 1>::Constant(fit_point.rows(), 1, translate_out.z());

              fit_point=fit_point*transform_out.transpose()+trans_vec;
            }

//            if (typeid(data_type)==typeid(float))
//            {
//              std::cout.flush();
//              eli::test::octave_start(1);
//              for (i=0; i<fit_point.rows(); ++i)
//              {
//                std::string name("upt"); name+=std::to_string(i);
//                eli::test::octave_print(1, fit_point.row(i), name);
//              }
//
//              eli::test::octave_finish(1, true);
//              std::cout << std::endl << std::endl;
//            }

            // find actual leading edge (point furthest from trailing edge)
            data_type xmin(fit_point.row(0).x());
            bool manual_leading_edge, extra_points_lower;

            n_le_act=0;
            for (n=0; n<static_cast<index_type>(fit_point.rows()-1); ++n)
            {
              if (fit_point.row(n).x()<=xmin)
              {
                xmin=fit_point.row(n).x();
                n_le_act=n;
              }
            }

            // need to translate, rotate, and scale to get actual leading edge at zero
            manual_leading_edge=n_le_act!=n_le_des;
            extra_points_lower=n_le_act>n_le_des;
            actual_leading_edge_t=0;
            if (manual_leading_edge)
            {
              point_type le_act_pt;
              data_type tscale, ttheta;

              te_pt_ave=0.5*(fit_point.row(n_te[LOWER])+fit_point.row(n_te[UPPER]));
              le_act_pt=fit_point.row(n_le_act);
              ttheta=atan2(te_pt_ave.y()-le_act_pt.y(), te_pt_ave.x()-le_act_pt.x());
              tscale=(te_pt_ave-le_act_pt).norm();

              Eigen::Matrix<data_type, dim__, dim__> rot_mat;
              Eigen::Matrix<data_type, Eigen::Dynamic, dim__> trans_vec(fit_point.rows(), dim__);
              data_type ct(std::cos(-ttheta)), st(std::sin(-ttheta));

              rot_mat << ct/tscale, -st/tscale, 0,
                         st/tscale,  ct/tscale, 0,
                         0,                0,   1/tscale;
              trans_vec.col(0)=Eigen::Matrix<data_type, Eigen::Dynamic, 1>::Constant(fit_point.rows(), 1, -le_act_pt.x());
              trans_vec.col(1)=Eigen::Matrix<data_type, Eigen::Dynamic, 1>::Constant(fit_point.rows(), 1, -le_act_pt.y());
              trans_vec.col(2)=Eigen::Matrix<data_type, Eigen::Dynamic, 1>::Constant(fit_point.rows(), 1, -le_act_pt.z());

              fit_point=(fit_point+trans_vec)*rot_mat.transpose();

              transform_out=rot_mat*transform_out;
              translate_out=(translate_out-le_act_pt)*rot_mat.transpose();
              actual_leading_edge_t=((extra_points_lower)?(-1):(1))*fit_point.row(n_le_des).x();
//              std::cout <<   "actual_leading_edge_t=" << actual_leading_edge_t;
//              std::cout << "\tleading_edge_desired =" << fit_point.row(n_le_des);
//              std::cout << "\tleading_edge_desired =" << (fit_point.row(n_le_des)-translate_out)*transform_out.transpose().inverse();
//              std::cout << std::endl;
            }

            // extract trailing edge thicknesses and slopes for fitting
            data_type dte[2], te_slope[2];

            dte[LOWER]=fit_point.row(n_te[LOWER]).y();
            dte[UPPER]=fit_point.row(n_te[UPPER]).y();
            te_slope[LOWER]=(dte[LOWER]-fit_point.row(n_te[LOWER]+1).y())/(fit_point.row(n_te[LOWER]).x()-fit_point.row(n_te[LOWER]+1).x());
            te_slope[UPPER]=(dte[UPPER]-fit_point.row(n_te[UPPER]-1).y())/(fit_point.row(n_te[UPPER]).x()-fit_point.row(n_te[UPPER]-1).x());

//            std::cout << "theta_out=" << theta << " te_slope_upper=" << te_slope[UPPER] << " te_slope_lower=" << te_slope[LOWER] << std::endl;


            // Fit points to 1-D control points
            index_type ncp, deg[2], ile[2], ite[2];
            deg[LOWER] = get_lower_shape_degree();
            deg[UPPER] = get_upper_shape_degree();
            ile[LOWER] = 0;
            ile[UPPER] = deg[LOWER]+1;
            ite[LOWER] = deg[LOWER];
            ite[UPPER] = deg[LOWER]+deg[UPPER]+1;
            ncp        = ite[UPPER]+1;

            // * Make sure that the two Shape functions satisfy the curvature continuity relation P0,u=-P0,l
            // * (Optional) Make sure that the two Shape functions satisfy the trailing edge slope condition
            // * (Optional) Make sure that desired leading edge is on curve if not where actual leading edge is
            index_type neqcond, npts;
            neqcond=(approx_te_slope) ? (3) : (1); // if approximating trailing edge slopes then adds two equality conditions
            npts=fit_point.rows()-3; // do not use lead edge (1) or trailing edges (2)
            if (manual_leading_edge)
            {
              ++neqcond; // one more equality condition for going through desired leading edge point
              --npts;    // one less fitting point because it is an equality condition
            }

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
            if (manual_leading_edge)
            {
              index_type itemp, ind_offset;

              // if desired leading edge is on lower surface
              if (extra_points_lower)
              {
                itemp=LOWER;
                ind_offset=0;
              }
              // else desired leading edge is on uppers surface
              else
              {
                itemp=UPPER;
                ind_offset=deg[LOWER]+1;
              }

              data_type t(fit_point.row(n_le_des).x());
              Eigen::Matrix<data_type, Eigen::Dynamic, 1> co;

              eli::geom::utility::bezier_coefficient_factors(co, t, deg[itemp]);
              for (index_type j=0; j<=deg[itemp]; ++j)
              {
                B(3,j+ind_offset)=co(j);
              }
              d_rhs(3)=(fit_point.row(n_le_des).y()-t*dte[itemp])/(sqrt(t)*(1-t));

//              std::cout << "B.row(3)=" << B.row(3) << std::endl << "d_rhs(3)=" << d_rhs(3) << std::endl;
            }

            // fill the least squares constraints
            {
              A.setZero();
              b_rhs.setZero();
              for (i=1; i<static_cast<index_type>(fit_point.rows()-1); ++i)
              {
                if ((i!=n_le_act) && (i!=n_le_des))
                {
                  index_type itemp, ind_offset, row_offset;

                  if (i<n_le_act)
                  {
                    itemp=LOWER;
                    ind_offset=0;
                    row_offset=1;
                  }
                  else
                  {
                    itemp=UPPER;
                    ind_offset=deg[LOWER]+1;
                    row_offset=2;
                  }
                  if (manual_leading_edge && (i>n_le_des))
                  {
                    ++row_offset;
                  }

                  data_type t(fit_point.row(i).x());
                  Eigen::Matrix<data_type, Eigen::Dynamic, 1> co;

                  eli::geom::utility::bezier_coefficient_factors(co, t, deg[itemp]);
                  for (index_type j=0; j<=deg[itemp]; ++j)
                  {
                    A(i-row_offset,j+ind_offset)=co(j);
                  }

                  // hack for rare cases when still get negative t
                  if (t<=0)
                  {
//                    std::cout << "negative t at point " << i << " :" << fit_point.row(i) << std::endl;
//                    std::cout << "fit points:" << fit_point << std::endl;
                    A.row(i-row_offset)=A.row(i-row_offset-1);
                    b_rhs(i-row_offset)=b_rhs(i-row_offset-1);
                  }
                  else
                  {
                    b_rhs(i-row_offset)=(fit_point.row(i).y()-t*dte[itemp])/(sqrt(t)*(1-t));
//                    std::cout << "i=" << std::setw(3) << i << "\tt=" << t << "\tb_rhs(" << std::setw(3) << i << ")=" << b_rhs(i-row_offset) << std::endl;
                  }
                }
              }
            }

            eli::mutil::opt::least_squares_eqcon(cp, A, b_rhs, B, d_rhs);

//            std::cout << "A=" << A << std::endl;
//            {
//              Eigen::JacobiSVD<Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic>> svd(A);
//              std::cout << "Cond(A)=" << svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1) << std::endl;
//            }
//            std::cout << "b_rhs=" << b_rhs << std::endl;
//            std::cout << "B=" << B << std::endl;
//            {
//              Eigen::JacobiSVD<Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic>> svd(B);
//              std::cout << "Cond(B)=" << svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1) << std::endl;
//            }
//            std::cout << "d_rhs=" << d_rhs << std::endl;
//            std::cout << "cp=" << cp << std::endl;

            // Create CST airfoil from control points
            cst.resize(deg[UPPER], deg[LOWER]);
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
                cp_temp[0] << cp(i+deg[LOWER]+1);
                cst.set_upper_control_point(cp_temp[0], i);
              }

              // set the trailing edge thickness of CST airfoil
              cst.set_trailing_edge_thickness(dte[UPPER], -dte[LOWER]);

//              if (manual_leading_edge && (typeid(data_type)==typeid(float)))
//              {
//                std::cout.flush();
//                eli::test::octave_start(1);
//                // print out the upper surface points
//                for (size_t n=0; n<fit_point.rows(); ++n)
//                {
//                  std::string name("pt"); name+=std::to_string(n);
//                  eli::test::octave_print(1, fit_point.row(n), name);
//                }
//                eli::test::octave_print(1, cst, "cst");
//                eli::test::octave_finish(1, true);
//              }
            }

            return true;
          }

          virtual bool create(piecewise<bezier, data_type, dim__, tolerance_type> &pc) const
          {
            cst_airfoil_type cst;
            point_type translate_out;
            Eigen::Matrix<data_type, dim__, dim__> tranform_out;
            data_type actual_leading_edge_t;

            // create the CST airfoil
            if (!create(cst, tranform_out, translate_out, actual_leading_edge_t))
            {
              assert(false);
              return false;
            }

//            if (typeid(data_type)==typeid(float))
//            {
//              std::cout.flush();
//              eli::test::octave_start(1);
//              // print out the upper surface points
//              eli::test::octave_print(1, cst, "cst");
//              eli::test::octave_finish(1, true);
//            }

            // build piecewise bezier curve using CST Airfoil creator
            {
              typedef eli::geom::curve::piecewise_cst_airfoil_creator<data__, 3, tolerance_type> airfoil_creator_type;

              piecewise<bezier, data_type, dim__, tolerance_type> pc_temp;
              typename piecewise<bezier, data_type, dim__, tolerance_type>::curve_type curve_temp;
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
              rtn_flag=pcst.create(pc_temp);
              if (!rtn_flag)
              {
                return false;
              }

              // split based on desired leading edge & redo indexing
              tolerance_type tol;
              if (!tol.approximately_equal(actual_leading_edge_t,0))
              {
                data_type t_split, tsqrt(std::sqrt(std::abs(actual_leading_edge_t))), dt[3];

                // since actual Bezier parameterization of CST curve is x=t^2, need to take sqrt() of CST t
                if (actual_leading_edge_t<0)
                {
                  t_split=this->get_t0()+(1-tsqrt)*this->get_segment_dt(0);
                  dt[0]=this->get_segment_dt(0);
                  dt[1]=tsqrt*this->get_segment_dt(1);
                  dt[2]=(1-tsqrt)*this->get_segment_dt(1);
                }
                else
                {
                  t_split=this->get_t0()+this->get_segment_dt(0)+tsqrt*this->get_segment_dt(1);
                  dt[0]=(1-tsqrt)*this->get_segment_dt(0);
                  dt[1]=tsqrt*this->get_segment_dt(0);
                  dt[2]=this->get_segment_dt(1);
                }
//                std::cout << "split index=" << t_split << "\tf(t_split)=" << pc_temp.f(t_split) << std::endl;
                pc_temp.split(t_split);

                pc.set_t0(this->get_t0());
                pc_temp.get(curve_temp, 0);
                pc.push_back(curve_temp, dt[0]);
                pc_temp.get(curve_temp, 1);
                pc.push_back(curve_temp, dt[1]);
                pc_temp.get(curve_temp, 2);
                pc.push_back(curve_temp, dt[2]);
              }
              else
              {
                pc=pc_temp;
              }
            }

            // scale, rotate, & translate back to original
            // Need to use the translation, rotation, and scaling as member variables. Get
            // them from the get operators.
            pc.translate(-translate_out);
            pc.rotate(tranform_out.inverse());

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
