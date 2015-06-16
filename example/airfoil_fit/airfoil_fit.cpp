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

#include <cstdlib>   // EXIT_SUCCESS
#include <iostream>  // std::iostream
#include <algorithm> // std::sort, std::find, std::find_if
#include <utility>   // std::pair
#include <iomanip>
#include <string>
#include <vector>
#include <list>

#include "eli/constants/math.hpp"
#include "eli/geom/point/distance.hpp"
#include "eli/geom/general/continuity.hpp"

#include "eli/geom/curve/bezier.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_general_creator.hpp"

// FIX: Need to make this accessible to test and example or come up with another way of presenting
//      curve results.
#include "../test/include/octave_helpers.hpp"

namespace eli
{
  namespace test
  {
    template<typename data__>
    void octave_param_print(int figno, const eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> &pc,
                            const std::string &name="", bool show_joints=false)
    {
      typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
      typedef typename piecewise_curve_type::curve_type curve_type;
      typedef typename piecewise_curve_type::data_type data_type;
      typedef typename piecewise_curve_type::point_type point_type;
      typedef typename piecewise_curve_type::index_type index_type;

      std::string nm, tbuf, xbuf, ybuf, zbuf, xpbuf, ypbuf, zpbuf, xppbuf, yppbuf, zppbuf;

      index_type i, j;

      // FIX: Need to update to show the joints
      assert(!show_joints);

      // build name
      if (name=="")
      {
        nm=random_string(5);
      }
      else
      {
        nm=name;
      }

      // get the parameters for the segments
      const index_type npts(130);
      index_type nsegs(pc.number_segments()), nt;
      data_type dt;
      std::vector<data_type> tseg(nsegs+1);
      std::vector<index_type> npseg(nsegs);
      std::vector<data_type> t;

      pc.get_parameters(tseg.begin());
      dt=(tseg[nsegs]-tseg[0])/npts;

      calculate_t(t, npseg, tseg, dt);
      nt=t.size();

      // set the curve points
      tbuf=nm+"_t=[";
      xbuf=nm+"_curv_x=[";
      ybuf=nm+"_curv_y=[";
      zbuf=nm+"_curv_z=[";
      xpbuf=nm+"_curv_xp=[";
      ypbuf=nm+"_curv_yp=[";
      zpbuf=nm+"_curv_zp=[";
      xppbuf=nm+"_curv_xpp=[";
      yppbuf=nm+"_curv_ypp=[";
      zppbuf=nm+"_curv_zpp=[";
      for (i=0; i<nt; ++i)
      {
        point_type pt;

        tbuf+=std::to_string(t[i]);
        pt=pc.f(t[i]);
        xbuf+=std::to_string(pt.x());
        ybuf+=std::to_string(pt.y());
        zbuf+=std::to_string(pt.z());
        pt=pc.fp(t[i]);
        xpbuf+=std::to_string(pt.x());
        ypbuf+=std::to_string(pt.y());
        zpbuf+=std::to_string(pt.z());
        pt=pc.fpp(t[i]);
        xppbuf+=std::to_string(pt.x());
        yppbuf+=std::to_string(pt.y());
        zppbuf+=std::to_string(pt.z());
        if (i<(nt-1))
        {
          tbuf+=", ";
          xbuf+=", ";
          ybuf+=", ";
          zbuf+=", ";
          xpbuf+=", ";
          ypbuf+=", ";
          zpbuf+=", ";
          xppbuf+=", ";
          yppbuf+=", ";
          zppbuf+=", ";
        }
      }
      tbuf+="];";
      xbuf+="];";
      ybuf+="];";
      zbuf+="];";
      xpbuf+="];";
      ypbuf+="];";
      zpbuf+="];";
      xppbuf+="];";
      yppbuf+="];";
      zppbuf+="];";

      std::cout << "% parameter curve: " << nm << std::endl;
      std::cout << "figure(" << figno << ");" << std::endl;
      std::cout << tbuf << std::endl;
      std::cout << xbuf << std::endl;
      std::cout << ybuf << std::endl;
      std::cout << zbuf << std::endl;
      std::cout << xpbuf << std::endl;
      std::cout << ypbuf << std::endl;
      std::cout << zpbuf << std::endl;
      std::cout << xppbuf << std::endl;
      std::cout << yppbuf << std::endl;
      std::cout << zppbuf << std::endl;

      std::cout << "setenv('GNUTERM', 'x11');" << std::endl;
      std::cout << "subplot(3, 1, 1);" << std::endl;
      std::cout << "plot(" << nm << "_t, " << nm << "_curv_x, " << "'-g', "
                           << nm << "_t, " << nm << "_curv_y, " << "'-b', "
                           << nm << "_t, " << nm << "_curv_z, " << "'-r');" << std::endl;
      std::cout << "xlabel('t');" << std::endl;
      std::cout << "ylabel('x, y, z');" << std::endl;
      std::cout << "subplot(3, 1, 2);" << std::endl;
      std::cout << "plot(" << nm << "_t, " << nm << "_curv_xp, " << "'-g', "
                           << nm << "_t, " << nm << "_curv_yp, " << "'-b', "
                           << nm << "_t, " << nm << "_curv_zp, " << "'-r');" << std::endl;
      std::cout << "xlabel('t');" << std::endl;
      std::cout << "ylabel('x^{\\prime}, y^{\\prime}, z^{\\prime}');" << std::endl;
      std::cout << "subplot(3, 1, 3);" << std::endl;
      std::cout << "plot(" << nm << "_t, " << nm << "_curv_xpp, " << "'-g', "
                           << nm << "_t, " << nm << "_curv_ypp, " << "'-b', "
                           << nm << "_t, " << nm << "_curv_zpp, " << "'-r');" << std::endl;
      std::cout << "xlabel('t');" << std::endl;
      std::cout << "ylabel('x^{\\prime\\prime}, y^{\\prime\\prime}, z^{\\prime\\prime}');" << std::endl;
    }
  }
}

typedef double data__;
typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
typedef piecewise_curve_type::curve_type curve_type;
typedef piecewise_curve_type::point_type point_type;
typedef piecewise_curve_type::data_type data_type;
typedef piecewise_curve_type::index_type index_type;
typedef piecewise_curve_type::tolerance_type tolerance_type;
typedef eli::geom::curve::piecewise_general_creator<data_type, 3, tolerance_type> general_creator_type;
typedef general_creator_type::joint_data joint_data_type;
typedef general_creator_type::joint_continuity continuity_type;

void create_airfoil(std::vector<point_type, Eigen::aligned_allocator<point_type>> &pts)
{
  pts.resize(40);
  // lower surface
  pts[ 0] << 1.0000,-0.0020, 0;
  pts[ 1] << 0.9012,-0.0087, 0;
  pts[ 2] << 0.8078,-0.0247, 0;
  pts[ 3] << 0.7017,-0.0453, 0;
  pts[ 4] << 0.6178,-0.0572, 0;
  pts[ 5] << 0.4123,-0.0684, 0;
  pts[ 6] << 0.3545,-0.0678, 0;
  pts[ 7] << 0.2986,-0.0657, 0;
  pts[ 8] << 0.2453,-0.0621, 0;
  pts[ 9] << 0.1499,-0.0511, 0;
  pts[10] << 0.1029,-0.0441, 0;
  pts[11] << 0.0741,-0.0365, 0;
  pts[12] << 0.0451,-0.0284, 0;
  pts[13] << 0.0143,-0.0112, 0;
  pts[14] << 0.0076,-0.0116, 0;
  pts[15] << 0.0029,-0.0077, 0;
  pts[16] << 0.0000, 0.0000, 0;
  // upper surface
  pts[17] << 0.0013, 0.0030, 0;
  pts[18] << 0.0045, 0.0094, 0;
  pts[19] << 0.0096, 0.0138, 0;
  pts[20] << 0.0165, 0.0183, 0;
  pts[21] << 0.0252, 0.0228, 0;
  pts[22] << 0.0477, 0.0315, 0;
  pts[23] << 0.0767, 0.0397, 0;
  pts[24] << 0.1118, 0.0473, 0;
  pts[25] << 0.1524, 0.0539, 0;
  pts[26] << 0.1980, 0.0593, 0;
  pts[27] << 0.2979, 0.0636, 0;
  pts[28] << 0.3015, 0.0665, 0;
  pts[29] << 0.3578, 0.0680, 0;
  pts[30] << 0.4160, 0.0680, 0;
  pts[31] << 0.4455, 0.0675, 0;
  pts[32] << 0.5049, 0.0652, 0;
  pts[33] << 0.5930, 0.0585, 0;
  pts[34] << 0.6501, 0.0514, 0;
  pts[35] << 0.7050, 0.0416, 0;
  pts[36] << 0.7623, 0.0297, 0;
  pts[37] << 0.8168, 0.0221, 0;
  pts[38] << 0.9074, 0.0108, 0;
  pts[39] << 1.0000, 0.0050, 0;
}

int main(int /*argc*/, char * /*argv*/[])
{
  typedef std::pair<index_type, continuity_type> joint_type;

  enum airfoil_trailing_edge_cap
  {
    no_te_cap = 0,
    blunt_te_cap = 1,
    force_chord_blunt_te_cap = 2,
    force_thick_blunt_te_cap = 3,
    round_te_cap = 11,
    force_chord_round_te_cap = 12,
    force_thick_round_te_cap = 13,
    sharp_te_cap = 21,
    force_chord_sharp_te_cap = 22,
    force_thick_sharp_te_cap = 23
  };

  index_type i, j;
  tolerance_type tol;

  //
  // User supplied inputs
  //
  std::vector<point_type, Eigen::aligned_allocator<point_type>> pt;
  index_type npts;

  // user supplied airfoil points
  create_airfoil(pt);
  npts=static_cast<index_type>(pt.size());

  // user supplied force joint and force not joint
  std::vector<index_type> not_joint(0);
  std::vector<joint_type> is_joint(0);
  not_joint.push_back(28);
  not_joint.push_back(29);
//  is_joint[0]=std::make_pair(32, continuity_type::C2);

  // user supplied trailing edge treatment
  airfoil_trailing_edge_cap te_cap(blunt_te_cap);
  assert( (te_cap==no_te_cap) || (te_cap==blunt_te_cap) );

  //
  // Closed trailing edge processing
  //
  point_type pt0(pt[0]), ptN(pt[npts-1]), pt_te;
  bool te_points_open;

  // determine if the trailing edge points are open
  te_points_open=!tol.approximately_equal(pt0, ptN);

  // create closed trailing point
  // if non-blunt is desired, then after creating airfoil replace the first and last segments (the
  //   two trailing edge segments) with the desired shape)
  // TODO: This does not account for forced at particular chord location or particular thickness
  if (te_points_open && te_cap!=no_te_cap)
  {
    if (tol.approximately_equal(pt0.y(), 0))
    {
      pt.push_back(pt0);
    }
    else if (tol.approximately_equal(ptN.y(), 0))
    {
      pt.insert(pt.begin(), ptN);
    }
    else if (pt0.y()*ptN.y()<0)
    {
      data_type xval;

      xval=pt0.x()-pt0.y()*((ptN.x()-pt0.x())/(ptN.y()-pt0.y()));
      pt_te << xval, static_cast<data_type>(0), static_cast<data_type>(0);
      pt.insert(pt.begin(), pt_te);
      pt.push_back(pt_te);
    }
    else
    {
      pt_te=0.5*(pt0+ptN);
      pt.insert(pt.begin(), pt_te);
      pt.push_back(pt_te);
    }
  }
  npts=pt.size();

  //
  // Process points for joints
  //

  // build running length for points
  // find leading edge
  // find any breaks in airfoil and mark them as joints
  std::vector<data_type> length_point(npts);
  std::vector<data_type> angle_point(npts);
  std::vector<data_type> radius_point(npts);
  std::vector<data_type> radius_rate_point(npts);

  // calculate the running length of each point
  length_point[0]=static_cast<data_type>(0);
  angle_point[0]=static_cast<data_type>(1);
  radius_point[0]=static_cast<data_type>(1);
  radius_rate_point[0]=static_cast<data_type>(0);
  for (i=1; i<npts-1; ++i)
  {
    // calculate the length
    length_point[i]=length_point[i-1]+eli::geom::point::distance(pt[i-1], pt[i]);

    // calculate the panel angle
    {
      point_type p_i, p_im1;

      p_i=pt[i+1]-pt[i];
      p_im1=pt[i]-pt[i-1];
      p_i.normalize();
      p_im1.normalize();

      angle_point[i]=p_i.dot(p_im1);
    }

    // calculate the point radius and radius rate of change
    {
      data_type xa, xb, xc, xab, xac, xbc, ya, yb, yc, yab, yac, ybc;

      xa=pt[i-1].x();
      xb=pt[i  ].x();
      xc=pt[i+1].x();
      ya=pt[i-1].y();
      yb=pt[i  ].y();
      yc=pt[i+1].y();
      xab=xa-xb;
      xac=xa-xc;
      xbc=xb-xc;
      yab=ya-yb;
      yac=ya-yc;
      ybc=yb-yc;

      radius_point[i]=0.5*std::sqrt((xab*xab+yab*yab)*(xac*xac+yac*yac)*(xbc*xbc+ybc*ybc))/(xc*yab-xb*yac+xa*ybc);

      if (i>1)
      {
        data_type drds_f((radius_point[i]-radius_point[i-1])/(length_point[i]-length_point[i-1])), drds_b((radius_point[i-1]-radius_point[i-2])/(length_point[i-1]-length_point[i-2]));
        radius_rate_point[i-1]=(length_point[i-1]-length_point[i-2])*drds_f+(length_point[i]-length_point[i-1])*drds_b;
      }
    }
  }
  length_point[i]=length_point[i-1]+eli::geom::point::distance(pt[i-1], pt[i]);
  angle_point[i]=1;
  radius_point[0]=radius_point[1];
  radius_point[i]=radius_point[i-1];
  radius_rate_point[i]=static_cast<data_type>(0);
  {
    data_type drds_f((radius_point[i]-radius_point[i-1])/(length_point[i]-length_point[i-1])), drds_b((radius_point[i-1]-radius_point[i-2])/(length_point[i-1]-length_point[i-2]));
    radius_rate_point[i-1]=(length_point[i-1]-length_point[i-2])*drds_f+(length_point[i]-length_point[i-1])*drds_b;
  }

  // set the first and last point as a C0 joint
  typedef std::list<joint_type> joint_collection_type;
  joint_collection_type joint_point;
  joint_collection_type::iterator jit;
  index_type ile(-1);
  data_type dist_to_origin(static_cast<data_type>(1000));
  const data_type le_close(static_cast<data_type>(1e-4)), le_near(static_cast<data_type>(0.01));
  const data_type angle_point_max(static_cast<data_type>(0.8192)), angle_point_abs_max(0.2);
  const data_type radius_rate_max(static_cast<data_type>(1));

  joint_point.push_back(std::make_pair(0, continuity_type::C0));
  joint_point.push_back(std::make_pair(npts-1, continuity_type::C0));
  for (i=1; i<npts-1; ++i)
  {
    // look for leading edge
    if (pt[i].x()<le_close)
    {
      point_type orig; orig << 0, 0, 0;
      data_type temp_dist(eli::geom::point::distance(orig, pt[i]));

      // if distance from current point to leading edge location is less than tolerance
      // then have candidate for leading edge
      if ( (temp_dist<le_close) && (temp_dist<dist_to_origin) )
      {
        dist_to_origin=temp_dist;
        ile=i;
      }
    }

    // joint if angle turn is larger than 35 deg -> cos(35 deg)=0.8192 && not close to leading edge
    if ( (angle_point[i]<angle_point_abs_max)
      || ((angle_point[i]<angle_point_max) && (std::abs(radius_rate_point[i])>radius_rate_max) && (pt[i].x()>le_near)) )
    {
      joint_point.push_back(std::make_pair(i, continuity_type::C0));
    }
  }

  // if have a leading edge identified
  if (ile>0)
  {
    joint_collection_type::iterator it=std::find(joint_point.begin(), joint_point.end(), std::make_pair(ile, continuity_type::C0));

    if (it==joint_point.end())
      joint_point.push_back(std::make_pair(ile, continuity_type::C2));
    else
      it->second=continuity_type::C2;
  }

  // go through user input vectors to unset and then set joints
  for (i=0; i<not_joint.size(); ++i)
  {
    auto comp_if = [idx=not_joint[i]](const std::pair<index_type, continuity_type> &a) -> bool
    {
      return (a.first == idx);
    };

    if ( (not_joint[i]>0) && (not_joint[i]<(npts-1)) )
    {
      joint_point.erase(std::remove_if(joint_point.begin(), joint_point.end(), comp_if));
    }
  }
  for (i=0; i<is_joint.size(); ++i)
  {
    auto comp_if = [idx=is_joint[i].first](const joint_type &a) -> bool
    {
      return (a.first == idx);
    };

    if ( (is_joint[i].first>0) && (is_joint[i].first<(npts-1)) )
    {
      joint_collection_type::iterator it=std::find_if(joint_point.begin(), joint_point.end(), comp_if);

      if (it==joint_point.end())
        joint_point.push_back(is_joint[i]);
      else
        it->second=is_joint[i].second;
    }
  }

  // sort final joint list
  auto comp_lt = [](const joint_type &a, const joint_type &b) -> bool
  {
    return (a.first < b.first);
  };
  joint_point.sort(comp_lt);

  // write out info for debugging
//  jit=joint_point.begin();
//  for (i=0; i<npts; ++i)
//  {
//    std::cout << "pt[" << std::setw(2) << i << "] @ " << std::setw(9) << length_point[i];
//    std::cout << "  angle=" << std::setw(9) << (static_cast<data_type>(180)/eli::constants::math<data_type>::pi())*std::acos(angle_point[i]);
//    std::cout << "  radius=" << std::setw(11) << radius_point[i];
//    std::cout << "  d(radius)/ds=" << std::setw(11) << radius_rate_point[i];
//    if (i==jit->first)
//    {
//      std::cout << " and is joint " << std::setw(2) << j;
//      switch (jit->second)
//      {
//        case (continuity_type::C0):
//          std::cout << " type C0";
//          break;
//        case (continuity_type::C1):
//          std::cout << " type C1";
//          break;
//        case (continuity_type::C2):
//          std::cout << " type C2";
//          break;
//        default:
//          std::cout << "How did I get here?";
//          break;
//      }
//      ++jit;
//    }
//    std::cout << std::endl;
//  }
//
//  std::cout << "We are done here" << std::endl;

  // if capping trailing edge, need account for that in the parameterization
  // Capping options:
  //    * None - do nothing
  //    * Blunt - vertical curve if points specify an open trailing edge
  //    * Round - rounded curve if points specify an open trailing edge
  //    * Sharp - extend trailing edge panel until intersect if points specify an open trailing edge
  //    * Force @ chord - Blunt, Round, or Sharp at specified chord location (extend trailing edge panels if needed)
  //    * Force thickness - Blunt, Round, or Sharp at chord location that yields desired thickness
  // These options will determine where the first and last airfoil point needs to be located

  // should have:
  //   * list of joints (at least lower trailing edge [& mid-chord???], leading edge, and [upper mid-chord??? &] trailing edge)
  //   * list of points to fit
  //   * trailing edge points in correct position for capping
  // build airfoil piecewise curve

  // build and replace trailing edge cap as needed

  // construct the piecewise curve
  index_type nsegs(static_cast<index_type>(joint_point.size()-1));
  std::vector<typename general_creator_type::joint_data> joint(nsegs+1);
  std::vector<typename general_creator_type::index_type> max_degree(nsegs);
  std::vector<typename general_creator_type::fit_data> fit_point(nsegs);
  std::vector<data_type> t(nsegs+1);
  const index_type max_degree_param(4);

  // fill each joint and time
  const data_type t_scale(4/length_point.back());

  jit=joint_point.begin();
  for (i=0; i<=nsegs; ++i, ++jit)
  {
    t[i]=t_scale*length_point[jit->first];
    joint[i].set_f(pt[jit->first]);
    if (jit->second > eli::geom::general::continuity::C0)
    {
      joint[i].set_continuity(jit->second);
    }
    if (i==1)
    {
      j=jit->first;
      joint[i].set_right_fp(0.5*(pt[j+1]-pt[j])/(length_point[j+1]-length_point[j])/t_scale);
    }
    else if (i==nsegs-1)
    {
      j=jit->first;
      joint[i].set_left_fp(0.5*(pt[j]-pt[j-1])/(length_point[j]-length_point[j-1])/t_scale);
    }
  }

  // fill the fit points and max degree
  jit=joint_point.begin();
  for (j=0; j<nsegs; ++j, ++jit)
  {
    joint_collection_type::iterator jitn(jit); ++jitn;
    max_degree[j]=max_degree_param;
    for (i=jit->first+1; i<jitn->first; ++i)
    {
      fit_point[j].add_point(pt[i]);
    }
  }

  // debug printing of fit points
//  for (j=0; j<nsegs; ++j)
//  {
//    for (i=0; i<fit_point[j].npoints(); ++i)
//    {
//      std::cout << "segment " << std::setw(2) << j << " point(" << std::setw(2) << i << ")=" << fit_point[j].get_point(i) << std::endl;
//    }
//  }

  general_creator_type gc;
  piecewise_curve_type c;
  bool rtn_flag;


  // create curve
  rtn_flag=gc.set_conditions(joint, fit_point, max_degree, false);
  if (!rtn_flag)
  {
    std::cerr << "Could not set the airfoil curve conditions!" << std::endl;
    return EXIT_FAILURE;
  }
  gc.set_t0(t[0]);
  for (i=0; i<nsegs; ++i)
  {
    gc.set_segment_dt(t[i+1]-t[i], i);
  }
  rtn_flag=gc.create(c);
  if (!rtn_flag)
  {
    std::cerr << "Could not create the airfoil curve!" << std::endl;
    return EXIT_FAILURE;
  }

//    std::cout << "subplot(5, 1, 1);" << std::endl;
//    octave_print(std::cout, "af", pts);
//    octave_print(std::cout, "free", bez);
//    std::cout << "plot(af_x, af_y, 'ok', free_x, free_y, '-k');" << std::endl;
//    std::cout << "axis([0,1.1,-0.1,0.1]);" << std::endl;
//    std::cout << "title('" << d << "th Order Bezier - Free Fit');" << std::endl;
//    std::cout << "text(0.3, 0, '\\Sigma d^2=" << err << "');" << std::endl;
  std::cout.flush();
  eli::test::octave_start(1);
  for (i=0; i<npts; ++i)
  {
    char buf[3];
    std::string nm="fp";

    std::sprintf(buf, "%02d", i);
    nm=nm+buf;
    eli::test::octave_print(1, pt[i], nm);
  }
  eli::test::octave_print(1, c, "piecewise");
  eli::test::octave_finish(1);
  std::cout << "axis equal" << std::endl;

  std::cout.flush();
  eli::test::octave_start(2);
  eli::test::octave_param_print(2, c, "piecewise1");
  eli::test::octave_finish(2);

  return EXIT_SUCCESS;
}

