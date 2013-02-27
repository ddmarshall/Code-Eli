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

#ifndef eli_geom_surface_curvature_hpp
#define eli_geom_surface_curvature_hpp

#include <iomanip>
#include <limits> // numeric_limits

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      namespace internal
      {
        template<typename surface__>
        void calculate_surface_terms(typename surface__::data_type &e, typename surface__::data_type &f, typename surface__::data_type &g,
                                     typename surface__::data_type &ee, typename surface__::data_type &ff, typename surface__::data_type &gg,
                                     const surface__ &s, const typename surface__::data_type &u, const typename surface__::data_type &v)
        {
          typename surface__::point_type Su, Sv, Suu, Suv, Svv, n;

          // surfaces in two dimensions have no curvature
          if (Su.innerSize()==2)
          {
            e=0;
            f=0;
            g=0;
            ee=1;
            ff=1;
            gg=1;
            return;
          }
          else if (Su.innerSize()!=3)
          {
            // how generalize curvature to surfaces in hyper-space?
            assert(false);
            return;
          }

          // calculate the surface derivatives
          Su=s.f_u(u, v);
          Sv=s.f_v(u, v);
          Suu=s.f_uu(u, v);
          Suv=s.f_uv(u, v);
          Svv=s.f_vv(u, v);
          n=s.normal(u, v);

          // calculate terms
          e=n.dot(Suu);
          f=n.dot(Suv);
          g=n.dot(Svv);
          ee=Su.dot(Su);
          ff=Su.dot(Sv);
          gg=Sv.dot(Sv);
          assert((ee*gg-ff*ff)>=0);
          assert(std::abs(1-std::sqrt(ee*gg-ff*ff)/Su.cross(Sv).norm())<1e-4);
        }

        template<typename surface__>
        void principal_curvature_calc(typename surface__::data_type &kmax, typename surface__::data_type &kmin,
                                      typename surface__::point_type &kmax_dir, typename surface__::point_type &kmin_dir, typename surface__::point_type &n,
                                      const surface__ &s, const typename surface__::data_type &u, const typename surface__::data_type &v, bool calc_dir)
        {
          typename surface__::data_type H, K, tmp;

          // check to make sure have valid curve
          assert(s.degree_u()>0);
          assert(s.degree_v()>0);

          // check to make sure given valid parametric value
          assert((u>=0) && (u<=1));
          assert((v>=0) && (v<=1));

          // calculate the surface parameters
          typename surface__::data_type e, f, g, ee, ff, gg;

          internal::calculate_surface_terms(e, f, g, ee, ff, gg, s, u, v);

          // calculate the mean and gaussian curvatures
          tmp=ee*gg-ff*ff;
          if (tmp==0)
          {
            kmax=std::numeric_limits<typename surface__::data_type>::max();
            kmin=std::numeric_limits<typename surface__::data_type>::max();
            return;
          }
          else
          {
            H=(e*gg-2*f*ff+g*ee)/(2*tmp);
            K=(e*g-f*f)/tmp;
          }


          // calculate the principal curvatures
          tmp=std::sqrt(H*H-K);
          kmax=H+tmp;
          kmin=H-tmp;

          // check if round-off error has crept in
          typename surface__::data_type tol=std::numeric_limits<typename surface__::data_type>::epsilon();
#ifdef ELI_QD_FOUND
          if (typeid(typename surface__::data_type)==typeid(dd_real))
          {
            tol=1e-9*static_cast<typename surface__::data_type>(std::numeric_limits<double>::epsilon());
          }
          else if (typeid(typename surface__::data_type)==typeid(qd_real))
          {
            tol=1e-6*static_cast<typename surface__::data_type>(std::numeric_limits<dd_real>::epsilon());
          }
#endif
          if ((tmp<1e-5) && (tmp/std::max(e*e*gg*gg, std::max(f*f*ff*ff, g*g*ee*ee))<tol))
          {
            kmax=H;
            kmin=H;
            tmp=0;
          }

          // calculate the principal curvature directions
          if (calc_dir)
          {
            typename surface__::point_type Su, Sv;
            typename surface__::data_type lambda_max, lambda_min;

            Su=s.f_u(u, v);
            Sv=s.f_v(u, v);
            n=Su.cross(Sv);
            n.normalize();

            // handle umbilic point case specially
            if (tmp==0)
            {
              kmax_dir=Su;
              kmin_dir=n.cross(kmax_dir);
            }
            else
            {
              tmp=g-kmax*gg;
              if (tmp==0)
              {
                kmax_dir=Sv;
              }
              else
              {
                lambda_max=-(f-kmax*ff)/tmp;
                kmax_dir=Su+lambda_max*Sv;
              }

              tmp=g-kmin*gg;
              if (tmp==0)
              {
                kmin_dir=Sv;
              }
              else
              {
                lambda_min=-(f-kmin*ff)/tmp;
                kmin_dir=Su+lambda_min*Sv;
              }
            }
            kmax_dir.normalize();
            kmin_dir.normalize();

            // need to return right-handed system kmax_dir:kmin_dir:n so check direction of kmin_dir
            tmp=n.cross(kmax_dir).dot(kmin_dir);
            if (tmp<0)
            {
              kmin_dir*=-1;
            }
            assert(std::abs(1-std::abs(tmp))<0.001); // tmp1 should be very close to 1 or -1
          }
        }
      }

      template<typename surface__>
      void mean_curvature(typename surface__::data_type &k, const surface__ &s, const typename surface__::data_type &u, const typename surface__::data_type &v)
      {
        // check to make sure have valid curve
        assert(s.degree_u()>0);
        assert(s.degree_v()>0);

        // check to make sure given valid parametric value
        assert((u>=0) && (u<=1));
        assert((v>=0) && (v<=1));

        // calculate the surface parameters
        typename surface__::data_type e, f, g, ee, ff, gg;

        internal::calculate_surface_terms(e, f, g, ee, ff, gg, s, u, v);

        // calculate the curvature
        typename surface__::data_type tmp;
        tmp=ee*gg-ff*ff;
        if ((tmp)==0)
        {
          k=std::numeric_limits<typename surface__::data_type>::max();
          return;
        }

        k=(e*gg-2*f*ff+g*ee)/(2*tmp);
      }

      template<typename surface__>
      void gaussian_curvature(typename surface__::data_type &k, const surface__ &s, const typename surface__::data_type &u, const typename surface__::data_type &v)
      {
        // check to make sure have valid curve
        assert(s.degree_u()>0);
        assert(s.degree_v()>0);

        // check to make sure given valid parametric value
        assert((u>=0) && (u<=1));
        assert((v>=0) && (v<=1));

        // calculate the surface parameters
        typename surface__::data_type e, f, g, ee, ff, gg;

        internal::calculate_surface_terms(e, f, g, ee, ff, gg, s, u, v);

        // calculate the curvature
        typename surface__::data_type tmp;
        tmp=ee*gg-ff*ff;
        if (tmp==0)
        {
          k=std::numeric_limits<typename surface__::data_type>::max();
          return;
        }

        k=(e*g-f*f)/tmp;
    }

      template<typename surface__>
      void principal_curvature(typename surface__::data_type &kmax, typename surface__::data_type &kmin, const surface__ &s, const typename surface__::data_type &u, const typename surface__::data_type &v)
      {
        typename surface__::point_type kmax_dir, kmin_dir, n;
        internal::principal_curvature_calc(kmax, kmin, kmax_dir, kmin_dir, n, s, u, v, false);
      }

      template<typename surface__>
      void principal_curvature(typename surface__::data_type &kmax, typename surface__::data_type &kmin,
                               typename surface__::point_type &kmax_dir, typename surface__::point_type &kmin_dir, typename surface__::point_type &n,
                               const surface__ &s, const typename surface__::data_type &u, const typename surface__::data_type &v)
      {
        internal::principal_curvature_calc(kmax, kmin, kmax_dir, kmin_dir, n, s, u, v, true);
      }
    }
  }
}

#endif
