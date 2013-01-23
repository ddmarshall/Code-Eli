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

#ifndef eli_poly_root_closed_form_hpp
#define eli_poly_root_closed_form_hpp

#include <vector>
#include <limits>

#include "eli/constants/math.hpp"

namespace eli
{
  namespace poly
  {
    namespace root
    {
      template<typename data__>
      int closed_form_0(std::vector<data__> &root, const data__ a[])
      {
        root.resize(0);

        if (a[0]==0)
          return std::numeric_limits<int>::max();

        return 0;
      }

      template<typename data__>
      int closed_form_1(std::vector<data__> &root, const data__ a[])
      {
        if (a[1]==0)
          return closed_form_0(root, a);

        root.resize(1);
        root[0]=-a[0]/a[1];
        return 1;
      }

      template<typename data__>
      int closed_form_2(std::vector<data__> &root, const data__ a[])
      {
        if (a[2]==0)
          return closed_form_1(root, a);

        data__ desc=a[1]*a[1]-4*a[2]*a[0];
        if (desc<0)
        {
          root.resize(0);
          return 0;
        }
        root.resize(2);

        data__ sz=std::max(std::max(std::abs(a[0]), std::abs(a[1])), std::abs(a[2]));
        if (desc<=sz*std::numeric_limits<data__>::epsilon())
          desc=static_cast<data__>(0);
        if (desc==0)
        {
          root[0]=-a[1]/(2*a[2]);
          root[1]=root[0];
          return 2;
        }

        data__ term=-(a[1]+((a[1]<0) ? -1 : 1)*std::sqrt(desc))/2;
        root[0]=term/a[2];
        root[1]=a[0]/term;

        std::sort(root.begin(), root.end());
        return 2;
      }

      /** inspired by numerical recipes section 5.6
        */
      template<typename data__>
      int closed_form_3(std::vector<data__> &root, const data__ a[])
      {
        if (a[3]==0)
          return closed_form_2(root, a);

        data__ aa(a[2]/a[3]), bb(a[1]/a[3]), cc(a[0]/a[3]);

        data__ q(aa*aa/9-bb/3), r, q3, r2;
        q=aa*aa/9-bb/3;
        r=aa*q/3-aa*bb/18+cc/2;
        q3=q*q*q;
        r2=r*r;

        // unique real root and double real root or triple root
        if (r2==q3)
        {
          // triple real root
          if (aa*aa==3*bb)
          {
            root.resize(3);
            root[0]=-aa/3;
            root[1]=-aa/3;
            root[2]=-aa/3;

            return 3;
          }

          data__ cbrt_r(std::cbrt(r));

          root.resize(3);

          root[0]=-2*cbrt_r-aa/3;
          root[1]=cbrt_r-aa/3;
          root[2]=root[1];

          std::sort(root.begin(), root.end());
          return 3;
        }

        // three real roots
        if (r2<q3)
        {
          data__ sqrt_q(std::sqrt(q));
          data__ theta;

          theta=std::acos(r/std::sqrt(q3));
          root.resize(3);
          root[0]=-2*sqrt_q*std::cos(theta/3)-aa/3;
          root[1]=-2*sqrt_q*std::cos((theta+eli::constants::math<data__>::two_pi())/3)-aa/3;
          root[2]=-2*sqrt_q*std::cos((theta-eli::constants::math<data__>::two_pi())/3)-aa/3;

          std::sort(root.begin(), root.end());
          return 3;
        }

        // must have only one real root
        data__ newa;
        newa=-((r<0) ? -1 : 1)*std::cbrt(std::abs(r)+std::sqrt(r2-q3));
        root.resize(1);
        root[0]=(newa+q/newa)-aa/3;

        return 1;
      }

      template<typename data__>
      int closed_form_4(std::vector<data__> &root, const data__ aa[])
      {
//         data__ zero(0), one(1), two(2), three(3), four(4), six(6), ninety(90);

        if (aa[4]==0)
          return closed_form_3(root, aa);

        /** based on Shmakov IJPAM 71(2)
         */
        data__ a=aa[3]/aa[4], b=aa[2]/aa[4], c=aa[1]/aa[4], d=aa[0]/aa[4];
        data__ sz=std::max(std::max(std::abs(a), std::abs(b)), std::max(std::abs(c), std::abs(d)));
        data__ delta, DELTA;
//         std::cout << "\ta=" << a << "\tb=" << b << "\tc=" << c << "\td=" << d << std::endl;

        // NOTE: These can actually be set to a variety of values
//         data__ alpha(2), beta(3), gamma(0); // made up
//         data__ alpha(0), beta(0), gamma(0); // Ferrari-Lagrange
//         data__ alpha(-1.0/4.0), beta(1), gamma(-1.0/4.0); // Descartes-Euler-Cardano
        data__ alpha(0), beta(1/3), gamma(0); // new scheme
        delta=alpha*a*a+beta*b;
        DELTA=gamma*a;

        // set the cubic coefficients and solve cubic equation for y
        data__ as, bs, cs, ds, A, B, C, y;
        data__ ay[4];
        std::vector<data__> ysroot;
        as=a+4*DELTA;
        bs=b+(3*a+6*DELTA)*DELTA;
        cs=c+(2*b+(3*a+4*DELTA)*DELTA)*DELTA;
        ds=d+(c+(b+(a+DELTA)*DELTA)*DELTA)*DELTA;
        A=3*delta-bs;
        B=(3*delta-2*bs)*delta+as*cs-4*ds;
        C=((delta-bs)*delta+(as*cs-4*ds))*delta+4*bs*ds-as*as*ds-cs*cs;
        ay[0]=C; ay[1]=B; ay[2]=A; ay[3]=1;
        closed_form_3(ysroot, ay);
//         std::cout << "A=" << A << "\tB=" << B << "\tC=" << C << std::endl;
//         std::cout << "A=" << 0 << "\tB=" << a*c-b*b/3-4*d << "\tC=" << a*b*c/3-a*a*d-2*b*b*b/27-c*c+8*b*d/3 << std::endl;
//         for (size_t i=0; i<ysroot.size(); ++i)
//           std::cout << "\tys[" << i << "]=" << ysroot[i];
//         std::cout << std::endl;

        // get the largest root to try and minimize chance of complex roots below
//         std::cout << "nyroots=" << ysroot.size() << std::endl;
        y=ysroot.back()+delta;
        // NOTE: if ys is approximately delta, then can have loss of significance and
        //       really should chose a different delta.
        // Note: have to compensate for low precision trig function calls in getting root of cubic
//         if (typeid(data__)==typeid(long double))
//         {
//           if (std::abs(y)<=sz*sz*std::numeric_limits<double>::epsilon())
//             y=zero;
//         }
//         else
        {
          if (std::abs(y)<=sz*sz*std::numeric_limits<data__>::epsilon())
            y=static_cast<data__>(0);
        }

        // solve the two quadratic equations. Need to track the discriminants in case complex
        // values for g or h
        data__ gprime, hprime, gdes, hdes, sigma;
        bool gcomplex, hcomplex;

        gprime=a/2;
        gdes=a*a-4*b+4*y;
//         std::cout << "gdes=" << gdes << "\tepsilon=" << 3*sz*sz*std::numeric_limits<data__>::epsilon() << std::endl;
        // Note: have to compensate for low precision trig function calls in getting root of cubic
//         if (typeid(data__)==typeid(long double))
//         {
//           if (std::abs(gdes)<=30*sz*sz*std::numeric_limits<double>::epsilon())
//             gdes=zero;
//         }
//         else
        {
          if (std::abs(gdes)<=3*sz*sz*std::numeric_limits<data__>::epsilon())
            gdes=static_cast<data__>(0);
        }
        if (gdes<0)
        {
//           std::cout << "Complex g quadratic roots!" << std::endl;
          gcomplex=true;
          gdes=std::sqrt(-gdes)/2;
        }
        else
        {
          gcomplex=false;
          gdes=std::sqrt(gdes)/2;
        }

        hprime=y/2;
        hdes=y*y-4*d;
//         std::cout << "hdes=" << hdes << "\tepsilon=" << 2*sz*sz*std::numeric_limits<data__>::epsilon() << std::endl;
        // Note: have to compensate for low precision trig function calls in getting root of cubic
//         if (typeid(data__)==typeid(long double))
//         {
//           if (std::abs(hdes)<=20*sz*sz*std::numeric_limits<double>::epsilon())
//             hdes=zero;
//         }
//         else
        {
          if (std::abs(hdes)<=2*sz*sz*std::numeric_limits<data__>::epsilon())
            hdes=static_cast<data__>(0);
        }
        if (hdes<0)
        {
//           std::cout << "Complex h quadratic roots!" << std::endl;
          hcomplex=true;
          hdes=std::sqrt(-hdes)/2;
        }
        else
        {
          hcomplex=false;
          hdes=std::sqrt(hdes)/2;
        }

        // need to find correct sign for h discriminant (sigma) that solves original equations
//         std::cout << "gprime=" << gprime << "\tgdes=" << gdes << "\thprime=" << hprime << "\thdes=" << hdes << std::endl;
        if ((gdes!=0) && (hdes!=0))
        {
          // catch special case when c is zero
          if (c==0)
          {
            if (gprime*hprime*gdes*hdes>=0)
              sigma=1;
            else
              sigma=-1;
          }
          // catch special case when gprime or hprime are zero
          else if ((gprime==0) || (hprime==0))
          {
            if (gdes*hdes>0)
             sigma=-1;
            else
              sigma=1;
          }
          else
          {
            data__ sig_cal=(gprime*hprime-c/2)/(gdes*hdes);

            // sig_cal should be very close to -1 or +1
            assert((std::abs(sig_cal)>0.9) && (std::abs(sig_cal)<1.1));
            if (sig_cal>=0)
              sigma=1;
            else
              sigma=-1;
          }
//           std::cout << "sigma=" << sigma << std::endl;
//           std::cout << "sigma(1)=" << (gprime*hprime-c/2)/(gdes*hdes) << std::endl;
//           std::cout << "zero?=" << (gprime+gdes)*(hprime-hdes)+(gprime-gdes)*(hprime+hdes)-c << std::endl;
//           std::cout << "zero?=" << (gprime+gdes)*(hprime+hdes)+(gprime-gdes)*(hprime-hdes)-c << std::endl;
        }
        else
        {
          sigma=1;
        }

        // now calculate the 4 roots
        data__ x12des, x34des;
        x12des=gprime*gprime+(gcomplex ? -1 : 1)*gdes*gdes-4*hprime;
        x34des=x12des;
        if (gcomplex || hcomplex)
        {
          if (hcomplex && !gcomplex)
          {
            root.resize(0);
            return 0;
          }

          // should we ever get here?
          assert(false);
        }
        else
        {
          // no complex math needs to be done
          x12des+=2*gprime*gdes-4*sigma*hdes;
          // Note: have to compensate for low precision math in long double
//           if (typeid(data__)==typeid(long double))
//           {
//             if (std::abs(x12des)<=ninety*sz*sz*std::numeric_limits<double>::epsilon())
//               x12des=zero;
//           }
//           else
          {
            if (std::abs(x12des)<=90*sz*sz*std::numeric_limits<data__>::epsilon())
              x12des=static_cast<data__>(0);
          }
          x34des+=-2*gprime*gdes+4*sigma*hdes;
          // Note: have to compensate for low precision math in long double
//           if (typeid(data__)==typeid(long double))
//           {
//             if (std::abs(x34des)<=ninety*sz*sz*std::numeric_limits<double>::epsilon())
//               x34des=zero;
//           }
//           else
          {
            if (std::abs(x34des)<=90*sz*sz*std::numeric_limits<data__>::epsilon())
              x34des=static_cast<data__>(0);
          }
//           std::cout << "x12des=" << x12des << "\tx34des=" << x34des << "\tepsilon=" << 90*sz*sz*std::numeric_limits<data__>::epsilon() << std::endl;
          if (x12des>=0)
          {
            x12des=std::sqrt(x12des);
            root.push_back((-(gprime+gdes)+x12des)/2);
            root.push_back((-(gprime+gdes)-x12des)/2);
          }
          if (x34des>=0)
          {
            x34des=std::sqrt(x34des);
            root.push_back((-(gprime-gdes)+x34des)/2);
            root.push_back((-(gprime-gdes)-x34des)/2);
          }
        }

//         for (size_t i=0; i<root.size(); ++i)
//           std::cout << "\troot[" << i << "]=" << root[i];
//         std::cout << std::endl;
        std::sort(root.begin(), root.end());
        return root.size();
      }

      template<typename data__>
      int closed_form(std::vector<data__> &root, const data__ a[], size_t degree)
      {
        switch (degree)
        {
          case(0):
          {
            return closed_form_0(root, a);
            break;
          }
          case(1):
          {
            return closed_form_1(root, a);
            break;
          }
          case(2):
          {
            return closed_form_2(root, a);
            break;
          }
          case(3):
          {
            return closed_form_3(root, a);
            break;
          }
          case(4):
          {
            return closed_form_4(root, a);
            break;
          }
          default:
          {
            assert(false);
            return -1;
            break;
          }
        }
      }
    }
  }
}

#endif
