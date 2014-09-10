/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_curve_pseudo_cst_base_hpp
#define eli_geom_curve_pseudo_cst_base_hpp

#include "eli/code_eli.hpp"

#include <cmath>

#include "eli/util/tolerance.hpp"

#include "eli/geom/general/continuity.hpp"

#include "eli/geom/curve/pseudo/explicit_bezier.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace pseudo
      {
        template<typename data__, typename tol__=eli::util::tolerance<data__> >
        class cst_base
        {
          private:
            typedef pseudo::explicit_bezier<data__, tol__> shape_curve_type;

          public:
            typedef typename shape_curve_type::data_type data_type;
            typedef typename shape_curve_type::point_type point_type;
            typedef typename shape_curve_type::dimension_type dimension_type;
            typedef typename shape_curve_type::control_point_type control_point_type;
            typedef typename shape_curve_type::index_type index_type;
            typedef typename shape_curve_type::tolerance_type tolerance_type;

          public:
            cst_base() : N1(0.5), N2(1), delta_te(0), shape_function(1) {}
            cst_base(const index_type &dim) : N1(0.5), N2(1), delta_te(0), shape_function(dim) {}
            cst_base(const data_type &N1r, const data_type &N2r, const index_type &dim)
              : N1(N1r), N2(N2r), delta_te(0), shape_function(dim) {}
            cst_base(const cst_base<data_type, tolerance_type> &cst)
              : N1(cst.N1), N2(cst.N2), delta_te(cst.delta_te), shape_function(cst.shape_function) {}
            virtual ~cst_base() {}

            bool operator==(const cst_base<data_type, tolerance_type> &cst) const
            {
              if (this==&cst)
                return true;
              if ((N1!=cst.N1) || (N2!=cst.N2))
                return false;
              return (shape_function==cst.shape_function);
            }

            bool operator!=(const cst_base<data_type, tolerance_type> &cst) const
            {
              return !operator==(cst);
            }

            cst_base & operator=(const cst_base<data_type, tolerance_type> &cst)
            {
              if (this == &cst)
              {
                return (*this);
              }

              N1=cst.N1;
              N2=cst.N2;
              shape_function=cst.shape_function;

              return (*this);
            }

            static dimension_type dimension() {return 2;}

            void resize(const index_type &dim)
            {
              shape_function.resize(dim);
            }

            index_type degree() const
            {
              return shape_function.degree();
            }

            const data_type & get_trailing_edge_thickness() const {return delta_te;}
            void set_trailing_edge_thickness(const data_type &dte)
            {
              if (dte<=0)
                delta_te=0;
              else
                delta_te=dte;
            }

            void get_shape_parameters(data_type &N1_out, data_type &N2_out) const
            {
              N1_out=N1;
              N2_out=N2;
            }

            data_type get_t0() const {return static_cast<data_type>(0);}
            data_type get_tmax() const {return static_cast<data_type>(1);}

            bool open() const
            {
              return true;
            }
            bool closed() const
            {
              return false;
            }

            void degree_promote()
            {
              shape_function.degree_promote();
            }

            bool degree_demote(const geom::general::continuity &continuity_degree=geom::general::C0)
            {
              return shape_function.degree_demote(continuity_degree);
            }

            void set_control_point(const control_point_type &cp_in, const index_type &i)
            {
              shape_function.set_control_point(cp_in, i);
            }

            control_point_type get_control_point(const index_type &i) const
            {
              return shape_function.get_control_point(i);
            }

            point_type f(const data_type &t) const
            {
              point_type rtn;
              rtn = class_f(t)*shape_function.f(t);
              rtn(1)+=t*delta_te;
              return rtn;
            }

            point_type fp(const data_type &t) const
            {
              point_type rtn;
              rtn = class_fp(t)*shape_function.f(t)+class_f(t)*shape_function.fp(t);
              rtn(1)+=delta_te;
              return rtn;
            }

            point_type fpp(const data_type &t) const
            {
              return class_fpp(t)*shape_function.f(t)+2*class_fp(t)*shape_function.fp(t)+class_f(t)*shape_function.fpp(t);
            }

          protected:
            void set_shape_params(const data_type &N1_in, const data_type &N2_in)
            {
              // do some range checking
              if ((N1_in<0) || (N1_in>1))
              {
                assert(false);
                return;
              }
              if ((N2_in<0) || (N2_in>1))
              {
                assert(false);
                return;
              }

              N1=N1_in;
              N2=N2_in;
            }

          private:
            data_type N1, N2;
            data_type delta_te;
            shape_curve_type shape_function;

          private:
            data_type class_f(const data_type &t) const
            {
              const data_type one_half(static_cast<data_type>(0.5));

              // front and back are special
              if ((t==0) || (t==1))
              {
                return 0;
              }

              // special class functions
              if (N1==0)
              {
                if (N2==0)
                {
                  return 1;
                }
                if (N2==1)
                {
                  return 1-t;
                }
                if (N2==one_half)
                {
                  return std::sqrt(1-t);
                }

                return std::pow(1-t, N2);
              }
              if (N1==1)
              {
                if (N2==0)
                {
                  return t;
                }
                if (N2==1)
                {
                  return t*(1-t);
                }
                if (N2==one_half)
                {
                  return t*std::sqrt(1-t);
                }

                return t*std::pow(1-t, N2);
              }
              if (N1==one_half)
              {
                if (N2==0)
                {
                  return std::sqrt(t);
                }
                if (N2==1)
                {
                  return std::sqrt(t)*(1-t);
                }
                if (N2==one_half)
                {
                  return std::sqrt(t*(t-1));
                }

                return std::sqrt(t)*std::pow(1-t, N2);
              }
              if (N2==0)
              {
                return std::pow(t, N1);
              }
              if (N2==1)
              {
                return std::pow(t, N1)*(1-t);
              }
              if (N2==one_half)
              {
                return std::pow(t, N1)*std::sqrt(1-t);
              }

              return std::pow(t, N1)*std::pow(1-t, N2);
            }

            data_type class_fp(const data_type &t) const
            {
              const data_type one_half(static_cast<data_type>(0.5)), three_half(static_cast<data_type>(1.5));
              const data_type one_quarter(static_cast<data_type>(0.25));
              const data_type max_val(std::numeric_limits<data_type>::max());

              // front is special
              if (t==0)
              {
                if (N1==0)
                {
                  if (N2==0)
                  {
                    return 0;
                  }
                  if (N2==1)
                  {
                    return -N2;
                  }
                  if (N2==one_half)
                  {
                    return N2;
                  }

                  return N2;
                }
                if (N1==1)
                {
                  return N1;
                }
                if (N1==one_half)
                {
                  return max_val;
                }
                if (N2==0)
                {
                  return max_val;
                }
                if (N2==1)
                {
                  return max_val;
                }
                if (N2==one_half)
                {
                  return max_val;
                }

                return max_val;
              }

              // back is special
              if (t==1)
              {
                if (N1==0)
                {
                  if (N2==0)
                  {
                    return N2;
                  }
                  if (N2==1)
                  {
                    return -N2;
                  }
                  if (N2==one_half)
                  {
                    return max_val;
                  }

                  return max_val;
                }
                if (N1==1)
                {
                  if (N2==0)
                  {
                    return N1;
                  }
                  if (N2==1)
                  {
                    return -N1;
                  }
                  if (N2==one_half)
                  {
                    return -max_val;
                  }

                  return -max_val;
                }
                if (N1==one_half)
                {
                  if (N2==0)
                  {
                    return N1;
                  }
                  if (N2==1)
                  {
                    return -N2;
                  }
                  if (N2==one_half)
                  {
                    return -max_val;
                  }

                  return -max_val;
                }
                if (N2==0)
                {
                  return N1;
                }
                if (N2==1)
                {
                  return -N2;
                }
                if (N2==one_half)
                {
                  return -max_val;
                }

                return -max_val;
              }

              // special class functions
              if (N1==0)
              {
                if (N2==0)
                {
                  return 1;
                }
                if (N2==1)
                {
                  return -N2;
                }
                if (N2==one_half)
                {
                  return one_half/std::sqrt(1-t);
                }

                return N2*std::pow(1-t, N2-1);
              }
              if (N1==1)
              {
                if (N2==0)
                {
                  return N1;
                }
                if (N2==1)
                {
                  return 1-2*t;
                }
                if (N2==one_half)
                {
                  return (1-three_half*t)/std::sqrt(1-t);
                }
              }
              if (N1==one_half)
              {
                if (N2==0)
                {
                  return -one_quarter/(t*std::sqrt(t));
                }
                if (N2==1)
                {
                  return one_half*(1-3*t)/std::sqrt(t);
                }
                if (N2==one_half)
                {
                  return one_half*(1-2*t)/std::sqrt(t*(1-t));
                }

                return one_half*std::sqrt(t)*std::pow(1-t, N2-1)*(1-(1+2*N2)*t);
              }
              if (N2==0)
              {
                return N1*std::pow(t, N1-1);
              }
              if (N2==1)
              {
                return std::pow(t, N1-1)*(N1-(N1+1)*t);
              }
              if (N2==one_half)
              {
                return one_half*(2*N1-(1+2*N1)*t)*std::pow(t, N1-1)/std::sqrt(1-t);
              }

              return (N1-(N1+N2)*t)*std::pow(t, N1-1)*std::pow(1-t, N2-1);
            }

            data_type class_fpp(const data_type &t) const
            {
              const data_type one_half(static_cast<data_type>(0.5));
              const data_type one_quarter(static_cast<data_type>(0.25)), three_quarter(static_cast<data_type>(0.75));
              const data_type max_val(std::numeric_limits<data_type>::max());

              // front is special
              if (t==0)
              {
                if (N1==0)
                {
                  if (N2==0)
                  {
                    return 0;
                  }
                  if (N2==1)
                  {
                    return 0;
                  }
                  if (N2==one_half)
                  {
                    return -one_quarter;
                  }

                  return N2*(N2-1);
                }
                if (N1==1)
                {
                  if (N2==0)
                  {
                    return 0;
                  }
                  if (N2==1)
                  {
                    return -2;
                  }
                  if (N2==one_half)
                  {
                    return -1;
                  }

                  return -2*N2;
                }
                if (N1==one_half)
                {
                  return -max_val;
                }
                if (N2==0)
                {
                  return -max_val;
                }
                if (N2==1)
                {
                  return -max_val;
                }
                if (N2==one_half)
                {
                  return -max_val;
                }

                return -max_val;
              }

              // back is special
              if (t==1)
              {
                if (N1==0)
                {
                  if (N2==0)
                  {
                    return 0;
                  }
                  if (N2==1)
                  {
                    return 0;
                  }
                  if (N2==one_half)
                  {
                    return -max_val;
                  }

                  return -max_val;
                }
                if (N1==1)
                {
                  if (N2==0)
                  {
                    return 0;
                  }
                  if (N2==1)
                  {
                    return -2;
                  }
                  if (N2==one_half)
                  {
                    return -max_val;
                  }

                  return -max_val;
                }
                if (N1==one_half)
                {
                  if (N2==0)
                  {
                    return -one_quarter;
                  }
                  if (N2==1)
                  {
                    return -1;
                  }
                  if (N2==one_half)
                  {
                    return -max_val;
                  }

                  return -max_val;
                }
                if (N2==0)
                {
                  return N2*(N2-1);
                }
                if (N2==1)
                {
                  return -2*N1;
                }
                if (N2==one_half)
                {
                  return -max_val;
                }

                return -max_val;
              }

              // special class functions
              if (N1==0)
              {
                if (N2==0)
                {
                  return 0;
                }
                if (N2==1)
                {
                  return 0;
                }
                if (N2==one_half)
                {
                  return -one_quarter/std::sqrt(1-t)/(1-t);
                }

                return N2*(N2-1)*std::pow(1-t, N2-2);
              }
              if (N1==1)
              {
                if (N2==0)
                {
                  return 0;
                }
                if (N2==1)
                {
                  return -2;
                }
                if (N2==one_half)
                {
                  return (three_quarter*t-1)/((1-t)*std::sqrt(1-t));
                }

                return std::pow(1-t, N2-2)*(N2*((N2-1)*t-2));
              }
              if (N1==one_half)
              {
                if (N2==0)
                {
                  return -one_quarter/(t*std::sqrt(t));
                }
                if (N2==1)
                {
                  return -one_quarter*(1+3*t)/(t*std::sqrt(t));
                }
                if (N2==one_half)
                {
                  data_type tmp(t*(1-t));
                  return -one_quarter/(tmp*std::sqrt(tmp));
                }

                return one_quarter*std::pow(1-t, N2-2)*((2*N2-1)*t*((2*N2-1)*t-2)-1)/(t*std::sqrt(t));
              }
              if (N2==0)
              {
                return N1*(N1-1)*std::pow(t, N1-2);
              }
              if (N2==1)
              {
                return N1*(N1-1-(N1+1)*t)*std::pow(t, N1-2);
              }
              if (N2==one_half)
              {
                data_type N12(N1*N1);
                return one_quarter*(((4*N12-1)*t+4*(N1-2*N12))*t+4*(N12-N1))*std::pow(t, N1-2)/((1-t)*std::sqrt(1-t));
              }

              return (t*(t*(N1*(1-N1)+N2*(1-N2))+2*N1*(1-(N1+N2)))+N1*(1-N1))*std::pow(t, N1-2)*std::pow(1-t, N2-2);
            }
        };
      }
    }
  }
}

#endif
