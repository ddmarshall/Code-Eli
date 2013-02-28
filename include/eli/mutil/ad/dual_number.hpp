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

#ifndef eli_mutil_ad_dual_number_hpp
#define eli_mutil_ad_dual_number_hpp

//
// dual value type
//
namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      // forward declare dual
      template<typename data__, bool comp_real_only> class dual;

      namespace dual_number
      {
        //
        // constant type
        //
        template <typename data__>
        class constant
        {
          public:
            typedef data__ data_type;

          private:
            const data_type & val;

          public:
            constant(const data_type &v) : val(v) {}

            data_type real() const {return val;}
            data_type nonreal() const {return 0;}
        };

        //
        // expression type
        //
        template <typename T__, bool comp_real_only>
        class expression
        {
          private:
            const T__ expr;

          public:
            typedef typename T__::data_type data_type;

            explicit expression(const T__ &e) : expr(e) {}

            const data_type real() const {return expr.real();}
            const data_type nonreal() const {return expr.nonreal();}

            // operator==
            bool operator==(const data_type &v) const
            {
              eli::mutil::ad::dual<data_type, comp_real_only> d(this->real(), this->nonreal());
              return d==v;
            }
            bool operator==(const expression<T__, comp_real_only> &e) const
            {
              eli::mutil::ad::dual<data_type, comp_real_only> dl(this->real(), this->nonreal()), dr(e.real(), e.nonreal());
              return dl==dr;
            }

            // operator==
            bool operator!=(const data_type &v) const
            {
              return !operator==(v);
            }
            bool operator!=(const expression<T__, comp_real_only> &e) const
            {
              return !operator==(e);
            }

            // operator<=
            bool operator<=(const data_type &v) const
            {
              eli::mutil::ad::dual<data_type, comp_real_only> d(this->real(), this->nonreal());
              return d<=v;
            }
            bool operator<=(const expression<T__, comp_real_only> &e) const
            {
              eli::mutil::ad::dual<data_type, comp_real_only> dl(this->real(), this->nonreal()), dr(e.real(), e.nonreal());
              return dl<=dr;
            }

            // operator<
            bool operator<(const data_type &v) const
            {
              eli::mutil::ad::dual<data_type, comp_real_only> d(this->real(), this->nonreal());
              return d<v;
            }
            bool operator<(const expression<T__, comp_real_only> &e) const
            {
              eli::mutil::ad::dual<data_type, comp_real_only> dl(this->real(), this->nonreal()), dr(e.real(), e.nonreal());
              return dl<dr;
            }

            // operator>=
            bool operator>=(const data_type &v) const
            {
              eli::mutil::ad::dual<data_type, comp_real_only> d(this->real(), this->nonreal());
              return d>=v;
            }
            bool operator>=(const expression<T__, comp_real_only> &e) const
            {
              eli::mutil::ad::dual<data_type, comp_real_only> dl(this->real(), this->nonreal()), dr(e.real(), e.nonreal());
              return dl>=dr;
            }

            // operator>
            bool operator>(const data_type &v) const
            {
              eli::mutil::ad::dual<data_type, comp_real_only> d(this->real(), this->nonreal());
              return d>v;
            }
            bool operator>(const expression<T__, comp_real_only> &e) const
            {
              eli::mutil::ad::dual<data_type, comp_real_only> dl(this->real(), this->nonreal()), dr(e.real(), e.nonreal());
              return dl>dr;
            }
        };
      }
    }
  }
}

// comparison operators
#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4244)
#endif
template <typename data2__, typename data__, bool comp_real_only>
bool operator==(const data2__ &v, const eli::mutil::ad::dual_number::expression<data__, comp_real_only> &e)
{
  return e==v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator!=(const data2__ &v, const eli::mutil::ad::dual_number::expression<data__, comp_real_only> &e)
{
  return e!=v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator<=(const data2__ &v, const eli::mutil::ad::dual_number::expression<data__, comp_real_only> &e)
{
  return e>=v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator<(const data2__ &v, const eli::mutil::ad::dual_number::expression<data__, comp_real_only> &e)
{
  return e>v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator>=(const data2__ &v, const eli::mutil::ad::dual_number::expression<data__, comp_real_only> &e)
{
  return e<=v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator>(const data2__ &v, const eli::mutil::ad::dual_number::expression<data__, comp_real_only> &e)
{
  return e<v;
}
#ifdef _MSC_VER
# pragma warning(pop)
#endif

namespace eli
{
  namespace mutil
  {
    namespace ad
    {
      template <typename data__, bool comp_real_only>
      class dual
      {
        public:
          typedef data__ data_type;

        private:
          data__ val[2];

        private:
          bool exact_equiv(const data__ &r, const data__ &nr) const
          {
            if (val[0]!=r)
              return false;
            if (val[1]!=nr)
              return false;
            return true;
          }
          bool nearly_equiv(const data__ &r, const data__ &nr, const data__ &eps) const
          {
            data__ rdif(val[0]-r), nrdif(val[1]-nr);
            if (std::abs(rdif)>eps)
              return false;
            if (std::abs(nrdif)>eps)
              return false;
            return true;
          }

          bool exact_less_than(const data__ &/*r*/, const data__ &/*nr*/) const
          {
            // not implemented
            assert(false);
            return false;
          }

          bool exact_less_than_equal(const data__ &/*r*/, const data__ &/*nr*/) const
          {
            // not implemented
            assert(false);
            return false;
          }

          bool exact_greater_than(const data__ &/*r*/, const data__ &/*nr*/) const
          {
            // not implemented
            assert(false);
            return false;
          }

          bool exact_greater_than_equal(const data__ &/*r*/, const data__ &/*nr*/) const
          {
            // not implemented
            assert(false);
            return false;
          }

          bool inexact_equiv(const data__ &r) const
          {
            return val[0]==r;
          }

          bool inexact_less_than(const data__ &r) const
          {
            return val[0]<r;
          }

          bool inexact_less_than_equal(const data__ &r) const
          {
            return val[0]<=r;
          }

          bool inexact_greater_than(const data__ &r) const
          {
            return val[0]>r;
          }

          bool inexact_greater_than_equal(const data__ &r) const
          {
            return val[0]>=r;
          }

        public:
          dual()
          {
  #ifdef DEBUG
            val[0]=static_cast<data_type>(0);
            val[1]=static_cast<data_type>(0);
  #endif
          }

          explicit dual(const data__ &v)
          {
            val[0]=v;
            val[1]=static_cast<data_type>(0);
          }

          dual(const data_type &v1, const data_type &v2)
          {
            val[0]=v1;
            val[1]=v2;
          }

          dual(const dual<data_type, comp_real_only> &d)
          {
            val[0]=d.real();
            val[1]=d.nonreal();
          }

          template<typename data2__>
          explicit dual(const dual_number::constant<data2__> &c)
          {
            val[0]=c.real();
            val[1]=c.nonreal();
          }

          template<typename data2__>
          explicit dual(const dual_number::expression<data2__, comp_real_only> &c)
          {
            val[0]=c.real();
            val[1]=c.nonreal();
          }

          ~dual() {}

          data_type real() const {return val[0];}
          data_type nonreal() const {return val[1];}
          void set_real(const data_type &r) {val[0]=r;}
          void set_nonreal(const data_type &nr) {val[1]=nr;}

          dual<data_type, comp_real_only> & operator=(const data_type &v)
          {
            val[0]=v;
            val[1]=static_cast<data_type>(0);
            return *this;
          }

          dual<data_type, comp_real_only> & operator=(const dual<data_type, comp_real_only> &d)
          {
            if (this!=&d)
            {
              val[0]=d.real();
              val[1]=d.nonreal();
            }

            return *this;
          }

          template <typename data2__>
          dual<data_type, comp_real_only> & operator=(const dual<data2__, comp_real_only> &d)
          {
            if (this!=&d)
            {
              val[0]=d.real();
              val[1]=d.nonreal();
            }

            return *this;
          }

          template <typename data2__>
          dual<data_type, comp_real_only> & operator=(const dual_number::constant<data2__> &c)
          {
            val[0]=c.real();
            return *this;
          }

          template <typename data2__>
          dual<data_type, comp_real_only> & operator=(const dual_number::expression<data2__, comp_real_only> &c)
          {
            val[0]=c.real();
            val[1]=c.nonreal();
            return *this;
          }

          // operator +=
          dual<data_type, comp_real_only> & operator+=(const data_type &d)
          {
            val[0]+=d;
            return *this;
          }
          dual<data_type, comp_real_only> & operator+=(const dual<data_type, comp_real_only> &d)
          {
            val[0]+=d.real();
            val[1]+=d.nonreal();
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator+=(const dual<data2__, comp_real_only> &d)
          {
            val[0]+=d.real();
            val[1]+=d.nonreal();
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator+=(const dual_number::constant<data2__> &c)
          {
            val[0]+=c.real();
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator+=(const dual_number::expression<data2__, comp_real_only> &e)
          {
            val[0]+=e.real();
            val[1]+=e.nonreal();
            return *this;
          }

          // operator -=
          dual<data_type, comp_real_only> & operator-=(const data_type &d)
          {
            val[0]-=d;
            return *this;
          }
          dual<data_type, comp_real_only> & operator-=(const dual<data_type, comp_real_only> &d)
          {
            val[0]-=d.real();
            val[1]-=d.nonreal();
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator-=(const dual<data2__, comp_real_only> &d)
          {
            val[0]-=d.real();
            val[1]-=d.nonreal();
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator-=(const dual_number::constant<data2__> &c)
          {
            val[0]-=c.real();
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator-=(const dual_number::expression<data2__, comp_real_only> &e)
          {
            val[0]-=e.real();
            val[1]-=e.nonreal();
            return *this;
          }

          // operator *=
          dual<data_type, comp_real_only> & operator*=(const data_type &d)
          {
            val[0]*=d;
            val[1]*=d;
            return *this;
          }
          dual<data_type, comp_real_only> & operator*=(const dual<data_type, comp_real_only> &d)
          {
            data_type val0(val[0]);

            val[0]=val0*d.real();
            val[1]=val0*d.nonreal()+val[1]*d.real();
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator*=(const dual<data2__, comp_real_only> &d)
          {
            data_type val0(val[0]);

            val[0]=val0*d.real();
            val[1]=val0*d.nonreal()+val[1]*d.real();
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator*=(const dual_number::constant<data2__> &c)
          {
            val[0]*=c.real();
            val[1]*=c.real();
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator*=(const dual_number::expression<data2__, comp_real_only> &e)
          {
            data_type val0(val[0]);

            val[0]=val0*e.real();
            val[1]=val0*e.nonreal()+val[1]*e.real();
            return *this;
          }

          // operator /=
          dual<data_type, comp_real_only> & operator/=(const data_type &d)
          {
            val[0]/=d;
            val[1]/=d;
            return *this;
          }
          dual<data_type, comp_real_only> & operator/=(const dual<data_type, comp_real_only> &d)
          {
            data_type val0(val[0]);
            data_type c(d.real());

            val[0]=val0/d.real();
            val[1]=(val[1]*c-val0*d.nonreal())/c/c;
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator/=(const dual<data2__, comp_real_only> &d)
          {
            data_type val0(val[0]);
            data2__ c(d.real());

            val[0]=val0/d.real();
            val[1]=(val[1]*c-val0*d.nonreal())/c/c;
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator/=(const dual_number::constant<data2__> &c)
          {
            val[0]/=c.real();
            val[1]/=c.real();
            return *this;
          }
          template <typename data2__>
          dual<data_type, comp_real_only> & operator/=(const dual_number::expression<data2__, comp_real_only> &e)
          {
            data_type val0(val[0]);
            typename dual_number::expression<data2__, comp_real_only>::data_type c(e.real());

            val[0]=val0/e.real();
            val[1]=(val[1]*c-val0*e.nonreal())/c/c;
            return *this;
          }

          bool exact(const data_type &d) const
          {
            return exact_equiv(d, 0);
          }
          bool exact(const dual<data_type, comp_real_only> &d) const
          {
            return exact_equiv(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool exact(const dual<data2__, comp_real_only> &d) const
          {
            return exact_equiv(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool exact(const dual_number::constant<data2__> &c) const
          {
            return exact_equiv(c.real(), 0);
          }
          template <typename data2__>
          bool exact(const dual_number::expression<data2__, comp_real_only> &e) const
          {
            return exact_equiv(e.real(), e.nonreal());
          }

          bool nearly(const data_type &d, const data_type &eps) const
          {
            return nearly_equiv(d, 0, eps);
          }
          bool nearly(const dual<data_type, comp_real_only> &d, const data_type &eps) const
          {
            return nearly_equiv(d.real(), d.nonreal(), eps);
          }
          template <typename data2__>
          bool nearly(const dual<data2__, comp_real_only> &d, const data_type &eps) const
          {
            return nearly_equiv(d.real(), d.nonreal(), eps);
          }
          template <typename data2__>
          bool nearly(const dual_number::constant<data2__> &c, const data_type &eps) const
          {
            return nearly_equiv(c.real(), 0, eps);
          }
          template <typename data2__>
          bool nearly(const dual_number::expression<data2__, comp_real_only> &e, const data_type &eps) const
          {
            return nearly_equiv(e.real(), e.nonreal(), eps);
          }

          // operator ==
          bool operator==(const data_type &d) const
          {
            if (comp_real_only)
              return inexact_equiv(d);
            else
              return exact_equiv(d, 0);
          }
          bool operator==(const dual<data_type, comp_real_only> &d) const
          {
            if (comp_real_only)
              return inexact_equiv(d.real());
            else
              return exact_equiv(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool operator==(const dual<data2__, comp_real_only> &d) const
          {
            if (comp_real_only)
              return inexact_equiv(d.real());
            else
              return exact_equiv(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool operator==(const dual_number::constant<data2__> &c) const
          {
            return operator==(c.real());
          }
          template <typename data2__>
          bool operator==(const dual_number::expression<data2__, comp_real_only> &e) const
          {
            if (comp_real_only)
              return inexact_equiv(e.real());
            else
              return exact_equiv(e.real(), e.nonreal());
          }

          // operator !=
          bool operator!=(const data_type &d) const
          {
            return !operator==(d);
          }
          bool operator!=(const dual<data_type, comp_real_only> &d) const
          {
            return !operator==(d);
          }
          template <typename data2__>
          bool operator!=(const dual<data2__, comp_real_only> &d) const
          {
            return !operator==(d);
          }
          template <typename data2__>
          bool operator!=(const dual_number::constant<data2__> &c) const
          {
            return !operator==(c);
          }
          template <typename data2__>
          bool operator!=(const dual_number::expression<data2__, comp_real_only> &e) const
          {
            return !operator==(e);
          }

          // operator <=
          bool operator<=(const data_type &d) const
          {
            if (comp_real_only)
              return inexact_less_than_equal(d);
            else
              return exact_less_than_equal(d, 0);
          }
          bool operator<=(const dual<data_type, comp_real_only> &d) const
          {
            if (comp_real_only)
              return inexact_less_than_equal(d.real());
            else
              return exact_less_than_equal(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool operator<=(const dual<data2__, comp_real_only> &d) const
          {
            if (comp_real_only)
              return inexact_less_than_equal(d.real());
            else
              return exact_less_than_equal(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool operator<=(const dual_number::constant<data2__> &c) const
          {
            return operator<=(c.real());
          }
          template <typename data2__>
          bool operator<=(const dual_number::expression<data2__, comp_real_only> &e) const
          {
            if (comp_real_only)
              return inexact_less_than_equal(e.real());
            else
              return exact_less_than_equal(e.real(), e.nonreal());
          }

          // operator <
          bool operator<(const data_type &d) const
          {
            if (comp_real_only)
              return inexact_less_than(d);
            else
              return exact_less_than(d, 0);
          }
          bool operator<(const dual<data_type, comp_real_only> &d) const
          {
            if (comp_real_only)
              return inexact_less_than(d.real());
            else
              return exact_less_than(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool operator<(const dual<data2__, comp_real_only> &d) const
          {
            if (comp_real_only)
              return inexact_less_than(d.real());
            else
              return exact_less_than(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool operator<(const dual_number::constant<data2__> &c) const
          {
            return operator<(c.real());
          }
          template <typename data2__>
          bool operator<(const dual_number::expression<data2__, comp_real_only> &e) const
          {
            if (comp_real_only)
              return inexact_less_than(e.real());
            else
              return exact_less_than(e.real(), e.nonreal());
          }

          // operator >=
          bool operator>=(const data_type &d) const
          {
            if (comp_real_only)
              return inexact_greater_than_equal(d);
            else
              return exact_greater_than_equal(d, 0);
          }
          bool operator>=(const dual<data_type, comp_real_only> &d) const
          {
            if (comp_real_only)
              return inexact_greater_than_equal(d.real());
            else
              return exact_greater_than_equal(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool operator>=(const dual<data2__, comp_real_only> &d) const
          {
            if (comp_real_only)
              return inexact_greater_than_equal(d.real());
            else
              return exact_greater_than_equal(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool operator>=(const dual_number::constant<data2__> &c) const
          {
            return operator>=(c.real());
          }
          template <typename data2__>
          bool operator>=(const dual_number::expression<data2__, comp_real_only> &e) const
          {
            if (comp_real_only)
              return inexact_greater_than_equal(e.real());
            else
              return exact_greater_than_equal(e.real(), e.nonreal());
          }

          // operator >
          bool operator>(const data_type &d) const
          {
            if (comp_real_only)
              return inexact_greater_than(d);
            else
              return exact_greater_than(d, 0);
          }
          bool operator>(const dual<data_type, comp_real_only> &d) const
          {
            if (comp_real_only)
              return inexact_greater_than(d.real());
            else
              return exact_greater_than(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool operator>(const dual<data2__, comp_real_only> &d) const
          {
            if (comp_real_only)
              return inexact_greater_than(d.real());
            else
              return exact_greater_than(d.real(), d.nonreal());
          }
          template <typename data2__>
          bool operator>(const dual_number::constant<data2__> &c) const
          {
            return operator>(c.real());
          }
          template <typename data2__>
          bool operator>(const dual_number::expression<data2__, comp_real_only> &e) const
          {
            if (comp_real_only)
              return inexact_greater_than(e.real());
            else
              return exact_greater_than(e.real(), e.nonreal());
          }

          void print(std::ostream &str) const
          {
            str << val[0] << " + " << val[1] << "ε";
          }

          // TODO: Implement this
          void input(std::istream &/*str*/)
          {
            assert(false);
          }

          // TODO: Implement this
          void write(std::ostream &/*str*/) const
          {
            assert(false);
          }

          // TODO: Implement this
          void read(std::istream &/*str*/)
          {
            assert(false);
          }
      };
    }
  }
}

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4244)
#endif

// comparison operators
template <typename data2__, typename data__, bool comp_real_only>
bool operator==(const data2__ &v, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d==v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator!=(const data2__ &v, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d!=v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator<=(const data2__ &v, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d>=v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator<(const data2__ &v, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d>v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator>=(const data2__ &v, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d<=v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator>(const data2__ &v, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d<v;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator==(const eli::mutil::ad::dual_number::constant<data2__> &c, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d==c;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator!=(const eli::mutil::ad::dual_number::constant<data2__> &c, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d!=c;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator<=(const eli::mutil::ad::dual_number::constant<data2__> &c, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d>=c;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator<(const eli::mutil::ad::dual_number::constant<data2__> &c, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d>c;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator>=(const eli::mutil::ad::dual_number::constant<data2__> &c, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d<=c;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator>(const eli::mutil::ad::dual_number::constant<data2__> &c, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d<c;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator==(const eli::mutil::ad::dual_number::expression<data2__, comp_real_only> &e, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d==e;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator!=(const eli::mutil::ad::dual_number::expression<data2__, comp_real_only> &e, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d!=e;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator<=(const eli::mutil::ad::dual_number::expression<data2__, comp_real_only> &e, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d>=e;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator<(const eli::mutil::ad::dual_number::expression<data2__, comp_real_only> &e, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d>e;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator>=(const eli::mutil::ad::dual_number::expression<data2__, comp_real_only> &e, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d<=e;
}
template <typename data2__, typename data__, bool comp_real_only>
bool operator>(const eli::mutil::ad::dual_number::expression<data2__, comp_real_only> &e, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  return d<e;
}
#ifdef _MSC_VER
# pragma warning(pop)
#endif

// I/O functions
template <typename data__, bool comp_real_only>
std::ostream & operator<<(std::ostream &ostr, const eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  d.print(ostr);
  return ostr;
}

template <typename data__, bool comp_real_only>
std::istream & operator>>(std::istream &istr, eli::mutil::ad::dual<data__, comp_real_only> &d)
{
  d.input(istr);
  return istr;
}
#endif
