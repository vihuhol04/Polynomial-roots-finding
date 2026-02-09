#ifndef INC_COMPLEXPRECISION
#define INC_COMPLEXPRECISION

/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2002-2024
 *                       Henrik Vestermark
 *                       Denmark
 *
 *                       All Rights Reserved
 *
 *   This source file is subject to the terms and conditions of the
 *   Henrik Vestermark License Agreement which restricts the manner
 *   in which it may be used.
 *   Mail: hve@hvks.com
 *
 *******************************************************************************
*/

/*
 *******************************************************************************
 *
 *
 * Module name     :   complexprecision.h
 * Module ID Nbr   :   
 * Description     :   Arbitrary complex precision class
 *                     Actually it a general complex class that works with both
 *                     standard types like int, float, double and float_precision
 *                     or int_precision or for that matter with any other types
 * --------------------------------------------------------------------------
 * Change Record   :   
 *
 * Version	Author/Date		Description of changes
 * -------  -----------		----------------------
 * 01.01	HVE/030331		Initial release
 * 01.02    hve/060203		Minor declaration bug fixed
 * 01.03    HVE/060217		Error in the formula for exp and log corrected
 * 01.04	HVE/14-Feb-2018	Change log(z). instead of using atan(imag()/real()) it has been replaced with atan2(imag(),real()) which is more correct
 * 01.05	HVE/05-Mar-2018	Added trigonometic functions for complex arguments
 * 01.06	HVE/07-Jul-2019	Included iostream header for increased portability
 * 01.07	HVE/22-Mar-2021 Updated License Info
 * 01.08    HVE/3-Oct-2023  Added conj() and norm() as function to be consistence with <complex>
 * 01.09    HVE-28-Mar-2024 added header file <stdexcept>
 * 01.10    HVE/26-Nov-2024 Fixed an error in cin >> cfp where the interim variable didnt carry the same precision as cfp
 * 01.11    HVE/9-Dec-2024  Added the precision() method for getting and setting the precision of the complex number
 * 
 * End of Change Record
 * --------------------------------------------------------------------------
*/

/* define version string */
static char _VC_[] = "@(#)complexprecision.h 01.11 -- Copyright (C) Henrik Vestermark";

#include "fprecision.h" // Required: complexprecision depends on float_precision, _PI, _float_table
#include <iostream>
#include <stdexcept>

//static_assert(__cplusplus >= 201703L, "The complexprecision.h code requires c++17 or higher.");

// Complex Precision template class
template<class _Ty> class complex_precision {
   _Ty re, im;
   public:
      typedef _Ty value_type;

      // Check that this class is only for float, double or float_precision
      static_assert(
          std::is_floating_point<_Ty>::value || std::is_same<_Ty, float_precision>::value,
          "_Ty must be one of the float types: float, double, long double or float_precision"
          );

      // constructor
      complex_precision( const _Ty& r = _Ty(0), const _Ty& i = _Ty(0)) : re(r), im(i) {}
      // constructor for any other type to _Ty
      template<class _X> complex_precision( const complex_precision<_X>& a ) : re((_Ty)a.real()), im((_Ty)a.imag()) {}
      
      // Coordinate functions
      _Ty real() const { return re; }
      _Ty imag() const { return im; }
      _Ty real( const _Ty& r )   { return ( re = r ); }
	  _Ty real( const complex_precision<_Ty>& r )   { return ( re = r.real() ); }
      _Ty imag( const _Ty& i )   { return ( im = i ); }
	  _Ty imag( const complex_precision<_Ty>& i )   { return ( im = i.imag() ); }
      _Ty norm() const { return re * re + im * im; }
      _Ty abs () const { _Ty a(re), b(im), c(re);  
                       a = (_Ty)re < _Ty(0) ? -re : re;
                       b = (_Ty)im < _Ty(0) ? -im : im;
                       if( a >= b )
                          {
                          c = im / re;
                          return a * sqrt( _Ty(1) + c * c );
                          }
                       else
                          {
                          c = re / im;
                          return b * sqrt( _Ty(1) + c * c );
                          }
                       }
      _Ty arg() const  { return atan2( im, re ); }
      complex_precision<_Ty> conj() const { complex_precision<_Ty> x(*this); x.real( re ); x.imag( -im ); return x; }
      _Ty *ref_real()   { return &re; }
      _Ty* ref_imag() { return &im; }
      size_t precision() const;			// Return the precision (number of decimal digits after the . )
      size_t precision(const size_t p);	// Set and return the new precision of the number
      std::string toString() {
          return "(" + re.toString() + "," + im.toString() + ")";
      }

      // Essential operators
      complex_precision<_Ty>& operator= ( const complex_precision<_Ty>& x )   { re = x.real(); im = x.imag(); return *this; }
      complex_precision<_Ty>& operator+=( const complex_precision<_Ty>& x )   { re += x.real(); im += x.imag(); return *this; }
      complex_precision<_Ty>& operator-=( const complex_precision<_Ty>& x )   { re -= x.real(); im -= x.imag(); return *this; }
      complex_precision<_Ty>& operator*=( const complex_precision<_Ty>& x )   { _Ty w(x.real()); w = re * x.real() - im * x.imag(); im = re * x.imag() + im * x.real(); re = w; return *this; }
      complex_precision<_Ty>& operator/=( const complex_precision<_Ty>& x );  // Too big to have here
     
      //class divide_by_zero {};
      class divide_by_zero : public std::logic_error {
      public:
          explicit divide_by_zero(const std::string& message)
              : std::logic_error(message) {}
      };
   };


/////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Input/output e.g. cout and cin
// 
////////////////////////////////////////////////////////////////////////////////////////////////

// The cin operator
//
template<class _Ty> std::ostream& operator<<( std::ostream& strm, const complex_precision<_Ty>& a )
	{ return strm << "(" << a.real() << "," << a.imag() << ")"; }

// The cout operator
template<class _Ty> std::istream& operator>>( std::istream& strm, complex_precision<_Ty>& c ) 
   {
   _Ty re(c.real()), im(c.imag());  // Ensure re & im is initialize with the same precision as c 
   char ch;
   re = _Ty(0);
   im = _Ty(0);
   // Parse the input
   if (strm >> ch)
   {
       if (ch == '(')
       {
           // Case: (re) or (re,im)
           if (strm >> re >> ch)
           {
               if (ch == ',')
               {
                   // Case: (re,im)
                   if (strm >> im >> ch && ch == ')')
                   {
                       c = complex_precision<_Ty>(re, im);
                       return strm;
                   }
               }
               else if (ch == ')')
               {
                   // Case: (re)
                   c = complex_precision<_Ty>(re, im);
                   return strm;
               }
           }
       }
       else
       {
           // Case: re (single number)
           strm.putback(ch);
           if (strm >> re)
           {
               c = complex_precision<_Ty>(re, im);
               return strm;
           }
       }
   }

   // If we reach here, input was invalid
   strm.setstate(std::ios::failbit);
   return strm;
   }

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Methods.
//  precision()     // Get precision
//  precision(prec) // Set precision
// 
// Other methods is done directly in the class definition
//
//////////////////////////////////////////////////////////////////////////////////////////////////

// return current precision of the left interval.
// default types like float and double just return the constant 9 and 17.
template<class _Ty> inline size_t complex_precision<_Ty>::precision() const
{
    if constexpr (std::is_same_v < _Ty, float>)
        return 9;
    else
        if constexpr (std::is_same_v < _Ty, double>)
            return 17;
        else
        {
            static_assert(std::is_same_v<_Ty, float_precision >, "Unsupported type for precision handling");
            return re.precision();
        }
    return 0;
}

// Set precision of both left and right interval. Only for float_precision type
// float and double types cant be changed.
template<class _Ty> inline size_t complex_precision<_Ty>::precision(const size_t p)
{
    if constexpr (std::is_same_v < _Ty, float_precision>)
    {	// Set both left ad right to the same precision
        re.precision(p);
        im.precision(p);
        return p;
    }
    else
        static_assert(std::is_same_v<_Ty, float_precision >, "Setting precision is only supported for float_precision type");

}


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Ovwerloading arithmetic & Boolean operators
//
//////////////////////////////////////////////////////////////////////////////////////////////////

// Arithmetic
template<class _Ty> complex_precision<_Ty> operator+( const complex_precision<_Ty>&, const complex_precision<_Ty>& );  // Binary
template<class _Ty> complex_precision<_Ty> operator+( const complex_precision<_Ty>& );                                 // Unary
template<class _Ty> complex_precision<_Ty> operator-( const complex_precision<_Ty>&, const complex_precision<_Ty>& );  // Binary
template<class _Ty> complex_precision<_Ty> operator-( const complex_precision<_Ty>& );                                 // Unary
template<class _Ty> complex_precision<_Ty> operator*( const complex_precision<_Ty>&, const complex_precision<_Ty>& );  // Binary
template<class _Ty> complex_precision<_Ty> operator/( const complex_precision<_Ty>&, const complex_precision<_Ty>& );  // Binary
                                                                                                                       
// Boolean Comparision Operators
template<class _Ty> bool operator==( const complex_precision<_Ty>&, const complex_precision<_Ty>& );
template<class _Ty> bool operator!=( const complex_precision<_Ty>&, const complex_precision<_Ty>& );


// Essential operators
//   a/=b
//
template<class _Ty> complex_precision<_Ty>& complex_precision<_Ty>::operator/=( const complex_precision<_Ty>& y )
   {
   if( y.real() == (_Ty)0 && y.imag() == (_Ty)0 )
      { throw complex_precision<_Ty>::divide_by_zero("Complex division by zero, is not allowed"); }

   if( ( y.real() < (_Ty)0 ? -y.real() : y.real() ) >= ( y.imag() < (_Ty)0 ? -y.imag() : y.imag() ) )
      {
      _Ty t(y.imag() / y.real() );        // Force same precision as y 
      _Ty t2(y.real() + y.imag() * t );   // Force same precision as y 
      _Ty t3(re + im * t );               // Force same precision as y 
      im -= re * t;
      im /= t2;
      re = t3 / t2;
      }
   else
      {
      _Ty t(y.real() / y.imag() );        // Force same precision as y 
      _Ty t2(y.real() * t + y.imag() );   // Force same precision as y 
      _Ty t3(re * t + im );               // Force same precision as y 
      im *= t;
      im -= re;
      im /= t2;
      re = t3 / t2;
      }

   return *this;
   }


// Result of a + b
//
template<class _Ty> complex_precision<_Ty> operator+( const complex_precision<_Ty>& x, const complex_precision<_Ty>& y )
   {
   complex_precision<_Ty> c(x);

   c += y;
   return c;
   }

// Unary + operator
//
template<class _Ty> complex_precision<_Ty> operator+( const complex_precision<_Ty>& x )
   {
   return x;
   }

// Result of a-b
//
template<class _Ty> complex_precision<_Ty> operator-( const complex_precision<_Ty>& x, const complex_precision<_Ty>& y )
   {
   complex_precision<_Ty> c(x);

   c -= y;
   return c;
   }

// Unary - operator
//
template<class _Ty> complex_precision<_Ty> operator-( const complex_precision<_Ty>& x )
   {
   complex_precision<_Ty> c(x);

   c.real( -c.real() );
   c.imag( -c.imag() );
   return c;
   }

// Result of a * b
//
template<class _Ty> complex_precision<_Ty> operator*( const complex_precision<_Ty>& x, const complex_precision<_Ty>& y )
   {
   complex_precision<_Ty> c(x);

   c *= y;
   return c;
   }

// result of a / b
//
template<class _Ty> complex_precision<_Ty> operator/( const complex_precision<_Ty>& x, const complex_precision<_Ty>& y )
   {
   complex_precision<_Ty> c(x);

   c /= y;
   return c;
   }

// Compare two complex number for equality
//
template<class _Ty> bool operator==( const complex_precision<_Ty>& a, const complex_precision<_Ty>& b )
   {
   return a.real() == b.real() && a.imag() == b.imag();
   }

// Compare two complex number for inequality
//
template<class _Ty> bool operator!=( const complex_precision<_Ty>& a, const complex_precision<_Ty>& b )
   {
   return a.real() != b.real() || a.imag() != b.imag();
   }


// Returns the absolute value of the complex number
//  As a function
template<class _Ty> _Ty abs(const complex_precision<_Ty>& x )
   {
   return x.abs();
   }

// Returns the norm of the complex number
//  As a function
template<class _Ty> _Ty norm(const complex_precision<_Ty>& x)
{
    return x.norm();
}

// Returns the complex conjugated conplex number
// As a function
template<class _Ty> complex_precision<_Ty> conj(const complex_precision<_Ty>& x)
{
    return x.conj();
}

// Conversion from Polar to Cartesian coordinates
//
template<class _Ty> complex_precision<_Ty> polarToComplex(_Ty magn, _Ty theta = (_Ty)0)
{
    if (theta == (_Ty)0)
        return complex_precision<_Ty>(magn);// sine cos(0)==1 and sin(0)==0

    return complex_precision<_Ty>(magn * cos(theta), magn * sin(theta));
}

// Conversion from Cartesian to Polar coordinates
//
template<class _Ty> std::pair<_Ty, _Ty> complexToPolar(_Ty x, _Ty y)
{
    _Ty magnitude(sqrt(x * x + y * y));
    _Ty angle(atan2(y, x)); // Angle in radians
    return std::make_pair(magnitude, angle);
}

// Calculate squareroot of a complex number
//
template<class _Ty> complex_precision<_Ty> sqrt( const complex_precision<_Ty> x )
   {
   _Ty w(x.real()), c(x.real()), d(x.real());  // Force the local variables to the same precisions as x (Only for float_precision)

   if( x.real() == _Ty(0) && x.imag() == _Ty(0) )
      w = _Ty(0);
   else
      {
      c = x.real() < _Ty(0) ? -x.real() : x.real();
      d = x.imag() < _Ty(0) ? -x.imag() : x.imag();
      if( c < d )
         {
         _Ty t(c / d);
         if( t < _Ty(0) ) t = -t;
         w = sqrt( d ) * sqrt( ( t + sqrt( _Ty(1) + t * t ) ) / _Ty(2) );
         }
      else
         {
         _Ty t(d / c);
         w = sqrt( c ) * sqrt( ( _Ty(1) + sqrt( _Ty(1) + t * t ) ) / _Ty(2) );
         }
      }

   if(  w == _Ty(0) )
      return complex_precision<_Ty>( w );
   
   if( x.real() >= _Ty(0) )
      return complex_precision<_Ty>( w, x.imag() / ( _Ty(2) * w ) );

   if( x.imag() >= _Ty(0) )
      return complex_precision<_Ty>( d / ( _Ty(2) * w ), w );

   return complex_precision<_Ty>( d / ( _Ty(2) * w ), -w );
   }

// Calculate log() of a complex number
//
template<class _Ty> complex_precision<_Ty> log( const complex_precision<_Ty> x )
   {
   return complex_precision<_Ty>( log( x.abs() ), atan2( x.imag(), x.real() ) );
   }

// Calcuate log10() of a complex number
//
template<class _Ty> complex_precision<_Ty> log10( const complex_precision<_Ty> x )
   {
   return log( x ) / complex_precision<_Ty>( log( _Ty(10) ), 0 );
   }

// Calcuate exp() of a complex number
//
template<class _Ty> complex_precision<_Ty> exp( const complex_precision<_Ty> x )
   {
   return complex_precision<_Ty>( exp( x.real() ) * cos( x.imag() ), exp( x.real() ) * sin( x.imag() ) );
   }

// Calculate pow() of a complex number raised to the power of a complex number
//
template<class _Ty> complex_precision<_Ty> pow( const complex_precision<_Ty> x, const complex_precision<_Ty> y )
   {
   complex_precision<_Ty> z(x);     // Force same precision as x
   z = log( x );
   z *= y;
   z = exp( z );
   return z;
   }

/****************************************************************************************************************
*
* Trigonometric functions for complex arguments
*	sin(x)
*	cos(x)
*	tan(x)
*	asin(x)
*	acos(x)
*	atan(x)
*
*****************************************************************************************************************
*/

// General template class for complex_precision<_Ty> sin
//
template<class _Ty> complex_precision<_Ty> sin(const complex_precision<_Ty> x)
{
	return complex_precision<_Ty>(sin(x.real())*cosh(x.imag()),cos(x.real()) * sinh(x.imag()));
}

// General template class for complex_precision<_Ty> cos
//
template<class _Ty> complex_precision<_Ty> cos(const complex_precision<_Ty> x)
{
	return complex_precision<_Ty>(cos(x.real())*cosh(x.imag()),-sin(x.real())*sinh(x.imag()));
}

// General template class for complex_precision<_Ty> tan
//
template<class _Ty> complex_precision<_Ty> tan(const complex_precision<_Ty> x)
{
	_Ty z(x.real());  // Force same precision
	
	z=cos((_Ty)2 * x.real()) + cosh((_Ty)2 * x.imag());
	
	return complex_precision<_Ty>(sin(2*x.real())/z,sinh(2*x.imag())/z);
}

// General template class for complex_precision<_Ty> asin
//
template<class _Ty> complex_precision<_Ty> asin(const complex_precision<_Ty> x)
{
	complex_precision<_Ty> z(x);		// Force same precision
	z = sqrt(complex_precision<_Ty>(1, 0) - x*x);
	z += x * complex_precision<_Ty>(0, 1);
	z = log(z);
	return complex_precision<_Ty>(z * complex_precision<_Ty>(0,-1));
}

// Specialization for acos for float_precision
//
inline complex_precision<float_precision> acos(const complex_precision<float_precision> x)
{
	complex_precision<float_precision> z(x);		// Force same precision
	z = sqrt(complex_precision<float_precision>(1, 0) - x*x);
	z += x * complex_precision<float_precision>(0, 1);
	z = log(z)* complex_precision<float_precision>(0, 1);
	return complex_precision<float_precision>(z + complex_precision<float_precision>(_float_table(_PI, (x.real()).precision())/(float_precision)2, 0));
}

// General template class for complex_precision<_Ty> acos
//
template<class _Ty> complex_precision<_Ty> acos(const complex_precision<_Ty> x)
{
	complex_precision<_Ty> z(x);		// Force same precision
	z = sqrt(complex_precision<_Ty>(1, 0) - x*x);
	z += x * complex_precision<_Ty>(0, 1);
	z = log(z)* complex_precision<_Ty>(0,1);
	return complex_precision<_Ty>( z + complex_precision<_Ty>(3.14159265358979323846/2, 0));
}

// General template class for complex_precision<_Ty> atan
//
template<class _Ty> complex_precision<_Ty> atan(const complex_precision<_Ty> x)
{
	complex_precision<_Ty> z(x), z1(x), z2(x);		// Force same precision
	z = x*complex_precision<_Ty>(0, 1);
	z1 = log(complex_precision<_Ty>(1,0)-z);
	z2 = log(z + complex_precision<_Ty>(1,0));
	z = z1 - z2;
	return complex_precision<_Ty>(z * complex_precision<_Ty>(0, 0.5));
}


/****************************************************************************************************************
*
* Hyperbolic functions for complex arguments
*	sinh(x)
*	cosh(x)
*	tanh(x)
*	asinh(x)
*	acosh(x)
*	atanh(x)
*
****************************************************************************************************************
*/

// General template class for complex_precision<_Ty> sinh
//
template<class _Ty> complex_precision<_Ty> sinh(const complex_precision<_Ty> x)
	{
	return complex_precision<_Ty>(sinh(x.real())*cos(x.imag()), cosh(x.real()) * sin(x.imag()));
	}

// General template class for complex_precision<_Ty> cosh
//
template<class _Ty> complex_precision<_Ty> cosh(const complex_precision<_Ty> x)
	{
	return complex_precision<_Ty>(cosh(x.real())*cos(x.imag()), sinh(x.real())*sin(x.imag()));
	}

// General template class for complex_precision<_Ty> tanh
//
template<class _Ty> complex_precision<_Ty> tanh(const complex_precision<_Ty> x)
	{
	_Ty z(x.real());  // Force same precision

	z = cosh(_Ty(2) * x.real()) + cos(_Ty(2) * x.imag());

	return complex_precision<_Ty>(sinh(2 * x.real()) / z, sin(2 * x.imag()) / z);
	}

// General template class for complex_precision<_Ty> asinh
//
template<class _Ty> complex_precision<_Ty> asinh(const complex_precision<_Ty> x)
	{
	complex_precision<_Ty> z(x);		// Force same precision
	z = sqrt(complex_precision<_Ty>(1, 0) + x*x);
	z += x;
	z = log(z);
	return complex_precision<_Ty>(z);
	}


// General template class for complex_precision<_Ty> acosh
//
//var return new Complex.log(Complex.add(x, Complex.mul(s1, s2))); }
template<class _Ty> complex_precision<_Ty> acosh(const complex_precision<_Ty> x)
	{
	complex_precision<_Ty> z1(x), z2(x);		// Force same precision
	
	z1 = sqrt(x - complex_precision<_Ty>(1, 0));
	z2 = sqrt(x + complex_precision<_Ty>(1, 0));

	return complex_precision<_Ty>( log(x + z1 * z2 ) );
	}

// General template class for complex_precision<_Ty> atanh
//
template<class _Ty> complex_precision<_Ty> atanh(const complex_precision<_Ty> x)
	{
	complex_precision<_Ty> z1(x), z2(x);		// Force same precision
	
	if (abs(x.real()) == (_Ty)1 && x.imag() == _Ty(0))
		return complex_precision<_Ty>(x.real()*_Ty(0), 0);  // Not correct

	z1 = log( complex_precision<_Ty>(1, 0) - x );
	z2 = log( complex_precision<_Ty>(1, 0) + x );
	return complex_precision<_Ty>(( z2 - z1 ) * complex_precision<_Ty>(0.5, 0 ));
	}

#endif
