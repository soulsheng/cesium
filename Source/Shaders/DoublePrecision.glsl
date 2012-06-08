/**
 * extPrecision.cg -- a set of Cg routines for doing extended precision
 *   floating point operations using double-floats stored in vec2 variables 
 *   (for real) or vec4 (for complex), each with 32-bit (24-bit mantissa)
 *   to give ~48 bits with the single-precision float exponent
 *
 * A. Thall
 * Dec. 6, 2004
 * 
 * This code copyright 2004, 2006, 2010  Andrew Thall
 * 
 * This code is offered as is for experimental and research purposes.  Its correctness
 *   is highly dependent on the nuances of the floating point behavior on any particular
 *   graphics card.  To understand the issues involved, see the Hida paper on extended
 *   precision for CPU code.  This author (Thall) makes no guarantees that any of this
 *   code will work as advertised.  Feel free to swipe from this code, but if you copy
 *   verbatim, I'd appreciate you crediting me with my work and would appreciate also
 *   hearing about your application!
 *
 * Many of these routines are based on similar routines in the QD library of Hida et al.,
 *   themselves derived from the doubledouble library of Briggs.  The author would also
 *   like to express thanks to Dr. David Bremer (LLNL) and Dr. Jeremy Meredith (ORNL)
 *   for discussions of their experiments in extended precision on GPUs.
 * The author of the current code asserts no ownership of algorithms presented herein,
 *   and has tried to cite attributions where appropropriate.
 * See the Hida papers on the QD library for more information on extended precision in general.
 * External documentation and information on extended precision for GPUs is available as
 *    a technical report http://andrewthall.org/papers/df64_qf128.pdf .
 *
 * If you're interested in extended precision for modern GPUs (2010) under CUDA and such,
 *    you should look at a more modern treatment such as "Supporting Extended Precision on
 *    Graphics Processors" by M. Lu, 2010.  This discusses double-doubles and quad-doubles
 *    implemented using CUDA, which has much better IEEE compliance and programmability,
 *    compared to my poor double-floats in Cg.
 */

const vec2 agi_df64E = vec2(2.718281745910645, 8.254840366817007e-008);
const vec2 agi_df64Log2 = vec2(0.6931471824645996, -1.904654212125934e-009);
const vec2 agi_df64Pi = vec2(3.141592741012573, -8.742277657347586e-008);
const vec2 agi_df64TwoPi = vec2(6.283185482025147, -1.748455531469517e-007);
const vec2 agi_df64PiOverTwo = vec2(1.570796370506287, -4.371138828673793e-008);
const vec2 agi_df64PiOverSixteen = vec2(0.1963495463132858, -5.463923535842241e-009);
const vec4 agi_df64SinCosOfPiOverSixteen = vec4(0.1950903236865997, -1.670471538872675e-009, 0.9807852506637573, 2.973947310636049e-008);
const vec4 agi_df64SinCosOfTwoPiOverSixteen = vec4(0.3826834261417389, 6.223350723644217e-009, 0.9238795042037964, 2.830748968563057e-008);
const vec4 agi_df64SinCosOfThreePiOverSixteen = vec4(0.5555702447891235, -1.176952135750753e-008, 0.8314695954322815, 1.687026340846387e-008);
const vec4 agi_df64SinCosOfFourPiOverSixteen = vec4(0.7071067690849304, 1.210161748588234e-008, 0.7071067690849304, 1.210161748588234e-008);

/**
 * the following splits a 24-bit IEEE floating point mantissa+E
 * into two numbers hi and low such that a = hi + low.  hi
 * will contain the first 12 bits, and low will contain the lower
 * order 12 bits.
 */

vec2 agi_split(float a) 
{
    float SPLIT = 4097.0; // (1 << 12) + 1;
    float t = a*SPLIT;
    float a_hi = t - (t - a);
    float a_lo = a - a_hi;
    
    return vec2(a_hi, a_lo);
}

/**
 * simulateously split the real and imaginary component
 * of a pair of floats into two
 */
vec4 agi_splitComplex(vec2 c) 
{
    float SPLIT = 4097.0; // (1 << 12) + 1;
    vec2 t = c*SPLIT; 
    vec2 c_hi = t - (t - c);
    vec2 c_lo = c - c_hi;
    return vec4(c_hi.x, c_lo.x, c_hi.y, c_lo.y);
}

vec2 agi_quickTwoSum(float a, float b) 
{
    float s = a + b;
    float e = b - (s - a);
    return vec2(s, e);
}

/**
 * does a quick sum of the high-order real and imaginary
 * components of a_ri and b_ri
 */
vec4 agi_quickTwoSumComplex(vec2 a_ri, vec2 b_ri) 
{

    vec2 s = a_ri + b_ri;
    vec2 e = b_ri - (s - a_ri);
    return vec4(s.x, e.x, s.y, e.y);
}

vec2 agi_twoSum(float a, float b) 
{
    float s = a + b;
    float v = s - a;
    float e = (a - (s - v)) + (b - v);
    return vec2(s, e);
}

vec4 agi_twoSumComplex(vec2 a_ri, vec2 b_ri) 
{
    vec2 s = a_ri + b_ri;
    vec2 v = s - a_ri;
    // xxAT Should be a 1.0* here on s?
    vec2 e = (a_ri - (s - v)) + (b_ri - v);
    return vec4(s.x, e.x, s.y, e.y);
}

vec2 agi_twoSubtract(float a, float b) 
{
    float s = a - b;
    float v = s - a;
    float err = (a - (s - v)) - (b + v);
    return vec2(s, err);
}

vec4 agi_twoSubtractComplex(vec2 a_ri, vec2 b_ri) 
{
    vec2 s = a_ri - b_ri;
    vec2 v = s - a_ri;
    vec2 err = (a_ri - (s - v)) - (b_ri + v);
    return vec4(s.x, err.x, s.y, err.y);
}

vec2 agi_twoProduct(float a, float b) 
{
    float p = a*b;
    vec2 aS = agi_split(a);
    vec2 bS = agi_split(b);
    float err = ((aS.x*bS.x - p) + aS.x*bS.y + aS.y*bS.x) + aS.y*bS.y;
    return vec2(p, err);
}

/**
 * faster agi_twoProduct that uses single call of split to complex version
 *   (which should run same speed as simple
 */
vec2 agi_twoProductFast(vec2 ab) 
{
    float p = ab.x*ab.y;
    vec4 S = agi_splitComplex(ab);
    float err = ((S.x*S.z - p) + S.x*S.w + S.y*S.z) + S.y*S.w;
    return vec2(p, err);
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
vec2 agi_twoSquare(float a) 
{
  float p = a*a;
  vec2 hilo = agi_split(a);
  float err = ((hilo.x*hilo.x - p) + 2.0*hilo.x*hilo.y) + hilo.y*hilo.y;
  return vec2(p, err);
}

/**
 * Computes fl(ar*ar), err(ar, ar), and fl(ai, ai) and err(ai*ai) for each component
 *    of a complex float (ar, ai)
 */
vec4 agi_twoSquareComplex(vec2 a_ri) 
{
    vec2 q = a_ri*a_ri;
    vec4 hilo = agi_splitComplex(a_ri);
    vec2 err = ((hilo.xz*hilo.xz - q) + 2.0*hilo.xz*hilo.yw) + hilo.yw*hilo.yw; 
    return vec4(q.x, err.x, q.y, err.y);
}
 
/* double-float * double-float */
vec2 agi_df64Multiply(vec2 a, vec2 b) 
{
    vec2 p;
    p = agi_twoProduct(a.x, b.x);
    p.y += a.x * b.y;
    p.y += a.y * b.x;
    p = agi_quickTwoSum(p.x, p.y);
    return p;
    /*
// Faster agi_twoProduct()
    p = agi_twoProductFast(vec2(a.x, b.x));
    p.y += dot(a, b.yx);    
    p = agi_quickTwoSum(p.x, p.y);
    return p;
    */
}

/* double-float * float */
vec2 agi_df64Multiply(vec2 a, float b) 
{
    vec2 p;
    p = agi_twoProduct(a.x, b);
//  p = agi_twoProductFast(vec2(a.x, b));
    p.y += (a.y * b);
    p = agi_quickTwoSum(p.x, p.y);
    return p;
}

vec2 agi_df64Multiply(float a, vec2 b) 
{
    return agi_df64Multiply(b, a);
}

/**
 * since a 2x mult just changes the exponent
 *   this should work for mult by any power of 2
 */
vec2 agi_df64MultFastX2(vec2 a) 
{
    return a*2.0;
}

vec2 agi_df64DivideFastX2(vec2 a) 
{
    return a * 0.5;
}

vec2 agi_df64MultiplyPower2(vec2 a, int b) 
{
    return a*float(b);
}

vec2 agi_df64DividePower2(vec2 a, int b) 
{
    return a/float(b);
}

/**
 * df64 + float
 */
vec2 agi_df64Add(vec2 a, float b) 
{
    vec2 s;
    s = agi_twoSum(a.x, b);
    s.y += a.y;
    s = agi_quickTwoSum(s.x, s.y);
    return s;
}

/**
 * float + df64
 */
vec2 agi_df64Add(float a, vec2 b) 
{
    vec2 s;
    s = agi_twoSum(b.x, a);
    s.y += b.y;
    s = agi_quickTwoSum(s.x, s.y);
    return s;
}

/** same as above but uses agi_twoSumComplex() to perform both
 *   twoSum() ops at the same time
 */
vec2 agi_df64Add(vec2 a, vec2 b) 
{
    vec4 st;
    st = agi_twoSumComplex(a, b);
    st.y += st.z;
    st.xy = agi_quickTwoSum(st.x, st.y);
    st.y += st.w;
    st.xy = agi_quickTwoSum(st.x, st.y);
    return st.xy;
}

/**
 * dfReal - float
 */
vec2 agi_df64Subtract(vec2 a, float b) 
{
    vec2 s;
    s = agi_twoSubtract(a.x, b);
    s.y += a.y;
    s = agi_quickTwoSum(s.x, s.y);
    return s;
}

/**
 * float - dfReal
 */
vec2 agi_df64Subtract(float a, vec2 b) 
{
    vec2 s;
    s = agi_twoSubtract(a, b.x);
    s.y -= b.y;
    s = agi_quickTwoSum(s.x, s.y);
    return s;
}

/**
 * the above _diff method can be improved using the agi_twoSubtractComplex()
 * to do the two agi_twoSubtract() calls simultaneously
 */
vec2 agi_df64Subtract(vec2 a, vec2 b) 
{
    vec4 st;
    st = agi_twoSubtractComplex(a, b);
    st.y += st.z;
    st.xy = agi_quickTwoSum(st.x, st.y);
    st.y += st.w;
    st.xy = agi_quickTwoSum(st.x, st.y);
    return st.xy;
}

/**
 * eq() & neq() -- equality tests
 */
// dfReal == float
bool agi_df64Equal(vec2 a, float b) 
{
  return (a.x == b && a.y == 0.0);
}

// float == dfReal
bool agi_df64Equal(float a, vec2 b) 
{
  return (a == b.x && b.y == 0.0);
}

/* double-float == double-float */
bool agi_df64Equal(vec2 a, vec2 b) 
{
  return (a.x == b.x && a.y == b.y);
}

// dfReal != float
bool agi_df64NotEqual(vec2 a, float b) 
{
  return (a.x != b || a.y != 0.0);
}

// float != dfReal
bool agi_df64NotEqual(float a, vec2 b) 
{
  return (a != b.x || b.y != 0.0);
}

/* double-float != double-float */
bool agi_df64NotEqual(vec2 a, vec2 b) 
{
  return (a.x != b.x || a.y != b.y);
}

/**
 * lt() & leq() -- less-than tests
 */
// dfReal < float
bool agi_df64LessThan(vec2 a, float b) 
{
  return (a.x < b || (a.x == b && a.y < 0.0));
}

// float < dfReal
bool agi_df64LessThan(float a, vec2 b) 
{
  return (a < b.x || (a == b.x && b.y > 0.0));
}

/* double-float < double-float */
bool agi_df64LessThan(vec2 a, vec2 b) 
{
  return (a.x < b.x || (a.x == b.x && a.y < b.y));
}

// dfReal <= float
bool agi_df64LessThanOrEqual(vec2 a, float b) 
{
  return (a.x < b || (a.x == b && a.y <= 0.0));
}

// float <= dfReal
bool agi_df64LessThanOrEqual(float a, vec2 b) 
{
    return (a < b.x || (a == b.x && b.y >= 0.0));
}

// dfReal <= dfReal
bool agi_df64LessThanOrEqual(vec2 a, vec2 b) 
{
  return (a.x < b.x || (a.x == b.x && a.y <= b.y));
}

/**
 * gt() & geq() -- greater-than tests
 */
// dfReal > float
bool agi_df64GreaterThan(vec2 a, float b) 
{
  return (a.x > b || (a.x == b && a.y > 0.0));
}

// float > dfReal
bool agi_df64GreaterThan(float a, vec2 b) 
{
  return (a > b.x || (a == b.x && b.y < 0.0));
}

/* double-float > double-float */
bool agi_df64GreaterThan(vec2 a, vec2 b) 
{
  return (a.x > b.x || (a.x == b.x && a.y > b.y));
}

// dfReal >= float
bool agi_df64GreaterThanOrEqual(vec2 a, float b) 
{
  return (a.x > b || (a.x == b && a.y >= 0.0));
}

// float >= dfReal
bool agi_df64GreaterThanOrEqual(float a, vec2 b) 
{
    return (a > b.x || (a == b.x && b.y <= 0.0));
}

// dfReal >= dfReal
bool agi_df64GreaterThanOrEqual(vec2 a, vec2 b) 
{
  return (a.x > b.x || (a.x == b.x && a.y >= b.y));
}

/*********** Squaring and higher powers **********/
vec2 agi_df64Square(vec2 a) 
{
  vec2 p;
  vec2 s;
  p = agi_twoSquare(a.x);
  p.y += 2.0 * a.x * a.y;
  p.y += a.y * a.y;
  s = agi_quickTwoSum(p.x, p.y);
  return s;
}

vec2 agi_df64Square(float a) 
{
    vec2 p1;
    p1 = agi_twoSquare(a);
    return agi_quickTwoSum(p1.x, p1.y);
}

/**
 * sqrt() -- this uses "Karp's Trick" according to Hida,
 *   which is just to use the single-precision sqrt as
 *   an approximation and does a Newton iteration using
 *   as few full-precision operations as possible.  See
 *   [Karp 93]
 *
 * "sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)
 *
 *   The approximation is accurate to twice the accuracy of x.
 *   Also, the multiplication (a*x) and [-]*x can be done with
 *   only half the precision."
 * 
 * @param A must be non-negative
 */
vec2 agi_df64Sqrt(vec2 A) 
{
//  float xn = (A.x == 0.0) ? 0.0 : rsqrt(A.x);

    float xn = inversesqrt(A.x);
    float yn = A.x*xn;
    vec2 ynsqr = agi_df64Square(yn);

    float diff = (agi_df64Subtract(A, ynsqr)).x;
    vec2 prodTerm = agi_twoProduct(xn, diff) * 0.5;
    
    return agi_df64Add(yn, prodTerm);
}

/**
 * div() -- this also uses "Karp's Trick", using the
 *   single-precision recip() as an approximation.
 *   Hida's work was hard to decipher, due to operator
 *   overloading. Will implement their version later
 *   and compare it to this for speed and accuracy.
 *   For details, see [Karp 93].
 *
 *   xn = recip(A)
 *   yn = B*xn
 *   "div(B, A) = yn + xn*(B - A*yn)   (approx)
 *
 *   Since Newton's is quadratic convergent, should get double
 *      prec. from single-prec. approx
 *
 *      xn is single-prec.
 *      yn is single-prec. approx by prod. of two single-prec.
 *      A*yn is double*single
 *      B - Ayn is double-double
 *      xn*(...) multiplies two single-prec and gives double.
 * 
 * @param A must be non-zero
 */
vec2 agi_df64Divide(vec2 B, vec2 A) 
{
    float xn = 1.0/A.x;

    float yn = B.x*xn;
    float diffTerm = (agi_df64Subtract(B, agi_df64Multiply(A, yn))).x;
    vec2 prodTerm = agi_twoProduct(xn, diffTerm);

    return agi_df64Add(yn, prodTerm);
}

/**
 * similar to above, but A only specified to f32 precision
 */
vec2 agi_df64Divide(vec2 B, float A) 
{
    float xn = 1.0/A;
    float yn = B.x*xn;
    float diffTerm = (agi_df64Subtract(B, agi_twoProduct(A, yn))).x;
    vec2 prodTerm = agi_twoProduct(xn, diffTerm);
    
    return agi_df64Add(yn, prodTerm);

}

/**
 * similar to above, but B only specified to f32 precision
 */
vec2 agi_df64Divide(float B, vec2 A) 
{
    float xn = 1.0/A.x;
    float yn = B*xn;
    float diffTerm = (agi_df64Subtract(B, agi_df64Multiply(A, yn))).x;
    vec2 prodTerm = agi_twoProduct(xn, diffTerm);
    
    return agi_df64Add(yn, prodTerm);

}

/**
 * both values are only f32
 */
 vec2 agi_df64Divide(float B, float A) 
 {
    float xn = 1.0/A;
    float yn = B*xn;
    float diffTerm = (agi_df64Subtract(B, agi_twoProduct(A, yn))).x;
    vec2 prodTerm = agi_twoProduct(xn, diffTerm);
    
    return agi_df64Add(yn, prodTerm);
}

/**
 * fabs is an easy one for df64 on the GPU
 *   --multiply a by sign of high-order float
 */
vec2 agi_df64FAbs(vec2 a) 
{
    return sign(a.x)*a;
}

/**
 * if hi-float is integer already, add floor(lo)
 *    else lo = 0
 */
vec2 agi_df64Floor(vec2 a) 
{
    /**
 * This will be faster on a vector processor, the below on a scalar
 *
    vec2 outVal = floor(a);
    if (outVal.x == a.x)
        outVal = agi_quickTwoSum(outVal.x, outVal.y);
    else
        outVal.y = 0;
 */
    vec2 outVal;
    outVal.x = floor(a.x);
    if (outVal.x != a.x)
        outVal.y = 0.0;
    else
        outVal = agi_quickTwoSum(a.x, floor(a.y));
    return outVal;
}


/**
 * if hi-float is integer already, add ceil(lo)
 *    else lo = 0
 */
vec2 agi_df64Ceil(vec2 a) 
{
/**
 * This will be faster on a vector processor, the below on a scalar
 *
    vec2 outVal = ceil(a);
    if (outVal.x == a.x)
        outVal = agi_quickTwoSum(outVal.x, outVal.y);
    else
        outVal.y = 0;
 */
    vec2 outVal;
    outVal.x = ceil(a.x);
    if (outVal.x != a.x)
        outVal.y = 0.0;
    else
        outVal = agi_quickTwoSum(a.x, ceil(a.y));
    return outVal;
}

/**
 * nint() -- nearest int rounding code.  Adapted from Yozo Hida's code.
 *    This seems simple and easy, but may be inaccurate for values
 *    very near i + 0.5.
 */
 vec2 agi_df64NearestInt(vec2 a) 
 {
    return agi_df64Floor(agi_df64Add(a, 0.5));
 }

/**
 * df64_rem -- compute quotient and rem satisfying 
 *     a = quot*b + rem,  |rem| <= b/2,
 * @return rem, save quot in out var
 * Note that this does *not* return fmod(a, b),
 *     since |quot*b| may be > |a|
 * NOTE:  this needs to be done in extended precision;
 *  --it assumes values are exact, but, as below, we need higher
 *    than df64 for our M_PI values to get correct product and
 *    remainders for the trig. range-reductions.
 */

// Like Miller's modr
// a=n*b+rem, |rem|<=b/2, exact result.
vec2 agi_df64Remainder(vec2 a, vec2 b, out vec2 quot) 
{
  vec2 temp;
  temp = agi_df64Divide(a, b);
  //float n = round(temp.x); //TODO What kind of rounding?
  float n = floor(temp.x + 0.5);
  temp = agi_twoProduct(n, b.x);
  vec2 rem = vec2(a.x, 0.0);
  temp = agi_df64Subtract(rem, temp);
  rem = vec2(a.y, 0.0);
  temp = agi_df64Add(rem, temp);
  rem = agi_df64Multiply(n, vec2(b.y, 0.0));
  rem = agi_df64Subtract(temp, rem);
  quot = vec2(n, 0.0);
  
  return rem;
}

/**
 * pre: a > 0, n integer, there may be something wrong here xxAT
 */
vec2 agi_df64NPow(vec2 a, int n) 
{
    /* if n = 0, return 1 unless a = 0, when Nan */
    /* for now xxAT don't bother checking a */
    vec2 outVal;
    int nmag = n >= 0 ? n : -n;

    vec2 r = a;
    
    if (int(mod(float(nmag), 2.0)) == 1)
        outVal = a;
    else
        outVal = vec2(1.0, 0.0);
        
    nmag = nmag/2;
    
    // TODO: Is this enough iterations?
    for (int i = 0; i < 1000; ++i)
    {
        if (nmag <= 0)
            break;
        
        r = agi_df64Square(r); 
        if (int(mod(float(nmag), 2.0)) == 1)
            outVal = agi_df64Multiply(outVal, r);
        nmag = nmag/2;
    }

    if (n < 0)
        outVal = agi_df64Divide(vec2(1.0, 0.0), outVal);

    return outVal;
}   
    
vec2 agi_df64NPow(float a, int n) 
{
    return agi_df64NPow(vec2(a, 0.0), n);
}

/**
 * agi_df64NPow2To6() -- compute a^64 = a^(2^6) power
 */
vec2 agi_df64NPow2To6(vec2 a) 
{
    vec2 outVal = a;
    
    outVal = agi_df64Square(outVal);
    outVal = agi_df64Square(outVal);
    outVal = agi_df64Square(outVal);
    outVal = agi_df64Square(outVal);
    outVal = agi_df64Square(outVal);
    outVal = agi_df64Square(outVal);

    return outVal;
}

/**
 * agi_df64NPow4() -- compute a^4
 */
vec2 agi_df64NPow4(vec2 a) 
{
    vec2 outVal = a;
    outVal = agi_df64Square(outVal);
    outVal = agi_df64Square(outVal);

    return outVal;
}

/**
 * expTAYLOR() -- compute Taylor series approx to e^a
 *   for |a| < 0.000542. [(ln 2)/2/64].  Should converge
 *   very rapidly.
 */
vec2 agi_df64ExpTaylor(vec2 a) 
{
    float thresh = 1.0e-20*exp(a.x);
    vec2 t;  /* Term being added. */
    vec2 p;  /* Current power of a. */
    vec2 f;  /* Denominator. */
    vec2 s;  /* Current partial sum. */
    vec2 x;  /* = -sqr(a) */
    float m;
    
    s = agi_df64Add(1.0, a);  // first two terms
    p = agi_df64Square(a);
    m = 2.0;
    f = vec2(2.0, 0.0);
    t = p/2.0;
    
    // TODO: Is this enough iterations?
    for (int i = 0; i < 1000; ++i)
    {
        if (abs(t.x) < thresh)
            break;
        
        s = agi_df64Add(s, t);
        p = agi_df64Multiply(p, a);
        m += 1.0;
        f = agi_df64Multiply(f, m);
        t = agi_df64Divide(p, f);
    }

    return agi_df64Add(s, t);
}

vec2 agi_df64Exp(vec2 a) 
{
  /* (by Hida):
     Strategy:  We first reduce the size of x by noting that
     
          exp(kr + m) = exp(m) * exp(r)^k

     Thus by choosing m to be a multiple of log(2) closest
     to x, we can make |kr| <= log(2) / 2 = 0.3466.  Now
     we can set k = 64, so that |r| <= 0.000542.  Then

          exp(x) = exp(kr + s log 2) = (2^s) * [exp(r)]^64

     Then exp(r) is evaluated using the familiar Taylor series.
     Reducing the argument substantially speeds up the convergence.
  */  
    vec2 outVal = vec2(0.0, 0.0);
    vec2 r;
    vec2 rem, df_z;
    
    // RANGE REDUCTION IS SCREWED --- HOW TO FIX??? xxAT
    int lgk = 2;//(int) 1;    // using 1, 2 gives NaNs, but better log(x)
    int k = 4; //(int) 2;
    /**
     * return +INF or 0 if exponent is out of range
     * (just let outVal.x = exp(a.hi)
     */
    
    if ((a.x >= 88.7) || (a.x <= -87.3)) 
//      outVal.x = exp(a.x);
        outVal.x = 0.0;
    else if (a.x == 0.0)
        outVal.x = 1.0;
    else if (agi_df64Equal(a, 1.0))
        outVal = agi_df64E;
    else if (abs(a.x) < 0.000542)
        outVal = agi_df64ExpTaylor(a);
    else {
        // ~150 fps
        // choice of test is 8: main_EXP_df64
        // EPS = 7.105427e-015 and numCount = 262144
        // For whole image, min error    = -18.486298 ulps
        //          max error    = 16.943582 ulps
        //          mean error   = -0.590102 ulps
        //          mean |error| = 2.199226 ulps
        //          r.m.s error  = 3.030591 ulps
        // with 0 NaNs or Infs
        r = a/float(k);
        r = agi_df64ExpTaylor(r);

        outVal = agi_df64NPow4(r);
    }
    /* else {
        // ~360 fps
        //  choice of test is 8: main_EXP_df64
        //  EPS = 7.105427e-015 and numCount = 262144
        //  For whole image, min error    = -90.356737 ulps
        //             max error    = 53.146445 ulps
        //             mean error   = -4.933327 ulps
        //             mean |error| = 13.192848 ulps
        //             r.m.s error  = 18.129237 ulps
        //  with 0 NaNs or Infs
 
        k = 64;
        lgk = 6;
        rem = agi_df64Remainder(a, agi_df64Log2, df_z);
        int z = (int) df_z.x;

        // since k is 2^6, can simply divide the df64 number pairwise
        r = rem/k;
        r = agi_df64ExpTaylor(r);
        r = df64_nTo2ToThe6(r);
        
        outVal = agi_df64Multiply(r, pow(2.0, z));
    }
    */
    return outVal;
}

vec2 agi_df64Exp(float a) {
    return agi_df64Exp(vec2(a, 0.0));
}

/**
 * As Hida et al., use a Newton iteration for log_e(x)
 */
vec2 agi_df64Log(vec2 a) 
{
      /* Strategy.  The Taylor series for log converges much more
     slowly than that of exp, due to the lack of the factorial
     term in the denominator.  Hence this routine instead tries
     to determine the root of the function

         f(x) = exp(x) - a

     using Newton iteration.  The iteration is given by

         x' = x - f(x)/f'(x) 
            = x - (1 - a * exp(-x))
            = x + a * exp(-x) - 1.
       
     Only one iteration is needed, since Newton's iteration
     approximately doubles the number of digits per iteration. */

    vec2 xi = vec2(0.0, 0.0);
    
    if (!agi_df64Equal(a, 1.0)) {
        if (a.x <= 0.0)
            xi = vec2(log(a.x));   // return NaN
        else {
            xi.x = log(a.x);
            xi = agi_df64Add(agi_df64Add(xi, agi_df64Multiply(agi_df64Exp(-xi), a)), -1.0);
        //  xi = agi_df64Add(agi_df64Add(xi, agi_df64Multiply(agi_df64Exp(-xi), a)), -1.0);
        //  xi = agi_df64Add(agi_df64Add(xi, agi_df64Multiply(agi_df64Exp(-xi), a)), -1.0);
        }
    }
    return xi;
}

/**
 * nthRoot() -- return the positive nth-root of input using
 *   Newton iteration
 * From Hida (xxAT),
       Strategy:  Use Newton's iteration to solve
          1/(x^n) - a = 0
       Newton iteration becomes
          x' = x + x * (1 - a * x^n) / n
       Since Newton's iteration converges quadratically, 
       we only need to perform it twice.

 * xxAT This can probably be done more efficiently as per agi_df64Sqrt().
 * @param x must be non-negative
 */
vec2 agi_df64NRoot(vec2 A, int n) 
{
    vec2 outVal;
    
    if (int(mod(float(n), 2.0)) == 0 && A.x <= 0.0)
        outVal = vec2(sqrt(-1.0), 0.0);
    else if (n == 0)
        outVal = vec2(1.0, 0.0);
    else if (n == 1)
        outVal = A;
    else if (n == 2)
        outVal = agi_df64Sqrt(A);
    else {
        vec2 r = agi_df64FAbs(A);
        vec2 x = agi_df64Exp(-log(r.x)/float(n));

        outVal = x;
        vec2 prodTerm = agi_df64Multiply(r, agi_df64NPow(x, n));
        vec2 sumExpr = agi_df64Divide(agi_df64Subtract(1.0, prodTerm), float(n));
        x = agi_df64Add(x, agi_df64Multiply(x, sumExpr));
        
        if (A.x < 0.0)
            x = -x;
        outVal = agi_df64Divide(1.0, x);
    }
    return outVal;      
}

/**
 * sincosTAYLOR(a) -- computes sin(a) == .xy
 *                             cos(a) == .zw
 *   using Taylor series, assuming fabs(a) < PI/32
 * @return vec4(df64:sin(a), df64:cos(a))
 */
vec4 agi_df64SinCosTaylor(vec2 a) {

    float thresh = 1.0e-20 * abs(a.x);
    vec2 t;  /* Term being added. */
    vec2 p;  /* Current power of a. */
    vec2 f;  /* Denominator. */
    vec2 s;  /* Current partial sum. */
    vec2 x;  /* = -sqr(a) */
    float m;

    vec2 sin_a, cos_a;
    if (a.x == 0.0) {
        sin_a = vec2(0.0, 0.0);
        cos_a = vec2(1.0, 0.0);
    }
    else {
        x = -agi_df64Square(a);
        s = a;
        p = a;
        m = 1.0;
        f = vec2(1.0, 0.0);
        
        // TODO: Is this enough iterations?
        for (int i = 0; i < 1000; ++i)
        {
            p = agi_df64Multiply(p, x);
            m += 2.0;
            f = agi_df64Multiply(f, m*(m-1.0));
            t = agi_df64Divide(p, f);
            s = agi_df64Add(s, t);
            
            if (abs(t.x) < thresh)
                break;
        }

        sin_a = s;
        cos_a = agi_df64Sqrt(agi_df64Add(1.0, -agi_df64Square(s)));
    }
    return vec4(sin_a, cos_a);
}

/**
 * sin(a) -- use Hida's method to do argument reduction and compute sin(x)
 * NOTE:  this really needs to be done in extended precision;
 *  --it assumes values are exact, but, as below, we need higher
 *    than df64 for our M_PI values to get correct product and
 *    remainders for the trig. range-reductions.
 *  Current values can be anywhere from great to really bad,
 *    whereas we'd like 0.5 ulps for our trig functions.
*/
vec2 agi_df64Sin(vec2 a) {  

    /* Strategy.  To compute sin(x), we choose integers a, b so that
    x = s + a * (pi/2) + b * (pi/16)
    and |s| <= pi/32.  Using the fact that 
    sin(pi/16) = 0.5 * sqrt(2 - sqrt(2 + sqrt(2)))
    we can compute sin(x) from sin(s), cos(s).  This greatly 
    increases the convergence of the sine Taylor series. 
    */
    vec2 outVal = vec2(0.0, 0.0);
    vec2 r, t, s, c, j, k;
    vec4 sin_t_cos_t;
    vec4 uv;
    
    if (a.x != 0.0) {
        /* First reduce modulo 2*pi so that |r| <= pi. */
        r = agi_df64Remainder(a, agi_df64TwoPi, j);  // xxAT better pass in 2*PI as const

        /* Now reduce by modulo pi/2 and then by pi/16 so that
        we obtain numbers a, b, and t. */
        t = agi_df64Remainder(r, agi_df64PiOverTwo, j);
        //float abs_j = round(abs(j.x)); // TODO: which rounding method?
        float abs_j = floor(abs(j.x) + 0.5);
        t = agi_df64Remainder(t, agi_df64PiOverSixteen, k);
        //float abs_k = round(abs(k.x)); // TODO: which rounding method?
        float abs_k = floor(abs(k.x) + 0.5);

        if (abs_j > 2.0 || abs_k > 4.0 )
            outVal = vec2(sqrt(-1.0), 0.0);
        else {

            sin_t_cos_t = agi_df64SinCosTaylor(t);

            if (abs_k == 0.0) {
                s = sin_t_cos_t.xy;
                c = sin_t_cos_t.zw;
            }
            else {
                if (abs_k == 1.0)
                    uv = agi_df64SinCosOfPiOverSixteen.zwxy;
                else if (abs_k == 2.0)
                    uv = agi_df64SinCosOfTwoPiOverSixteen.zwxy;
                else if (abs_k == 3.0)
                    uv = agi_df64SinCosOfThreePiOverSixteen.zwxy;
                else 
                    uv = agi_df64SinCosOfFourPiOverSixteen.zwxy;

                float signK = sign(k.x);
                
                s = agi_df64Add(agi_df64Multiply(uv.xy, sin_t_cos_t.xy),  signK*agi_df64Multiply(uv.zw, sin_t_cos_t.zw));
                c = agi_df64Add(agi_df64Multiply(uv.xy, sin_t_cos_t.zw), -signK*agi_df64Multiply(uv.zw, sin_t_cos_t.xy));
            }

            if (abs_j == 0.0)
                outVal = s;
            else if (j.x == 1.0)
                outVal = c;
            else if (j.x == -1.0)
                outVal = -c;
            else
                outVal = -s;
        }
    }
    return outVal;;
}

/**
 * compute both sine and cosine of a -- virtually identical to agi_df64Sin()
 */
vec4 agi_df64SinCos(vec2 a) {  
    /* Strategy.  To compute sin(x), we choose integers a, b so that
    x = s + a * (pi/2) + b * (pi/16)
    and |s| <= pi/32.  Using the fact that 
    sin(pi/16) = 0.5 * sqrt(2 - sqrt(2 + sqrt(2)))
    we can compute sin(x) from sin(s), cos(s).  This greatly 
    increases the convergence of the sine Taylor series. 
    */
    
    vec4 outVal = vec4(0.0, 0.0, 1.0, 0.0);
    vec2 r, t, s, c, j, k, dummy;
    vec4 sin_t_cos_t;
    vec4 uv;
    
    if (a.x != 0.0) {
        
        /* First reduce modulo 2*pi so that |r| <= pi. */
        r = agi_df64Remainder(a, agi_df64TwoPi, dummy);  // xxAT better pass in 2*PI as const

        /* Now reduce by modulo pi/2 and then by pi/16 so that
        we obtain numbers a, b, and t. */
        t = agi_df64Remainder(r, agi_df64PiOverTwo, j);
        //float abs_j = round(abs(j.x)); // TODO: which rounding method?
        float abs_j = floor(abs(j.x) + 0.5);
        t = agi_df64Remainder(t, agi_df64PiOverSixteen, k);
        //float abs_k = round(abs(k.x)); // TODO: which rounding method?
        float abs_k = floor(abs(k.x) + 0.5);

        if (abs_j > 2.0 || abs_k > 4.0 ) {
            float i = sqrt(-1.0);
            outVal = vec4(i, 0.0, i, 0.0);
        }
        else {
            sin_t_cos_t = agi_df64SinCosTaylor(t);

            if (abs_k == 0.0) {
                s = sin_t_cos_t.xy;
                c = sin_t_cos_t.zw;
            }
            else {
                if (abs_k == 1.0)
                    uv = agi_df64SinCosOfPiOverSixteen.zwxy;
                else if (abs_k == 2.0)
                    uv = agi_df64SinCosOfTwoPiOverSixteen.zwxy;
                else if (abs_k == 3.0)
                    uv = agi_df64SinCosOfThreePiOverSixteen.zwxy;
                else 
                    uv = agi_df64SinCosOfFourPiOverSixteen.zwxy;

                float signK = sign(k.x);
                
                s = agi_df64Add(agi_df64Multiply(uv.xy, sin_t_cos_t.xy),  signK*agi_df64Multiply(uv.zw, sin_t_cos_t.zw));
                c = agi_df64Add(agi_df64Multiply(uv.xy, sin_t_cos_t.zw), -signK*agi_df64Multiply(uv.zw, sin_t_cos_t.xy));
            }

            if (abs_j == 0.0)
                outVal = vec4(s, c);
            else if (j.x == 1.0)
                outVal = vec4(c, -s);
            else if (j.x == -1.0)
                outVal = vec4(-c, s);
            else
                outVal = vec4(-s, -c);
        }
    }
    
    return outVal;;
}

// Alternate sin function that uses the taylor series expansion of sine around zero.
vec2 agi_df64Sin2(vec2 a)
{
    vec2 coef1 = vec2(0.1666666716337204, -4.967053879312289e-009);
    vec2 coef2 = vec2(0.008333333767950535, -4.34617203337595e-010);
    vec2 coef3 = vec2(0.0001984127011382952, -2.725596874933456e-012);
    vec2 coef4 = vec2(2.755731884462875e-006, 3.793571224297229e-014);
    vec2 coef5 = vec2(2.50521079436794e-008, 4.417623044648367e-016);
    vec2 coef6 = vec2(1.605904437207428e-010, -5.352526511562726e-018);
    vec2 coef7 = vec2(7.647163609812713e-013, 1.220071047117829e-020);

    vec2 a2 = agi_df64Multiply(a, a);

    vec2 result = agi_df64Subtract(agi_df64Multiply(agi_df64Add(agi_df64Multiply(a2, -coef7), coef6), a2), coef5);
    result = agi_df64Subtract(agi_df64Multiply(agi_df64Add(agi_df64Multiply(result, a2), coef4), a2), coef3);
    result = agi_df64Subtract(agi_df64Multiply(agi_df64Add(agi_df64Multiply(result, a2), coef2), a2), coef1);
    result = agi_df64Multiply(agi_df64Add(agi_df64Multiply(result, a2), vec2(1.0, 0.0)), a);

    return result;
}

/**
 * cos(x) -- this is easy.  sincos() is almost identical in runtime to sin();
 *   any transformation of the sin's result will be more costly than simply
 *   running sinoos()
 */
vec2 agi_df64Cos(vec2 a) {
    return agi_df64SinCos(a).zw;
}

/**
 * Complex number routines, including
 *
 *    agi_cdf64Add(vec4, vec4)
 *    agi_cdf64Subtract(vec4, vec4)
 *    agi_cdf64Multiply(vec4, vec4)
 *    agi_cdf64ComplexConjugate(vec4)
 *    agi_cdf64Divide(vec4, vec4)
 *    agi_cdf64SquareComponentwise(vec4) -- returns componentwise sqr,
 *       i.e., sqrc(vec4(a, b)) = sqrc(a + bi) = (a*a + b*b*i)
 */

/**
 * The following code was ported to GLSL but hasn't been tested.
 
// does a df add of two complex df numbers,
//   where the .xy is the df real components 
//             .zw is the df complex components
vec4 agi_cdf64Add(vec4 a_ri, vec4 b_ri) 
{
    vec4 s, t;
    s = agi_twoSumComplex(a_ri.xz, b_ri.xz);
    t = agi_twoSumComplex(a_ri.yw, b_ri.yw);
    s.yw += t.xz;
    s = agi_quickTwoSumComplex(s.xz, s.yw);
    s.yw += t.yw;
    s = agi_quickTwoSumComplex(s.xz, s.yw);
    return s;
}

vec4 agi_cdf64Subtract(vec4 a_ri, vec4 b_ri) 
{
    vec4 s, t;
    s = agi_twoSubtractComplex(a_ri.xz, b_ri.xz);
    t = agi_twoSubtractComplex(a_ri.yw, b_ri.yw);
    s.yw += t.xz;
    s = agi_quickTwoSumComplex(s.xz, s.yw);
    s.yw += t.yw;
    s = agi_quickTwoSumComplex(s.xz, s.yw);
    return s;
}

// xxAT -- a true complex multiply.  Check for accuracy and efficiency.
vec4 agi_cdf64Multiply(vec4 a, vec4 b) 
{
    return vec4(agi_df64Subtract(agi_df64Multiply(a.xy, b.xy), agi_df64Multiply(a.zw, b.zw)),
                  agi_df64Add(agi_df64Multiply(a.xy, b.zw), agi_df64Multiply(a.zw, b.xy)));
}

// xxAT -- returns the complex conjugate
vec4 agi_cdf64ComplexConjugate(vec4 a) 
{
    return vec4(a.xy, -a.zw);
}

// squares the real and imaginary terms of a complex number
//   individually (i.e., (r, i) --> (r*r, i*i)
vec4 agi_cdf64SquareComponentwise(vec4 a) 
{
    vec4 p;
    vec4 s;
    p = agi_twoSquareComplex(a.xz);
    p.yw += 2.0*a.xz*a.yw;
    p.yw += a.yw*a.yw;
    s = agi_quickTwoSumComplex(p.xz, p.yw);
    return s;
}

// computes the square of the complex number (a + bi)
//   as (a^2 - b^2) + 2abi
vec4 agi_cdf64Square(vec4 a) 
{
    vec4 zed;
    vec4 aCompSqr = agi_cdf64SquareComponentwise(a);
    zed.xy = agi_df64Subtract(aCompSqr.xy, aCompSqr.zw);
    zed.zw = agi_df64Multiply(a.xy, a.zw);
    zed.zw *= 2.0;
    return zed;
}

// xxAT -- complex division
// @param A must be non-zero
vec4 agi_cdf64Divide(vec4 B, vec4 A) 
{
    vec4 subprod = agi_cdf64Multiply(B, agi_cdf64ComplexConjugate(A));
    vec4 normA = agi_cdf64SquareComponentwise(A);
    vec2 denom = agi_df64Add(normA.xy, normA.zw);

    return vec4(agi_df64Divide(subprod.xy, denom), agi_df64Divide(subprod.zw, denom));
}

// double-float * double-float
vec2 agi_dfRealMultiply(vec2 a, vec2 b) 
{
    vec2 p;

    //    p = agi_twoProduct(a.x, b.x);
    //    p.y += a.x * b.y;
    //    p.y += a.y * b.x;
    //    p = agi_quickTwoSum(p.x, p.y);
    //    return p;
    
// Faster agi_twoProduct()
    p = agi_twoProductFast(vec2(a.x, b.x));
    p.y += dot(a.xy, b.yx); 
    p = agi_quickTwoSum(p.x, p.y);
    return p;
}

vec2 agi_dfRealAdd(vec2 a, vec2 b) 
{
    vec4 st;
    st = agi_twoSumComplex(a, b);
    st.y += st.z;
    st.xy = agi_quickTwoSum(st.x, st.y);
    st.y += st.w;
    st.xy = agi_quickTwoSum(st.x, st.y);
    return st.xy;
}

vec2 agi_dfRealSubtract(vec2 a, vec2 b) 
{
    vec4 st;
    st = agi_twoSubtractComplex(a, b);
    st.y += st.z;
    st.xy = agi_quickTwoSum(st.x, st.y);
    st.y += st.w;
    st.xy = agi_quickTwoSum(st.x, st.y);
    return st.xy;
}

vec4 agi_dfCompAdd(vec4 a_ri, vec4 b_ri) 
{
    vec4 s, t;
    s = agi_twoSumComplex(a_ri.xz, b_ri.xz);
    t = agi_twoSumComplex(a_ri.yw, b_ri.yw);
    s.yw += t.xz;
    s = agi_quickTwoSumComplex(s.xz, s.yw);
    s.yw += t.yw;
    s = agi_quickTwoSumComplex(s.xz, s.yw);
    return s;
}

vec4 agi_dfCompSquareComponentwise(vec4 a) 
{
    vec4 p;
    vec4 s;
    p = agi_twoSquareComplex(a.xz);
    p.yw += 2.0*a.xz*a.yw;
    p.yw += a.yw*a.yw;
    s = agi_quickTwoSumComplex(p.xz, p.yw);
    return s;
}*/
