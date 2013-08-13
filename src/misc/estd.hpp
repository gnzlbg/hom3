// // -----------------------------------------------------------------------------
// //! Enhanced C++11:
// //! @brief: Functionality omitted in the c++11 std results in weird idioms.
// // -----------------------------------------------------------------------------
// // Includes:
// #include <memory>
// #include <list>
// #include <forward_list>
// #include <vector>
// #include <deque>
// #include <map>
// #include <unordered_map>
// #include <set>
// #include <unordered_set>
// #include <algorithm>
// // -----------------------------------------------------------------------------

// //! Puts the enhancements in the estd namespace
// namespace estd {

// // -----------------------------------------------------------------------------
// //! Idioms:
// // -----------------------------------------------------------------------------

// //! erase_remove_if:
// //! complexity: O(n)
// //! Supports:
// //! sequence containers: vector, list, deque, forward_list
// //! associative containers: map, unordered_map, set

// namespace erase_remove_if_detail_ { // erase_remove_if implementation details

// //! container tags:
// struct contiguous_tag {};
// struct list_tag {};
// struct associative_tag {};

// //! container traits:
// template<class T> struct traits {
// typedef typename T::erase_remove_if::unknown_type::error type; };
// template<class T,class U> struct traits<std::vector<T,U>>
// {typedef contiguous_tag type;};
// template<class T,class U> struct traits<std::deque<T,U>>
// {typedef contiguous_tag type;};
// template<class T,class U> struct traits<std::list<T,U>>
// {typedef list_tag type;};
// template<class T,class U> struct traits<std::forward_list<T,U>>
// {typedef list_tag type;};
// template<class T,class U,class V,class W> struct traits<std::map<T,U,V,W>>
// {typedef associative_tag type;};
// template<class T,class U,class V,class W> struct traits<std::unordered_map<T,U,V,W>>
// {typedef associative_tag type;};
// template <class T,class U,class V,class W> struct traits<std::multimap<T,U,V,W>>
// {typedef associative_tag type;};
// template <class T,class U,class V,class W> struct traits<std::unordered_multimap<T,U,V,W>>
// {typedef associative_tag type;};
// template<class T,class U,class V> struct traits<std::set<T,U,V>>
// {typedef associative_tag type;};
// template <class T,class U,class V> struct traits<std::unordered_set<T,U,V>>
// {typedef associative_tag type;};
// template <class T,class U,class V> struct traits<std::multiset<T,U,V>>
// {typedef associative_tag type;};
// template <class T,class U,class V> struct traits<std::unordered_multiset<T,U,V>>
// {typedef associative_tag type;};

// //! implementation:
// template<class T,class Predicate> void impl_(T& c,Predicate&& p,contiguous_tag){
//   c.erase(std::remove_if(std::begin(c),std::end(c),p),std::end(c));}

// template<class T,class Predicate> void impl_(T& c,Predicate&& p,list_tag){
//   c.remove_if(p);}

// template<class T,class Predicate> void impl_(T& c,Predicate&& p,associative_tag){
//   auto b = std::begin(c), e = std::end(c);
//   for(; b != e; )
//     if (p(*b)) c.erase(b++); //incr. it before deleting to avoid invalidation
//     else ++b;}

// } // erase_remove_if_detail_

// //! erase_remove_if interface:
// template<class T, class Predicate>
// void erase_remove_if(T& c, Predicate&& p){
//   typedef typename erase_remove_if_detail_::traits<T>::type container_type;
//   erase_remove_if_detail_::impl_(c,std::forward<Predicate>(p),container_type());
// }

// // -----------------------------------------------------------------------------
// //! Iterators:
// // -----------------------------------------------------------------------------

// //! Missing non-member functions:
// template<class T> auto cbegin (T& t)->decltype(t.cbegin()) {return t.cbegin() ;}
// template<class T> auto rbegin (T& t)->decltype(t.rbegin()) {return t.rbegin() ;}
// template<class T> auto crbegin(T& t)->decltype(t.crbegin()){return t.crbegin();}
// template<class T> auto cend   (T& t)->decltype(t.cend())   {return t.cend()   ;}
// template<class T> auto rend   (T& t)->decltype(t.rend())   {return t.rend()   ;}
// template<class T> auto crend  (T& t)->decltype(t.crend())  {return t.crend()  ;}

// // -----------------------------------------------------------------------------
// //! Memory:
// // -----------------------------------------------------------------------------

// //! Analog to make_shared:
// template<class T, class... Args>
// std::unique_ptr<T> make_unique(Args&&... args) {
//     return std::unique_ptr<T>( new T(std::forward<Args>(args)...) );
// }

// // -----------------------------------------------------------------------------
// //! Numeric:
// // -----------------------------------------------------------------------------


// //! Compile-time numeric functions:

// namespace ct_numeric_detail_ {
// template<class T> constexpr T factImpl_(const T n){
//   return n==0 ? 1 : n*factImpl_(n-1);}
// constexpr uint64_t powImpl_(uint64_t e,uint64_t b){
//   return e==0 ? 1 : b*powImpl_(e-1,b);}
// constexpr uint64_t logImpl_(uint64_t n,uint64_t b){
//   return n==0 ? 0 : ( n==1 ? 0 : 1 + logImpl_(n/b,b) ); }
// }

// //! Static Log_Base_(N)
// template <uint64_t Base, uint64_t N> struct log {
//   // Overflows for N > 2^64, but you can't test N if it overflows
//   static const uint64_t value = ct_numeric_detail_::logImpl_(N,Base);
// };

// template <uint64_t Base, uint64_t Exp> struct pow {
//   static_assert(Exp < ((64-1)/( log<Base, 2>::value ) ),"POW: INTEGER OVERFLOW");
//   static const uint64_t value = ct_numeric_detail_::powImpl_(Exp,Base);
// };

// //! Static Factorial: N!
// template <uint64_t N,class T = void> struct factorial {
//   static_assert( (N < 1755),"FACTORIAL: LONG DOUBLE OVERFLOW");
//   static constexpr long double value
//   = ct_numeric_detail_::factImpl_(static_cast<long double>(N));
// };

// //! Static Factorial: N!
// template <uint64_t N> struct factorial<N,uint64_t> {
//   static_assert( (N < 21) ,"FACTORIAL: INTEGER OVERFLOW");
//   static const uint64_t value = ct_numeric_detail_::factImpl_(N);
// };

// constexpr double sqrt_helper(const double x, const double g) {
//   return abs(g-x/g) < tol ? g : sqrt_helper(x,(g+x/g)/2.0);
// }

// constexpr double sqrt(const double x) { return sqrt_helper(x,1.0); }

// //! Division:
// template <class T, class U>
// constexpr auto div(T a, U b) -> decltype(a/b) { return a / b; }

// //! Multiplication:
// template <class T, class U>
// constexpr auto mult(T a, U b) -> decltype(a*b) { return a * b; }

// //! Power:
// template <class T>
// constexpr T powt(T a, uint64_t b) { return b==0 ? 1 : a*powt(a,b-1); }

// #define M_PI   3.141592653589793
// #define M_PI_2 1.570796326794897
// #define M_E    2.718281828459045

// constexpr double tol = 0.001;

// constexpr double abs(const double x)    { return x < 0.0 ? -x : x; }

// constexpr double square(const double x) { return x*x; }

// constexpr double sqrt_helper(const double x, const double g) {
//   return abs(g-x/g) < tol ? g : sqrt_helper(x,(g+x/g)/2.0);
// }

// constexpr double sqrt(const double x) { return sqrt_helper(x,1.0); }

// constexpr double cube(const double x) { return x*x*x; }

// // Based on the triple-angle formula: sin 3x = 3 sin x - 4 sin ^3 x
// constexpr
// double sin_helper(const double x) {
//   return x < tol ? x : 3*(sin_helper(x/3.0)) - 4*cube(sin_helper(x/3.0));
// }

// constexpr
// double sin(const double x) {
//   return sin_helper(x < 0 ? -x+M_PI : x);
// }

// //sinh 3x = 3 sinh x + 4 sinh ^3 x
// constexpr
// double sinh_helper(const double x) {
//   return x < tol ? x : 3*(sinh_helper(x/3.0)) + 4*cube(sinh_helper(x/3.0));
// }

// //sinh 3x = 3 sinh x + 4 sinh ^3 x
// constexpr
// double sinh(const double x) {
//   return x < 0 ? -sinh_helper(-x) : sinh_helper(x);
// }

// constexpr double cos (const double x) { return sin(M_PI_2 - x); }

// constexpr double cosh(const double x) { return sqrt(1.0 + square(sinh(x))); }

// constexpr
// double pow(double base, int exponent) {
//   return exponent <  0 ? 1.0 / pow(base,-exponent) :
//          exponent == 0 ? 1.                        :
//          exponent == 1 ? base                      :
//          base * pow(base,exponent-1);
// }

// // atan formulae from http://mathonweb.com/algebra_e-book.htm
// // x - x^3/3 + x^5/5 - x^7/7+x^9/9  etc.
// constexpr
// double atan_poly_helper(const double res,  const double num1,
//                         const double den1, const double delta) {
//   return res < tol ? res :
//          res + atan_poly_helper((num1*delta)/(den1+2.)-num1/den1,
//                                  num1*delta*delta,den1+4.,delta);
// }

// constexpr
// double atan_poly(const double x) {
//   return x + atan_poly_helper(pow(x,5)/5.-pow(x,3)/3., pow(x,7), 7., x*x);
// }

// // Define an M_PI_6? Define a root 3?
// constexpr
// double atan_identity(const double x) {
//   return x <= (2. - sqrt(3.)) ? atan_poly(x) :
//          (M_PI_2 / 3.) + atan_poly((sqrt(3.)*x-1)/(sqrt(3.)+x));
// }

// constexpr
// double atan_cmplmntry(const double x) {
//   return (x < 1) ? atan_identity(x) : M_PI_2 - atan_identity(1/x);
// }

// constexpr
// double atan(const double x) {
//   return (x >= 0) ? atan_cmplmntry(x) : -atan_cmplmntry(-x);
// }

// constexpr
// double atan2(const double y, const double x) {
//   return           x >  0 ? atan(y/x)        : 
//          y >= 0 && x <  0 ? atan(y/x) + M_PI :
//          y <  0 && x <  0 ? atan(y/x) - M_PI :
//          y >  0 && x == 0 ?  M_PI_2          :
//          y <  0 && x == 0 ? -M_PI_2          : 0;   // 0 == undefined
// }

// constexpr
// double nearest(double x) {
//   return (x-0.5) > (int)x ? (int)(x+0.5) : (int)x;
// }

// constexpr
// double fraction(double x) {
//   return (x-0.5) > (int)x ? -(((double)(int)(x+0.5))-x) : x-((double)(int)(x));
// }

// constexpr
// double exp_helper(const double r) {
//   return 1.0 + r + pow(r,2)/2.0   + pow(r,3)/6.0   +
//                    pow(r,4)/24.0  + pow(r,5)/120.0 +
//                    pow(r,6)/720.0 + pow(r,7)/5040.0;
// }

// // exp(x) = e^n . e^r (where n is an integer, and -0.5 > r < 0.5
// // exp(r) = e^r = 1 + r + r^2/2 + r^3/6 + r^4/24 + r^5/120
// constexpr
// double exp(const double x) {
//   return pow(M_E,nearest(x)) * exp_helper(fraction(x));
// }

// constexpr
// double mantissa(const double x) {
//   return x >= 10.0 ? mantissa(x *  0.1) :
//          x <  1.0  ? mantissa(x * 10.0) :
//          x;
// }

// // log(m) = log(sqrt(m)^2) = 2 x log( sqrt(m) )
// // log(x) = log(m x 10^p) = 0.86858896 ln( sqrt(m) ) + p
// constexpr
// int exponent_helper(const double x, const int e) {
//   return x >= 10.0 ? exponent_helper(x *  0.1, e+1) :
//          x <  1.0  ? exponent_helper(x * 10.0, e-1) :
//          e;
// }

// constexpr
// int exponent(const double x) { return exponent_helper(x,0); }

// constexpr
// double log_helper2(const double y) {
//   return 2.0 * (y + pow(y,3)/3.0 + pow(y,5)/5.0 +
//                     pow(y,7)/7.0 + pow(y,9)/9.0 + pow(y,11)/11.0);
// }

// // log in the range 1-sqrt(10)
// constexpr
// double log_helper(const double x) { return log_helper2((x-1.0)/(x+1.0)); }

// // n.b. log 10 is 2.3025851
// // n.b. log m = log (sqrt(m)^2) = 2 * log sqrt(m)
// constexpr
// double log(const double x) {
//   return x == 0 ? -std::numeric_limits<double>::infinity() :
//          x <  0 ?  std::numeric_limits<double>::quiet_NaN() :
//          2.0 * log_helper(sqrt(mantissa(x))) + 2.3025851 * exponent(x);
// }

// namespace traits {

// template <typename Condition, typename T = void>
// using EnableIf = typename std::enable_if<Condition::value, T>::type;

// template <typename Condition, typename T = void>
// using DisableIf = typename std::enable_if<!Condition::value, T>::type;

// }


// }
