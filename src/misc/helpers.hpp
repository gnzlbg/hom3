// #ifndef HELPERS_HPP_
// #define HELPERS_HPP_
// ////////////////////////////////////////////////////////////////////////////////
// /// Includes:
// #include "types.hpp"
// #include "error.hpp"
// ////////////////////////////////////////////////////////////////////////////////
// namespace hom3 {
// ////////////////////////////////////////////////////////////////////////////////

// namespace helpers {

// template <typename T> struct identity { using type = T; }; 

// template <class _Compare, class _ForwardIterator, class _Tp>
// _ForwardIterator
// __lower_bound(_ForwardIterator __first, _ForwardIterator __last, const _Tp& __value_, _Compare __comp)
// {
//   typedef typename std::iterator_traits<_ForwardIterator>::difference_type difference_type;
//     difference_type __len = std::distance(__first, __last);
//     // std::cerr << "lower_bound len: " << __len << "\n";
//     // std::cerr << "value: " << __value_(0) << " " << __value_(1) << " " << __value_(2) << "\n";
//     // auto init_pos = __first;
//     // auto end_pos = __last;
//     while (__len != 0)
//     {
//         difference_type __l2 = __len / 2;
//         _ForwardIterator __m = __first;
//         std::advance(__m, __l2);
//         // std::cerr << "c m: " << (__m - init_pos) << " | " << (*__m)(0) << " "
//         //           << (*__m)(1) << " " << (*__m)(2) << "\n";
//         if (__comp(*__m, __value_))
//         {
//           // std::cerr << "cmp true: value at " << (__m - init_pos) << " = "
//           //           << (*__m)(0) << " " << (*__m)(1) << " " << (*__m)(2)
//           //           << " < value: " << __value_(0) << " " << __value_(1) << " " << __value_(2) << "\n";
//             __first = ++__m;
//             __len -= __l2 + 1;
//         }
//         else {
//           // std::cerr << "cmp false: value at " << (__m - init_pos) << " = "
//           //           << (*__m)(0) << " " << (*__m)(1) << " " << (*__m)(2)
//           //           << " >= value: " << __value_(0) << " " << __value_(1) << " " << __value_(2) << "\n";
          
//             __len = __l2;

//         }
//     }
//     return __first;
// }

// }

// ////////////////////////////////////////////////////////////////////////////////
// } // hom3 namespace
// ////////////////////////////////////////////////////////////////////////////////
// #endif
