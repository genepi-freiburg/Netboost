/**
 * Compiler/preprocessor flags must be set in Makefile (could be set dynamically
 * according to current hardware if used via rcppSource() -- this can be later
 * reimplemented at least as hardware check on installing machine).
 *
 * Flags:
 *  - NO_AVX   -- disable AVX code
 *  - NO_FMA   -- disable FMA code  (completely removed anyway atm)
 *
 * Flags with autotools:
 *  - HAVE_AVX
 *  - HAVE_FMA
 *
 * Intrinsics are only used if both AVX/FMA are available in compiler, else
 * compilation may fail.
 *
 * Code currently only tested on GCC 4.7 - 7.0 & Clang 3.7
 *
 * Portability:
 *  - "__restrict" used. Non-standard, but supported by GCC, Clang, MSVC, Intel
 *
 * Optimization:
 * - Loop unrolling for inner loops are left to the compiler
 *
 * Notes:
 * - Armadillo version is about three times slower.
 * - RcppParallel version still in debugging.
 * - FMA doesn't bring real advantage to AVX (few percent and
 *   even being slower in many cases -- tested on different Xeons).
 *
 * Jun 2016, jo
 * R package Dec 2016, jo
 */

// [[Rcpp::interfaces(r)]]

#include <Rcpp.h>

#include <cstddef>
#include <memory>
#include <vector>
#include <limits>

// Base flags from autotools.
#include "config.h"

// If PACKAGE_NAME is defined, complete autotools build is assumed.
#ifdef PACKAGE_NAME
  #ifndef HAVE_AVX
    #define NO_AVX
  #endif
  #ifndef HAVE_FMA3
    #define NO_FMA
  #endif
// "Traditional" build, using NO_AVX/NO_FMA flags.
#else
  #ifdef NO_AVX
    #undef HAVE_AVX
  #endif
  #ifdef NO_FMA
    #undef HAVE_FMA3
  #endif
#endif

// GCC includes all AVX/SSE intrinsics
// @TODO Recheck MSVC & Intel for their names.
#if defined(HAVE_AVX) && defined(HAVE_FMA3)
  #include <immintrin.h>
#endif

using namespace Rcpp;
using namespace std;

constexpr bool DEBUG = false;
constexpr int AVX_UNROLL = 4;

// Align data by x-bits (log2(byteborder)).
// Min are 32-bytes to use aligned AVX commands, 64-bytes for better cache-line
// usage.
constexpr int ALIGN_BITS = 6;

constexpr int ALIGN_BYTES = 1 << ALIGN_BITS;

class Boosting 
{
private:
  // Raw data holder for both data matrix and cov-cache matrix.
  // Must both be copies, as columns must be aligned (for use of fast AVX
  // reads and better cache-line usage).
  // Both are only used in the constructor, later only the column pointer
  // arrays are used (data_cols, Z_cols).
  std::shared_ptr<double> data_heap;
  std::shared_ptr<double> z_means_heap;
  
  std::vector<double *> data_cols;  //!< Pointers to columns in data
  std::vector<double *> Z_cols;     //!< Pointers to columns in z_means

  const unsigned int stepno;        //!< Boosting steps
  const double nu;                  //!< Constant (update)
  const size_t ncol;                //!< Input matrix columns
  const size_t nrow;                //!< Input matrix row
  
  // Lambda/function for arrays manual delete is required (for smart pointers
  // to allocated arrays).
  const std::function<void(double *)> heap_deleter = [](double *p) {
    delete[] p;
  };
  
  /**
   * @brief Sum of vector multiplication (inplace pointers, no overlap, non-AVX)
   *
   * May be unrolled by compiler depending on compilation settings.
   *
   * @param[in] a  Reference (pointer) to array a
   * @param[in] b  Reference (pointer) to array b
   * @param[in] len  Size of arrays/vectors
   * @return Sum of a*b
   */
  inline double vecmul_sum_inplace(double* &a, double* &b, const size_t len) {
    const double *end = a + len;
    
    double sum = 0;
    
    while (a < end)
      sum += *a++ * *b++;
    
    return sum;
  }
  
  /**
   * @brief Sum of vector multiplication (copy pointers, no overlap)
   *
   * May be unrolled by compiler depending on compilation settings.
   *
   * @param[in] a  Pointer to array a
   * @param[in] b  Pointer to array b
   * @param[in] len  Size of arrays/vectors
   * @return Sum of a*b
   */
  inline double vecmul_sum(double* __restrict a, double* __restrict b, const size_t len) {
    const double *end = a + len;
    double sum = 0;
    
    while (a < end)
      sum += *a++ * *b++;
    
    return sum;
  }
  
  /**
   * @brief Sum of vector multiplication using AVX (sum(a*b)). General usable.
   *
   * Fast vector sum function using the AVX vector unit. Compilation must be
   * including AVX support, as intrinsics are used (which basically transform
   * into straight assembler).
   *
   * Analogous to R: sum(a*b) with a and b being vectors.
   *
   * This version is a general purpose one, which means it has no constraints
   * on the data organisation in the input arrays.
   *
   * May be unrolled by compiler depending on compilation settings.
   *
   * *NO* security checks (like length(a) != length(b))
   *
   * @param[in]  a  Array/Vector (pointer changed, data const)
   * @param[in]  b  Array/Vector (pointer changed, data const)
   * @param[in]  len  Length of vectors
   * @return Sum of vector multiplication of a and b
   */
  inline double vecmul_sum_avx(const double* __restrict a, const double* __restrict b,
                               const size_t len) {
#ifdef HAVE_AVX
    double avx_sum = 0;
    
    // AVX only comes into account if more than 4 doubles are around.
    if ((len / AVX_UNROLL) > 0 ) {
      // Clear result register (YMM)
      __m256d sum = _mm256_setzero_pd();
      __m256d x, y;
      
      for (auto i = (len / AVX_UNROLL); i > 0 ; i--, a += AVX_UNROLL, b += AVX_UNROLL) {
        // Fetch 32 bytes from both vectors. Read unaligned (loadu)
        x = _mm256_loadu_pd(a);          // x{1..4}
        y = _mm256_loadu_pd(b);          // y{1..4}
        
        x = _mm256_mul_pd(x, y);         // x = x*y
        sum = _mm256_add_pd(x, sum);     // sum = x*y + sum
        
        // Honestly, unrolling has no real benefit here (even tested with much larger
        // unrolling values).
        if(AVX_UNROLL > 4) {
          x = _mm256_loadu_pd(a+4);
          y = _mm256_loadu_pd(b+4);
          
          x = _mm256_mul_pd(x, y);
          sum = _mm256_add_pd(x, sum);
        }
      }
      
      // Sum all four doubles in register 'sum' into a single value.
      __m256d hsum = _mm256_add_pd(sum, _mm256_permute2f128_pd(sum, sum, 0x1));
      
      // Copy lowest register double into memory.
      _mm_store_sd(&avx_sum, _mm_hadd_pd(_mm256_castpd256_pd128(hsum), _mm256_castpd256_pd128(hsum)));
    }
    
    // Processing the leftover doubles.
    for (auto i = len % AVX_UNROLL; i > 0; --i)
      avx_sum += *a++ * *b++;
    
    return avx_sum;
#else
    Rf_error("Calling of AVX function without compiler support defined (AVX usage disabled or not available).");
#endif
  }

  /**
   * @brief Faster vector multiplication using AVX (sum(a*b)) for aligned data.
   *
   * Faster than normal AVX version (10-50%, depending what is cached), but with
   * strict constraints:
   * - Data (both pointers!) must be aligned to (at least) 32-byte boundary
   * - Data after a+len and b+len must be 0 until next 32-byte boundary
   *
   * If contraints are not fulfilled, program crashes (unaligned) or returns
   * wrong results (non zero filled).
   *
   * Both contraints allow to read aligned (unaligned often two cache lines must
   * be filled for a single 32-byte read) and to skip adding of the rest (by
   * simply adding over the end).
   *
   * May be unrolled by compiler depending on compilation settings.
   *
   * @param[in] ar  Vector/array aligned to 32-byte boundary, zero filled end
   * @param[in] br  Vector/array aligned to 32-byte boundary, zero filled end
   * @param[in] len  Amount of elements in both vectors
   * @return Sum of a*b
   */
  //  double vecmul_sum_avx_aligned(double *a, double *b, const size_t len) __attribute__ ((__option__("avx")));
  double vecmul_sum_avx_aligned(double * __restrict ar, double * __restrict br,
                                const size_t len) {
#ifdef HAVE_AVX
    __m256d sum = _mm256_setzero_pd();
    
    // Using GCC the compiler can be informed data is aligned (may help on auto
    // unrolling of loop).
#ifdef __GNUG__
    auto a = (double *) __builtin_assume_aligned(ar, ALIGN_BYTES);
    auto b = (double *) __builtin_assume_aligned(br, ALIGN_BYTES);
#else
    auto a = ar;
    auto b = ar;
#endif    
    
    //  double *end = a + (len / 4) + ((len % 4) != 0);
    
    if (true) {
      size_t add = 0;
      
      // Note: unrolled version is actually slower. Simply because loop control
      // is executed in parallel, as well as loading.
      while (add < len) {
        __m256d x = _mm256_load_pd(a + add);    // Aligned read 32-bytes
        __m256d y = _mm256_load_pd(b + add);
        __m256d z = _mm256_mul_pd(x, y);  // z = x * y
        sum = _mm256_add_pd(z, sum);      // sum = sum + z
        add += 4;
        
        /* Translates to (a + add is implicit), GCC 4.8, -O2:
         vmovapd (%rdi,%rax,8), %ymm1
         vmulpd  (%rsi,%rax,8), %ymm1, %ymm1
         addq    $4, %rax
         cmpq    %rax, %rdx
         vaddpd  %ymm0, %ymm1, %ymm0
         ja      .L5
         */
      }
    } else {
      double *end = a + len;
      while (a < end) {
        __m256d x = _mm256_load_pd(a);    // Aligned read 32-bytes
        __m256d y = _mm256_load_pd(b);
        __m256d z = _mm256_mul_pd(x, y);  // z = x * y
        sum = _mm256_add_pd(z, sum);      // sum = sum + z
        a += 4;
        b += 4;
        
        /* Translates to:
         .L5:
         vmovapd (%rdi), %ymm1
         addq    $32, %rdi
         addq    $32, %rsi
         vmulpd  -32(%rsi), %ymm1, %ymm1
         cmpq    %rdi, %rax
         vaddpd  %ymm0, %ymm1, %ymm0
         ja      .L5
         */
      }
    }
    
    // Transfer (AVX) YMM register to memory -- depending on compiler, it may be
    // returned in (SSE) XMM register (e.g. GCC does return single doubles in a
    // SSE registers).
    double sum_memory;
    
    __m256d hsum = _mm256_add_pd(sum, _mm256_permute2f128_pd(sum, sum, 0x1));
    _mm_store_sd(&sum_memory, _mm_hadd_pd(_mm256_castpd256_pd128(hsum), _mm256_castpd256_pd128(hsum)));
    
    return sum_memory;
#else
    Rf_error("Calling of AVX function without compiler support or AVX disabled.");
#endif
  }
  
  /**
   * @brief Align pointer to margin with zeroing gaps.
   *
   * Small helper to align a pointer to a given margin (margin must be power of 2).
   * Beside that, all space outside of the margins is filled with 0 (any type).
   * In a perfect world this would be implemented as allocator for the STL
   * containers, but this is a lazy world...
   *
   * @tparam T         Type of pointer to be aligned
   * @tparam ALIGNBIT  Border to be aligned (2^ALIGNBITS) or log2(ALIGNBYTES)
   * @param[in] src    Pointer to be aligned.
   * @return Aligned   pointer (>= src. If > src, area between src and ret is zeroed)
   */
  template <typename T, int ALIGNBIT>
  T* align(T* src) {
    // As logical operations on pointers are supported by even less compilers.
    size_t pointer_int = (size_t) src;
    
    // Amount of misaligned bytes.
    size_t misaligned_bytes = pointer_int & ((1<<ALIGNBIT) - 1);
    
    // printf("DO: %lli MISALIGN: %lli (mask: %d). ADDING %d ZEROS\n", pointer_int, misaligned_bytes,
    // 	 ((1<<ALIGNBIT) - 1), (((1<<ALIGNBIT) - misaligned_bytes) / sizeof(T)));
    
    // If pointer was unaligned, zero space until aligned.
    if (misaligned_bytes > 0) {
      for (T i = 0; i < (((1<<ALIGNBIT) - misaligned_bytes) / sizeof(T)); i++) {
        *src++ = 0;
      }
    }
    
    return src;
  }
  
public:
  /**
   * Constructor: set up used data structures for speeding up the calculation
   * of the covariance.
   * Basically this requires two full copies of rdata (would be only one if
   * there would be a possibility to get the data pointers from the Rcpp data
   * structures).
   *
   * As the constructor is only called once, this is not required to be critically
   * optimized.
   */
  Boosting(const NumericMatrix &rdata, unsigned int stepno_ = 20, double nu_ = 0.1) :
        stepno(stepno_),
        nu(nu_),
        ncol(rdata.ncol()),
        nrow(rdata.nrow()) {
    // In the copies space must be reserved for pointer alignments (of the column
    // vectors).
    const size_t size_aligned_copy = ncol * (nrow + (ALIGN_BYTES/sizeof(double)) + 1);
    
    if (DEBUG)
      Rcout << "ALIGNED DOUBLES: " << size_aligned_copy << endl;
    
    data_heap.reset(new double[size_aligned_copy], heap_deleter);
    z_means_heap.reset(new double[size_aligned_copy], heap_deleter);
    
    // Column pointers: setup for columns. Faster than ncol times push_back
    // (this is assumed, not measured).
    data_cols.resize(ncol);
    Z_cols.resize(ncol);
    
    // Aligned starting pointer to data copies.
    auto seq_data = align<double, ALIGN_BITS>(data_heap.get());
    auto seq_z_means = align<double, ALIGN_BITS>(z_means_heap.get());
    
    //    Rprintf( "STARTING POINTER: %p - %p\n", seq_data, seq_z_means);
    
    // Copy data with aligned columns.
    // Build Z = data - colMeans(data). Store aligned.
    for (size_t col=0; col < ncol; col++) {
      data_cols[col] = seq_data;
      Z_cols[col] = seq_z_means;
      
      // Copy rdata[,col] to current running pointer.
      std::copy(rdata.column(col).begin(), rdata.column(col).end(), seq_data);
      
      double col_sum = 0;           // Column sum for mean
      
      for (size_t row=0; row < nrow; row++)
        col_sum += *seq_data++;
      
      // Mean
      double col_mean = col_sum / nrow;
      
      auto tmp_seq_data = data_cols[col];
      
      // Z[,col] = data[,col] - mean(data[,col])
      for (size_t row=0; row < nrow; row++)
        *seq_z_means++ = *tmp_seq_data++ - col_mean;
      
      if (DEBUG)
        Rprintf( "  DATA POINTER REAL LENGTH: %li / %li dbl (matrix: %d -- max %d)\n",
                 (size_t)seq_data - (size_t)data_heap.get(),
                 ((size_t)seq_data - (size_t)data_heap.get()) / sizeof(double),
                 (col+1) * nrow,
                 (col+1) * (nrow + 256/sizeof(double))
        );
      
      // Realign pointers before copy next column.
      seq_data = align<double, ALIGN_BITS>(seq_data);
      seq_z_means = align<double, ALIGN_BITS>(seq_z_means);
    }
    
    if (DEBUG) {
      Rprintf( "END POINTER: %p - %p\n", seq_data, seq_z_means);
      Rprintf( "  DATA END POINTER SHOULD BE max: %p (REAL LENGTH: %li / %li dbl)\n",
               data_heap.get() + size_aligned_copy,
               (size_t)seq_data - (size_t)data_heap.get(),
               ((size_t)seq_data - (size_t)data_heap.get()) / sizeof(double)
      );
      Rprintf( "  zMEANS END POINTER SHOULD BE: %p (DIFF: %li bytes / %li dbl)\n",
               z_means_heap.get() + size_aligned_copy,
               (size_t)seq_z_means - (size_t)z_means_heap.get(),
               ((size_t)seq_z_means - (size_t)z_means_heap.get()) / sizeof(double)
      );
    }
    
    // for (auto col = 0; col < ncol; col++) {
    //   Rprintf("COL %d. DATA: %p. MEANS: %p\n",
    // 	      col, data_cols[col], Z_cols[col]);
    //   Rcout << "DATA: " << endl;
    //   for (auto row = 0; row < nrow; row++)
    // 	Rcout << data_cols[col][row] << endl;
    //   Rcout << "MEANS: " << endl;
    //   for (auto row = 0; row < nrow; row++)
    // 	Rcout << Z_cols[col][row] << endl;
    // }
    
    // // Convert data to flat vector.
    // // Note: this is a copy, therefore original data can be accessed with 'data'
    // data = Rcpp::as<std::vector<double>>(rdata);  // NumericMatrix => std::vector
    
    // // This is not working...
    // Z_means = data;
    
    // double *p = (double *) rdata.begin();
    // Rprintf( "CHECK: %p, VALUE: %f -- MAT: %f\n", p, *p, rdata(0,0));
    // Rprintf( "  -- CHECK 10: %p, VALUE: %f -- MAT: %f\n", p + 10, *(p+10), rdata(10,0));
    
    // // Pointer to data[0,0]
    // auto data_pointer = data.data();
    // auto seq = Z_means.data();
    
    // // // Build Z = data - colMeans(data)
    // for (auto col=0; col < ncol; col++) {
    //   // En-passant build column pointer on data (so no own loop is required)
    //   data_cols[col] = &data_pointer[col * nrow];
    
    //   // Save pointer to colum start
    //   Z_cols[col] = seq;
    
    //   double col_sum = 0;           // Column sum for mean
    
    //   for (auto row=0; row < nrow; row++)
    // 	col_sum += *seq++;
    
    //   // Mean
    //   double col_mean = col_sum / nrow;
    
    //   auto s = seq;
    
    //   // Z[,col] - mean(Z[,col])
    //   while (s > Z_cols[col])
    //   	*(--s) -= col_mean;
    // }
  }
  
  ~Boosting() {
    // Should not be required, as no memory allocated with new/malloc
    // But as those pointers are pointing inside of another structure, ensure
    // they are not freed.
    for (size_t col = 0; col < ncol; col++ ) {
      data_cols[col] = nullptr;
      Z_cols[col] = nullptr;
    }
    
    if (DEBUG)
      Rcout << "Destruct Boosting object." << endl;
  }
  
  /**
   * @brief Boosting.
   *
   * Boosting implementation.
   * Templated to have a single main function, which can be compiled to usenetboost:::
   * different backends for the sum of a vector multiplication, which takes
   * most of the time here.
   *
   * @tparam AVX  Use AVX unit for vector-mul-sum
   * @tparam FMA  Use FMA unit for vector-mul-sum (currently unsupported)
   * @tparam ALIGNED Use functions for aligned loading
   * @param[in] col_y  Column (in original data-matrix!) to processes
   * @return Vector with chosen element indizes (max. size stepno)
   */
  template <bool AVX, bool FMA, bool ALIGNED>
  IntegerVector filter_column_gen(size_t col_y) {
    // Alloc unibeta vector on the heap. Smart pointer for automatic free.
    // Advantage over std::vector: memory is not zeroed on alloc (we fill it
    // anyway).
    std::shared_ptr<double> unibeta_ptr(new double[ncol - 1], heap_deleter);

    // Now, as we have a "cleaner", we continue C style.
    double *unibeta = unibeta_ptr.get();
    
    size_t offset = 0;

    // R-Code:
    //   X = data[,-col_y]
    //   y = data[,col_y]
    //   unibeta = colSums(X*y)
    //
    // NOTE: matrix X is only "created" by skipping column from data.
    auto y = data_cols[col_y];
    
    for (size_t col = 0; col < ncol; col++) {
      // Skip column if it is y.
      if (col == col_y) {
        offset = 1;
      }
      else {
        double col_sum;
        
        // Choose specific routine (this cascade is completely optimized away
        // as parameters are from the template).
        // (Could also be done via template specialisation)
        if (AVX) {
          if (ALIGNED)
            col_sum = vecmul_sum_avx_aligned(data_cols[col], y, nrow);
          else
            col_sum = vecmul_sum_avx(data_cols[col], y, nrow);
        } else {   // Fallback: X86 version
          col_sum = vecmul_sum(data_cols[col], y, nrow);
        }
        
        unibeta[col - offset] = col_sum;
      }
    }
    
    // Alias reference: reuse already allocated vector unibeta as actual_nom.
    auto &actual_nom = unibeta;

    double actual_update = 0;
    size_t actual_sel = 0;
    
    // Index of selected betas (for faster search and delete)
    // (As beta ist a very sparse vector, allocating and clearing it and
    // searching for non-0 values later would be overkill compared to this
    // approach in which only those values which are non-0 are stored).
    // KEY: Index in beta vector (length(ncol - 1)!)
    // VALUE: Actual beta-value. May be updated several times
    std::map<decltype(col_y),double> beta_selected;

    for (unsigned int bstep=0; bstep < stepno; bstep++) {
      if (bstep > 0) {
        auto offset = (actual_sel >= col_y) ? 1 : 0;
        
        auto yy = Z_cols[actual_sel + offset];
        
        // R: actual.nom <- actual.nom - actual.update*cov(yy,x)
        for (size_t col = 0, offset = 0; col < (ncol - 0); col++) {
          if (col == col_y) {
            offset = 1;
            continue;
          }
          
          double cov_col;
          
          if (AVX) {
            if (ALIGNED)
              cov_col = vecmul_sum_avx_aligned(Z_cols[col], yy, nrow);
            else
              cov_col = vecmul_sum_avx(Z_cols[col], yy, nrow);
          } else {
            cov_col = vecmul_sum(Z_cols[col], yy, nrow);
          }
          
          // Finish cov() calculation with 1/(n-1) * sum...
          cov_col /= (nrow - 1);
          
          actual_nom[col - offset] -= (actual_update * cov_col);
        }
      }
      
      /**
       * R:
       *      actual.score <- actual.nom^2
       *      actual.sel <- which.max(actual.score)
       */

      // Smallest number possible on type double
      double max_score = std::numeric_limits<double>::min();
      
      //      auto seq = actual_nom.data();
      auto seq = actual_nom;
      
      // Note: search is in actual_nom, which has length ncol - 1,
      // so indizes are not equal to those in data for index > col_y
      for (size_t col = 0; col < ncol - 1; col++, seq++) {
        double seq_pow2 = *seq * *seq;
        
        // ^2 is actually faster than using abs()
        if (seq_pow2 > max_score) {
          max_score = seq_pow2;
          actual_sel = col;     // Save index of maximum value.
        }
      }
      
      /**
       * R:
       *      actual.update <- nu*actual.nom[actual.sel]
       *      actual.beta[actual.sel] <- actual.beta[actual.sel] + actual.update
       */
      actual_update = nu * actual_nom[actual_sel];
      
      // Store those betas which are != 0. As stored in map, if several positions
      // are selected multiple times, their value is updated.
      // [Note: accessing with []-operator on a non-existing element creates it
      // with default value (0 for double)]
      beta_selected[actual_sel] += actual_update;
      
      if (DEBUG)
        Rprintf( "C STEP: %d (ACTUAL_SEL: %d, UPDATE: %f)\n   ",
                 bstep + 1, actual_sel + 1, actual_update);
    }
    
    // Rcout << "SELECTED INDIZES IN BETA: " << endl;
    // for (auto i: beta_selected)
    //   Rprintf("INDEX: %d - VALUE: %f\n", i, beta[i]);
    
    // Result set: selected indizes (as indizes in original data matrix) and
    // also directly usable in R.
    IntegerVector selected_pairs;
    
    for (auto &bs: beta_selected) {
      auto elem = bs.first;
      
      // Index to index in original data-array
      // Also directly change to 1-based indizes.
      elem += (elem >= col_y) ? 2 : 1;
      
      selected_pairs.push_back(elem);
    }
    
    return selected_pairs;
  }
  
  size_t get_ncol() {
    return this->ncol;
  }
  
  size_t get_nrow() {
    return this->nrow;
  }
};

// Global object to hold optimized setup (which is done once).
Boosting *boost = nullptr;

// Should be replaced with function reference to templated method.
int mode = 0;

//' @title Initialise boosting with chosen accelerator hardware (x86, AVX, FMA)
//' 
//' @param data Matrix
//' @param stepno Amount of steps
//' @param mode_ Accelerator mode (0: x86, 1: FMA, 2: AVX)
//' @return none
// [[Rcpp::export(name=cpp_filter_base)]]
void filter_base(const NumericMatrix &data, unsigned int stepno = 20,
                 int mode_ = 2) {
  // If AVX is switched off, mode is always 0.
#ifndef HAVE_AVX
  mode = 0;
#else
  mode = mode_;
#endif
  
  if (boost == nullptr)
    boost = new Boosting(data, stepno);
}

//' @title Boosting cleanup (required to free memory)
//' 
//' @return none
// [[Rcpp::export(name=cpp_filter_end)]]
void filter_end() {
  if (boost != nullptr)
    delete boost;
  
  boost = nullptr;
}

//' @title Single boosting step
//' 
//' @details Must be initialised before using @see{filter_base}
//' 
//' @param col_y Row in data matrix
//' @return integer vector
// [[Rcpp::export(name=cpp_filter_step)]]
IntegerVector rcpp_filter_step(size_t col_y) {
  if (boost == nullptr)
    Rf_error("boost not initialized. Call rcpp_filter_base first");
  
  // AVX, no-FMA, aligned
  if (mode == 2)
    return boost->filter_column_gen<true,false,true>(col_y-1);
  // no-AVX, no-FMA, aligned
  else
    return boost->filter_column_gen<false,false,true>(col_y-1);
}
