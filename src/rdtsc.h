/**
 * \file rdtsc.h
 * \brief Portable rdtsc/rdtscp wrapper for both C and C++
 *
 * \details Small header-only library to access RDTSCP or RDTSC cycle counter.
 *
 * _Note:_
 * this only accesses the counter for the current core, so a process should
 * pin itself to a single core during measurement (so the operating system
 * does not switch to another core during runtime). On Linux e.g. taskset can
 * be used.
 *
 * _Note:_ this is single-core only. For multi-core measurements, use PCM (e.g.
 * with libs from Intel or Agner). But those require drivers or root access,
 * which the basic cycle counters do not.
 *
 * C requirement:       C99 (portable types. Using GCC: use --std=gnu99)
 * C++ requirement:     C++11 (portable types)
 *
 * Supported compilers: GNU, Clang
 * Supported untested:  ICC, MSVC
 *
 * By default RDTSCP is used, this can be switched back to RDTSC by defining
 * symbol USE_RDTSCP.
 *
 * No live probe is done whether the CPU supports RDTSCP (by cpuid) to keep
 * overhead small. All recent Intel and AMD CPUs support both commands.
 *
 * \author    Jochen Knaus
 * \version   0.2
 * \date      2016-10-11
 * \warning   May not compile on unsupported compilers
 * \copyright GPL 2
 */

#ifndef TICTOC__RDTSC_H
#define TICTOC__RDTSC_H

// Define to use RDTSCP instead of RDTSC instruction.
// #define USE_RDTSCP  1

// Only makes sense on i386/x86 platforms.
#if !( (defined(i386) || defined(__i386) || defined(__i386__)) || \
       defined(__x86_64) || defined(__x86_64__) || \
       defined(__amd64) || defined(__amd64__))
  #warning "Only i386/x86 platforms supported."
#endif

// Alternatively uint64_t can be defined manually.
#if defined(__cplusplus)
  // Requires C++11 (before C++11 is depending on compiler)
  #include <cstdint>
#else
  // Requires C99
  #include <stdint.h>
#endif

/*! @brief Get current state from the cyclecounters (via rdtscp/rdtsc)
 *
 * @return Current state of the cycle counter
 */
inline uint64_t cycles() 
{
  // GCC compiler collection.
  // Also works with Clang.
  // Also ICC/Linux is compatible with this (@todo untested)
  #if defined(__GNUC__) || defined(__GNUG__) || defined(__clang__) || \
      (defined(__linux__) && defined(__INTELCOMPILER))
    uint64_t cycles;

    // Probe 32/64 bit on compilation (as constant it's removed during compil.)
    if (sizeof(void*) == 8) {
      // Using RDTSC (first clear pipeline, then probe)
      #ifdef USE_RDTSC
        // rdtsc must not be optimized away, as rdtsc is side effect.
        // Flushing pipeline according to recommendation in Intels performance
        // guide.
        asm volatile(
          "cpuid\n\t"           // flush instruction pipeline
          "rdtsc\n\t"           // edx (high) - eax (low) 
          "shl  $32, %%rdx\n\t" // edx = edx << 32
          "or   %%rdx, %0"      // eax = eax | ( edx << 32 )
          : "=a" (cycles)       // result in eax
          :
          : "rbx", "rcx", "rdx" // cpuid trashes a-d (a used as result)
        );
      #else
        // RDTSCP (instruction waits itself, but cpuid afterwards to prevent
        // Out-of-order execution of later instructions)
        asm volatile (
          "rdtscp\n\t"          // edx (high) - eax (low)
          "shl  $32, %%rdx\n\t" // edx = edx << 32
          "mov  %%rax, %0\n\t"
          "or   %%rdx, %0\n\t"  // eax = eax | ( edx << 32 )
          "cpuid\n\t"           // prevent out of order (and trashes a-d)
          : "=r" (cycles)       // result in normal reg
          :
          : "rax", "rbx", "rcx", "rdx" // cpuid trashes a-d (a used as result)
        );
      #endif
      // 32-bit versions.
    } else {
      uint32_t high, low;

      #ifdef USE_RDTSC
        // rdtsc must not be optimized away, as rdtsc is side effect.
        // Flushing pipeline according to recommendation in Intels performance
        // guide.
        asm volatile (
          "cpuid\n\t"           // flush instruction pipeline
          "rdtsc\n\t"           // edx (high) - eax (low) 
          : "=d" (high), "=a" (low)       // result in eax assigned to var
          :
          : "ebx", "ecx"        // cpuid trashes a-d (a/d used as result)
        );
      #else
        asm volatile (
          "rdtscp\n\t"          // edx (high) - eax (low)
          "mov  %%edx, %0\n\t"
          "mov  %%eax, %1\n\t"
          "cpuid\n\t"           // prevent out of order (and trashes a-d)
          : "=r" (high), "=r" (low)       // split result in two 32-bit registers
          :
          : "eax", "ebx", "ecx", "edx"    // trashes a-d
        );
      #endif

      cycles = ((uint64_t) high << 32) | (uint64_t)low;
    }

    return cycles;
  #else
    // Both VC and Intel codes untested(!)
    // With MSVC counter can be accessed via intrinsic (since V-Stuido 2005)
    #ifdef _MSC_VER
      // Hopefully compiler does not shuffle the intrinsics
      #include <intrin.h>

      #ifdef USE_RDTSC
        __cpuid();           // Serialization still required.
        return __rdtsc();    // uint64_t
      #else
        uint64_t cycles = __rdtscp();
        __cpuid();
        return cycles;
      #endif
    #else
      // On Linux, ICC should compile GCC code above.
      // @todo Check if header files are required and this compiles.
      // @todo Check if cpuid is used (or implicitely added), no info found yet
      #ifdef __INTELCOMPILER
        #ifdef USE_RDTSC
          return _rdtsc();
        #else
          return _rdtscp();
	#endif
      #else
        #error "Dead end in compiler selection building rdtsc()"
      #endif
    #endif
  #endif
}

#endif
