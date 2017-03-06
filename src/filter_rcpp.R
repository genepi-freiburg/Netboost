#' Setup: probing hardware capabilities (SSE 4.2 and AVX) and then compile the
#' main program with according feature sets.
#' Therefore two programs are compiled:
#' - CPU probe (rcpp_test_cpu.cpp): test CPU for specs
#' - Boosting (rcpp_boosting.cpp) : "main" program
#'
#' @param path  String path to C++ files
#' @param optim Logical: activate aggressive optimization ("-O3 -funroll-loops")
#' @return Mode flags with AVX/SSE4.x usable (and already compiled)
setup_rcpp_filter <- function(path=".", optim=TRUE) {
    # Save current C++ compiler flags.
    cxx.flags <- Sys.getenv("PKG_CXXFLAGS")

    ## Test current CPU if SSE 4 and AVX units are available: compile a tiny
    ## testing program:
    ## Test is heavily compiler dependent (currently only GCC/Clang/MSVC)
    ## => check C source if there are problems.
    Sys.setenv("PKG_CXXFLAGS" = "-mno-avx -mno-sse4.2 -mno-sse4.1")
    sourceCpp(file.path(path, "rcpp_test_cpu.cpp"))

    ## Compiler flags to build main program: always C++11
    cxxflags <- c("-std=c++11")

    ## With O3 GCC starts to include auto-vectorisation on it's own (which does
    ## not really works in the given program), loop unrolling works fine, but
    ## benefit is not too big (data loading outweights all).
    if (optim) {
        cxxflags <- c(cxxflags, "-O3", "-funroll-loops")
    }

    use_avx <- FALSE
    
    ## If no SSE 4.1&4.2 is around, AVX is automatically skipped.
    if (!has_sse4_2()) {
        ## Preprocessor flag NO_AVX is required to exclude explicit AVX
        ## code during compilation (as probing for machine flags with the
        ## preprocessor does not seem to be possible).
        cxxflags <- c(cxxflags, "-mno-sse4.2", "-mno-avx", "-mno-fma", "-DNO_AVX")
        print("No SSE 4.2 detected. Skipping both SSE and AVX (slower)")
    } else {
        cxxflags <- c(cxxflags, "-msse4.2")

        if (!has_avx()) {
            cxxflags <- c(cxxflags, "-mno-avx", "-mno-fma", "-DNO_AVX")
            print("SSE 4.2 present, but no AVX detected. Fastmode disabled.")
        } else {
            use_avx <- TRUE
            cxxflags <- c(cxxflags, "-mavx")
            print("Very good, AVX is usable.")
        }
    }

##    print(paste(cxxflags))

    ## Setup compiler for main program: C++11 and detected hardware
    ## Problem: SSE/AVX must be enabled on compile => required to
    ## separate the C++ files.
    Sys.setenv("PKG_CXXFLAGS" = paste(cxxflags, sep="", collapse=" "))
    sourceCpp(file.path(path, "rcpp_boosting.cpp"))

    ## Restore compiler settings.
    Sys.setenv("PKG_CXXFLAGS" = cxx.flags)

    return (use_avx);
}

#' Boosting via C++ function. Parallelisation by R-package parallel with forking
#' (overhead of this method does not fall into account as single steps are ~10s).
#'
#' @param datan     Dataset
#' @param stepno    Integer amount of boosting steps
#' @param until     Stop at index/column (if 0: iterate through all columns)
#' @param progress  Integer. If > 0, print progress after X steps (mind: parallel!)
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential)
#' @param mode      Integer. Mode flags received by setup_rcpp_filter()
#' @return List (for each column vector with items received by boosting)
all_rcpp_para <- function(datan, stepno = 20, until = 0,
                          progress=100,
                          cores=2, mode=0) {
    if (until == 0)
        until = ncol(datan);

    ## Initialize data structures for optimized boosting (once)
    rcpp_filter_base(datan, stepno, mode=mode);

    ## Parallelization "conventional" via mclapply. Not really accountable overhead,
    ## as single calls take ~10 seconds.
    if (cores > 1) {
        print(paste("Parallel version:", cores, "cores"))
        
        ret <- mclapply(seq(1, until),
                        function(x) {
                            if ((progress > 0) && ((x-1) %% progress == 0)) {
                                print(sprintf("idx: %d (%.1f%%) - %s", x, x * 100 / until, date()))
                            }

                            rcpp_filter_step(x)
                        },
                        mc.cores=cores)
    } else {   ## Sequential function for debugging.
        print(paste("Sequential version"))

        ret <- lapply(seq(1, until),
                      function(x) {
                          if (progress && ((x-1) %% progress == 0)) {
                              print(timestamp(quiet=TRUE))
                              print(sprintf("idx: %d (%.1f%%)", x, x * 100 / until))
                          }
                          
                          rcpp_filter_step(x)
                      })
    }    

    ## Important!: stop (free memory, else suitable memory is still blocked)
    rcpp_filter_end();

    return(ret)
}

#' Transform result from C++ function to list of data frames format used by
#' Pascal.
#' (C++ returns list with each element containing the vector of pairing elements
#' with the element of the given list index).
#' 
#' @param list.data List. Result from all_rcpp_para()
#' @return List of data.frames (each with two columns with selected pairs)
transform_result <- function(list.data=list()) {
    result <- lapply(1:length(list.data), function(x) {
        return(as.data.frame(cbind(as.integer(list.data[[x]]),
                                   as.integer(rep(x,length(list.data[[x]]))))))
    })

    return(result)
}

library(Rcpp)
library(parallel)

load("../data.Rdata")

## Setup and compile C++ function (according to current CPU features)
mode <- setup_rcpp_filter(path=".", optim=TRUE)

## Boosting: must be run after setup_rcpp_filter()
result <- all_rcpp_para(data, stepno=20, cores=24, mode=mode, progress=1000)

## Transform results to the format of the original R program
result_format_pascal <- transform_result(result)

save(file="./boost_cpp_results.Rdata", result)
