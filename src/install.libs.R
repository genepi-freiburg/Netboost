#' Build of the mcupgma executables inplace and copy them to temp build.
#' Basically, the build using make is not perfect here, but also has
#' advantages: no own Makefile for the Rcpp parts is required (default is
#' safer). As mcupgma has a fixed Makefile without dependencies, but given
#' options, basically a rewrite would be required to inject the R compiler
#' settings.
#' This leaves "make" as only hard dependency (and running *nix).
#' C/C++ compiler settings from OS/defaults.
#'
#' Oct 2017, jo

if (WINDOWS) stop("mcupgma build on Windows is currently not supported")

#print(paste("WINDOWS:", WINDOWS))
#print(paste("R_PACKAGE_NAME", R_PACKAGE_NAME))
print(paste("R_PACKAGE_SOURCE", R_PACKAGE_SOURCE))
print(paste("R_PACKAGE_DIR", R_PACKAGE_DIR))
print(paste("R_ARCH", R_ARCH))
print(paste("SHLIB_EXT", SHLIB_EXT))

## mcupgma lying in src folder
src_mcupgma  <- file.path(R_PACKAGE_SOURCE, 'src', 'mcupgma')

## Destination path of R package
dest_lib <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))

## MCUPGMA folder (binaries, scripts and additional files)
dest_mcupgma <- file.path(R_PACKAGE_DIR, 'mcupgma')

print(paste("INSTALL PATH: ",
            file.path(src_mcupgma, "scripts", "install_path.mk")))

# Original MCUPGMA Makefile: create installation path as separate file with
# only export of variable INSTALL_PATH.
system2("echo", paste0("export INSTALL_PATH := ",
                       dest_mcupgma,
                       " > ",
                       file.path(src_mcupgma, "scripts", "install_path.mk")))

## Build mcupgma package (package builds itself in own folder)
system2(command="make",                   ## Clean
        args=c(paste0("-C ", src_mcupgma),
               "-f Makefile", 
               "clean"))

system2(command="make",                   ## Build all
        args=c(paste0("-C ", src_mcupgma),
               "-f Makefile", 
               "all"))

## TODO Add checks.

print(paste("SRC: ", src_mcupgma))
print(paste("DEST:", dest_mcupgma))

dir.create(dest_mcupgma, recursive = TRUE, showWarnings = TRUE)

## Executables created in two directories (clustering_round and
## clustering_utils) plus some Shell/Perl scripts in directory
## scripts
for (path in c(file.path('clustering_round', 'bin'),
               file.path('clustering_util', 'bin'),
               file.path('scripts'))) {
  copy_path <- file.path(src_mcupgma, path, '*')
  files <- Sys.glob(copy_path)

  ## Very simply checks (to extend later)
  if (length(files) == 0)
    warning(paste("MCUPGMA build: no files in directory: ",
                  copy_path,
                  " (build failed)"))

  ## Any empty files?
  if (!all(sapply(files, function(x) file.info(x)$size > 0)))
    warning(paste("MCUPGMA build: empty files in directory: ",
                  copy_path,
                  " (build failed)"))
    
  print(paste("COPY FROM:", copy_path, "TO:", dest_mcupgma, "FILES:"))
  print(files)
  file.copy(files, dest_mcupgma, overwrite = TRUE)
}

files <- Sys.glob(paste0(dest_mcupgma, "/*"))
print(paste("HAVE:", files))

## Copy own (non-mcugpma) libraries (basically only netboost.so) to dest.
## Basically the code from R manual.
dir.create(dest_lib, recursive = TRUE, showWarnings = FALSE)
files <- Sys.glob(paste0("*", SHLIB_EXT))
file.copy(files, dest_lib, overwrite = TRUE)
if(file.exists("symbols.rds"))
  file.copy("symbols.rds", dest_lib, overwrite = TRUE)
