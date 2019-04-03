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

# Internal helper: create visible error messages.
nbErrorMsg <- function(msg="") {
  row <- paste0(rep("*", 70), collapse = "")

  return(paste0(c("", row, msg, row, ""), collapse = "\n", sep = ""))
}

# Internal helper: check if a given program is installed
# (simple check: calling creates an error)
nbCheckPrg <- function(prg, txt) {
  # Check if program is present (all feature --version or ignore it)
  program <- tryCatch(system2(prg, "--version",
                              stdout = TRUE, stderr = TRUE),
                      error=function(err) NULL)
  
  if (is.null(program))
    return(paste(txt, " not installed (call '", prg, "')", sep=""))
  else
    return(NULL)
}

if (WINDOWS)
  stop(nbErrorMsg("mcupgma build on Windows is currently not supported"))

# Check if required prequisites for mcupgma are all available
# (as in SystemRequirements in DESCRIPTION + others)
cnt <- c(nbCheckPrg("make", "(GNU) make"),
         nbCheckPrg("bash", "Bash shell"),
         nbCheckPrg("perl", "Perl interpreter"),
         nbCheckPrg("gzip", "GZip packer"),
         nbCheckPrg("echo", "echo shell command"))

# Any missings?
if (length(cnt) > 0) {
  # List formating could be done more elegant...
  cnt[1] <- paste0("- ", cnt[1])
  stop(nbErrorMsg(gsub("\n$", "", perl = TRUE,
                       paste0(cnt, sep = "\n", collapse = "- "))))
}


#print(paste("R_PACKAGE_NAME", R_PACKAGE_NAME))
print(paste("R_PACKAGE_SOURCE", R_PACKAGE_SOURCE))
print(paste("R_PACKAGE_DIR", R_PACKAGE_DIR))
print(paste("R_ARCH", R_ARCH))
print(paste("SHLIB_EXT", SHLIB_EXT))

## mcupgma lying in src folder
src_mcupgma  <- file.path(R_PACKAGE_SOURCE, 'src', 'mcupgma')

# Destination path of R package
dest_lib <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))

## MCUPGMA folder (binaries, scripts and additional files)
dest_mcupgma <- file.path(R_PACKAGE_DIR, 'mcupgma')

print(paste("INSTALL PATH: ",
            file.path(src_mcupgma, "scripts", "install_path.mk")))

# Original MCUPGMA Makefile: create installation path as separate file with
# only export of variable INSTALL_PATH.
# @TODO Replace with R file writing.
system2("echo", paste0("export INSTALL_PATH := ",
                       dest_mcupgma,
                       " > ",
                       file.path(src_mcupgma, "scripts", "install_path.mk")))

## Build mcupgma package (package builds itself in own folder)
retMake1 <- system2(command = "make",                ## Clean
                    args = c(paste0("-C ", src_mcupgma),
                             "-f Makefile", 
                             "clean"))

## Stop on error (default)
if (retMake1 != 0)
  stop(nbErrorMsg("Build mcupgma: make clean failed"))

## Build all.
retMake2 <- system2(command = "make",
                    args = c(paste0("-C ", src_mcupgma),
                             "-f Makefile", 
                             "all"))

if (retMake2 != 0)
  stop(nbErrorMsg("Build mcupgma: make all failed"))

## TODO Add checks.

print(paste("SRC: ", src_mcupgma))
print(paste("DEST:", dest_mcupgma))

if (!dir.create(dest_mcupgma, recursive = TRUE, showWarnings = TRUE))
  stop(nbErrorMsg(paste("Unable to create:", dest_mcupgma)))

## Executables created in two directories (clustering_round and
## clustering_utils) plus some Shell/Perl scripts in directory
## scripts.
## (Many file probes here, incl. dublicates, but less than 20
## files so no problem)
for (path in c(file.path('clustering_round', 'bin'),
               file.path('clustering_util', 'bin'),
               file.path('scripts'))) {
  ## Source folder
  copy_path <- file.path(src_mcupgma, path, '*')
  
  ## All files in folder
  files <- Sys.glob(copy_path)

  ## Skip all sub-directories, all files in flat folder.
  files <- files[which(sapply(files, function(x) file.info(x)$isdir) != TRUE)]
  
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

  print(paste("INSTALL FROM:", copy_path, "TO:", dest_mcupgma, "FILES:"))
##  print(files)
  ## Write all files. (return is logical for each element)
  retCopy <- file.copy(files, dest_mcupgma, overwrite = TRUE)

  if (!all(retCopy)) {
    broken <- files[which(retCopy == FALSE)]
    print("FILES UNABLE TO COPY:")
    print(broken)
    stop(nbErrorMsg("Errors during mcupgma file copies (see list above)"))
  }
}

# @TODO Unused?
files <- Sys.glob(paste0(dest_mcupgma, "/*"))

# This fix should be in Makevars, but as Gnu make syntax only allowed on
# Windows, it hastingly added here...
# @TODO Remove
sapply(Sys.glob(paste0("*", SHLIB_EXT)), function(file) {
    if (file != tolower(file)) {
#        print(paste("Create:", tolower(file), "from", file))
        file.copy(c(file), c(tolower(file)), overwrite = TRUE)
    }
})

## Copy own (non-mcugpma) libraries (basically only netboost.so) to dest.
## Basically the code from R manual
dir.create(dest_lib, recursive = TRUE, showWarnings = FALSE)
files <- Sys.glob(paste0("*", SHLIB_EXT))
file.copy(files, dest_lib, overwrite = TRUE)
print(files)

if(file.exists("symbols.rds"))
  file.copy("symbols.rds", dest_lib, overwrite = TRUE)
