###########################################################################
## Makefile generated for component 'navigation_codegen_entry'. 
## 
## Makefile     : navigation_codegen_entry_rtw.mk
## Generated on : Sat May 16 19:50:16 2026
## Final product: .\navigation_codegen_entry.lib
## Product type : static-library
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile
# COMPILER_COMMAND_FILE   Compiler command listing model reference header paths
# CMD_FILE                Command file
# MODELLIB                Static library target

PRODUCT_NAME              = navigation_codegen_entry
MAKEFILE                  = navigation_codegen_entry_rtw.mk
MATLAB_ROOT               = C:\PROGRA~1\MATLAB\R2025b
MATLAB_BIN                = C:\PROGRA~1\MATLAB\R2025b\bin
MATLAB_ARCH_BIN           = $(MATLAB_BIN)\win64
START_DIR                 = C:\Users\trist\Documents\GitHub\simulink-canards
TGT_FCN_LIB               = ISO_C
SOLVER_OBJ                = 
CLASSIC_INTERFACE         = 0
MODEL_HAS_DYNAMICALLY_LOADED_SFCNS = 
RELATIVE_PATH_TO_ANCHOR   = ..\..\..
COMPILER_COMMAND_FILE     = navigation_codegen_entry_rtw_comp.rsp
CMD_FILE                  = navigation_codegen_entry_rtw.rsp
C_STANDARD_OPTS           = 
CPP_STANDARD_OPTS         = 
NODEBUG                   = 1
MODELLIB                  = navigation_codegen_entry.lib

###########################################################################
## TOOLCHAIN SPECIFICATIONS
###########################################################################

# Toolchain Name:          Microsoft Visual C++ 2022 v17.0 | nmake (64-bit Windows)
# Supported Version(s):    17.0
# ToolchainInfo Version:   2025b
# Specification Revision:  1.0
# 
#-------------------------------------------
# Macros assumed to be defined elsewhere
#-------------------------------------------

# C_STANDARD_OPTS
# CPP_STANDARD_OPTS
# NODEBUG
# cvarsdll
# cvarsmt
# conlibsmt
# ldebug
# conflags
# cflags

#-----------
# MACROS
#-----------

MW_EXTERNLIB_DIR    = $(MATLAB_ROOT)\extern\lib\win64\microsoft
MW_LIB_DIR          = $(MATLAB_ROOT)\lib\win64
CPU                 = AMD64
APPVER              = 5.02
CVARSFLAG           = $(cvarsmt)
CFLAGS_ADDITIONAL   = -D_CRT_SECURE_NO_WARNINGS
CPPFLAGS_ADDITIONAL = -EHs -D_CRT_SECURE_NO_WARNINGS /wd4251 /Zc:__cplusplus
LIBS_TOOLCHAIN      = $(conlibs)

TOOLCHAIN_SRCS = 
TOOLCHAIN_INCS = 
TOOLCHAIN_LIBS = 

#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: Microsoft Visual C Compiler
CC = cl

# Linker: Microsoft Visual C Linker
LD = link

# C++ Compiler: Microsoft Visual C++ Compiler
CPP = cl

# C++ Linker: Microsoft Visual C++ Linker
CPP_LD = link

# Archiver: Microsoft Visual C/C++ Archiver
AR = lib

# MEX Tool: MEX Tool
MEX_PATH = $(MATLAB_ARCH_BIN)
MEX = "$(MEX_PATH)\mex"

# Download: Download
DOWNLOAD =

# Execute: Execute
EXECUTE = $(PRODUCT)

# Builder: NMAKE Utility
MAKE = nmake


#-------------------------
# Directives/Utilities
#-------------------------

CDEBUG              = -Zi
C_OUTPUT_FLAG       = -Fo
LDDEBUG             = /DEBUG
OUTPUT_FLAG         = -out:
CPPDEBUG            = -Zi
CPP_OUTPUT_FLAG     = -Fo
CPPLDDEBUG          = /DEBUG
OUTPUT_FLAG         = -out:
ARDEBUG             =
STATICLIB_OUTPUT_FLAG = -out:
MEX_DEBUG           = -g
RM                  = @del
ECHO                = @echo
MV                  = @ren
RUN                 = @cmd /C

#--------------------------------------
# "Faster Runs" Build Configuration
#--------------------------------------

ARFLAGS              = /nologo
CFLAGS               = $(cflags) $(CVARSFLAG) $(CFLAGS_ADDITIONAL) $(C_STANDARD_OPTS) \
                       /O2 /Oy-
CPPFLAGS             = /TP $(cflags) $(CVARSFLAG) $(CPPFLAGS_ADDITIONAL) $(CPP_STANDARD_OPTS) \
                       /O2 /Oy-
CPP_LDFLAGS          = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN)
CPP_SHAREDLIB_LDFLAGS  = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN) \
                         -dll -def:$(DEF_FILE)
DOWNLOAD_FLAGS       =
EXECUTE_FLAGS        =
LDFLAGS              = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN)
MEX_CPPFLAGS         =
MEX_CPPLDFLAGS       =
MEX_CFLAGS           =
MEX_LDFLAGS          =
MAKE_FLAGS           = -f $(MAKEFILE)
SHAREDLIB_LDFLAGS    = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN) \
                       -dll -def:$(DEF_FILE)



###########################################################################
## OUTPUT INFO
###########################################################################

PRODUCT = .\navigation_codegen_entry.lib
PRODUCT_TYPE = "static-library"
BUILD_TYPE = "Static Library"

###########################################################################
## INCLUDE PATHS
###########################################################################

INCLUDES_BUILDINFO = 

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_CUSTOM = 
DEFINES_STANDARD = -DMODEL=navigation_codegen_entry

DEFINES = $(DEFINES_CUSTOM) $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = $(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry_data.c $(START_DIR)\codegen\lib\navigation_codegen_entry\rt_nonfinite.c $(START_DIR)\codegen\lib\navigation_codegen_entry\rtGetNaN.c $(START_DIR)\codegen\lib\navigation_codegen_entry\rtGetInf.c $(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry_initialize.c $(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry_terminate.c $(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry.c $(START_DIR)\codegen\lib\navigation_codegen_entry\norm.c $(START_DIR)\codegen\lib\navigation_codegen_entry\eye.c $(START_DIR)\codegen\lib\navigation_codegen_entry\airdata_atmos.c $(START_DIR)\codegen\lib\navigation_codegen_entry\mpower.c $(START_DIR)\codegen\lib\navigation_codegen_entry\atan2.c $(START_DIR)\codegen\lib\navigation_codegen_entry\inv.c $(START_DIR)\codegen\lib\navigation_codegen_entry\ekf_correct.c $(START_DIR)\codegen\lib\navigation_codegen_entry\svd.c $(START_DIR)\codegen\lib\navigation_codegen_entry\xzlangeM.c $(START_DIR)\codegen\lib\navigation_codegen_entry\xnrm2.c $(START_DIR)\codegen\lib\navigation_codegen_entry\pad_filter.c $(START_DIR)\codegen\lib\navigation_codegen_entry\dynamics.c $(START_DIR)\codegen\lib\navigation_codegen_entry\dynamics_jacobian.c $(START_DIR)\codegen\lib\navigation_codegen_entry\xrotg.c $(START_DIR)\codegen\lib\navigation_codegen_entry\xzlascl.c

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = navigation_codegen_entry_data.obj rt_nonfinite.obj rtGetNaN.obj rtGetInf.obj navigation_codegen_entry_initialize.obj navigation_codegen_entry_terminate.obj navigation_codegen_entry.obj norm.obj eye.obj airdata_atmos.obj mpower.obj atan2.obj inv.obj ekf_correct.obj svd.obj xzlangeM.obj xnrm2.obj pad_filter.obj dynamics.obj dynamics_jacobian.obj xrotg.obj xzlascl.obj

ALL_OBJS = $(OBJS)

###########################################################################
## PREBUILT OBJECT FILES
###########################################################################

PREBUILT_OBJS = 

###########################################################################
## LIBRARIES
###########################################################################

LIBS = 

###########################################################################
## SYSTEM LIBRARIES
###########################################################################

SYSTEM_LIBS = /LIBPATH:"$(MATLAB_ROOT)\bin\win64" "$(MATLAB_ROOT)\bin\win64\libiomp5md.lib"

###########################################################################
## ADDITIONAL TOOLCHAIN FLAGS
###########################################################################

#---------------
# C Compiler
#---------------

CFLAGS_ = /source-charset:utf-8
CFLAGS_OPTS = /openmp /wd4101
CFLAGS_BASIC = $(DEFINES) @$(COMPILER_COMMAND_FILE)

CFLAGS = $(CFLAGS) $(CFLAGS_) $(CFLAGS_OPTS) $(CFLAGS_BASIC)

#-----------------
# C++ Compiler
#-----------------

CPPFLAGS_ = /source-charset:utf-8
CPPFLAGS_OPTS = /openmp /wd4101
CPPFLAGS_BASIC = $(DEFINES) @$(COMPILER_COMMAND_FILE)

CPPFLAGS = $(CPPFLAGS) $(CPPFLAGS_) $(CPPFLAGS_OPTS) $(CPPFLAGS_BASIC)

#---------------
# C++ Linker
#---------------

CPP_LDFLAGS_ = /nodefaultlib:vcomp  

CPP_LDFLAGS = $(CPP_LDFLAGS) $(CPP_LDFLAGS_)

#------------------------------
# C++ Shared Library Linker
#------------------------------

CPP_SHAREDLIB_LDFLAGS_ = /nodefaultlib:vcomp  

CPP_SHAREDLIB_LDFLAGS = $(CPP_SHAREDLIB_LDFLAGS) $(CPP_SHAREDLIB_LDFLAGS_)

#-----------
# Linker
#-----------

LDFLAGS_ = /nodefaultlib:vcomp  

LDFLAGS = $(LDFLAGS) $(LDFLAGS_)

#--------------------------
# Shared Library Linker
#--------------------------

SHAREDLIB_LDFLAGS_ = /nodefaultlib:vcomp  

SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS) $(SHAREDLIB_LDFLAGS_)

###########################################################################
## INLINED COMMANDS
###########################################################################


!include $(MATLAB_ROOT)\rtw\c\tools\vcdefs.mak


###########################################################################
## PHONY TARGETS
###########################################################################

.PHONY : all build clean info prebuild download execute set_environment_variables


all : build
	@cmd /C "@echo ### Successfully generated all binary outputs."


build : set_environment_variables prebuild $(PRODUCT)


prebuild : 


download : $(PRODUCT)


execute : download


set_environment_variables : 
	@set INCLUDE=$(INCLUDES);$(INCLUDE)
	@set LIB=$(LIB)


###########################################################################
## FINAL TARGET
###########################################################################

#---------------------------------
# Create a static library         
#---------------------------------

$(PRODUCT) : $(OBJS) $(PREBUILT_OBJS)
	@cmd /C "@echo ### Creating static library "$(PRODUCT)" ..."
	$(AR) $(ARFLAGS) -out:$(PRODUCT) @$(CMD_FILE)
	@cmd /C "@echo ### Created: $(PRODUCT)"


###########################################################################
## INTERMEDIATE TARGETS
###########################################################################

#---------------------
# SOURCE-TO-OBJECT
#---------------------

.c.obj:
	$(CC) $(CFLAGS) -Fo"$@" "$<"


.cpp.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


.cc.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


.cxx.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(RELATIVE_PATH_TO_ANCHOR)}.c.obj:
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(RELATIVE_PATH_TO_ANCHOR)}.cpp.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(RELATIVE_PATH_TO_ANCHOR)}.cc.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(RELATIVE_PATH_TO_ANCHOR)}.cxx.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(START_DIR)\codegen\lib\navigation_codegen_entry}.c.obj:
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(START_DIR)\codegen\lib\navigation_codegen_entry}.cpp.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(START_DIR)\codegen\lib\navigation_codegen_entry}.cc.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(START_DIR)\codegen\lib\navigation_codegen_entry}.cxx.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(START_DIR)}.c.obj:
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(START_DIR)}.cpp.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(START_DIR)}.cc.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(START_DIR)}.cxx.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


navigation_codegen_entry_data.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry_data.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry_data.c"


rt_nonfinite.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\rt_nonfinite.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\rt_nonfinite.c"


rtGetNaN.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\rtGetNaN.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\rtGetNaN.c"


rtGetInf.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\rtGetInf.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\rtGetInf.c"


navigation_codegen_entry_initialize.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry_initialize.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry_initialize.c"


navigation_codegen_entry_terminate.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry_terminate.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry_terminate.c"


navigation_codegen_entry.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\navigation_codegen_entry.c"


norm.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\norm.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\norm.c"


eye.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\eye.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\eye.c"


airdata_atmos.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\airdata_atmos.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\airdata_atmos.c"


mpower.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\mpower.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\mpower.c"


atan2.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\atan2.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\atan2.c"


inv.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\inv.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\inv.c"


ekf_correct.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\ekf_correct.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\ekf_correct.c"


svd.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\svd.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\svd.c"


xzlangeM.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\xzlangeM.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\xzlangeM.c"


xnrm2.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\xnrm2.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\xnrm2.c"


pad_filter.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\pad_filter.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\pad_filter.c"


dynamics.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\dynamics.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\dynamics.c"


dynamics_jacobian.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\dynamics_jacobian.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\dynamics_jacobian.c"


xrotg.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\xrotg.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\xrotg.c"


xzlascl.obj : "$(START_DIR)\codegen\lib\navigation_codegen_entry\xzlascl.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\navigation_codegen_entry\xzlascl.c"


###########################################################################
## DEPENDENCIES
###########################################################################

$(ALL_OBJS) : rtw_proj.tmw $(COMPILER_COMMAND_FILE) $(MAKEFILE)


###########################################################################
## MISCELLANEOUS TARGETS
###########################################################################

info : 
	@cmd /C "@echo ### PRODUCT = $(PRODUCT)"
	@cmd /C "@echo ### PRODUCT_TYPE = $(PRODUCT_TYPE)"
	@cmd /C "@echo ### BUILD_TYPE = $(BUILD_TYPE)"
	@cmd /C "@echo ### INCLUDES = $(INCLUDES)"
	@cmd /C "@echo ### DEFINES = $(DEFINES)"
	@cmd /C "@echo ### ALL_SRCS = $(ALL_SRCS)"
	@cmd /C "@echo ### ALL_OBJS = $(ALL_OBJS)"
	@cmd /C "@echo ### LIBS = $(LIBS)"
	@cmd /C "@echo ### MODELREF_LIBS = $(MODELREF_LIBS)"
	@cmd /C "@echo ### SYSTEM_LIBS = $(SYSTEM_LIBS)"
	@cmd /C "@echo ### TOOLCHAIN_LIBS = $(TOOLCHAIN_LIBS)"
	@cmd /C "@echo ### CFLAGS = $(CFLAGS)"
	@cmd /C "@echo ### LDFLAGS = $(LDFLAGS)"
	@cmd /C "@echo ### SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS)"
	@cmd /C "@echo ### CPPFLAGS = $(CPPFLAGS)"
	@cmd /C "@echo ### CPP_LDFLAGS = $(CPP_LDFLAGS)"
	@cmd /C "@echo ### CPP_SHAREDLIB_LDFLAGS = $(CPP_SHAREDLIB_LDFLAGS)"
	@cmd /C "@echo ### ARFLAGS = $(ARFLAGS)"
	@cmd /C "@echo ### MEX_CFLAGS = $(MEX_CFLAGS)"
	@cmd /C "@echo ### MEX_CPPFLAGS = $(MEX_CPPFLAGS)"
	@cmd /C "@echo ### MEX_LDFLAGS = $(MEX_LDFLAGS)"
	@cmd /C "@echo ### MEX_CPPLDFLAGS = $(MEX_CPPLDFLAGS)"
	@cmd /C "@echo ### DOWNLOAD_FLAGS = $(DOWNLOAD_FLAGS)"
	@cmd /C "@echo ### EXECUTE_FLAGS = $(EXECUTE_FLAGS)"
	@cmd /C "@echo ### MAKE_FLAGS = $(MAKE_FLAGS)"


clean : 
	$(ECHO) "### Deleting all derived files ..."
	@if exist $(PRODUCT) $(RM) $(PRODUCT)
	$(RM) $(ALL_OBJS)
	$(ECHO) "### Deleted all derived files."


