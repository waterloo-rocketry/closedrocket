###########################################################################
## Makefile generated for component 'GNC_codegen'. 
## 
## Makefile     : GNC_codegen_rtw.mk
## Generated on : Fri Jun 05 20:32:15 2026
## Final product: ./GNC_codegen.a
## Product type : static-library
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile
# MODELLIB                Static library target

PRODUCT_NAME              = GNC_codegen
MAKEFILE                  = GNC_codegen_rtw.mk
MATLAB_ROOT               = /Applications/MATLAB_R2025b.app
MATLAB_BIN                = /Applications/MATLAB_R2025b.app/bin
MATLAB_ARCH_BIN           = $(MATLAB_BIN)/maca64
START_DIR                 = /Users/shiminghe/Desktop/Waterloo_Rocketry/closedrocket/embedded-coder
TGT_FCN_LIB               = ISO_C
SOLVER_OBJ                = 
CLASSIC_INTERFACE         = 0
MODEL_HAS_DYNAMICALLY_LOADED_SFCNS = 
RELATIVE_PATH_TO_ANCHOR   = ../../..
C_STANDARD_OPTS           = -fno-common -fexceptions -fwrapv
CPP_STANDARD_OPTS         = -std=c++14 -fno-common -fexceptions -fwrapv
MODELLIB                  = GNC_codegen.a

###########################################################################
## TOOLCHAIN SPECIFICATIONS
###########################################################################

# Toolchain Name:          Xcode with Clang | gmake (64-bit Mac)
# Supported Version(s):    
# ToolchainInfo Version:   2025b
# Specification Revision:  1.0
# 
#-------------------------------------------
# Macros assumed to be defined elsewhere
#-------------------------------------------

# C_STANDARD_OPTS
# CPP_STANDARD_OPTS

#-----------
# MACROS
#-----------

ARCHS                               = arm64
XCODE_SDK_VER                       = $(shell perl $(MATLAB_ROOT)/rtw/c/tools/macsdkver.pl)
XCODE_SDK                           = MacOSX$(XCODE_SDK_VER).sdk
XCODE_DEVEL_DIR                     = $(shell xcode-select -print-path)
XCODE_SDK_ROOT                      = $(XCODE_DEVEL_DIR)/Platforms/MacOSX.platform/Developer/SDKs/$(XCODE_SDK)
XCODE_LD_VERSION                    = $(shell xcrun ld -v 2>&1 >/dev/null | head -n 1)
XCODE_LD_VERSION_IS_LD64            = $(findstring ld64-,$(XCODE_LD_VERSION))
XCODE_DYLD_NO_WARN_DUPLICATE_LDFLAG = -Wl,-no_warn_duplicate_libraries
NO_WARN_DUPLICATE_LIBRARIES         = $(if $(XCODE_LD_VERSION_IS_LD64),,$(XCODE_DYLD_NO_WARN_DUPLICATE_LDFLAG))

TOOLCHAIN_SRCS = 
TOOLCHAIN_INCS = 
TOOLCHAIN_LIBS = 

#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: Clang C Compiler
CC = xcrun clang

# Linker: Clang Linker
LD = xcrun clang++

# C++ Compiler: Clang C++ Compiler
CPP = xcrun clang++

# C++ Linker: Clang C++ Linker
CPP_LD = xcrun clang++

# Archiver: Clang Archiver
AR = xcrun ar

# MEX Tool: MEX Tool
MEX_PATH = $(MATLAB_ARCH_BIN)
MEX = "$(MEX_PATH)/mex"

# Download: Download
DOWNLOAD =

# Execute: Execute
EXECUTE = $(PRODUCT)

# Builder: GMAKE Utility
MAKE_PATH = %MATLAB%/bin/maca64
MAKE = "$(MAKE_PATH)/gmake"


#-------------------------
# Directives/Utilities
#-------------------------

CDEBUG              = -g
C_OUTPUT_FLAG       = -o
LDDEBUG             = -g
OUTPUT_FLAG         = -o
CPPDEBUG            = -g
CPP_OUTPUT_FLAG     = -o
CPPLDDEBUG          = -g
OUTPUT_FLAG         = -o
ARDEBUG             =
STATICLIB_OUTPUT_FLAG =
MEX_DEBUG           = -g
RM                  = @rm -f
ECHO                = @echo
MV                  = @mv
RUN                 =

#--------------------------------------
# "Faster Runs" Build Configuration
#--------------------------------------

ARFLAGS              = ruvs
CFLAGS               = -c -isysroot $(XCODE_SDK_ROOT) -arch $(ARCHS) $(C_STANDARD_OPTS) -mmacosx-version-min=12.0 \
                       -O3 -ffp-contract=off
CPPFLAGS             = -c -isysroot $(XCODE_SDK_ROOT) -arch $(ARCHS) $(CPP_STANDARD_OPTS) -mmacosx-version-min=12.0 \
                       -O3 -ffp-contract=off
CPP_LDFLAGS          = -arch $(ARCHS) -isysroot $(XCODE_SDK_ROOT) -Wl,-rpath,@executable_path -Wl,-rpath,@executable_path/$(RELATIVE_PATH_TO_ANCHOR) $(NO_WARN_DUPLICATE_LIBRARIES)
CPP_SHAREDLIB_LDFLAGS  = -arch $(ARCHS)  \
                         -dynamiclib -install_name @rpath/$(notdir $(PRODUCT)) -isysroot $(XCODE_SDK_ROOT) \
                         -Wl,$(LD_NAMESPACE) $(LD_UNDEFS) \
                         $(NO_WARN_DUPLICATE_LIBRARIES)
DOWNLOAD_FLAGS       =
EXECUTE_FLAGS        =
LDFLAGS              = -arch $(ARCHS) -isysroot $(XCODE_SDK_ROOT) -Wl,-rpath,@executable_path -Wl,-rpath,@executable_path/$(RELATIVE_PATH_TO_ANCHOR) $(NO_WARN_DUPLICATE_LIBRARIES)
MEX_CPPFLAGS         =
MEX_CPPLDFLAGS       =
MEX_CFLAGS           =
MEX_LDFLAGS          =
MAKE_FLAGS           = -j $(MAX_MAKE_JOBS) -l $(MAX_MAKE_LOAD_AVG) -f $(MAKEFILE)
SHAREDLIB_LDFLAGS    = -arch $(ARCHS)  \
                       -dynamiclib -install_name @rpath/$(notdir $(PRODUCT)) -isysroot $(XCODE_SDK_ROOT) \
                       -Wl,$(LD_NAMESPACE) $(LD_UNDEFS) \
                       $(NO_WARN_DUPLICATE_LIBRARIES)



###########################################################################
## OUTPUT INFO
###########################################################################

PRODUCT = ./GNC_codegen.a
PRODUCT_TYPE = "static-library"
BUILD_TYPE = "Static Library"

###########################################################################
## INCLUDE PATHS
###########################################################################

INCLUDES_BUILDINFO = -I$(START_DIR)/codegen/lib/GNC_codegen -I$(START_DIR) -I$(MATLAB_ROOT)/extern/include

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_CUSTOM = 
DEFINES_STANDARD = -DMODEL=GNC_codegen

DEFINES = $(DEFINES_CUSTOM) $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = $(START_DIR)/codegen/lib/GNC_codegen/GNC_codegen_data.c $(START_DIR)/codegen/lib/GNC_codegen/rt_nonfinite.c $(START_DIR)/codegen/lib/GNC_codegen/rtGetNaN.c $(START_DIR)/codegen/lib/GNC_codegen/rtGetInf.c $(START_DIR)/codegen/lib/GNC_codegen/GNC_codegen_initialize.c $(START_DIR)/codegen/lib/GNC_codegen/GNC_codegen_terminate.c $(START_DIR)/codegen/lib/GNC_codegen/controller_codegen_entry.c $(START_DIR)/codegen/lib/GNC_codegen/eye.c $(START_DIR)/codegen/lib/GNC_codegen/navigation_codegen_entry.c $(START_DIR)/codegen/lib/GNC_codegen/norm.c $(START_DIR)/codegen/lib/GNC_codegen/airdata_atmos.c $(START_DIR)/codegen/lib/GNC_codegen/mpower.c $(START_DIR)/codegen/lib/GNC_codegen/atan2.c $(START_DIR)/codegen/lib/GNC_codegen/inv.c $(START_DIR)/codegen/lib/GNC_codegen/ekf_correct.c $(START_DIR)/codegen/lib/GNC_codegen/svd.c $(START_DIR)/codegen/lib/GNC_codegen/xzlangeM.c $(START_DIR)/codegen/lib/GNC_codegen/xnrm2.c $(START_DIR)/codegen/lib/GNC_codegen/pad_filter.c $(START_DIR)/codegen/lib/GNC_codegen/dynamics.c $(START_DIR)/codegen/lib/GNC_codegen/dynamics_jacobian.c $(START_DIR)/codegen/lib/GNC_codegen/diag.c $(START_DIR)/codegen/lib/GNC_codegen/xrotg.c $(START_DIR)/codegen/lib/GNC_codegen/xzlascl.c

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = GNC_codegen_data.o rt_nonfinite.o rtGetNaN.o rtGetInf.o GNC_codegen_initialize.o GNC_codegen_terminate.o controller_codegen_entry.o eye.o navigation_codegen_entry.o norm.o airdata_atmos.o mpower.o atan2.o inv.o ekf_correct.o svd.o xzlangeM.o xnrm2.o pad_filter.o dynamics.o dynamics_jacobian.o diag.o xrotg.o xzlascl.o

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

SYSTEM_LIBS =  -lm

###########################################################################
## ADDITIONAL TOOLCHAIN FLAGS
###########################################################################

#---------------
# C Compiler
#---------------

CFLAGS_BASIC = $(DEFINES) $(INCLUDES)

CFLAGS += $(CFLAGS_BASIC)

#-----------------
# C++ Compiler
#-----------------

CPPFLAGS_BASIC = $(DEFINES) $(INCLUDES)

CPPFLAGS += $(CPPFLAGS_BASIC)

###########################################################################
## INLINED COMMANDS
###########################################################################

###########################################################################
## PHONY TARGETS
###########################################################################

.PHONY : all build clean info prebuild download execute


all : build
	@echo "### Successfully generated all binary outputs."


build : prebuild $(PRODUCT)


prebuild : 


download : $(PRODUCT)


execute : download


###########################################################################
## FINAL TARGET
###########################################################################

#---------------------------------
# Create a static library         
#---------------------------------

$(PRODUCT) : $(OBJS) $(PREBUILT_OBJS)
	@echo "### Creating static library "$(PRODUCT)" ..."
	$(AR) $(ARFLAGS)  $(PRODUCT) $(OBJS)
	@echo "### Created: $(PRODUCT)"


###########################################################################
## INTERMEDIATE TARGETS
###########################################################################

#---------------------
# SOURCE-TO-OBJECT
#---------------------

%.o : %.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : %.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : %.cc
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : %.CC
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : %.cp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : %.CPP
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : %.c++
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : %.C++
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : %.cxx
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : %.CXX
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : %.C
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.cc
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.CC
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.cp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.CPP
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.c++
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.C++
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.cxx
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.CXX
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.C
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.cc
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.CC
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.cp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.CPP
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.c++
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.C++
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.cxx
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.CXX
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/codegen/lib/GNC_codegen/%.C
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.cc
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.CC
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.cp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.CPP
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.c++
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.C++
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.cxx
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.CXX
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.C
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


GNC_codegen_data.o : $(START_DIR)/codegen/lib/GNC_codegen/GNC_codegen_data.c
	$(CC) $(CFLAGS) -o "$@" "$<"


rt_nonfinite.o : $(START_DIR)/codegen/lib/GNC_codegen/rt_nonfinite.c
	$(CC) $(CFLAGS) -o "$@" "$<"


rtGetNaN.o : $(START_DIR)/codegen/lib/GNC_codegen/rtGetNaN.c
	$(CC) $(CFLAGS) -o "$@" "$<"


rtGetInf.o : $(START_DIR)/codegen/lib/GNC_codegen/rtGetInf.c
	$(CC) $(CFLAGS) -o "$@" "$<"


GNC_codegen_initialize.o : $(START_DIR)/codegen/lib/GNC_codegen/GNC_codegen_initialize.c
	$(CC) $(CFLAGS) -o "$@" "$<"


GNC_codegen_terminate.o : $(START_DIR)/codegen/lib/GNC_codegen/GNC_codegen_terminate.c
	$(CC) $(CFLAGS) -o "$@" "$<"


controller_codegen_entry.o : $(START_DIR)/codegen/lib/GNC_codegen/controller_codegen_entry.c
	$(CC) $(CFLAGS) -o "$@" "$<"


eye.o : $(START_DIR)/codegen/lib/GNC_codegen/eye.c
	$(CC) $(CFLAGS) -o "$@" "$<"


navigation_codegen_entry.o : $(START_DIR)/codegen/lib/GNC_codegen/navigation_codegen_entry.c
	$(CC) $(CFLAGS) -o "$@" "$<"


norm.o : $(START_DIR)/codegen/lib/GNC_codegen/norm.c
	$(CC) $(CFLAGS) -o "$@" "$<"


airdata_atmos.o : $(START_DIR)/codegen/lib/GNC_codegen/airdata_atmos.c
	$(CC) $(CFLAGS) -o "$@" "$<"


mpower.o : $(START_DIR)/codegen/lib/GNC_codegen/mpower.c
	$(CC) $(CFLAGS) -o "$@" "$<"


atan2.o : $(START_DIR)/codegen/lib/GNC_codegen/atan2.c
	$(CC) $(CFLAGS) -o "$@" "$<"


inv.o : $(START_DIR)/codegen/lib/GNC_codegen/inv.c
	$(CC) $(CFLAGS) -o "$@" "$<"


ekf_correct.o : $(START_DIR)/codegen/lib/GNC_codegen/ekf_correct.c
	$(CC) $(CFLAGS) -o "$@" "$<"


svd.o : $(START_DIR)/codegen/lib/GNC_codegen/svd.c
	$(CC) $(CFLAGS) -o "$@" "$<"


xzlangeM.o : $(START_DIR)/codegen/lib/GNC_codegen/xzlangeM.c
	$(CC) $(CFLAGS) -o "$@" "$<"


xnrm2.o : $(START_DIR)/codegen/lib/GNC_codegen/xnrm2.c
	$(CC) $(CFLAGS) -o "$@" "$<"


pad_filter.o : $(START_DIR)/codegen/lib/GNC_codegen/pad_filter.c
	$(CC) $(CFLAGS) -o "$@" "$<"


dynamics.o : $(START_DIR)/codegen/lib/GNC_codegen/dynamics.c
	$(CC) $(CFLAGS) -o "$@" "$<"


dynamics_jacobian.o : $(START_DIR)/codegen/lib/GNC_codegen/dynamics_jacobian.c
	$(CC) $(CFLAGS) -o "$@" "$<"


diag.o : $(START_DIR)/codegen/lib/GNC_codegen/diag.c
	$(CC) $(CFLAGS) -o "$@" "$<"


xrotg.o : $(START_DIR)/codegen/lib/GNC_codegen/xrotg.c
	$(CC) $(CFLAGS) -o "$@" "$<"


xzlascl.o : $(START_DIR)/codegen/lib/GNC_codegen/xzlascl.c
	$(CC) $(CFLAGS) -o "$@" "$<"


###########################################################################
## DEPENDENCIES
###########################################################################

$(ALL_OBJS) : rtw_proj.tmw $(MAKEFILE)


###########################################################################
## MISCELLANEOUS TARGETS
###########################################################################

info : 
	@echo "### PRODUCT = $(PRODUCT)"
	@echo "### PRODUCT_TYPE = $(PRODUCT_TYPE)"
	@echo "### BUILD_TYPE = $(BUILD_TYPE)"
	@echo "### INCLUDES = $(INCLUDES)"
	@echo "### DEFINES = $(DEFINES)"
	@echo "### ALL_SRCS = $(ALL_SRCS)"
	@echo "### ALL_OBJS = $(ALL_OBJS)"
	@echo "### LIBS = $(LIBS)"
	@echo "### MODELREF_LIBS = $(MODELREF_LIBS)"
	@echo "### SYSTEM_LIBS = $(SYSTEM_LIBS)"
	@echo "### TOOLCHAIN_LIBS = $(TOOLCHAIN_LIBS)"
	@echo "### CFLAGS = $(CFLAGS)"
	@echo "### LDFLAGS = $(LDFLAGS)"
	@echo "### SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS)"
	@echo "### CPPFLAGS = $(CPPFLAGS)"
	@echo "### CPP_LDFLAGS = $(CPP_LDFLAGS)"
	@echo "### CPP_SHAREDLIB_LDFLAGS = $(CPP_SHAREDLIB_LDFLAGS)"
	@echo "### ARFLAGS = $(ARFLAGS)"
	@echo "### MEX_CFLAGS = $(MEX_CFLAGS)"
	@echo "### MEX_CPPFLAGS = $(MEX_CPPFLAGS)"
	@echo "### MEX_LDFLAGS = $(MEX_LDFLAGS)"
	@echo "### MEX_CPPLDFLAGS = $(MEX_CPPLDFLAGS)"
	@echo "### DOWNLOAD_FLAGS = $(DOWNLOAD_FLAGS)"
	@echo "### EXECUTE_FLAGS = $(EXECUTE_FLAGS)"
	@echo "### MAKE_FLAGS = $(MAKE_FLAGS)"


clean : 
	$(ECHO) "### Deleting all derived files ..."
	$(RM) $(PRODUCT)
	$(RM) $(ALL_OBJS)
	$(ECHO) "### Deleted all derived files."


