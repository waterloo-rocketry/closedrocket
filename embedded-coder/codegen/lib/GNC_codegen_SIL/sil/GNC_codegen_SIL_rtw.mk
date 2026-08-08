###########################################################################
## Makefile generated for component 'GNC_codegen_SIL'. 
## 
## Makefile     : GNC_codegen_SIL_rtw.mk
## Generated on : Sat Aug 08 13:00:24 2026
## Final product: .\GNC_codegen_SIL.exe
## Product type : executable
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile
# MODELREF_LINK_RSPFILE   Linker command listing model reference link objects
# COMPILER_COMMAND_FILE   Compiler command listing model reference header paths
# CMD_FILE                Command file

PRODUCT_NAME              = GNC_codegen_SIL
MAKEFILE                  = GNC_codegen_SIL_rtw.mk
MATLAB_ROOT               = C:\PROGRA~1\MATLAB\R2025b
MATLAB_BIN                = C:\PROGRA~1\MATLAB\R2025b\bin
MATLAB_ARCH_BIN           = $(MATLAB_BIN)\win64
START_DIR                 = C:\Users\trist\Documents\GitHub\simulink-canards\embedded-coder
TGT_FCN_LIB               = ISO_C
SOLVER_OBJ                = 
CLASSIC_INTERFACE         = 0
MODEL_HAS_DYNAMICALLY_LOADED_SFCNS = 
RELATIVE_PATH_TO_ANCHOR   = ..\..\..\..
MODELREF_LINK_RSPFILE     = GNC_codegen_SIL_rtw_ref.rsp
COMPILER_COMMAND_FILE     = GNC_codegen_SIL_rtw_comp.rsp
CMD_FILE                  = GNC_codegen_SIL_rtw.rsp
C_STANDARD_OPTS           = 
CPP_STANDARD_OPTS         = 
NODEBUG                   = 1

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

PRODUCT = .\GNC_codegen_SIL.exe
PRODUCT_TYPE = "executable"
BUILD_TYPE = "Executable"

###########################################################################
## INCLUDE PATHS
###########################################################################

INCLUDES_BUILDINFO = 

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_ = -DCODER_ASSUMPTIONS_ENABLED=1 -DXIL_SIGNAL_HANDLER=1 -DSIL_DISABLE_SUBNORMAL_SUPPORT=0
DEFINES_CUSTOM = 
DEFINES_OPTS = -DRTIOSTREAM_RX_BUFFER_BYTE_SIZE=50000 -DRTIOSTREAM_TX_BUFFER_BYTE_SIZE=50000 -DCODE_INSTRUMENTATION_ENABLED=1 -DMEM_UNIT_BYTES=1 -DMemUnit_T=uint8_T
DEFINES_STANDARD = -DMODEL=GNC_codegen_SIL

DEFINES = $(DEFINES_) $(DEFINES_CUSTOM) $(DEFINES_OPTS) $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = $(START_DIR)\codegen\lib\GNC_codegen_SIL\sil\xil_interface.c $(START_DIR)\codegen\lib\GNC_codegen_SIL\sil\sil_main.c $(MATLAB_ROOT)\toolbox\coder\rtiostream\src\rtiostreamtcpip\rtiostream_tcpip.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_interface_lib.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_data_stream.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_services.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xilcomms_rtiostream.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_rtiostream.c $(MATLAB_ROOT)\toolbox\coder\rtiostream\src\utils\rtiostream_utils.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\target_io.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\coder_assumptions_app.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\coder_assumptions_data_stream.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\coder_assumptions_rtiostream.c $(START_DIR)\codegen\lib\GNC_codegen_SIL\sil\xil_instrumentation.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\codeinstr_data_stream.c $(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\codeinstr_rtiostream.c $(MATLAB_ROOT)\toolbox\coder\profile\src\host_timer_x86.c $(START_DIR)\codegen\lib\GNC_codegen_SIL\target\_coder_GNC_codegen_SIL_target.c

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = xil_interface.obj sil_main.obj rtiostream_tcpip.obj xil_interface_lib.obj xil_data_stream.obj xil_services.obj xilcomms_rtiostream.obj xil_rtiostream.obj rtiostream_utils.obj target_io.obj coder_assumptions_app.obj coder_assumptions_data_stream.obj coder_assumptions_rtiostream.obj xil_instrumentation.obj codeinstr_data_stream.obj codeinstr_rtiostream.obj host_timer_x86.obj _coder_GNC_codegen_SIL_target.obj

ALL_OBJS = $(OBJS)

###########################################################################
## PREBUILT OBJECT FILES
###########################################################################

PREBUILT_OBJS = 

###########################################################################
## LIBRARIES
###########################################################################

MODELREF_LIBS = C:\Users\trist\Documents\GitHub\simulink-canards\embedded-coder\codegen\lib\GNC_codegen_SIL\instrumented\GNC_codegen_SIL.lib

LIBS = $(START_DIR)\codegen\lib\GNC_codegen_SIL\coderassumptions\lib\GNC_codegen_SIL_ca.lib

###########################################################################
## SYSTEM LIBRARIES
###########################################################################

SYSTEM_LIBS = 

###########################################################################
## ADDITIONAL TOOLCHAIN FLAGS
###########################################################################

#---------------
# C Compiler
#---------------

CFLAGS_ = /source-charset:utf-8
CFLAGS_BASIC = $(DEFINES) @$(COMPILER_COMMAND_FILE)

CFLAGS = $(CFLAGS) $(CFLAGS_) $(CFLAGS_BASIC)

#-----------------
# C++ Compiler
#-----------------

CPPFLAGS_ = /source-charset:utf-8
CPPFLAGS_BASIC = $(DEFINES) @$(COMPILER_COMMAND_FILE)

CPPFLAGS = $(CPPFLAGS) $(CPPFLAGS_) $(CPPFLAGS_BASIC)

###########################################################################
## INLINED COMMANDS
###########################################################################


!include $(MATLAB_ROOT)\rtw\c\tools\vcdefs.mak


###########################################################################
## PHONY TARGETS
###########################################################################

.PHONY : all build buildobj clean info prebuild download execute set_environment_variables


all : build
	@cmd /C "@echo ### Successfully generated all binary outputs."


build : set_environment_variables prebuild $(PRODUCT)


buildobj : set_environment_variables prebuild $(OBJS) $(PREBUILT_OBJS) $(LIBS)
	@cmd /C "@echo ### Successfully generated all binary outputs."


prebuild : 


download : $(PRODUCT)


execute : download
	@cmd /C "@echo ### Invoking postbuild tool "Execute" ..."
	$(EXECUTE) $(EXECUTE_FLAGS)
	@cmd /C "@echo ### Done invoking postbuild tool."


set_environment_variables : 
	@set INCLUDE=$(INCLUDES);$(INCLUDE)
	@set LIB=$(LIB)


###########################################################################
## FINAL TARGET
###########################################################################

#-------------------------------------------
# Create a standalone executable            
#-------------------------------------------

$(PRODUCT) : $(OBJS) $(PREBUILT_OBJS) $(MODELREF_LIBS) $(LIBS)
	@cmd /C "@echo ### Creating standalone executable "$(PRODUCT)" ..."
	$(LD) $(LDFLAGS) -out:$(PRODUCT) @$(CMD_FILE) @$(MODELREF_LINK_RSPFILE) $(LIBS) $(SYSTEM_LIBS) $(TOOLCHAIN_LIBS)
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


{$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c}.c.obj:
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c}.cpp.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c}.cc.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c}.cxx.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\toolbox\coder\profile\src}.c.obj:
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\toolbox\coder\profile\src}.cpp.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\toolbox\coder\profile\src}.cc.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\toolbox\coder\profile\src}.cxx.obj:
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


xil_interface.obj : "$(START_DIR)\codegen\lib\GNC_codegen_SIL\sil\xil_interface.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\GNC_codegen_SIL\sil\xil_interface.c"


sil_main.obj : "$(START_DIR)\codegen\lib\GNC_codegen_SIL\sil\sil_main.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\GNC_codegen_SIL\sil\sil_main.c"


rtiostream_tcpip.obj : "$(MATLAB_ROOT)\toolbox\coder\rtiostream\src\rtiostreamtcpip\rtiostream_tcpip.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\coder\rtiostream\src\rtiostreamtcpip\rtiostream_tcpip.c"


xil_interface_lib.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_interface_lib.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_interface_lib.c"


xil_data_stream.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_data_stream.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_data_stream.c"


xil_services.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_services.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_services.c"


xilcomms_rtiostream.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xilcomms_rtiostream.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xilcomms_rtiostream.c"


xil_rtiostream.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_rtiostream.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\xil_rtiostream.c"


rtiostream_utils.obj : "$(MATLAB_ROOT)\toolbox\coder\rtiostream\src\utils\rtiostream_utils.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\coder\rtiostream\src\utils\rtiostream_utils.c"


target_io.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\target_io.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\target_io.c"


coder_assumptions_app.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\coder_assumptions_app.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\coder_assumptions_app.c"


coder_assumptions_data_stream.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\coder_assumptions_data_stream.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\coder_assumptions_data_stream.c"


coder_assumptions_rtiostream.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\coder_assumptions_rtiostream.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\coder_assumptions_rtiostream.c"


xil_instrumentation.obj : "$(START_DIR)\codegen\lib\GNC_codegen_SIL\sil\xil_instrumentation.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\GNC_codegen_SIL\sil\xil_instrumentation.c"


codeinstr_data_stream.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\codeinstr_data_stream.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\codeinstr_data_stream.c"


codeinstr_rtiostream.obj : "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\codeinstr_rtiostream.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\rtw\targets\pil\c\codeinstr_rtiostream.c"


host_timer_x86.obj : "$(MATLAB_ROOT)\toolbox\coder\profile\src\host_timer_x86.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(MATLAB_ROOT)\toolbox\coder\profile\src\host_timer_x86.c"


_coder_GNC_codegen_SIL_target.obj : "$(START_DIR)\codegen\lib\GNC_codegen_SIL\target\_coder_GNC_codegen_SIL_target.c"
	$(CC) $(CFLAGS) -Fo"$@" "$(START_DIR)\codegen\lib\GNC_codegen_SIL\target\_coder_GNC_codegen_SIL_target.c"


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


