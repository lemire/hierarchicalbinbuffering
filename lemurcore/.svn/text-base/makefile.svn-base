#
.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h
#
#

DEBUGFLAGS =  -g3 -Wall -O2 -fpermissive 
#
%.o : %.cpp %.h common.h molasseexception.h
	echo ****compiling $<
	g++ -c $(DEBUGFLAGS)  $< -I../function

#####################
LemurCoreObjects = mathutil.o chunkeddatacube.o  datacube.o  ramcube.o molasseexception.o  optimizationoptions.o  pyrps.o fileutil.o
#
lemurcore.a : $(LemurCoreObjects)
	ar rs lemurcore.a $?  

all: lemurcore.a tags
	
tags:  
	ls ../* -1d | xargs ctags -aR

.PHONY : clean tags

clean :
	rm -f  *.o
	rm -f	 *.a
	rm -f tags

