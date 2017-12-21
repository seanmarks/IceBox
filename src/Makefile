# NOTES:
# -g	Extra symbolic debugging info for use with the gdb debugger
# -o	Specify output file name
# -c	Compile into object file
# -Wall	Print all warning messages
#
# $@	Name of the target
# $<	First dependency of the target
# $^	All dependencies of the target
#
# Linking occurs when the object files (.o) are put together into an executable (.exe)

# Output program (LINK_TARGET)
PROJECT=IceBox.exe
# Compiler
CC= g++
OBJECTS= main.o 
# Compiler flags (define CPLUSPLUS for xdrfile preprocessor commands and such)
CFLAGS= -DCPLUSPLUS -Wall -g -std=c++11 -O3 -c

# Linking
XDR_LIB=/usr/local/lib
XDR_INCL=/usr/local/include/xdrfile
LFLAGS= -I $(XDR_INCL) -L $(XDR_LIB) -lxdrfile -lm -g


# Make all
all : $(PROJECT)

# Linking
$(PROJECT) : $(OBJECTS)
	$(CC) $(LFLAGS) -o ../$@ $^

# Compiling
%.o : %.cpp	
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f $(OBJECTS) 
	rm -f ../$(PROJECT)

