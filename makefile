#make file

all:    #target name
	mpicc wsn.c utils.c alert.c balloon.c node.c -lm -lpthread -o wsn