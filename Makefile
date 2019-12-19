#############################################################################
# Makefile for building: Jeny_openMP
# Generated by qmake (3.1) (Qt 5.9.7)
# Project:  Jeny_openMP.pro
# Template: app

# #############################################################################

# MAKEFILE      = Makefile

# ####### Compiler, tools and options

# CC            = gcc
# CXX           = g++
# DEFINES       = -DQT_DEPRECATED_WARNINGS -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB
# CFLAGS        = -pipe -O2 -Wall -W -D_REENTRANT -fPIC $(DEFINES)

# CFLAGS += $(shell mpicc -showme:compile)
# LDFLAGS += $(shell mpicc -showme:link)

all:
	mpiCC main.cpp utils.h transpose_openmp.h transpose_MPI.h -o main -fopenmp -O2 -Wall -std=c++11

only_main: 
	mpiCC main.cpp -o only_main -fopenmp -O2 -Wall -std=c++11
