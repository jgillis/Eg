cmake_minimum_required(VERSION 2.8)
project (HELLO NONE) 

MESSAGE(STATUS "\n\nhello world\n\n")

FIND_PROGRAM(OCTAVE_EXE octave
PATHS "C:/Octave/*/bin"
)

MESSAGE(STATUS "Octave : <${OCTAVE_EXE}>")
