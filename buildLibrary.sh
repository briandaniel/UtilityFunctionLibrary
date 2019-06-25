#!/bin/bash


# Root folder
PROGRAM_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"



echo "  Program will be installed in $PROGRAM_ROOT/Install"

# Create install director if it does not exist
if [ ! -d $PROGRAM_ROOT/Install ]
then

  mkdir -p $PROGRAM_ROOT/Install

  echo "  Install directory did not exist. Install directory has now been created."

else

  echo "  Install directory already exists. "
 
fi

  echo -e "\n  Running CMake configuration script: "

  cd $PROGRAM_ROOT/Build

  echo -e "  $PROGRAM_ROOT/Build/configure-utility-function-library-cpp \n \n"

  ./configure-utility-function-library-cpp

  echo -e " \n !! If the CMake script compiled successfully go to the build folder and run make. \n"
