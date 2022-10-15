#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#ifndef infos_h
#define infos_h

#define OPEN_OUTPUT 0
#define CLOSE_OUTPUT 1
#define SINGLE_OUTPUT 2
#define ADD_OUTPUT 3

#define numberOfChannels 8
#define HALF_LEN_X 300 // [mm]
#define HALF_LEN_Y 22.5 // [mm]

#define ORESET       "\033[0m"
#define OBOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define OBOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define OBOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define OBOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define OBOLDWHITE   "\033[1m\033[37m"      /* Bold White */
#define OBOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */


void PrintColor(std::string message, std::string color){
  std::cout << color << message << ORESET << std::endl;
}

#endif