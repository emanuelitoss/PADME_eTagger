#include <iostream>
using std::cout;
using std::endl;

#ifndef libraryPrint_h
#define libraryPrint_h

// output colors
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