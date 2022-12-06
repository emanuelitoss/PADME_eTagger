#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

#ifndef infos_h
#define infos_h

#define OPEN_OUTPUT 0
#define CLOSE_OUTPUT 1
#define SINGLE_OUTPUT 2
#define ADD_OUTPUT 3

#define DX 0
#define SX 1

#define numberOfChannels 8
#define HALF_LEN_X 300 // [mm]
#define HALF_LEN_Y 22.5 // [mm]
#define max_time 10. //[ns] (range for plots)
#define min_time -2.0 //[ns] (range for plots)
#define OFFSET 8 //[ns] >> I do not have times >~ 8 ns

#define OPTION_Q_DIFFERENCE 0
#define OPTION_Q_DOS 1
#define OPTION_DELTA_POSITION 0
#define OPTION_POSITION 1

#define ORESET   "\033[0m"
#define OBLACK   "\033[30m"      /* Black */
#define ORED     "\033[31m"      /* Red */
#define OGREEN   "\033[32m"      /* Green */
#define OYELLOW  "\033[33m"      /* Yellow */
#define OBLUE    "\033[34m"      /* Blue */
#define OMAGENTA "\033[35m"      /* Magenta */
#define OCYAN    "\033[36m"      /* Cyan */
#define OWHITE   "\033[37m"      /* White */
#define OBOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define OBOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define OBOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define OBOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define OBOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define OBOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define OBOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define OBOLDWHITE   "\033[1m\033[37m"      /* Bold White */

void PrintColor(std::string message, std::string color){
  std::cout << color << message << ORESET << std::endl;
}

double_t Mean(vector <double_t> vec){
    
    double_t mean = 0;
    int num = vec.size();

    for (int i = 0; i<num; ++i) mean += vec[i];

    mean/=num;
    return mean;
}

double_t StdDeviation(vector <double_t> vec){
    
    double_t err = 0, mean = Mean(vec);
    int num = vec.size();

    for (int i = 0; i<num; ++i) err += (vec[i]*vec[i]);

    err /= (num-1);
    err -= mean*mean*num/(num-1);
    err = sqrt(err);

    return err;

}

#endif