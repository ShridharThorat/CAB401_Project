//
//  cvtree_parallel.hpp
//  cvtree_project
//
//  Created by Shridhar Thorat on 8/9/2023.
//

#ifndef cvtree_parallel_hpp
#define cvtree_parallel_hpp

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <fstream>
#include "/opt/homebrew/Cellar/libomp/16.0.6/include/omp.h"
#include <vector>
#include <thread>

//#include <os/signpost.h>
//#include <os/log.h>

// Macros
#define encode(ch) code[ch-'A']
#define LEN 6
#define AA_NUMBER 20
#define EPSILON 1e-010

// const
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};


#endif /* cvtree_parallel_hpp */
