#include "GenHam.h"

// 3x3x2 lattice with periodic boundary conditions
//----------------------------------------------------------
void GENHAM::Bonds_18p(){

  Bond.resize(Nsite,4);
  //each site has either 1 or two bonds
  //
  //   1           6   4  
  //   |     or     \ /
  //  (0)           (1)
  //                
  // and an entry of -99 indicates the absence of the 2nd bond
  
//  Bond = 0, 1, -99,   0,
//         1, 10, 6,    0,
//         2, 3, -99,   0,
//         3, 6, 8,     0,
//         4, 5, -99,   0,
//         5, 8, 10,    0, 
//         6, 7, -99,   0,
//         7, 16, 12,   1,
//         8, 9, -99,   0,
//         9, 12, 14,   1,
//         10, 11, -99, 0,
//         11, 14, 16,  1,
//         12, 13, -99, 0,
//         13, 4, 0,    0,
//         14, 15, -99, 0,
//         15, 0, 2,    0,
//         16, 17, -99, 0,
//         17, 2, 4,    0;

//the LAST entry of 1 indicates that the NONZERO 2nd bond is FM

//Below is the same lattice with OPEN boundary conditions

  Bond = 0, 1, -99,    0,
         1, 10, 6,     0,
         2, 3, -99,    0,
         3, 6, 8,      0,
         4, 5, -99,    0,
         5, 8, 10,     0, 
         6, 7, -99,    0,
         7, 16, 12,    1,
         8, 9, -99,    0,
         9, 12, 14,    1,
         10, 11, -99,  0,
         11, 14, 16,   1,
         12, 13, -99,  0,
         13, -99, -99, 0, //OBC
         14, 15, -99,  0,
         15, -99, -99, 0,
         16, 17, -99,  0,
         17, -99, -99, 0;


}//MakeBonds

