#include "GenHam.h"

// 4x3x2 lattice with periodic boundary conditions
//----------------------------------------------------------
void GENHAM::Bonds_24(){

  Bond.resize(Nsite,4);
  //each site has either 1 or two bonds
  //
  //   1          14   8  
  //   |     or     \ /
  //  (0)           (1)
  //                
  // and an entry of -99 indicates the absence of the 2nd bond
  
  Bond = 0, 1, -99,   0,
         1, 14, 8,    0,
         2, 3, -99,   0,
         3, 8, 10,    0,
         4, 5, -99,   0,
         5, 10, 12,   0, 
         6, 7, -99,   0,
         7, 12, 14,   0,

         8, 9, -99,   0,
         9, 22, 16,   1,
         10, 11, -99, 0,
         11, 16, 18,  1,
         12, 13, -99, 0,
         13, 18, 20,  1,
         14, 15, -99, 0,
         15, 20, 22,  1,

         16, 17, -99, 0,
         17, 6, 0,    0,
         18, 19, -99, 0,
         19, 0, 2,    0,
         20, 21, -99, 0,
         21, 2, 4,    0, 
         22, 23, -99, 0,
         23, 4, 6,    0;


//the LAST entry of 1 indicates that the NONZERO 2nd bond is FM



}//MakeBonds

