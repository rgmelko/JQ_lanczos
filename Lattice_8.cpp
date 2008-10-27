#include "GenHam.h"

// Square (4x4) conventional 16 site lattice
//----------------------------------------------------------
void GENHAM::Bonds_8(){

  Bond.resize(Nsite,4);
  //each site has either 1 or two bonds
  //
  //   1           6   4  
  //   |     or     \ /
  //  (0)           (1)
  //                
  // and an entry of -99 indicates the absence of the 2nd bond
  //
  Bond = 0, 1, -99,  0,
         1, 6, 4,    0,
         2, 3, -99,  0,
         3, 4, 6,    0,
         4, 5, -99,  0,
         5, 2, 0,    1, //FM
         6, 7, -99,  0,
         7, 0, 2,    1; //FM

//the LAST entry of 1 indicates that the NONZERO 2nd bond is FM



}//MakeBonds

