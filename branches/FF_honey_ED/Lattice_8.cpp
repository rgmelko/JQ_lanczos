#include "GenHam.h"

// Square (4x4) conventional 16 site lattice
//----------------------------------------------------------
void GENHAM::Bonds_8(){

  Bond.resize(Nsite,3);
  //horizontal    vertical
  Bond = 0, 1, 5,    //   This is all you need for the Heisenberg model
         1, 0, 1,   //   Each row indexes two sites associated with ll
         2, 3, 6,   //   
         3, 0, 7,   //       1
         4, 5, 8,   //       |           e.g. site 0 related to site 1 and 4
         5, 6, 9,   //     ( 0 ) 
         6, 7, 10,  
         7, 4, 11,  //     5   7
         8, 4, 11,  //
         9, 4, 11,
         10, 4, 11,
         11, 4, 11;



}//MakeBonds

