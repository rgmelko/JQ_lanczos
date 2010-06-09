#include "GenHam.h"

// Square (4x4) conventional 16 site lattice
//----------------------------------------------------------
void GENHAM::Bonds_16B(){

  Bond.resize(Nsite,3);
  //horizontal    vertical
  Bond = 0, 1, 4,   //   This is all you need for the Heisenberg model
         1, 2, 5,   //   Each row indexes two sites associated with ll
         2, 3, 6,   //   
         3, 0, 7,   //     4
         4, 5, 8,   //     |           e.g. site 0 related to site 1 and 4
         5, 6, 9,   //   ( 0 ) - 1
         6, 7, 10,
         7, 4, 11,
         8, 9, 12,
         9, 10, 13,
         10, 11, 14,
         11, 8, 15,
         12, 13, 0,
         13, 14, 1,
         14, 15, 2,
         15, 12, 3;


  OtherTwoX.resize(Nsite,4);
  OtherTwoX = 4, 5, 12, 13, //   In this case, four adjacant sites on a plaquette 
              5, 6, 13, 14, //   are associated with a X bond (indexed above) 
              6, 7, 14, 15, //     
              7, 4, 15, 12, //    4     5 
              8, 9, 0, 1,   //    |     |
              9, 10, 1, 2,  //    0-(0)-1     <-Bond 0 
              10, 11, 2, 3, //    |     |
              11, 8, 3, 0,  //    12    13
              12, 13, 4, 5,
              13, 14, 5, 6,
              14, 15, 6, 7,
              15, 12, 7, 4,
              0, 1, 8, 9,
              1, 2, 9, 10,
              2, 3, 10, 11,
              3, 0, 11, 8;

  OtherTwoY.resize(Nsite,4);
  OtherTwoY = 3, 7, 1, 5,   //    In this case, four adjacant sites on a plaquette 
              0, 4, 2, 6,   //    are associated with a Y bond (indexed above) 
              1, 5, 3, 7,   //      
              2, 6, 0, 4,   //     7 -- 4 -- 5 
              7, 11, 5, 9,  //          |
              4, 8, 6, 10,  //         (0)     <-Bond 0 
              5, 9, 7, 11,  //          |
              6, 10, 4, 8,  //     3 -- 0 -- 1
              11, 15, 9, 13,
              8, 12, 10, 14,
              9, 13, 11, 15,
              10, 14, 8, 12,
              15, 3, 13, 1,
              12, 0, 14, 2,
              13, 1, 15, 3,
              14, 2, 12, 0;

  PlaqX.resize(Nsite,4);
  PlaqX = 0, 1, 5, 4,    //  Usual 4 site plaquette with spin flips on X bonds
          1, 2, 6, 5,    //  
          2, 3, 7, 6,    //   l - k
          3, 0, 4, 7,    //   |   |
          4, 5, 9, 8,    //   i - j
          5, 6, 10, 9,
          6, 7, 11, 10,
          7, 4, 8, 11,
          8, 9, 13, 12,
          9, 10, 14, 13,
          10, 11, 15, 14,
          11, 8, 12, 15,
          12, 13, 1, 0,
          13, 14, 2, 1,
          14, 15, 3, 2,
          15, 12, 0, 3;

  PlaqY.resize(Nsite,4);
  PlaqY = 1, 5, 4, 0,    //  Usual 4 site plaquette with spin flips on Y bonds
          2, 6, 5, 1,    //  
          3, 7, 6, 2,    //   k - j
          0, 4, 7, 3,    //   |   |
          5, 9, 8, 4,    //   l - i
          6, 10, 9, 5,
          7, 11, 10, 6,
          4, 8, 11, 7,
          9, 13, 12, 8,
          10, 14, 13, 9,
          11, 15, 14, 10,
          8, 12, 15, 11,
          13, 1, 0, 12,
          14, 2, 1, 13,
          15, 3, 2, 14,
          12, 0, 3, 15;

//  for (int i=0; i<Nsite; i++) {   
//    cout<<i<<" "<<OtherTwoY(i,0)<<" ";
//    cout<<OtherTwoY(i,1)<<" ";
//    cout<<OtherTwoY(i,2)<<" ";
//    cout<<OtherTwoY(i,3)<<"\n";
//  }


}//MakeBonds

