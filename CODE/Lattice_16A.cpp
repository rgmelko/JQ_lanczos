#include "GenHam.h"

// tilted Betts 16 site lattice
//----------------------------------------------------------
void GENHAM::Bonds_16A(){

  int i, tempI;
  Array<int,1> CheckSum(Nsite);

  Bond.resize(Nsite,3);
  //horizontal    vertical
  Bond =0, 1, 7,  // 0     //   This is all you need for the Heisenberg model
        1, 2, 4,  // 1     //   Each row indexes two sites associated with ll
        2, 3, 5,  // 2     //   
        3, 0, 6,  // 3     //     7
        7, 4, 11, // 4     //     |           e.g. site 0 related to site 1 and 4
        4, 5, 8,  // 5     //   ( 0 ) - 1
        5, 6, 9,  // 6
        6, 7, 10, // 7
        11, 8, 14, // 8
        8, 9, 15, // 9
        9, 10, 12,// 10
        10, 11, 13,// 11
        14, 15, 2, // 12
        15, 12, 3, // 13
        12, 13, 0, // 14
        13, 14, 1; // 15

//  CheckSum=0;
//  for (i=0;i<Nsite;i++){  //Checksum
//    CheckSum(Bond(i,0))++; tempI = CheckSum(Bond(i,0));
//    if (tempI > 2) cout<<"checksum error 0 "<<Bond(i,0)<<"\n";
//    CheckSum(Bond(i,1))++; tempI = CheckSum(Bond(i,1));
//    if (tempI > 2) cout<<"checksum error 1 "<<Bond(i,1)<<"\n";
//  }

  OtherTwoX.resize(Nsite,4);
  OtherTwoX = 7, 4, 12, 13, //   In this case, four adjacant sites on a plaquette 
              4, 5, 13, 14, //   are associated with a X bond (indexed above) 
              5, 6, 14, 15, //     
              6, 7, 15, 12, //    7     4 
              11, 8, 0, 1,  //    |     |
              8, 9, 1, 2,   //    0-(0)-1     <-Bond 0 
              9, 10, 2, 3,  //    |     |
              10, 11, 3, 0, //    12    13
              14, 15, 7, 4,
              15, 12, 4, 5,
              12, 13, 5, 6,
              13, 14, 6, 7,
              2, 3, 11, 8,
              3, 0, 8, 9,
              0, 1, 9, 10,
              1, 2, 10, 11;

//  CheckSum=0;
//  for (i=0;i<Nsite;i++){  //Checksum
//    CheckSum(OtherTwoX(i,0))++; tempI = CheckSum(OtherTwoX(i,0));
//    if (tempI > 4) cout<<"checksum error 0 "<<OtherTwoX(i,0)<<"\n";
//    CheckSum(OtherTwoX(i,1))++; tempI = CheckSum(OtherTwoX(i,1));
//    if (tempI > 4) cout<<"checksum error 1 "<<OtherTwoX(i,1)<<"\n";
//    CheckSum(OtherTwoX(i,2))++; tempI = CheckSum(OtherTwoX(i,2));
//    if (tempI > 4) cout<<"checksum error 2 "<<OtherTwoX(i,2)<<"\n";
//    CheckSum(OtherTwoX(i,3))++; tempI = CheckSum(OtherTwoX(i,3));
//    if (tempI > 4) cout<<"checksum error 3 "<<OtherTwoX(i,3)<<"\n";
//  }

  OtherTwoY.resize(Nsite,4);
  OtherTwoY = 3, 6, 1, 4,   //    In this case, four adjacant sites on a plaquette 
              0, 7, 2, 5,   //    are associated with a Y bond (indexed above) 
              1, 4, 3, 6,   //      
              2, 5, 0, 7,   //     6 -- 7 -- 4 
              6, 10, 4, 8,  //          |
              7, 11, 5, 9,  //         (0)     <-Bond 0 
              4, 8, 6, 10,  //          |
              5, 9, 7, 11,  //     3 -- 0 -- 1
              10, 13, 8, 15,
              11, 14, 9, 12,
              8, 15, 10, 13,
              9, 12, 11, 14,
              13, 1, 15, 3,
              14, 2, 12, 0,
              15, 3, 13, 1,
              12, 0, 14, 2;

  PlaqX.resize(Nsite,4);
  PlaqX = 0, 1, 4, 7,    //  Usual 4 site plaquette with spin flips on X bonds
          1, 2, 5, 4,    //  
          2, 3, 6, 5,    //   l - k
          3, 0, 7, 6,    //   |   |
          7, 4, 8, 11,   //   i - j
          4, 5, 9, 8,
          5, 6, 10, 9,
          6, 7, 11, 10,
          11, 8, 15, 14,
          8, 9, 12, 15,
          9, 10, 13, 12,
          10, 11, 14, 13,
          14, 15, 3, 2,
          15, 12, 0, 3,
          12, 13, 1, 0,
          13, 14, 2, 1;

  PlaqY.resize(Nsite,4);
  PlaqY = 1, 4, 7, 0,    //  Usual 4 site plaquette with spin flips on Y bonds
          2, 5, 4, 1,    //  
          3, 6, 5, 2,    //   k - j
          0, 7, 6, 3,    //   |   |
          4, 8, 11, 7,   //   l - i
          5, 9, 8, 4,
          6, 10, 9, 5,
          7, 11, 10, 6,
          8, 15, 14, 11,
          9, 12, 15, 8,
          10, 13, 12, 9,
          11, 14, 13, 10,
          15, 3, 2, 14,
          12, 0, 3, 15,
          13, 1, 0, 12,
          14, 2, 1, 13;

//  for (int i=0; i<Nsite; i++) {   
//    cout<<i<<" "<<OtherTwoY(i,0)<<" ";
//    cout<<OtherTwoY(i,1)<<" ";
//    cout<<OtherTwoY(i,2)<<" ";
//    cout<<OtherTwoY(i,3)<<"\n";
//  }


}//MakeBonds

