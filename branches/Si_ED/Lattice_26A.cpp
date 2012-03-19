#include "GenHam.h"

// Betts 26A site lattice
//----------------------------------------------------------
void GENHAM::Bonds_26A(){

  int i, tempI;
  Array<int,1> CheckSum(Nsite);

  //In square indexing: PlaqX can replace Bond()  (changed in trunk 607)
  PlaqX.resize(Nsite,4);
  PlaqX = 0, 19, 2, 1,    //  Usual 4 site plaquette with spin flips on X bonds
          1, 2, 6, 5,   
          2, 3, 7, 6,  
          3, 24, 8, 7, 
          4, 5, 10, 9,   
          5, 6, 11, 10,
          6, 7, 12, 11,
          7, 8, 13, 12,
          8, 0, 1, 13,
          9, 10, 16, 15,
          10, 11, 17, 16,
          11, 12, 18, 17,
          12, 13, 4, 18,
          13, 1, 5, 4,
          14, 15, 20, 19,
          15, 16, 21, 20,
          16, 17, 22, 21,
          17, 18, 23, 22,
          18, 4, 9, 23,
          19, 20, 3, 2,
          20, 21, 24, 3,
          21, 22, 25, 24,
          22, 23, 14, 25,
          23, 9, 15, 14,
          24, 25, 0, 8,
          25, 14, 19, 0;

  CheckSum=0;
  for (i=0;i<Nsite;i++){  //Checksum
    CheckSum(PlaqX(i,0))++; tempI = CheckSum(PlaqX(i,0));
    if (tempI > 4) cout<<"cHecksum error 0 "<<PlaqX(i,0)<<"\n";
    CheckSum(PlaqX(i,1))++; tempI = CheckSum(PlaqX(i,1));
    if (tempI > 4) cout<<"cHecksum error 1 "<<PlaqX(i,1)<<"\n";
    CheckSum(PlaqX(i,2))++; tempI = CheckSum(PlaqX(i,2));
    if (tempI > 4) cout<<"cHecksum error 2 "<<PlaqX(i,2)<<"\n";
    CheckSum(PlaqX(i,3))++; tempI = CheckSum(PlaqX(i,3));
    if (tempI > 4) cout<<"cHecksum error 3 "<<PlaqX(i,3)<<"\n";
  }


  OtherTwoX.resize(Nsite,4);                       

  OtherTwoX.resize(Nsite,4);
  OtherTwoX = 1, 2, 25, 14,   // 0
              5, 6, 0, 19,    // 1
              6, 7, 19, 20,   // 2       
              7, 8, 20, 21,   // 3 
              9, 10, 13, 1,   // 4
              10, 11, 1, 2,   // 5
              11, 12, 2, 3,   // 6
              12, 13, 3, 24,  // 7
              13, 1, 24, 25,  // 8
              15, 16, 4, 5,   // 9
              16, 17, 5, 6,   // 10
              17, 18, 6, 7,   // 11
              18, 4, 7, 8,    // 12
              4, 5, 8, 0,     // 13
              19, 20, 23, 9,  // 14
              20, 21, 9, 10,  // 15
              21, 22, 10, 11, // 16
              22, 23, 11, 12, // 17
              23, 9, 12, 13,  // 18
              2, 3, 14, 15,   // 19
              3, 24, 15, 16,  // 20
              24, 25, 16, 17, // 21
              25, 14, 17, 18, // 22
              14, 15, 18, 4,  // 23
              8, 0, 21, 22,   // 24
              0, 19, 22, 23;  // 25

  CheckSum=0;
  for (i=0;i<Nsite;i++){  //Checksum
    CheckSum(OtherTwoX(i,0))++; tempI = CheckSum(OtherTwoX(i,0));
    if (tempI > 4) cout<<"Checksum error 0 "<<OtherTwoX(i,0)<<"\n";
    CheckSum(OtherTwoX(i,1))++; tempI = CheckSum(OtherTwoX(i,1));
    if (tempI > 4) cout<<"Checksum error 1 "<<OtherTwoX(i,1)<<"\n";
    CheckSum(OtherTwoX(i,2))++; tempI = CheckSum(OtherTwoX(i,2));
    if (tempI > 4) cout<<"Checksum error 2 "<<OtherTwoX(i,2)<<"\n";
    CheckSum(OtherTwoX(i,3))++; tempI = CheckSum(OtherTwoX(i,3));
    if (tempI > 4) cout<<"Checksum error 3 "<<OtherTwoX(i,3)<<"\n";
  }

  OtherTwoY.resize(Nsite,4);
  OtherTwoY = 8, 13, 19, 2,    // 0
              13, 4, 2, 6,     // 1
              1, 5, 3, 7,      // 2       
              2, 6, 24, 8,     // 3   
              18, 23, 5, 10,   // 4    
              4, 9, 6, 11,     // 5  
              5, 10, 7, 12,    // 6   
              6, 11, 8, 13,    // 7
              7, 12, 0, 1,     // 8
              23, 14, 10, 16,   // 9  
              9, 15, 11, 17,   // 10 
              10, 16, 12, 18,  // 11  
              11, 17, 13, 4,   // 12 
              12, 18, 1, 5,    // 13
              25, 0, 15, 20,     // 14
              14, 19, 16, 21,   // 15 
              15, 20, 17, 22,   // 16 
              16, 21, 18, 23,    // 17
              17, 22, 4, 9,     // 18
              0, 1, 20, 3,     // 19
              19, 2, 21, 24,     // 20
              20, 3, 22, 25,     // 21
              21, 24, 23, 14,     // 22
              22, 25, 9, 15,     // 23
              3, 7, 25, 0,     // 24
              24, 8, 14, 19;    // 25

  CheckSum=0;
  for (i=0;i<Nsite;i++){  //Checksum
    CheckSum(OtherTwoY(i,0))++; tempI = CheckSum(OtherTwoY(i,0));
    if (tempI > 4) cout<<"Checksum error 0 "<<OtherTwoY(i,0)<<"\n";
    CheckSum(OtherTwoY(i,1))++; tempI = CheckSum(OtherTwoY(i,1));
    if (tempI > 4) cout<<"Checksum error 1 "<<OtherTwoY(i,1)<<"\n";
    CheckSum(OtherTwoY(i,2))++; tempI = CheckSum(OtherTwoY(i,2));
    if (tempI > 4) cout<<"Checksum error 2 "<<OtherTwoY(i,2)<<"\n";
    CheckSum(OtherTwoY(i,3))++; tempI = CheckSum(OtherTwoY(i,3));
    if (tempI > 4) cout<<"Checksum error 3 "<<OtherTwoY(i,3)<<"\n";
  }

//  for (int i=0; i<Nsite; i++) {   
//    cout<<i<<" "<<OtherTwoX(i,0)<<" ";
//    cout<<OtherTwoX(i,1)<<" ";
//    cout<<OtherTwoX(i,2)<<" ";
//    cout<<OtherTwoX(i,3)<<"\n";
//  }


}//MakeBonds

