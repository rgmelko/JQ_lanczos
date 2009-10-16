#include "GenHam.h"

// Betts 22A site lattice
//----------------------------------------------------------
void GENHAM::Bonds_22A(){

  int i, tempI;
  Array<int,1> CheckSum(Nsite);

  //In square indexing: PlaqX can replace Bond()  (changed in trunk 607)
  PlaqX.resize(Nsite,4);
  PlaqX = 0, 1, 6, 5,
          1, 2, 7, 6,   
          2, 3, 8, 7,  
          3, 4, 9, 8, 
          4, 5, 10, 9,   
          5, 6, 11, 10,
          6, 7, 12, 11,
          7, 8, 13, 12,
          8, 9, 14, 13,
          9, 10, 15, 14,
          10, 11, 16, 15,
          11, 12, 17, 16,
          12, 13, 18, 17,
          13, 14, 19, 18,
          14, 15, 20, 19,
          15, 16, 21, 20,
          16, 17, 0, 21,
          17, 18, 1, 0,
          18, 19, 2, 1,
          19, 20, 3, 2,
          20, 21, 4, 3,
          21, 0, 5, 4;

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
  //cout<<"Checksum 1"<<endl;


  OtherTwoX.resize(Nsite,4);                       

  OtherTwoX.resize(Nsite,4);
  OtherTwoX = 5, 6, 17, 18,   // 0
              6, 7, 18, 19,    // 1
              7, 8, 19, 20,   // 2       
              8, 9, 20, 21,   // 3 
              9, 10, 21, 0,   // 4
              10, 11, 0, 1,   // 5
              11, 12, 1, 2,   // 6
              12, 13, 2, 3,  // 7
              13, 14, 3, 4,  // 8
              14, 15, 4, 5,   // 9
              15, 16, 5, 6,   // 10
              16, 17, 6, 7,   // 11
              17, 18, 7, 8,    // 12
              18, 19, 8, 9,     // 13
              19, 20, 9, 10,   // 14
              20, 21, 10, 11,  // 15
              21, 0, 11, 12, // 16
              0, 1, 12, 13, // 17
              1, 2, 13, 14,  // 18
              2, 3, 14, 15,   // 19
              3, 4, 15, 16,   // 20
              4, 5, 16, 17;   // 21

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
  //cout<<"Checksum 2"<<endl;

  OtherTwoY.resize(Nsite,4);
  OtherTwoY = 21, 4, 1, 6,    // 0
              0, 5, 2, 7,     // 1
              1, 6, 3, 8,      // 2       
              2, 7, 4, 9,     // 3   
              3, 8, 5, 10,   // 4    
              4, 9, 6, 11,     // 5  
              5, 10, 7, 12,    // 6   
              6, 11, 8, 13,    // 7
              7, 12, 9, 14,     // 8
              8, 13, 10, 15,   // 9  
              9, 14, 11, 16,   // 10 
              10, 15, 12, 17,  // 11  
              11, 16, 13, 18,   // 12 
              12, 17, 14, 19,    // 13
              13, 18, 15, 20,     // 14
              14, 19, 16, 21,   // 15 
              15, 20, 17, 0,   // 16 
              16, 21, 18, 1,    // 17
              17, 0, 19, 2,     // 18
              18, 1, 20, 3,     // 19
              19, 2, 21, 4,     // 20
              20, 3, 0, 5;     // 21

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
  //cout<<"Checksum 3"<<endl;

//  for (int i=0; i<Nsite; i++) {   
//    cout<<i<<" "<<OtherTwoX(i,0)<<" ";
//    cout<<OtherTwoX(i,1)<<" ";
//    cout<<OtherTwoX(i,2)<<" ";
//    cout<<OtherTwoX(i,3)<<"\n";
//  }


}//MakeBonds

