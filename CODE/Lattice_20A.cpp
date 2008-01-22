#include "GenHam.h"

// Betts 20A site lattice
//----------------------------------------------------------
void GENHAM::Bonds_20A(){

  int i, tempI;
  Array<int,1> CheckSum(Nsite);

//  int nX = 5;
  
  //In square indexing: PlaqX can replace Bond()  (changed in trunk 607)
  PlaqX.resize(Nsite,4);
  PlaqX = 0, 14, 2, 1,    //  Usual 4 site plaquette with spin flips on X bonds
          1, 2, 6, 5,   
          2, 3, 7, 6,  
          3, 18, 8, 7, 
          4, 5, 10, 9,   
          5, 6, 11, 10,
          6, 7, 12, 11,
          7, 8, 13, 12,
          8, 0, 1, 13,
          9, 10, 15, 14,
          10, 11, 16, 15,
          11, 12, 17, 16,
          12, 13, 4, 17,
          13, 1, 5, 4,
          14, 15, 3, 2,
          15, 16, 18, 3,
          16, 17, 19, 18,
          17, 4, 9, 19,
          18, 19, 0, 8,
          19, 9, 14, 0;

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
  OtherTwoX = 1, 2, 19, 9,  //    In this case, four adjacant sites on a plaquette 
              5, 6, 0, 14,   //    are associated with a Y bond (indexed above) 
              6, 7, 14, 15,        
              7, 8, 15, 16,   
              9, 10, 13, 1,  
              10, 11, 1, 2,  
              11, 12, 2, 3,  
              12, 13, 3, 18,
              13, 1, 18, 19,
              14, 15, 4, 5,
              15, 16, 5, 6,
              16, 17, 6, 7,
              17, 4, 7, 8,
              4, 5, 8, 0,
              2, 3, 9, 10,
              3, 18, 10, 11,
              18, 19, 11, 12,
              19, 9, 12, 13,
              8, 0, 16, 17,
              0, 14, 17, 4;

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
  OtherTwoY = 8, 13, 14, 2,  //    In this case, four adjacant sites on a plaquette 
              13, 4, 2, 6,   //    are associated with a Y bond (indexed above) 
              1, 5, 3, 7,        
              2, 6, 18, 8,   
              17, 19, 5, 10,  
              4, 9, 6, 11,  
              5, 10, 7, 12,  
              6, 11, 8, 13,
              7, 12, 0, 1,
              19, 0, 10, 15,
              9, 14, 11, 16,
              10, 15, 12, 17,
              11, 16, 13, 4,
              12, 17, 1, 5,
              0, 1, 15, 3,
              14, 2, 16, 18,
              15, 3, 17, 19,
              16, 18, 4, 9,
              3, 7, 19, 0,
              18, 8, 9, 14;

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

