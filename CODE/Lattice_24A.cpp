#include "GenHam.h"

// Betts 18A site lattice
//----------------------------------------------------------
void GENHAM::Bonds_24A(){

  int i, tempI;
  Array<int,1> CheckSum(Nsite);

  int nX = 4;
  
  //In square indexing: PlaqX can replace Bond()  (changed in trunk 607)
  PlaqX.resize(Nsite,4);
  PlaqX = 21, 0, 3, 2,    //  Usual 4 site plaquette with spin flips on X bonds
          0, 5, 4, 3,     //  
          5, 6, 11, 4,   //   l - k
          6, 7, 12, 11,     //   |   |
          2, 3, 8, 7,     //   i - j
          3, 4, 9, 8,
          4, 11, 10, 9,
          11, 12, 17, 10,
          7, 8, 13, 12,
          8, 9, 14, 13,
          9, 10, 15, 14,
          10, 17, 16, 15,
          12, 13, 18, 17,
          13, 14, 19, 18,
          14, 15, 20, 19,
          15, 16, 21, 20,
          17, 18, 22, 16,
          18, 19, 23, 22,
          19, 20, 1, 23,
          20, 21, 2, 1,
          16, 22, 0, 21,
          22, 23, 5, 0,
          23, 1, 6, 5,
          1, 2, 7, 6;

  CheckSum=0;
  for (i=0;i<Nsite;i++){  //Checksum
    CheckSum(PlaqX(i,0))++; tempI = CheckSum(PlaqX(i,0));
    if (tempI > 4) cout<<"checksum error 0 "<<PlaqX(i,0)<<"\n";
    CheckSum(PlaqX(i,1))++; tempI = CheckSum(PlaqX(i,1));
    if (tempI > 4) cout<<"checksum error 1 "<<PlaqX(i,1)<<"\n";
    CheckSum(PlaqX(i,2))++; tempI = CheckSum(PlaqX(i,2));
    if (tempI > 4) cout<<"checksum error 1 "<<PlaqX(i,2)<<"\n";
    CheckSum(PlaqX(i,3))++; tempI = CheckSum(PlaqX(i,3));
    if (tempI > 4) cout<<"checksum error 1 "<<PlaqX(i,3)<<"\n";
  }


  int first2 = nX;
  int second2 = Nsite-nX;

  OtherTwoX.resize(Nsite,4);                       
  for (int ii=0; ii<Nsite; ii++){   //   In this case, four adjacant sites on a plaquette 
    OtherTwoX(ii,0) = PlaqX(first2,0); //   are associated with a X bond (indexed above) 
    OtherTwoX(ii,1) = PlaqX(first2,1);       //     
    OtherTwoX(ii,2) = PlaqX(second2,0);              //    1     2 
    OtherTwoX(ii,3) = PlaqX(second2,1);              //    |     |         
    first2++; if (first2 == Nsite) first2 = 0;       //   13-(0)-0     <-Bond 0    
    second2++; if (second2 == Nsite) second2 = 0;    //    |     |         
  }                                                  //    8    17


  OtherTwoY.resize(Nsite,4);
  OtherTwoY = 20, 1, 0, 3,  //    In this case, four adjacant sites on a plaquette 
              21, 2, 5, 4,   //    are associated with a Y bond (indexed above) 
              0, 3, 6, 11,   //      
              5, 4, 7, 12,   //    16 -- 1 -- 2 
              1, 6, 3, 8,    //          |
              2, 7, 4, 9,   //         (0)     <-Bond 0 
              3, 8, 11, 10,   //          |
              4, 9, 12, 17,  //    12 -- 13 -- 0
              6, 11, 8, 13,
              7, 12, 9, 14,
              8, 13, 10, 15,
              9, 14, 17, 16,
              11, 10, 13, 18,
              12, 17, 14, 19,
              13, 18, 15, 20,
              14, 19, 16, 21,
              10, 15, 18, 22,
              17, 16, 19, 23,
              18, 22, 20, 1,
              19, 23, 21, 2,
              15, 20, 22, 0,
              16, 21, 23, 5,
              22, 0, 1, 6,
              23, 5, 2, 7;
  CheckSum=0;
  for (i=0;i<Nsite;i++){  //Checksum
    CheckSum(OtherTwoY(i,0))++; tempI = CheckSum(OtherTwoY(i,0));
    if (tempI > 4) cout<<"checksum error 0 "<<OtherTwoY(i,0)<<"\n";
    CheckSum(OtherTwoY(i,1))++; tempI = CheckSum(OtherTwoY(i,1));
    if (tempI > 4) cout<<"checksum error 1 "<<OtherTwoY(i,1)<<"\n";
    CheckSum(OtherTwoY(i,2))++; tempI = CheckSum(OtherTwoY(i,2));
    if (tempI > 4) cout<<"checksum error 1 "<<OtherTwoY(i,2)<<"\n";
    CheckSum(OtherTwoY(i,3))++; tempI = CheckSum(OtherTwoY(i,3));
    if (tempI > 4) cout<<"checksum error 1 "<<OtherTwoY(i,3)<<"\n";
  }

//  for (int i=0; i<Nsite; i++) {   
//    cout<<i<<" "<<OtherTwoX(i,0)<<" ";
//    cout<<OtherTwoX(i,1)<<" ";
//    cout<<OtherTwoX(i,2)<<" ";
//    cout<<OtherTwoX(i,3)<<"\n";
//  }


}//MakeBonds

