#include "GenHam.h"

// Mishmash's 24M lattice to test DMRG
//----------------------------------------------------------
void GENHAM::Bonds_24M(){

  int i, tempI;
  Array<int,1> CheckSum(Nsite);

  int nX = 3;
  
  //In square indexing: PlaqX can replace Bond()  (changed in trunk 607)
  PlaqX.resize(Nsite,4);
  PlaqX = 0, 1, 4, 3,    //  Usual 4 site plaquette with spin flips on X bonds
          1, 2, 5, 4,     //  
          2, 0, 3, 5,     //   l - k
          3, 4, 7, 6,     //   |   |
          4, 5, 8, 7,     //   i - j
          5, 3, 6, 8,
          6, 7, 10, 9,
          7, 8, 11, 10,
          8, 6, 9, 11,
          9, 10, 13, 12,
          10, 11, 14, 13,
          11, 9, 12, 14,
          12, 13, 16, 15,
          13, 14, 17, 16,
          14, 12, 15, 17,
          15, 16, 19, 18,
          16, 17, 20, 19,
          17, 15, 18, 20,
          18, 19, 22, 21,
          19, 20, 23, 22,
          20, 18, 21, 23,
          21, 22, 1, 0,
          22, 23, 2, 1,
          23, 21, 0, 2;

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
  OtherTwoY = 2, 5, 1, 4,  //    In this case, four adjacant sites on a plaquette 
              0, 3, 2, 5,   //    are associated with a Y bond (indexed above) 
              1, 4, 0, 3,   //      

              5, 8, 4, 7,   //    16 -- 1 -- 2 
              3, 6, 5, 8, 
              4, 7, 3, 6, 

              8, 11, 7, 10, 
              6, 9, 8, 11, 
              7, 10, 6, 9,

              11, 14, 10, 13,
              9, 12, 11, 14, 
              10, 13, 9, 12, 

              14, 17, 13, 16,
              12, 15, 14, 17,
              13, 16, 12, 15, 

              17, 20, 16, 19, 
              15, 18, 17, 20,
              16, 19, 15, 18,

              20, 23, 19, 22, 
              18, 21, 20, 23, 
              19, 22, 18, 21,

              23, 2, 22, 1,
              21, 0, 23, 2, 
              22, 1, 21, 0; 

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
//    cout<<i<<" "<<OtherTwoY(i,0)<<" ";
//    cout<<OtherTwoY(i,1)<<" ";
//    cout<<OtherTwoY(i,2)<<" ";
//    cout<<OtherTwoY(i,3)<<"\n";
//  }


}//MakeBonds

