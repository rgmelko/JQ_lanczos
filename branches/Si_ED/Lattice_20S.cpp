#include "GenHam.h"

// Betts 18A site lattice
//----------------------------------------------------------
void GENHAM::Bonds_20S(){

  int i, tempI;
  Array<int,1> CheckSum(Nsite);

  int nX = 4;
  
  //In square indexing: PlaqX can replace Bond()  (changed in trunk 607)
  PlaqX.resize(Nsite,4);
  PlaqX = 0, 1, 5, 4,    //  Usual 4 site plaquette with spin flips on X bonds
          1, 2, 6, 5,   
          2, 3, 7, 6,  
          3, 0, 4, 7, 
          4, 5, 9, 8,   
          5, 6, 10, 9,
          6, 7, 11, 10,
          7, 4, 8, 11,
          8, 9, 13, 12,
          9, 10, 14, 13,
          10, 11, 15, 14,
          11, 8, 12, 15,
          12, 13, 17, 16,
          13, 14, 18, 17,
          14, 15, 19, 18,
          15, 12, 16, 19,
          16, 17, 1, 0,
          17, 18, 2, 1,
          18, 19, 3, 2,
          19, 16, 0, 3;

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
  OtherTwoY = 3, 7, 1, 5,  //    In this case, four adjacant sites on a plaquette 
              0, 4, 2, 6,   //    are associated with a Y bond (indexed above) 
              1, 5, 3, 7,   //      
              2, 6, 0, 4,   //    16 -- 1 -- 2 
              7, 11, 5, 9,    //          |
              4, 8, 6, 10,   //         (0)     <-Bond 0 
              5, 9, 7, 11,   //          |
              6, 10, 4, 8,  //    12 -- 13 -- 0
              11, 15, 9, 13,
              8, 12, 10, 14,
              9, 13, 11, 15,
              10, 14, 8, 12,
              15, 19, 13, 17,
              12, 16, 14, 18,
              13, 17, 15, 19,
              14, 18, 12, 16,
              19, 3, 17, 1,
              16, 0, 18, 2,
              17, 1, 19, 3,
              18, 2, 16, 0;
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

