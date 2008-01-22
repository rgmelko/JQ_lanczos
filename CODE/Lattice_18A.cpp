#include "GenHam.h"

// Betts 18A site lattice
//----------------------------------------------------------
void GENHAM::Bonds_18A(){

  int i, tempI;
  Array<int,1> CheckSum(Nsite);

  int nX = 3;
  
  //In square indexing: PlaqX can replace Bond()  (changed in trunk 607)
  PlaqX.resize(Nsite,4);
  PlaqX = 13, 0, 2, 1,    //  Usual 4 site plaquette with spin flips on X bonds
          0, 9, 3, 2,     //  
          9, 10, 14, 3,   //   l - k
          1, 2, 6, 5,     //   |   |
          2, 3, 7, 6,     //   i - j
          3, 14, 8, 7,
          5, 6, 11, 10,
          6, 7, 12, 11,
          7, 8, 13, 12,

          10, 11, 15, 14,
          11, 12, 16, 15,
          12, 13, 1, 16,
          14, 15, 17, 8,
          15, 16, 4, 17,
          16, 1, 5, 4,
          8, 17, 0, 13,
          17, 4, 9, 0,
          4, 5, 10, 9;

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
  OtherTwoY = 12, 16, 0, 2,  //    In this case, four adjacant sites on a plaquette 
              13, 1, 9, 3,   //    are associated with a Y bond (indexed above) 
              0, 2, 10, 14,   //      
              16, 4, 2, 6,   //    16 -- 1 -- 2 
              1, 5, 3, 7,    //          |
              2, 6, 14, 8,   //         (0)     <-Bond 0 
              4, 9, 6, 11,   //          |
              5, 10, 7, 12,  //    12 -- 13 -- 0
              6, 11, 8, 13,
              9, 3, 11, 15,
              10, 14, 12, 16,
              11, 15, 13, 1,
              3, 7, 15, 17,
              14, 8, 16, 4,
              15, 17, 1, 5,
              7, 12, 17, 0,
              8, 13, 4, 9,
              17, 0, 5, 10;

//  for (int i=0; i<Nsite; i++) {   
//    cout<<i<<" "<<OtherTwoY(i,0)<<" ";
//    cout<<OtherTwoY(i,1)<<" ";
//    cout<<OtherTwoY(i,2)<<" ";
//    cout<<OtherTwoY(i,3)<<"\n";
//  }


}//MakeBonds

