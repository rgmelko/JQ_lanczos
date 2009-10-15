#include "GenHam.h"

// 30 site rectangular (5x6) lattice
//----------------------------------------------------------
void GENHAM::Bonds_30R(){

  int i, tempI;
  Array<int,1> CheckSum(Nsite);

  int nX = 5;
  
  //In square indexing: PlaqX can replace Bond()  (changed in trunk 607)
  PlaqX.resize(Nsite,4);
  PlaqX = 0, 1, 6, 5,   
          1, 2, 7, 6,   
          2, 3, 8, 7,   
          3, 4, 9, 8,   
          4, 0, 5, 9,   

          5, 6, 11, 10,
          6, 7, 12, 11,
          7, 8, 13, 12,
          8, 9, 14, 13,
          9, 5, 10, 14,

          10, 11, 16, 15,
          11, 12, 17, 16,
          12, 13, 18, 17,
          13, 14, 19, 18,
          14, 10, 15, 19,

          15, 16, 21, 20,
          16, 17, 22, 21,
          17, 18, 23, 22,
          18, 19, 24, 23,
          19, 15, 20, 24,

          20, 21, 26, 25,
          21, 22, 27, 26,
          22, 23, 28, 27,
          23, 24, 29, 28,
          24, 20, 25, 29,

          25, 26, 1, 0,
          26, 27, 2, 1,
          27, 28, 3, 2,
          28, 29, 4, 3,
          29, 25, 0, 4;

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
  OtherTwoY = 4, 9, 1, 6,  //    In this case, four adjacant sites on a plaquette 
              0, 5, 2, 7,   //    are associated with a Y bond (indexed above) 
              1, 6, 3, 8, 
              2, 7, 4, 9, 
              3, 8, 0, 5, 
              
              9, 14, 6, 11, 
              5, 10, 7, 12, 
              6, 11, 8, 13, 
              7, 12, 9, 14, 
              8, 13, 5, 10, 

              14, 19, 11, 16, 
              10, 15, 12, 17, 
              11, 16, 13, 18, 
              12, 17, 14, 19, 
              13, 18, 10, 15, 

              19, 24, 16, 21, 
              15, 20, 17, 22, 
              16, 21, 18, 23, 
              17, 22, 19, 24, 
              18, 23, 15, 20, 

              24, 29, 21, 26, 
              20, 25, 22, 27, 
              21, 26, 23, 28, 
              22, 27, 24, 29, 
              23, 28, 20, 25, 

              29, 4, 26, 1, 
              25, 0, 27, 2, 
              26, 1, 28, 3, 
              27, 2, 29, 4, 
              28, 3, 25, 0; 
             
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

