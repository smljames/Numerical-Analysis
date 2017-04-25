#include "MAT.h"
#include "VEC.h"
#include <iostream>
using namespace std;

/*
1 2 3 |1  14   
5 2 4 |2  21
3 1 1 |3  8


*/

int main(void){
    MAT A(3), G(3);
    VEC x(3), y(3), b(3);

    A[0][0] = 1;
    A[0][1] = 2;
    A[0][2] = 3;
    A[1][0] = 5;
    A[1][1] = 2;
    A[1][2] = 4;
    A[2][0] = 3;
    A[2][1] = 1;
    A[2][2] = 1;

    G = luFact(A);
    b[0] = 14;
    b[1] = 21;
    b[2] = 8;

    y = fwdSubs(G, b);
    x = bckSubs(G, y);
    
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            cout << G[i][j] << " ";
        }
        cout << endl;
    }

    for(int i=0; i<3; i++){
        cout << x[i] << endl;
    }
}
