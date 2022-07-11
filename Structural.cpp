#include <bits/stdc++.h>
using namespace std;

// Global array of all member stiffness matrix accessible to all functions
float m_g_stiff[100][6][6];  // m_g_stiff[ith member][][] denotes a 2-D Element-stiffness Matrix
float g_stiff[100][100];     // Global Stiffness Matrix
double D[100];               // For storing the displacement
double F[100];               // For storing the forces at ith node
double Q[100];               // For storing the forces of ith member


// Function to compute each member stiffness matrix
void member_stiff_calc(float x_near, float y_near, float x_far, float y_far, int n){
    float A;                  
    float E;                
    double L;

    float x1 = x_near; float y1 = y_near;
    float x2 = x_far; float y2 = y_far;

    fflush(stdin);

    L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    cout << "L = " << L << endl;

    cout << "Enter area of cross section in sq.m: ";
    cin >> A;
    cout << "Enter E value in GPA: ";
    cin >> E;

    // // Uncommment the following line to use E in GPA
    // E = E*pow(10,9);


    float lx = (x2-x1)/L;     // cos0
    float ly = (y2-y1)/L;     // sin0
    float K = A*E/L;          // stiffness of the truss element

    cout << "lx= " << lx << endl;
    cout << "ly= " << ly << endl;
    cout << "AE/L= " << K << endl;


    // *************** ELEMENT STIFFNESS MATRIX ******************//
    // [ (cos0)^2       cos0.sin0       -(cos0)^2        -cos0.sin0]
    // [ cos0.sin0      (sin0)^2        -cos0.sin0        -(sin0)^2]
    // [ -(cos0)^2      -cos0.sin0       (cos0)^2         cos0.sin0]
    // [ -cos0.sin0     -(sin0)^2        cos0.sin0        (sin0)^2 ]
    m_g_stiff[n][1][1] = m_g_stiff[n][3][3] = lx*lx; 
    m_g_stiff[n][2][2] = m_g_stiff[n][4][4] = ly*ly;
    m_g_stiff[n][1][2] = m_g_stiff[n][2][1] = m_g_stiff[n][3][4] = m_g_stiff[n][4][3] = lx*ly;
    m_g_stiff[n][1][3] = m_g_stiff[n][3][1] = -1*m_g_stiff[n][1][1];
    m_g_stiff[n][2][4] = m_g_stiff[n][4][2] = -1*m_g_stiff[n][2][2];
    m_g_stiff[n][1][4] = m_g_stiff[n][2][3] = m_g_stiff[n][3][2] = m_g_stiff[n][4][1] = -1*m_g_stiff[n][1][2];


    // Storing lx*AE/l and ly*AE/l to compute member forces
    m_g_stiff[n][1][5] = lx*K;
    m_g_stiff[n][2][5] = ly*K;


    // Calculating Element-stiffness Matrix for the nth member 
    for(int i=1; i<=4; i++){
        for(int j=1; j<=4; j++){
            m_g_stiff[n][i][j] = m_g_stiff[n][i][j] * K;
        }
    }

}

// Function to assemble all member stiffness matrices
void assemble(int node, int no_of_member){
    int ndof = node*2;

    for(int i=0; i<ndof; i++){
        for(int j=0; j<ndof; j++){
            g_stiff[i][j] = 0;

            //Loop to search all member matrices
            for(int k=0; k<no_of_member; k++){
                for(int l=1; l<=4; l++){
                    if(m_g_stiff[k][l][0] == (i+1)){
                        for(int m=1; m<=4; m++){
                            if(m_g_stiff[k][0][m] == (j+1)){
                                g_stiff[i][j] = g_stiff[i][j] + m_g_stiff[k][l][m];
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}


// Function to solve for all the displacements and forces
void d_solve(int n_uk, int n_tot){
    // Link for the Inversion of matrix:
    // Inverse Matrix Using Gauss-Jordan / Row Reduction: https://youtu.be/cJg2AuSFdjw
    // Linear Algebra in C++ - Part 2a - Compute Matrix Inverse (Theory): https://youtu.be/wOlG_fnd3v8?t=793
    

    // since g_stiff[50][50], for Gauss-Jordan Elimination method, we need double space 
    // since we are augmenting two matrix [A | I] so matrix[100][100], 
    float matrix[100][100], inv[100][100], t, a;
    
    // creating a temp Matrix for augmentation
    int i, j, k;
    for(i=0; i<n_uk; i++){
        for(j=0; j<n_uk; j++){
            matrix[i][j] = g_stiff[i][j];
        }
    }

    // Augmenting the matrix as [A | I] to find [I | (A)^-1] using Gauss-Jordan Elimination
    for(i=0; i<n_uk; i++){
        for(j=n_uk; j<2*n_uk; j++){
            if(i==(j-n_uk))
                matrix[i][j] = 1.0;
            else
                matrix[i][j] = 0.0;
        }
    }

    // converting [A | I] to find [I | (A)^-1]
    for(i=0; i<n_uk; i++){
        t = matrix[i][i];

        for(j=i; j<2*n_uk; j++){
            matrix[i][j] = matrix[i][j]/t;
        }

        for(j=0; j<n_uk; j++){
            if(i != j){
                t = matrix[j][i];
                for(k=0; k<2*n_uk; k++){
                    matrix[j][k]=matrix[j][k]-t*matrix[i][k];
                }
            }
        }
    }


    printf("The inverse matrix is: \n");
    for(i=0; i<n_uk; i++){
        for(j=n_uk, k=0; j<2*n_uk; j++,k++){
            cout << matrix[i][j];
            inv[i][k] = matrix[i][j];
            printf("\t");
        }
        printf("\n");
    }


    for(i=0; i<n_uk; i++){
        D[i] = 0.0;
        for(int j=0; j<n_uk; j++){
            D[i] = D[i] + (inv[i][j]*F[j]);
        }
    }

    for(; i<n_tot; i++){
        F[i] = 0.0;
        for(int j=0; j<n_uk; j++){
            F[i] = F[i] + (g_stiff[i][j]*D[j]);
        }
    }
}

//Function to find each member force
double member_force(int n){
    // Qf(local) = AE/L * [-lx  -ly  lx  ly]  *  Transpose( [Dxi  Dyi  Dxj  Dyj] )
    float v1 = m_g_stiff[n][1][5];
    float v2 = m_g_stiff[n][2][5];

    double temp1[4], temp2[4];

    temp1[0] = -1*v1;
    temp1[1] = -1*v2;
    temp1[2] = v1;
    temp1[3] = v2;

    int index[4];
    index[0] = m_g_stiff[n][0][1];
    index[1] = m_g_stiff[n][0][2];
    index[2] = m_g_stiff[n][0][3];
    index[3] = m_g_stiff[n][0][4];

    temp2[0] = D[index[0]-1];
    temp2[1] = D[index[1]-1];
    temp2[2] = D[index[2]-1];
    temp2[3] = D[index[3]-1];

    // Qf(local) = AE/L * [-lx  -ly  lx  ly]  *  Transpose( [Dxi  Dyi  Dxj  Dyj] )
    double sum = 0.0;
    for(int i=0; i<4; i++){
        sum = sum + (temp1[i]*temp2[i]);
    }
    return sum;
}

void Line_Seperator(){
    cout << endl; 
    for(int i=0; i<60; i++){
        cout<<"*"; 
    }
    cout << endl;
}

int main(){
    float coord[50][2];          // for taking coordinate {x, y} of 50 nodes
    int nDOF[50][2];             // for taking the label of force in x and y direction of 50 nodes

    Line_Seperator();
    cout << "Welcome to the truss structure analyzer using the Stiffness method" << endl;
    cout << "This is a general purpose program to analyze simple trusses having upto 50 nodes and 100 members" << endl;
    Line_Seperator();

    int nnode, nmember;
    cout << "Enter the number of nodes in the structure: ";
    cin >> nnode;
    cout << "Enter the number of members in the structure: ";
    cin >> nmember;

    Line_Seperator();
    cout << "Enter the coordinates of the node in the order you would number it." << endl;
    cout << "Also enter the DOF(Degree of Freedom) number associated with that node" << endl;

    //Loop to input the coordinates
    for(int i=0; i<nnode; i++){
        cout << "NODE " << i+1 << endl;
        cout << "x-coordinate: ";
        cin >> coord[i][0];

        cout << "y-coordinate: ";
        cin >> coord[i][1];

        cout << "Label of force in X direction: ";
        cin >> nDOF[i][0];

        cout << "Label of force in Y direction: ";
        cin >> nDOF[i][1];

        cout << endl;
    }

    Line_Seperator();

    //Loop to input the member conenctivity and values of A,E,L for each member
    cout << "Now specify the connectivity between node" << endl;
    cout << "Please follow the same node numbering as before." << endl;
    cout << "Enter the node number and the required values" << endl;

    for(int i=0; i<nmember; i++){
        cout << "MEMBER " << i+1 << endl;
        
        int near, far;
        cout << "Enter near node number: ";
        cin >> near;
        cout << "Enter far node number: ";
        cin >> far;

        float XN,XF,YN,YF;
        XN = coord[near-1][0];
        XF = coord[far-1][0];
        YN = coord[near-1][1];
        YF = coord[far-1][1];

        //assigning correct DOFs for ease of assembly later
        m_g_stiff[i][0][1] = m_g_stiff[i][1][0] = nDOF[near-1][0];
        m_g_stiff[i][0][2] = m_g_stiff[i][2][0] = nDOF[near-1][1];
        m_g_stiff[i][0][3] = m_g_stiff[i][3][0] = nDOF[far-1][0];
        m_g_stiff[i][0][4] = m_g_stiff[i][4][0] = nDOF[far-1][1];

        //Function to compute the member stiffness matrix in gloabal form
        member_stiff_calc(XN, YN, XF, YF, i);

        cout<<"Member stiffness matrix is: \n";
        for(int num1=0; num1<=4; num1++){
            for(int num2=0; num2<=4; num2++){
                cout << setw(11) << m_g_stiff[i][num1][num2] << "  ";
            }
            cout << endl;
        }
        cout << endl << endl;
    }


    Line_Seperator();

    // Assembling all the Element-Stiffness matrix to obtain Global Stiffness Matrix
    assemble(nnode, nmember);

    // Printing Global Stiffness Matrix
    cout << "Global Stiffness Matrix is: " << endl;
    for(int i=0; i<nnode*2; i++){
        for(int j=0; j<nnode*2; j++){
            cout << setw(11) << g_stiff[i][j] << "   ";
        }
        cout << endl;
    }


    //Enter constraints and conditons
    int n_uk, n_k;  // uk-Unknown & k-Known
    cout << "Number of known forces/Unknown Displacements: ";
    cin >> n_uk;
    n_k = (nnode*2) - n_uk;

    cout << "Enter all the known forces at nodes(In order of DOF defined earlier)" << endl;

    int num;
    for(num=1; num<=n_uk; num++){
        cout << "F" << num << "=";
        cin >> F[num-1];
    }

    cout<<"Enter all the known displacements/constraints/support settlements at the node\n";
    for(; num<=nnode*2; num++){
        cout << "D" << num << "=";
        cin >> D[num-1];
    }


    // Function to solve for all the displacements and forces
    d_solve(n_uk, nnode*2);

    ofstream myfile;
    myfile.open("truss_result.csv");


    Line_Seperator();

    cout<<"\nALL RESULTS WILL BE SAVED IN TRUSS_RESULT.CSV(Excel File) BY THIS PROGRAM AUTOMATICALLY\n";

    Line_Seperator();


    myfile << "Nodal Force " << "," << "in N" << endl;
    cout << "All nodal forces are: " << endl;
    for(int i=0; i<nnode*2; i++){
        cout << "F" << i+1 << " = " << F[i] << endl;
        myfile << "F" << i+1 << "," << F[i] << endl;
    }

    myfile << "Nodal Displacement " << "," << "in m" << endl;
    cout << "All nodal displacements are: " << endl;
    for(int i=0; i<nnode*2; i++){
        cout << "D" << i+1 << " = " << D[i] << endl;
        myfile << "D" << i+1 << "," << " | " << D[i] << endl;
    }

    myfile << "Member forces " << "," << "in N" << endl;
    myfile << "Positive values indicate tension and negative compression" << endl;

    cout << "Now member forces are: \n";
    for(int i=0; i<nmember; i++){
            Q[i] = member_force(i);
            cout << "Q" << i+1 << " | " << Q[i];
            cout << endl;
            myfile << "Q" << i+1 << "," << Q[i] << endl;
    }

    myfile.close();
    return 0;
}