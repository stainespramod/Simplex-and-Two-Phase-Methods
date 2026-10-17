#include <iostream>
#include <vector>
#include <algorithm>
#include<fstream>
using namespace std;

// note: in certain cases, you may not get the right answer because of precision errors with long double

/*

format of input file:
Enter the total number of variables in the problem(including slack variables).
Enter the coefficients of each of the variables of the objective function to be maximised.
Enter the number of constraint equations m.
For each of the m equations, follow these steps:
1.Enter the coefficients of each of the variables.
2.Enter the value of the constant.

*/

void print_tableau(ofstream &outputFile,vector<vector<long double>> &A,vector<long double> &xB, vector<long double> &yB, vector<long double> &redcost, long double z_max, int n, vector<long double> &z){
    outputFile<<"Final tableau:\n";
    int i=0;
    outputFile<<"cB        yB        xB       ";
    for(int i=1;i<=n;i++) outputFile<<"y"<<i<<"       ";
    outputFile<<"\n";
    for(vector<long double> a:A){
        outputFile<<fixed<<setprecision(2)<<z[yB[i]-1]<<"      "<<"y"<<fixed<<setprecision(0)<<yB[i]<<"     "<<fixed<<setprecision(2)<<xB[i]<<"     ";
        for(long double x:a) outputFile<<fixed<<setprecision(2)<<x<<"     ";
        outputFile<<"\n";
        i++;
    }
    outputFile<<"          z      "<<fixed<<setprecision(2)<<z_max<<"     ";
    for(long double x:redcost) outputFile<<fixed<<setprecision(2)<<x<<"     ";
    outputFile<<"\n";
}

void alternative_optima(ofstream &outputFile,vector<vector<long double>> &A, vector<long double> &xB,vector<long double> &redcost, int m, int n, int z_max, vector<long double> &yB, int pivot_column){
    vector<long double> minratio_test=xB,sol1_yB=yB,sol1_xB=xB;
    long double min_ratio=INT_MAX,min_redcost=*min_element(redcost.begin(),redcost.end());
    int pivot_row;
    for(int i=0;i<m;i++){
        if(A[i][pivot_column]>0&&minratio_test[i]/A[i][pivot_column]<min_ratio){
            min_ratio=minratio_test[i]/A[i][pivot_column];
            pivot_row=i;
        }
    }
    long double pivot_element=A[pivot_row][pivot_column];
    for(int i=0;i<n;i++) A[pivot_row][i]/=pivot_element;
    xB[pivot_row]/=pivot_element;
    pivot_element=A[pivot_row][pivot_column];
    for(int i=0;i<m;i++){
        if(i==pivot_row) continue;
        long double factor=A[i][pivot_column]/pivot_element;
        for(int j=0;j<n;j++) A[i][j]-=A[pivot_row][j]*factor;
        xB[i]-=xB[pivot_row]*factor;
    }
    long double factor=redcost[pivot_column]/pivot_element;
    for(int j=0;j<n;j++) redcost[j]-=A[pivot_row][j]*factor;
    z_max-=xB[pivot_row]*factor;
    yB[pivot_row]=pivot_column+1;
    outputFile<<"Maximum achievable value = "<<z_max<<"\n";
    outputFile<<"The basic variables can assume any values on the line segment between:\n";
    for(int i=0;i<m;i++) outputFile<<"y"<<sol1_yB[i]<<" : "<<sol1_xB[i]<<" ";
    outputFile<<"and\n";
    for(int i=0;i<m;i++) outputFile<<"y"<<yB[i]<<" : "<<xB[i]<<" ";
    outputFile<<"\n\n";
}

void unbounded_case(ofstream &outputFile){
    outputFile<<"This is a case of unbounded solution\n\n";
}

void base_case(ofstream &outputFile, vector<vector<long double>> &A, vector<long double> &xB,vector<long double> &redcost, int m, int n, int z_max, vector<long double> &yB, vector<long double> &z){
    vector<long double> minratio_test=xB;
    long double min_ratio=INT_MAX,pivot_row,pivot_column,min_redcost=*min_element(redcost.begin(),redcost.end());
    if(min_redcost==0){
        for(int i=0;i<n;i++){
            if(redcost[i]==0){
                int flag=0;
                for(int j=0;j<m;j++){
                    if(i==yB[j]-1){
                        flag=1;
                        break;
                    }
                }
                if(flag) continue;
                alternative_optima(outputFile,A,xB,redcost,m,n,z_max,yB,i);
                print_tableau(outputFile,A,xB,yB,redcost,z_max,n,z);
                return;
            }
        }
    }
    if(min_redcost>=0){
        outputFile<<"Maximum achievable value = "<<z_max<<"\n";
        outputFile<<"The values of the basic variables are:\n";
        for(int i=0;i<m;i++) outputFile<<"y"<<yB[i]<<" : "<<xB[i]<<" ";
        outputFile<<"\n\n";
        print_tableau(outputFile,A,xB,yB,redcost,z_max,n,z);
        return;
    }
    for(int i=0;i<n;i++){
        if(redcost[i]==min_redcost){
            pivot_column=i;
            int flag=0;
            for(int j=0;j<m;j++){
                if(A[j][pivot_column]>0){
                    flag=1;
                    break;
                }
                else continue;
            }
            if(!flag){
                unbounded_case(outputFile);
                print_tableau(outputFile,A,xB,yB,redcost,z_max,n,z);
                return;
            }
            break;
        }
    }
    for(int i=0;i<m;i++){
        if(A[i][pivot_column]>0&&minratio_test[i]/A[i][pivot_column]<min_ratio){
            min_ratio=minratio_test[i]/A[i][pivot_column];
            pivot_row=i;
        }
    }
    long double pivot_element=A[pivot_row][pivot_column];
    for(int i=0;i<n;i++) A[pivot_row][i]/=pivot_element;
    xB[pivot_row]/=pivot_element;
    pivot_element=A[pivot_row][pivot_column];
    for(int i=0;i<m;i++){
        if(i==pivot_row) continue;
        long double factor=A[i][pivot_column]/pivot_element;
        for(int j=0;j<n;j++) A[i][j]-=A[pivot_row][j]*factor;
        xB[i]-=xB[pivot_row]*factor;
    }
    long double factor=redcost[pivot_column]/pivot_element;
    for(int j=0;j<n;j++) redcost[j]-=A[pivot_row][j]*factor;
    z_max-=xB[pivot_row]*factor;
    yB[pivot_row]=pivot_column+1;
    base_case(outputFile,A,xB,redcost,m,n,z_max,yB,z);
}

void matrix_build(ifstream &inputFile,vector<vector<long double>> &A, vector<long double> &b, int m, int n_max){
    for(int i=0;i<m;i++){
        for(int j=0;j<n_max;j++){
            long double x;
            inputFile>>x;
            A[i].push_back(x);
        }
        long double x;
        inputFile>>x;
        b.push_back(x);
    }
}

int main(){
    ifstream inputFile("input_simplex3.txt");
    ofstream outputFile("output_simplex3.txt");
    int m,n_max;
    inputFile>>n_max;
    vector<long double> z(n_max);
    for(long double &x:z) inputFile>>x;
    inputFile>>m;
    vector<vector<long double>> A(m);
    vector<long double> b,xB,yB,redcost(n_max);
    for(int i=1;i<=m;i++) yB.push_back(n_max-m+i);
    matrix_build(inputFile,A,b,m,n_max);
    xB=b;
    for(int i=0;i<n_max;i++) redcost[i]=-z[i];
    base_case(outputFile,A,xB,redcost,m,n_max,0,yB,z);
    inputFile.close();
    outputFile.close();
}