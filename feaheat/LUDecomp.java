package feaheat;
import java.lang.*;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.lang.Math;

public class LUDecomp {
  double LUDecomp(int N, double M[][], int mi, int mj, int side) {
    int n = ( N + 1 ) * ( N + 1 );
    int count = N + 2;
    double L[][] = new double[n][n]; //Lower matrix
    double U[][] = new double[n][n]; //Upper matrix
    double O[][] = new double[n][n]; //Output matrix

    for (int j = 0; j < count; ++j) {
      U[0][j] = M[0][j];
    }
    for (int i = 0; i < count; ++i) {
      L[i][0] = M[i][0] / U[0][0];
    }
    for (int i = 0; i < n; ++i) {
      L[i][i] = 1;
    }
    for (int i = 1; i < n; ++i) {
      if ( ( count + i ) < n ) {
        for (int j = i; j < ( count + i ); ++j) {
          double U_sum = 0;
          for (int k = 1; k < i + 1; ++k) {
            U_sum = U_sum + ( L[i][k - 1] * U[k - 1][j] );
          }
          U[i][j] = M[i][j] - U_sum;
        }
        for (int j = i + 1; j < ( count + i ); ++j) {
          double L_sum = 0;
          for (int k = 1; k < i + 1; ++k) {
            L_sum = L_sum + ( L[j][k - 1] * U[k - 1][i] );
          }
          L[j][i] = ( M[j][i] - L_sum ) / U[i][i];
        }
      }
      else {
        for (int j = i; j < n; ++j) {
          double U_sum = 0;
          for (int k = 1; k < i + 1; ++k) {
            U_sum = U_sum + ( L[i][k - 1] * U[k - 1][j] );
          }
          U[i][j] = M[i][j] - U_sum;
        }
        for (int j = i + 1; j < n; ++j) {
          double L_sum = 0;
          for (int k = 1; k < i + 1; ++k) {
            L_sum = L_sum + ( L[j][k - 1] * U[k - 1][i] );
          }
          L[j][i] = ( M[j][i] - L_sum ) / U[i][i];
        }
      }
    }

    switch (side) {
      case 1:
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
            O[i][j] = L[i][j];
          }
        }
        break;
      case 2:
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
            O[i][j] = U[i][j];
          }
        }
        break;
      default:
        break;
    }

    return O[mi][mj];
  }
}
