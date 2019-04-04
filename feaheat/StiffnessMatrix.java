package feaheat;
import java.lang.*;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.lang.Math;

public class StiffnessMatrix {
  double StiffnessMatrix(int N, double t, double x_size, double y_size, int numEle, double k_xx, double k_yy, int K_list[], double K_val[][], int mi, int mj) {

    //--------------------------------------------------------------------------
    //Definitions
    int E_l = 2 * N * N, w = 2 * N;
    int n = ( N + 1 ) * ( N + 1 );
    int sqn = N + 1;
    double gamma1 = (double) x_size / ( 2 * y_size );
    double gamma2 = (double) y_size / ( 2 * x_size );
    double y_sizerec = (double) 1 / y_size;
    double kxb = 1 + ( 0.5 * ( ( gamma2 * k_xx ) + ( gamma1 * k_yy ) - 2 ) );
    double KK[][] = new double[E_l][2];
    double W[][][] = new double[sqn][sqn][3];
    double M[][] = new double[n][n];

    //--------------------------------------------------------------------------
    //K matrix
    for (int i = 0; i < E_l; ++i) {
      KK[i][0] = k_xx;
      KK[i][1] = k_yy;
    }

    for (int i = 0; i < numEle; ++i) {
      KK[K_list[i]][0] = K_val[i][0];
      KK[K_list[i]][1] = K_val[i][1];
    }

    //Matrix generation
    W[0][0][0] = 2 * kxb * t;
    W[N][N][0] = 2 * kxb * t;
    W[0][0][1] = 4 * kxb * t;
    W[N][N][1] = 4 * kxb * t;
    W[0][0][2] = - 0.5 * t * x_size * y_sizerec * k_yy;
    W[N][N][2] = - 0.5 * t * x_size * y_sizerec * k_yy;

    if (N >= 2) {
      for (int i = 1; i < N; ++i) {
        W[i][i][0] = 4 * kxb * t;
        W[i][i][1] = 8 * kxb * t;
        W[i][i][2] = - 1 * t * x_size * y_sizerec * k_yy;
      }
    }

    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < 2; ++j) {
        W[i][i + 1][j] = - ( W[0][0][j] + W[j][j][2] );
        W[i + 1][i][j] = - ( W[0][0][j] + W[j][j][2] );
      }
    }

    for (int k = sqn; k < n - sqn; k += sqn) {
      for (int i = 0; i < sqn; ++i) {
        for (int j = 0; j < sqn; ++j) {
          M[i + k][j + k] = W[i][j][1];
        }
      }
    }

    for (int i = 0; i < sqn; ++i) {
      for (int j = 0; j < sqn; ++j) {
        M[i][j] = W[i][j][0];
        M[i + n - sqn][j + n - sqn] = W[i][j][0];
      }
    }

    for (int k = 0; k < n - sqn; k += sqn) {
      for (int i = 0; i < sqn; ++i) {
        for (int j = 0; j < sqn; ++j) {
          M[i + k][sqn + j + k] = W[i][j][2];
          M[sqn + i + k][j + k] = W[i][j][2];
        }
      }
    }

    int q1 = 0, q5 = 0;
    for (int i = 0; i < N; ++i) {
      int q6 = 0;
      for (int j = 0; j < w; ++j) {
        int q2 = q1 + sqn, q3 = q1 + N + 2, q4 = q1 - 1;
        if ( getParity( j + q6 + 1 ) == true ) {
          //Odd parity
          M[q1 + q5][q1 + q5] = M[q1 + q5][q1 + q5] + ( 0.5 * t * KK[j + ( w * i )][1] ) - 0.5;
          M[q2 + q5][q2 + q5] = M[q2 + q5][q2 + q5] + ( 0.5 * t * KK[j + ( w * i )][1] ) - 0.5;
          M[q1 + q5][q2 + q5] = M[q1 + q5][q2 + q5] - ( 0.5 * t * KK[j + ( w * i )][1] ) + 0.5;
          M[q2 + q5][q1 + q5] = M[q2 + q5][q1 + q5] - ( 0.5 * t * KK[j + ( w * i )][1] ) + 0.5;
          M[q2 + q5][q2 + q5] = M[q2 + q5][q2 + q5] + ( 0.5 * t * KK[j + ( w * i )][0] ) - 0.5;
          M[q3 + q5][q3 + q5] = M[q3 + q5][q3 + q5] + ( 0.5 * t * KK[j + ( w * i )][0] ) - 0.5;
          M[q2 + q5][q3 + q5] = M[q2 + q5][q3 + q5] - ( 0.5 * t * KK[j + ( w * i )][0] ) + 0.5;
          M[q3 + q5][q2 + q5] = M[q3 + q5][q2 + q5] - ( 0.5 * t * KK[j + ( w * i )][0] ) + 0.5;

          q1 ++;
          q6 ++;
        }
        else if ( getParity( j + q6 + 1 ) == false ) {
          //Even parity
          M[q1 + q5][q1 + q5] = M[q1 + q5][q1 + q5] + ( 0.5 * t * KK[j + ( w * i )][1] ) - 0.5;
          M[q2 + q5][q2 + q5] = M[q2 + q5][q2 + q5] + ( 0.5 * t * KK[j + ( w * i )][1] ) - 0.5;
          M[q1 + q5][q2 + q5] = M[q1 + q5][q2 + q5] - ( 0.5 * t * KK[j + ( w * i )][1] ) + 0.5;
          M[q2 + q5][q1 + q5] = M[q2 + q5][q1 + q5] - ( 0.5 * t * KK[j + ( w * i )][1] ) + 0.5;
          M[q4 + q5][q4 + q5] = M[q4 + q5][q4 + q5] + ( 0.5 * t * KK[j + ( w * i )][0] ) - 0.5;
          M[q1 + q5][q1 + q5] = M[q1 + q5][q1 + q5] + ( 0.5 * t * KK[j + ( w * i )][0] ) - 0.5;
          M[q4 + q5][q1 + q5] = M[q4 + q5][q1 + q5] - ( 0.5 * t * KK[j + ( w * i )][0] ) + 0.5;
          M[q1 + q5][q4 + q5] = M[q1 + q5][q4 + q5] - ( 0.5 * t * KK[j + ( w * i )][0] ) + 0.5;
        }
        else {}
      }
      q5 ++;
    }

    return M[mi][mj];
  }
  private static boolean getParity(int number) {
    boolean parity = false;
    while(number != 0) {
      parity = !parity;
      number = number & ( number - 1 );
    }
    return parity;
  }
}
