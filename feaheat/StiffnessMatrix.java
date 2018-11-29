package feaheat;
import java.lang.*;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.lang.Math;

public class StiffnessMatrix {
  double StiffnessMatrix(int N, double t, double x_size, double y_size, double k_xx, double k_yy, int mi, int mj) {

    //--------------------------------------------------------------------------
    //Definitions
    int sqn = N + 1;
    int n = sqn * sqn;
    double x_ele = x_size / N, y_ele = y_size / N;
    int E_l = 2 * N * N;
    int w = 2 * N;
    double A = x_size * y_size;
    double A_ele = A / E_l; //Area per element
    double z = 2 * A_ele;
    double x[][] = new double[3][E_l];
    double y[][] = new double[3][E_l];
    int nd[][] = new int[3][E_l];
    double ap[] = new double[E_l];
    double bp[] = new double[E_l];
    double cp[] = new double[E_l];
    double dp[] = new double[E_l];
    double ep[] = new double[E_l];
    double fp[] = new double[E_l];
    double D[][] = new double[2][2];
    double M[][] = new double[n][n];
    double W[][][] = new double[3][3][E_l];
    double I[][][] = new double[n][n][E_l];
    double B[][][] = new double[2][3][E_l];
    double BT[][][] = new double[3][2][E_l];
    double DB[][][] = new double[2][3][E_l];

    //--------------------------------------------------------------------------
    //Vectors
    for (int b = 0; b < 3; ++b) {
      for (int a = 0; a < E_l; ++a) {
        x[b][a] = xVector(N, x_ele, E_l, a, b);
      }
    }
    for (int b = 0; b < 3; ++b) {
      for (int a = 0; a < E_l; ++a) {
        y[b][a] = yVector(N, y_ele, E_l, a, b);
      }
    }
    for (int b = 0; b < 3; ++b) {
      for (int a = 0; a < E_l; ++a) {
        nd[b][a] = ndVector(N, E_l, a, b);
      }
    }

    //--------------------------------------------------------------------------
    //Matrix Operations
    //B and BT matrix elements
    for (int i = 0; i < E_l; ++i) {
      ap[i] = y[1][i] - y[2][i];
      bp[i] = y[2][i] - y[0][i];
      cp[i] = y[0][i] - y[1][i];
      dp[i] = x[2][i] - x[1][i];
      ep[i] = x[0][i] - x[2][i];
      fp[i] = x[1][i] - x[0][i];
    }

    //D matrix
    D[0][0] = k_xx;
    D[1][1] = k_yy;

    //B and BT matrices
    for (int i = 0; i < E_l; ++i) {
      B[0][0][i] = ap[i] / z;
      B[0][1][i] = bp[i] / z;
      B[0][2][i] = cp[i] / z;
      B[1][0][i] = dp[i] / z;
      B[1][1][i] = ep[i] / z;
      B[1][2][i] = fp[i] / z;

      BT[0][0][i] = ap[i] / z;
      BT[1][0][i] = bp[i] / z;
      BT[2][0][i] = cp[i] / z;
      BT[0][1][i] = dp[i] / z;
      BT[1][1][i] = ep[i] / z;
      BT[2][1][i] = fp[i] / z;
    }

    //D by B matrix
    for (int i = 0; i < E_l; ++i) {
      DB[0][0][i] = D[0][0] * B[0][0][i] * t * A_ele;
      DB[0][1][i] = D[0][0] * B[0][1][i] * t * A_ele;
      DB[0][2][i] = D[0][0] * B[0][2][i] * t * A_ele;
      DB[1][0][i] = D[1][1] * B[1][0][i] * t * A_ele;
      DB[1][1][i] = D[1][1] * B[1][1][i] * t * A_ele;
      DB[1][2][i] = D[1][1] * B[1][2][i] * t * A_ele;
    }

    //Computing BT * D * B
    for (int k = 0; k < E_l; ++k) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          W[i][j][k] = ( BT[i][0][k] * DB[0][j][k] ) + ( BT[i][1][k] * DB[1][j][k] );
        }
      }
    }

    //Indexing matrix
    for (int k = 0; k < E_l; ++k) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          I[nd[i][k] - 1][nd[j][k] - 1][k] = W[i][j][k];
        }
      }
    }

    //Stiffness matrix
    for (int k = 0; k < E_l; ++k) {
      for (int i = 0; i < n; ++i) { //I0
        for (int j = 0; j < n; ++j) {
          M[i][j] = M[i][j] + I[i][j][k];
        }
      }
    }

    return M[mi][mj];
  }
  private static double xVector(int N, double x_ele, int E_l, int a, int b) {
    //Defines exact position of all x components of all nodes
    double x[][] = new double[3][E_l];
    int w = E_l / N;
    x[0][0] = 0;
    x[0][1] = 0;
    x[1][0] = 0;
    x[1][1] = x_ele;
    x[2][0] = x_ele;
    x[2][1] = x_ele;

    if (E_l > 2) {
      for (int k = 0; k < 3; ++k) {
        for (int i = 2; i < w; i += 2) {
          x[k][i] = x[k][i - 2] + x_ele;
          x[k][i + 1] = x[k][i - 1] + x_ele;
        }
        for (int i = 1; i < N; ++i) {
          for (int j = 0; j < w; ++j) {
            x[k][( i * w ) + j] = x[k][j];
          }
        }
      }
    }

    return x[b][a];
  }
  private static double yVector(int N, double y_ele, int E_l, int a, int b) {
    //Defines exact position of all y components of all nodes
    double y[][] = new double[3][E_l];
    int w = E_l / N;
    y[0][0] = 0;
    y[0][1] = 0;
    y[1][0] = y_ele;
    y[1][1] = 0;
    y[2][0] = y_ele;
    y[2][1] = y_ele;

    if (E_l > 2) {
      for (int k = 0; k < 3; ++k) {
        for (int i = 2; i < w; i += 2) {
          y[k][i] = y[k][i - 2];
          y[k][i + 1] = y[k][i - 1];
        }
        for (int i = 1; i < N; ++i) {
          for (int j = 0; j < w; ++j) {
            y[k][( i * w ) + j] = y[k][( ( i - 1 ) * w ) + j] + y_ele;
          }
        }
      }
    }

    return y[b][a];
  }
  private static int ndVector(int N, int E_l, int a, int b) {
    //Defines which nodes correspond to each induvidual element
    int nd[][] = new int[3][E_l];
    int w = E_l / N;
    int sqn = N + 1;
    nd[0][0] = 1;
    nd[0][1] = 1;
    nd[1][0] = N + 2;
    nd[1][1] = 2;
    nd[2][0] = N + 3;
    nd[2][1] = N + 3;

    if (E_l > 2) {
      for (int k = 0; k < 3; ++k) {
        for (int i = 2; i < w; i += 2) {
          nd[k][i] = nd[k][i - 2] + 1;
          nd[k][i + 1] = nd[k][i - 1] + 1;
        }
        for (int i = 1; i < N; ++i) {
          for (int j = 0; j < w; ++j) {
            nd[k][( i * w ) + j] = nd[k][( ( i - 1 ) * w ) + j] + sqn;
          }
        }
      }
    }

    return nd[b][a];
  }
}
