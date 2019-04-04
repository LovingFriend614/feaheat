package feaheat;
import java.lang.*;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.lang.Math;

public class BoundaryMatrix {
  double BoundaryMatrix(int N, double t, double x_size, double y_size, int numEle, double k_xx, double k_yy, int K_list[], double K_val[][], int mi, int mj) {

    //--------------------------------------------------------------------------
    //Definitions
    int sqn = N + 1;
    int n = 2 * sqn;
    double x_ele = x_size, y_ele = y_size / N;
    int w = 2 * N;
    double A = x_size * y_size;
    double A_ele = A / w; //Area per element
    double z = 2 * A_ele;
    double x[][] = new double[3][w];
    double y[][] = new double[3][w];
    int nd[][] = new int[3][w];
    double ap[] = new double[w];
    double bp[] = new double[w];
    double cp[] = new double[w];
    double dp[] = new double[w];
    double ep[] = new double[w];
    double fp[] = new double[w];
    double KK[][] = new double[w][2];
    double M[][] = new double[n][n];
    double W[][][] = new double[3][3][w];
    double I[][][] = new double[n][n][w];
    double B[][][] = new double[2][3][w];
    double BT[][][] = new double[3][2][w];
    double DB[][][] = new double[2][3][w];

    //--------------------------------------------------------------------------
    //Vectors
    for (int b = 0; b < 3; ++b) {
      for (int a = 0; a < w; ++a) {
        x[b][a] = xVector(x_ele, w, a, b);
      }
    }
    for (int b = 0; b < 3; ++b) {
      for (int a = 0; a < w; ++a) {
        y[b][a] = yVector(y_ele, w, a, b);
      }
    }
    for (int b = 0; b < 3; ++b) {
      for (int a = 0; a < w; ++a) {
        nd[b][a] = ndVector(w, a, b);
      }
    }

    //--------------------------------------------------------------------------
    //Matrix Operations
    //B and BT matrix elements
    for (int i = 0; i < w; ++i) {
      ap[i] = y[1][i] - y[2][i];
      bp[i] = y[2][i] - y[0][i];
      cp[i] = y[0][i] - y[1][i];
      dp[i] = x[2][i] - x[1][i];
      ep[i] = x[0][i] - x[2][i];
      fp[i] = x[1][i] - x[0][i];
    }

    //K matrix
    for (int i = 0; i < w; ++i) {
      KK[i][0] = k_xx;
      KK[i][1] = k_yy;
    }

    /*for (int i = 0; i < numEle; ++i) {
      KK[K_list[i]][0] = K_val[i][0];
      KK[K_list[i]][1] = K_val[i][1];
    }*/

    //B and BT matrices
    for (int i = 0; i < w; ++i) {
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
    for (int i = 0; i < w; ++i) {
      DB[0][0][i] = KK[i][0] * B[0][0][i] * t * A_ele;
      DB[0][1][i] = KK[i][0] * B[0][1][i] * t * A_ele;
      DB[0][2][i] = KK[i][0] * B[0][2][i] * t * A_ele;
      DB[1][0][i] = KK[i][1] * B[1][0][i] * t * A_ele;
      DB[1][1][i] = KK[i][1] * B[1][1][i] * t * A_ele;
      DB[1][2][i] = KK[i][1] * B[1][2][i] * t * A_ele;
    }

    //Computing BT * D * B
    for (int k = 0; k < w; ++k) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          W[i][j][k] = ( BT[i][0][k] * DB[0][j][k] ) + ( BT[i][1][k] * DB[1][j][k] );
        }
      }
    }

    //Indexing matrix
    for (int k = 0; k < w; ++k) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          I[nd[i][k] - 1][nd[j][k] - 1][k] = W[i][j][k];
        }
      }
    }

    //Boundary Stiffness matrix
    for (int k = 0; k < w; ++k) {
      for (int i = 0; i < n; ++i) { //I0
        for (int j = 0; j < n; ++j) {
          M[i][j] = M[i][j] + I[i][j][k];
        }
      }
    }

    return M[mi][mj];
  }
  private static double xVector(double x_ele, int w, int a, int b) {
    //Defines exact position of all x components of all nodes
    double x[][] = new double[3][w];
    x[0][0] = 0;
    x[0][1] = 0;
    x[1][0] = 0;
    x[1][1] = x_ele;
    x[2][0] = x_ele;
    x[2][1] = x_ele;

    if (w > 2) {
      for (int k = 0; k < 3; ++k) {
        for (int i = 2; i < w; i += 2) {
          x[k][i] = x[k][i - 2];
          x[k][i + 1] = x[k][i - 1];
        }
      }
    }

    return x[b][a];
  }
  private static double yVector(double y_ele, int w, int a, int b) {
    //Defines exact position of all y components of all nodes
    double y[][] = new double[3][w];
    y[0][0] = 0;
    y[0][1] = 0;
    y[1][0] = y_ele;
    y[1][1] = 0;
    y[2][0] = y_ele;
    y[2][1] = y_ele;

    if (w > 2) {
      for (int k = 0; k < 3; ++k) {
        for (int i = 2; i < w; i += 2) {
          y[k][i] = y[k][i - 2] + y_ele;
          y[k][i + 1] = y[k][i - 1] + y_ele;
        }
      }
    }

    return y[b][a];
  }
  private static int ndVector(int w, int a, int b) {
    //Defines which nodes correspond to each induvidual element
    int nd[][] = new int[3][w];
    nd[0][0] = 1;
    nd[1][0] = 3;
    nd[2][0] = 4;
    nd[0][1] = 1;
    nd[1][1] = 2;
    nd[2][1] = 4;

    if (w > 2) {
      for (int i = 2; i < w; i += 2) {
        nd[0][i] = nd[0][i - 2] + 2;
        nd[1][i] = nd[1][i - 2] + 2;
        nd[2][i] = nd[2][i - 2] + 2;
        nd[0][i + 1] = nd[0][i - 1] + 2;
        nd[1][i + 1] = nd[1][i - 1] + 2;
        nd[2][i + 1] = nd[2][i - 1] + 2;
      }
    }

    return nd[b][a];
  }
}
