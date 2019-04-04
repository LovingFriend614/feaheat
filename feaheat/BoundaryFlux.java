package feaheat;
import java.lang.*;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.lang.Math;

public class BoundaryFlux {
  double BoundaryFlux(int N, double t, double x_size, double y_size, int numEle, double k_xx, double k_yy, int K_list[], double K_val[][], double T[], int fi) {

    BoundaryMatrix bm = new BoundaryMatrix();
    int sqn = N + 1;
    int n = 2 * sqn;
    double M[][] =new double[n][n];
    double f[] = new double[n];

    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        M[i][j] = bm.BoundaryMatrix(N, t, x_size, y_size, numEle, k_xx, k_yy, K_list, K_val, i, j);
      }
    }
    for (int i = 0; i < n; ++i) {
      double f_sum = 0;
      for (int j = 0; j < n; ++j) {
        f_sum = f_sum + ( M[i][j] * T[j] );
      }
      f[i] = f_sum;
    }

    return f[fi];
  }
}
