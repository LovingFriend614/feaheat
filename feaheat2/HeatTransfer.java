package feaheat2;
import java.lang.*;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.lang.Math;

public class HeatTransfer {
  public static void main(String argv[]) throws FileNotFoundException {

    File inputFile = new File("InputFile.txt");
    Scanner input = new Scanner(inputFile);
    Scanner fileInput = new Scanner(System.in);
    System.out.println("Name of output file to be generated: ");
    String fileName = fileInput.next();
    PrintWriter outputFile = new PrintWriter(fileName + ".csv");
    BoundaryFlux flux = new BoundaryFlux();
    StiffnessMatrix smatrix = new StiffnessMatrix();
    long timeInitial = System.nanoTime();

    int N = (int) input.nextInt();
    int n = ( N + 1 ) * ( N + 1 );
    int qn = (int) 2 + ( 2 * N );
    int sqn = N + 1;
    int nn = 2 * sqn;
    int E_l = 2 * N * N;
    int w = E_l / N;
    int p = n - 1;
    int pp = nn - 1;
    int count = N + 2;

    double Tb[] = new double[nn]; //Boundary temperature
    double Tm[] = new double[sqn]; //Measured temperature
    double T[] = new double[n]; //Full temperature vector
    double fb[] = new double[nn]; //Boundary flux
    double f[] = new double[n]; //Flux along all sample
    double M[][] = new double[n][n]; //Stiffness matrix
    double L[][] = new double[n][n]; //Lower triangular matrix
    double U[][] = new double[n][n]; //Upper triangular matrix
    double Y[] = new double[n];

    double x_size = (double) input.nextDouble();
    double y_size = (double) input.nextDouble();
    double t = (double) input.nextDouble();
    double k_xx = (double) input.nextDouble();
    double k_yy = (double) input.nextDouble();
    double Ts = (double) input.nextDouble();

    for (int i = 0; i < sqn; ++i) {
      Tm[i] = (double) input.nextDouble();
    }

    int numEle = input.nextInt();
    int K_list[] = new int[numEle];
    double K_val[][] = new double[numEle][2];

    if ( numEle != 0 ) {
      for (int i = 0; i < numEle; ++i) {
        K_list[i] = input.nextInt() - 1;
      }
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < numEle; ++j) {
          K_val[j][i] = input.nextDouble();
        }
      }
    }
    else {
      for (int i = 0; i < numEle; ++i) {
        K_list[i] = 0;
      }
      for (int i = 0; i < numEle; ++i) {
        K_val[i][0] = k_xx;
      }
      for (int i = 0; i < numEle; ++i) {
        K_val[i][1] = k_yy;
      }
    }

    //--------------------------------------------------------------------------
    //Boundary temperature vector
    int Tnd = 0;
    for (int i = 0; i < nn; i += 2) {
      Tb[i] = Ts;
      Tb[i + 1] = Tm[Tnd];
      Tnd = Tnd + 1;
    }

    //--------------------------------------------------------------------------
    //Boundary flux vector
    long timeBFi = System.nanoTime();
    for (int i = 0; i < nn; ++i) {
      fb[i] = flux.BoundaryFlux(N, t, x_size, y_size, numEle, k_xx, k_yy, K_list, K_val, Tb, i);
    }
    long timeBFf = System.nanoTime();
    double boundaryFluxTime = (double) ( timeBFf - timeBFi ) / 1000000000;

    //--------------------------------------------------------------------------
    //Flux vector expanded for full mesh
    f[0] = fb[0];
    f[p] = fb[pp];

    int fnd = 1;
    for (int i = N; i < p; i += sqn) {
      f[i] = fb[fnd];
      f[i + 1] = fb[fnd + 1];
      fnd = fnd + 2;
    }

    //--------------------------------------------------------------------------
    //Stiffness matrix and LU matrices
    long timeMatrixi = System.nanoTime();
    long stiffnessi = System.nanoTime();
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        M[i][j] = smatrix.StiffnessMatrix(N, t, x_size, y_size, numEle, k_xx, k_yy, K_list, K_val, i, j);
      }
    }
    long stiffnessf = System.nanoTime();
    long lui = System.nanoTime();

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

    long luf = System.nanoTime();
    long timeMatrixf = System.nanoTime();
    double matrixTime = (double) ( timeMatrixf - timeMatrixi ) / 1000000000;
    double stiffnessTime = (double) ( stiffnessf - stiffnessi ) / 1000000000;
    double luTime = (double) ( luf - lui ) / 1000000000;
    //--------------------------------------------------------------------------
    //Calculating Y vector for LU
    Y[0] = f[0];

    for (int i = 1; i < n; ++i) {
      double y_sum = 0;
      for (int j = 0; j < i; ++j) {
        y_sum = y_sum + ( L[i][j] * Y[j] );
      }
      Y[i] = f[i] - y_sum;
    }

    //--------------------------------------------------------------------------
    //Calculating temperature vector with U matrix
    T[p] = Tb[pp];

    for (int i = p - 1; i >= 0; --i) {
      double T_sum = 0;
      for (int j = i + 1; j < n; ++j) {
        T_sum = T_sum + ( U[i][j] * T[j] );
      }
      T[i] = ( Y[i] - T_sum ) / U[i][i];
    }

    //--------------------------------------------------------------------------
    //Print temperature to .csv file
    /*for (int i = 0; i < n; ++i) { //Matrix print M
      for (int j = 0; j < n; ++j) {
        outputFile.print(round(M[i][j], 2) + "\t");
      }
      outputFile.print("\n");
    }
    outputFile.print("\n");
    for (int i = 0; i < n; ++i) { //Matrix print U
      for (int j = 0; j < n; ++j) {
        outputFile.print(round(U[i][j], 2) + "\t");
      }
      outputFile.print("\n");
    }
    outputFile.print("\n");
    for (int i = 0; i < n; ++i) { //Matrix print L
      for (int j = 0; j < n; ++j) {
        outputFile.print(round(L[i][j], 2) + "\t");
      }
      outputFile.print("\n");
    }*/

    int Tout = 0;
    for (int i = 0; i < sqn; ++i) {
      for (int j = 0; j < sqn; ++j) {
        outputFile.print(round(T[j + Tout], 2) + "\t");
      }
      Tout = Tout + sqn;
      outputFile.print("\n");
    }
    outputFile.close();

    long timeFinal = System.nanoTime();
    double genTime = (double) ( timeFinal - timeInitial ) / 1000000000;
    System.out.println("Boundary matrix generated in " + boundaryFluxTime + " seconds");
    System.out.println("Stiffness matrix and LU generated in " + matrixTime + " seconds");
    System.out.println("Stiffness matrix generated in " + stiffnessTime + " seconds");
    System.out.println("LU generated in " + luTime + " seconds");
    System.out.println("File generated in " + genTime + " seconds");
  }
  private static double round(double value, int precision) {
    int scale = (int) Math.pow(10, precision);
    return (double) Math.round(value * scale) / scale;
  }
  private static double LUDecomp(int N, double M[][], int mi, int mj, int side) {
    int n = ( N + 1 ) * ( N + 1 );
    double L[][] = new double[n][n];
    double U[][] = new double[n][n];
    double O[][] = new double[n][n];

    for (int j = 0; j < n; ++j) {
      U[0][j] = M[0][j];
    }
    for (int i = 0; i < n; ++i) {
      L[i][0] = M[i][0] / U[0][0];
    }
    for (int i = 0; i < n; ++i) {
      L[i][i] = 1;
    }
    for (int i = 1; i < n; ++i) {
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
