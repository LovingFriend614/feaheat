package feaheat;
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
    LUDecomp lu = new LUDecomp();
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
      fb[i] = flux.BoundaryFlux(N, t, x_size, y_size, k_xx, k_yy, Tb, i);
    }
    long timeBFf = System.nanoTime();

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
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        M[i][j] = smatrix.StiffnessMatrix(N, t, x_size, y_size, k_xx, k_yy, i, j);
      }
    }
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        L[i][j] = lu.LUDecomp(N, M, i, j, 1);
      }
    }
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        U[i][j] = lu.LUDecomp(N, M, i, j, 2);
      }
    }
    long timeMatrixf = System.nanoTime();

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
    System.out.println("File generated in " + genTime + " seconds");
  }
  private static double round(double value, int precision) {
    int scale = (int) Math.pow(10, precision);
    return (double) Math.round(value * scale) / scale;
  }
}

