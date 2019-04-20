package feaheat2;
import java.lang.*;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.Random;
import java.lang.Math;

public class FEAGen {
  double[][] FEAGen(int N, double KK[][]) {
    int E_l = 2 * N * N;
    int ni = 4000; //Initial size of gene pool.
    int ng = 2000; //Size of gene pool (even).
    int nr = ng / 2; //Number of parents or kids.
    int nu = nr / 2;
    int ne = 4; //Number of elite states..
    double mp = 5.0 / 100; //Mutation percentage.
    int count_max = 1500; //Max number of generations.
    boolean c[][][] = new boolean[ng][E_l][2]; //Chromosomes.
    int f[] = new int[ng]; //Cost function.
    Random rnd = new Random();
    boolean d1[][][] = new boolean[ni][E_l][2];
    boolean r[][] = new boolean[E_l][2]; //Random values for vector
    boolean kbool[][] = new boolean[E_l][2]; //Thermal cond in bool.
    int e1[] = new int[ni]; //Cost of each member of initial pool
    int index1[] = new int[ni];
    int index[] = new int[ng];
    boolean d[][][] = new boolean[ng][E_l][2];
    int e[] = new int[ng];

    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < E_l; ++i) {
        if ( KK[i][j] < 20 ) {
          kbool[i][j] = false;
        }
        if ( KK[i][j] > 20 ) {
          kbool[i][j] = true;
        }
      }
    }

    for (int k = 0; k < ni; ++k) {
      for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < E_l; ++i) {
          r[i][j] = rnd.nextBoolean();
          d1[k][i][j] = r[i][j];
        }
      }
      e1[k] = (int) cost(N, r, kbool);
      index1[k] = k;
    }

    e1 = sort(N, e1, index1);
    for (int k = 0; k < ng; ++k) {
      f[k] = e1[k];
      for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < E_l; ++i) {
          c[k][i][j] = d1[index1[k]][i][j];
        }
      }
    }

    int count = 0;
    while ( count < count_max ) {
      ++count;
      for (int k = 0; k < ng; ++k) {
        for (int j = 0; j < 2; ++j) {
          for (int i = 0; i < E_l; ++i) {
            d[k][i][j] = c[k][i][j];
          }
        }
        e[k] = f[k];
        index[k] = k;
      }

      index = shuffle(index);
      int g = 0;
      for (int k = 0; k <  nr; ++k) {
        if ( e[index[g]] < e[index[g + 1]] ) {
          for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < E_l; ++i) {
              c[k][i][j] = d[index[g]][i][j];
            }
          }
          f[k] = e[index[g]];
        }
        else {
          for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < E_l; ++i) {
              c[k][i][j] = d[index[g + 1]][i][j];
            }
          }
          f[k] = e[index[g + 1]];
        }
        g += 2;
      }

      g = 0;
      for (int k = nr; k < ( nr + nu ); ++k) {
        int nx = 1 + (int) ( E_l * rnd.nextDouble() );
        for (int j = 0; j < 2; ++j) {
          for (int i = 0; i < nx; ++i) {
            c[k][i][j] = c[g][i][j];
            c[k + nu][i][j] = c[g + 1][i][j];
          }
        }
        for (int j = 0; j < 2; ++j) {
          for (int i = nx; i < E_l; ++i) {
            c[k][i][j] = c[g + 1][i][j];
            c[k + nu][i][j] = c[g][i][j];
          }
        }
        g += 2;
      }

      for (int k = 0; k < ng; ++k) {
        for (int j = 0; j < 2; ++j) {
          for (int i = 0; i < E_l; ++i) {
            d[k][i][j] = c[k][i][j];
            r[i][j] = d[k][i][j];
          }
        }

        e[k] = cost(N, r, kbool);
        index[k] = k;
      }

      e = sort(N, e, index);
      for (int k = 0; k < ng; ++k){
        f[k] = e[k];
        for (int j = 0; j < 2; ++j) {
          for (int i = 0; i < E_l; ++i) {
            c[k][i][j] = d[index[k]][i][j];
          }
        }
      }

      //Mutation in elite configuration.
      for (int k = 0; k < ne; ++k) {
        for (int j = 0; j < 2; ++j) {
          for (int i = 0; i < E_l; ++i) {
            r[i][j] = c[k][i][j];
          }
        }

        int mb = (int) ( E_l * mp + 1 );
        for (int j = 0; j < 2; ++j) {
          for (int i = 0; i < mb; ++i) {
            int ib = (int) ( E_l * rnd.nextDouble() );
            r[ib][j] = !r[ib][j];
          }
        }

        for (int i = 0; i < ng; ++i) {
          e[i] = cost(N, r, kbool);
        }
        if ( e[k] < f[k] ) {
          for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < E_l; ++i) {
              c[k][i][j] = r[i][j];
            }
          }
          f[k] = e[k];
        }
      }

      //Mutation in other configurations.
      int mmax = (int) ( ( ng - ne ) * E_l * mp + 1 );
      for (int k = 0; k < mmax; ++k) {
        for (int j = 0; j < 2; ++j) {
          int ig = (int) ( ( ( ng - ne ) * rnd.nextDouble() ) + ne );
          int ib = (int) ( E_l * rnd.nextDouble() );
          c[ig][ib][j] = !c[ig][ib][j];
        }
      }

      for (int k = 0; k < ng; ++k) {
        for (int j = 0; j < 2; ++j) {
          for (int i = 0; i < E_l; ++i) {
            d[k][i][j] = c[k][i][j];
            r[i][j] = d[k][i][j];
          }
        }
        index[k] = k;
      }

      e = sort(N, e, index);
      for (int k = 0; k < ng; ++k){
        f[k] = e[k];
        for (int j = 0; j < 2; ++j) {
          for (int i = 0; i < E_l; ++i) {
            c[k][i][j] = d[index[k]][i][j];
          }
        }
      }

      System.out.println(count + " out of " + count_max);
    }
    double[][] KKgena = export(N, c);

    return KKgena;
  }

  //----------------------------------------------------------------------------
  //Method to evaluate the cost for a given variable array.
  private static int cost(int N, boolean r[][], boolean kbool[][]) {
    int E_l = 2 * N * N;
    int g = 0;

    for (int i = 0; i < E_l; ++i) {
      if ( r[i][0] == kbool[i][0] ) {
        g += 1;
      }
      else if ( r[i][1] == kbool[i][1] ) {
        g += 1;
      }
    }

    return g;
  }

  //----------------------------------------------------------------------------
  //Method to sort array x[] from lowest to highest.
  private static int[] sort(int N, int x[], int index[]) {
    int m = x.length;

    for (int i = 0; i < m; ++i) {
      for (int j = ( i + 1 ); j < m; ++j) {
        if ( x[i] > x[j] ) {
          int xtmp = x[i];
          x[i] = x[j];
          x[j] = xtmp;

          int itmp = index[i];
          index[i] = index[j];
          index[j] = itmp;
        }
      }
    }

    return x;
  }

  //----------------------------------------------------------------------------
  //Method to shuffle the index array.
  private static int[] shuffle(int index[]) {
    int m = index.length;
    Random rnd = new Random();

    for (int i = 0; i < m; ++i) {
      int j = (int) ( m * rnd.nextDouble() );
      if ( j != i ) {
        int itmp = index[i];
        index[i] = index[j];
        index[j] = itmp;
      }
    }

    return index;
  }

  //----------------------------------------------------------------------------
  //Method for decoding resulting kbool vector.
  private static double[][] export(int N, boolean c[][][]) {
    int E_l = 2 * N * N;
    double KKgen[][] = new double[E_l][2]; //Best chromosome.

    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < E_l; ++i) {
        if ( c[0][i][j] == true ) {
          KKgen[i][j] = 1;
        }
        else if ( c[0][i][j] == false ) {
          KKgen[i][j] = 100;
        }
      }
    }

    return KKgen;
  }
}
