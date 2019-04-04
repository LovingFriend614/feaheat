package feaheat2;
import java.lang.*;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.lang.Math;

public class FEAGen {
  double FEAGen(int N, double KK[][], int mi, int mj) {
    int ni = 256; //Initial size of gene pool.
    int ng = 128; //Size of gene pool (even).
    int nb = 64; //Number of bits used in encoding.
    int nc = 60; //Number of charges.
    int nv = 2 * N * N; //Number of variables.
    int nd = nv * nb; //Number of digits in chromosome.
    int nr = ng / 2; //Number of parents or kids.
    int ne = 1; //Number of elite states..
    int gmax = 500; //Max number of generations.
    double pm = 1.0 / 100; //Mutation percentage.
    double pi = (double) Math.PI;
    boolean c[][] = new boolean[ng][nd]; //Chromosomes.
    double f[] = new double[ng]; //Cost function.

    initiate();
    int g = 0;
    boolean contd = true;
    //ProgressBar pb = new ProgressBar("GenAlg", gmax);
    //pb.start();
    while ( contd && g < gmax ) {
      ++g;
      select();
      cross();
      rank();
      mutate();
      //pb.step();
      System.out.println(g + " out of " + gmax);
    }
    //pb.stop();
    export();
  }

  //----------------------------------------------------------------------------
  //Method to initialize by creating zeroth gen of gene pool
  private static void initiate() {
    Random rnd = new Random();
    boolean d[][] = new boolean[ni][nd];
    boolean w[] = new boolean[nd]; //Values from encoding
    double r[] = new double[nv];
    double e[] = new double[ni];
    int index[] = new int[ni];

    for (int i = 0; i < ni; ++i) {
      for (int j = 0; j < nv; ++j) {
        r[j] = rnd.nextDouble();
      }
      e[i] = cost(r);
      index[i] = i;
      w = encode(r,nb);
      for (int j = 0; j < nd; ++j) {
        d[i][j] = w[j];
      }
    }

    sort(e,index);
    for (int i = 0; i < ng; ++i) {
      f[i] = e[i];
      for (int j = 0; j < nd; ++j) {
        c[i][j] = d[index[i]][j];
      }
    }
  }

  //----------------------------------------------------------------------------
  //Method to evaluate the cost for a given variable array.
  private static double cost(double r[]) {
    double g = 0;
    int k = 0;
    double V_0 = 1.09e3;
    double r_0 = 0.330;
    double q_e = 1;
    double epsilon_0 = 8.85e-12;
    double ri[] = new double[nc];
    double theta[] = new double[nc];
    double phi[] = new double[nc];
    double k_e = - ( q_e * q_e ) / ( 4 * pi * epsilon_0 );

    for (int i = 0; i < nc; ++i) {
      ri[i] = r_0 * (double) r[k];
      theta[i] = pi * r[k + 1];
      phi[i] = 2 * pi * r[k + 2];
      k += 3;
    }

    for (int i = 0; i < ( nc - 1 ); ++i) {
      double xi = ri[i] * Math.sin(theta[i]) * Math.cos(phi[i]);
      double yi = ri[i] * Math.sin(theta[i]) * Math.sin(phi[i]);
      double zi = ri[i] * Math.cos(theta[i]);
      for ( int j = ( i + 1 ); j < nc; ++j) {
        double dx = xi - ri[j] * Math.sin(theta[j]) * Math.cos(phi[j]);
        double dy = yi - ri[j] * Math.sin(theta[j]) * Math.sin(phi[j]);
        double dz = zi - ri[j] * Math.cos(theta[j]);
        double rf = Math.sqrt( ( dx * dx ) + ( dy * dy ) + ( dz * dz ) );
        g += ( k_e / rf ) + ( V_0 * Math.exp( - rf / r_0 ) );
      }
    }
    return g;
  }

  //----------------------------------------------------------------------------
  //Method to encode real array into binary array
  private static boolean[] encode(double r[], int m) {
    int n = r.length;
    boolean w[] = new boolean[ n * m ];

    for (int i = 0; i < n; ++i) {
      double sum = r[i];
      w[ i * m ] = false;
      if ( (int) ( 0.5 + sum ) == 1 ) {
        w[ i * m ] = true;
      }
      double d = 2;

      for (int j = 1; j < m; ++j) {
        if ( w[ i * m + j - 1 ] == true ) {
          sum -= 1 / d;
        }
        w[ i * m + j ] = false;
        if ( (int) (0.5 + d * sum ) == 1 ) {
          w[ i * m + j ] = true;
        }
        d *= 2;
      }
    }
    return w;
  }

  //----------------------------------------------------------------------------
  //Method to sort array x[] from lowest to highest.
  private static void sort(double x[], int index[]) {
    int m = x.length;

    for (int i = 0; i < m; ++i) {
      for (int j = ( i + 1 ); j < m; ++j) {
        if ( x[i] > x[j] ) {
          double xtmp = x[i];
          x[i] = x[j];
          x[j] = xtmp;

          int itmp = index[i];
          index[i] = index[j];
          index[j] = itmp;
        }
      }
    }
  }

  //----------------------------------------------------------------------------
  //Method to run tournaments in selecting parents.
  private static void select() {
    int index[] = new int[ng];
    boolean d[][] = new boolean[ng][nd];
    double e[] = new double[ng];

    for (int i = 0; i < ng; ++i) {
      for (int j = 0; j < nd; ++j) {
        d[i][j] = c[i][j];
      }
      e[i] = f[i];
      index[i] = i;
    }

    shuffle(index);
    int k = 0;
    for (int i = 0; i <  nr; ++i) {
      if ( e[index[k]] < e[index[ k + 1 ]] ) {
        for (int j = 0; j < nd; ++j) {
          c[i][j] = d[index[k]][j];
        }
        f[i] = e[index[k]];
      }
      else {
        for (int j = 0; j < nd; ++j) {
          c[i][j] = d[index[ k + 1 ]][j];
        }
        f[i] = e[index[ k + 1 ]];
      }
      k += 2;
    }
  }

  //----------------------------------------------------------------------------
  //Method to shuffle the index array.
  private static void shuffle(int index[]) {
    int k = index.length;
    Random rnd = new Random();
    for (int i = 0; i < k; ++i) {
      int j = (int) ( k * rnd.nextDouble() );
      if ( j != i ) {
        int itmp = index[i];
        index[i] = index[j];
        index[j] = itmp;
      }
    }
  }

  //----------------------------------------------------------------------------
  //Method to perform single-point crossover operations.
  private static void cross() {
    int k = 0;
    int nu = nr / 2;
    Random rnd = new Random();

    for (int i = nr; i < ( nr + nu ); ++i) {
      int nx = 1 + (int) ( nd * rnd.nextDouble() );
      for (int j = 0; j < nx; ++j) {
        c[i][j] = c[k][j];
        c[i + nu][j] = c[k + 1][j];
      }
      for (int j = nx; j < nd; ++j) {
        c[i][j] = c[k + 1][j];
        c[i + nu][j] = c[k][j];
      }
      k += 2;
    }
  }

  //----------------------------------------------------------------------------
  //Method to rank chromosomes in the population.
  private static void rank() {
    boolean d[][] = new boolean[ng][nd];
    double r[] = new double[nv];
    double e[] = new double[ng];
    int index[] = new int[ng];

    for (int i = 0; i < ng; ++i) {
      for (int j = 0; j < nd; ++j) {
        d[i][j] = c[i][j];
      }

      r = decode(d[i],nb);
      e[i] = cost(r);
      index[i] = i;
    }

    sort(e,index);
    for (int i = 0; i < ng; ++i){
      f[i] = e[i];
      for (int j = 0; j < nd; ++j) {
        c[i][j] = d[index[i]][j];
      }
    }
  }

  //----------------------------------------------------------------------------
  //Method to array an array.
  private static double[] decode(boolean w[], int m) {
    int n = w.length/m;
    double r[] = new double[n];
    for (int i = 0; i < n; ++i) {
      double d = 2;
      double sum = 0;
      for (int j = 0; j < m; ++j) {
        if ( w[i * m + j] == true ) {
          sum += 1/d;
        }
        d *= 2;
      }
      r[i] = sum + 1/d;
    }
    return r;
  }

  //----------------------------------------------------------------------------
  //Method to mutate a percentage of bits in chromosomes (not best one).
  private static void mutate() {
    Random rnd = new Random();
    double r[] = new double[nv];
    boolean w[] = new boolean[nd];

    //Mutation in elite configuration.
    for (int i = 0; i < ne; ++i) {
      for (int j = 0; j < nd; ++j) {
        w[j] = c[i][j];
      }

      int mb = (int) ( nd * pm + 1 );
      for (int j = 0; j < mb; ++j) {
        int ib = (int) ( nd * rnd.nextDouble() );
        w[ib] = !w[ib];
      }

      r = decode(w,nb);
      double e = cost(r);
      if ( e < f[i] ) {
        for (int j = 0; j < nd; ++j) {
          c[i][j] = w[j];
        }
        f[i] = e;
      }
    }

    //Mutation in other configurations.
    int mmax = (int) ( ( ng - ne ) * nd * pm + 1 );
    for (int i = 0; i < mmax; ++i) {
      int ig = (int) ( ( ( ng - ne ) * rnd.nextDouble() ) + ne );
      int ib = (int) ( nd * rnd.nextDouble() );
      c[ig][ib] = !c[ig][ib];
    }

    //Rank the chromosomes.
    rank();
  }

  //----------------------------------------------------------------------------
  //Method to save coordinates and output cost function of best chromosome.
  private static void export() throws FileNotFoundException {
    int k = 0;
    double r_0 = 0.330;
    double r[] = new double[nv];
    boolean w[] = new boolean[nd];
    double ri[] = new double[nc];
    double theta[] = new double[nc];
    double phi[] = new double[nc];
    double x[] = new double[nc];
    double y[] = new double[nc];
    double z[] = new double[nc];
    double rescale = Math.sqrt( nc / 4.0 );

    // Write out the coordinates of the best chromosome
    PrintWriter q = new PrintWriter
      (new FileOutputStream("genalgdata.txt"), true);
    for (int i = 0; i < nd; ++i) {
      w[i] = c[0][i];
    }

    r = decode(w,nb);
    for (int i = 0; i < nc; ++i) {
      ri[i] = r_0 * r[k];
      theta[i] = pi * r[k + 1];
      phi[i] = 2 * pi * r[k + 2];
      x[i] = rescale * ri[i] * Math.sin(theta[i]) * Math.cos(phi[i]);
      y[i] = rescale * ri[i] * Math.sin(theta[i]) * Math.sin(phi[i]);
      z[i] = rescale * ri[i] * Math.cos(theta[i]);
      q.println(x[i] + "  " + y[i] + "  " + z[i]);
      k += 3;
    }

    // Print out the energy of the best chromosome
    for (int i = 0; i < ng; ++i) {
      System.out.println("Length: " + f[i]);
    }
  }
}
