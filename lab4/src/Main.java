import java.io.FileNotFoundException;
import java.io.PrintStream;

public class Main
{
    public static void main(String[] args)
    {
        Task.GlobalRelaxation();
        Task.LocalRelaxation();
    }
}

class Task
{
    private static final int epsilon = 1;
    private static final double dt = 0.1;
    private static final double dx = dt, dy = dt;
    private static final int nx = 150;
    private static final int ny = 100;
    private static final double xmax = dx * nx;
    private static final double ymax = dy * ny;
    private static final double v1 = 10;
    private static final double v2 = 0;
    private static final double sigmax = 0.1*xmax;
    private static final double sigmay = 0.1*ymax;
    private static final double TOL = 1e-8;

    public static void GlobalRelaxation()
    {
        try
        {
            var console = System.out;
            double[] omegaG = {0.6, 1.0};
            for(var omega : omegaG)
            {
                var f = new PrintStream("glob_s_"+omega+".txt");
                double[][] Vn = new double[nx+1][ny+1];
                double[][] Vs = new double[nx+1][ny+1];

                // INITIALIZATION
                for(int j = 0; j <= nx; j++)
                {
                    Vn[j][0] = v1;
                    Vs[j][0] = v1;
                }

                // ITERATION
                int iter = 1;
                double S = 7500.0;
                System.setOut(f);
                System.out.println(iter + " " + S);
                System.setOut(console);
                iter++;
                double Sp;
                do
                {
                    // Vn
                    for(int i = 1; i <= nx-1; i++)
                    {
                        for(int j = 1; j <= ny-1; j++)
                        {
                            Vn[i][j] = 0.25 * (Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] + ((dt*dt)/epsilon)*ro(i, j));
                        }
                    }
                    // Border conditions
                    for(int j = 1; j <= ny - 1; j++)
                    {
                        Vn[nx][j] = Vn[nx-1][j];
                        Vn[0][j] = Vn[1][j];
                    }
                    // Vs
                    for(int i = 0; i <= nx; i++)
                    {
                        for(int j = 0; j <= ny; j++)
                        {
                            Vs[i][j] = (1.0 - omega)*Vs[i][j] + omega*Vn[i][j];
                        }
                    }

                    // Stop condition
                    Sp = S;
                    S = stop(Vn, nx, ny);

                    // Printing & Debug
                    System.setOut(f);
                    System.out.println(iter + " " + S);
                    System.setOut(console);
                    debugIterPrint(iter, 1000, S, Sp);

                    // iterator incrementation
                    iter++;
                }
                while(Math.abs((S - Sp) / Sp) > TOL);

                System.out.println("Finished with "+iter+" iterations.");
                saveGlobalToFile("glob_err_"+omega+".txt", "glob_v_"+omega+".txt", Vn);
            }
        }
        catch(FileNotFoundException ignored)
        { }
    }

    public static void LocalRelaxation()
    {
        try
        {
            var console = System.out;
            double[] omegaL = {1.0, 1.4, 1.8, 1.9};
            for(var omega : omegaL)
            {
                var f = new PrintStream("loc_s_"+omega+".txt");
                double[][] V = new double[nx+1][ny+1];

                // INITIALIZATION
                for(int j = 0; j <= nx; j++)
                {
                    V[j][0] = v1;
                }

                // ITERATION
                int iter = 1;
                double S = 7500.0;
                System.setOut(f);
                System.out.println(iter + " " + S);
                System.setOut(console);
                iter++;
                double Sp;
                do
                {
                    // Vn
                    for(int i = 1; i <= nx-1; i++)
                    {
                        for(int j = 1; j <= ny-1; j++)
                        {
                            V[i][j] = (1.0 - omega) * V[i][j] +
                                      ((0.25 * omega) *
                                      (V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1] +
                                       (dt*dt) / epsilon * ro(i, j)));
                        }
                    }
                    // Border conditions
                    for(int j = 1; j <= ny-1; j++)
                    {
                        V[nx][j] = V[nx-1][j];
                        V[0][j] = V[1][j];
                    }

                    // Stop condition
                    Sp = S;
                    S = stop(V, nx, ny);

                    // Printing & Debug
                    System.setOut(f);
                    System.out.println(iter + " " + S);
                    System.setOut(console);
                    debugIterPrint(iter, 1000, S, Sp);

                    // iterator incrementation
                    iter++;
                }
                while(Math.abs((S - Sp) / Sp) > TOL);
                System.out.println("Finished with "+iter+" iterations.");
            }
        }
        catch(FileNotFoundException ignored)
        { }
    }

    private static double ro1(double i, double j)
    {
        return Math.exp((-1.0*Math.pow(i - 0.35*xmax, 2)/(sigmax*sigmax)) -
                        (Math.pow(j - 0.5*ymax, 2)/(sigmay*sigmay)));
    }

    private static double ro2(double i, double j)
    {
        return -1.0 * Math.exp((-1.0 * Math.pow(i - 0.65*xmax, 2)/(sigmax*sigmax)) -
                         (Math.pow(j - 0.5*ymax, 2)/(sigmay*sigmay)));
    }

    private static double ro(double x, double y)
    {
         return ro1(x*dx, y*dy) + ro2(x*dx, y*dy);
    }

    private static double stop(double[][] V, int nx, int ny)
    {
        double S = 0.0;
        for(int i = 0; i <= nx - 1; i++)
        {
            for(int j = 0; j <= ny - 1; j++)
            {
                S += dt*dt *
                    ((0.5*Math.pow((V[i+1][j]-V[i][j])/dt, 2)) +
                    (0.5*Math.pow((V[i][j+1]-V[i][j])/dt, 2)) -
                    (ro(i, j)*V[i][j]));
            }
        }
        return S;
    }

    private static void debugIterPrint(int iter, int step, double S, double Sprev)
    {
        if(iter%step == 0)
        {
            System.out.println(" -- Iteration "+iter+":");
            System.out.println(" ---- S = "+S+", Sprev = "+Sprev+",");
            System.out.println(" ---- STOP = "+Math.abs((S - Sprev)/Sprev)+", TOL = "+TOL);
        }
    }

    private static void saveGlobalToFile(String filename, String filename2, double[][] V)
    {
        var console = System.out;
        try
        {
            var f1 = new PrintStream(filename);
            var f2 = new PrintStream(filename2);
            System.setOut(f1);
            for (int i = 0; i <= nx; i++)
            {
                for (int j = 0; j <= ny; j++)
                {
                    double err;
                    if(i==0 || j==0 || i==nx || j==ny)
                    {
                        err = 0;
                    }
                    else
                    {
                        err = (((V[i+1][j] - 2.0*V[i][j] + V[i-1][j]) / (dt*dt)) +
                              ((V[i][j+1] - 2.0*V[i][j] + V[i][j-1]) / (dt*dt))) +
                                     (ro(i, j) / epsilon);
                    }
                    System.out.println(i*dt+" "+j*dt+" "+err);
                }
            }
            System.setOut(f2);
            for(int i = 0; i <= nx; i++)
            {
                for(int j = 0; j <= ny; j++)
                {
                    System.out.println(i*dt+" "+j*dt+" "+V[i][j]);
                }
            }
        }
        catch(FileNotFoundException ignored)
        { }
        System.setOut(console);
    }
}
