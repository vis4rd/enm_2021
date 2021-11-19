import java.io.PrintStream;
import java.io.File;
import java.io.FileNotFoundException;

public class Zadanie
{
    public static final double x_0 = 0.01;
    public static final double v_0 = 0.0;
    public static final double dt_0 = 1.0;
    public static final double S = 0.75;
    public static final int p = 2;
    public static final double t_max = 40.0;
    public static final double alpha = 5.0;
    public static final double[] TOL = {0.01, 0.00001};

    private static double[] calcTrapezy(double x, double v, double dt)
    {
        double delta = Math.pow(10, -10);
        double xp = x;
        double vp = v;

        while(true)
        {
            double F = F(x, v, xp, vp, dt);
            double G = G(x, v, xp, vp, dt);

            double a11 = 1.0;
            double a12 = -dt/2.0;
            double a21 = -(dt/2.0) * (-2.0 * alpha * x * v - 1.0);
            double a22 = 1.0 - (dt/2.0) * alpha * (1.0 - x*x);

            double dx = (-F * a22 - ((-G)*a12)) / (a11*a22 - a12*a21);
            double dv = (a11 * (-G) - a21*(-F)) / (a11*a22 - a12*a21);

            x += dx;
            v += dv;
            if(Math.abs(dx) < delta && Math.abs(dv) < delta) { break; }
        }
        double[] result = {x, v};
        return result;
    }

    public static void Trapezy() throws FileNotFoundException
    {      
        PrintStream[] ps = new PrintStream[2];
        ps[0] = new PrintStream(new File("trapez.txt"));
        ps[1] = new PrintStream(new File("trapez2.txt"));
        PrintStream console = System.out;

        for(int i = 0; i < 2; i++)
        {
            var ct = TOL[i];
            var o = ps[i];
            
            double x = x_0;
            double v = v_0;
            double dt = dt_0;
            double t = 0.0;
            
            while(t < t_max)
            {
                var first = calcTrapezy(x, v, dt);
                var second = calcTrapezy(first[0], first[1], dt);
                var doubl = calcTrapezy(x, v, 2*dt);

                double err_x = (second[0] - doubl[0]) / (Math.pow(2, p) - 1.0);
                double err_v = (second[1] - doubl[1]) / (Math.pow(2, p) - 1.0);

                double err_max = Math.max(Math.abs(err_x), Math.abs(err_v));

                if(err_max < ct)
                {
                    t += (2.0 * dt);
                    x = second[0];
                    v = second[1];
                    System.setOut(o);
                    System.out.println(t + " " + dt + " " + x + " " + v);
                    System.setOut(console);
                }
                dt *= Math.pow((S * ct) / err_max, (1.0 / (p + 1.0)));         
            }
        }
    }

    private static double[] calcRK2(double x, double v, double dt)
    {
        double k1x = f(v);
        double k1v = g(x, v);
        double k2x = f(v) + dt*k1v;
        double k2v = alpha * (1.0 - Math.pow(x + dt*k1x, 2)) * (v + dt*k1v) - (x + dt*k1x);

        double[] result = {x + dt/2.0*(k1x + k2x), v + dt/2.0*(k1v + k2v)};
        return result;
    }

    public static void RK2() throws FileNotFoundException
    {      
        PrintStream[] ps = new PrintStream[2];
        ps[0] = new PrintStream(new File("RK2.txt"));
        ps[1] = new PrintStream(new File("RK22.txt"));
        PrintStream console = System.out;

        for(int i = 0; i < 2; i++)
        {
            var ct = TOL[i];
            var o = ps[i];
            
            double x = x_0;
            double v = v_0;
            double dt = dt_0;
            double t = 0.0;
            
            while(t < t_max)
            {
                var first = calcRK2(x, v, dt);
                var second = calcRK2(first[0], first[1], dt);
                var doubl = calcRK2(x, v, 2*dt);

                double err_x = (second[0] - doubl[0]) / (Math.pow(2, p) - 1.0);
                double err_v = (second[1] - doubl[1]) / (Math.pow(2, p) - 1.0);

                double err_max = Math.max(Math.abs(err_x), Math.abs(err_v));

                if(err_max < ct)
                {
                    t += (2.0 * dt);
                    x = second[0];
                    v = second[1];
                    System.setOut(o);
                    System.out.println(t + " " + dt + " " + x + " " + v);
                    System.setOut(console);
                }
                dt *= Math.pow((S * ct) / err_max, (1.0 / (p + 1.0)));         
            }
        }
    }

    private static double f(double v) { return v; }
    private static double g(double x, double v) { return alpha*(1.0 - x*x)*v - x; }

    private static double F(double x, double v, double xp, double vp, double dt) { return x - xp - (dt/2.0)*(f(vp) + f(v)); }
    private static double G(double x, double v, double xp, double vp, double dt) { return v - vp - (dt/2.0)*(g(xp, vp) + g(x, v)); }
}

