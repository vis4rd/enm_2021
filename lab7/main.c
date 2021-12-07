#include <stdio.h>
#include <math.h>
#include <string.h>

static const double DT = 0.01;
static const double DX = DT;
static const double DY = DT;
static const double RHO = 1;
static const double MU = 1;
static const unsigned NX = 200u;
static const unsigned NY = 90u;
static const unsigned I1 = 50u;
static const unsigned J1 = 55u;
static const unsigned IT_MAX = 20000u;

double get_x(const unsigned i);
double get_y(const unsigned j);
double get_q_out(const double q_in, const double *y);
void bc_psi(double psi[NX+1][NY+1], const double *y, const double q_in, const double q_out);
double bigGamma(double psi[NX+1][NY+1], double zeta[NX+1][NY+1]);
void relaxation(double psi[NX+1][NY+1],
    double zeta[NX+1][NY+1],
    double v_x[NX+1][NY+1],
    double v_y[NX+1][NY+1],
    int q_in);

int main(void)
{
    double psi[NX+1][NY+1];
    double zeta[NX+1][NY+1];
    double v_x[NX+1][NY+1];
    double v_y[NX+1][NY+1];
    relaxation(psi, zeta, v_x, v_y, -1000);
    relaxation(psi, zeta, v_x, v_y, -4000);
    relaxation(psi, zeta, v_x, v_y, 4000);    
}

double get_x(const unsigned i)
{
    return DT * i;
}

double get_y(const unsigned j)
{
    return DT * j;
}

double get_q_out(const double q_in, const double *y)
{
    return q_in *
        ( ( pow(y[NY], 3) - pow(y[J1], 3) - 
            3.0 * y[J1] * pow(y[NY], 2) +
            3.0 * pow(y[J1], 2) * y[NY] )
            / pow(y[NY], 3) );
}

void bc_psi(double psi[NX+1][NY+1], const double *y, const double q_in, const double q_out)
{
    for(unsigned j = J1; j <= NY; j++)  // BORDER A
        psi[0][j] = (q_in / (2.0*MU)) * ((pow(y[j], 3) / 3.0) - (pow(y[j], 2) / 2.0) * (y[J1] + y[NY]) + (y[j] * y[J1] * y[NY]));

    for(unsigned j = 0u; j <= NY; j++)  // BORDER C
        psi[NX][j] = (q_out / (2.0*MU)) * ((pow(y[j], 3) / 3.0) - (pow(y[j], 2) / 2.0 * y[NY])) + ((q_in * pow(y[J1], 2) * (-y[J1] + 3.0*y[NY])) / (12.0 * MU));

    for(unsigned i = 1u; i <= NX - 1u; i++)  // BORDER B
        psi[i][NY] = psi[0][NY];

    for(unsigned i = I1; i <= NX - 1u; i++)  // BORDER D
        psi[i][0] = psi[0][J1];

    for(unsigned j = 1u; j <= J1; j++)  // BORDER E
        psi[I1][j] = psi[0][J1];

    for(unsigned i = 1u; i <= I1; i++)  // BORDER F
        psi[i][J1] = psi[0][J1];
}

void bc_zeta(double zeta[NX+1][NY+1], const double *y, double psi[NX+1][NY+1], const double q_in, const double q_out)
{
    for(unsigned j = J1; j <= NY; j++)  // BORDER A
        zeta[0][j] = (q_in / (2.0*MU)) * ((2.0 * y[j]) - y[J1] - y[NY]);

    for(unsigned j = 0u; j <= NY; j++)  // BORDER C
        zeta[NX][j] = (q_out / (2.0*MU)) * ((2.0 * y[j]) - y[NY]);

    for(unsigned i = 1u; i <= NX - 1u; i++)  // BORDER B
        zeta[i][NY] = (2.0 / pow(DT, 2)) * (psi[i][NY-1] - psi[i][NY]);

    for(unsigned i = I1 + 1u; i <= NX - 1; i++)  // BORDER D
        zeta[i][0] = (2.0 / pow(DT, 2)) * (psi[i][1] - psi[i][0]);

    for(unsigned j = 1u; j <= J1 - 1u; j++)  // BORDER E
        zeta[I1][j] = (2.0 / pow(DT, 2)) * (psi[I1+1][j] - psi[I1][j]);

    for(unsigned i = 1u; i <= I1; i++)  // BORDER F
        zeta[i][J1] = (2.0 / pow(DT, 2)) * (psi[i][J1+1] - psi[i][J1]);

    zeta[I1][J1] = 0.5 * (zeta[I1-1][J1] + zeta[I1][J1-1]);  // POINT E/F
}

double bigGamma(double psi[NX+1][NY+1], double zeta[NX+1][NY+1])
{
    const unsigned J2 = J1 + 2;
    double gamma = 0.0;
    for(unsigned i = 1u; i <= NX - 1u; i++)
        gamma += (psi[i+1][J2] + psi[i-1][J2] + psi[i][J2+1] + psi[i][J2-1] - (4.0*psi[i][J2]) - DT*DT*zeta[i][J2]);
    return gamma;
}

void relaxation(double psi[NX+1][NY+1], double zeta[NX+1][NY+1], double v_x[NX+1][NY+1], double v_y[NX+1][NY+1], int q_in)
{
    for(unsigned i = 0u; i <= NX; i++)
    {
        for(unsigned j = 0u; j <= NY; j++)
        {
            psi[i][j] = 0;
            zeta[i][j] = 0;
            v_x[i][j] = 0;
            v_y[i][j] = 0;
        }
    }

    double y[NY+1];
    for(unsigned j = 0u; j <= NY; j++)
        y[j] = get_y(j);

    const double q_out = get_q_out(q_in, y);

    bc_psi(psi, y, q_in, q_out);
    
    double omega = 0.0;
    for(unsigned it = 1u; it <= IT_MAX; it++)
    {
        if(it < 2000u) { omega = 0.0; }
        else { omega = 1.0; }

        for(unsigned i = 1u; i <= NX - 1u; i++)
        {
            for(unsigned j = 1u; j <= NY - 1u; j++)
            {
                // A: i = 0, J1 < j < NY
                // B: i = *, j = NY
                // C: i = NX, j = *
                // D: I1 < i < NX, j = *
                // E: i = I1, 0 < j < J1
                // F: 0 < i < I1, j = J1

                if(i > I1 || j > J1)
                {
                    psi[i][j] = 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - DT*DT*zeta[i][j]);
                    zeta[i][j] = 0.25 * (zeta[i+1][j] + zeta[i-1][j] + zeta[i][j+1] + zeta[i][j-1])
                        - (omega * RHO / (16.0 * MU))
                        * ( ((psi[i][j+1] - psi[i][j-1]) * (zeta[i+1][j] - zeta[i-1][j]))
                        - ((psi[i+1][j] - psi[i-1][j]) * (zeta[i][j+1] - zeta[i][j-1])) );
                }
                v_x[i][j] = (psi[i][j+1] - psi[i][j])/DY;
                v_y[i][j] = -(psi[i+1][j] - psi[i][j])/DX;
            }
        }
        bc_zeta(zeta, y, psi, q_in, q_out);
        double gamma = bigGamma(psi, zeta);
    }
    
    const int N = snprintf(NULL, 0, "%d", q_in);
    const int N2 = snprintf(NULL, 0, "results_.txt");
    char filename[N+N2+1];
    sprintf(filename, "results_%d.txt", q_in);
    FILE *file = fopen(filename, "w");
    for(unsigned i = 0u; i <= NX; i++)
    {
        for(unsigned j = 0u; j <= NY; j++)
        {
            if( !(i > I1 || j > J1) )
                fprintf(file, "NaN NaN NaN NaN\n");
            else
                fprintf(file, "%lf %lf %lf %lf\n", psi[i][j], zeta[i][j], v_x[i][j], v_y[i][j]);
        }
    }
    fclose(file);
}
