#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#define active 1

static const double DELTA = 0.01;
static const unsigned NX = 400u;
static const unsigned NY = 90u;
static const unsigned I1 = 200u;
static const unsigned I2 = 210u;
static const unsigned J1 = 50u;
static const double SIGMA = 10 * DELTA;
static const double XA = 0.45;
static const double YA = 0.45;
static const double ITMAX = 10000;
static const double K[5] = { ITMAX/5.0, 2.0/5.0*ITMAX, 3.0/5.0*ITMAX, 4.0/5.0*ITMAX, ITMAX };
static const double M_PI = 3.1415926535897932;

double getX(const unsigned i);
double getY(const unsigned j);
double getDT(double vx[NX+1][NY+1], double vy[NX+1][NY+1]);
void crank_nicolson(double u[NX+1][NY+1], double up[NX+1][NY+1], double vx[NX+1][NY+1], double vy[NX+1][NY+1], const double D);
void speed_vector(double vx[NX+1][NY+1], double vy[NX+1][NY+1]);
void starting_condition(double u[NX+1][NY+1]);
void advection_diffussion(const double D);
double getC(double up[NX+1][NY+1]);
double avgX(double up[NX+1][NY+1]);
int doubleMatToFile(const char* filename, double mat[NX+1][NY+1]);
void debug_log(const char* output, ...);

int main(void)
{
    advection_diffussion(0.0);
    advection_diffussion(0.1);
}

double getX(const unsigned i) { return DELTA * i; }
double getY(const unsigned j) { return DELTA * j; }
double getDT(double vx[NX+1][NY+1], double vy[NX+1][NY+1])
{
    double max = 0.0;
    for (unsigned i = 0; i <= NX; i++)
        for (unsigned j = 0; j <= NY; j++)
            if(max < sqrt(pow(vx[i][j], 2) + pow(vy[i][j], 2)))
                max = sqrt(pow(vx[i][j], 2) + pow(vy[i][j], 2));
    return DELTA / (4.0 * max);
}

void crank_nicolson(
    double u[NX+1][NY+1],
    double up[NX+1][NY+1],
    double vx[NX+1][NY+1],
    double vy[NX+1][NY+1],
    const double D)
{
    double dt = getDT(vx, vy);
    memcpy(u, up, (NX+1)*(NY+1)*sizeof(double));
    for(unsigned micro = 1u; micro <= 20; micro++)
    {
        for(unsigned i = 0u; i <= NX; i++)
        {
            for(unsigned j = 1u; j<= NY - 1u; j++)
            {
                if( (i >= I1) && (i <= I2) && (j <= J1) )
                {
                    continue;
                }
                else if(i == 0u)
                {
                    u[i][j] = (1.0 / (1.0 + 2.0*D*dt/(DELTA*DELTA)))
                        * ( up[i][j]
                        - 0.5*dt*vx[i][j] 
                            * ( ((up[i+1][j] - up[NX][j]) / (2.0*DELTA))
                                + (u[i+1][j] - u[NX][j])/(2.0*DELTA) )
                        - 0.5*dt*vy[i][j]
                            * ( ((up[i][j+1] - up[i][j-1]) / (2.0*DELTA))
                                + (u[i][j+1] - u[i][j-1])/(2.0*DELTA) )
                        + 0.5*dt*D
                            * ( ((up[i+1][j] + up[NX][j] + up[i][j+1]+up[i][j-1]-4.0*up[i][j])/(DELTA*DELTA))
                                + (u[i+1][j]+u[NX][j]+u[i][j+1]+u[i][j-1])/(DELTA*DELTA) ) );
                }
                else if(i == NX)
                {
                    u[i][j] = (1.0 / (1.0 + 2.0*D*dt/(DELTA*DELTA)))
                        * ( up[i][j]
                        - 0.5*dt*vx[i][j] 
                            * ( ((up[0][j] - up[i-1][j]) / (2.0*DELTA))
                                + (u[0][j] - u[i-1][j])/(2.0*DELTA) )
                        - 0.5*dt*vy[i][j]
                            * ( ((up[i][j+1] - up[i][j-1]) / (2.0*DELTA))
                                + (u[i][j+1] - u[i][j-1])/(2.0*DELTA) )
                        + 0.5*dt*D
                            * ( ((up[0][j] + up[i-1][j] + up[i][j+1]+up[i][j-1]-4.0*up[i][j])/(DELTA*DELTA))
                                + (u[0][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])/(DELTA*DELTA) ) );
                }
                else
                {
                    u[i][j] = (1.0 / (1.0 + 2.0*D*dt/(DELTA*DELTA)))
                        * ( up[i][j]
                        - 0.5*dt*vx[i][j] 
                            * ( ((up[i+1][j] - up[i-1][j]) / (2.0*DELTA))
                                + (u[i+1][j] - u[i-1][j])/(2.0*DELTA) )
                        - 0.5*dt*vy[i][j]
                            * ( ((up[i][j+1] - up[i][j-1]) / (2.0*DELTA))
                                + (u[i][j+1] - u[i][j-1])/(2.0*DELTA) )
                        + 0.5*dt*D
                            * ( ((up[i+1][j] + up[i-1][j] + up[i][j+1]+up[i][j-1]-4.0*up[i][j])/(DELTA*DELTA))
                                + (u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])/(DELTA*DELTA) ) );
                }
            }
        }
    }
}

void speed_vector(double vx[NX+1][NY+1], double vy[NX+1][NY+1])
{
    double psi[NX+1][NY+1];
    unsigned temp = 0;
    FILE *file = fopen("../field_stream.data", "r");
    if(!file) { debug_log("Error opening the .txt file.\n"); exit(0); }
    for(unsigned i = 0u; i <= NX; i++)
    {
        for(unsigned j = 0u; j <= NY; j++)
        {
            if(fscanf(file, "%6u%8u%14lf", &temp, &temp, &psi[i][j]) == 0)
            {
                debug_log("Error reading .txt file.\n");
                exit(0);
            }
        }
    }
    fclose(file);

    for(unsigned i = 1u; i <= NX - 1u; i++)
    {
        for(unsigned j = 1u; j <= NY - 1u; j++)
        {
            vx[i][j] = (psi[i][j+1] - psi[i][j-1]) / (2.0*DELTA);
            vy[i][j] = -(psi[i+1][j] - psi[i-1][j]) / (2.0*DELTA);
        }
    }
    for(unsigned i = I1; i <= I2; i++)
    {
        for(unsigned j = 0u; j <= J1; j++)
        {
            vx[i][j] = 0.0;
            vy[i][j] = 0.0;
        }
    }
    for(unsigned i = 1u; i <= NX - 1u; i++)
    {
        vx[i][0] = 0.0;
        vy[i][NY] = 0.0;
    }
    for(unsigned j = 0u; j <= NY; j++)
    {
        vx[0][j] = vx[1][j];
        vx[NX][j] = vx[NX-1][j];
    }
    debug_log("Speed area calculated.\n");
}

void starting_condition(double u[NX+1][NY+1])
{
    for(unsigned i = 0u; i <= NX; i++)
    {
        for(unsigned j = 0u; j <= NY; j++)
        {
            u[i][j] = 1.0 / (2.0*M_PI*SIGMA*SIGMA)
                * exp(-(pow(getX(i) - XA, 2) + pow(getY(j) - YA, 2)) / (2.0*SIGMA*SIGMA));
        }
    }
    debug_log("Starting condition calculated.\n");
}

void advection_diffussion(const double D)
{
    debug_log("Starting calculation for D = %lf\n", D);
    double u0[NX+1][NY+1];
    double u1[NX+1][NY+1];
    double vx[NX+1][NY+1];
    double vy[NX+1][NY+1];

    speed_vector(vx, vy);
    starting_condition(u0);

    unsigned kter = 0;
    char file_c[32]; sprintf(file_c, "c_avgx_D_%3.1lf.txt", D);
    FILE *file = fopen(file_c, "w");
    for(unsigned iter = 1u; iter <= ITMAX; iter++)
    {
        crank_nicolson(u1, u0, vx, vy, D);
        memcpy(u0, u1, (NX+1)*(NY+1)*sizeof(double));

        fprintf(file, "%lf %lf\n", getC(u1), avgX(u1));
        if(iter == K[kter] && kter <= 5)
        {
            char filename[32];
            sprintf(filename, "u_t_%u_D_%3.1lf.txt", kter, D);
            doubleMatToFile(filename, u1);
            debug_log(" - %u iterations calculated, continuing...\n", iter);
            kter++;
        }
    }
    fclose(file);
    char file_vx[32]; sprintf(file_vx, "vx_D_%3.1lf.txt", D);
    char file_vy[32]; sprintf(file_vy, "vy_D_%3.1lf.txt", D);
    doubleMatToFile(file_vx, vx);
    doubleMatToFile(file_vy, vy);
    debug_log("Flow calculated.\n\n");
}

double getC(double up[NX+1][NY+1])
{
    double c = 0.0;
    for(unsigned i = 1u; i <= NX - 1u; i++)
    {
        for(unsigned j = 1u; j <= NY - 1; j++)
        {
            c += (up[i][j] * DELTA * DELTA);
        }
    }
    return c;
}

double avgX(double up[NX+1][NY+1])
{
    double x = 0.0;
    for(unsigned i = 1u; i <= NX - 1u; i++)
    {
        for(unsigned j = 1u; j <= NY - 1; j++)
        {
            x += getX(i) * up[i][j] * DELTA * DELTA;
        }
    }
    return x;
}

int doubleMatToFile(const char* filename, double mat[NX+1][NY+1])
{
    FILE* fout = fopen(filename, "w");
    if(!fout) { return 0; }

    for(unsigned i = 0u; i <= NX; i++)
    {
        for(unsigned j = 0u; j <= NY; j++)
        {
            if(fprintf(fout, "%lf\n", mat[i][j]) == 0)
            {
                fclose(fout);
                return 0;
            }
        }
    }

    fclose(fout);
    return 1;
}

void debug_log(const char* output, ...)
{
	if(active)
	{
		va_list ptr;
		va_start(ptr, output);
		vprintf(output, ptr);
		va_end(ptr);
	}
	return;
}
