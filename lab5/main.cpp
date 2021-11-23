#include <cmath>
#include <iostream>

#include <FileUtils.hpp>

// Probably works only with GCC 10

constexpr static double dt = 0.2;
constexpr static double dx = dt;
constexpr static double dy = dt;
constexpr static unsigned nx = 128;
constexpr static unsigned ny = 128;
constexpr static double xmax = dx*nx;
constexpr static double ymax = dy*ny;
constexpr static double TOL = 1e-8;

constexpr static void multigrid();
constexpr static void setupBorderConditions(double (&V)[nx+1][ny+1]);
constexpr static void discretization(double (&V)[nx+1][ny+1], const unsigned K);
constexpr static void fineGrid(double (&V)[nx+1][ny+1], const unsigned K);
constexpr static double stop(const double (&V)[nx+1][ny+1], const unsigned K);

int main()
{
    multigrid();
    return 0;
}

constexpr static void multigrid()
{
    double V[nx+1][ny+1] = {0};
    setupBorderConditions(V);

    unsigned iter = 1;
    for(const auto K : {16u, 8u, 4u, 2u, 1u})
    {
        fu::removeFile("grid_k"+std::to_string(K)+"_s.txt");
        double S = 0;
        double Sp;
        do
        {
            discretization(V, K);
            
            Sp = S;
            S = stop(V, K);  // CALCULATING STOP CONDITION

            fu::impl::appendToFile(std::ofstream{"grid_k"+std::to_string(K)+"_s.txt", std::ios::app}, iter, S);
            iter++;
        }
        while(std::fabs((S - Sp)/Sp) >= TOL);

        fu::saveMatricesToFile("grid_k"+std::to_string(K)+"_v.txt", V);

        if(K >= 2)  // GRID FINING
        {
            fineGrid(V, K);
            setupBorderConditions(V);
        }
    }
}

constexpr static void setupBorderConditions(double (&V)[nx+1][ny+1])
{
    for(unsigned i{0}; i <= nx; i++)
    {
        V[0][i] = std::sin(M_PI * i*dy / ymax);
        V[i][ny] = -std::sin(2.0 * M_PI * i*dx / xmax);
        V[nx][i] = std::sin(M_PI * i*dy / ymax);
        V[i][0] = std::sin(2.0 * M_PI * i*dx / xmax);
    }
}

constexpr static void discretization(double (&V)[nx+1][ny+1], const unsigned K)
{
    for(auto i{K}; i <= nx - K; i+=K)
    {
        for(auto j{K}; j <= ny - K; j+=K)
        {
            V[i][j] = 0.25 * (V[i+K][j] + V[i-K][j] + V[i][j+K] + V[i][j-K]);
        }
    }
}

constexpr static void fineGrid(double (&V)[nx+1][ny+1], const unsigned K)
{
    for(unsigned i{0}; i <= nx - K; i+=K)
    {
        for(unsigned j{0}; j <= ny - K; j+=K)
        {
            V[i+K/2][j+K/2] = 0.25 * (V[i][j] + V[i+K][j] + V[i][j+K] + V[i+K][j+K]);
            V[i+K][j+K/2] = 0.5 * (V[i+K][j] + V[i+K][j+K]);
            V[i+K/2][j+K] = 0.5 * (V[i][j+K] + V[i+K][j+K]);
            V[i+K/2][j] = 0.5 * (V[i][j] + V[i+K][j]);
            V[i][j+K/2] = 0.5 * (V[i][j] + V[i][j+K]);
        }
    }
}

constexpr static double stop(const double (&V)[nx+1][ny+1], const unsigned K)
{
    double kdt = K*dx;
    double S = 0;
    for(unsigned i{0}; i <= nx-K; i+=K)
    {
        for(unsigned j{0}; j <= ny-K; j+=K)
        {
            S += ( ((kdt*kdt)/2.0) *
            (std::pow((V[i+K][j]-V[i][j])/(2*kdt) + (V[i+K][j+K] - V[i][j+K])/(2*kdt), 2) +
            std::pow((V[i][j+K]-V[i][j])/(2*kdt) + (V[i+K][j+K]-V[i+K][j])/(2*kdt), 2)) );
        }
    }
    return S;
}