#include <cmath>
#include <iostream>
#include <fstream>

extern "C" {
#include <mgmres.h>
}

constexpr static double dt = 0.1;
constexpr static double dx = dt;
constexpr static double dy = dt;
template <unsigned nx, unsigned ny> constexpr unsigned Nl = (nx+1) * (ny+1);

template <
    unsigned nx,
    unsigned ny,
    int eps_1 = 1,
    int eps_2 = 1,
    int V1 = 10,
    int V2 = -10,
    int V3 = 10,
    int V4 = -10,
    unsigned task = 0u>
static void sparse();

template <unsigned nx, unsigned ny>
constexpr static auto get_ij(const unsigned l);

template <
    unsigned nx,
    unsigned ny,
    int V1 = 10,
    int V2 = -10,
    int V3 = 10,
    int V4 = -10,
    unsigned task = 0u>
constexpr static void fillCSRDBC(
    double (&b)[Nl<nx, ny>],
    double (&a)[5u * Nl<nx, ny>],
    int (&ia)[Nl<nx, ny> + 1u],
    int (&ja)[5u * Nl<nx, ny>],
    const int (&eps)[Nl<nx, ny>]);

template <unsigned nx, unsigned ny>
constexpr static void fillEps_l(
    int (&eps_l)[Nl<nx, ny>],
    const int &eps_1,
    const int &eps_2);

template <typename T, unsigned nx, unsigned ny>
constexpr static void saveToFile(const T &arr, const unsigned N, const char *prefix, const char *suffix = "A");

template <typename T, unsigned nx, unsigned ny>
constexpr static void saveVToFile(const T &V, const char *filename);

template <unsigned nx, unsigned ny, unsigned task = 0u>
constexpr static auto get_rhos(const unsigned i, const unsigned j);

int main()
{
    // 3
    sparse<4, 4, 1, 1, 10, -10, 10, -10, 3u>();

    // 5: a, b, c
    sparse<50, 50>();
    sparse<100, 100>();
    sparse<200, 200>();

    // 6: a, b, c
    sparse<100, 100, 1, 1, 0, 0, 0, 0, 6u>();
    sparse<100, 100, 1, 2, 0, 0, 0, 0, 6u>();
    sparse<100, 100, 1, 10, 0, 0, 0, 0, 6u>();
    return 0;
}

template <unsigned nx, unsigned ny, int eps_1, int eps_2, int V1, int V2, int V3, int V4, unsigned task>
static void sparse()
{
    // decl
    double V[Nl<nx, ny>] = {0.0};
    double b[Nl<nx, ny>] = {0.0};
    double a[5u * Nl<nx, ny>] = {0.0};
    int ia[Nl<nx, ny> + 1u] = {0};
    int ja[5u * Nl<nx, ny>] = {0};
    int eps[Nl<nx, ny>] = {0};

    // init
    std::fill(std::begin(ia), std::end(ia), -1);
    fillEps_l<nx, ny>(eps, eps_1, eps_2);
    fillCSRDBC<nx, ny, V1, V2, V3, V4, task>( b, a, ia, ja, eps);

    if constexpr(task == 3u)
    {
        // save
        saveToFile<decltype(b), nx, ny>(b, Nl<nx, ny>, "init", "B");
        saveToFile<decltype(a), nx, ny>(a, (Nl<nx, ny>) * 5u, "init", "A");
    }
    else
    {
        // calc
        pmgmres_ilu_cr(Nl<nx, ny>, ia[Nl<nx, ny>], ia, ja, a, V, b, 500, 500, 1e-8, 1e-8);

        // save
        saveVToFile<decltype(V), nx, ny>(V, (
            "V_"+std::to_string(nx)+
            "_eps1_"+std::to_string(eps_1)+
            "_eps2_"+std::to_string(eps_2)+
            "_V_"+std::to_string(V1)+
            "_data.txt"
            ).c_str());
    }   
}

template <unsigned nx, unsigned ny>
constexpr static unsigned get_l(const unsigned i, const unsigned j)
{
    return i + j * (nx + 1);
}

template <unsigned nx, unsigned ny>
constexpr static auto get_ij(const unsigned l)
{
    return std::make_pair(l - std::floor(l / (nx + 1.0)) * (nx + 1.0), std::floor(l / (nx + 1.0)));
}

template <unsigned nx, unsigned ny, int V1, int V2, int V3, int V4, unsigned task>
constexpr static void fillCSRDBC(
    double (&b)[Nl<nx, ny>],
    double (&a)[5u * Nl<nx, ny>],
    int (&ia)[Nl<nx, ny> + 1u],
    int (&ja)[5u * Nl<nx, ny>],
    const int (&eps)[Nl<nx, ny>])
{
    int k = -1;
    for(unsigned l{0}; l < Nl<nx, ny>; l++)
    {
        int border = 0;
        double vb = 0.0;
        auto [i, j] = get_ij<nx, ny>(l);
        auto [rho1, rho2] = get_rhos<nx, ny, task>(i, j);
        if(i == 0)
        {
            border = 1;
            vb = V1;
        }
        if(j == ny)
        {
            border = 1;
            vb = V2;
        }
        if(i == nx)
        {
            border = 1.0;
            vb = V3;
        }
        if(j == 0)
        {
            border = 1;
            vb = V4;
        }

        b[l] = - (rho1 + rho2);
        if(border == 1)
        {
            b[l] = vb;
        }

        ia[l] = -1.0;
        if((static_cast<int>(l) - static_cast<int>(nx) - 1 >= 0) && (border == 0))
        {
            k++;
            if(ia[l] < 0)
            {
                ia[l] = k;
            }
            a[k] = eps[l]/(dt*dt);
            ja[k] = l - nx - 1;
        }
        if((static_cast<int>(l) - 1 >= 0) && (border == 0))
        {
            k++;
            if(ia[l] < 0)
            {
                ia[l] = k;
            }
            a[k] = eps[l]/(dt*dt);
            ja[k] = l - 1;
        }

        //
        k++;
        if(ia[l] < 0.0)
        {
            ia[l] = k;
        }
        if(border == 0)
        {
            a[k] = (-(2.0*eps[l]+eps[l+1]+eps[l+nx+1])/(dt*dt));
        }
        else
        {
            a[k] = 1.0;
        }
        ja[k] = l;
        //

        if((l < Nl<nx, ny>) && (border == 0))
        {
            k++;
            a[k] = eps[l+1]/(dt*dt);
            ja[k] = l+1;
        }
        if((l < Nl<nx, ny> - nx - 1) && (border == 0))
        {
            k++;
            a[k] = eps[l+nx+1]/(dt*dt);
            ja[k] = l + nx + 1;
        }
    }
    ia[Nl<nx, ny>] = k+1;
}

template <unsigned nx, unsigned ny>
constexpr static void fillEps_l(int (&eps_l)[Nl<nx, ny>], const int &eps_1, const int &eps_2)
{
    for(unsigned l{0}; l < Nl<nx, ny>; l++)
    {
        auto [i,j] = get_ij<nx, ny>(l);
        if(i <= (nx/2u))
        {
            eps_l[l] = eps_1;
        }
        else
        {
            eps_l[l] = eps_2;
        }
    }
}

template <typename T, unsigned nx, unsigned ny>
constexpr static void saveToFile(const T &arr, const unsigned N, const char *prefix, const char *suffix)
{
    std::ofstream file(std::string(prefix)+std::string(suffix)+".txt");
    for(unsigned l{0}; l < N; l++)
    {
        if(std::fabs(arr[l]) > 0)
        {
            auto [i, j] = get_ij<nx, ny>(l);
            file << l << " " << i << " " << j << " " << arr[l] << std::endl;
        }
    }
    file.close();
}

template <typename T, unsigned nx, unsigned ny>
constexpr static void saveVToFile(const T &V, const char *filename)
{
    std::ofstream file(filename);
    for(unsigned l{0}; l < Nl<nx, ny>; l++)
    {
        auto [i, j] = get_ij<nx, ny>(l);
        file << l << " " << i << " " << j << " " << V[l] << std::endl;
    }
    file.close();
}

template <unsigned nx, unsigned ny, unsigned task>
constexpr static auto get_rhos(const unsigned i, const unsigned j)
{
    if constexpr(task == 6u)
    {
        double x = dt * i;
        double y = dt * j;
        double xmax = dt * nx;
        double ymax = dt * ny;
        double sigma = xmax / 10.0;
        double rho1 = std::exp(-std::pow((x - 0.25*xmax)/sigma, 2) - std::pow((y - 0.5*ymax)/sigma, 2));
        double rho2 = -1.0 * std::exp(-std::pow((x - 0.75*xmax)/sigma, 2) - std::pow((y - 0.5*ymax)/sigma, 2));
        return std::make_pair(rho1, rho2);
    }
    else
    {
        return std::make_pair(0.0, 0.0);
    }
}
