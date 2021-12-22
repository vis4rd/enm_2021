#include <fstream>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

// constants
constexpr static std::size_t NX = 40;
constexpr static std::size_t NY = 40;
constexpr static std::size_t N = ((NX + 1) * (NY + 1));
constexpr static std::size_t IT_MAX = 2000;
constexpr static double DELTA = 1;
constexpr static double D_T = 1;
constexpr static double T_A = 40;
constexpr static double T_B = 0;
constexpr static double T_C = 30;
constexpr static double T_D = 0;
constexpr static double K_B = 0.1;
constexpr static double K_D = 0.6;

// abbreviations
auto &gms = gsl_matrix_set;
auto &gmg = gsl_matrix_get;
auto &gvs = gsl_vector_set;
auto &gvg = gsl_vector_get;

// function declarations
auto get_l(std::size_t i, std::size_t j);
void heat_diffusion();
void cranck_nicolson(gsl_matrix *A, gsl_matrix *B, gsl_vector* c, gsl_vector* T);
double diff(gsl_vector* T, std::size_t l);

// main function
auto main() -> int
{
	heat_diffusion();
	return 0;
}

auto get_l(std::size_t i, std::size_t j)
{
	return i + j*(NX+1u);
}

void heat_diffusion()
{
	int signum;

	gsl_matrix *A = gsl_matrix_calloc(N, N);
	gsl_matrix *B = gsl_matrix_calloc(N, N);
	gsl_vector *c = gsl_vector_calloc(N);
	gsl_vector *d = gsl_vector_calloc(N);
    gsl_vector *T = gsl_vector_calloc(N);
    gsl_permutation* p = gsl_permutation_calloc(N);

	cranck_nicolson(A, B, c, T);
	gsl_linalg_LU_decomp(A, p, &signum);

	for(std::size_t it{0}; it <= IT_MAX; it++)
	{
		gsl_blas_dgemv(CblasNoTrans, 1.0, B, T, 0.0, d);  // B * T
        gsl_blas_daxpy(1.0, c, d);
        gsl_linalg_LU_solve(A, p, d, T);

        if (it == 99 
		|| it == 199
        || it == 499
        || it == 999 
        || it == 1999)
		{
			std::stringstream fpath; fpath << "../results/T_" << it+1 << ".txt";
			std::stringstream fpath2; fpath2 << "../results/grad2T_" << it+1 << ".txt";
			
			std::ofstream file1(fpath.str().c_str());
			std::ofstream file2(fpath2.str().c_str());

			for(std::size_t i{1}; i < NX; i++)
			{
				for(std::size_t j{1}; j < NY; j++)
				{
					auto l = get_l(i, j);
					file1 << gvg(T, l) << " ";
					file2 << diff(T, l) << " ";
				}
				file1 << "\n";
				file2 << "\n";
			}
			file1.close();
			file2.close();
        }
	}

	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_vector_free(c);
	gsl_vector_free(d);
	gsl_vector_free(T);
	gsl_permutation_free(p);
}

void cranck_nicolson(gsl_matrix *A, gsl_matrix *B, gsl_vector* c,  gsl_vector* T)
{
	auto val = D_T / (2.0*std::pow(DELTA, 2));
	for(std::size_t i{1}; i < NX; i++)
	{
		for(std::size_t j{1}; j < NY; j++)
		{
			auto l = get_l(i, j);

			gms(A, l, l-NX-1, val);
            gms(A, l, l-1, val);
            gms(A, l, l+1, val);
            gms(A, l, l+NX+1, val);

            gms(A, l, l, -4*val-1);

            gms(B, l, l-NX-1, -val);
            gms(B, l, l-1, -val);
            gms(B, l, l+1, -val);
            gms(B, l, l+NX+1, -val);

            gms(B, l, l, 4*val-1);
		}
	}

	// left, right border
	val = 1.0;
	for(std::size_t i{0}; i <= NX; i += NX)
	{
		for(std::size_t j{0}; j <= NY; j++)
		{
			auto l = get_l(i, j);
			gms(A, l, l, val);
            gms(B, l, l, val);
            gvs(c, l, 0.0);
		}
	}

	// upper border at n+1
	val = 1.0 / (K_B * DELTA);
	for(std::size_t i{1}; i < NX; i++)
	{
		auto l = get_l(i, NY);
		gms(A, l, l-NX-1, -val);
		gms(A, l, l, 1.0 + val);
		gvs(c, l, T_B);
		for (std::size_t j{0}; j < N; j++)
		{
			gms(B, l, j, 0.0);
		}
	}

	// lower border at n+1
	for(std::size_t i{1}; i < NX; i++)
	{
		auto l = get_l(i, 0);
		gms(A, l, l, 1+val);
		gms(A, l, l+NX+1, -val);
		gvs(c, l, T_D);
		for(std::size_t j{0}; j < N; j++)
		{
			gms(B, l, j, 0.0);
		}
	}

	// vector T
	for(std::size_t j{0}; j <= NY; j++)
	{
		auto l = get_l(0, j);  // left border
		gvs(T, l, T_A);
		
		l = get_l(NX, j);  // right border
		gvs(T, l, T_C);
		
		for (std::size_t i{1}; i < NX; i++)
		{
			l = get_l(i, j);  // the rest
			gvs(T, l, 0.0);
		}
	}
}

double diff(gsl_vector* T, std::size_t l)
{
    return ((gvg(T, l+1)    - 2.0*gvg(T, l) + gvg(T, l-1))    / std::pow(DELTA, 2))
		+  ((gvg(T, l+NX+1) - 2.0*gvg(T, l) + gvg(T, l-NX-1)) / std::pow(DELTA, 2));
}
