#include <iostream>
#include <vector>
#include <complex>

const double pi = std::acos(-1);
const std::complex<double> j(0, 1);

std::vector<std::complex<double>> fft(const std::vector<double>& x)
{
    const int N = x.size();
    std::vector<std::complex<double>> F(N);

    // base case
    if(N == 2)
    {
        F[0] = x[0] + x[1];
        F[1] = x[0] - x[1];
    }
    else
    {
        // populate W_N^k for all k 0..N/2-1
        std::vector<std::complex<double>> W(N/2);
        const auto& two_pi_over_N = 2.*pi/N;
        for(int k = 0; k < N/2; ++k)
            W[k] = std::exp(two_pi_over_N*k*-j);

        // partion x into even and odd indices
        std::vector<double> evens(N/2);
        std::vector<double> odds(N/2);
        for(int n = 0; n < N/2; ++n)
        {
            evens[n] = x[2*n];
            odds[n] = x[2*n+1];
        }

        // divide
        const auto& G = fft(evens);
        const auto& H = fft(odds);

        // conquer
        for(int k = 0; k < N/2; ++k)
        {
            F[k] = G[k] + W[k]*H[k];
            F[k+N/2] = G[k] - W[k]*H[k];
        }
    }

    return F;
}


int main()
{
    std::vector<double> x{1., 2., 3., 4.};
    const auto& F = fft(x);
    for(const auto& e : x)
        std::cout << e << " ";
    std::cout << std::endl;
    for(const auto& f : F)
        std::cout << f << " ";
    std::cout << std::endl;
}
