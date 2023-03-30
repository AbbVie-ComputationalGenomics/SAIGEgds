// ===========================================================
//
// vectorization.h: optimization with vectorization
//
// Copyright (C) 2019-2020    Xiuwen Zheng / AbbVie-ComputationalGenomics
//
// This file is part of SAIGEgds.
//
// SAIGEgds is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// SAIGEgds is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public
// License along with SAIGEgds.
// If not, see <http://www.gnu.org/licenses/>.

#ifndef COREARRAY_NO_COMPILER_OPTIMIZE
#if defined(__clang__) && !defined(__APPLE__)
    #pragma clang optimize on
#endif
#if defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=4))
    #pragma GCC optimize("O3")
    #define MATH_OFAST    __attribute__((optimize("Ofast")))
#endif
#endif

#ifndef MATH_OFAST
#   define MATH_OFAST
#endif


// Function multiversioning (requiring target_clones)
#if (defined(__GNUC__) && (__GNUC__ >= 6) && !defined(COREARRAY_NO_TARGET)) && !defined(__AVX512F__) && !defined(_WIN32) && !defined(__APPLE__) && !defined(__arm__) && !defined(__aarch64__)
#   if defined(__x86_64__) || defined(__i386__)
#       define COREARRAY_HAVE_TARGET
#       define COREARRAY_TARGET(opt)    __attribute__((target(opt)))
#       define COREARRAY_HAVE_TARGET_CLONES
#       define COREARRAY_TARGET_CLONES    \
            __attribute__((target_clones("avx512f","avx2","avx","sse3","sse2","default")))
#   endif
//#elif defined(__clang__)  // not support
//#   define COREARRAY_HAVE_TARGET
//#   define COREARRAY_TARGET(opt)           __attribute__((target(opt)))
//#   define COREARRAY_TARGET_CLONES(opt)    __attribute__((target_clones(opt)))
#else
#   define COREARRAY_TARGET(opt)
#   define COREARRAY_TARGET_CLONES
#endif


#include <stddef.h>


namespace vectorization
{
	/// get mean and sd
	void f64_mean_sd(const double x[], size_t n, double &mean, double &sd);
	/// get max and min
	void f64_maxmin(const double x[], size_t n, double &max, double &min);
	/// get max, min, median
	void f64_medmaxmin(const double x[], size_t n, double &med, double &max, double &min);

	/// return allele frequency and impute genotype using the mean
	void f64_af_ac_impute(double *ds, size_t n, double &AF, double &AC,
		int &Num, int buf_idx[]);
	/// get the index of each nonzero value in x and return the number of nonzeros
	size_t f64_nonzero_index(size_t n, const double *x, int *i);

	/// y[i] += x
	void f64_add(size_t n, double x, double *y);
	/// y[i] += x[i]
	void f64_add(size_t n, const double *x, double *y);
	/// y[i] = x - y[i]
	void f64_sub(size_t n, double x, double *y);
	/// y[i] *= x
	void f64_mul(size_t n, double x, double *y);
	/// sum_i x[i]*y[i]
	double f64_dot(size_t n, const double *x, const double *y);
	/// sum_i x[i]
	double f64_sum(size_t n, const double *x);
	/// x[i] = x[i] / sum_i x[i] (excluding not finite numbers)
	void f64_normalize(size_t n, double *x);

	/// out1 = sum_i x[i]*y[i], out2 = sum_i y[i]*y[i]
	void f64_dot_sp(size_t n, const double *x, const double *y,
		double &out1, double &out2);
	/// out1 = sum_i x1[i]*y[i], out2 = sum_i x2[i]*y[i]*y[i]
	void f64_dot_sp2(size_t n, const double *x1, const double *x2,
		const double *y, double &out1, double &out2);

	/// vec(p_m) = mat(x_{m*n}) * vec(y_n), y is a sparse vector
	void f64_mul_mat_vec(size_t n, size_t m, const double *x,
		const double *y, double *p);
	/// vec(p_m) = mat(x_{m*n}) * vec(y_n), y is a sparse vector with indices
	void f64_mul_mat_vec_sp(size_t n_idx, const int *idx, size_t m,
		const double *x, const double *y, double *p);
	/// vec(p_n) = t(mat(x_{m*n})) * vec(y_m), with a subset
	void f64_mul_mat_vec_sub(size_t n, const int *idx, size_t m,
		const double *x, const double *y, double *p);
	/// vec(p_n) = vec(x_n) - t(mat(y_{m*n})) * vec(z_m)
	void f64_sub_mul_mat_vec(size_t n, size_t m,
		const double *x, const double *y, const double *z, double *p);
	/// t(vec(y)) * mat(x) * vec(y)
	double f64_sum_mat_vec(size_t n, const double *x, const double *y);
}
