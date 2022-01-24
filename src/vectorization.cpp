// ===========================================================
//
// vectorization.cpp: optimization with vectorization
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

#include "vectorization.h"
#include <Rdefines.h>
#include <R.h>


using namespace std;


// ========================================================================= //

/// SIMD version
extern "C" SEXP saige_simd_version()
{
#ifdef COREARRAY_HAVE_TARGET_CLONES
	const bool fmv = true;
	const bool avx512f = __builtin_cpu_supports("avx512f") != 0;
	const bool avx2 = __builtin_cpu_supports("avx2") != 0;
	const bool avx  = __builtin_cpu_supports("avx")  != 0;
	const bool sse3 = __builtin_cpu_supports("sse3") != 0;
	const bool sse2 = __builtin_cpu_supports("sse2") != 0;
#else
	const bool fmv = false;
#ifdef __AVX512F__
	const bool avx512f = true;
#else
	const bool avx512f = false;
#endif
#ifdef __AVX2__
	const bool avx2 = true;
#else
	const bool avx2 = false;
#endif
#ifdef __AVX__
	const bool avx = true;
#else
	const bool avx = false;
#endif
#ifdef __SSE3__
	const bool sse3 = true;
#else
	const bool sse3 = false;
#endif
#ifdef __SSE2__
	const bool sse2 = true;
#else
	const bool sse2 = false;
#endif
#endif
	char buffer[256], *p=buffer;
	if (avx512f)
	{
		strcpy(p, "AVX512F"); p += strlen(p);
	} else if (avx2)
	{
		strcpy(p, "AVX2"); p += strlen(p);
	} else if (avx)
	{
		strcpy(p, "AVX"); p += strlen(p);
	} else if (sse3)
	{
		strcpy(p, "SSE3"); p += strlen(p);
	} else if (sse2)
	{
		strcpy(p, "SSE2"); p += strlen(p);
	} else {
		strcpy(p, "basic x86"); p += strlen(p);
	}
	if (fmv)
	{
		strcpy(p, " (FMV)"); p += strlen(p);
	}
	return mkString(buffer);
}


// ========================================================================= //

namespace vectorization
{

/// get mean and sd
COREARRAY_TARGET_CLONES
	void f64_mean_sd(const double x[], size_t n, double &mean, double &sd)
{
	size_t m=0;
	double sum = 0, sum2 = 0;
	for (size_t i=0; i < n; i++)
	{
		const double v = x[i];
		if (isfinite(v))
			{ m++; sum += v; sum2 += v*v; }
	}
	mean = sd = R_NaN;
	if (m > 0)
	{
		mean = sum / m;
		if (m > 1) sd = sqrt((sum2 - sum*sum/m) / (m - 1));
	}
}

/// get max and min
void f64_maxmin(const double x[], size_t n, double &max, double &min)
{
	double vmax = -INFINITY, vmin = INFINITY;
	for (size_t i=0; i < n; i++)
	{
		const double v = x[i];
		if (isfinite(v))
		{
			if (v > vmax) vmax = v;
			if (v < vmin) vmin = v;
		}
	}
	if (!isfinite(vmax)) vmax = R_NaN;
	if (!isfinite(vmin)) vmin = R_NaN;
	max = vmax; min = vmin;
}

/// get max, min, median
void f64_medmaxmin(const double x[], size_t n, double &med, double &max, double &min)
{
	double vmax = -INFINITY, vmin = INFINITY;
	size_t num = 0;
	for (size_t i=0; i < n; i++)
	{
		const double v = x[i];
		if (isfinite(v))
		{
			if (v > vmax) vmax = v;
			if (v < vmin) vmin = v;
			num ++;
		}
	}
	if (!isfinite(vmax)) vmax = R_NaN;
	if (!isfinite(vmin)) vmin = R_NaN;
	max = vmax; min = vmin;
	// find median
	if (num > 0)
	{
		const size_t i1_med = (num - 1) / 2;
		const size_t i2_med = num / 2;
		double v1 = R_NaN, v2 = R_NaN;
		vmin = INFINITY; num = 0;
		for (size_t i=0; i < n; i++)
		{
			const double v = x[i];
			if (isfinite(v))
			{
				if (v <= vmin)
				{
					vmin = v;
					if (num == i1_med) v1 = v;
					if (num == i2_med) v2 = v;
					if (num > i2_med) break;
					num ++;
				}
			}
		}
		med = (v1 + v2) * 0.5;
	} else
		med = R_NaN;
}

/// return allele frequency and impute genotype using the mean
COREARRAY_TARGET_CLONES
	void f64_af_ac_impute(double *ds, size_t n, double &AF, double &AC, int &Num, int buf_idx[])
{
	double sum = 0;
	int num = 0, *pIdx = buf_idx;
	for (size_t i=0; i < n; i++)
	{
		if (isfinite(ds[i]))
			{ sum += ds[i]; num ++; }
		else
			*pIdx++ = i;
	}
	AF = (num > 0) ? (sum/(2*num)) : R_NaN;
	AC = sum; Num = num;
	if (num < (int)n)
	{
		double d = AF * 2;
		for (; buf_idx < pIdx; ) ds[*buf_idx++] = d;
	}
}


/// get the index of each nonzero value in x and return the number of nonzeros
COREARRAY_TARGET_CLONES size_t f64_nonzero_index(size_t n, const double *x, int *i)
{
	size_t n_i = 0;
	for (size_t j=0; j < n; j++)
		if (x[j] != 0) i[n_i++] = j;
	return n_i;
}


/// y[i] += x
COREARRAY_TARGET_CLONES MATH_OFAST void f64_add(size_t n, double x, double *y)
{
	for (size_t i=0; i < n; i++) y[i] += x;
}


/// y[i] += x[i]
COREARRAY_TARGET_CLONES MATH_OFAST void f64_add(size_t n, const double *x, double *y)
{
	for (size_t i=0; i < n; i++) y[i] += x[i];
}


/// y[i] = x - y[i]
COREARRAY_TARGET_CLONES MATH_OFAST void f64_sub(size_t n, double x, double *y)
{
	for (size_t i=0; i < n; i++) y[i] = x - y[i];
}


/// y[i] = x * y[i]
COREARRAY_TARGET_CLONES MATH_OFAST void f64_mul(size_t n, double x, double *y)
{
	for (size_t i=0; i < n; i++) y[i] *= x;
}


/// sum_i x[i]*y[i]
COREARRAY_TARGET_CLONES MATH_OFAST
	double f64_dot(size_t n, const double *x, const double *y)
{
	double sum = 0;
	for (size_t i=0; i < n; i++) sum += x[i] * y[i];
	return sum;
}


/// sum_i x[i]
COREARRAY_TARGET_CLONES MATH_OFAST double f64_sum(size_t n, const double *x)
{
	double sum = 0;
	for (size_t i=0; i < n; i++) sum += x[i];
	return sum;
}


/// x[i] = x[i] / sum_i x[i] (excluding not finite numbers)
COREARRAY_TARGET_CLONES void f64_normalize(size_t n, double *x)
{
	double sum = 0;
	for (size_t i=0; i < n; i++)
		if (isfinite(x[i])) sum += x[i];
	if (sum > 0)
	{
		sum = 1 / sum;
		for (size_t i=0; i < n; i++)
			if (isfinite(x[i])) x[i] *= sum;
	}
}


/// out1 = sum_i x[i]*y[i], out2 = sum_i y[i]*y[i]
COREARRAY_TARGET_CLONES MATH_OFAST
	void f64_dot_sp(size_t n, const double *x, const double *y, double &out1, double &out2)
{
	double sum1=0, sum2=0;
	for (size_t i=0; i < n; i++)
	{
		sum1 += x[i] * y[i];
		sum2 += y[i] * y[i];
	}
	out1 = sum1; out2 = sum2;
}


/// out1 = sum_i x1[i]*y[i], out2 = sum_i x2[i]*y[i]*y[i]
COREARRAY_TARGET_CLONES MATH_OFAST
	void f64_dot_sp2(size_t n, const double *x1, const double *x2, const double *y, double &out1, double &out2)
{
	double sum1=0, sum2=0;
	for (size_t i=0; i < n; i++)
	{
		sum1 += x1[i] * y[i];
		sum2 += x2[i] * y[i] * y[i];
	}
	out1 = sum1; out2 = sum2;
}


/// vec(p_m) = mat(x_{m*n}) * vec(y_n)
COREARRAY_TARGET_CLONES MATH_OFAST
	void f64_mul_mat_vec(size_t n, size_t m, const double *x, const double *y, double *p)
{
	memset(p, 0, sizeof(double)*m);
	switch (m)
	{
	case 1:
		for (size_t k=0; k < n; k++) p[0] += y[k] * x[k];
		break;
	case 2:
		for (size_t k=0; k < n; k++, x+=2)
		{
			const double alpha = y[k];
			if (alpha != 0)
				{ p[0] += alpha * x[0];  p[1] += alpha * x[1]; }
		}
		break;
	case 3:
		for (size_t k=0; k < n; k++, x+=3)
		{
			const double alpha = y[k];
			if (alpha != 0)
			{
				p[0] += alpha * x[0];  p[1] += alpha * x[1];
				p[2] += alpha * x[2];
			}
		}
		break;
	case 4:
		for (size_t k=0; k < n; k++, x+=4)
		{
			const double alpha = y[k];
			if (alpha != 0)
			{
				p[0] += alpha * x[0];  p[1] += alpha * x[1];
				p[2] += alpha * x[2];  p[3] += alpha * x[3];
			}
		}
		break;
	case 5:
		for (size_t k=0; k < n; k++, x+=5)
		{
			const double alpha = y[k];
			if (alpha != 0)
			{
				p[0] += alpha * x[0];  p[1] += alpha * x[1];
				p[2] += alpha * x[2];  p[3] += alpha * x[3];
				p[4] += alpha * x[4];
			}
		}
		break;
	case 6:
		for (size_t k=0; k < n; k++, x+=6)
		{
			const double alpha = y[k];
			if (alpha != 0)
			{
				p[0] += alpha * x[0];  p[1] += alpha * x[1];
				p[2] += alpha * x[2];  p[3] += alpha * x[3];
				p[4] += alpha * x[4];  p[5] += alpha * x[5];
			}
		}
		break;
	case 7:
		for (size_t k=0; k < n; k++, x+=7)
		{
			const double alpha = y[k];
			if (alpha != 0)
			{
				p[0] += alpha * x[0];  p[1] += alpha * x[1];
				p[2] += alpha * x[2];  p[3] += alpha * x[3];
				p[4] += alpha * x[4];  p[5] += alpha * x[5];
				p[6] += alpha * x[6];
			}
		}
		break;
	case 8:
		for (size_t k=0; k < n; k++, x+=8)
		{
			const double alpha = y[k];
			if (alpha != 0)
			{
				p[0] += alpha * x[0];  p[1] += alpha * x[1];
				p[2] += alpha * x[2];  p[3] += alpha * x[3];
				p[4] += alpha * x[4];  p[5] += alpha * x[5];
				p[6] += alpha * x[6];  p[7] += alpha * x[7];
			}
		}
		break;
	default:
		for (size_t k=0; k < n; k++, x+=m)
		{
			double alpha = y[k];
			if (alpha != 0) // if sparse vector
				for (size_t i=0; i < m; i++) p[i] += alpha * x[i];
		}
	}
}


/// vec(p_m) = mat(x_{m*n}) * vec(y_n), y is a sparse vector with indices
COREARRAY_TARGET_CLONES MATH_OFAST
	void f64_mul_mat_vec_sp(size_t n_idx, const int *idx, size_t m, const double *x, const double *y, double *p)
{
	memset(p, 0, sizeof(double)*m);
	switch (m)
	{
	case 1:
		for (size_t k=0; k < n_idx; k++)
		{
			const size_t i = idx[k];
			p[0] += y[i] * x[i];
		}
		break;
	case 2:
		for (size_t k=0; k < n_idx; k++)
		{
			const size_t i = idx[k];
			const double alpha = y[i], *xx = &x[2*i];
			p[0] += alpha * xx[0];  p[1] += alpha * xx[1];
		}
		break;
	case 3:
		for (size_t k=0; k < n_idx; k++)
		{
			const size_t i = idx[k];
			const double alpha = y[i], *xx = &x[3*i];
			p[0] += alpha * xx[0];  p[1] += alpha * xx[1];
			p[2] += alpha * xx[2];
		}
		break;
	case 4:
		for (size_t k=0; k < n_idx; k++)
		{
			const size_t i = idx[k];
			const double alpha = y[i], *xx = &x[4*i];
			p[0] += alpha * xx[0];  p[1] += alpha * xx[1];
			p[2] += alpha * xx[2];  p[3] += alpha * xx[3];
		}
		break;
	case 5:
		for (size_t k=0; k < n_idx; k++)
		{
			const size_t i = idx[k];
			const double alpha = y[i], *xx = &x[5*i];
			p[0] += alpha * xx[0];  p[1] += alpha * xx[1];
			p[2] += alpha * xx[2];  p[3] += alpha * xx[3];
			p[4] += alpha * xx[4];
		}
		break;
	case 6:
		for (size_t k=0; k < n_idx; k++)
		{
			const size_t i = idx[k];
			const double alpha = y[i], *xx = &x[6*i];
			p[0] += alpha * xx[0];  p[1] += alpha * xx[1];
			p[2] += alpha * xx[2];  p[3] += alpha * xx[3];
			p[4] += alpha * xx[4];  p[5] += alpha * xx[5];
		}
		break;
	case 7:
		for (size_t k=0; k < n_idx; k++)
		{
			const size_t i = idx[k];
			const double alpha = y[i], *xx = &x[7*i];
			p[0] += alpha * xx[0];  p[1] += alpha * xx[1];
			p[2] += alpha * xx[2];  p[3] += alpha * xx[3];
			p[4] += alpha * xx[4];  p[5] += alpha * xx[5];
			p[6] += alpha * xx[6];
		}
		break;
	case 8:
		for (size_t k=0; k < n_idx; k++)
		{
			const size_t i = idx[k];
			const double alpha = y[i], *xx = &x[8*i];
			p[0] += alpha * xx[0];  p[1] += alpha * xx[1];
			p[2] += alpha * xx[2];  p[3] += alpha * xx[3];
			p[4] += alpha * xx[4];  p[5] += alpha * xx[5];
			p[6] += alpha * xx[6];  p[7] += alpha * xx[7];
		}
		break;
	default:
		for (size_t k=0; k < n_idx; k++)
		{
			const size_t i = idx[k];
			const double alpha = y[i], *xx = &x[m*i];
			for (size_t j=0; j < m; j++) p[j] += alpha * xx[j];
		}
	}
}


/// vec(p_n) = t(mat(x_{m*n})) * vec(y_m), with a subset
COREARRAY_TARGET_CLONES MATH_OFAST
	void f64_mul_mat_vec_sub(size_t n, const int *idx, size_t m, const double *x, const double *y, double *p)
{
	for (size_t i=0; i < n; i++)
	{
		size_t k = idx[i];
		const double *xx = &x[m*k];
		double sum = 0;
		for (size_t j=0; j < m; j++) sum += y[j] * xx[j];
		p[i] = sum;
	}
}


/// vec(p_n) = vec(x_n) - t(mat(y_{m*n})) * vec(z_m)
COREARRAY_TARGET_CLONES MATH_OFAST
	void f64_sub_mul_mat_vec(size_t n, size_t m, const double *x, const double *y, const double *z, double *p)
{
	switch (m)
	{
	case 1:
		for (size_t i=0; i < n; i++) p[i] = x[i] - z[0]*y[i];
		break;
	case 2:
		for (size_t i=0; i < n; i++, y+=2)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1]);
		break;
	case 3:
		for (size_t i=0; i < n; i++, y+=3)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2]);
		break;
	case 4:
		for (size_t i=0; i < n; i++, y+=4)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2] + z[3]*y[3]);
		break;
	case 5:
		for (size_t i=0; i < n; i++, y+=5)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2] + z[3]*y[3] +
				z[4]*y[4]);
		break;
	case 6:
		for (size_t i=0; i < n; i++, y+=6)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2] + z[3]*y[3] +
				z[4]*y[4] + z[5]*y[5]);
		break;
	case 7:
		for (size_t i=0; i < n; i++, y+=7)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2] + z[3]*y[3] +
				z[4]*y[4] + z[5]*y[5] + z[6]*y[6]);
		break;
	case 8:
		for (size_t i=0; i < n; i++, y+=8)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2] + z[3]*y[3] +
				z[4]*y[4] + z[5]*y[5] + z[6]*y[6] + z[7]*y[7]);
		break;
	default:
		for (; n > 0; n--)
		{
			double sum = 0;
			for (size_t i=0; i < m; i++) sum += y[i] * z[i];
			y += m;
			*p++ = (*x++) - sum;
		}
	}
}


/// t(vec(y)) * mat(x) * vec(y)
COREARRAY_TARGET_CLONES MATH_OFAST
	double f64_sum_mat_vec(size_t n, const double *x, const double *y)
{
	double sum = 0;
	for (size_t i=0; i < n; i++)
	{
		const double *xx = &x[n*i], a = y[i];
		for (size_t j=0; j < n; j++)
			sum += a * y[j] * xx[j];
	}
	return sum;
}

}
