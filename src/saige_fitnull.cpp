// ===========================================================
//
// saige_fitnull.cpp: C++ implementation of fitting the null model
//
// Copyright (C) 2019-2020    Xiuwen Zheng / AbbVie-ComputationalGenomics
//
// This file is part of SAIGEgds. It was created based on the original SAIGE
// C++ and R codes in the SAIGE package. Compared with the original SAIGE,
// all single-precision floating-point numbers are changed to double precision,
// and a more efficient algorithm with packed 2-bit and sparse genotypes is
// implemented to calculate the cross product of genetic relationship matrix.
//
// SAIGEgds is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as published
// by the Free Software Foundation.
//
// SAIGEgds is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with SAIGEgds.
// If not, see <http://www.gnu.org/licenses/>.

#include "vectorization.h"
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/parallel_for.h>
#include <vector>
#include <algorithm>


using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
using namespace vectorization;


// ========================================================================= //
// Define Intel TBB macro with a thread index starting from 0
// requiring C++11

#if RCPP_PARALLEL_USE_TBB

#define PARALLEL_HEAD(SIZE, balancing)    \
	tbb::parallel_for(  \
		balancing ? tbb::blocked_range<size_t>(0, SIZE) :  \
		tbb::blocked_range<size_t>(0, SIZE, SIZE/NumThreads + (SIZE % NumThreads ? 1:0)),  \
		[&](const tbb::blocked_range<size_t> &r)  \
	{  \
		const int th_idx = tbb::this_task_arena::current_thread_index();  \
		if (th_idx < 0 || th_idx >= NumThreads)  \
			throw std::invalid_argument( \
				"Invalid tbb::this_task_arena::current_thread_index()!");

#define PARALLEL_FOR(i, SIZE, balancing)    \
		PARALLEL_HEAD(SIZE, balancing)    \
		for (size_t i=r.begin(); i < r.end(); i++)

#define PARALLEL_RANGE(st, ed, SIZE, balancing)    \
		PARALLEL_HEAD(SIZE, balancing)    \
		const size_t st = r.begin(), ed = r.end();

#define PARALLEL_END    });

#else

#define PARALLEL_FOR(i, SIZE, balancing)    \
	{  \
		const int th_idx = 0;  \
		for (size_t i=0; i < SIZE; i++)
#define PARALLEL_RANGE(st, ed, SIZE, balancing)    \
	{  \
		const int th_idx = 0;  \
		const size_t st = 0, ed = SIZE;
#define PARALLEL_END    }

#endif



// ========================================================================= //
// R functions

inline static NumericVector SEXP_VEC(const dvec &x)
{
	return NumericVector(x.begin(), x.end());
}

inline static void set_seed(unsigned int seed)
{
	Environment base_env("package:base");
	Function set_seed_r = base_env["set.seed"];
	set_seed_r(seed);
}


// ========================================================================= //
// Store 2-bit packed SNP genotypes

#define BYTE    unsigned char

static int NumThreads = 0;  //< the number of threads

static BYTE *Geno_PackedRaw = NULL;    //< the pointer to the 2-bit packed genotypes
static size_t Geno_NumSamp = 0;        //< the number of samples
static size_t Geno_PackedNumSamp = 0;  //< the number of bytes for packed samples
static size_t Geno_NumVariant = 0;     //< the number of variants

static double *buf_std_geno = NULL;   //< a 4-by-n_variant look-up matrix
static double *buf_diag_grm = NULL;   //< n_samp-length, sigma_i = sum_j adj.g[i,j]^2
static double *buf_crossprod = NULL;  //< nThread-by-n_samp matrix

// Internal 2-bit genotype lookup tables
static bool lookup_has_init = false;
static BYTE num_valid[256], num_sum[256];

static void init_lookup_table()
{
	if (!lookup_has_init)
	{
		for (int i=0; i < 256; i++)
		{
			int b0 = i & 0x03;
			int b1 = (i >> 2) & 0x03;
			int b2 = (i >> 4) & 0x03;
			int b3 = (i >> 6) & 0x03;
			num_valid[i] = (b0<3 ? 1:0) + (b1<3 ? 1:0) + (b2<3 ? 1:0) + (b3<3 ? 1:0);
			num_sum[i] = (b0<3 ? b0:0) + (b1<3 ? b1:0) + (b2<3 ? b2:0) + (b3<3 ? b3:0);
		}
		lookup_has_init = true;
	}
}

inline static double sq(double v) { return v*v; }
inline static double ds_nan(BYTE g) { return (g < 3) ? g : R_NaN; }


/// Store SNP genotype in the 2-bit packed format
RcppExport SEXP saige_store_2b_geno(SEXP rawgeno, SEXP num_samp,
	SEXP r_buf_geno, SEXP r_buf_sigma, SEXP r_buf_crossprod)
{
BEGIN_RCPP

	// initialize the basic genotype variables
	RawMatrix RawGeno(rawgeno);
	Geno_PackedRaw = (BYTE*)&RawGeno[0];
	Geno_NumSamp = Rf_asInteger(num_samp);
	Geno_PackedNumSamp = RawGeno.nrow();
	Geno_NumVariant = RawGeno.ncol();

	// set the buffer for get_crossprod_b_grm()
	NumericMatrix mat(r_buf_crossprod);
	buf_crossprod = REAL(r_buf_crossprod);
	NumThreads = mat.ncol();
	if (NumThreads > (int)Geno_NumSamp) NumThreads = Geno_NumSamp;
	if (NumThreads > (int)Geno_NumVariant) NumThreads = Geno_NumVariant;
	if (NumThreads < 1) NumThreads = 1;

	// build the look-up table of standardized genotypes
	init_lookup_table();
	buf_std_geno = REAL(r_buf_geno);
	PARALLEL_FOR(i, Geno_NumVariant, true)
	{
		BYTE *g = Geno_PackedRaw + Geno_PackedNumSamp*i;
		// calculate allele frequency
		int n_valid=0, sum=0;
		for (size_t j=0; j < Geno_PackedNumSamp; j++)
		{
			n_valid += num_valid[g[j]];
			sum += num_sum[g[j]];
		}
		double af = double(sum) / (2*n_valid);
		double inv = 1 / sqrt(2*af*(1-af));
		if (!R_FINITE(af) || !R_FINITE(inv))
			af = inv = 0;
		double *p = &buf_std_geno[4*i];
		p[0] = (0 - 2*af) * inv; p[1] = (1 - 2*af) * inv;
		p[2] = (2 - 2*af) * inv; p[3] = 0;
	}
	PARALLEL_END

	// calculate diag(grm)
	buf_diag_grm = REAL(r_buf_sigma);
	memset(buf_diag_grm, 0, sizeof(double)*Geno_NumSamp);
	for (size_t i=0; i < Geno_NumVariant; i++)
	{
		BYTE *g = Geno_PackedRaw + Geno_PackedNumSamp*i;
		const double *base = buf_std_geno + 4*i;
		size_t n = Geno_NumSamp;
		double *p = buf_diag_grm;
		for (; n >= 4; n-=4, p+=4)
		{
			BYTE gg = *g++;
			p[0] += sq(base[gg & 0x03]);
			p[1] += sq(base[(gg >> 2) & 0x03]);
			p[2] += sq(base[(gg >> 4) & 0x03]);
			p[3] += sq(base[gg >> 6]);
		}
		for (BYTE gg = (n>0 ? *g : 0); n > 0; n--)
		{
			(*p++) += sq(base[gg & 0x03]);
			gg >>= 2;
		}
	}
	f64_mul(Geno_NumSamp, 1.0 / Geno_NumVariant, buf_diag_grm);

END_RCPP
}


// ========================================================================= //
// Store sparse structure of SNP genotypes

/// List vector of sparse structure of genotypes, with a format:
/// integer vector: n1, n2, n3, n1-length int vector, n2-length, n3-length
static SEXP Geno_Sparse = NULL;

/// Get sparse structure of genotypes
RcppExport SEXP saige_get_sparse(SEXP geno, SEXP buffer)
{
BEGIN_RCPP
	const size_t num = Rf_length(geno);
	BYTE *gs = (BYTE*)RAW(geno);
	// determine whether need to flip
	int n=0, sum=0;
	for (size_t i=0; i < num; i++)
		if (gs[i] < 3) { sum += gs[i]; n++; }
	if (sum > n)
	{
		// flip
		for (size_t i=0; i < num; i++)
			if (gs[i] < 3) gs[i] = 2 - gs[i];
	}
	// get sparse structure
	int *base = INTEGER(buffer), *p = base+3;
	int &n1 = base[0], &n2 = base[1], &n3 = base[2];
	n1 = n2 = n3 = 0;
	// genotype: 1
	for (size_t i=0; i < num; i++)
		if (gs[i]==1) { *p++ = i; n1++; }
	// genotype: 2
	for (size_t i=0; i < num; i++)
		if (gs[i]==2) { *p++ = i; n2++; }
	// genotype: missing
	for (size_t i=0; i < num; i++)
		if (gs[i]==3) { *p++ = i; n3++; }
	// output
	return IntegerVector(base, p);
END_RCPP
}


/// Initialize sparse structure of genotypes
RcppExport SEXP saige_store_sp_geno(SEXP sp_geno_list, SEXP num_samp,
	SEXP r_buf_geno, SEXP r_buf_sigma, SEXP r_buf_crossprod)
{
BEGIN_RCPP

	// initialize the basic genotype variables
	Geno_PackedRaw = NULL;
	Geno_Sparse = sp_geno_list;
	Geno_NumSamp = Rf_asInteger(num_samp);
	Geno_NumVariant = Rf_length(sp_geno_list);

	// set the buffer for get_crossprod_b_grm()
	NumericMatrix mat(r_buf_crossprod);
	buf_crossprod = REAL(r_buf_crossprod);
	NumThreads = mat.ncol();
	if (NumThreads > (int)Geno_NumSamp) NumThreads = Geno_NumSamp;
	if (NumThreads > (int)Geno_NumVariant) NumThreads = Geno_NumVariant;
	if (NumThreads < 1) NumThreads = 1;

	// build the look-up table of standardized genotypes
	buf_std_geno = REAL(r_buf_geno);
	PARALLEL_FOR(i, Geno_NumVariant, true)
	{
		int *pg = INTEGER(VECTOR_ELT(Geno_Sparse, i));
		int n_valid = Geno_NumSamp - pg[2];
		int sum = pg[0] + 2*pg[1];
		double af = double(sum) / (2*n_valid);
		double inv = 1 / sqrt(2*af*(1-af));
		if (!R_FINITE(af) || !R_FINITE(inv))
			af = inv = 0;
		double *p = &buf_std_geno[4*i];
		p[0] = (0 - 2*af) * inv; p[1] = (1 - 2*af) * inv;
		p[2] = (2 - 2*af) * inv; p[3] = 0;
		p[1] -= p[0]; p[2] -= p[0]; p[3] -= p[0]; // adjustment
	}
	PARALLEL_END

	// calculate diag(grm)
	buf_diag_grm = REAL(r_buf_sigma);
	memset(buf_diag_grm, 0, sizeof(double)*Geno_NumSamp);
	double adj_g0 = 0, v;
	for (size_t i=0; i < Geno_NumVariant; i++)
	{
		const double *p = &buf_std_geno[4*i];
		const int *pg = INTEGER(VECTOR_ELT(Geno_Sparse, i));
		const int *ii = pg + 3;
		// g0^2
		double g0_2 = sq(p[0]); adj_g0 += g0_2;
		// g1^2
		v = sq(p[1] + p[0]) - g0_2;
		for (int k=0; k < pg[0]; k++) buf_diag_grm[*ii++] += v;
		// g2^2
		v = sq(p[2] + p[0]) - g0_2;
		for (int k=0; k < pg[1]; k++) buf_diag_grm[*ii++] += v;
		// g3^2
		v = sq(p[3] + p[0]) - g0_2;
		for (int k=0; k < pg[2]; k++) buf_diag_grm[*ii++] += v;
	}
	f64_add(Geno_NumSamp, adj_g0, buf_diag_grm);
	f64_mul(Geno_NumSamp, 1.0 / Geno_NumVariant, buf_diag_grm);

END_RCPP
}


// ========================================================================= //

/// Get a dosage vector for a specific SNP
static void get_geno_ds(int snp_idx, dvec &ds)
{
	ds.resize(Geno_NumSamp);
	if (!Geno_PackedRaw)
	{
		const int *pg = INTEGER(VECTOR_ELT(Geno_Sparse, snp_idx));
		const int *i = pg + 3;
		// g0
		memset(&ds[0], 0, sizeof(double)*Geno_NumSamp);
		// g1
		for (int k=0; k < pg[0]; k++) ds[*i++] = 1;
		// g2
		for (int k=0; k < pg[1]; k++) ds[*i++] = 2;
		// g3 -- missing
		for (int k=0; k < pg[2]; k++) ds[*i++] = R_NaN;
	} else {
		const BYTE *g = Geno_PackedRaw + Geno_PackedNumSamp * snp_idx;
		double *p = &ds[0];
		size_t n = Geno_NumSamp;
		for (; n >= 4; n-=4, p+=4)
		{
			BYTE gg = *g++;
			p[0] = ds_nan(gg & 0x03);
			p[1] = ds_nan((gg >> 2) & 0x03);
			p[2] = ds_nan((gg >> 4) & 0x03);
			p[3] = ds_nan(gg >> 6);
		}
		for (BYTE gg = (n>0 ? *g : 0); n > 0; n--)
		{
			*p++ = ds_nan(gg & 0x03);
			gg >>= 2;
		}
	}
}


// ========================================================================= //

/// Cross-product of standardized genotypes (G) and a numeric vector
/// Input: b (n_samp-length)
/// Output: out_b (n_samp-length) = GRM * b = G' G b
static COREARRAY_TARGET_CLONES MATH_OFAST
	void get_crossprod_b_grm(const dcolvec &b, dvec &out_b)
{
	// initialize
	memset(buf_crossprod, 0, sizeof(double)*Geno_NumSamp*NumThreads);
	double sum_b = f64_sum(Geno_NumSamp, &b[0]);
	dvec sum_cp_g0(NumThreads);
	sum_cp_g0.zeros();

	// crossprod with b
	if (!Geno_PackedRaw)
	{
		// sparse genotypes
		PARALLEL_FOR(i, Geno_NumVariant, true)
		{
			const double *p = &buf_std_geno[4*i], *pb = &b[0];
			const int *pg = INTEGER(VECTOR_ELT(Geno_Sparse, i)), *ii = pg + 3;
			// g0 * b
			double dot = sum_b * p[0];
			// g1 * b
			for (int k=0; k < pg[0]; k++) dot += p[1] * pb[*ii++];
			// g2 * b
			for (int k=0; k < pg[1]; k++) dot += p[2] * pb[*ii++];
			// g3 * b
			for (int k=0; k < pg[2]; k++) dot += p[3] * pb[*ii++];

			// update buf_crossprod += dot .* std.geno
			double v, *pbb = buf_crossprod + Geno_NumSamp * th_idx;
			ii = pg + 3;
			// g0 * dot
			sum_cp_g0[th_idx] += dot * p[0];
			// g1 * dot
			v = dot * p[1];
			for (int k=0; k < pg[0]; k++) pbb[*ii++] += v;
			// g2 * dot
			v = dot * p[2];
			for (int k=0; k < pg[1]; k++) pbb[*ii++] += v;
			// g3 * dot
			v = dot * p[3];
			for (int k=0; k < pg[2]; k++) pbb[*ii++] += v;
		}
		PARALLEL_END
	} else {
		// dense packed genotypes
		PARALLEL_FOR(i, Geno_NumVariant, true)
		{
			const BYTE *g = Geno_PackedRaw + Geno_PackedNumSamp*i;
			const double *base = buf_std_geno + 4*i, *pb = &b[0];

			// get dot = sum(std.geno .* b)
			double dot = 0;
			size_t n = Geno_NumSamp;
			for (; n >= 4; n-=4, pb+=4)
			{
				BYTE gg = *g++;
				dot += base[gg & 0x03] * pb[0] + base[(gg >> 2) & 0x03] * pb[1] +
					base[(gg >> 4) & 0x03] * pb[2] + base[gg >> 6] * pb[3];
			}
			for (BYTE gg = (n>0 ? *g : 0); n > 0; n--)
			{
				dot += base[gg & 0x03] * (*pb++);
				gg >>= 2;
			}

			// update buf_crossprod += dot .* std.geno
			double *pbb = buf_crossprod + Geno_NumSamp * th_idx;
			g = Geno_PackedRaw + Geno_PackedNumSamp*i;
			n = Geno_NumSamp;
			for (; n >= 4; n-=4, pbb+=4)
			{
				BYTE gg = *g++;
				pbb[0] += dot * base[gg & 0x03];
				pbb[1] += dot * base[(gg >> 2) & 0x03];
				pbb[2] += dot * base[(gg >> 4) & 0x03];
				pbb[3] += dot * base[gg >> 6];
			}
			for (BYTE gg = (n>0 ? *g : 0); n > 0; n--)
			{
				(*pbb++) += dot * base[gg & 0x03];
				gg >>= 2;
			}
		}
		PARALLEL_END
	}

	// normalize out_b
	out_b.resize(Geno_NumSamp);
	double sum_g0 = sum(sum_cp_g0);
	PARALLEL_RANGE(st, ed, Geno_NumSamp, false)
	{
		size_t len = ed - st;
		const double *s = buf_crossprod + st;
		double *p = &out_b[st];
		memset(p, 0, sizeof(double)*len);
		for (int i=0; i < NumThreads; i++, s += Geno_NumSamp)
			f64_add(len, s, p);
		if (!Geno_PackedRaw)
			f64_add(len, sum_g0, p);
		f64_mul(len, 1.0/Geno_NumVariant, p);
	}
	PARALLEL_END
}

  
/// Diagonal Sigma = tau[0] * diag(1/W) + tau[1] * diag(GRM)
/// Input: w, tau
/// Output: out_sigma
static COREARRAY_TARGET_CLONES
	void get_diag_sigma(const dvec& w, const dvec& tau, dvec &out_sigma)
{
	out_sigma.resize(Geno_NumSamp);
	PARALLEL_RANGE(st, ed, Geno_NumSamp, false)
	{
		const double tau0 = tau[0], tau1 = tau[1];
		double *out = &out_sigma[0];
		for (size_t i=st; i < ed; i++)
		{
			double v = tau0 / w[i] + tau1 * buf_diag_grm[i];
			if (v < 1e-4) v = 1e-4;
			out[i] = v;
		}
	}
	PARALLEL_END
}


/// Diagonal Sigma = tau[0] * b * diag(1/W) + tau[1] * diag(GRM, b)
/// Input: w, tau
/// Output: out_sigma
static COREARRAY_TARGET_CLONES
	dvec get_crossprod(const dcolvec &b, const dvec& w, const dvec& tau)
{
	const double tau0 = tau[0], tau1 = tau[1];
	if (tau1 == 0)
	{
		return(tau0 * (b % (1/w)));
	} else {
		dvec out_b;
		get_crossprod_b_grm(b, out_b);
		return(tau0 * (b % (1/w)) + tau1 * out_b);
	}
}


/// PCG algorithm for diagonal of Sigma = tau[0] * diag(1/W) + tau[1] * GRM
/// Input: w, tau, b, maxiterPCG, tolPCG
static COREARRAY_TARGET_CLONES
	dvec PCG_diag_sigma(const dvec &w, const dvec &tau, const dvec &b,
		int maxiterPCG, double tolPCG)
{
	dvec r = b, r1, minv;
	get_diag_sigma(w, tau, minv);
	minv = 1 / minv;

	dvec z = minv % r, z1;
	dvec p = z;
	dvec x(Geno_NumSamp);
	x.zeros();

	int iter = 0;
	while ((iter < maxiterPCG) && (sum(r % r) > tolPCG))
	{
		iter = iter + 1;
		dvec Ap = get_crossprod(p, w, tau);
		double a = sum(r % z) / sum(p % Ap);
		x += a * p;
		r1 = r - a * Ap;
		z1 = minv % r1;

		double bet = sum(z1 % r1) / sum(z % r);
		p = z1 + bet*p;
		z = z1;
		r = r1;
	}

	if (iter >= maxiterPCG)
		Rprintf("PCG does not converge (may need to increase 'maxiter').\n");

	return(x);
}


/// Calculate the coefficient of variation for mean of a vector
inline static double calcCV(const dvec &x)
{
	double x_mean = mean(x);
	double x_sd = stddev(x);
	return(x_sd / (x_mean * int(x.n_elem)));
}


/// Calculate the trace of matrix for binary outcomes
static double get_trace(const dmat &Sigma_iX, const dmat& X, const dvec& w,
	const dvec& tau, const dmat& cov, int nrun, int maxiterPCG, double tolPCG,
	double traceCVcutoff, int seed)
{
	set_seed(seed);
	dmat Sigma_iXt = Sigma_iX.t();
	dvec Sigma_iu;  
	dcolvec Pu;
	dvec Au, u;

	int nrunStart = 0;
	int nrunEnd = nrun;
	double traceCV = traceCVcutoff + 0.1;
	dvec buf(nrun);
	buf.zeros();

	// Hutchinson's randomized trace estimator
	while (traceCV > traceCVcutoff)
	{
		for (int i=nrunStart; i < nrunEnd; i++)
		{
			// buf[i] = Pu Au = (u' P) (G' G u) = u' P G' G u
			u = 2 * as<dvec>(rbinom(Geno_NumSamp, 1, 0.5)) - 1;
			Sigma_iu = PCG_diag_sigma(w, tau, u, maxiterPCG, tolPCG);
			Pu = Sigma_iu - Sigma_iX * (cov *  (Sigma_iXt * u));
			get_crossprod_b_grm(u, Au);
			buf[i] = dot(Au, Pu);
		}
		traceCV = calcCV(buf);
		if (traceCV > traceCVcutoff)
		{
			nrunStart = nrunEnd;
			nrunEnd = nrunEnd + 10;
			buf.resize(nrunEnd);
			Rprintf("CV for trace random estimator using %d runs is %g > %g\n",
				nrun, traceCV, traceCVcutoff);
			Rprintf("try %d runs ...\n", nrunEnd);
		}
	}

	return(mean(buf));
}


/// Calculate the trace of matrix for quantitative outcomes
static void get_trace_q(const dmat &Sigma_iX, const dmat& X, const dvec& w,
	const dvec& tau, const dmat& cov, int nrun, int maxiterPCG, double tolPCG,
	double traceCVcutoff, int seed, double &outTrace0, double &outTrace1)
{
	set_seed(seed);
	dmat Sigma_iXt = Sigma_iX.t();
	dvec Sigma_iu;  
	dcolvec Pu;
	dvec Au, u;

	int nrunStart = 0;
	int nrunEnd = nrun;
	double traceCV, traceCV0;
	traceCV = traceCV0 = traceCVcutoff + 0.1;
	dvec buf(nrun), buf0(nrun);
	buf.zeros(); buf0.zeros();

	// Hutchinson's randomized trace estimator
	while ((traceCV > traceCVcutoff) || (traceCV0 > traceCVcutoff))
	{
		for (int i=nrunStart; i < nrunEnd; i++)
		{
			// buf[i] = Pu Au = (u' P) (G' G u) = u' P G' G u
			u = 2 * as<dvec>(rbinom(Geno_NumSamp, 1, 0.5)) - 1;
			Sigma_iu = PCG_diag_sigma(w, tau, u, maxiterPCG, tolPCG);
			Pu = Sigma_iu - Sigma_iX * (cov *  (Sigma_iXt * u));
			get_crossprod_b_grm(u, Au);
			buf[i]  = dot(Au, Pu);
			buf0[i] = dot(u, Pu);
		}
		traceCV  = calcCV(buf);
		traceCV0 = calcCV(buf0);
		if ((traceCV > traceCVcutoff) || (traceCV0 > traceCVcutoff))
		{
			nrunStart = nrunEnd;
			nrunEnd = nrunEnd + 10;
			buf.resize(nrunEnd);
			buf0.resize(nrunEnd);
			Rprintf("CV for trace random estimator using %d runs is %g > %g\n",
				nrun, traceCV, traceCVcutoff);
			Rprintf("try %d runs ...\n", nrunEnd);
		}
	}

	outTrace0 = mean(buf0);
	outTrace1 = mean(buf);
}


/// matrix inverse
inline static dmat mat_inv(const dmat &m)
{
	dmat rv, xs = symmatu(m);
	if (!auxlib::inv_sympd(rv, xs))
	{
		// xs is singular or not positive definite (possibly due to rounding error)
		// try inv(), if inv() still fails throw an exception and stop fitting
		Rprintf("Warning: arma::inv_sympd(), matrix is singular or not positive definite, use arma::inv() instead.\n");
		rv = inv(xs);
	}
	return rv;
}


/// Calculate fixed and random effect coefficients given Y, X, w, tau
/// Input:  Y, X, w, tau, maxiterPCG, tolPCG
/// Output: Sigma_iY, Sigma_iX, cov, alpha, eta
static void get_coeff_w(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &tau, int maxiterPCG, double tolPCG,
	dvec &Sigma_iY, dmat &Sigma_iX, dmat &cov, dvec &alpha, dvec &eta)
{
	int n_col_X = X.n_cols;
	Sigma_iY = PCG_diag_sigma(w, tau, Y, maxiterPCG, tolPCG);
	// Sigma_iX = (X' Sigma^-1)'
	Sigma_iX.resize(Geno_NumSamp, n_col_X);
	dvec xv_i;
	for(int i = 0; i < n_col_X; i++)
	{
		xv_i = X.col(i);
		Sigma_iX.col(i) = PCG_diag_sigma(w, tau, xv_i, maxiterPCG, tolPCG);
	}
	// cov = (X' Sigma^-1 X)^-1
	cov = mat_inv(X.t() * Sigma_iX);
	// alpha = (X' Sigma^-1 X)^-1 X' Sigma^-1 Y
	alpha = cov * (Sigma_iX.t() * Y);
	eta = Y - tau[0] * (Sigma_iY - Sigma_iX * alpha) / w;
}


///
static dmat get_sigma_X(dvec &w, dvec &tau, dmat &X, int maxiterPCG,
	double tolPCG)
{
	int ncol = X.n_cols;
	dmat Sigma_iX1(Geno_NumSamp, ncol);
	for(int i = 0; i < ncol; i++)
		Sigma_iX1.col(i) = PCG_diag_sigma(w, tau, X.col(i), maxiterPCG, tolPCG);
	return(Sigma_iX1);
}


// ========================================================================= //

/// Calculate fixed and random effect coefficients given Y, X, tau
/// Input:  Y, X, tau, ...
/// Output: alpha, eta, W, ...
static void get_coeff(const dvec &y, const dmat &X, const dvec &tau,
	const List &family, const dvec &alpha0, const dvec &eta0,
	const dvec &offset, int maxiterPCG, int maxiter, double tolPCG,
	bool verbose,
	dvec &Y, dvec &mu, dvec &alpha, dvec &eta, dvec &W, dmat &cov,
	dvec &Sigma_iY, dmat &Sigma_iX)
{
	// initialize
	const double tol_coef = 0.1;
	Function fc_linkinv = wrap(family["linkinv"]);
	Function fc_mu_eta = wrap(family["mu.eta"]);
	Function fc_variance = wrap(family["variance"]);

	mu = as<dvec>(fc_linkinv(eta0));
	dvec mu_eta = as<dvec>(fc_mu_eta(eta0));
	Y = eta0 - offset + (y - mu)/mu_eta;
	W = (mu_eta % mu_eta) / as<dvec>(fc_variance(mu));

	// iterate ...
	dvec a0 = alpha0;
	for (int i=0; i < maxiter; i++)
	{
		get_coeff_w(Y, X, W, tau, maxiterPCG, tolPCG, Sigma_iY, Sigma_iX,
			cov, alpha, eta);

		eta += offset;
		mu = as<dvec>(fc_linkinv(eta));
		mu_eta = as<dvec>(fc_mu_eta(eta));
		Y = eta - offset + (y - mu)/mu_eta;
		W = (mu_eta % mu_eta) / as<dvec>(fc_variance(mu));

		if (max(abs(alpha - a0)/(abs(alpha) + abs(a0) + tol_coef)) < tol_coef)
			break;
		a0 = alpha;
	}
}


// Get average information (AI) for binary outcomes
static void get_AI_score(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &tau, const dvec &Sigma_iY, const dmat &Sigma_iX, const dmat &cov,
	int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff, int seed,
	double &outYPAPY, double &outTrace, double &outAI)
{
	dmat Sigma_iXt = Sigma_iX.t();
	dvec PY = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Y));
	dvec APY;
	get_crossprod_b_grm(PY, APY);
	// output
	outYPAPY = dot(PY, APY);
	outTrace = get_trace(Sigma_iX, X, w, tau, cov, nrun, maxiterPCG, tolPCG,
		traceCVcutoff, seed);
	dvec PAPY_1 = PCG_diag_sigma(w, tau, APY, maxiterPCG, tolPCG);
	dvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
	outAI = dot(APY, PAPY);
}

// Get average information (AI) for quantitative outcomes
static void get_AI_score_q(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &tau, const dvec &Sigma_iY, const dmat &Sigma_iX, const dmat &cov,
	int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff, int seed,
	double outYPAPY[], double outTrace[], dmat &outAI)
{
	dmat Sigma_iXt = Sigma_iX.t();
	dvec PY = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Y));
	dvec A0PY = PY;
	dvec APY;
	get_crossprod_b_grm(PY, APY);
	// get YPAPY
	outYPAPY[0] = dot(PY, APY);   // YPAPY
	outYPAPY[1] = dot(PY, A0PY);  // YPA0PY
	// get Trace
	get_trace_q(Sigma_iX, X, w, tau, cov, nrun, maxiterPCG, tolPCG,
		traceCVcutoff, seed, outTrace[0], outTrace[1]);
	// get AI
	dmat AI(2,2);
	dvec PA0PY_1 = PCG_diag_sigma(w, tau, A0PY, maxiterPCG, tolPCG);
	dvec PA0PY = PA0PY_1 - Sigma_iX * (cov * (Sigma_iXt * PA0PY_1));
	AI(0,0) = dot(A0PY, PA0PY);
	dvec PAPY_1 = PCG_diag_sigma(w, tau, APY, maxiterPCG, tolPCG);
	dvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
	AI(1,1) = dot(APY, PAPY);
	AI(1,0) = AI(0,1) = dot(A0PY, PAPY);
	outAI = AI;
}


// Update tau for binary outcomes
static dvec fitglmmaiRPCG(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &in_tau, const dvec &Sigma_iY, const dmat &Sigma_iX,
	const dmat &cov, int nrun, int maxiterPCG, double tolPCG, double tol,
	double traceCVcutoff, int seed)
{
	double YPAPY, Trace, AI;
	get_AI_score(Y, X, w, in_tau, Sigma_iY, Sigma_iX, cov, nrun,
		maxiterPCG, tolPCG, traceCVcutoff, seed,
		YPAPY, Trace, AI);
  	double score = YPAPY - Trace;
  	double Dtau = score / AI;
  	dvec tau = in_tau;
  	dvec tau0 = in_tau;
  	tau[1] = tau0[1] + Dtau;

  	for(size_t i=0; i < tau.n_elem; i++)
		if (tau[i] < tol) tau[i] = 0;

  	double step = 1.0;
  	while (tau[1] < 0.0)
  	{
		step *= 0.5;
		tau[1] = tau0[1] + step * Dtau;
  	}

  	for(size_t i=0; i < tau.n_elem; i++)
		if (tau[i] < tol) tau[i] = 0;

  	return tau;
}

// Update tau for quantitative outcomes
static dvec fitglmmaiRPCG_q(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &in_tau, const dvec &Sigma_iY, const dmat &Sigma_iX,
	const dmat &cov, int nrun, int maxiterPCG, double tolPCG, double tol,
	double traceCVcutoff, int seed)
{
	uvec zero_v = (in_tau < tol);
	double YPAPY[2], Trace[2];
	dmat AI;
	get_AI_score_q(Y, X, w, in_tau, Sigma_iY, Sigma_iX, cov, nrun,
		maxiterPCG, tolPCG, traceCVcutoff, seed,
		YPAPY, Trace, AI);
	dvec score(2);
	score[0] = YPAPY[1] - Trace[0];
	score[1] = YPAPY[0] - Trace[1];
	dvec Dtau = solve(AI, score);

	dvec tau0 = in_tau;
	dvec tau = tau0 + Dtau;
	tau.elem( find(zero_v % (tau < tol)) ).zeros();

	double step = 1.0;
	while (tau[0] < 0.0 || tau[1]  < 0.0)
	{
		step *= 0.5;
		tau = tau0 + step * Dtau;
		tau.elem( find(zero_v % (tau < tol)) ).zeros();
	}
	tau.elem( find(tau < tol) ).zeros();

	return tau;
}


// ========================================================================= //

/// Print a numeric vector
inline static void print_vec(const char *s, dvec &x, bool nl=true)
{
	Rprintf("%s(", s);
	for (size_t i=0; i < x.n_elem; i++)
	{
		if (i > 0) Rprintf(", ");
		Rprintf("%0.7g", x[i]);
	}
	if (nl) Rprintf(")\n"); else Rprintf(")");
}


// Fitting the null model with binary outcomes
RcppExport SEXP saige_fit_AI_PCG_binary(SEXP r_fit0, SEXP r_X, SEXP r_tau,
	SEXP r_param)
{
BEGIN_RCPP

	// parameters for fitting the model
	List param(r_param);
	const double tol = Rf_asReal(param["tol"]);
	const double tol_inv_2 = 1 / (tol*tol);
	const double tolPCG = Rf_asReal(param["tolPCG"]);
	const int seed = Rf_asInteger(param["seed"]);
	const int maxiter = Rf_asInteger(param["maxiter"]);
	const int maxiterPCG = Rf_asInteger(param["maxiterPCG"]);
	const int nrun = Rf_asInteger(param["nrun"]);
	const double traceCVcutoff = Rf_asReal(param["traceCVcutoff"]);
	const bool verbose = Rf_asLogical(param["verbose"])==TRUE;

	List fit0(r_fit0);
	dvec y = as<dvec>(fit0["y"]);
	dmat X = as<dmat>(r_X);
	dvec offset(y.size());
	if (Rf_isNull(fit0["offset"]))
		offset.zeros();
	else
		offset = as<dvec>(fit0["offset"]);

	List family = fit0["family"];
	Function fc_mu_eta = wrap(family["mu.eta"]);
	dvec eta = as<dvec>(fit0["linear.predictors"]);
	dvec eta0 = eta;
	dvec mu = as<dvec>(fit0["fitted.values"]);
	dvec mu_eta = as<dvec>(fc_mu_eta(eta0));
	dvec Y = eta - offset + (y - mu) / mu_eta;
	dvec alpha0 = as<dvec>(fit0["coefficients"]);
	dvec alpha = alpha0;
	dmat cov;

	dvec tau = as<dvec>(r_tau);
	dvec tau0 = tau;

	dvec re_Y, re_mu, re_alpha, re_eta, re_W, re_Sigma_iY;
	dmat re_cov, re_Sigma_iX;
	get_coeff(y, X, tau, family, alpha0, eta0, offset, maxiterPCG, maxiter,
		tolPCG, verbose,
		re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX);

	double YPAPY, Trace, AI;
	get_AI_score(re_Y, X, re_W, tau, re_Sigma_iY, re_Sigma_iX, re_cov, nrun,
		maxiterPCG, tolPCG, traceCVcutoff, seed,
		YPAPY, Trace, AI);

	tau[1] = std::max(0.0, tau0[1] + tau0[1]*tau0[1]*(YPAPY - Trace)/y.size());
	if (verbose)
	{
		Rprintf("Initial variance component estimates, tau:\n");
		Rprintf("    Sigma_E: %g, Sigma_G: %g\n", tau[0], tau[1]);
	}

	int iter = 1;
	for (; iter <= maxiter; iter++)
	{
		if (verbose)
		{
			Rprintf("Iteration %d:\n", iter);
			print_vec("    tau: ", tau);
			print_vec("    fixed coeff: ", alpha);
		}

		alpha0 = re_alpha;
		tau0 = tau;
		eta0 = eta;

		// find the next tau
		for (int itry=1; itry <= 11; itry++)
		{
			get_coeff(y, X, tau0, family, alpha0, eta0, offset, maxiterPCG,
				maxiter, tolPCG, verbose,
				re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX);
			tau = fitglmmaiRPCG(re_Y, X, re_W, tau0, re_Sigma_iY, re_Sigma_iX,
				re_cov, nrun, maxiterPCG, tolPCG, tol, traceCVcutoff, seed);
			if (max(tau) > tol_inv_2)
			{
				if (itry <= 10)
				{
					tau0[1] *= 0.5;
					if (verbose)
					{
						print_vec("    tau: ", tau, false);
						Rprintf(", large variance estimate observed, retry (%d) ...\n", itry);
						print_vec("    set new tau: ", tau0);
					}
					continue;
				} else {
					if (verbose) print_vec("tau: ", tau);
					throw std::overflow_error(
					"Large variance estimate observed in the iterations, model not converged!");
				}
			}
			break;
		}

		cov = re_cov; alpha = re_alpha; eta = re_eta;
		Y = re_Y; mu = re_mu;

		if (tau[1] == 0) break;
		if (max(abs(tau-tau0) / (abs(tau)+abs(tau0)+tol)) < tol) break;
	}

	get_coeff(y, X, tau, family, alpha0, eta0, offset, maxiterPCG, maxiter,
		tolPCG, verbose,
		re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX);
	cov = re_cov; alpha = re_alpha; eta = re_eta;
	Y = re_Y; mu = re_mu;

	if (verbose)
	{
		print_vec("Final tau: " , tau);
		print_vec("    fixed coeff: ", alpha);
	}

	return List::create(
		_["coefficients"] = SEXP_VEC(alpha),
		_["tau"] = SEXP_VEC(tau),
		_["linear.predictors"] = SEXP_VEC(eta),
		_["fitted.values"] = SEXP_VEC(mu),
		_["residuals"] = SEXP_VEC(y - mu),
		_["cov"] = cov,
		_["converged"] = bool(iter <= maxiter));

END_RCPP
}


// Fitting the null model with quantitative outcomes
RcppExport SEXP saige_fit_AI_PCG_quant(SEXP r_fit0, SEXP r_X, SEXP r_tau,
	SEXP r_param)
{
BEGIN_RCPP

	// parameters for fitting the model
	List param(r_param);
	const double tol = Rf_asReal(param["tol"]);
	const double tol_inv_2 = 1 / (tol*tol);
	const double tolPCG = Rf_asReal(param["tolPCG"]);
	const int seed = Rf_asInteger(param["seed"]);
	const int maxiter = Rf_asInteger(param["maxiter"]);
	const int maxiterPCG = Rf_asInteger(param["maxiterPCG"]);
	const int nrun = Rf_asInteger(param["nrun"]);
	const double traceCVcutoff = Rf_asReal(param["traceCVcutoff"]);
	const bool verbose = Rf_asLogical(param["verbose"])==TRUE;

	List fit0(r_fit0);
	dvec y = as<dvec>(fit0["y"]);
	const int n = y.size();
	dmat X = as<dmat>(r_X);
	dvec offset(y.size());
	if (Rf_isNull(fit0["offset"]))
		offset.zeros();
	else
		offset = as<dvec>(fit0["offset"]);

	List family = fit0["family"];
	Function fc_mu_eta = wrap(family["mu.eta"]);
	dvec eta = as<dvec>(fit0["linear.predictors"]);
	dvec eta0 = eta;
	dvec mu = as<dvec>(fit0["fitted.values"]);
	dvec mu_eta = as<dvec>(fc_mu_eta(eta0));
	dvec Y = eta - offset + (y - mu) / mu_eta;
	dvec alpha0 = as<dvec>(fit0["coefficients"]);
	dvec alpha = alpha0;
	dmat cov;

	dvec tau = as<dvec>(r_tau);
	dvec tau0 = tau;

	dvec re_Y, re_mu, re_alpha, re_eta, re_W, re_Sigma_iY;
	dmat re_cov, re_Sigma_iX;
	get_coeff(y, X, tau, family, alpha0, eta0, offset, maxiterPCG, maxiter,
		tolPCG, verbose,
		re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX);

	if (verbose)
	{
		Rprintf("Initial variance component estimates, tau:\n");
		Rprintf("    Sigma_E: %g, Sigma_G: %g\n", tau[0], tau[1]);
	}

	double YPAPY[2], Trace[2];
	dmat AI;
	get_AI_score_q(re_Y, X, re_W, tau, re_Sigma_iY, re_Sigma_iX, re_cov, nrun,
		maxiterPCG, tolPCG, traceCVcutoff, seed,
		YPAPY, Trace, AI);

	tau[0] = std::max(0.0, tau0[0] + tau0[0]*tau0[0]*(YPAPY[1] - Trace[0])/n);
	tau[1] = std::max(0.0, tau0[1] + tau0[1]*tau0[1]*(YPAPY[0] - Trace[1])/n);

	int iter = 1;
	for (; iter <= maxiter; iter++)
	{
		if (verbose)
		{
			Rprintf("Iteration %d:\n", iter);
			print_vec("    tau: ", tau);
			print_vec("    fixed coeff: ", alpha);
		}

		alpha0 = re_alpha;
		tau0 = tau;
		eta0 = eta;

		// find the next tau
		for (int itry=1; itry <= 11; itry++)
		{
			get_coeff(y, X, tau0, family, alpha0, eta0, offset, maxiterPCG,
				maxiter, tolPCG, verbose,
				re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX);
			tau = fitglmmaiRPCG_q(re_Y, X, re_W, tau0, re_Sigma_iY,
				re_Sigma_iX, re_cov, nrun, maxiterPCG, tolPCG,
				tol, traceCVcutoff, seed);
			if (max(tau) > tol_inv_2)
			{
				if (verbose) print_vec("tau: ", tau, false);
				if (itry <= 10)
				{
					tau0[1] *= 0.5;
					if (verbose)
					{
						Rprintf(", large variance estimate observed, retry (%d) ...\n", itry);
						print_vec("    set new tau: ", tau0);
					}
					continue;
				} else {
					if (verbose) Rprintf("\n");
					throw std::overflow_error(
					"Large variance estimate observed in the iterations, model not converged!");
				}
			}
			break;
		}

		cov = re_cov; alpha = re_alpha; eta = re_eta;
		Y = re_Y; mu = re_mu;

		if (tau[0] <= 0)
		{
			print_vec("    tau: ", tau);
			throw std::overflow_error("Sigma_E = 0, model not converged!");
		}
		if (max(abs(tau-tau0) / (abs(tau)+abs(tau0)+tol)) < tol) break;
	}

	get_coeff(y, X, tau, family, alpha0, eta0, offset, maxiterPCG, maxiter,
		tolPCG, verbose,
		re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX);
	cov = re_cov; alpha = re_alpha; eta = re_eta;
	Y = re_Y; mu = re_mu;

	if (verbose)
	{
		print_vec("Final tau: " , tau);
		print_vec("    fixed coeff: ", alpha);
	}

	return List::create(
		_["coefficients"] = SEXP_VEC(alpha),
		_["tau"] = SEXP_VEC(tau),
		_["linear.predictors"] = SEXP_VEC(eta),
		_["fitted.values"] = SEXP_VEC(mu),
		_["residuals"] = SEXP_VEC(y - mu),
		_["cov"] = cov,
		_["converged"] = bool(iter <= maxiter));

END_RCPP
}



// ========================================================================= //

/// Calculate variance ratio for binary outcomes
RcppExport SEXP saige_calc_var_ratio_binary(SEXP r_fit0, SEXP r_glmm,
	SEXP r_noK, SEXP r_param, SEXP r_marker_list)
{
BEGIN_RCPP

	List fit0(r_fit0);
	List glmm(r_glmm);
	List obj_noK(r_noK);
	List param(r_param);
	IntegerVector rand_index(r_marker_list);

	// parameters for fitting the model
	const double tolPCG = Rf_asReal(param["tolPCG"]);
	const int maxiterPCG = Rf_asInteger(param["maxiterPCG"]);
	const double ratioCVcutoff = Rf_asReal(param["ratioCVcutoff"]);
	int num_marker = Rf_asInteger(param["num.marker"]);
	const bool verbose = Rf_asLogical(param["verbose"])==TRUE;

	List family = fit0["family"];
	Function fc_mu_eta = wrap(family["mu.eta"]);
	Function fc_variance = wrap(family["variance"]);

	dvec eta = as<dvec>(fit0["linear.predictors"]);
	dvec mu = as<dvec>(fit0["fitted.values"]);
	dvec mu_eta = as<dvec>(fc_mu_eta(eta));
	dvec W = (mu_eta % mu_eta) / as<dvec>(fc_variance(mu));
	dvec tau = as<dvec>(glmm["tau"]);
	dmat X1 = as<dmat>(obj_noK["X1"]);
	dmat Sigma_iX = get_sigma_X(W, tau, X1, maxiterPCG, tolPCG);

	dvec y = as<dvec>(fit0["y"]);
	dmat noK_XXVX_inv = as<dmat>(obj_noK["XXVX_inv"]);
	dmat noK_XV = as<dmat>(obj_noK["XV"]);

	double ratioCV = ratioCVcutoff + 0.1;
	int num_tested = 0, snp_idx = 0;
	const int num_rand_snp = rand_index.length();

	dvec G0(Geno_NumSamp);
	vector<int> buf_idx(Geno_NumSamp);
	vector<int> lst_idx;
	vector<double> lst_maf, lst_mac, lst_var1, lst_var2, lst_ratio;

	while (ratioCV > ratioCVcutoff && snp_idx < num_rand_snp)
	{
		while (num_tested < num_marker && snp_idx < num_rand_snp)
		{
			const int i_snp = rand_index[snp_idx++];
			get_geno_ds(i_snp - 1, G0);

			double AF, AC;
			int Num;
			f64_af_ac_impute(&G0[0], Geno_NumSamp, AF, AC, Num, &buf_idx[0]);
			if (AF > 0.5)
			{
				f64_sub(Geno_NumSamp, 2, &G0[0]);
				AC = 2*Num - AC;
				AF = 1 - AF;
			}
			if (AC <= 20) continue;  // suggested by the paper

			// adjusted genotypes
			dvec G = G0 - noK_XXVX_inv * (noK_XV * G0);
			dvec g = G / sqrt(AC);
			dvec Sigma_iG = PCG_diag_sigma(W, tau, G, maxiterPCG, tolPCG);
			dvec adj = Sigma_iX * mat_inv(X1.t() * Sigma_iX) * X1.t() * Sigma_iG;

			double var1 = (sum(G % Sigma_iG) - sum(G % adj)) / AC;
			double var2 = sum(mu % (1 - mu) % g % g);
			double ratio = var1 / var2;

			num_tested ++;
			lst_idx.push_back(i_snp);
			lst_maf.push_back(AF);
			lst_mac.push_back(AC);
			lst_var1.push_back(var1);
			lst_var2.push_back(var2);
			lst_ratio.push_back(ratio);
			if (verbose)
			{
				Rprintf("%6d, maf: %0.4f, mac: %g,\tratio: %0.4f (var1: %.3g, var2: %.3g)\n",
					num_tested, AF, AC, ratio, var1, var2);
			}
		}

		ratioCV = calcCV(lst_ratio);
		if (ratioCV > ratioCVcutoff)
		{
			if (verbose)
			{
				Rprintf(
					"CV for variance ratio estimate using %d markers is %g > ratioCVcutoff (%g), try more markers ...\n",
					num_marker, ratioCV, ratioCVcutoff);
			}
			num_marker += 10;
		}
	}

	return DataFrame::create(
		_["id"] = lst_idx,    _["maf"] = lst_maf,   _["mac"] = lst_mac,
		_["var1"] = lst_var1, _["var2"] = lst_var2, _["ratio"] = lst_ratio);

END_RCPP
}


/// Calculate variance ratio for quantitative outcomes
RcppExport SEXP saige_calc_var_ratio_quant(SEXP r_fit0, SEXP r_glmm,
	SEXP r_noK, SEXP r_param, SEXP r_marker_list)
{
BEGIN_RCPP

	List fit0(r_fit0);
	List glmm(r_glmm);
	List obj_noK(r_noK);
	List param(r_param);
	IntegerVector rand_index(r_marker_list);

	// parameters for fitting the model
	const double tolPCG = Rf_asReal(param["tolPCG"]);
	const int maxiterPCG = Rf_asInteger(param["maxiterPCG"]);
	const double ratioCVcutoff = Rf_asReal(param["ratioCVcutoff"]);
	int num_marker = Rf_asInteger(param["num.marker"]);
	const bool verbose = Rf_asLogical(param["verbose"])==TRUE;

	List family = fit0["family"];
	Function fc_mu_eta = wrap(family["mu.eta"]);
	Function fc_variance = wrap(family["variance"]);

	dvec eta = as<dvec>(fit0["linear.predictors"]);
	dvec mu = as<dvec>(fit0["fitted.values"]);
	dvec mu_eta = as<dvec>(fc_mu_eta(eta));
	dvec W = (mu_eta % mu_eta) / as<dvec>(fc_variance(mu));
	dvec tau = as<dvec>(glmm["tau"]);
	dmat X1 = as<dmat>(obj_noK["X1"]);
	dmat Sigma_iX = get_sigma_X(W, tau, X1, maxiterPCG, tolPCG);

	dvec y = as<dvec>(fit0["y"]);
	dmat noK_XXVX_inv = as<dmat>(obj_noK["XXVX_inv"]);
	dmat noK_XV = as<dmat>(obj_noK["XV"]);

	double ratioCV = ratioCVcutoff + 0.1;
	int num_tested = 0, snp_idx = 0;
	const int num_rand_snp = rand_index.length();

	dvec G0(Geno_NumSamp);
	vector<int> buf_idx(Geno_NumSamp);
	vector<int> lst_idx;
	vector<double> lst_maf, lst_mac, lst_var1, lst_var2, lst_ratio;

	while (ratioCV > ratioCVcutoff && snp_idx < num_rand_snp)
	{
		while (num_tested < num_marker && snp_idx < num_rand_snp)
		{
			const int i_snp = rand_index[snp_idx++];
			get_geno_ds(i_snp - 1, G0);

			double AF, AC;
			int Num;
			f64_af_ac_impute(&G0[0], Geno_NumSamp, AF, AC, Num, &buf_idx[0]);
			if (AF > 0.5)
			{
				f64_sub(Geno_NumSamp, 2, &G0[0]);
				AC = 2*Num - AC;
				AF = 1 - AF;
			}
			if (AC <= 20) continue;  // suggested by the paper

			// adjusted genotypes
			dvec G = G0 - noK_XXVX_inv * (noK_XV * G0);
			dvec g = G / sqrt(AC);
			dvec Sigma_iG = PCG_diag_sigma(W, tau, G, maxiterPCG, tolPCG);
			dvec adj = Sigma_iX * mat_inv(X1.t() * Sigma_iX) * X1.t() * Sigma_iG;

			double var1 = (sum(G % Sigma_iG) - sum(G % adj)) / AC;
			double var2 = sum(g % g);
			double ratio = var1 / var2;

			num_tested ++;
			lst_idx.push_back(i_snp);
			lst_maf.push_back(AF);
			lst_mac.push_back(AC);
			lst_var1.push_back(var1);
			lst_var2.push_back(var2);
			lst_ratio.push_back(ratio);
			if (verbose)
			{
				Rprintf("%6d, maf: %0.4f, mac: %g,\tratio: %0.4f (var1: %.3g, var2: %.3g)\n",
					num_tested, AF, AC, ratio, var1, var2);
			}
		}

		ratioCV = calcCV(lst_ratio);
		if (ratioCV > ratioCVcutoff)
		{
			if (verbose)
			{
				Rprintf(
					"CV for variance ratio estimate using %d markers is %g > ratioCVcutoff (%g)\n",
					num_marker, ratioCV, ratioCVcutoff);
			}
			num_marker += 10;
			if (verbose) Rprintf("try %d markers ...\n", num_marker);
		}
	}

	return DataFrame::create(
		_["id"] = lst_idx,    _["maf"] = lst_maf,   _["mac"] = lst_mac,
		_["var1"] = lst_var1, _["var2"] = lst_var2, _["ratio"] = lst_ratio);

END_RCPP
}
