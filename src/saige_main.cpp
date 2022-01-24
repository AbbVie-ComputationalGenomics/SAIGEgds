// ===========================================================
//
// saige_main.cpp: SAIGE association analysis
//
// Copyright (C) 2019-2021    Xiuwen Zheng / AbbVie-ComputationalGenomics
//
// This file is part of SAIGEgds.
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

/* To ensure atanpi, cospi, sinpi, tanpi are defined */
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
#   define __STDC_WANT_IEC_60559_FUNCS_EXT__    1
#endif

#include <Rconfig.h>
#include "vectorization.h"
#include <Rcpp.h>
#include <algorithm>
#include <Rdefines.h>
#include <math.h>
#include <Rmath.h>


using namespace Rcpp;
using namespace vectorization;



// ========================================================================= //
// internal functions

/// SPAtest
extern "C" double Saddle_Prob(double q, double m1, double var1, size_t n_g,
	const double mu[], const double g[], double cutoff, bool &converged, double *p_noadj);

extern "C" double Saddle_Prob_Fast(double q, double m1, double var1, size_t n_g,
	const double mu[], const double g[], size_t n_nonzero, const int nonzero_idx[],
	double cutoff, bool &converged, double buf_spa[], double *p_noadj);

/// square
inline double sq(double v) { return v*v; }
/// minor allele frequency
inline double MAF(double af) { return std::min(af, 1-af); }


// ========================================================================= //

static double threshold_maf = 0;  //< the threshold of MAF filter
static double threshold_mac = 0;  //< the threshold of MAC filter
static double threshold_missing = 1;      //< the threshold of missing proportion per variant
static double threshold_pval_spa = 0.05;  //< the threshold of p-value filter for SPA

static int mod_NSamp = 0;   //< the number of samples
static int mod_NCoeff = 0;  //< the number of beta coefficients

static double *mod_tau = NULL;           //< variance components: tau[0], tau[1]
static double *mod_y = NULL;             //< a n_samp-length vector
static double *mod_mu = NULL;            //< a n_samp-length vector
static double *mod_y_mu = NULL;          //< a n_samp-length vector, y-mu
static double *mod_mu2 = NULL;           //< a n_samp-length vector, mu*(1-mu)
static double *mod_t_XXVX_inv = NULL;    //< a K-by-n_samp matrix
static double *mod_XV = NULL;            //< a K-by-n_samp matrix
static double *mod_t_XVX_inv_XV = NULL;  //< a K-by-n_samp matrix
static double *mod_XVX = NULL;           //< a K-by-K matrix
static double *mod_t_X = NULL;           //< a K-by-n_samp matrix
static double *mod_S_a = NULL;           //< a K-length vector

static double mod_varRatio = 0;

static double *buf_dosage = NULL;    //< temporary buffer for real dosages
static double *buf_coeff = NULL;     // beta coefficients
static double *buf_adj_g = NULL;     //< genotype after adjusting for fixed effects
static int *buf_index = NULL;
static double *buf_B = NULL;
static double *buf_g_tilde = NULL;   //< est. adj. genotypes
static double *buf_X1 = NULL;        //< ncol(X1)
static double *buf_spa = NULL;       //< buffer for SPA calculation

// aggregate tests
static double threshold_summac = 0;  //< the threshold of weighted sum MAC
static int num_wbeta = 0;            //< # of beta parameters
static double *buf_wbeta = NULL;     //< beta parameters
static int num_unitsz = 0;           //< max unit size
static double *buf_unitsz = NULL;    //< length() = max unit size
static double threshold_acatv_mac = 0;  //< the threshold of weighted sum MAC

#define IDX_i    buf_index[i]


/// initialize internal parameters from the model object
RcppExport SEXP saige_score_test_init(SEXP model)
{
BEGIN_RCPP
	List M(model);
	// threshold setting
	threshold_maf = Rf_asReal(M["maf"]);
	if (!R_FINITE(threshold_maf)) threshold_maf = -1;
	threshold_mac = Rf_asReal(M["mac"]);
	if (!R_FINITE(threshold_mac)) threshold_mac = -1;
	threshold_missing = Rf_asReal(M["missing"]);
	if (!R_FINITE(threshold_missing)) threshold_missing = 1;
	threshold_pval_spa = Rf_asReal(M["spa.pval"]);
	if (!R_FINITE(threshold_pval_spa)) threshold_pval_spa = 0.05;
	// model parameters
	mod_NSamp = Rf_length(M["y"]);
	mod_NCoeff = NumericMatrix(wrap(M["XV"])).nrow();
	mod_tau = REAL(M["tau"]);
	mod_y = REAL(M["y"]);
	mod_mu = REAL(M["mu"]);
	mod_y_mu = REAL(M["y_mu"]);
	mod_mu2 = REAL(M["mu2"]);
	mod_t_XXVX_inv = REAL(M["t_XXVX_inv"]);
	mod_XV = REAL(M["XV"]);
	mod_t_XVX_inv_XV = REAL(M["t_XVX_inv_XV"]);
	mod_XVX = REAL(M["XVX"]);
	mod_t_X = REAL(M["t_X"]);
	mod_S_a = REAL(M["S_a"]);
	mod_varRatio = Rf_asReal(M["var.ratio"]);
	// buffer
	buf_dosage = REAL(M["buf_dosage"]);
	buf_coeff = REAL(M["buf_coeff"]);
	buf_adj_g = REAL(M["buf_adj_g"]);
	buf_index = INTEGER(M["buf_index"]);
	buf_B = REAL(M["buf_B"]);
	buf_g_tilde = REAL(M["buf_g_tilde"]);
	buf_X1 = REAL(M["buf_X1"]);
	buf_spa = REAL(M["buf_spa"]);
	// buffer for aggregate tests
	threshold_summac = Rf_asReal(M["summac"]);
	if (!R_FINITE(threshold_summac)) threshold_summac = -1;
	threshold_acatv_mac = Rf_asReal(M["acatv_mac"]);
	if (!R_FINITE(threshold_acatv_mac)) threshold_acatv_mac = 10;
	num_wbeta = Rf_length(M["buf_wbeta"]) / 2;  // # of columns
	buf_wbeta = REAL(M["buf_wbeta"]);
	num_unitsz = Rf_length(M["buf_unitsz"]) / 5;
	buf_unitsz = REAL(M["buf_unitsz"]);
END_RCPP
}



// ========================================================================= //

static const char *ERR_DS_TYPE = "Invalid type of dosages.";
static const char *ERR_DS_LEN  = "Invalid length of dosages: %d.";
static const char *ERR_DS_MAT  = "Input dosage should be a matrix.";
static const char *ERR_DS_DIM  = "Invalid dimension of dosages: %dx%d.";

/// get numeric dosages from an R object
static double *get_ds(SEXP ds, R_xlen_t n, R_xlen_t start)
{
	const R_xlen_t ntot = Rf_xlength(ds);
	if (start < 0) start = 0;
	if (n - start > ntot)
		Rf_error(ERR_DS_LEN, ntot);

	const double NaN = R_NaN;
	R_xlen_t i = 0;
	switch (TYPEOF(ds))
	{
	case REALSXP:
		return REAL(ds) + start;
	case INTSXP:
		for (const int *p = INTEGER(ds) + start; i < n; i++)
			buf_dosage[i] = (p[i] != NA_INTEGER) ? p[i] : NaN;
		return buf_dosage;
	case RAWSXP:
		for (const Rbyte *p = RAW(ds) + start; i < n; i++)
			buf_dosage[i] = (p[i] != Rbyte(0xFF)) ? p[i] : NaN;
		return buf_dosage;
	}
	Rf_error(ERR_DS_TYPE);
	return NULL;
}

/// single variant test for quantitative outcomes
static bool single_test_quant(size_t num_samp, double G[], double &oAF,
	double &omac, int &onum, double &obeta, double &oSE, double &opval)
{
	// calc allele freq, and impute geno using the mean
	double AF, AC;
	int Num;
	f64_af_ac_impute(&G[0], num_samp, AF, AC, Num, buf_index);

	const double maf = std::min(AF, 1 - AF);
	const double mac = std::min(AC, 2*Num - AC);
	const double missing = double(num_samp - Num) / num_samp;
	if ((Num > 0) && (maf > 0) && (maf >= threshold_maf) &&
		(mac >= threshold_mac) && (missing <= threshold_missing))
	{
		bool minus = (AF > 0.5);
		if (minus) f64_sub(mod_NSamp, 2, &G[0]);

		double pval, beta;
		size_t n_nonzero = 0;
		const double inv_sqrt_mac = 1.0 / sqrt(mac);
		const double inv_mac = 1.0 / mac;

		const bool is_sparse = maf < 0.05;
		if (is_sparse)
		{
			// get the number of nonzeros and the nonzero indices
			n_nonzero = f64_nonzero_index(mod_NSamp, &G[0], buf_index);
			// buf_coeff = XVX_inv_XV * G
			f64_mul_mat_vec_sp(n_nonzero, buf_index, mod_NCoeff,
				mod_t_XVX_inv_XV, &G[0], buf_coeff);
			// buf_B = t(X) * buf_coeff
			f64_mul_mat_vec_sub(n_nonzero, buf_index, mod_NCoeff, mod_t_X,
				buf_coeff, buf_B);
			// g_tilde = G - B
			for (size_t i=0; i < n_nonzero; i++)
				buf_g_tilde[i] = G[IDX_i] - buf_B[i];
			// var2 = t(buf_coeff) %*% XVX %*% buf_coeff - sum(B^2) + sum(g_tilde^2)
			double var2 = f64_sum_mat_vec(mod_NCoeff, mod_XVX, buf_coeff);
			for (size_t i=0; i < n_nonzero; i++)
				var2 += sq(buf_g_tilde[i]) - sq(buf_B[i]);
			double var1 = var2 * inv_mac * mod_varRatio;
			// S1 = sum(y_mu .* g_tilde)
			double S1 = 0;
			for (size_t i=0; i < n_nonzero; i++)
				S1 += mod_y_mu[IDX_i] * buf_g_tilde[i];
			// buf_X1 = t(X1) * (y-mu)
			f64_mul_mat_vec_sp(n_nonzero, buf_index, mod_NCoeff, mod_t_X,
				mod_y_mu, buf_X1);
			// S2 = sum((buf_X1 - mod_S_a) .* buf_coeff)
			double S2 = 0;
			for (int i=0; i < mod_NCoeff; i++)
				S2 += (buf_X1[i] - mod_S_a[i]) * buf_coeff[i];
			//
			double Tstat = (S1 + S2) * inv_sqrt_mac / mod_tau[0];
			pval = ::Rf_pchisq(Tstat*Tstat/var1, 1, FALSE, FALSE);
			beta = Tstat / var1 * inv_sqrt_mac;

		} else {
			// calculate adjusted genotypes: adj_g = G - XXVX_inv * (XV * G)
			// buf_coeff = XV * G
			f64_mul_mat_vec(mod_NSamp, mod_NCoeff, mod_XV, &G[0], buf_coeff);
			// buf_adj_g = G - XXVX_inv * buf_coeff
			f64_sub_mul_mat_vec(mod_NSamp, mod_NCoeff, &G[0], mod_t_XXVX_inv,
				buf_coeff, buf_adj_g);

			// inner product
			double S, var;
			// S = sum((y - mu) .* buf_adj_g
			// var = sum(buf_adj_g .* buf_adj_g)
			f64_dot_sp(mod_NSamp, mod_y_mu, buf_adj_g, S, var);
			double Tstat = S * inv_sqrt_mac / mod_tau[0];
			var *= inv_mac * mod_varRatio;

			// p-value and beta
			pval = ::Rf_pchisq(Tstat*Tstat/var, 1, FALSE, FALSE);
			beta = Tstat / var * inv_sqrt_mac;
		}

		if (minus) beta = -beta;
		double SE = fabs(beta/::Rf_qnorm5(pval/2, 0, 1, TRUE, FALSE));

		// output
		oAF = AF; omac = mac; onum = Num;
		obeta = beta; oSE = SE; opval = pval;
		return true;
	} else
		return false;
}

/// single variant test for binary outcomes
static bool single_test_bin(size_t num_samp, double G[],
	double &oAF, double &omac, int &onum, double &obeta, double &oSE,
	double &opval, double &opval_noadj, bool &oconverged)
{
	// calc allele freq, and impute geno using the mean
	double AF, AC;
	int Num;
	f64_af_ac_impute(&G[0], num_samp, AF, AC, Num, buf_index);

	const double maf = std::min(AF, 1 - AF);
	const double mac = std::min(AC, 2*Num - AC);
	const double missing = double(num_samp - Num) / num_samp;
	if ((Num > 0) && (maf > 0) && (maf >= threshold_maf) &&
		(mac >= threshold_mac) && (missing <= threshold_missing))
	{
		bool minus = (AF > 0.5);
		if (minus) f64_sub(mod_NSamp, 2, &G[0]);

		double pval_noadj, beta;
		size_t n_nonzero = 0;
		const bool is_sparse = maf < 0.05;
		if (is_sparse)
		{
			// get the number of nonzeros and the nonzero indices
			n_nonzero = f64_nonzero_index(mod_NSamp, &G[0], buf_index);
			// buf_coeff = XVX_inv_XV * G
			f64_mul_mat_vec_sp(n_nonzero, buf_index, mod_NCoeff,
				mod_t_XVX_inv_XV, &G[0], buf_coeff);
			// buf_B = t(X) * buf_coeff
			f64_mul_mat_vec_sub(n_nonzero, buf_index, mod_NCoeff, mod_t_X,
				buf_coeff, buf_B);
			// g_tilde = G - B
			for (size_t i=0; i < n_nonzero; i++)
				buf_g_tilde[i] = G[IDX_i] - buf_B[i];
			// var2 = t(buf_coeff) %*% XVX %*% buf_coeff - sum(B^2 .* mu2) + sum(g_tilde^2 .* mu2)
			double var2 = f64_sum_mat_vec(mod_NCoeff, mod_XVX, buf_coeff);
			for (size_t i=0; i < n_nonzero; i++)
				var2 += (sq(buf_g_tilde[i]) - sq(buf_B[i])) * mod_mu2[IDX_i];
			double var1 = var2 * mod_varRatio;
			// S1 = sum(y_mu .* g_tilde)
			double S1 = 0;
			for (size_t i=0; i < n_nonzero; i++)
				S1 += mod_y_mu[IDX_i] * buf_g_tilde[i];
			// buf_X1 = t(X1) * (y-mu)
			f64_mul_mat_vec_sp(n_nonzero, buf_index, mod_NCoeff, mod_t_X,
				mod_y_mu, buf_X1);
			// S2 = sum((buf_X1 - mod_S_a) .* buf_coeff)
			double S2 = 0;
			for (int i=0; i < mod_NCoeff; i++)
				S2 += (buf_X1[i] - mod_S_a[i]) * buf_coeff[i];
			//
			double S = S1 + S2;
			pval_noadj = ::Rf_pchisq(S*S/var1, 1, FALSE, FALSE);
			beta = S / var1;

		} else {
			// calculate adjusted genotypes: adj_g = G - XXVX_inv * (XV * G)
			// buf_coeff = XV * G
			f64_mul_mat_vec(mod_NSamp, mod_NCoeff, mod_XV, &G[0], buf_coeff);
			// buf_adj_g = G - XXVX_inv * buf_coeff
			f64_sub_mul_mat_vec(mod_NSamp, mod_NCoeff, &G[0], mod_t_XXVX_inv,
				buf_coeff, buf_adj_g);
			// inner product
			double S, var;
			// S = sum((y - mu) .* buf_adj_g)
			// var = sum(mu*(1-mu) .* buf_adj_g .* buf_adj_g)
			f64_dot_sp2(mod_NSamp, mod_y_mu, mod_mu2, buf_adj_g, S, var);
			var *= mod_varRatio;
			// p-value and beta
			pval_noadj = ::Rf_pchisq(S*S/var, 1, FALSE, FALSE);
			beta = S / var;
		}

		double pval = pval_noadj;
		bool converged = R_FINITE(pval_noadj) != 0;

		// need further SPAtest or not?
		if (converged && (pval_noadj <= threshold_pval_spa))
		{
			// calculate adjusted genotypes
			if (is_sparse)
			{
				// need adjusted genotypes, adj_g = G - XXVX_inv * (XV * G)
				// buf_coeff = XV * G
				f64_mul_mat_vec_sp(n_nonzero, buf_index, mod_NCoeff,
					mod_XV, &G[0], buf_coeff);
				// buf_adj_g = G - XXVX_inv * buf_coeff
				f64_sub_mul_mat_vec(mod_NSamp, mod_NCoeff, &G[0], mod_t_XXVX_inv,
					buf_coeff, buf_adj_g);
			}
			// minor allele count
			double AC2 = minus ? (2*Num - AC) : AC;
			// adj_g = adj_g / sqrt(AC2)
			f64_mul(mod_NSamp, 1/sqrt(AC2), buf_adj_g);
			// q = sum(y .* adj_g)
			double q = f64_dot(mod_NSamp, mod_y, buf_adj_g);
			double m1, var2;
			// m1 = sum(mu .* adj_g)
			// var2 = sum(mu*(1-mu) .* adj_g .* adj_g)
			f64_dot_sp2(mod_NSamp, mod_mu, mod_mu2, buf_adj_g, m1, var2);
			double var1 = var2 * mod_varRatio;
			double Tstat = q - m1;
			double qtilde = Tstat/sqrt(var1) * sqrt(var2) + m1;

			// need the number of nonzeros and the nonzero indices
			if (!is_sparse)
				n_nonzero = f64_nonzero_index(mod_NSamp, &G[0], buf_index);

			// call Saddle_Prob in SPAtest
			pval = Saddle_Prob_Fast(qtilde, m1, var2, mod_NSamp, mod_mu,
				buf_adj_g, n_nonzero, buf_index, 2, converged, buf_spa, NULL);
			if (pval==0 && pval_noadj>0)
				{ pval = pval_noadj; converged = false; }

			// effect size
			beta = (Tstat / var1) / sqrt(AC2);
		}

		if (minus) beta = -beta;
		double SE = fabs(beta / ::Rf_qnorm5(pval/2, 0, 1, TRUE, FALSE));

		// output
		oAF = AF; omac = mac; onum = Num; obeta = beta; oSE = SE;
		opval = pval; opval_noadj = pval_noadj;
		oconverged = converged;
		return true;
	} else
		return false;
}


// ====================================

/// calculate single-variant p-values for quantitative outcomes
RcppExport SEXP saige_score_test_quant(SEXP dosage)
{
BEGIN_RCPP
	// get numeric dosage
	if (Rf_xlength(dosage) != mod_NSamp)
		Rf_error(ERR_DS_LEN, Rf_xlength(dosage));
	double *G = get_ds(dosage, mod_NSamp, 0);

	int num = 0;
	double AF, mac, beta, SE, pval;
	AF = mac = beta = SE = pval = R_NaN;

	if (single_test_quant(mod_NSamp, G, AF, mac, num, beta, SE, pval))
	{
		NumericVector ans(6);
		ans[0] = AF;    ans[1] = mac;   ans[2] = num;
		ans[3] = beta;  ans[4] = SE;    ans[5] = pval;
		return ans;
	} else
		return R_NilValue;
END_RCPP
}

/// calculate single-variant p-values for binary outcomes
RcppExport SEXP saige_score_test_bin(SEXP dosage)
{
BEGIN_RCPP
	// get numeric dosage
	if (Rf_xlength(dosage) != mod_NSamp)
		Rf_error(ERR_DS_LEN, Rf_xlength(dosage));
	double *G = get_ds(dosage, mod_NSamp, 0);

	int num = 0;
	bool converged = 0;
	double AF, mac, beta, SE, pval, pval_noadj;
	AF = mac = beta = SE = pval = pval_noadj = R_NaN;

	if (single_test_bin(mod_NSamp, G, AF, mac, num, beta, SE, pval, pval_noadj,
		converged))
	{
		NumericVector ans(8);
		ans[0] = AF;    ans[1] = mac;   ans[2] = num;
		ans[3] = beta;  ans[4] = SE;    ans[5] = pval;
		ans[6] = pval_noadj;
		ans[7] = converged ? 1 : 0;
		return ans;
	} else
		return R_NilValue;
END_RCPP
}


// ========================================================================= //
// Aggregate Tests

// get the number of columns
static int ds_mat_dim(SEXP dosage)
{
	if (Rf_isNull(dosage)) return 0;
	// get nrow and ncol
	SEXP dm = GET_DIM(dosage);
	if (Rf_isNull(dm) || Rf_length(dm) != 2)
		Rf_error(ERR_DS_MAT);
	const int *p = INTEGER(dm);
	if (mod_NSamp != p[0])
		Rf_error(ERR_DS_DIM, p[0], p[1]);
	const int n_snp = p[1];
	if (n_snp > num_unitsz)
		Rf_error(ERR_DS_DIM, p[0], p[1]);
	return n_snp;
}

static void ds_mat_mafmac(double out_maf[], double out_mac[], SEXP ds_mat,
	size_t n_samp, size_t n_snp)
{
	size_t i = 0;
	switch (TYPEOF(ds_mat))
	{
	case RAWSXP:
		for (const Rbyte *p = RAW(ds_mat); i < n_snp; i++, p+=n_samp)
		{
			int n=0, s=0;
			for (size_t j=0; j < n_samp; j++)
				if (p[j] != Rbyte(0xFF)) { n++; s+=p[j]; }
			out_maf[i] = (n > 0) ? MAF((double)s / (2*n)) : R_NaN;
			out_mac[i] = std::min(s, 2*n-s);
		}
		break;
	case INTSXP:
		for (const int *p = INTEGER(ds_mat); i < n_snp; i++, p+=n_samp)
		{
			int n=0, s=0;
			for (size_t j=0; j < n_samp; j++)
				if (p[j] != NA_INTEGER) { n++; s+=p[j]; }
			out_maf[i] = (n > 0) ? MAF((double)s / (2*n)) : R_NaN;
			out_mac[i] = std::min(s, 2*n-s);
		}
		break;
	case REALSXP:
		for (const double *p = REAL(ds_mat); i < n_snp; i++, p+=n_samp)
		{
			int n=0; double s=0;
			for (size_t j=0; j < n_samp; j++)
				if (R_FINITE(p[j])) { n++; s+=p[j]; }
			out_maf[i] = (n > 0) ? MAF(s / (2*n)) : R_NaN;
			out_mac[i] = std::min(s, 2*n-s);
		}
		break;
	default:
		if (!Rf_isNull(ds_mat)) Rf_error(ERR_DS_TYPE);
	}
}

static void ds_mat_burden(double out_ds[], SEXP ds_mat, size_t n_samp,
	size_t n_snp, const double weight[])
{
	// initialize
	memset(out_ds, 0, sizeof(double)*n_samp);
	// data type
	size_t i = 0;
	switch (TYPEOF(ds_mat))
	{
	case RAWSXP:
		for (const Rbyte *s = RAW(ds_mat); i < n_snp; i++, s+=n_samp)
		{
			if (R_FINITE(weight[i]))
			{
				// get average
				int n=0, sum=0;
				for (size_t j=0; j < n_samp; j++)
					if (s[j] != Rbyte(0xFF)) { n++; sum += s[j]; }
				// imputed by the mean
				double *p = out_ds, w = weight[i], m = (double)sum / n;
				if (sum <= n)
				{
					for (size_t j=0; j < n_samp; j++)
						p[j] += (s[j] != Rbyte(0xFF)) ? (s[j] * w) : (m * w);
				} else {
					// flip
					m = 2 - m;
					for (size_t j=0; j < n_samp; j++)
						p[j] += (s[j] != Rbyte(0xFF)) ? ((2-s[j]) * w) : (m * w);
				}
			}
		}
		break;
	case INTSXP:
		for (const int *s = INTEGER(ds_mat); i < n_snp; i++, s+=n_samp)
		{
			if (R_FINITE(weight[i]))
			{
				// get average
				int n=0, sum=0;
				for (size_t j=0; j < n_samp; j++)
					if (s[j] != NA_INTEGER) { n++; sum += s[j]; }
				// imputed by the mean
				double *p = out_ds, w = weight[i], m = (double)sum / n;
				if (sum <= n)
				{
					for (size_t j=0; j < n_samp; j++)
						p[j] += (s[j] != NA_INTEGER) ? (s[j] * w) : (m * w);
				} else {
					// flip
					m = 2 - m;
					for (size_t j=0; j < n_samp; j++)
						p[j] += (s[j] != NA_INTEGER) ? ((2-s[j]) * w) : (m * w);
				}
			}
		}
		break;
	case REALSXP:
		for (const double *s = REAL(ds_mat); i < n_snp; i++, s+=n_samp)
		{
			if (R_FINITE(weight[i]))
			{
				// get average
				int n=0, sum=0;
				for (size_t j=0; j < n_samp; j++)
					if (R_FINITE(s[j])) { n++; sum += s[j]; }
				// imputed by the mean
				double *p = out_ds, w = weight[i], m = (double)sum / n;
				if (sum <= n)
				{
					for (size_t j=0; j < n_samp; j++)
						p[j] += R_FINITE(s[j]) ? (s[j] * w) : (m * w);
				} else {
					// flip
					m = 2 - m;
					for (size_t j=0; j < n_samp; j++)
						p[j] += R_FINITE(s[j]) ? ((2-s[j]) * w) : (m * w);
				}
			}
		}
		break;
	default:
		if (!Rf_isNull(ds_mat)) Rf_error(ERR_DS_TYPE);
	}
}

// ====================================

/// calculate burden p-values for binary outcomes
RcppExport SEXP saige_burden_test_bin(SEXP dosage)
{
BEGIN_RCPP
	// get # of SNPs
	const int n_snp = ds_mat_dim(dosage);
	// get MAF and MAC
	double *maf = buf_unitsz;
	double *mac = buf_unitsz + num_unitsz;
	double *ws  = buf_unitsz + 2*num_unitsz;
	ds_mat_mafmac(maf, mac, dosage, mod_NSamp, n_snp);

	// summarize maf & mac
	NumericVector ans(8 + 6*num_wbeta);
	f64_mean_sd(maf, n_snp, ans[0], ans[1]);
	f64_maxmin(maf, n_snp, ans[3], ans[2]);
	f64_mean_sd(mac, n_snp, ans[4], ans[5]);
	f64_maxmin(mac, n_snp, ans[7], ans[6]);

	// for each beta weight
	for (int i=0; i < num_wbeta; i++)
	{
		// get weights
		const double b1 = buf_wbeta[2*i + 0];
		const double b2 = buf_wbeta[2*i + 1];
		for (int j=0; j < n_snp; j++)
			ws[j] = Rf_dbeta(maf[j], b1, b2, FALSE);
		f64_normalize(n_snp, ws);

		// burden with weights
		double *G = buf_dosage;
		ds_mat_burden(G, dosage, mod_NSamp, n_snp, ws);
		double summac = f64_sum(mod_NSamp, G) * n_snp;  // sum of mac

		bool converged = 0;
		double AF, mac, beta, SE, pval, pval_noadj;
		AF = mac = beta = SE = pval = pval_noadj = R_NaN;

		if (summac >= threshold_summac && summac > 0)
		{
			// single-variant calculation
			int num = 0;
			single_test_bin(mod_NSamp, G, AF, mac, num, beta, SE, pval,
				pval_noadj, converged);
		}

		// set the output
		const int st = 8 + 6*i;
		ans[0+st] = summac; ans[1+st] = beta;   ans[2+st] = SE;
		ans[3+st] = pval;   ans[4+st] = pval_noadj;
		ans[5+st] = converged ? 1 : 0;
	}

	// output
	return ans;
END_RCPP
}

/// calculate burden p-values for quantitative outcomes
RcppExport SEXP saige_burden_test_quant(SEXP dosage)
{
BEGIN_RCPP
	// get # of SNPs
	const int n_snp = ds_mat_dim(dosage);
	// get MAF and MAC
	double *maf = buf_unitsz;
	double *mac = buf_unitsz + num_unitsz;
	double *ws  = buf_unitsz + 2*num_unitsz;
	ds_mat_mafmac(maf, mac, dosage, mod_NSamp, n_snp);

	// summarize maf & mac
	NumericVector ans(8 + 4*num_wbeta);
	f64_mean_sd(maf, n_snp, ans[0], ans[1]);
	f64_maxmin(maf, n_snp, ans[3], ans[2]);
	f64_mean_sd(mac, n_snp, ans[4], ans[5]);
	f64_maxmin(mac, n_snp, ans[7], ans[6]);

	// for each beta weight
	for (int i=0; i < num_wbeta; i++)
	{
		// get weights
		const double b1 = buf_wbeta[2*i + 0];
		const double b2 = buf_wbeta[2*i + 1];
		for (int j=0; j < n_snp; j++)
			ws[j] = Rf_dbeta(maf[j], b1, b2, FALSE);
		f64_normalize(n_snp, ws);

		// burden with weights
		double *G = buf_dosage;
		ds_mat_burden(G, dosage, mod_NSamp, n_snp, ws);
		double summac = f64_sum(mod_NSamp, G) * n_snp;  // sum of mac

		double AF, mac, beta, SE, pval;
		AF = mac = beta = SE = pval = R_NaN;

		if (summac >= threshold_summac && summac > 0)
		{
			// single-variant calculation
			int num;
			single_test_quant(mod_NSamp, G, AF, mac, num, beta, SE, pval);
		}

		// set the output
		const int st = 8 + 4*i;
		ans[0+st] = summac;  ans[1+st] = beta;
		ans[2+st] = SE;      ans[3+st] = pval;
	}

	// output
	return ans;
END_RCPP
}


// ====================================

static double acat_pval(R_xlen_t n, const double pval[], const double w[],
	bool throw_error);

/// calculate ACAT-V p-values for binary outcomes
RcppExport SEXP saige_acatv_test_bin(SEXP dosage)
{
BEGIN_RCPP
	// get # of SNPs
	const int n_snp = ds_mat_dim(dosage);
	// get MAF and MAC
	double *maf = buf_unitsz;
	double *mac = buf_unitsz + num_unitsz;
	ds_mat_mafmac(maf, mac, dosage, mod_NSamp, n_snp);

	// summarize maf & mac
	NumericVector ans(10 + num_wbeta*4);
	f64_mean_sd(maf, n_snp, ans[0], ans[1]);
	f64_maxmin(maf, n_snp, ans[3], ans[2]);
	f64_mean_sd(mac, n_snp, ans[4], ans[5]);
	f64_maxmin(mac, n_snp, ans[7], ans[6]);

	// buffer
	double *w_burden = buf_unitsz + 2*num_unitsz;
	double *w_pval = buf_unitsz + 3*num_unitsz;
	double *pvals  = buf_unitsz + 4*num_unitsz;

	// for each beta weight
	for (int w_i=0; w_i < num_wbeta; w_i++)
	{
		// get weights
		const double b1 = buf_wbeta[2*w_i + 0];
		const double b2 = buf_wbeta[2*w_i + 1];
		int n_single = 0;  // # of SNPs for single variant tests
		int n_burden = 0;  // # of SNPs for burden test
		double summaf_burden = 0;  // sum of MAF for burden test
		// for each SNP
		for (R_xlen_t j=0; j < n_snp; j++)
		{
			if (mac[j] >= threshold_acatv_mac)
			{
				// single variant tests
				double *G = get_ds(dosage, mod_NSamp, mod_NSamp*j);
				// run
				int num = 0;
				bool converged = 0;
				double AF, mac, beta, SE, pval, pval_noadj;
				AF = mac = beta = SE = pval = pval_noadj = R_NaN;
				single_test_bin(mod_NSamp, G, AF, mac, num, beta, SE,
					pval, pval_noadj, converged);
				// save
				const double p = maf[j];
				w_pval[n_single] = sq(Rf_dbeta(p, b1, b2, FALSE)) * p * (1-p);
				pvals[n_single] = pval;  // p-value for the SNP
				n_single ++;
				w_burden[j] = R_NaN;  // no weight for burden
			} else {
				// burden test for very rare variants
				n_burden ++;
				summaf_burden += maf[j];
				w_burden[j] = Rf_dbeta(maf[j], b1, b2, FALSE);
			}
		}
		// if has burden
		if (n_burden > 0)
		{
			f64_normalize(n_snp, w_burden);
			// burden with weights
			double *G = buf_dosage;
			ds_mat_burden(G, dosage, mod_NSamp, n_snp, w_burden);
			double summac = f64_sum(mod_NSamp, G) * n_snp;  // sum of mac
			if (summac >= threshold_summac && summac > 0)
			{
				bool converged = 0;
				double AF, mac, beta, SE, pval, pval_noadj;
				AF = mac = beta = SE = pval = pval_noadj = R_NaN;
				int num = 0;
				single_test_bin(mod_NSamp, G, AF, mac, num, beta, SE, pval,
					pval_noadj, converged);
				if (R_FINITE(pval))
				{
					const double p = summaf_burden / n_burden;
					w_pval[n_single] = sq(Rf_dbeta(p, b1, b2, FALSE)) * p * (1-p);
					pvals[n_single] = pval;
					n_single ++;
				}
			}
		}

		// set the output
		if (w_i == 0)
		{
			ans[8] = n_single - n_burden;
			ans[9] = n_burden;
		}
		const int st = 10 + w_i * 4;
		ans[st + 0] =
			(n_single > 0) ? acat_pval(n_single, pvals, w_pval, false) : R_NaN;
		f64_medmaxmin(pvals, n_single, ans[st+1], ans[st+3], ans[st+2]);
	}

	// output
	return ans;
END_RCPP
}

/// calculate burden p-values for quantitative outcomes
RcppExport SEXP saige_acatv_test_quant(SEXP dosage)
{
BEGIN_RCPP
	Rf_error("'saige_acatv_test_quant' not implemented.");
	return R_NilValue;
END_RCPP
}


// ====================================

/// calculate ACAT-O p-values for binary outcomes
RcppExport SEXP saige_acato_test_bin(SEXP dosage)
{
BEGIN_RCPP

	// get # of SNPs
	const int n_snp = ds_mat_dim(dosage);
	// get MAF and MAC
	double *maf = buf_unitsz;
	double *mac = buf_unitsz + num_unitsz;
	ds_mat_mafmac(maf, mac, dosage, mod_NSamp, n_snp);

	// summarize maf & mac
	NumericVector ans(8 + 1 + 2*num_wbeta);
	f64_mean_sd(maf, n_snp, ans[0], ans[1]);
	f64_maxmin(maf, n_snp, ans[3], ans[2]);
	f64_mean_sd(mac, n_snp, ans[4], ans[5]);
	f64_maxmin(mac, n_snp, ans[7], ans[6]);

	// buffer
	double *w_burden = buf_unitsz + 2*num_unitsz;
	double *w_pval = buf_unitsz + 3*num_unitsz;
	double *pvals  = buf_unitsz + 4*num_unitsz;

	// burden tests
	for (int w_i=0; w_i < num_wbeta; w_i++)
	{
		// get weights
		const double b1 = buf_wbeta[2*w_i + 0];
		const double b2 = buf_wbeta[2*w_i + 1];
		for (int j=0; j < n_snp; j++)
			w_burden[j] = Rf_dbeta(maf[j], b1, b2, FALSE);
		f64_normalize(n_snp, w_burden);

		// burden with weights
		double *G = buf_dosage;
		ds_mat_burden(G, dosage, mod_NSamp, n_snp, w_burden);
		double summac = f64_sum(mod_NSamp, G) * n_snp;  // sum of mac

		double pval = R_NaN;
		if (summac >= threshold_summac && summac > 0)
		{
			// single-variant calculation
			int num = 0;
			bool converged = 0;
			double AF, mac, beta, SE, pval_noadj;
			single_test_bin(mod_NSamp, G, AF, mac, num, beta, SE, pval,
				pval_noadj, converged);
		}
		ans[9 + 2*w_i] = pval;
	}

	// ACAT-V tests
	for (int w_i=0; w_i < num_wbeta; w_i++)
	{
		// get weights
		const double b1 = buf_wbeta[2*w_i + 0];
		const double b2 = buf_wbeta[2*w_i + 1];
		int n_single = 0;  // # of SNPs for single variant tests
		int n_burden = 0;  // # of SNPs for burden test
		double summaf_burden = 0;  // sum of MAF for burden test
		// for each SNP
		for (R_xlen_t j=0; j < n_snp; j++)
		{
			if (mac[j] >= threshold_acatv_mac)
			{
				// single variant tests
				double *G = get_ds(dosage, mod_NSamp, mod_NSamp*j);
				// run
				int num = 0;
				bool converged = 0;
				double AF, mac, beta, SE, pval, pval_noadj;
				AF = mac = beta = SE = pval = pval_noadj = R_NaN;
				single_test_bin(mod_NSamp, G, AF, mac, num, beta, SE,
					pval, pval_noadj, converged);
				// save
				const double p = maf[j];
				w_pval[n_single] = sq(Rf_dbeta(p, b1, b2, FALSE)) * p * (1-p);
				pvals[n_single] = pval;  // p-value for the SNP
				n_single ++;
				w_burden[j] = R_NaN;  // no weight for burden
			} else {
				// burden test for very rare variants
				n_burden ++;
				summaf_burden += maf[j];
				w_burden[j] = Rf_dbeta(maf[j], b1, b2, FALSE);
			}
		}

		// if has burden
		if (n_burden > 0)
		{
			f64_normalize(n_snp, w_burden);
			// burden with weights
			double *G = buf_dosage;
			ds_mat_burden(G, dosage, mod_NSamp, n_snp, w_burden);
			double summac = f64_sum(mod_NSamp, G) * n_snp;  // sum of mac
			if (summac >= threshold_summac && summac > 0)
			{
				bool converged = 0;
				double AF, mac, beta, SE, pval, pval_noadj;
				AF = mac = beta = SE = pval = pval_noadj = R_NaN;
				int num = 0;
				single_test_bin(mod_NSamp, G, AF, mac, num, beta, SE, pval,
					pval_noadj, converged);
				if (R_FINITE(pval))
				{
					const double p = summaf_burden / n_burden;
					w_pval[n_single] = sq(Rf_dbeta(p, b1, b2, FALSE)) * p * (1-p);
					pvals[n_single] = pval;
					n_single ++;
				}
			}
		}
		ans[10 + 2*w_i] =
			(n_single > 0) ? acat_pval(n_single, pvals, w_pval, false) : R_NaN;
	}

	// combined p-value
	const int n_pval = 2*num_wbeta;
	double *ws = buf_unitsz;
	if (n_pval > 5*num_unitsz)
		ws = REAL(NEW_NUMERIC(n_pval));
	for (int i=0; i < n_pval; i++) ws[i] = 1;
	ans[8] = acat_pval(n_pval, &ans[9], ws, false);

	// output
	return ans;
END_RCPP
}

/// calculate burden p-values for quantitative outcomes
RcppExport SEXP saige_acato_test_quant(SEXP dosage)
{
BEGIN_RCPP
	Rf_error("'saige_acato_test_quant' not implemented.");
	return R_NilValue;
END_RCPP
}



// ========================================================================= //

#ifdef HAVE_ATANPI
extern double atanpi(double x);
#else
inline static double atanpi(double x) { return atan(x) / M_PI; }
#endif

static const double ROUND_ZERO = 1e-300;
static const double ROUND_ONE  = 1 - 1e-16;

/// p-value from ACAT combination method
static double acat_pval(R_xlen_t n, const double pval[], const double w[],
	bool throw_error)
{
	// get the weight sum
	double sumw = 0;
	for (R_xlen_t i=0; i < n; i++)
		if (R_FINITE(pval[i]) && R_FINITE(w[i])) sumw += w[i];
	if (sumw <= 0)
	{
		if (throw_error)
			Rf_error("the sum of weights should be > 0.");
		else
			return R_NaN;
	}
	// get statistic
	double Tstat = 0;
	for (R_xlen_t i=0; i < n; i++)
	{
		double p = pval[i];
		if (R_FINITE(p) && R_FINITE(w[i]))
		{
			// check p-value
			if (p < 0 || p > 1)
			{
				if (throw_error)
					Rf_error("Invalid input p-value: %g.", p);
				else
					return R_NaN;
			}
			if (p < ROUND_ZERO)
				p = ROUND_ZERO;  // almost the smallest number > 0
			else if (p > ROUND_ONE)
				p = ROUND_ONE;  // almost the closest number around 1
			// calc stat
			if (p >= 1e-15)
			{
				Tstat += w[i] * tanpi(0.5 - p);
			} else {
				// based on the taylor series expansion,
				//   Series[tan[(1/2-x)*pi], {x, 0, 5}]
				//   1/(pi*x) - pi/3*x - pi^3/45*x^3 - 2*pi^5/945*x^5 + O(x^7)
				Tstat += w[i] / p / M_PI;
			}
		}
	}
	Tstat /= sumw;
	// get p-value from Tstat, and return
	if (Tstat <= 5e+14)
		return 0.5 - atanpi(Tstat);
	else
		return 1.0 / Tstat * M_1_PI;
}

RcppExport SEXP saige_acat_p(SEXP pval, SEXP weight)
{
	const R_xlen_t n = Rf_xlength(pval);
	// check pval
	if (n <= 0)
		Rf_error("the number of p-values should be > 0.");
	else if (n == 1)
		return pval;
	// check weight
	if (Rf_isNull(weight))
	{
		weight = NEW_NUMERIC(n);
		double *w = REAL(weight);
		for (R_xlen_t i=0; i < n; i++) w[i] = 1;
	}
	// check
	if (n != Rf_xlength(weight))
		Rf_error("weights should have the same length as p-values.");
	if (TYPEOF(pval) != REALSXP)
		Rf_error("p-values should be numeric.");
	if (TYPEOF(weight) != REALSXP)
		Rf_error("weights should be numeric.");
	// calculate
	double v = acat_pval(n, REAL(pval), REAL(weight), true);
	// output
	return Rf_ScalarReal(v);
}


// ========================================================================= //

RcppExport SEXP saige_simd_version();
RcppExport SEXP saige_store_2b_geno(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP saige_store_sp_geno(SEXP, SEXP, SEXP, SEXP, SEXP);

/// initialize the package
RcppExport void R_init_SAIGEgds(DllInfo *info)
{
	#define CALL(name, num)	   { #name, (DL_FUNC)&name, num }

	static R_CallMethodDef callMethods[] =
	{
		CALL(saige_score_test_init, 1),
		CALL(saige_simd_version, 0),
		CALL(saige_store_2b_geno, 5),
		CALL(saige_store_sp_geno, 5),
		CALL(saige_burden_test_bin, 1),
		CALL(saige_burden_test_quant, 1),
		CALL(saige_acatv_test_bin, 1),
		CALL(saige_acato_test_bin, 1),
		CALL(saige_acat_p, 2),
		{ NULL, NULL, 0 }
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
