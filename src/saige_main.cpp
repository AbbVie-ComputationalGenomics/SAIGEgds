// ===========================================================
//
// saige_main.cpp: SAIGE association analysis
//
// Copyright (C) 2019-2020    Xiuwen Zheng / AbbVie-ComputationalGenomics
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

#include "vectorization.h"
#include <Rcpp.h>
#include <algorithm>


using namespace Rcpp;
using namespace vectorization;



// ========================================================================= //
// internal functions

/// SPAtest
extern "C" double Saddle_Prob(double q, double m1, double var1, size_t n_g,
	const double mu[], const double g[], double cutoff, bool &converged);

extern "C" double Saddle_Prob_Fast(double q, double m1, double var1, size_t n_g,
	const double mu[], const double g[], size_t n_nonzero, const int nonzero_idx[],
	double cutoff, bool &converged, double buf_spa[]);

/// square
inline double sq(double v) { return v*v; }


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

static double *buf_dosage = NULL;    //< temporary buffer for real dosages if the input is not real numbers
static double *buf_coeff = NULL;
static double *buf_adj_g = NULL;     //< genotype after adjusting for fixed effects
static int *buf_index = NULL;
static double *buf_B = NULL;
static double *buf_g_tilde = NULL;   //<
static double *buf_spa = NULL;       //< the buffer for SPA calculation
static double *buf_tmp = NULL;

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
	buf_spa = REAL(M["buf_spa"]);
	buf_tmp = REAL(M["buf_tmp"]);
END_RCPP
}



// ========================================================================= //

static double *get_real_dosage(SEXP dosage, size_t &out_n_samp)
{
	const size_t num_samp = Rf_length(dosage);
	if (size_t(mod_NSamp) != num_samp)
		throw std::invalid_argument("Invalid length of dosages.");
	out_n_samp = num_samp;

	int *I;
	unsigned char *p;
	const double NaN = R_NaN;

	switch (TYPEOF(dosage))
	{
		case REALSXP:
			return REAL(dosage);
		case INTSXP:
			I = INTEGER(dosage);
			for (size_t i=0; i < num_samp; i++)
				buf_dosage[i] = (I[i] != NA_INTEGER) ? I[i] : NaN;
			return buf_dosage;
		case RAWSXP:
			p = (unsigned char *)RAW(dosage);
			for (size_t i=0; i < num_samp; i++)
				buf_dosage[i] = (p[i] != 0xFF) ? p[i] : NaN;
			return buf_dosage;
	}
	throw std::invalid_argument(
		"Invalid data type for dosages (should be one of RAW, INT or REAL).");
}


/// calculate p-values for quantitative outcome
RcppExport COREARRAY_TARGET_CLONES SEXP saige_score_test_quant(SEXP dosage)
{
BEGIN_RCPP

	// dosages and imputed
	size_t num_samp = Rf_length(dosage);
	double *G = get_real_dosage(dosage, num_samp);

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
			// buf_tmp = t(X1) * (y-mu)
			f64_mul_mat_vec_sp(n_nonzero, buf_index, mod_NCoeff, mod_t_X,
				mod_y_mu, buf_tmp);
			// S2 = sum((buf_tmp - mod_S_a) .* buf_coeff)
			double S2 = 0;
			for (int i=0; i < mod_NCoeff; i++)
				S2 += (buf_tmp[i] - mod_S_a[i]) * buf_coeff[i];
			//
			double Tstat = (S1 + S2) * inv_sqrt_mac / mod_tau[0];
			pval = ::Rf_pchisq(Tstat*Tstat/var1, 1, FALSE, FALSE);
			beta = Tstat / var1 * inv_sqrt_mac;

		} else {
			// adj_g = G - XXVX_inv * (XV * G), adjusted genotypes
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

		NumericVector ans(6);
		ans[0] = AF;    ans[1] = mac;   ans[2] = Num;
		ans[3] = beta;  ans[4] = SE;    ans[5] = pval;
		return ans;
	} else {
		return R_NilValue;
	}

END_RCPP
}


/// calculate p-values for binary outcome
RcppExport COREARRAY_TARGET_CLONES SEXP saige_score_test_bin(SEXP dosage)
{
BEGIN_RCPP

	// dosages and imputed
	size_t num_samp = Rf_length(dosage);
	double *G = get_real_dosage(dosage, num_samp);

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
			// buf_tmp = t(X1) * (y-mu)
			f64_mul_mat_vec_sp(n_nonzero, buf_index, mod_NCoeff, mod_t_X,
				mod_y_mu, buf_tmp);
			// S2 = sum((buf_tmp - mod_S_a) .* buf_coeff)
			double S2 = 0;
			for (int i=0; i < mod_NCoeff; i++)
				S2 += (buf_tmp[i] - mod_S_a[i]) * buf_coeff[i];
			//
			double S = S1 + S2;
			pval_noadj = ::Rf_pchisq(S*S/var1, 1, FALSE, FALSE);
			beta = S / var1;

		} else {
			// adj_g = G - XXVX_inv * (XV * G), adjusted genotypes
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
				buf_adj_g, n_nonzero, buf_index, 2, converged, buf_spa);
			if (pval==0 && pval_noadj>0)
				{ pval = pval_noadj; converged = false; }

			// effect size
			beta = (Tstat / var1) / sqrt(AC2);
		}

		if (minus) beta = -beta;
		double SE = fabs(beta / ::Rf_qnorm5(pval/2, 0, 1, TRUE, FALSE));

		NumericVector ans(8);
		ans[0] = AF;    ans[1] = mac;   ans[2] = Num;
		ans[3] = beta;  ans[4] = SE;    ans[5] = pval;
		ans[6] = pval_noadj;
		ans[7] = converged ? 1 : 0;
		return ans;
	} else {
		return R_NilValue;
	}

END_RCPP
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
		{ NULL, NULL, 0 }
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
