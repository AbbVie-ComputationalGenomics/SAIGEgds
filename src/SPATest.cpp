// ===========================================================
//
// SPATest.cpp: C implementation of part of the R SPATest package
//
// Copyright (C) 2019    Xiuwen Zheng / AbbVie-ComputationalGenomics
//
// This file is part of SAIGEgds. It was created based on the R codes in the
// SPAtest package.
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
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>
#include <cfloat>
#include <cmath>


inline static double sq(double v) { return v*v; }

// inline static int sign(double v) { return (v>0) ? 1 : ((v<0) ? -1 : 0); }


inline static COREARRAY_TARGET_CLONES
	double Korg(double t, size_t n_g, const double mu[], const double g[])
{
	double sum = 0;
	for (size_t i=0; i < n_g; i++)
	{
		double m_i = mu[i];
		sum += log(1 - m_i + m_i * exp(g[i] * t));
	}
	return sum;
}


inline static COREARRAY_TARGET_CLONES
	double K1_adj(double t, size_t n_g, const double mu[], const double g[], double q)
{
	double sum = 0;
	for (size_t i=0; i < n_g; i++)
	{
		double m_i=mu[i], g_i=g[i];
		sum += m_i * g_i / ((1-m_i) * exp(-g_i * t) + m_i);
	}
	return sum - q;
}


inline static COREARRAY_TARGET_CLONES
	double K2(double t, size_t n_g, const double mu[], const double g[])
{
	double sum = 0;
	for (size_t i=0; i < n_g; i++)
	{
		double m_i=mu[i], one_m_i=1-m_i;
		double g_i=g[i], exp_i=exp(-g_i*t);
		double v = (one_m_i * m_i * g_i * g_i * exp_i) / sq(one_m_i * exp_i + m_i);
		if (R_FINITE(v)) sum += v;
	}
	return sum;
}


// .Machine$double.eps^0.25
static const double root_tol = sqrt(sqrt(DBL_EPSILON));
static const int MaxNumIter = 1000;


inline static void COREARRAY_TARGET_CLONES
	getroot_K1(const double g_pos, const double g_neg,
		double &root, int &n_iter, bool &converged, double init, size_t n_g,
		const double mu[], const double g[], double q, double tol=root_tol,
		int maxiter=MaxNumIter)
{
	if (q>=g_pos || q<=g_neg)
	{
		root = R_PosInf; n_iter = 0;
		converged = true;

	} else {
		double t = root = init;
		double K1_eval = K1_adj(t, n_g, mu, g, q);
		double prevJump = R_PosInf;
		converged = false;
		for (int n_iter=1; n_iter <= maxiter; n_iter++)
		{
			double K2_eval = K2(t, n_g, mu, g);
			double tnew = t - K1_eval/K2_eval;
			if (!R_FINITE(tnew))
				break;
			if (fabs(tnew - t) < tol)
			{
				converged = true;
				break;
			}
			double newK1 = K1_adj(tnew, n_g, mu, g, q);
			if (sign(K1_eval) != sign(newK1))
			{
				if (fabs(tnew - t) > prevJump-tol)
				{
					tnew = t + sign(newK1 - K1_eval) * prevJump * 0.5;
					newK1 = K1_adj(tnew, n_g, mu, g, q);
					prevJump *= 0.5;
				} else {
					prevJump = fabs(tnew - t);
				}
			}
			root = t = tnew;
			K1_eval = newK1;
		}
	}
}

inline static void COREARRAY_TARGET_CLONES
	getroot_K1_fast(const double g_pos, const double g_neg,
		double &root, int &n_iter, bool &converged, double init, size_t n_nonzero,
		const double mu[], const double g[], double q, double NAmu, double NAsigma,
		double tol=root_tol, int maxiter=MaxNumIter)
{
	if (q>=g_pos || q<=g_neg)
	{
		root = R_PosInf; n_iter = 0;
		converged = true;

	} else {
		double t = root = init;
		double K1_eval = K1_adj(t, n_nonzero, mu, g, q) + NAmu + NAsigma * t;
		double prevJump = R_PosInf;
		converged = false;
		for (int n_iter=1; n_iter <= maxiter; n_iter++)
		{
			double K2_eval = K2(t, n_nonzero, mu, g) + NAsigma;
			double tnew = t - K1_eval/K2_eval;
			if (!R_FINITE(tnew))
				break;
			if (fabs(tnew - t) < tol)
			{
				converged = true;
				break;
			}
			double newK1 = K1_adj(tnew, n_nonzero, mu, g, q) + NAmu + NAsigma * tnew;
			if (sign(K1_eval) != sign(newK1))
			{
				if (fabs(tnew - t) > prevJump-tol)
				{
					tnew = t + sign(newK1 - K1_eval) * prevJump * 0.5;
					newK1 = K1_adj(tnew, n_nonzero, mu, g, q) + NAmu + NAsigma * tnew;
					prevJump *= 0.5;
				} else {
					prevJump = fabs(tnew - t);
				}
			}
			root = t = tnew;
			K1_eval = newK1;
		}
	}
}


inline static double COREARRAY_TARGET_CLONES
	Get_Saddle_Prob(double zeta, size_t n_g, const double mu[], const double g[],
		double q)
{
	double k1 = Korg(zeta, n_g, mu, g);
	double k2 = K2(zeta, n_g, mu, g);
	double pval = 0;
	if (R_FINITE(k1) && R_FINITE(k1))
	{
		double temp1 = zeta * q - k1;
		double w = sign(zeta) * sqrt(2 * temp1);
		double v = zeta * sqrt(k2);
		double Z_test = w + 1/w * log(v/w);
		if (Z_test > 0)
			pval = ::Rf_pnorm5(Z_test, 0, 1, FALSE, FALSE);
		else
			pval = - ::Rf_pnorm5(Z_test, 0, 1, TRUE, FALSE);
	}
	return pval;
}

inline static double COREARRAY_TARGET_CLONES
	Get_Saddle_Prob_fast(double zeta, size_t n_nonzero, const double mu[],
		const double g[], double q, double NAmu, double NAsigma)
{
	double k1 = Korg(zeta, n_nonzero, mu, g) + NAmu * zeta + 0.5 * NAsigma * zeta * zeta;
	double k2 = K2(zeta, n_nonzero, mu, g) + NAsigma;
	double pval = 0;
	if (R_FINITE(k1) && R_FINITE(k1))
	{
		double temp1 = zeta * q - k1;
		double w = sign(zeta) * sqrt(2 * temp1);
		double v = zeta * sqrt(k2);
		double Z_test = w + 1/w * log(v/w);
		if (Z_test > 0)
			pval = ::Rf_pnorm5(Z_test, 0, 1, FALSE, FALSE);
		else
			pval = - ::Rf_pnorm5(Z_test, 0, 1, TRUE, FALSE);
	}
	return pval;
}



// m1 <- sum(mu * g),  var1 <- sum(mu * (1-mu) * g^2)
extern "C" double COREARRAY_TARGET_CLONES
	Saddle_Prob(double q, double m1, double var1, size_t n_g, const double mu[],
		const double g[], double cutoff, bool &converged)
{
	double s = q - m1;
	double qinv = -s + m1;
	double pval_noadj = Rf_pchisq(s*s/var1, 1, FALSE, FALSE);
	double pval;
	double g_pos=0, g_neg=0;
	double init=false;

	while (true)
	{
		converged = true;
		if (cutoff < 0.1) cutoff = 0.1;

		if (fabs(q - m1)/sqrt(var1) < cutoff)
		{
			pval = pval_noadj;
		} else {
			// need initializing
			if (!init)
			{
				init = true;
				for (size_t i=0; i < n_g; i++)
				{
					double v = g[i];
					if (v > 0) g_pos += v; else g_neg += v;
				}
			}
			//
			double uni1_root, uni2_root;
			int n_iter1, n_iter2;
			bool conv1, conv2;
			getroot_K1(g_pos, g_neg, uni1_root, n_iter1, conv1, 0, n_g, mu, g, q);
			getroot_K1(g_pos, g_neg, uni2_root, n_iter2, conv2, 0, n_g, mu, g, qinv);
			if (conv1 && conv2)
			{
				double p1 = Get_Saddle_Prob(uni1_root, n_g, mu, g, q);
				double p2 = Get_Saddle_Prob(uni2_root, n_g, mu, g, qinv);
				pval = fabs(p1) + fabs(p2);
			} else {
				pval = pval_noadj;
				converged = false;
				break;
			}
		}

		if (pval!=0 && pval_noadj/pval>1000)
			cutoff *= 2;
		else
			break;
	}
	return pval;
}


// m1 <- sum(mu * g),  var1 <- sum(mu * (1-mu) * g^2)
extern "C" double COREARRAY_TARGET_CLONES
	Saddle_Prob_Fast(double q, double m1, double var1, size_t n_g, const double mu[],
		const double g[], size_t n_nonzero, const int nonzero_idx[], double cutoff,
		bool &converged, double buf_spa[])
{
	double s = q - m1;
	double qinv = -s + m1;
	double pval_noadj = Rf_pchisq(s*s/var1, 1, FALSE, FALSE);
	double pval;
	double NAmu, NAsigma;
	double g_pos=0, g_neg=0;
	double init=false;

	while (true)
	{
		converged = true;
		if (cutoff < 0.1) cutoff = 0.1;

		if (fabs(q - m1)/sqrt(var1) < cutoff)
		{
			pval = pval_noadj;
		} else {
			// need initializing
			if (!init)
			{
				init = true;
				// get g_pos, g_neg
				for (size_t i=0; i < n_g; i++)
				{
					double v = g[i];
					if (v > 0) g_pos += v; else g_neg += v;
				}
				// re-save mu and g to buf_spa, calculate NAmu and NAsigma
				NAmu = m1, NAsigma = var1;
				for (size_t i=0; i < n_nonzero; i++)
				{
					size_t k = nonzero_idx[i];
					double g_k, mu_k;
					buf_spa[i] = g_k = g[k];
					buf_spa[i + n_nonzero] = mu_k = mu[k];
					NAmu -= g_k * mu_k;
					NAsigma -= g_k * g_k * mu_k * (1 - mu_k);
				}
				g = &buf_spa[0]; mu = &buf_spa[n_nonzero];
			}
			//
			double uni1_root, uni2_root;
			int n_iter1, n_iter2;
			bool conv1, conv2;
			getroot_K1_fast(g_pos, g_neg, uni1_root, n_iter1, conv1, 0, n_nonzero,
				mu, g, q, NAmu, NAsigma);
			getroot_K1_fast(g_pos, g_neg, uni2_root, n_iter2, conv2, 0, n_nonzero,
				mu, g, qinv, NAmu, NAsigma);
			if (conv1 && conv2)
			{
				double p1 = Get_Saddle_Prob_fast(uni1_root, n_nonzero, mu, g, q,
					NAmu, NAsigma);
				double p2 = Get_Saddle_Prob_fast(uni2_root, n_nonzero, mu, g, qinv,
					NAmu, NAsigma);
				pval = fabs(p1) + fabs(p2);
			} else {
				pval = pval_noadj;
				converged = false;
				break;
			}
		}

		if (pval!=0 && pval_noadj/pval>1000)
			cutoff *= 2;
		else
			break;
	}
	return pval;
}
