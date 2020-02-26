// ===========================================================
//
// SPATest.cpp: C implementation of part of the R SPATest package
//
// Copyright (C) 2019-2020    Xiuwen Zheng / AbbVie-ComputationalGenomics
//
// This file is part of SAIGEgds. It was created based on the R codes in the
// SPAtest package with the reference:
//     A Fast and Accurate Algorithm to Test for Binary Phenotypes and Its
//     Application to PheWAS. Dey R, Schmidt EM, Abecasis GR, Lee S.
//     Am J Hum Genet. 2017 Jul 6;101(1):37-49. doi:10.1016/j.ajhg.2017.05.014.
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


// Cumulant-generating function (CGF) of the score statistic: K(t)

/// The log term of CGF
inline static COREARRAY_TARGET_CLONES MATH_OFAST  // auto-vectorize if possible
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


/// The first term of the first-order derivative of CGF
inline static COREARRAY_TARGET_CLONES MATH_OFAST  // auto-vectorize if possible
	double K1_adj(double t, size_t n_g, const double mu[], const double g[],
		double q)
{
	double sum = 0;
	for (size_t i=0; i < n_g; i++)
	{
		double m_i=mu[i], g_i=g[i];
		sum += m_i * g_i / ((1-m_i) * exp(-g_i * t) + m_i);
	}
	return sum - q;
}


/// The second-order derivative of CGF
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


/// Root-finding algorithm for a full saddle-point method
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


/// Root-finding algorithm for a partially normal approximation approach
inline static void COREARRAY_TARGET_CLONES
	getroot_K1_fast(const double g_pos, const double g_neg, double &root,
		int &n_iter, bool &converged, double init, size_t n_nonzero,
		const double mu[], const double g[], double q, double NAmu,
		double NAsigma, double tol=root_tol, int maxiter=MaxNumIter)
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
			double newK1 = K1_adj(tnew, n_nonzero, mu, g, q) +
				NAmu + NAsigma * tnew;
			if (sign(K1_eval) != sign(newK1))
			{
				if (fabs(tnew - t) > prevJump-tol)
				{
					tnew = t + sign(newK1 - K1_eval) * prevJump * 0.5;
					newK1 = K1_adj(tnew, n_nonzero, mu, g, q) +
						NAmu + NAsigma * tnew;
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


/// Get a p-value from a full saddle-point method
inline static double COREARRAY_TARGET_CLONES
	get_saddle_prob(double t, size_t n_g, const double mu[], const double g[],
		double q)
{
	if (!R_FINITE(t)) return 0;
	double K  = Korg(t, n_g, mu, g);
	double k2 = K2(t, n_g, mu, g);
	double pval = 0;
	if (R_FINITE(K) && R_FINITE(k2))
	{
		double w = sign(t) * sqrt(2 * (t * q - K));
		double v = t * sqrt(k2);
		double z = w + log(v/w) / w;
		if (z > 0)
			pval = ::Rf_pnorm5(z, 0, 1, FALSE, FALSE);
		else
			pval = - ::Rf_pnorm5(z, 0, 1, TRUE, FALSE);
	}
	return pval;
}


/// Get a p-value from a partially normal approximation approach
inline static double COREARRAY_TARGET_CLONES
	get_saddle_prob_fast(double t, size_t n_nonzero, const double mu[],
		const double g[], double q, double NAmu, double NAsigma)
{
	if (!R_FINITE(t)) return 0;
	double K  = Korg(t, n_nonzero, mu, g) + NAmu * t + 0.5 * NAsigma * t * t;
	double k2 = K2(t, n_nonzero, mu, g) + NAsigma;
	double pval = 0;
	if (R_FINITE(K) && R_FINITE(k2))
	{
		double w = sign(t) * sqrt(2 * (t * q - K));
		double v = t * sqrt(k2);
		double z = w + log(v/w) / w;
		if (z > 0)
			pval = ::Rf_pnorm5(z, 0, 1, FALSE, FALSE);
		else
			pval = - ::Rf_pnorm5(z, 0, 1, TRUE, FALSE);
	}
	return pval;
}



/// Get a p-value from the CGF of a full saddle-point method
/// Input: m1 <- sum(mu * g), var1 <- sum(mu * (1-mu) * g^2)
/// Output: p-value
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
			double root1, root2;  // the roots of equation
			int ni1, ni2;  // the number of iterations
			bool conv1, conv2;  // whether the algorithm converges or not
			getroot_K1(g_pos, g_neg, root1, ni1, conv1, 0, n_g, mu, g, q);
			getroot_K1(g_pos, g_neg, root2, ni2, conv2, 0, n_g, mu, g, qinv);
			if (conv1 && conv2)
			{
				double p1 = get_saddle_prob(root1, n_g, mu, g, q);
				double p2 = get_saddle_prob(root2, n_g, mu, g, qinv);
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


/// Get a p-value from faster calculation of the CGF by a partially
///     normal approximation approach
/// Input: m1 <- sum(mu * g), var1 <- sum(mu * (1-mu) * g^2)
/// Output: p-value
extern "C" double COREARRAY_TARGET_CLONES
	Saddle_Prob_Fast(double q, double m1, double var1, size_t n_g,
		const double mu[], const double g[], size_t n_nonzero,
		const int nonzero_idx[], double cutoff, bool &converged,
		double buf_spa[])
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
			double root1, root2;  // the roots of equation
			int ni1, ni2;  // the number of iterations
			bool conv1, conv2;  // whether the algorithm converges or not
			getroot_K1_fast(g_pos, g_neg, root1, ni1, conv1, 0, n_nonzero,
				mu, g, q, NAmu, NAsigma);
			getroot_K1_fast(g_pos, g_neg, root2, ni2, conv2, 0, n_nonzero,
				mu, g, qinv, NAmu, NAsigma);
			if (conv1 && conv2)
			{
				double p1 = get_saddle_prob_fast(root1, n_nonzero, mu, g, q,
					NAmu, NAsigma);
				double p2 = get_saddle_prob_fast(root2, n_nonzero, mu, g, qinv,
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
