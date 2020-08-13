// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   18:09
// Email:eddy_lpopc@163.com
#include "LpGuessChecker.h"
#include <strstream>
#include <vector>
namespace Lpopc
{
	void LpGuessChecker::GetGuess(shared_ptr<OptimalProblem> opt_problem, shared_ptr<RPMGenerator>& rpm, shared_ptr<LpCalculateData>& calcul_data)
	{
		LP_DBG_START_FUN("LpGuessChecker::GetGuess()")
		static int igrid = 1;
		vector<shared_ptr<ps>> Ps(calcul_data->numphases_);
		int nlpGuessNum = 0;
		std::vector<vec> opt_nlpguess(calcul_data->numphases_);
		for (int iphase = 0; iphase<calcul_data->numphases_; iphase++)
		{
			shared_ptr<Phase> currentphase(opt_problem->GetPhase(iphase));
			int nstates = calcul_data->SIZES_[iphase][0] ;
			int ncontrol = calcul_data->SIZES_[iphase][1];
			int nparameters = calcul_data->SIZES_[iphase][2];
			int npaths = calcul_data->SIZES_[iphase][3];
			int nevents = calcul_data->SIZES_[iphase][4];
			std::vector<double> tGuess = currentphase->GetTimeGuess();
			std::vector<std::vector<double>> xGuess = currentphase->GetStateGuess();
			std::vector<std::vector<double>> uGuess = currentphase->GetControlGuess();
			std::vector<double> pGuess = currentphase->GetparameterGuess();

			if (ncontrol != 0){
				uGuess = currentphase->GetControlGuess();
			}
			//std::vector<double> pGuess;
			if (nparameters != 0){
				pGuess = currentphase->GetparameterGuess();
			}
			//Check guess for proper format
			if (igrid == 1){
				if (tGuess.size()<2){
					std::string errmsg1, errmsg2;
					errmsg1 = "Guess  must have a least two points in Phase:";
					std::strstream ss;
					ss << iphase + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_GUESSCHECKER_ERROR, errmsg1 + errmsg2);
				}
				else {
					double tem = 0; bool unique = true;
					tem = tGuess[0];
					for (int i = 1; i<tGuess.size(); i++)
					{
						if (tGuess[i] == tem) unique = false;
					}
					if (!unique){
						std::string errmsg1, errmsg2;
						errmsg1 = "Guess for time  does not contain unique valuesin phase ";
						std::strstream ss;
						ss << iphase + 1; ss >> errmsg2;
						LP_THROW_EXCEPTION(LP_GUESSCHECKER_ERROR, errmsg1 + errmsg2);
					}
				}//tguess
				if (xGuess.size() != nstates){
					std::string errmsg1, errmsg2;
					errmsg1 = "Number of states in guess does not match limits in phase ";
					std::strstream ss;
					ss << iphase + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_GUESSCHECKER_ERROR, errmsg1 + errmsg2);
				}
				else{
					for (int i = 1; i<xGuess.size(); i++)
					{
						if (xGuess[i].size() != tGuess.size()){
							std::string errmsg1, errmsg2;
							errmsg1 = "Size of state does not match size of time in guess of phase ";
							std::strstream ss;
							ss << iphase + 1; ss >> errmsg2;
							LP_THROW_EXCEPTION(LP_GUESSCHECKER_ERROR, errmsg1 + errmsg2);
						}
					}
				}//xguess
				if (uGuess.size() != ncontrol){
					std::string errmsg1, errmsg2;
					errmsg1 = "Number of controls in guess does not match limits in phase ";
					std::strstream ss;
					ss << iphase + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_GUESSCHECKER_ERROR, errmsg1 + errmsg2);
				}
				else{
					for (int i = 1; i<uGuess.size(); i++)
					{
						if (uGuess[i].size() != tGuess.size()){
							std::string errmsg1, errmsg2;
							errmsg1 = "Size of controls does not match size of time in guess of phase ";
							std::strstream ss;
							ss << iphase + 1; ss >> errmsg2;
							LP_THROW_EXCEPTION(LP_GUESSCHECKER_ERROR, errmsg1 + errmsg2);
						}
					}
				}//uguess

				if (pGuess.size() != nparameters){
					std::string errmsg1, errmsg2;
					errmsg1 = "Number of parameters in guess does not match limits in phase ";
					std::strstream ss;
					ss << iphase + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_GUESSCHECKER_ERROR, errmsg1 + errmsg2);
				}//pguess
			}//if int igrid;

			std::vector<lp_index> nodesPerInterval = currentphase->GetNodesPerInterval();
			std::vector<double> meshPoints = currentphase->GetMeshPoints();
			size_t numseg = nodesPerInterval.size();
			double t0Guess = tGuess[0];
			double tfGuess = tGuess[tGuess.size() - 1];

			rpm->initialize(numseg, meshPoints, nodesPerInterval);
			Ps[iphase].reset(new ps());
			Ps[iphase]->Points = rpm->RPMpoints();
			Ps[iphase]->Weights = rpm->RPMweights();
			Ps[iphase]->D = rpm->DifferentiationMatrix();
			Ps[iphase]->Diag = rpm->DifferentiationMatrixDiag();
			Ps[iphase]->Doffdiag = rpm->DifferentiationMatrixOffDiag();
			//vector<double>  tau_plus_ends(Ps[iphase]->Points.GetLength()+1);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//	double* rpmpoints=Ps[iphase]->Points->Values();
			/*	for(int itor=0;itor<tau_plus_ends.size()-1;itor++){
			tau_plus_ends[itor]=rpmpoints[itor];
			}*/
			//	tau_plus_ends[tau_plus_ends.size()-1]=1.0;
			vec tau_plus_ends(Ps[iphase]->Points.n_elem + 1, 1);
			tau_plus_ends.rows(0, tau_plus_ends.n_elem - 2) = Ps[iphase]->Points;
			tau_plus_ends(tau_plus_ends.n_elem - 1) = 1.0;
			vec tauGuess(tGuess.size(), 1);
			//	vector<double>tauGuess(tGuess.size());
			for (int itor = 0; itor<tauGuess.n_elem; itor++){
				tauGuess(itor) = 2 * (tGuess[itor] - t0Guess) / (tfGuess - t0Guess) - 1;
			}

			// 2.18
			vector<vec >xtauGuess(xGuess.size());
			vector<vec >utauGuess(uGuess.size());
			int totalxguess = 0, totaluguess = 0;
			for (int itor = 0; itor<xGuess.size(); itor++){
				vec xtauGuessi(tau_plus_ends.n_elem);
				vec xGuessVector(xGuess[itor]);
				xtauGuess[itor] = (xtauGuessi);
				spline_interpolation(xtauGuess[itor], tau_plus_ends, tauGuess, xGuessVector);
				totalxguess += xtauGuess[itor].n_elem;
			}
			for (int itor = 0; itor<uGuess.size(); itor++){
				vec utauGuessi(rpm->RPMpoints().n_elem, 1);
				utauGuess[itor] = (utauGuessi);
				vec uGuessVector(uGuess[itor]);

				spline_interpolation(utauGuess[itor], Ps[iphase]->Points, tauGuess, uGuessVector);
				totaluguess += utauGuess[itor].n_elem;
			}
			//nlpGuess{iphase,1} = [xinterp(:); uinterp(:); t0Guess; tfGuess; pGuess];
			/*		SmartPtr<DenseVectorSpace> nlpguessi_space=new DenseVectorSpace(totalxguess+totaluguess+2+nparameters);
			SmartPtr<DenseVector> nlpguessi=nlpguessi_space->MakeNewDenseVector();*/
			vec nlpguessi(totalxguess + totaluguess + 2 + nparameters);

			//double* nlpguessi_value=nlpguessi->Values();
			int rowshift = 0;
			for (int itor = 0; itor<nstates; itor++){
				nlpguessi.subvec(rowshift, rowshift + xtauGuess[itor].n_elem - 1) =
					xtauGuess[itor];
				rowshift += xtauGuess[itor].n_elem;
				/*	for(int index=0;index<xtauGuess[0].size();index++){
				nlpguessi_value[itor+index]=xtauGuess[itor][index];
				}*/
			}
			for (int itor = 0; itor<ncontrol; itor++){
				nlpguessi.subvec(rowshift, rowshift + utauGuess[itor].n_elem - 1) =
					utauGuess[itor];
				rowshift += utauGuess[itor].n_elem;
				/*			for(int index=0;index<utauGuess[0].size();index++){
				nlpguessi_value[itor+index]=utauGuess[itor][index];
				}*/
			}
			nlpguessi(rowshift) = t0Guess;
			rowshift++;
			nlpguessi(rowshift) = tfGuess;
			rowshift++;
			double* nlpguessi_value = nlpguessi.memptr();
			for (int itor = 0; itor<nparameters; itor++){

				nlpguessi_value[totalxguess + totaluguess + 2 + itor] = pGuess[itor];
			}
			opt_nlpguess[iphase]=(nlpguessi);
			//opt_nlpguess[iphase].print("stateguess");
			nlpGuessNum += nlpguessi.n_elem;
		}//for iphase
		calcul_data->nlpGuess = opt_nlpguess;
		calcul_data->PS = Ps;
		calcul_data->nlpGuessVector = zeros(nlpGuessNum, 1);
		int index_shift = 0;
		for (int i = 0; i < calcul_data->nlpGuess.size(); i++)
		{
			calcul_data->nlpGuessVector.subvec(index_shift, index_shift + calcul_data->nlpGuess[i].n_elem - 1) = calcul_data->nlpGuess[i];
			index_shift += calcul_data->nlpGuess[i].n_elem;
		}

	}

	

	void LpGuessChecker::spline_second_derivative(double *x, double *y, int n, double *d2y)
	{
		LP_DBG_START_FUN("spline_second_derivative")
		int i, j;
		double *c = d2y;
		double hi, him1, alphai, li = 0;

		double* mu = new double[n];
		double*  z = new double[n];
		mu[0] = 0.0;
		z[0] = 0.0;
		for (i = 1; i < n - 1; i++){
			him1 = x[i] - x[i - 1];
			hi = x[i + 1] - x[i];
			alphai = 3.0 / hi*(y[i + 1] - y[i]) - 3.0 / him1*(y[i] - y[i - 1]);
			li = 2 * (x[i + 1] - x[i - 1]) - him1*mu[i - 1];
			mu[i] = hi / li;
			z[i] = (alphai - him1*z[i - 1]) / li;
		}
		c[n - 1] = 0.0;
		for (j = n - 2; j >= 0; j--) {
			c[j] = z[j] - mu[j] * c[j + 1];
		}
		for (j = 1; j < n - 1; j++) {
			c[j] = 2 * c[j];
		}
		delete[] mu;
		delete[] z;
	}

	void LpGuessChecker::spline_interpolation(double* y, double& x, double* xdata, double* ydata, int n)
	{
		LP_DBG_START_FUN("spline_interpolation")
		int kleft, kright, k;
		double h, A, B, C, D;
		double *d2y = new double[n];
		kleft = 1;

		spline_second_derivative(xdata, ydata, n, d2y);

		kright = n;
		while (kright - kleft > 1) {
			k = (int)((kright + kleft) / 2);
			if (xdata[k - 1] > x) kright = k;
			else kleft = k;
		}
		h = xdata[kright - 1] - xdata[kleft - 1];
		if (h == 0.0){
			std::string errmsg1, errmsg2;
			errmsg1 = "Bad xdata input to routine spline_interpolation()";

			LP_THROW_EXCEPTION(SPLINE_ERROR, errmsg1);
		}
		A = (xdata[kright - 1] - x) / h;
		B = (x - xdata[kleft - 1]) / h;
		C = (pow(A, 3) - A)*(h*h) / 6.0;
		D = (pow(B, 3) - B)*(h*h) / 6.0;
		// Evaluate the cubic spline polynomial
		*y = A*ydata[kleft - 1] + B*ydata[kright - 1] + C*d2y[kleft - 1] + D*d2y[kright - 1];
		delete[]d2y;
	}

	void LpGuessChecker::spline_interpolation(vec& y, vec& x, vec& xdata, vec& ydata)
	{
		LP_DBG_START_FUN("spline_interpolation")
		LP_ASSERT_EXCEPTION(ydata.n_elem == xdata.n_elem, SPLINE_ERROR, "ydata and xdata element numbers are different");

		y = zeros(x.n_rows);
		double* xPtr = x.memptr();
		double* yPtr = y.memptr();
		double* xdataPtr = xdata.memptr();
		double* ydataPtr = ydata.memptr();
		for (int i = 0; i < x.n_elem; i++)
		{
			spline_interpolation(yPtr + i, xPtr[i], xdataPtr, ydataPtr, xdata.n_elem);
		}
	}

 void LpGuessChecker::spline_interpolation(double* y, double x, vec& xdata, vec& ydata)
{
	LP_DBG_START_FUN("spline_interpolation")
	LP_ASSERT_EXCEPTION(ydata.n_elem == xdata.n_elem, SPLINE_ERROR, "ydata and xdata element numbers are different");
	double* xdataPtr = xdata.memptr();
	double* ydataPtr = ydata.memptr();
	spline_interpolation(y, x, xdataPtr, ydataPtr, xdata.n_elem);

}

}//end of namespace Lpopc