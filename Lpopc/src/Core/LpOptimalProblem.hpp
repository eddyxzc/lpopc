// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:29
// Email:eddy_lpopc@163.com
#ifndef  LPOPTIMALPROBLEM_HPP_
#define  LPOPTIMALPROBLEM_HPP_

#include "LpConf.h"
#include "LpFunctionWrapper.h"
#include "LpException.hpp"
#include<vector>
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(LP_OTIMALPROBLEM_ERRROR);

	struct 	Limit{
		Limit(double const state0,double const state,double const statef){
			this->state[0]=state0;
			this->state[1]=state;
			this->state[2]=statef;
		};
		void GetLimit(double& state0,double& state,double& statef){
			state0=this->state[0];
			state=this->state[1];
			statef=this->state[2];
		};
		double state[3];
	};
	class Phase{
	public:
		Phase();
		Phase(int phase_index,int statenum,int controlnum,int parameternum,
			int pathnum,int eventnum):
		phase_index_(phase_index),statenum_(statenum)
			,controlnum_(controlnum),parameternum_(parameternum),pathnum_(pathnum),eventnum_(eventnum){
		hasduration_=false;};

		 ~Phase(){};

		void get_optimal_info(int& state_number,int& control_number,int& parameter_number,int& path_number,int& event_number)
		{
			state_number=statenum_;
			control_number=controlnum_;
			parameter_number=parameternum_;
			path_number=pathnum_;
			event_number=eventnum_;
		};

		void SetTimeMin(double t0,double tf){
			vtimemin.reset(new Limit(t0,0,tf));

		};
		void GetTimeMin(double &t0,double &tf){
			t0=vtimemin->state[0];
			tf=vtimemin->state[2];
		};
		void SetTimeMax(double t0,double tf){
			vtimemax.reset(new Limit(t0,0,tf));

		};
		void GetTimeMax(double &t0,double &tf){
			t0=vtimemax->state[0];
			tf=vtimemax->state[2];
		};
		void SetStateMin(double state0,double state,double statef){
			shared_ptr<Limit> statemin (new Limit(state0,state,statef));
			vstatemin.push_back(statemin);
		};
		std::vector<shared_ptr<Limit>>& GetstateMin(){
			return vstatemin;
		};
		void SetStateMax(double state0,double state,double statef){
			shared_ptr<Limit> statemax (new Limit(state0,state,statef));
			vstatemax.push_back(statemax);
		};
		std::vector<shared_ptr<Limit>>& GetstateMax(){
			return vstatemax;
		};
		void SetcontrolMin(double controlmin){
			vcontrolmin.push_back(controlmin);
		};
		std::vector<double>& GetcontrolMin(){
			return vcontrolmin;
		};

		void SetcontrolMax(double controlmax){
			vcontrolmax.push_back(controlmax);
		};
		std::vector<double>& GetcontrolMax(){
			return vcontrolmax;
		};

		void SetparameterlMin(double paramin){
			vparametermin.push_back(paramin);
		};
		std::vector<double>& GetparameterMin(){
			return vparametermin;
		};
		void SetparameterMax(double paramax){
			vparametermax.push_back(paramax);
		};
		std::vector<double>& GetparameterMax(){
			return vparametermax;
		};
		void SetpathMin(double pathmin){
			vpathmin.push_back(pathmin);
		};
		std::vector<double>& GetpathMin(){
			return vpathmin;
		};
		void SetpathMax(double pathmax){
			vpathmax.push_back(pathmax);
		};
		std::vector<double>& GetpathMax(){
			return vpathmax;
		};
		void SeteventMin(double eventmin){
			veventmin.push_back(eventmin);
		};
		std::vector<double>& GeteventMin(){
			return veventmin;
		};
		void SeteventMax(double eventmax){
			veventmax.push_back(eventmax);
		};
		std::vector<double>& GeteventMax(){
			return veventmax;
		};
		void SetDuration(double durationmin,double durationmax){
			vduration.reset( new Limit(durationmin,0,durationmax));
			hasduration_=true;
		};


		void SetTimeGuess(double guess){
			vtimeguess.push_back(guess);
		};
		std::vector<double>& GetTimeGuess(){
			return vtimeguess;
		};

		void SetStateGuess(int stateindex, double guess){
			if(vstateguess.size()>=stateindex){
				vstateguess[stateindex-1].push_back(guess);
			}else if(stateindex==(vstateguess.size()+1)){
				std::vector<double> stateguess;
				stateguess.push_back(guess);
				vstateguess.push_back(stateguess);
			}
		};
		std::vector<std::vector<double>>& GetStateGuess(){
			return vstateguess;
		};
		void SetControlGuess(int controlindex,double guess){
			if(vcontrolguess.size()>=controlindex){
				vcontrolguess[controlindex-1].push_back(guess);
			}else if(controlindex==(vcontrolguess.size()+1)){
				std::vector<double> controlguess;
				controlguess.push_back(guess);
				vcontrolguess.push_back(controlguess);
			}
		};
		std::vector<std::vector<double>>& GetControlGuess(){
			return vcontrolguess;
		};
		void SetparameterGuess(double guess){
			vparameterguess.push_back(guess);
		};
		std::vector<double>& GetparameterGuess(){
			return vparameterguess;
		};

		void SetMeshPoints(double meshpoint){
			meshpoints.push_back(meshpoint);
		};
		std::vector<double>& GetMeshPoints()
		{
			return meshpoints;
		};
		void SetNodesPerInterval(int nodes)
		{
			nodesperinterval.push_back(nodes);
		};
		std::vector<lp_index>& GetNodesPerInterval(){
			return nodesperinterval;
		};
		void SetTotalNodes(int tem )
		{
			nodes=tem;
		};
		int GetTotalNodes(){
			return nodes;
		};
		void Getduration(double &min,double &max){
			min=vduration->state[0];
			max=vduration->state[2];
		};
		bool HasDuration(){
			if (hasduration_){
				return true;
			}else{
				return false;
			}
		};
	private:
		Phase(const Phase&);
		Phase& operator=(const Phase&);
		const int phase_index_;
		const int statenum_;
		const int controlnum_;
		const int parameternum_;
		const int pathnum_;
		const int eventnum_;

		std::vector<shared_ptr<Limit>> vstatemin;
		std::vector<double> vcontrolmin;
		std::vector<double> vparametermin;
		std::vector<double> vpathmin;
		std::vector<double> veventmin;
		shared_ptr<Limit> vtimemin;
		shared_ptr<Limit> vduration;
		bool hasduration_;

		bool hasdurationmin;
		std::vector<shared_ptr<Limit>> vstatemax;
		std::vector<double> vcontrolmax;
		std::vector<double> vparametermax;
		std::vector<double> vpathmax;
		std::vector<double> veventmax;
		shared_ptr<Limit> vtimemax;
		std::vector<double> vtimeguess;
		std::vector<std::vector<double>> vstateguess;
		std::vector<std::vector<double>> vcontrolguess;
		std::vector<double> vparameterguess;

		std::vector <double> meshpoints;
		std::vector <lp_index> nodesperinterval;
		 int nodes;
	};

	class Linkage
	{
	public:
		Linkage(int ipair,int left,int right):pairindex(ipair),leftphase(left),rightphase(right)
		{};
		Linkage();
		void SetLinkMin(double min)
		{
			linkmin.push_back(min);
		}
		void SetLinkMax(double max)
		{
			linkmax.push_back(max);
		}
		std::vector<double>& GetLinkageMin()
		{
			return linkmin;
		};
		std::vector<double>& GetLinkageMax()
		{
			return linkmax;
		};
		int LeftPhase(){
			return leftphase-1;
		};
		int RightPhase(){
			return rightphase-1;
		};
		 ~Linkage(){};
	private:
		Linkage(const Linkage&);
		Linkage& operator=(const Linkage&);
		const int pairindex;
		int leftphase;
		int rightphase;
		std::vector<double> linkmin;
		std::vector<double> linkmax;
	};

	class OptimalProblem 
	{
	public:
		OptimalProblem(int numphase,int numlinkage,shared_ptr<FunctionWrapper> userfun):numphase_(numphase),
			numlink_(numlinkage),userfunction_(userfun)
		{};

		virtual~ OptimalProblem(){
		};
		void AddPhase(shared_ptr<Phase>& addingPhase){
			Phases_.push_back(addingPhase);
		};
		void AddLinkage(shared_ptr<Linkage>& addingLinkage){
			Linkage_.push_back(addingLinkage);
		};

		shared_ptr<Phase>& GetPhase(const size_t phaseindex){
			LP_ASSERT_EXCEPTION((phaseindex<numphase_)&&(phaseindex>=0),LP_OTIMALPROBLEM_ERRROR,
				"The phase index is out of rang in Function 'GetPhase' ");
			return Phases_[phaseindex];

		};

		shared_ptr<Linkage>& GetLinkage(const size_t linkindex){
			LP_ASSERT_EXCEPTION((linkindex < numphase_) && (linkindex >= 0), LP_OTIMALPROBLEM_ERRROR,
				"The linkage index is out of rang in Function 'GetLinkage' ");
			return Linkage_[linkindex];

		};
		int GetPhaseNum(){
			return numphase_;
		};
		int GetLinkageNum(){
			return numlink_;
		};
		shared_ptr<FunctionWrapper>& GetOpimalProblemFuns(){return userfunction_;};
	private:
		OptimalProblem();
		void operator=(const OptimalProblem&);
		OptimalProblem(const OptimalProblem&);
		const int numphase_;
		const int numlink_;
		shared_ptr<FunctionWrapper> userfunction_;
		std::vector<shared_ptr<Phase> >Phases_;
		std::vector<shared_ptr<Linkage>> Linkage_;
	};

}//end of Lpopc
#endif