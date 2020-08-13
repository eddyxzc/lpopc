#ifndef LPLIUHPMESHREFINEALG_HPP
#define LPLIUHPMESHREFINEALG_HPP
#include "LpConf.h"
#include "LpMeshRefineImpletation.hpp"
#include "LpNLPWrapper.hpp"
#include "LpCalculateData.hpp"
#include <map>
namespace Lpopc
{
	//class PhMeshRefineAlg :public MeshRefineImpl
	LP_DECLARE_EXCEPTION(LIU_HP_ERROR);
	enum MeshOperation
	{
		////1.reducing the degree of polynomial approximation
		//2.merging mesh intervals
		NOT_SATISFIED,//not satisfied
		SATISFIED,
		REDUCED,
		MERGED//already merged,can't merge with next mesh
		
	};
	class LiuHpMeshRefineAlg : public MeshRefineImpl
	{
	public:
		LiuHpMeshRefineAlg(int Nmax, double ratio_R, double tol,
			shared_ptr<FunctionWrapper> Funs,
			shared_ptr<LpCalculateData> Data,
			shared_ptr<LpReporter>reporter)
			:Nmax_(Nmax), ratio_R_(ratio_R), mesh_tol_(tol),
			Funs_(Funs), Data_(Data), reporter_(reporter)
		{
			mesh_index_ = 0;
		}
		virtual ~LiuHpMeshRefineAlg() {};

		virtual bool RefineMesh(shared_ptr<OptimalProblem>optpro, std::vector<shared_ptr<LpMesh>>& tempMeshVector);
		static void GetLagrangeInterpCoefficientsImpl(int N,mat& alj);
	private:
		LiuHpMeshRefineAlg(const LiuHpMeshRefineAlg&);
		void operator= (const LiuHpMeshRefineAlg &);

		 void GetLagrangeInterpPowerCoefficients(int N,mat& alj);
		
		void ModifySegment(shared_ptr<OptimalProblem> optpro, size_t iphase, size_t segindex, mat seg_error,uvec state_index, shared_ptr<LpMesh> newmesh)const;
		uword Dividing_mesh(size_t iphase, size_t segindex, uvec state_index, double seg_error, shared_ptr<OptimalProblem> optpro);
		uword Increasing_N(size_t iphase, size_t segindex, uvec state_index, double seg_error, shared_ptr<OptimalProblem> optpro);
		void Reducing_N(size_t iphase, size_t segindex, uvec state_index, rowvec betai, uword& return_N, shared_ptr<OptimalProblem> optpro);
		bool Merging_mesh(size_t iphase, size_t segindex,
			std::vector<shared_ptr<LpMesh>>& tempmesh, std::vector<mat> tempstate, uword N_k, uword N_km1);
		/**
		return true when the function is smooth in this interval
		*/
		bool CanWeIncreaseN(size_t iphase, size_t segindex, uvec state_index);
		void calculate2nd_derive(vec t, mat x, vec & interp_tau, mat &derive2nd);
		/// input a vec X,return power series coefficients
		static void CalculateDi(vec x, vec& D_i);
		double mesh_tol_;
		double ratio_R_;
		size_t Nmax_;
		size_t mesh_index_;
		static std::map<int, mat> alj_map_;//coefficients of power series map to N(N is the number of LGR points)
		shared_ptr<FunctionWrapper> Funs_;
		shared_ptr<LpCalculateData> Data_;
		shared_ptr<LpReporter> reporter_;

		std::vector<std::vector<shared_ptr<LpMesh>>> mesh_history_;
		std::vector<std::vector<mat>> state_history_;
		std::vector<std::vector<vec>> mesh_points_history_;
		std::vector < std::vector<std::vector<MeshOperation>>> tag_history_;

	};

	
}// namespace Lpopc
#endif // !LPLIUHPMESHREFINEALG_HPP

