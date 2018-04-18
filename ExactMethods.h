#ifndef EXACTMETHODS
#define EXACTMETHODS
#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>
#include "Structures.h"
#include "LocalMethods2D.h"
//#include "gurobi_c++.h"

namespace ExactMethods{

	/****************************************************************************************/
	template <class T> Assignment PA3_Permutations(T &graph){
		Assignment a;
		if(graph.getNumberOfDimensions() == 3 && graph.getDimension(0) == graph.getDimension(1) && graph.getDimension(1) == graph.getDimension(2)){
			int n = graph.getDimension(0);

			std::vector<dTuple> assignment(n, dTuple(3, 0));

			int *d2, *d3, best_cost = INT_MAX;
			d2 = new int[n], d3 = new int[n];

			for(int i = 0; i < n; i++)
				assignment[i][0] = i, d2[i] = i, d3[i] = i;

			do{
				do{
					int current_cost = 0;
					for(int i = 0; i < n; i++)
						current_cost += graph.getEdgeCost(i, d2[i], d3[i]);
					if(current_cost < best_cost){
						best_cost = current_cost;
						for(int i = 0; i < n; i++)
							assignment[i][1] = d2[i], assignment[i][2] = d3[i];
					}
				}while(std::next_permutation(d3, d3+n));
			}while(std::next_permutation(d2, d2+n));
			a.matching = assignment, a.cost = best_cost;
			delete[] d2, delete[] d3;
		}
		else{
			std::cout << "Error: This instance can not be solved by this method (PA3_DynamicProgramming)\n";
			exit(0);
		}
		return a;
	}

	/*******************************************************************************************/
	template <class T> int PA3_bitmask(int mask1, int mask2, const int n, int**& dp, T &graph){

		if(dp[mask1][mask2] < 0){
			if(mask1 == 0 || mask2 == 0)
				return dp[mask1][mask2] = 0;

			int cont = 0, mincost = INT_MAX;
			for(int i = 0; i < n; i++)
				cont += (mask1&(1<<i)) == 0 ? 0 : 1;

			for(int i = 0; i < n; i++)
				if(mask1&(1<<i))
					for(int j = 0; j < n; j++)
						if(mask2&(1<<j))
							mincost = std::min(PA3_bitmask(mask1^(1<<i), mask2^(1<<j), n, dp, graph)+graph.getEdgeCost(cont-1, i, j), mincost);
			return dp[mask1][mask2] = mincost;
		}
		return dp[mask1][mask2];
	}

	template <class T> Assignment PA3_getSolution(int mask1, int mask2, const int n, int**& dp, T &graph){
//		std::cout << "PAM_getSolution\n";
		bool salir;
		Assignment assignment;
		assignment.cost = dp[mask1][mask2];
		dTuple tuple(3, 0);
		for(int cont = n-1; cont >= 0; cont--){
//			printf("%x %x %x\n", cont, mask1, mask2);
			salir = false;
			for(int i = 0; i < n && !salir; i++){
				if(mask1&(1<<i))
					for(int j = 0; j < n; j++)
						if(mask2&(1<<j)){
							if(dp[mask1][mask2] == (dp[mask1^(1<<i)][mask2^(1<<j)]+graph.getEdgeCost(cont, i, j))){
								salir = true;
								mask1 = mask1^(1<<i), mask2 = mask2^(1<<j);
								tuple[0] = cont, tuple[1] = i, tuple[2] = j;
								assignment.matching.push_back(tuple);
								break;
							}
						}
			}
		}
		std::sort(assignment.matching.begin(), assignment.matching.end());
		return assignment;
	}

	template <class T> Assignment PA3_DynamicProgramming(T &graph){
//		std::cout << "Running PA3_DynamicProgramming\n";
		Assignment a;
		if(graph.getNumberOfDimensions() == 3 && graph.getDimension(0) == graph.getDimension(1) && graph.getDimension(1) == graph.getDimension(2)){
			int **dp, lim = 1<<graph.getDimension(0);
			dp = new int*[lim];
			for(int i = 0; i < lim; i++){
				dp[i] = new int[lim];
				for(int j = 0; j < lim; j++)
					dp[i][j] = -1;
			}
			ExactMethods::PA3_bitmask(lim-1, lim-1, graph.getDimension(0), dp, graph);
			a = PA3_getSolution(lim-1, lim-1, graph.getDimension(0), dp, graph);
//			std::cout << dp[lim-1][lim-1]<<"\n";
			for(int i = 0; i < lim; i++)
				delete[] dp[i];
			delete[] dp;
		}
		else{
			std::cout << "Error: This instance can not be solved by this method (PA3_DynamicProgramming)\n";
			exit(0);
		}
		return a;
	}

	/**********************************************************************************************/

	template <class T> Assignment PA3_PermutationsHungarian(T &graph){
		std::cout << "Running PA3_PermutationsHungarian\n";
		Assignment answer(INT_MAX), aux;
		if(graph.getNumberOfDimensions() == 3 && graph.getDimension(0) == graph.getDimension(1) && graph.getDimension(1) == graph.getDimension(2)){
			int n = graph.getDimension(0);
			aux.matching = std::vector<dTuple>(n, dTuple(3, 0));
			int *d2 = new int[n];

			for(int i = 0; i < n; i++)
				aux.matching[i][0] = aux.matching[i][2] = i, d2[i] = i;

			do{
				aux.cost = 0;
				for(int i = 0; i < n; i++){
					aux.matching[i][1] = d2[i];
					aux.cost += graph.getEdgeCost(i, d2[i], i);
				}

				LocalMethods2D::RoundLocal2Dim(1, 2, aux, graph);

				if(aux.cost	< answer.cost)
					answer = aux;
			}while(std::next_permutation(d2, d2+n));
			delete[] d2;
		}
		else{
			std::cout << "Error: This instance can not be solved by this method (PA3_PermutationsHungarian)\n";
			exit(0);
		}
		return answer;
	}


	/**************************************************************************************/

	template <class T> Assignment PA3_PermutationsHungarianFast(T &graph){
		std::cout << "Running PA3_PermutationsHungarianFast\n";
		Assignment answer(INT_MAX);
		if(graph.getNumberOfDimensions() == 3 && graph.getDimension(0) == graph.getDimension(1) && graph.getDimension(1) == graph.getDimension(2)){
			int n = graph.getDimension(0);
			std::vector<dTuple> assignment(n, dTuple(3, 0));
			int *d2 = new int[n];

			for(int i = 0; i < n; i++)
				assignment[i][0] = i, d2[i] = i;

			Vertex *G = new Vertex[n];
			for(int x = 0; x < n; x++)
				G[x].Allocate(n);

			do{
				for(int i = 0; i < n; i++)
					G[i].length = 0;
				for(int i = 0; i < n; i++)
					for(int j = 0; j < n; j++)
						G[d2[i]].add(j, graph.getEdgeCost(i, d2[i], j));

				std::vector<std::pair<std::pair<int, int>, int> > matching = BasicMethods::SolveAP(G, n, n, 3);
				int currentCost = 0;
				for(int i = 0; i < matching.size(); i++)
					currentCost += matching[i].second;

				if(currentCost	< answer.cost){
					answer.cost = currentCost;
					std::sort(matching.begin(), matching.end());
					for(int i = 0; i < n; i++)
						assignment[i][1] = d2[i], assignment[i][2] = matching[d2[i]].first.second;
					answer.matching = assignment;
				}
			}while(std::next_permutation(d2, d2+n));
			delete[] d2;
			delete[] G;
		}
		else{
			std::cout << "Error: This instance can not be solved by this method (PA3_PermutationsHungarianFast)\n";
			exit(0);
		}
		return answer;
	}

	/*****************************************************************************************************/

/*		void AllocateMemoryGRBVar(GRBVar ****p, const int n, GRBModel &model){
		*p = new GRBVar**[n];
		for(int i = 0; i < n; i++){
			(*p)[i] = new GRBVar*[n];
			for(int j = 0; j < n; j++){
				(*p)[i][j] = new GRBVar[n];
				for(int k = 0; k < n; k++)
					(*p)[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}
		// Integrate new variables
		model.update();
		return;
	}

	void DeallocateMemoryGRBVar(GRBVar ****p, const int n){
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++)
				delete[] (*p)[i][j];
			delete[] (*p)[i];
		}
		delete[] *p;
		return;
	}

	template <class T> GRBLinExpr GenerateObjetiveFunction(const int n, GRBVar ***p, T &graph){
		GRBLinExpr expr;
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				for(int k = 0; k < n; k++)
					expr += graph.getEdgeCost(i, j, k)*p[i][j][k];
		return expr;
	}

	std::vector<dTuple> getGurobiSolution(const int n, GRBVar ***p, GRBModel &model){
		std::vector<dTuple> assignment(n, dTuple(3, 0));
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				for(int k = 0; k < n; k++)
					if(p[i][j][k].get(GRB_DoubleAttr_X) >= 1.0)
						assignment[i][0] = i, assignment[i][1] = j, assignment[i][2] = k;
//		std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << "\n";
		return assignment;
	}
	//*/

	template <class T> Assignment PA3_GurobiSolver(T &graph){
//		std::cout << "Running PA3_Gurobi_solver\n";
		Assignment assignment;
/*		if(graph.getNumberOfDimensions() == 3 && graph.getDimension(0) == graph.getDimension(1) && graph.getDimension(1) == graph.getDimension(2)){
			int n = graph.getDimension(0);
			try {
				static GRBEnv env = GRBEnv();
				env.set(GRB_IntParam_LogToConsole, 0);
				GRBModel model = GRBModel(env);

				// Create variables
				GRBVar ***p;
				AllocateMemoryGRBVar(&p, n, model);
				GRBLinExpr objetiveFuntion = GenerateObjetiveFunction(n, p, graph);

				// Set objective
				model.setObjective(objetiveFuntion, GRB_MINIMIZE);

				// Add constraints:
				for(int i = 0; i < n; i++){
					GRBLinExpr expr;
					for(int j = 0; j < n; j++)
						for(int k = 0; k < n; k++)
							expr += p[i][j][k];
					model.addConstr(expr == 1);
				}

				for(int j = 0; j < n; j++){
					GRBLinExpr expr;
					for(int i = 0; i < n; i++)
						for(int k = 0; k < n; k++)
							expr += p[i][j][k];
					model.addConstr(expr == 1);
				}

				for(int k = 0; k < n; k++){
					GRBLinExpr expr;
					for(int i = 0; i < n; i++)
						for(int j = 0; j < n; j++)
							expr += p[i][j][k];
					model.addConstr(expr == 1);
				}

				// Optimize model
				model.optimize();
				assignment.matching = getGurobiSolution(n, p, model);
				assignment.cost = (int)model.get(GRB_DoubleAttr_ObjVal);

				DeallocateMemoryGRBVar(&p, n);
			} catch(GRBException e) {
				std::cout << "Error code = " << e.getErrorCode() << "\n";
				std::cout << e.getMessage() << "\n";
			} catch(...) {
				std::cout << "Exception during optimization\n";
			}
		}
		else{
			std::cout << "Error: This instance can not be solved by this method (PA3_Gurobi_solver)\n";
			exit(0);
		}//*/
		return assignment;
	}


	//
/*	void AllocateMemoryGRBVar(GRBVar *****p, const int n, GRBModel &model){
		*p = new GRBVar***[n];
		for(int i = 0; i < n; i++){
			(*p)[i] = new GRBVar**[n];
			for(int j = 0; j < n; j++){
				(*p)[i][j] = new GRBVar*[n];
				for(int k = 0; k < n; k++){
					(*p)[i][j][k] = new GRBVar[n];
					for(int a = 0; a < n; a++)
						(*p)[i][j][k][a] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				}
			}
		}
		// Integrate new variables
		model.update();
		return;
	}

	void DeallocateMemoryGRBVar(GRBVar *****p, const int n){
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++)
					delete[] (*p)[i][j][k];
				delete[] (*p)[i][j];
			}
			delete[] (*p)[i];
		}
		delete[] *p;
		return;
	}

	template <class T> GRBLinExpr GenerateObjetiveFunction(const int n, GRBVar ****p, T &graph){
		GRBLinExpr expr;
		for(int i = 0; i < n; i++)
			for(int j = 0; j <  n; j++)
				for(int k = 0; k < n; k++)
					for(int a = 0; a < n; a++)
						expr += graph.getEdgeCost(i, j, k, a)*p[i][j][k][a];
		return expr;
	}

	std::vector<dTuple> getGurobiSolution(const int n, GRBVar ****p, GRBModel &model){
		std::vector<dTuple> assignment(n, dTuple(4, 0));
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				for(int k = 0; k < n; k++)
					for(int a = 0; a < n; a++)
						if(p[i][j][k][a].get(GRB_DoubleAttr_X) >= 1.0)
							assignment[i][0] = i, assignment[i][1] = j, assignment[i][2] = k, assignment[i][3] = a;
		std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << "\n";
		return assignment;
	}//*/

	template <class T> Assignment PA4_GurobiSolver(T &graph){
		std::cout << "Running PA4_Gurobi_solver\n";
		Assignment assignment;
/*		if(graph.getNumberOfDimensions() == 4 && graph.getDimension(0) == graph.getDimension(1) && graph.getDimension(1) == graph.getDimension(2)
			&& graph.getDimension(2) == graph.getDimension(3)){
			int n = graph.getDimension(0);
			try {
				GRBEnv env = GRBEnv();
				env.set(GRB_IntParam_LogToConsole, 0);
				GRBModel model = GRBModel(env);

				// Create variables
				GRBVar ****p;
				AllocateMemoryGRBVar(&p, n, model);

				GRBLinExpr objetiveFuntion = GenerateObjetiveFunction(n, p, graph);

				model.setObjective(objetiveFuntion, GRB_MINIMIZE);

				// Add constraints:
				for(int i = 0; i < n; i++){
					GRBLinExpr expr;
					for(int j = 0; j < n; j++)
						for(int k = 0; k < n; k++)
							for(int a = 0; a < n; a++)
								expr += p[i][j][k][a];
					model.addConstr(expr == 1);
				}

				for(int j = 0; j < n; j++){
					GRBLinExpr expr;
					for(int i = 0; i < n; i++)
						for(int k = 0; k < n; k++)
							for(int a = 0; a < n; a++)
								expr += p[i][j][k][a];
					model.addConstr(expr == 1);
				}

				for(int k = 0; k < n; k++){
					GRBLinExpr expr;
					for(int i = 0; i < n; i++)
						for(int j = 0; j < n; j++)
							for(int a = 0; a < n; a++)
								expr += p[i][j][k][a];
					model.addConstr(expr == 1);
				}

				for(int a = 0; a < n; a++){
					GRBLinExpr expr;
					for(int i = 0; i < n; i++)
						for(int j = 0; j < n; j++)
							for(int k = 0; k < n; k++)
								expr += p[i][j][k][a];
					model.addConstr(expr == 1);
				}

				// Optimize model
				model.optimize();

				assignment.matching = getGurobiSolution(n, p, model);
				assignment.cost = (int)model.get(GRB_DoubleAttr_ObjVal);

				DeallocateMemoryGRBVar(&p, n);
			} catch(GRBException e) {
				std::cout << "Error code = " << e.getErrorCode() << "\n";
				std::cout << e.getMessage() << "\n";
			} catch(...) {
				std::cout << "Exception during optimization\n";
			}
		}//*/
		return assignment;
	}


/*	void AllocateMemoryGRBVar(GRBVar **p, const std::vector<int> &dims, GRBModel &model){
		int size = 1;
		for(int i = 0; i < dims.size(); i++)
			size *= dims[i];
		*p = new GRBVar[size];
		for(int i = 0; i < size; i++)
			(*p)[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		// Integrate new variables
		model.update();
		return;
	}

	void DeallocateMemoryGRBVar(GRBVar **p){
		delete[] *p;
		return;
	}

	template <class T> void GenerateObjetiveFunction(const int s, const std::vector<int> &dims, GRBVar *p, const int index,
		const int prod, GRBLinExpr &expr, T &graph){
		if(s == 0){
			expr += graph.getEdgeCost(index)*p[index];
			return;
		}
		for(int i = 0; i < dims[s-1]; i++)
			GenerateObjetiveFunction(s-1, dims, p, i*prod+index, prod*dims[s-1], expr, graph);
		return;
	}

	dTuple getTuple(int index, std::vector<int> dims){
		dTuple e(dims.size(), 0);
		for(int i = (int)dims.size()-1; i >= 0; i--){
			e[i] = index%dims[i];
			index /= dims[i];
		}
		return e;
	}

	std::vector<dTuple> getGurobiSolution(std::vector<int> dims, GRBVar *p, GRBModel &model){
		std::vector<dTuple> assignment;
		int size = 1;
		for(int i = 0; i < dims.size(); i++)
			size *= dims[i];
		for(int index = 0; index < size; index++)
			if(p[index].get(GRB_DoubleAttr_X) >= 1.0)
				assignment.push_back(getTuple(index, dims));
//		std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << "\n";
		return assignment;
	}

	void addRestriction(const int s, const int fixed, const int restriction, const int index, const int prod,
		const std::vector<int> &dims, GRBLinExpr &expr, GRBVar *p){
		if(s == 0)
			expr += p[index];
		else{
			if(fixed != s){
				for(int i = 0; i < dims[s-1]; i++)
					addRestriction(s-1, fixed, restriction, i*prod+index, prod*dims[s-1], dims, expr, p);
			}
			else
				addRestriction(s-1, fixed, restriction, restriction*prod+index, prod*dims[s-1], dims, expr, p);
		}
		return;
	}//*/

	template <class T> Assignment PAM_GurobiSolver(T &graph){
//		std::cout << "Running PAM_Gurobi_solver************\n";
		Assignment assignment;
/*		try {
			static GRBEnv env = GRBEnv();
//			env.set(GRB_IntParam_OutputFlag, 1);
			env.set(GRB_IntParam_LogToConsole, 0);
			GRBModel model = GRBModel(env);

			// Create variables
			GRBVar *p;
			AllocateMemoryGRBVar(&p, graph.getDimensions(), model);

			GRBLinExpr objetiveFuntion;
			GenerateObjetiveFunction(graph.getNumberOfDimensions(), graph.getDimensions(), p, 0, 1, objetiveFuntion, graph);
			model.setObjective(objetiveFuntion, GRB_MINIMIZE);

			// Add constraints:
			for(int s = 0; s < graph.getNumberOfDimensions(); s++)
				for(int k = 0; k < graph.getDimension(s); k++){
					GRBLinExpr expr;
					addRestriction(graph.getNumberOfDimensions(), s+1, k, 0, 1, graph.getDimensions(), expr, p);
					model.addConstr(expr == 1);
				}
			// Optimize model
			model.optimize();

			assignment.matching = getGurobiSolution(graph.getDimensions(), p, model);
			assignment.cost = (int)model.get(GRB_DoubleAttr_ObjVal);
			DeallocateMemoryGRBVar(&p);
		} catch(GRBException e) {
			std::cout << "Error code = " << e.getErrorCode() << "\n";
			std::cout << e.getMessage() << "\n";
		} catch(...) {
			std::cout << "Exception during optimization\n";
		}//*/
		return assignment;
	}
//*/
	template <class T> std::pair<std::pair<double, double>, std::pair<int, int> > ExactMethods(const int exactMethodType,
		Assignment &assignment, T &graph){
		clock_t comienzo = clock();
		std::cout << "Start Exact Method\n";
		switch(exactMethodType){
			case 71:{
				assignment = PA3_Permutations(graph);
			}break;
			case 72:{
				assignment = PA3_DynamicProgramming(graph);
			}break;
			case 73:{
				assignment = PA3_PermutationsHungarianFast(graph);
			}break;//*/
			case 74:{
				assignment = PA3_GurobiSolver(graph);
			}break;
			case 75:{
				assignment = PA4_GurobiSolver(graph);
			}break;
			case 76:{
				assignment = PAM_GurobiSolver(graph);
			}break;//*/
			default:break;
		}
		std::cout << "\nFinal cost = " << assignment.cost << "\n";
		std::sort(assignment.matching.begin(), assignment.matching.end());
		Output:: PrintSolution(assignment.matching);
		std::cout << "End ExactMethod\n";
		double final_time = ((clock()-comienzo)/(double)CLOCKS_PER_SEC);
		std::cout << "The running time was "<< final_time <<" seconds.\n";
		return std::make_pair(std::make_pair(assignment.cost, final_time), std::make_pair(1, assignment.cost));
	}
};

#endif
