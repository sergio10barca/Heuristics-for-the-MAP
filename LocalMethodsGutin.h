#ifndef LOCALMETHODSGUTIN
#define LOCALMETHODSGUTIN
#include <vector>
#include <algorithm>
#include <utility>
#include "Output.h"
#include "BasicMethods.h"

namespace LocalMethodsGutin{

	/*
		Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on:
	*/
	template <class T> void UpdateSolution(const bool *F, std::vector<std::pair<std::pair<int, int>, int> > &matching, Assignment &solution, T &graph){
		int n = (int)matching.size(), s = (int)graph.getNumberOfDimensions(), &cost = solution.cost;
		cost = 0;
		std::vector<dTuple> aux = solution.matching;
		for(int i = 0; i < n; i++){
			for(int d = 0; d < s; d++)
				if(F[d])
					solution.matching[matching[i].first.first][d] = aux[matching[i].first.second][d];
			cost += matching[i].second;
		}
		std::sort(solution.matching.begin(), solution.matching.end());
		BasicMethods::verify(solution, graph, "UpdateSolution");
		return;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on:*/
	template <class T> bool DimensionwiseVariation(Assignment &solution, T &graph, int apAlgorithm = 3){
//		std::cout << "Running DimensionwiseVariation\n";
		int n = graph.getDimension(0), previousCost = solution.cost;
		bool *F = NULL;
		F = new bool[graph.getNumberOfDimensions()];
		for(int i = 0; i < graph.getNumberOfDimensions(); i++)
			F[i] = false;
		for(int K = 0; K < graph.getNumberOfDimensions(); K++){
			Vertex *G = NULL;
			F[K] = true;
			BasicMethods::GetBipartiteGraph(F, G, graph.getNumberOfDimensions(), n, solution.matching, graph);
			std::vector<std::pair<std::pair<int, int>, int> > matching = BasicMethods::SolveAP(G, n, n, apAlgorithm);
			UpdateSolution(F, matching, solution, graph);
			F[K] = false;
			delete[] G;
		}
		delete[] F;
		return previousCost > solution.cost;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on:	*/
	template <class T> bool MultiDimensionwiseVariation(Assignment &solution, T &graph){
//		std::cout << "Running MultiDimensionwiseVariation\n";
		int n = graph.getDimension(0), lim = (1<<graph.getNumberOfDimensions())-1, previousCost = solution.cost;
		bool *F = NULL, *visited;
		F = new bool[graph.getNumberOfDimensions()];
		for(int i = 0; i < graph.getNumberOfDimensions(); i++)
			F[i] = false;
		visited = new bool[lim];
		for(int i = 0; i < lim; i++)
			visited[i] = true;
		for(int K = 1; K < lim; K++){
			int notK = (~K)&(lim-1);
			if(visited[K]){
				visited[K] = visited[notK] = false;
				Vertex *G = NULL;
				for(int d = 0, j = K; d < graph.getNumberOfDimensions(); j>>=1, d++)
					F[d] = j&1;

				BasicMethods::GetBipartiteGraph(F, G, graph.getNumberOfDimensions(), n, solution.matching, graph);
				std::vector<std::pair<std::pair<int, int>, int> > matching = BasicMethods::SolveAP(G, n, n, 3);
				UpdateSolution(F, matching, solution, graph);
				delete[] G;
			}
		}
		delete[] F, delete[] visited;
		return previousCost > solution.cost;
	}

	/*
		Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:
	*/
	template <class T> void Variable2Opt(const int s, const int *indices, int &bestCost, Assignment &solution, dTuple tuples[], T &graph){
		if(s == 1){
			int cost = graph.getEdgeCost(solution.matching[indices[0]])+graph.getEdgeCost(solution.matching[indices[1]]);
			if(cost < bestCost)
				bestCost = cost, tuples[0] = solution.matching[indices[0]], tuples[1] = solution.matching[indices[1]];
			return;
		}
		else{
			Variable2Opt(s-1, indices, bestCost, solution, tuples, graph);
			std::swap(solution.matching[indices[0]][s-1], solution.matching[indices[1]][s-1]);
			Variable2Opt(s-1, indices, bestCost, solution, tuples, graph);
			std::swap(solution.matching[indices[0]][s-1], solution.matching[indices[1]][s-1]);
		}
		return;
	}

	/*
		Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:
	*/
	template <class T> void Variable3Opt(const int s, const int *indices, int &bestCost, Assignment &solution, dTuple tuples[], T &graph){
		if(s == 1){
//			Output::PrintSolution(solution, graph);
			int cost = graph.getEdgeCost(solution.matching[indices[0]]) + graph.getEdgeCost(solution.matching[indices[1]])
				+ graph.getEdgeCost(solution.matching[indices[2]]);
			if(cost < bestCost)
				bestCost = cost, tuples[0] = solution.matching[indices[0]], tuples[1] = solution.matching[indices[1]],
				tuples[2] = solution.matching[indices[2]];
			return;
		}
		else{
			Variable3Opt(s-1, indices, bestCost, solution, tuples, graph);
			std::swap(solution.matching[indices[0]][s-1], solution.matching[indices[1]][s-1]);
			Variable3Opt(s-1, indices, bestCost, solution, tuples, graph);
			std::swap(solution.matching[indices[1]][s-1], solution.matching[indices[2]][s-1]);
			Variable3Opt(s-1, indices, bestCost, solution, tuples, graph);
			std::swap(solution.matching[indices[0]][s-1], solution.matching[indices[1]][s-1]);
			Variable3Opt(s-1, indices, bestCost, solution, tuples, graph);
			std::swap(solution.matching[indices[1]][s-1], solution.matching[indices[2]][s-1]);
			Variable3Opt(s-1, indices, bestCost, solution, tuples, graph);
			std::swap(solution.matching[indices[0]][s-1], solution.matching[indices[1]][s-1]);
			Variable3Opt(s-1, indices, bestCost, solution, tuples, graph);
			std::swap(solution.matching[indices[1]][s-1], solution.matching[indices[2]][s-1]);
		}
		return;
	}

	/*
		Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on:
	*/
	template <class T> void VariableKOpt(const int s, const int K, const int *indices, int &bestCost, Assignment &solution,
		dTuple tuples[], T &graph){
		if(s == 1){
			int cost = 0;
			for(int i = 0; i < K; i++)//This takes O(n)
				cost += graph.getEdgeCost(solution.matching[indices[i]]);
			if(cost < bestCost){
				bestCost = cost;
				for(int i = 0; i < K; i++)
					tuples[i] = solution.matching[indices[i]];
			}
			return;
		}
		else{
			int *v = new int[K];
			for(int i = 0; i < K; i++)
				v[i] = solution.matching[indices[i]][s-1];
			std::sort(v, v+K);
			do{
				for(int i = 0; i < K; i++)
					solution.matching[indices[i]][s-1] = v[i];
				VariableKOpt(s-1, K, indices, bestCost, solution, tuples, graph);
			}while(std::next_permutation(v, v+K));
			delete[] v;
		}
		return;
	}

	/*
		Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on:
	*/
	template <class T> bool TwoOpt(Assignment &solution, T &graph){
//		std::cout << "Running TwoOptGutin\n";
		int n = graph.getDimension(0), indices[2], previousCost = solution.cost;
		for(int i = 0; i < n; i++){
			indices[0] = i;
			for(int j = i+1; j < n; j++){
				indices[1] = j;
//				std::cout << "Iteration " << i << " " << j << "\n";
				int cost = graph.getEdgeCost(solution.matching[i])+graph.getEdgeCost(solution.matching[j]), bestCost;
				bestCost = cost;
				dTuple tuples[2] = {solution.matching[i], solution.matching[j]};
				Variable2Opt(graph.getNumberOfDimensions(), indices, bestCost, solution, tuples, graph);
//				VariableKOpt(graph.getNumberOfDimensions(), 2, indices, bestCost, solution, tuples, graph);
				solution.matching[i] = tuples[0], solution.matching[j] = tuples[1];
				solution.cost += (bestCost-cost);//*/
				BasicMethods::verify(solution, graph, "TwoOptGutin");
			}
		}
		return previousCost > solution.cost;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on: July, 2016	*/
	template <class T> bool ThreeOpt(Assignment &solution, T &graph){
//		std::cout << "Running ThreeOpt\n";
		int n = graph.getDimension(0), indices[3], previousCost = solution.cost;
		for(int i = 0; i < n; i++){
			indices[0] = i;
			for(int j = i+1; j < n; j++){
				indices[1] = j;
				for(int k = j+1; k < n; k++){
					indices[2] = k;
					int cost = graph.getEdgeCost(solution.matching[i])+graph.getEdgeCost(solution.matching[j])+graph.getEdgeCost(solution.matching[k]);
					int bestCost = cost;
					dTuple tuples[3] = {solution.matching[i], solution.matching[j], solution.matching[k]};
					Variable3Opt(graph.getNumberOfDimensions(), indices, bestCost, solution, tuples, graph);
//					VariableKOpt(graph.getNumberOfDimensions(), 3, indices, bestCost, solution, tuples, graph);
					solution.matching[i] = tuples[0], solution.matching[j] = tuples[1], solution.matching[k] = tuples[2];
					solution.cost += (bestCost-cost);//*/

					BasicMethods::verify(solution, graph, "ThreeOpt");
				}
			}
		}
		return previousCost > solution.cost;
	}

	template <class T> void KOpt(const int index, const int prev, const int K, const int n, int *indices, Assignment &solution, T &graph){
		if(index == K){
			dTuple *tuples = new dTuple[K];
			int cost = 0, bestCost;
			for(int i = 0; i < K; i++)
				cost += graph.getEdgeCost(tuples[i] = solution.matching[indices[i]]);
			bestCost = cost;
			VariableKOpt(graph.getNumberOfDimensions(), K, indices, bestCost, solution, tuples, graph);
			//solution = bestSolution;
			for(int i = 0; i < K; i++)
				solution.matching[ indices[i] ] = tuples[i];
			solution.cost += (bestCost - cost);
			delete[] tuples;
		}
		else{
			for(int i = prev; i < n; i++){
				indices[index] = i;
				KOpt(index+1, i+1, K, n, indices, solution, graph);
			}
		}
		return;
	}

	template <class T> bool VarOpt(const int K, Assignment &solution, T &graph){
//		std::cout << "Running "<<K<<"-Opt\n";
		int n = graph.getDimension(0), previousCost = solution.cost, *indices;
		indices = new int[K];
		KOpt(0, 0, K, n, indices, solution, graph);
		BasicMethods::verify(solution, graph, "K-OptGutin");
		delete[] indices;
		return previousCost > solution.cost;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: July, 2016
		Last Update on:	*/
	template <class T> bool CompleteDimensionwiseVariation(Assignment &solution, T &graph){
//		std::cout << "Running Complete DimensionwiseVariation\n";
		int previousCost;
		do{
			previousCost = solution.cost;
			DimensionwiseVariation(solution, graph);
		}while(previousCost > solution.cost);
		return false;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: July, 2016
		Last Update on:	*/
	template <class T> bool CompleteMultiDimensionwiseVariation(Assignment &solution, T &graph){
//		std::cout << "Running Complete MultiDimensionwiseVariation\n";
		int previousCost;
		do{
			previousCost = solution.cost;
			MultiDimensionwiseVariation(solution, graph);
		}while(previousCost > solution.cost);
		return false;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: July, 2016
		Last Update on:	*/
	template <class T> bool CompleteTwoOpt(Assignment &solution, T &graph){
//		std::cout << "Running Complete Two Opt\n";
		int previousCost;
		do{
			previousCost = solution.cost;
			TwoOpt(solution, graph);
		}while(previousCost > solution.cost);
		return false;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: July, 2016
		Last Update on:	*/
	template <class T> bool CompleteThreeOpt(Assignment &solution, T &graph){
//		std::cout << "Running Complete Three Opt\n";
		int previousCost;
		do{
			previousCost = solution.cost;
			ThreeOpt(solution, graph);
		}while(previousCost > solution.cost);
		return false;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: July, 2016
		Last Update on:	*/
	template <class T> bool CompleteKOpt(const int K, Assignment &solution, T &graph){
//		std::cout << "Running Complete "<<K<<" Opt\n";
		int previousCost;
		do{
			previousCost = solution.cost;
			VarOpt(K, solution, graph);
		}while(previousCost > solution.cost);
		return false;
	}

	template <class T> void CombinationsFInterchange(int n, int k, dTuple &answer, dTuple &u, dTuple &w, dTuple &minimum, T &graph){
		if(n == 0 || k == 0){
			if(k > 0)
				return;
			for(int i = 0; i < n; i++)
				answer[i] = u[i];
			if(minimum.size() == 0 || graph.getEdgeCost(answer) < graph.getEdgeCost(minimum))
				minimum = answer;
			return;
		}
		answer[n-1] = w[n-1];//if j \in F
		CombinationsFInterchange(n-1, k-1, answer, u, w, minimum, graph);
		if(n-1 >= k){
			answer[n-1] = u[n-1];//if j not \in F
			CombinationsFInterchange(n-1, k, answer, u, w, minimum, graph);
		}
		return;
	}

	template<class T> dTuple CalculateArgMin(dTuple u, dTuple w, T &graph){
//		std::cout << "CalculateArgMin \n";
		int s = (int)u.size(), F = s/2;
		dTuple minimum, current(s, 0);
		for(int i = 1; i <= F; i++)
			CombinationsFInterchange(s, i, current, u, w, minimum, graph);
		return minimum;
	}

	dTuple CalculateVComplement(dTuple &eC, dTuple &eM, dTuple &v){
		dTuple vBar = v;
		for(int i = 0; i < v.size(); i++)
			vBar[i] = v[i] == eC[i] ? eM[i] : eC[i];
		return vBar;
	}

	template <class T> bool VOpt(Assignment &solution, T &graph){
//		std::cout << "Running VOpt\n";
		int n = (int)solution.matching.size(), previousCost = solution.cost;
		bool *available;
		//for every c = 1 ... n. I have 0-indexed vectors
		available = new bool[n];
		//for every c 0 1, 2, ..., n
		for(int c = 0; c < n; c++){
			//Initialize the total gain G = 0, the best assignment Abest = A, and a set of available vector indices I = {1, 2, .., n}\{c}
			int m = 0, Gain = 0;
			Assignment ABest = solution;
			for(int i = 0; i < n; i++)
				available[i] = true;
			available[c] = false;
			for(int availablesCont = n-1; availablesCont > 0; availablesCont--){
				//set m = argmin
				dTuple vMin, vBar, &eC = solution.matching[c];
				for(int i = 0; i < n; i++){
					if(available[i]){
						dTuple current = CalculateArgMin(eC, solution.matching[i], graph);
						if(vMin.size() == 0 || graph.getEdgeCost(current) < graph.getEdgeCost(vMin)){
							m = i;
							vMin = current;
						}
					}
				}
				vBar = CalculateVComplement(eC, solution.matching[m], vMin);
				Gain = Gain + graph.getEdgeCost(eC) - graph.getEdgeCost(vMin);
//				std::cout << "Gain = "<< Gain << "\n";
				if(Gain <= 0)
					break;
				//Replace eM with v and eC with vBar. Mark eM as unavailable.
				solution.cost -= graph.getEdgeCost(solution.matching[m]);
				solution.cost -= graph.getEdgeCost(eC);
				solution.matching[m] = vMin;
				eC = vBar;
				solution.cost += graph.getEdgeCost(solution.matching[m]);
				solution.cost += graph.getEdgeCost(eC);
				available[m] = false;
//				std::cout <<"costs = "<< solution.cost <<" "<< ABest.cost<< "\n";
				if(solution.cost < ABest.cost)
					ABest = solution;
			}
/*			std::cout << "final costs "  << ABest.cost << " " <<solution.cost<< "\n";//*/
			solution = ABest;
		}
		delete[] available;
		BasicMethods::verify(solution, graph, "VarOptGutin");
		return previousCost > solution.cost;
	}


	template <class T> bool MDVH2_3OPT(Assignment &assignment, T &graph){
//		std::cout << "Running VOpt\n";
		int prev_cost = assignment.cost;
		do{
			CompleteMultiDimensionwiseVariation(assignment, graph);
			if(prev_cost <= assignment.cost){
				break;
			}
			prev_cost = assignment.cost;
			CompleteKOpt(3, assignment, graph);
			if(prev_cost <= assignment.cost){
				break;
			}
			prev_cost = assignment.cost;
		}while(true);
		return true;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: July, 2016
		Last Update on:	*/
	template <class T> std::pair<std::pair<double, double>, std::pair<int, int> > IteratedLocalSearchGK(
		Assignment &solution, T &graph,
		const int heuristic, const int iterations, const bool seconds, const int K){
		std::cout << "Running Multi Start for Methods provided by Gutin and Karapetyan\n";
		clock_t comienzo = clock();
		double avg_time = 0, avg_solution = 0;
		int iteration = 0;
		Assignment assignment = solution;
		for(iteration = 0; (!seconds && iteration < iterations) ||
			(seconds && (int)BasicMethods::GetRunningTime(comienzo) < iterations); iteration++){
			switch(heuristic){
				case 86:{
					CompleteDimensionwiseVariation(assignment, graph);
				}break;
				case 87:{
					CompleteMultiDimensionwiseVariation(assignment, graph);
				}break;
				case 88:{
					CompleteTwoOpt(assignment, graph);
				}break;
				case 89:{
					CompleteThreeOpt(assignment, graph);
				}break;
				case 90:{
					CompleteKOpt(K, assignment, graph);
				}break;
				default:break;
			}
//			std::cout <<"Iter "<<iteration << " "<<assignment.cost << " " << solution.cost <<"\n";
			avg_solution += assignment.cost;
			if(solution.cost > assignment.cost){
				solution = assignment;
				std::cout <<"Iter = "<<iteration << " before "<<assignment.cost 
					<< " new " << solution.cost<<" time "<<(int)BasicMethods::GetRunningTime(comienzo) <<"\n";
			}

			assignment = solution;
			for(int dim = 1; dim < graph.getNumberOfDimensions(); dim++)
				BasicMethods::Perturbation(assignment, 2, dim);
			assignment.cost = BasicMethods::CalculateCost(assignment.matching, graph);
		}
		avg_time = BasicMethods::GetRunningTime(comienzo) / (double) iteration;
		avg_solution /= (double) iteration;
		return std::make_pair(std::make_pair(avg_solution, avg_time), std::make_pair(iteration, solution.cost));
	}


	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: July, 2016
		Last Update on:	*/
	template <class T> std::pair<std::pair<double, double>, std::pair<int, int> > MultiStartGeneric(
		Assignment &solution, T &graph,
		const int heuristic, const int iterations, const bool seconds, const int K){
		std::cout << "Running Multi Start for Methods provided by Gutin and Karapetyan\n";
		clock_t comienzo = clock();
		double avg_time = 0, avg_solution = 0;
		int iteration = 0;
		for(iteration = 0; (!seconds && iteration < iterations) ||
			(seconds && (int)BasicMethods::GetRunningTime(comienzo) < iterations); iteration++){
			Assignment assignment;
			assignment.matching = BasicMethods::RandomizedSolution(graph.getDimensions());
			assignment.cost     = BasicMethods::CalculateCost(assignment.matching, graph);
			switch(heuristic){
				case 57:{
					CompleteDimensionwiseVariation(assignment, graph);
					CompleteKOpt(2, assignment, graph);
				}break;
				case 58:{
					CompleteMultiDimensionwiseVariation(assignment, graph);
					CompleteKOpt(2, assignment, graph);
				}break;
				case 59:{
					CompleteMultiDimensionwiseVariation(assignment, graph);
					VOpt(assignment, graph);
				}break;
				case 60:{
					MDVH2_3OPT(assignment, graph);
				}break;
				case 66:{
					CompleteDimensionwiseVariation(assignment, graph);
				}break;
				case 67:{
					CompleteMultiDimensionwiseVariation(assignment, graph);
				}break;
				case 68:{
					CompleteTwoOpt(assignment, graph);
				}break;
				case 69:{
					CompleteThreeOpt(assignment, graph);
				}break;
				case 70:{
					CompleteKOpt(K, assignment, graph);
				}break;
				default:break;
			}
//			std::cout <<"Iter "<<iteration << " "<<assignment.cost << " " << solution.cost <<"\n";
			avg_solution += assignment.cost;
			if(solution.cost > assignment.cost){
				solution = assignment;
				//std::cout <<"Iter "<<iteration << " "<<assignment.cost << " " << solution.cost <<"\n";
			}
		}
		avg_time = BasicMethods::GetRunningTime(comienzo) / (double) iteration;
		avg_solution /= (double) iteration;
		return std::make_pair(std::make_pair(avg_solution, avg_time), std::make_pair(iteration, solution.cost));
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on: July, 2016*/
	template <class T> std::pair<std::pair<double, double>, std::pair<int, int> > LocalMethodsGutin(const int localSearchType, const int iterations,
		const bool seconds, const int K_value, Assignment &assignment, T &graph){
		clock_t comienzo = clock();
		std::cout << "Start LocalSearchGutin\n";
		bool success = true;
		int successIterations = 0;
		int widthIter = BasicMethods::getWidth(iterations), widthCost = BasicMethods::getWidth(assignment.cost);
		std::pair<std::pair<double, double>, std::pair<int, int> > summary;
		for(int count = 1; (!seconds && count <= iterations) ||
			( seconds && (int)BasicMethods::GetRunningTime(comienzo) < iterations ); count++){
			success = false;
			switch(localSearchType){
				case 51:{
					success = DimensionwiseVariation(assignment, graph);
				}break;
				case 52:{
					success = MultiDimensionwiseVariation(assignment, graph);
				}break;
				case 53:{
					success = TwoOpt(assignment, graph);
				}break;
				case 54:{
					success = ThreeOpt(assignment, graph);
				}break;
				case 55:{
					success = VarOpt(K_value, assignment, graph);
				}break;
				case 56:{
					success = VOpt(assignment, graph);
				}break;
				case 61:{
					success = CompleteDimensionwiseVariation(assignment, graph);
				}break;
				case 62:{
					success = CompleteMultiDimensionwiseVariation(assignment, graph);
				}break;
				case 63:{
					success = CompleteTwoOpt(assignment, graph);
				}break;
				case 64:{
					success = CompleteThreeOpt(assignment, graph);
				}break;
				case 57: case 58: case 59: case 60: case 66: case 67: case 68: case 69: case 70:{
					summary = MultiStartGeneric(assignment, graph, localSearchType, iterations, seconds, K_value);
					count = iterations;
				}break;
				case 86: case 87: case 88: case 89: case 90:{
					IteratedLocalSearchGK(assignment, graph, localSearchType, iterations, seconds, K_value);
					count = iterations;
				}break;
				default:{
					count = iterations;
				}break;
			}
			if(success){
				printf("(%*d %*d %*d %lf)\n",
					widthIter, count, 2, localSearchType, widthCost, assignment.cost, BasicMethods::GetRunningTime(comienzo));
			}//*/
			else{
				break;
			}
		}
		return summary;
	}

};

#endif
