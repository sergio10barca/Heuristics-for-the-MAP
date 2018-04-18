#ifndef LOCALMETHODS3D
#define LOCALMETHODS3D
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>
#include "Structures.h"
#include "Output.h"
#include "BasicMethods.h"
#include "ExactMethods.h"
#include "LocalMethodsGutin.h"

namespace LocalMethods3D{

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on: */
	template <class T> void UpdateSolution(const int *F, Assignment &sol3D, Assignment &solution, CompleteHyperGraph &graph3D, T &graph){
		int n = (int)sol3D.matching.size(), s = graph.getNumberOfDimensions(), &cost = solution.cost;
		std::vector<dTuple> aux = solution.matching;
		cost = 0;
		for(int i = 0; i < n; i++){
			for(int d = 0; d < s; d++)
				solution.matching[i][d] = aux[sol3D.matching[i][F[d]]][d];
			cost += graph3D.getEdgeCost(sol3D.matching[i]);
		}
		std::sort(solution.matching.begin(), solution.matching.end());
		BasicMethods::verify(solution, graph, "UpdateSolution");
		return;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on: */
	template <class T> void UpdateSolution3D(const int *F, const int *vertices, Assignment &sol3D, Assignment &solution,
		CompleteHyperGraph &graph3D, T &graph){
		int n = (int)sol3D.matching.size(), s = graph.getNumberOfDimensions();
		std::vector<dTuple> aux = solution.matching;
		for(int i = 0; i < n; i++)
			solution.cost -= graph.getEdgeCost(solution.matching[vertices[i]] );

		for(int i = 0; i < n; i++){
			for(int d = 0; d < s; d++)
				solution.matching[vertices[i]][d] = aux[ vertices[sol3D.matching[i][F[d]]] ][d];
			solution.cost += graph3D.getEdgeCost(sol3D.matching[i]);
		}

		BasicMethods::verify(solution, graph, "UpdateSolution");
		return;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:	*/
	template <class T> void UpdateSolutionMD(const int *dimensions, const int *hedges, Assignment &solMD, Assignment &solution,
		CompleteHyperGraph &graphMD, T &graph){
		int n = (int)solMD.matching.size(), s = graph.getNumberOfDimensions();
		std::vector<dTuple> aux = solution.matching;
		for(int i = 0; i < n; i++)
			solution.cost -= graph.getEdgeCost(solution.matching[hedges[i]] );

		for(int i = 0; i < n; i++){
			for(int d = 0; d < s; d++)
				if(dimensions[d] == 1)
					solution.matching[hedges[i]][d] = aux[ hedges[solMD.matching[i][d]] ][d];
			solution.cost += graphMD.getEdgeCost(solMD.matching[i]);
		}

		BasicMethods::verify(solution, graph, "UpdateSolution");
		return;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on:	*/
	template <class T> bool DimensionwiseVariation3D(Assignment &solution, T &graph){
//		std::cout << "Running DimensionwiseVariation3D\n";
		int n = graph.getDimension(0), previousCost = solution.cost, *F = new int[graph.getNumberOfDimensions()];
		for(int i = 0; i < graph.getNumberOfDimensions(); i++)
			F[i] = 2;
		for(int dim1 = 0; dim1 < graph.getNumberOfDimensions(); dim1++){
			F[dim1] = 0;
			for(int dim2 = dim1+1; dim2 < graph.getNumberOfDimensions(); dim2++){
				CompleteHyperGraph graph3D(3, n);
				F[dim2] = 1;
/*				for(int i = 0; i < graph.getNumberOfDimensions(); i++)
					printf("%d%c", F[i], i+1 == graph.getNumberOfDimensions() ? '\n' : ' ');//*/
				BasicMethods::GetTripartiteGraph(F, graph3D, graph.getNumberOfDimensions(), n, solution.matching, graph);
//				std::cout << "generate graph "<<dim1 <<" "<<dim2 << "\n";
				Assignment answer = ExactMethods::PAM_GurobiSolver(graph3D);
//				Output::PrintSolution(answer, graph3D);
				UpdateSolution(F, answer, solution, graph3D, graph);
				F[dim2] = 2;
				if(graph.getNumberOfDimensions() == 3){
					delete[] F;
					return false;
				}
			}
			F[dim1] = 2;
		}
		delete[] F;
		return previousCost > solution.cost;
	}

	template <class T> void GenerateSomeVersusRest(int index, int prev_dim, int K, int *F, Assignment &solution, T &graph){
		if(index == K){
			int n = solution.matching.size();
			CompleteHyperGraph graphKD(K+1, n);
			int *indices = new int[K+2];
/*			for(int i = 0; i < graph.getNumberOfDimensions(); i++)
				printf("%d%c", F[i], i+1 == graph.getNumberOfDimensions() ? '\n' : ' ');//*/
			BasicMethods::GetKPartiteGraph(F, graphKD, graph.getNumberOfDimensions(), n, solution.matching, graph,
				0, K+1, indices);//*/
			Assignment answer = ExactMethods::PAM_GurobiSolver(graphKD);
			UpdateSolution(F, answer, solution, graphKD, graph);
		}
		else{
			for(int dim = prev_dim + 1; dim < graph.getNumberOfDimensions(); dim++){
				F[dim] = index;
				GenerateSomeVersusRest(index + 1, dim, K, F, solution, graph);
				F[dim] = K;
			}
		}
		return;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on: */
	int getIndex(std::vector<int> F, const int base){
		int value = 0;
		for(int i = 0; i < F.size(); i++)
			value = (value*base) + F[i];
		return value;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:	*/
	int getSDValue(std::vector<int> F, const int sets){
		int *cont = new int[sets];
		for(int i = 0; i < sets; i++)
			cont[i] = 0;
		for(int i = 0; i < (int)F.size(); i++)
			cont[F[i]]++;
		for(int i = 0; i < sets; i++)
			if(cont[i] == 0)
				return -1;
		return getIndex(F, sets);
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:*/
	template <class T> void Generate3DCombinations(const int max_dim, const int index, const int next,
		std::vector<int> F, std::set<int> &visited, Assignment &solution, T &graph){
		if(index == 0){
			std::reverse(F.begin(), F.end());
			int combination = getSDValue(F, 3);
			if(visited.find(combination) == visited.end()){
				visited.insert(combination);
/*				for(int i = 0; i < F.size(); i++)
					printf("%d%c", F[i], i+1 == F.size() ? '\n' : ' ');//*/
				int n = (int)solution.matching.size();
				CompleteHyperGraph graph3D(3, n);
				int *F1 = new int[F.size()], *indices = new int[F.size()];
				for(int i = 0; i < F.size(); i++)
					F1[i] = F[i];
				BasicMethods::GetKPartiteGraph(F1, graph3D, graph.getNumberOfDimensions(), n, solution.matching, graph,
					0, 3, indices);//*/
//				BasicMethods::GetTripartiteGraph(F1, graph3D, graph.getNumberOfDimensions(), n, solution.matching, graph);
				Assignment answer = ExactMethods::PAM_GurobiSolver(graph3D);
				UpdateSolution(F1, answer, solution, graph3D, graph);
				delete[] F1;
				delete[] indices;
			}
		}
		else{
			for(int i = 0; i <= std::min(next, max_dim); i++){
				F[index-1] = i;
				Generate3DCombinations(max_dim, index-1, std::max(i+1, next), F, visited, solution, graph);
			}
		}
		return;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:	*/
	template <class T> bool MultiDimensionwiseVariation3D(Assignment &solution, T &graph){
//		std::cout << "Running MultiDimensionwiseVariation3D\n";
		int previousCost = solution.cost, K = 3;
		std::vector<int> F(graph.getNumberOfDimensions(), 0);
		std::set<int> visited;
		visited.insert(-1);
		Generate3DCombinations(K-1, graph.getNumberOfDimensions()-1, 1, F, visited, solution, graph);
//		std::cout << "Nodes = " << visited.size() << "\n";
		return previousCost > solution.cost && graph.getNumberOfDimensions() != 3;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:*/
	template <class T> void GenerateKDCombinations(const int max_dim, const int index, const int next,
		std::vector<int> F, std::set<int> &visited, Assignment &solution, T &graph){
		if(index == 0){
			std::reverse(F.begin(), F.end());
			int K = graph.getNumberOfDimensions();
			int combination = getSDValue(F, K - 1);
			if(visited.find(combination) == visited.end()){
				visited.insert(combination);
/*				for(int i = 0; i < F.size(); i++)
					printf("%d%c", F[i], i + 1 == F.size() ? '\n' : ' ');//*/
				int n = (int)solution.matching.size();
				CompleteHyperGraph graphKD(K - 1, n);
				int *F1 = new int[F.size()], *indices = new int[F.size()];
				for(int i = 0; i < F.size(); i++)
					F1[i] = F[i];
				BasicMethods::GetKPartiteGraph(F1, graphKD, graph.getNumberOfDimensions(), n, solution.matching, graph,
					0, K - 1, indices);//*/
				Assignment answer = ExactMethods::PAM_GurobiSolver(graphKD);
				UpdateSolution(F1, answer, solution, graphKD, graph);
				delete[] F1;
				delete[] indices;
			}
		}
		else{
			for(int i = 0; i <= std::min(next, max_dim); i++){
				F[index-1] = i;
				GenerateKDCombinations(max_dim, index-1, std::max(i+1, next), F, visited, solution, graph);
			}
		}
		return;
	}


	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: July, 2016
		Last Update on:	*/
	template <class T> bool CompleteDimensionwiseVariation3D(Assignment &solution, T &graph){
//		std::cout << "Running Complete DimensionwiseVariation3D\n";
		int previousCost;
		do{
			previousCost = solution.cost;
			DimensionwiseVariation3D(solution, graph);
//			std::cout << previousCost << " " << solution.cost << "\n";
		}while(previousCost > solution.cost && graph.getNumberOfDimensions() != 3);
		return false;
	}


	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: July, 2016
		Last Update on:	*/
	template <class T> bool CompleteMultiDimensionwiseVariation3D(Assignment &solution, T &graph){
//		std::cout << "Running CompleteMultiDimensionwiseVariation3D\n";
		int previousCost;
		do{
			previousCost = solution.cost;
			MultiDimensionwiseVariation3D(solution, graph);
//			std::cout << previousCost << " " << solution.cost << "\n";
		}while(previousCost > solution.cost && graph.getNumberOfDimensions() != 3);
		return false;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on:
	*/
	template <class T> void VariableKOpt(const int dim, const int s, const int K, const int *indices, const int n, const int index,
		Assignment &solution, CompleteHyperGraph &graphMD, T &graph){
		if(dim == s){
			graphMD.addEdge(graph.getEdgeCost(index));
		}
		else{
			for(int i = 0; i < K; i++)
				VariableKOpt(dim+1, s, K, indices, n, (index * n) + solution.matching[indices[i]][dim], solution, graphMD, graph);
		}
		return;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:	*/
	template <class T> void KOpt(const int index, const int prev, const int K, const int n, int *indices, Assignment &solution,
		CompleteHyperGraph &graphMD, T &graph){
		if(index == K){
			int cost = 0;
			for(int i = 0; i < K; i++)
				cost += graph.getEdgeCost(solution.matching[indices[i]]);
			graphMD.setDimensions(graph.getNumberOfDimensions(), K);
			VariableKOpt(0, graph.getNumberOfDimensions(), K, indices, n, 0, solution, graphMD, graph);
			Assignment A = ExactMethods::PAM_GurobiSolver(graphMD);
			Assignment aux = solution;
			for(int i = 0; i < K; i++)
				for(int j = 0; j < graph.getNumberOfDimensions(); j++)
					solution.matching[indices[i]][j] = aux.matching[indices[A.matching[i][j]]][j];
			solution.cost += (A.cost - cost);
		}
		else{
			for(int i = prev; i < n; i++){
				indices[index] = i;
				KOpt(index+1, i+1, K, n, indices, solution, graphMD, graph);
			}
		}
		return;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on: */
	template <class T> bool VarOptGurobi(const int K, Assignment &solution, T &graph){
//		std::cout << "Running "<<K<<"-OptGurobi\n";
		int n = graph.getDimension(0), previousCost = solution.cost, *indices;
		indices = new int[K];
		CompleteHyperGraph graphMD(graph.getNumberOfDimensions(), K);
		KOpt(0, 0, K, n, indices, solution, graphMD, graph);
		BasicMethods::verify(solution, graph, "K-OptGutin");
		delete[] indices;
		return previousCost > solution.cost;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on: */
	template <class T> bool CompleteKOptGurobi(const int K, Assignment &solution, T &graph){
		std::cout << "Running Complete"<<K<<"OptGurobi\n";
		int previousCost;
		do{
			previousCost = solution.cost;
			VarOptGurobi(K, solution, graph);
		}while(previousCost > solution.cost);
		return true;
	}

	/*
		Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:
	*/
	template <class T> void VariableDimensions(const int dim, const int s, const int K, const int fixedHEdge,
		const int *indices, const int *dimensions, const int n, const int index,
		Assignment &solution, CompleteHyperGraph &graphMD, T &graph){
		if(dim == s){

			graphMD.addEdge(graph.getEdgeCost(index));
//			dTuple tuple = graph.getEdge(index);
//			Output::PrintDTuple(tuple);
//			std::cout << graph.getEdgeCost(index) << std::endl;
		}
		else{
			if(dimensions[dim] == 1){
				for(int i = 0; i < K; i++)
					VariableDimensions(dim+1, s, K, fixedHEdge, indices, dimensions, n, (index * n) + solution.matching[indices[i]][dim],
					solution, graphMD, graph);
			}
			else{
				VariableDimensions(dim+1, s, K, fixedHEdge, indices, dimensions, n, (index * n) + solution.matching[fixedHEdge][dim],
					solution, graphMD, graph);
			}
		}
		return;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: July, 2016
		Last Update on:	*/
	template <class T> std::pair< std::pair<double, double>, std::pair<int, int> > MultiStartGeneric(
		Assignment &solution, T &graph,
		const int heuristic, const int iterations, const bool seconds, const int K){
		clock_t comienzo = clock();
		std::cout << "Running Multi Start Dimension Wise Variation\n";
		double avg_time = 0, avg_solution = 0;
		int iteration = 0;
		for(iteration = 0; (!seconds && iteration < iterations) ||
			( seconds && (int)BasicMethods::GetRunningTime(comienzo) < iterations ); iteration++){
			Assignment assignment;
			assignment.matching = BasicMethods::RandomizedSolution(graph.getDimensions());
			assignment.cost     = BasicMethods::CalculateCost(assignment.matching, graph);
			switch(heuristic){
				case 11:{
					CompleteMultiDimensionwiseVariation3D(assignment, graph);
					LocalMethodsGutin::CompleteKOpt(2, assignment, graph);
				}break;
				case 12:{
					int cost = assignment.cost;
					do{
						CompleteMultiDimensionwiseVariation3D(assignment, graph);
						if(cost <= assignment.cost){
							break;
						}
						std::cout << "sigue 1----------\n";
						cost = assignment.cost;
						LocalMethodsGutin::CompleteKOpt(3, assignment, graph);
						if(cost <= assignment.cost){
							break;
						}
						std::cout << "sigue 2----------\n";
						cost = assignment.cost;
					}while(true);
				}break;
				case 46:{
					CompleteDimensionwiseVariation3D(assignment, graph);
				}break;
				case 47:{
					CompleteMultiDimensionwiseVariation3D(assignment, graph);
				}break;
				case 48:{
					CompleteKOptGurobi(2, assignment, graph);
				}break;
				case 49:{
					CompleteKOptGurobi(3, assignment, graph);
				}break;
				case 50:{
					CompleteKOptGurobi(K, assignment, graph);
				}break;
				default:break;
			}
			avg_solution += assignment.cost;
//			std::cout <<"Iter "<< iteration << " "<<assignment.cost << " " << solution.cost <<"\n";
			if(solution.cost > assignment.cost){
				solution = assignment;
//				std::cout <<"Iter "<< iteration << " "<<assignment.cost << " " << solution.cost <<"\n";
			}
		}
		avg_time = BasicMethods::GetRunningTime(comienzo) / (double) iteration;
		avg_solution /= (double) iteration;
		return std::make_pair(std::make_pair(avg_solution, avg_time), std::make_pair(iteration, solution.cost));
	}
};

#endif

