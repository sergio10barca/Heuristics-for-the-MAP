#ifndef LOCALMETHODS2D
#define LOCALMETHODS2D
#include <algorithm>
#include <cstdlib>
#include <map>
#include <vector>
#include "Structures.h"
#include "Output.h"

namespace LocalMethods2D{
	/*
		Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:
	*/
	template <class T> bool RoundLocal2Dim(const int dimA, const int dimB, Assignment &solution, T &graph){
//		cout << "Start RoundLocal" << endl;
		int numberWomen = 0, numberMen = 0, previousCost = solution.cost;
		Vertex *G = NULL;

		BasicMethods::GetBipartiteGraph(dimA, dimB, G, numberWomen, numberMen, solution.matching, graph);
		std::vector<std::pair<std::pair<int, int>, int> > matching;

		matching = BasicMethods::SolveAP(G, numberWomen, numberMen, 3);
		solution.cost = 0;
		for(unsigned int i = 0; i < matching.size(); i++)
			for(unsigned int j = 0; j < solution.matching.size(); j++){
				if(solution.matching[j][dimA] == matching[i].first.first){
					solution.matching[j][dimB] = matching[i].first.second;
					solution.cost += matching[i].second;
					break;
				}				
			}
//		std::cout << "The matching cost from hungarian was = "<<solution.cost << " "<<numberWomen<<"\n";
		delete[] G;
		BasicMethods::verify(solution, graph, "Error en RoundLocal");
		return previousCost > solution.cost;
	}

	/*
		Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: May, 2015
		Last Update on:
	*/
	template <class T> bool RoundLocal2DimSubSet(const int dimB, Assignment &solution, T &graph){
//		std::cout << "Start RoundLocal2DimSubSet";
		int nVertices = (rand()%(solution.matching.size()-1))+2, *permutation = new int[solution.matching.size()], *original_vertices;
		int previousCost = solution.cost;
//		printf("tamano solution = %d -> tamano subset = %d dims %d %d\n", solution.matching.size(), nVertices, dimA, dimB);

		for(unsigned int i = 0; i < solution.matching.size(); i++)
			permutation[i] = i;
		std::random_shuffle(permutation, permutation+solution.matching.size());

		original_vertices = new int[nVertices];
		std::sort(permutation, permutation+nVertices);

		for(int i = 0; i < nVertices; i++){
			original_vertices[i] = solution.matching[permutation[i]][dimB];
//			printf("(%d,%d)%c", permutation[i], original_vertices[i], i == nVertices-1 ? '\n' : ' ');
		}

		Vertex *G = NULL;

		BasicMethods::GetBipartiteSubGraph(dimB, G, nVertices, permutation, original_vertices, solution.matching, graph);
		std::vector<std::pair<std::pair<int, int>, int> > matching;
		matching = BasicMethods::SolveAP(G, nVertices, nVertices, 3);

//		printf("Solution matching %d \n", solution.cost);
		for(unsigned int i = 0; i < matching.size(); i++){
//			printf("%d %d %d\n", permutation[matching[i].first.first], matching[i].first.first, matching[i].first.second);
			solution.cost -= graph.getEdgeCost(solution.matching[permutation[matching[i].first.first]]);
			solution.matching[permutation[matching[i].first.first]][dimB] = original_vertices[matching[i].first.second];
			solution.cost += matching[i].second;
		}
//		printf("Solution new matching cost %d \n", solution.cost);
//		cout << "The matching cost from hungarian was = "<<solution.cost << " "<<numberWomen<<endl;
		delete[] G;
		delete[] permutation;
		delete[] original_vertices;

		BasicMethods::verify(solution, graph, "Error en RoundLocal2DimSubSet");
		return previousCost > solution.cost;
	}

	/* */
	template <class T> bool RoundLocal2DimAtRandom(Assignment &solution, T &graph){
//		cout << "Start RoundLocalAtRandom" << endl;

		BasicMethods::validateAssignment(solution.matching);
		if(!BasicMethods::validateAssignmentSize(solution.matching, 2))
			return false;

		if(graph.getNumberOfDimensions() < 2){
			std::cout << "The method cannot be applied.\n";
			exit(0);
		}

		//chosening the dimensions to fix
		int fixedDimensionA, fixedDimensionB;
		fixedDimensionA = rand()%graph.getNumberOfDimensions();
		while((fixedDimensionB = rand()%graph.getNumberOfDimensions()) == fixedDimensionA);
		
		return RoundLocal2Dim(std::min(fixedDimensionA, fixedDimensionB), std::max(fixedDimensionA, fixedDimensionB), solution, graph);
	}

	template <class T> bool MultiRoundLocal2Dim(Assignment &solution, T &graph){
//		cout << "Start MultiRoundLocal" << endl;

		BasicMethods::validateAssignment(solution.matching);
		if(!BasicMethods::validateAssignmentSize(solution.matching, 2))
			return false;

		//chosening the dimensions to fix
		if(graph.getNumberOfDimensions() < 2){
			std::cout << "The method cannot be applied.\n";
			return false;
		}

		int cont;
		bool success = true;
		for(cont = 0; success; cont++){
			success = false;
			for(int dimA = 0; dimA < graph.getNumberOfDimensions(); dimA++)
				for(int dimB = dimA+1; dimB < graph.getNumberOfDimensions(); dimB++){
//					cout << dimA << " " << dimB << endl;
					success |= RoundLocal2Dim(dimA, dimB, solution, graph);
				}
		}
		return cont > 1;
	}

	template <class T> bool MultiRoundLocal2DimAtRandom(Assignment &solution, T &graph){
//		cout << "Start MultiRoundLocalAtRandom" << endl;

		BasicMethods::validateAssignment(solution.matching);
		if(!BasicMethods::validateAssignmentSize(solution.matching, 2))
			return false;

		//chosening the dimensions to fix
		if(graph.getNumberOfDimensions() < 2){
			std::cout << "The method cannot be applied.\n";
			return false;
		}

		int cont;
		bool success = true;
		std::vector<std::pair<int, int> > dims;
		for(int dimA = 0; dimA < graph.getNumberOfDimensions(); dimA++)
			for(int dimB = dimA+1; dimB < graph.getNumberOfDimensions(); dimB++)
				dims.push_back(std::make_pair (dimA, dimB));
		for(cont = 0; success; cont++){
			success = false;
			std::random_shuffle(dims.begin(), dims.end());
			for(unsigned int i = 0; i < dims.size(); i++){
//				std::cout << dims[i].first << " + " << dims[i].second << "\n";
				success |= RoundLocal2Dim(dims[i].first, dims[i].second, solution, graph);
			}
		}
		return false;
	}

};

#endif