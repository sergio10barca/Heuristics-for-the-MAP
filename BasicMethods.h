#ifndef BASICMETHODS
#define BASICMETHODS
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <set>
#include "Output.h"
#include "LinearAssignment/StructureAP.h"
#include "LinearAssignment/HungarianMethod.h"
#include "LinearAssignment/FlowAssign.h"
#include "LinearAssignment/Auction.h"

namespace BasicMethods{
	double GetRunningTime(clock_t start){
		return ((clock()-start)/(double)CLOCKS_PER_SEC);
	}

	int getWidth(int value){
		int width = 0;
		do{
			width++;
			value/=10;
		}while(value>0);
		return width;
	}

	bool validateAssignment(std::vector<dTuple> &assignment){
		if(assignment.size() <= 0){
			std::cout << "Error: The assignment should not be null.\n";
			exit(0);
			//return false;
		}
		return true;
	}

	bool validateAssignmentSize(std::vector<dTuple> &assignment, unsigned const int assignmentSize){
		if(assignment.size() < assignmentSize){
			std::cout << "Error: The assignment should be at least of size "<< assignmentSize << " in order to apply this method.\n";
			return false;
		}
		return true;
	}

	template <class T> edgeCost CalculateCost(const std::vector<dTuple> &assignment, T &G){
		edgeCost assignmentCost = 0;
		for(unsigned int i = 0; i < assignment.size(); i++)
			assignmentCost += G.getEdgeCost(assignment[i]);
		return assignmentCost;
	}

	template <class T> bool verify(Assignment &solution, T &graph, const char* message){
		
		for(int dim = 0; dim < solution.matching[0].size(); dim++){
			std::vector<bool> used(solution.matching.size(), false);
			for(int j = 0; j < solution.matching.size(); j++){
				if(used[solution.matching[j][dim]] == true){
					std::cout << "Dimension " << dim << "\n";
					std::cout << "The vertex " << solution.matching[j][dim] << " is repeated\n";
					exit(0);
				}
				used[solution.matching[j][dim]] = true;
			}
		}
		if(CalculateCost(solution.matching, graph)  != solution.cost){
			printf("Error by calculating the new cost %s algorithm.\n", message);
			exit(0);
		}
		return true;
	}

	template <class T> void updateCostGraph(bool increase, int left, int right, Assignment &solution, T &G){
		if(increase){
			for(int i = left; i <= right; i++)
				solution.cost = solution.cost + G.getEdgeCost(solution.matching[i]);
		}
		else{
			for(int i = left; i <= right; i++)
				solution.cost = solution.cost - G.getEdgeCost(solution.matching[i]);
		}
		return;
	}

	std::vector<dTuple> RandomizedSolution(const std::vector<int> &dimensionSize){
//		std::cout << "Start RandomizedSolution\n";
		int minimumSize = INT_MAX, dimensions = (int)dimensionSize.size();
		for(int i = 0; i < dimensions; i++)
			minimumSize = std::min(minimumSize, dimensionSize[i]);

		std::vector<dTuple> assignment(minimumSize, dTuple(dimensions, 0));
		for(int i = 0; i < dimensions; i++){
			std::vector<int> permutation(dimensionSize[i], 0);
			for(int j = 0; j < dimensionSize[i]; j++)
				permutation[j] = j;
			std::random_shuffle(permutation.begin(), permutation.end());
			for(int j = 0; j < minimumSize; j++)
				assignment[j][i] = permutation[j];
		}
		std::sort(assignment.begin(), assignment.end());
//		Output::PrintSolution(assignment);//*/
//		std::cout << "End RandomizedSolution\n";
		return assignment;
	}

	std::vector<std::pair<std::pair<int, int>, int> > SolveAP(Vertex*& G, const int numberWomen, const int numberMen, const int option){
		switch(option){
			case 1:{
				Hungarian H;
		        return H.HungarianMethod(numberWomen, numberWomen, numberMen, G);//
			}break;
			case 2:{
				FlowAssign Fa;
				return Fa.Flow_assign(numberWomen, numberWomen, numberMen, G, 4);//
			}break;
			case 3:{
		        Auction A;
		        return A.AuctionAlgorithm(numberWomen, numberMen, G);//
			}break;
			default:break;
		}
		std::vector<std::pair<std::pair<int, int>, int> > v;
		return v;
	}

	bool AssignmentCompare(std::vector<dTuple> &a, std::vector<dTuple> &b){
		bool first = a == b, second;
		unsigned int i;
		for(i = 0; i < a.size() && i < b.size() && a[i] == b[i]; i++);
		if(i == a.size() && i == b.size())
			second = true;
		else
			second = false;
		if(first != second){
			std::cout << "Comparing error:"<<first <<" "<<second << "\n";
			exit(0);
		}
		return first;
	}

	template <class T> void GetBipartiteGraph(int dimA, int dimB, Vertex*& G, int &numberWomen, int &numberMen, std::vector<dTuple> &assignment,
		T &graph){
		numberWomen = graph.getDimension(dimA), numberMen =  graph.getDimension(dimB);

		G = new Vertex[numberWomen];
		for(int x = 0; x < numberWomen; x++)
			G[x].Allocate(numberMen);
//		printf("%d %d\n", numberWomen,  numberMen);

		for(unsigned int i = 0; i < assignment.size(); i++){
			int previousVertex = assignment[i][dimB];

			for(int vertex = 0; vertex < graph.getDimension(dimB); vertex++){
				assignment[i][dimB] = vertex;
				if(graph.isEdge(assignment[i]))
					G[assignment[i][dimA]].add(vertex, graph.getEdgeCost(assignment[i]));
			}
//			printf("-- %d %d %d\n",i, assignment[i][dimA], G[assignment[i][dimA]].length);
			assignment[i][dimB] = previousVertex;
		}
		return;
	}

	/*	Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: April, 2015
		Last Update on:	*/
	template <class T> void GetBipartiteGraph(const bool *F, Vertex*& G, const int s, const int n,
		std::vector<dTuple> &assignment, T &graph){
		G = new Vertex[n];
		for(int x = 0; x < n; x++)
			G[x].Allocate(n);
		dTuple v(s, 0);
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
//				printf("G[%d][%d]\n", i, j);
				for(int d = 0; d < s; d++)
					v[d] = F[d] ? assignment[j][d] : assignment[i][d];
				if(graph.isEdge(v))
					G[i].add(j, graph.getEdgeCost(v));
			}
		}
		return;
	}

	template <class T> void GetBipartiteSubGraph(const int dimB, Vertex*& G, const int nVertices,
		const int *p, const int *vertices, std::vector<dTuple> &assignment, T &graph){

		G = new Vertex[nVertices];
		for(int x = 0; x < nVertices; x++)
			G[x].Allocate(nVertices);

		for(int i = 0; i < nVertices; i++){
			const int &index = p[i];
			int previousVertex = assignment[index][dimB];

			for(int j = 0; j < nVertices; j++){
				assignment[index][dimB] = vertices[j];//vertices[j] is the actual vertex but j is the renamed
				if(graph.isEdge(assignment[index]))
					G[i].add(j, graph.getEdgeCost(assignment[index]));
			}
//			printf("-- %d\n", G[assignment[i][dimA]].length);
			assignment[index][dimB] = previousVertex;
		}
		return;
	}

	template <class T> void GetKPartiteGraph(const int *F, CompleteHyperGraph &G,  const int s, const int n,
		std::vector<dTuple> &assignment, T &graph, const int dimension_index, const int k, int *indices){
//		printf("%d %d %d\n", dimension_index, indices[dimension_index], k); 
		if(dimension_index == k){
			dTuple v(s, 0);
			for(int d = 0; d < s; d++){
				for(int index_k = 0; index_k < k; index_k++){
					if(F[d] == index_k){
						v[d] = assignment[indices[index_k]][d];
						break;
					}
				}
			}
/*			for(int i = 0; i < v.size(); i++)
				printf("%d%c", v[i], i+1==v.size() ?'\n' : ' ');//*/
			if(graph.isEdge(v)){
				G.addEdge(graph.getEdgeCost(v));
			}
			else{
				std::cout << "We have not a k-partite complete graph\n";
				exit(0);
			}
		}
		else{
			for(int i = 0; i < n; i++){
				indices[dimension_index] = i;	
				BasicMethods::GetKPartiteGraph(F, G, s, n, assignment, graph, dimension_index+1, k, indices);
			}
		}
		return;
	}


	template <class T> void GetTripartiteGraph(const int *F, CompleteHyperGraph &G, const int s, const int n,
		std::vector<dTuple> &assignment, T &graph){
		dTuple v(s, 0);
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				for(int k = 0; k < n; k++){
					for(int d = 0; d < s; d++)
						v[d] = F[d] == 0 ? assignment[i][d] : (F[d] == 1 ? assignment[j][d] : assignment[k][d]);
	/*				for(int z = 0; z < v.size(); z++)
						printf("%d%c", v[z], z+1==v.size() ?'\n' : ' ');//*/

					if(graph.isEdge(v))
						G.addEdge(graph.getEdgeCost(v));
					else{
						std::cout << "We have not a tripartite complete graph\n";
						exit(0);
					}
				}
		return;
	}

	template <class T> void GetTripartiteSubGraph(const int *F, CompleteHyperGraph &G, const int s, const int n,
		const int *vertices, const std::vector<dTuple> &assignment, T &graph){
		dTuple v(s, 0);
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++){
					for(int d = 0; d < s; d++)
						v[d] = F[d] == 0 ? assignment[vertices[i]][d] : (F[d] == 1 ? assignment[vertices[j]][d] : assignment[vertices[k]][d]);
					if(graph.isEdge(v))
						G.addEdge(graph.getEdgeCost(v));
					else{
						std::cout << "We have not a tripartite complete graph\n";
						exit(0);
					}
				}
			}
		}
		return;
	}

	template <class T> void Get3DimSubGraph(const int *dim, const int nVertices, const int *p, const int *vertices_dim1, const int *vertices_dim2,
		CompleteHyperGraph &G, std::vector<dTuple> &assignment, T &graph){

		for(int i = 0; i < nVertices; i++){
			const int &index = p[i];
			int previousVertexDim1 = assignment[index][dim[1]];

			for(int j = 0; j < nVertices; j++){
				int previousVertexDim2 = assignment[index][dim[2]];
				assignment[index][dim[1]] = vertices_dim1[j];//vertices[j] is the actual vertex but j is the renamed
				for(int k = 0; k < nVertices; k++){
					assignment[index][dim[2]] = vertices_dim2[k];//vertices[j] is the actual vertex but j is the renamed
					if(graph.isEdge(assignment[index]))
						G.addEdge(graph.getEdgeCost(assignment[index]));
					else
						G.addEdge(INT_MAX);
				}
				assignment[index][dim[2]] = previousVertexDim2;
			}
			assignment[index][dim[1]] = previousVertexDim1;
		}	//*/
		return;
	}

	int next_power(const int N){
		int zeros = 0, i, n = N, pot = 1;
		for(i = 0; n > 0; n>>=1, i++, pot*=2){
			zeros += !(n&1);
		}
		if(zeros == i-1){
			return N;
		}
		return pot;
	}

	void Perturbation(Assignment &solution, const int size, const int dimension){
//		printf("Perturbation\n");
		int n = solution.matching.size();
		int *v;
		v = new int[n];
		for(int i = 0; i < n; i++){
			v[i] = i;
		}
		std::random_shuffle(v, v+n);
/*		for(int i = 0; i < size; i++)
			printf("%d%c", v[i], i+1 == size ? '\n' : ' ');//*/
		int first = solution.matching[ v[0] ][ dimension ];
		for(int i = 1; i < size; i++){
			solution.matching[ v[i-1] ][ dimension ] = solution.matching[ v[i] ][ dimension ];

		}
		solution.matching[ v[ size - 1 ] ][ dimension ] = first;
		delete []v;
//		printf("End Perturbation\n");
		return;
	}

};

#endif
