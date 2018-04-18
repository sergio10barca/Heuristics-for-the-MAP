#ifndef LOCALMETHODS
#define LOCALMETHODS
#include <algorithm>
#include <cstdlib>
#include <map>
#include <vector>
#include "Structures.h"
#include "Output.h"
#include "BasicMethods.h"
#include "LocalMethods2D.h"
#include "LocalMethods3D.h"
#include "LocalMethodsGutin.h"

namespace LocalMethods{

	void RightRotation(const int dim, const int startPosition, const int endPosition, Assignment &solution){
		const int lastVertex = solution.matching[endPosition][dim];
		for(int index = endPosition-1; index >= startPosition; index--){
			solution.matching[index+1][dim] = solution.matching[index][dim];
		}
		solution.matching[startPosition][dim] = lastVertex;
		return;
	}//*/

	void LeftRotation(const int dim, const int startPosition, const int endPosition, Assignment &solution){
		const int firstVertex = solution.matching[startPosition][dim];
		for(int index = startPosition+1; index <= endPosition; index++){
			solution.matching[index-1][dim] = solution.matching[index][dim];
		}
		solution.matching[endPosition][dim] = firstVertex;
		return;
	}//*/

	/*
	Created By: Sergio Perez
	Modified by: Sergio Perez
	Created on: May, 2014
	Last Update on:*/
	template <class T> bool CircularRotation(const int iter, Assignment &solution, T &graph){
//		std::cout << "Start CircularRotation\n";
		bool improved = true;
		int iterations = 1;
		while(improved == true && iterations <= iter){
			improved = false;
			iterations++;
			for(int chosenDim = 0; chosenDim < graph.getNumberOfDimensions(); chosenDim++){
				for(int startPosition = 0; startPosition < solution.matching.size(); startPosition++){
					for(int endPosition = startPosition+1; endPosition < solution.matching.size(); endPosition++){
						const int vertex_orig = solution.matching[startPosition][chosenDim];
						Assignment A = solution;
						do{
							BasicMethods::updateCostGraph(false, startPosition, endPosition, A, graph);
							LocalMethods::LeftRotation(chosenDim, startPosition, endPosition, A);
							BasicMethods::updateCostGraph(true, startPosition, endPosition, A, graph);

							if(A.cost < solution.cost){
								solution = A;
								improved = true;
							}
						}while(vertex_orig != A.matching[startPosition][chosenDim]);
					}
				}
			}
			BasicMethods::verify(solution, graph, "Error en CircularRotation");
//			std::cout << "Iteration " << (iterations-1) << " cost " << solution.cost << "\n";
		}
//		std::cout << "End CircularRotation "<<(iterations-1) <<"\n";
		return improved;
	}

	void ReverseVertices(const int dim, const int startPosition, const int endPosition, Assignment &solution){
		for(int left = startPosition, right = endPosition; left < right; left++, right--){
			std::swap(solution.matching[left][dim], solution.matching[right][dim]);
		}
		return;
	}

	/*
	Created By: Sergio Perez
	Modified by: Sergio Perez
	Created on: May, 2014
	Last Update on:*/
	template <class T> bool Inversion(const int iter, Assignment &solution, T &graph){
//		std::cout << "Start Inversion\n";
		BasicMethods::validateAssignment(solution.matching);
		if(!BasicMethods::validateAssignmentSize(solution.matching, 2))
			return false;

		bool improved = true;
		int iterations = 1;
		while(improved == true && iterations <= iter){
			improved = false;
			iterations++;
			for(int chosenDim = 0; chosenDim < graph.getNumberOfDimensions(); chosenDim++){
				for(int startPosition = 0; startPosition < solution.matching.size(); startPosition++){
					for(int endPosition = startPosition+1; endPosition < solution.matching.size(); endPosition++){
						edgeCost previousCost = solution.cost;

						BasicMethods::updateCostGraph(false, startPosition, endPosition, solution, graph);
						LocalMethods::ReverseVertices(chosenDim, startPosition, endPosition, solution);
						BasicMethods::updateCostGraph(true, startPosition, endPosition, solution, graph);
						if(solution.cost >= previousCost){
							BasicMethods::updateCostGraph(false, startPosition, endPosition, solution, graph);
							LocalMethods::ReverseVertices(chosenDim, startPosition, endPosition, solution);
							BasicMethods::updateCostGraph(true, startPosition, endPosition, solution, graph);
						}else{
							improved = true;
						}
					}
				}
			}
			BasicMethods::verify(solution, graph, "Error en Inversion");
//			std::cout << "Iteration " << (iterations-1) << " cost " << solution.cost << "\n";
		}
//		std::cout << "End Inversion "<<(iterations-1) <<"\n";
		return improved;
	}

	/*
	Created By: Sergio Perez
	Modified by: Sergio Perez
	Created on: May, 2014
	Last Update on://*/
	template <class T> bool TwoOpt(const int iter, Assignment &solution, T &graph){
//		std::cout << "Start TwoOpt\n";
		BasicMethods::validateAssignment(solution.matching);
		if(!BasicMethods::validateAssignmentSize(solution.matching, 2))
			return false;

		bool improved = true;
		int iterations = 1;
		while(improved == true && iterations <= iter){
			improved = false;
			iterations++;
			for(int chosenDim = 0; chosenDim < graph.getNumberOfDimensions(); chosenDim++){
				for(int startPosition = 0; startPosition < solution.matching.size(); startPosition++){
					for(int endPosition = startPosition+1; endPosition < solution.matching.size(); endPosition++){
						edgeCost previousCost = solution.cost;

						BasicMethods::updateCostGraph(false, startPosition, endPosition, solution, graph);
						std::swap(solution.matching[startPosition][chosenDim], solution.matching[endPosition][chosenDim]);
						BasicMethods::updateCostGraph(true, startPosition, endPosition, solution, graph);
						if(solution.cost >= previousCost){
							BasicMethods::updateCostGraph(false, startPosition, endPosition, solution, graph);
							std::swap(solution.matching[startPosition][chosenDim], solution.matching[endPosition][chosenDim]);
							BasicMethods::updateCostGraph(true, startPosition, endPosition, solution, graph);
						}else{
							improved = true;
						}
					}
				}
			}
			BasicMethods::verify(solution, graph, "Error en TwoOpt");
//			std::cout << "Iteration " << (iterations-1) << " cost " << solution.cost << "\n";
		}
//		std::cout << "End TwoOpt "<<(iterations-1) <<"\n";
		return improved;
	}

	/*
	Created By: Sergio Perez
	Modified by: Sergio Perez
	Created on: May, 2014
	Last Update on:
	This functions takes O(numberOfDimensions|solution|^2) time for independent costs//*/
	template <class T> bool Kvertex(const int iter, const int k, Assignment &solution, T &graph){
//		std::cout << "Start Kvertex "<< k <<"\n" ;
		BasicMethods::validateAssignment(solution.matching);
		if(!BasicMethods::validateAssignmentSize(solution.matching, 2))
			return false;

		bool improved = true;
		int iterations = 1;
		while(improved == true && iterations <= iter){
			improved = false;
			iterations++;
			for(int chosenDim = 1; chosenDim < graph.getNumberOfDimensions(); chosenDim++){
				for(int startPosition = 0; startPosition <= (solution.matching.size() - k); startPosition++){
					const int endPosition = startPosition + k - 1;
					edgeCost previousCost = solution.cost;

					std::vector<int> permutation(k, 0), best_perm;
					std::vector<int> original(endPosition - startPosition + 1, 0);
					for(int i = 0; i < permutation.size(); i++){
						permutation[i] = i;
						original[i] = solution.matching[ i + startPosition ][chosenDim];
					}
					best_perm = permutation;
					do{
						BasicMethods::updateCostGraph(false, startPosition, endPosition, solution, graph);
						for(int j = 0; j < original.size(); j++){
							solution.matching[ startPosition + j ][chosenDim] = original[ permutation[j] ];
						}
						BasicMethods::updateCostGraph(true, startPosition, endPosition, solution, graph);
						BasicMethods::verify(solution, graph, "Error en Kvertex");
						if(solution.cost < previousCost){
							best_perm = permutation;
							previousCost = solution.cost;
							improved = true;
						}
					}while(next_permutation(permutation.begin(), permutation.end()));
					BasicMethods::updateCostGraph(false, startPosition, endPosition, solution, graph);
					for(int j = 0; j < original.size(); j++){
						solution.matching[ startPosition + j ][chosenDim] = original[ best_perm[j] ];
					}
					BasicMethods::updateCostGraph(true, startPosition, endPosition, solution, graph);
				}
			}
			BasicMethods::verify(solution, graph, "Error en Kvertex");
//			std::cout << "Iteration " << (iterations-1) << " cost " << solution.cost << "\n";
		}
//		std::cout << "End kvertex "<<(iterations-1) <<"\n";
		return improved;
	}


	/*
	Created By: Sergio Perez
	Modified by: Sergio Perez
	Created on: May, 2015
	Last Update on:	*/
	template <class T> bool SelectAllBasic(const int iterations, const int K_value, Assignment &solution, T &graph){
//		std::cout << "Start BasicAllBasic\n";
		bool improved = true;
		for(int i = 1; improved == true && i <= iterations; i++){
			improved = false;
			improved |= TwoOpt(1, solution, graph);
			improved |= Inversion(1, solution, graph);
			improved |= CircularRotation(1, solution, graph);
			improved |= Kvertex(1, K_value, solution, graph);
//			std::cout << "iteration " << i << "\n";
		}
//		std::cout << "End BasicAllBasic\n";
		return improved;
	}

	template <class T> std::pair< std::pair<double, double>, std::pair<int, int> > MultiStartBasicLSH(
		Assignment &solution, T &graph,
		const int heuristic, const int iterations, const bool seconds, const int k_value){
		std::cout << "Running Multi Start Basic LSH\n";
		clock_t comienzo = clock();
		double avg_time = 0, avg_solution = 0;
		int iteration = 0;
		for(iteration = 0; (!seconds && iteration < iterations) ||
			(seconds && (int)BasicMethods::GetRunningTime(comienzo) < iterations); iteration++){
			Assignment assignment;
			assignment.matching = BasicMethods::RandomizedSolution(graph.getDimensions());
			assignment.cost     = BasicMethods::CalculateCost(assignment.matching, graph);
			switch(heuristic){
				case 6:{
					TwoOpt(10, assignment, graph);
				}break;
				case 7:{
					Inversion(10, assignment, graph);
				}break;
				case 8:{
					CircularRotation(10, assignment, graph);
				}break;
				case 9:{
					Kvertex(10, k_value, assignment, graph);
				}break;
				case 10:{
					SelectAllBasic(10, k_value, assignment, graph);
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

	template <class T> std::pair< std::pair<double, double>, std::pair<int, int> > IteratedLocalSearch(
		Assignment &solution, T &graph,
		const int heuristic, const int iterations, const bool seconds, const int k_value){
		std::cout << "Running Iterated Local Search\n";
		clock_t comienzo = clock();
		double avg_time = 0, avg_solution = 0;
		int iteration = 0;
		Assignment assignment = solution;
		for(iteration = 0; (!seconds && iteration < iterations) ||
			(seconds && (int)BasicMethods::GetRunningTime(comienzo) < iterations); iteration++){

			switch(heuristic){
				case 16:{
					TwoOpt(10, assignment, graph);
				}break;
				case 17:{
					Inversion(10, assignment, graph);
				}break;
				case 18:{
					CircularRotation(10, assignment, graph);
				}break;
				case 19:{
					Kvertex(10, k_value, assignment, graph);
				}break;
				case 20:{
					SelectAllBasic(10, k_value, assignment, graph);
				}break;
				default:break;
			}

//			printf("Iter %d %d %d - %d\n", iteration, solution.cost, assignment.cost, iterations);
			if(solution.cost > assignment.cost){
				std::cout <<"Iter = "<<iteration << " before "<<assignment.cost 
					<< " new " << solution.cost<<" time "<<(int)BasicMethods::GetRunningTime(comienzo) <<"\n";
				solution = assignment;
			}

			assignment = solution;
			for(int dim = 1; dim < graph.getNumberOfDimensions(); dim++)
				BasicMethods::Perturbation(assignment, 2, dim);
			assignment.cost = BasicMethods::CalculateCost(assignment.matching, graph);
//			printf("new cost = %d\n", assignment.cost);
		}
		avg_time = BasicMethods::GetRunningTime(comienzo) / (double) iteration;
		avg_solution /= (double) iteration;
		return std::make_pair(std::make_pair(avg_solution, avg_time), std::make_pair(iteration, solution.cost));
	}

	/*
	Created By: Sergio Perez
	Modified by: Sergio Perez
	Created on: February, 2015
	Last Update on:
	*/
	template <class T> bool MultiKSwap(const int K, Assignment &solution, T &graph){
//		std::cout << "Start MultiKSwap\n";
		BasicMethods::validateAssignment(solution.matching);
		int chosenDim = rand()%graph.getNumberOfDimensions(), *indices = new int[K], *vertices = new int[K];

		std::vector<int> p(solution.matching.size(), 0);
		for(unsigned int i = 0; i < solution.matching.size(); i++)
			p[i] = i;

		std::random_shuffle(p.begin(), p.end());

		for(int i = 0; i < std::min(K, (int)solution.matching.size()); i++){
			indices[i] = p[i];
			vertices[i] = solution.matching[indices[i]][chosenDim];
		}

		Assignment assignmentCopy = solution;
		edgeCost previousCost = solution.cost;

//		printf("chosendim %d %d %d %d\n", chosenDim, indices[0], indices[1], indices[2]);

		std::sort(vertices, vertices+K);
		do{

			for(int i = 0; i < K; i++)
				BasicMethods::updateCostGraph(false, indices[i], indices[i], assignmentCopy, graph);

			for(int i = 0; i < K; i++){
				assignmentCopy.matching[indices[i]][chosenDim] = vertices[i];
//				printf("%d%c", vertices[i], i == (K-1) ? '\n' : ' ');
			}

			for(int i = 0; i < K; i++)
				BasicMethods::updateCostGraph(true, indices[i], indices[i], assignmentCopy, graph);

			if(assignmentCopy.cost < solution.cost)
				solution = assignmentCopy;

		}while(std::next_permutation(vertices, vertices+K));

		delete[] indices;
		delete[] vertices;

		BasicMethods::verify(solution, graph, "MultiKSwap");

//		printf("End MultiKSwap %d\n", solution.cost < previousCost);
		return solution.cost < previousCost;
	}//*/


	/*
	Created By: Sergio Perez
	Modified by: Sergio Perez
	Created on: February, 2015
	Last Update on:
	*/

	template <class T> bool MultiKSwapPermutations(const int chosenDim, std::vector<int> &indices, Assignment &solution, T &graph){
//		printf("MultiKSwapPermutations\n");
		int K = (int)indices.size(), *vertices = new int[indices.size()];
		for(int i = 0; i < K; i++){
			vertices[i] = solution.matching[indices[i]][chosenDim];
		}

		Assignment assignmentCopy = solution;
		edgeCost previousCost = solution.cost;
		std::sort(vertices, vertices+K);
		do{

			for(int i = 0; i < K; i++)
				BasicMethods::updateCostGraph(false, indices[i], indices[i], assignmentCopy, graph);

			for(int i = 0; i < K; i++){
				assignmentCopy.matching[indices[i]][chosenDim] = vertices[i];
//				printf("%d%c", vertices[i], i == (K-1) ? '\n' : ' ');
			}

			for(int i = 0; i < K; i++)
				BasicMethods::updateCostGraph(true,  indices[i], indices[i], assignmentCopy, graph);

			if(assignmentCopy.cost < solution.cost){
				printf("mejoro %d %d\n", assignmentCopy.cost, solution.cost);
				solution = assignmentCopy;
			}

		}while(std::next_permutation(vertices, vertices+K));

		delete[] vertices;
		return solution.cost < previousCost;
	}

	/*
	Created By: Sergio Perez
	Modified by: Sergio Perez
	Created on: February, 2015
	Last Update on:
	*/
	template <class T> bool MultiKSwapCombinations(const int N, const int K, const int s, std::vector<int> &p, Assignment &solution, T &graph){
		if(K == 0){
			return MultiKSwapPermutations(s, p, solution, graph);
		}
		p[K-1] = N-1;
		bool answer = MultiKSwapCombinations(N-1, K-1, s, p, solution, graph);
		if(N > K){
			answer = answer | MultiKSwapCombinations(N-1, K, s, p, solution, graph);
		}
		return answer;
	}

		/*
		Created By: Sergio Perez
		Modified by: Sergio Perez
		Created on: February, 2015
		Last Update on:
	*/
	template <class T> bool MultiKSwapExhaustive(const int K, Assignment &solution, T &graph){
		std::cout << "Start MultiKSwapExhaustive "<<K<<"\n";
		BasicMethods::validateAssignment(solution.matching);
		bool answer = false;
		int n = (int)solution.matching.size();
		std::vector<int> p(std::min(K, n), 0);
//		Output::PrintSolution(solution, graph);
		std::sort(solution.matching.begin(), solution.matching.end());
		for(int s = 1; s < (int)graph.getNumberOfDimensions(); s++){
//			std::cout << "Start MultiKSwapCombinations\n";
			answer = answer | MultiKSwapCombinations((int)solution.matching.size(), std::min(K, n), s, p, solution, graph);
		}

		BasicMethods::verify(solution, graph, "MultiKSwapExhaustive");
		return answer;
	}//*/

	template <class T> void SelectDimensionWiseVariation(const int dvh_type,
		Assignment &assignment, T &graph){
			switch(dvh_type){
				case 101:{
					LocalMethodsGutin::DimensionwiseVariation(assignment, graph);
				}break;
				case 102:{
					LocalMethodsGutin::MultiDimensionwiseVariation(assignment, graph);
				}break;
				case 103:{
					LocalMethodsGutin::CompleteDimensionwiseVariation(assignment, graph);
				}break;
				case 104:{
					LocalMethodsGutin::CompleteMultiDimensionwiseVariation(assignment, graph);
				}break;
				case 105:{
					LocalMethodsGutin::MDVH2_3OPT(assignment, graph);
				}break;
				case 106:{
					LocalMethods3D::DimensionwiseVariation3D(assignment, graph);
				}break;
				case 107:{
					LocalMethods3D::MultiDimensionwiseVariation3D(assignment, graph);
				}break;//*/
				case 108:{
					LocalMethods3D::CompleteDimensionwiseVariation3D(assignment, graph);
				}break;
				case 109:{
					LocalMethods3D::CompleteMultiDimensionwiseVariation3D(assignment, graph);
				}break;//*/
				case 111:{
					int apAlgorithm = 1;
					LocalMethodsGutin::DimensionwiseVariation(assignment, graph, apAlgorithm);
				}break;//*/
				case 112:{
					int apAlgorithm = 2;
					LocalMethodsGutin::DimensionwiseVariation(assignment, graph, apAlgorithm);
				}break;//*/
				case 113:{
					int apAlgorithm = 3;
					LocalMethodsGutin::DimensionwiseVariation(assignment, graph, apAlgorithm);
				}break;//*/
				default:break;
			}
		return;
	}

	template <class T> int LocalSearch(const char* inputNameFile,
		const int localSearchType, const int iterations, const bool seconds,
		const int K_value, Assignment &assignment, T &graph){
		clock_t comienzo = clock();
		std::cout << "Start LocalSearch method "<< localSearchType << " iterations " << iterations << " K_value " << K_value << "\n";
		bool success = false;
		int successIterations = 0;
		std::pair<std::pair<double, double>, std::pair<int, int> > summary;
		Assignment mejor = assignment;
		std::cout << "Initial cost = " << assignment.cost << "\n";
		int widthIter = BasicMethods::getWidth(iterations), widthCost = BasicMethods::getWidth(assignment.cost);
		FILE *salida_iter;
		salida_iter = fopen("iterations_random.txt", "a");
		if(localSearchType <= 50){
			for(int count = 1; (!seconds && count <= iterations)
				|| (seconds && (int)BasicMethods::GetRunningTime(comienzo) < iterations) ; count++){
				switch(localSearchType){
					case 1:{
						TwoOpt(iterations, assignment, graph);
						count = iterations;
					}break;
					case 2:{
						Inversion(iterations, assignment, graph);
						count = iterations;
					}break;
					case 3:{
						CircularRotation(iterations, assignment, graph);
						count = iterations;
					}break;
					case 4:{
						Kvertex(iterations, K_value, assignment, graph);
						count = iterations;
					}break;
					case 5:{
						SelectAllBasic(iterations, K_value, assignment, graph);
						count = iterations;
					}break;
					case 6: case 7: case 8: case 9: case 10:{
						LocalMethods::MultiStartBasicLSH(assignment, graph, localSearchType, iterations, seconds, K_value);
						count = iterations;
					}break;
					case 16: case 17: case 18: case 19:case 20:{
						LocalMethods::IteratedLocalSearch(assignment, graph, localSearchType, iterations, seconds, K_value);
						count = iterations;
					}break;
					case 21:{
						success = MultiKSwap(rand()%5+2, assignment, graph);
						successIterations += success;
					}break;
					case 22:{
						success = false;
						for(int z = 2; z < 5; z++)
							success |= MultiKSwapExhaustive(z, assignment, graph);
						successIterations += success;
						if(!success)
							count = iterations;
					}break;//*/
					case 31:{
						success = LocalMethods3D::DimensionwiseVariation3D(assignment, graph);
						successIterations += success;
						if(!success)
							count = iterations;
					}break;
					case 32:{
						success = LocalMethods3D::MultiDimensionwiseVariation3D(assignment, graph);
						successIterations += success;
						if(!success)
							count = iterations;
					}break;
					case 37:{
						success = LocalMethods3D::VarOptGurobi(K_value, assignment, graph);
						successIterations += success;
						if(!success)
							count = iterations;
					}break;
					case 41:{
						success = LocalMethods3D::CompleteDimensionwiseVariation3D(assignment, graph);
					}break;
					case 42:{
						success = LocalMethods3D::CompleteMultiDimensionwiseVariation3D(assignment, graph);
					}break;
					case 11: case 12: case 46: case 47: case 48: case 49: case 50:{
						summary = LocalMethods3D::MultiStartGeneric(assignment, graph,
							localSearchType, iterations, seconds, K_value);
						count = iterations;//forcing to end;
					}break;//*/
					default:{
						printf("Ningun metodo empata con el criterio\n");
						count = iterations;//forcing to end;
					}break;
				}
				if(success){
					printf("(%*d %*d %*d %lf)\n",
						widthIter, count, widthCost, assignment.cost, widthCost, mejor.cost,  BasicMethods::GetRunningTime(comienzo));
				}//*/
				if(assignment.cost != BasicMethods::CalculateCost(assignment.matching, graph)){
					std::cout << "Hubo error en el metodo local "<<localSearchType << "\n";
					exit(0);
				}
				if(assignment.cost < mejor.cost)
					mejor = assignment;
			}
		}
		if((localSearchType >= 51 && localSearchType <= 70) || (80 <= localSearchType && localSearchType < 100)){
			summary = LocalMethodsGutin::LocalMethodsGutin(localSearchType, iterations,
				seconds, K_value, assignment, graph);
		}//*/
		if(70 < localSearchType && localSearchType < 80){
			summary = ExactMethods::ExactMethods(localSearchType, assignment, graph);
		}
		std::cout << "\nFinal cost = " << assignment.cost << "\n";
		if(assignment.cost > mejor.cost)
			assignment = mejor;
		std::cout << "Best final cost "<< assignment.cost << "\n";
		std::sort(assignment.matching.begin(), assignment.matching.end());
//		Output::PrintSolution(assignment);
		std::cout << "Number of success iterations " << successIterations<< "\n";
		std::cout << "End LocalSearch\n";
		std::cout << "The running time was "<<BasicMethods::GetRunningTime(comienzo)<<" seconds.\n";
		std::fprintf(salida_iter, "instance=%s;dim=%d;n=%d;heuristic=%d;iterations=%d;best_cost=%d;avg_cost=%0.1lf;time=%0.1lf\n",
			inputNameFile,	graph.getNumberOfDimensions(), graph.getDimension(0), localSearchType,
			summary.second.first, summary.second.second, summary.first.first, summary.first.second);
		fclose(salida_iter);
		return 1;
	}
};

#endif
