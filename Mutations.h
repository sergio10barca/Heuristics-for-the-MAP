#ifndef MUTATIONS
#define MUTATIONS
#include <vector>
#include <iostream>
#include <algorithm>
#include "Structures.h"
#include "BasicMethods.h"
#include "LocalMethods.h"

namespace Mutations{
	/*This functions takes O((s - 1) * n) time*/
	template <class T> void RandomSwapping(Assignment &solution, T &graph){
//		std::cout << "Start Mutation by RandomSwapping\n";
		BasicMethods::validateAssignment(solution.matching);

		/**************Performing some random swaps at each dimension *****************************/
		/*The dimension 0 should no be changed*/
		for(int dim = 1; dim < solution.matching[0].size(); dim++){
			int indexU = rand()%solution.matching.size(), indexV = rand()%solution.matching.size();

			BasicMethods::updateCostGraph(false, indexU, indexU, solution, graph);
			BasicMethods::updateCostGraph(false, indexV, indexV, solution, graph);

			std::swap(solution.matching[indexU][dim], solution.matching[indexV][dim]);

			BasicMethods::updateCostGraph(true, indexU, indexU, solution, graph);
			BasicMethods::updateCostGraph(true, indexV, indexV, solution, graph);
		}
		BasicMethods::verify(solution, graph, "Mutation by RandomSwapping");
		return;
	}

	/*This functions takes O((s - 1) * n) time*/
	template <class T> void RandomInversion(Assignment &solution, T &graph){
//		std::cout << "Start Mutation by RandomInversion\n";
		BasicMethods::validateAssignment(solution.matching);

		/*The dimension 0 should no be changed*/
		for(int dim = 1; dim < solution.matching[0].size(); dim++){
			int index[2] = {rand()%(int)solution.matching.size(), rand()%(int)solution.matching.size()};
			std::sort(index, index+2);

			BasicMethods::updateCostGraph(false, index[0], index[1], solution, graph);

			for(int i = index[0], j = index[1]; i < j; i++, j--){
				std::swap(solution.matching[index[0]][dim], solution.matching[index[1]][dim]);
			}
			BasicMethods::updateCostGraph(true, index[0], index[1], solution, graph);
		}
		BasicMethods::verify(solution, graph, "Mutation by RandomInversion");
		return;
	}

	/*This functions takes O((s - 1) * n) time*/
	template <class T> void RandomPermutation(Assignment &solution, T &graph){
//		std::cout << "Start Mutation by RandomPermutation\n";
		BasicMethods::validateAssignment(solution.matching);

		/*The dimension 0 should no be changed*/
		for(int dim = 1; dim < solution.matching[0].size(); dim++){
			int index[2] = {rand()%(int)solution.matching.size(), rand()%(int)solution.matching.size()};
			std::sort(index, index+2);
			std::vector<int> permutation(index[1]-index[0]+1);

			BasicMethods::updateCostGraph(false, index[0], index[1], solution, graph);

			for(int i = index[0]; i <= index[1]; i++){
				permutation[ i - index[0] ] = solution.matching[ i ][ dim ];
			}
			std::next_permutation(permutation.begin(), permutation.end());
			for(int i = index[0]; i <= index[1]; i++){
				solution.matching[ i ][ dim ] = permutation[ i - index[0] ];
			}

			BasicMethods::updateCostGraph(true, index[0], index[1], solution, graph);
		}
		BasicMethods::verify(solution, graph, "Mutation by RandomPermutation");
		return;
	}

	template <class T> void SelectedMutation(const int typeMutation,
		const double probMutation, const int lsh_type, std::vector<Assignment> &individual, T &graph){
//		std::cout << "Start SelectedMutation\n";
		for(unsigned int i = 0; i < individual.size(); i++){
			if(rand()%10000 < ((int)(probMutation*100))){
				Assignment aux = individual[i];
//				std::cout << "antes "<<individual[i].cost << "\n";
				switch(typeMutation){
					case 1:{
						Mutations::RandomSwapping(individual[i], graph);
					}break;
					case 2:{
						Mutations::RandomInversion(individual[i], graph);
					}break;
					case 3:{
						Mutations::RandomPermutation(individual[i], graph);
					}break;
					default:break;
				}
				LocalMethods::SelectDimensionWiseVariation(lsh_type, individual[i], graph);
//				std::cout << "nuevo despues "<<individual[i].cost << "\n";
			}
		}
//		std::cout << "End SelectedMutation\n";
		return;
	}

};
#endif
