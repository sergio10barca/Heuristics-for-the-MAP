#ifndef CROSSOVERS
#define CROSSOVERS
#include <vector>
#include <set>
#include <algorithm>
#include "LocalMethods.h"
#include "Output.h"

namespace Crossovers{

	bool checker(std::vector<dTuple> &individual){
		for(unsigned int i = 0; i < individual[0].size(); i++){
			std::set<int> q;
			for(unsigned int j = 0; j < individual.size(); j++){
				if(q.find(individual[j][i]) != q.end())
					return false;
				q.insert(individual[j][i]);
			}
		}
		return true;
	}

	template <class T> std::vector<Assignment > PartiallyMappedCrossover(Assignment &father1, Assignment &father2, T &G){
//		std::cout << "Start PartiallyMappedCrossover\n";
		std::vector<Assignment> offspring(2);
		offspring[0] = father1, offspring[1] = father2;

		int point1 = rand()%father1.matching.size(), point2 = rand()%father1.matching.size();
		if(point1 > point2)
			std::swap(point1, point2);
//		std::cout << "Crossover points " << point1 << " " << point2 << "\n";

		for(int i = point1; i <= point2; i++){
			std::swap(offspring[0].matching[i], offspring[1].matching[i]);
		}

		for(int noff = 0; noff < 2; noff++){//for each son
			for(int dim = 1; dim < G.getNumberOfDimensions(); dim++){
//				std::cout << "noff " << noff<< " dim " << dim << "\n";
				std::map<int, int> already;
				std::map<int, int>::iterator it;
				for(int i = point1; i <= point2; i++){
					already[offspring[noff].matching[i][dim]] = i;
				}

				for(int i = 0; i < offspring[noff].matching.size(); i++){
					if(i < point1 || i > point2){
						while((it = already.find(offspring[noff].matching[i][dim])) != already.end()){
							offspring[noff].matching[i][dim] = offspring[(noff+1)%2].matching[it->second][dim];
						}
					}
				}
			}
		}

		for(int i = 0; i < 2; i++){
			if(!checker(offspring[i].matching)){
				std::cout << "Error at offspring"<<i << "\n";
				exit(0);
			}
			offspring[i].cost = BasicMethods::CalculateCost(offspring[i].matching, G);
		}

/*		Output::PrintSolution(father1, G);
		Output::PrintSolution(father2, G);
		Output::PrintSolution(offspring[0], G);
		Output::PrintSolution(offspring[1], G);
		std::cout << "End PartiallyMappedCrossover\n";//*/
		return offspring;
	}

	template <class T> std::vector<Assignment > CycledCrossover(Assignment &father1, Assignment &father2, T &G){
//		std::cout << "Start CycledCrossover\n";
		std::vector<Assignment> offspring(2);

		offspring[0] = father1, offspring[1] = father2;

		const int point = rand()%father1.matching.size();
//		std::cout << "Cycled crossing point " << point << "\n";

		std::vector<int> indexes(offspring[0].matching.size());

		for(int dim = 1; dim < G.getNumberOfDimensions(); dim++){
			for(int index = 0; index < offspring[0].matching.size(); index++){
				indexes[ offspring[0].matching[index][dim] ] = index;
			}

			int index, start_vertex = offspring[0].matching[point][dim];
			for(index = point; start_vertex != offspring[1].matching[index][dim]; ){
				int next_point = offspring[1].matching[index][dim];
				std::swap(offspring[0].matching[index][dim], offspring[1].matching[index][dim]);
				index = indexes[ next_point ];
			}
			std::swap(offspring[0].matching[index][dim], offspring[1].matching[index][dim]);
		}

		for(int i = 0; i < 2; i++){
			if(!checker(offspring[i].matching)){
				std::cout << "Error at offspring"<<i << "\n";
				exit(0);
			}
			offspring[i].cost = BasicMethods::CalculateCost(offspring[i].matching, G);
		}

/*		Output::PrintSolution(father1, G);
		Output::PrintSolution(father2, G);
		Output::PrintSolution(offspring[0], G);
		Output::PrintSolution(offspring[1], G);
		std::cout << "End CycledCrossover\n";//*/
		return offspring;
	}

	template <class T> std::vector<Assignment > OrderCrossover(Assignment &father1, Assignment &father2, T &G){
//		std::cout << "Start OrderCrossover\n";
		std::vector<Assignment> offspring(2), fathers(2);
		fathers[0] = father1, fathers[1] = father2;
		offspring[0] = father1, offspring[1] = father2;

		int point[2];
		for(int i = 0; i < 2; i++)
			point[i] = rand()%father1.matching.size();
		std::sort(point, point+2);
//		std::cout << point[0]<<" "<<point[1] << "\n";

		std::vector<bool> flag(offspring[0].matching.size());
		std::vector<int>  vertices_order(offspring[0].matching.size());

		for(int noff = 0; noff < 2; noff++){//for each son
			for(int dim = 1; dim < G.getNumberOfDimensions(); dim++){
				int next = (noff+1)&1;
				for(int i = 0; i < flag.size(); i++){
					flag[offspring[noff].matching[i][dim]] = point[0] <= i && i <= point[1];//set true elements in point[0]...point[1]
					vertices_order[i] = fathers[next].matching[i][dim];
				}	
				int next_index_start = 0;
				for(int index = 0; index < offspring[noff].matching.size(); index++){
					if(point[0] <= index && index <= point[1]){
						continue;
					}
					while(next_index_start < flag.size() && flag[vertices_order[next_index_start]] == true){
						next_index_start++;
					}
					offspring[noff].matching[index][dim] = vertices_order[next_index_start];
					next_index_start++;
				}
			}
		}

		for(int i = 0; i < 2; i++){
			if(!checker(offspring[i].matching)){
				std::cout << "Error at offspring"<<i << "\n";
				exit(0);
			}
			offspring[i].cost = BasicMethods::CalculateCost(offspring[i].matching, G);
		}

/*		Output::PrintSolution(father1, G);
		Output::PrintSolution(father2, G);
		Output::PrintSolution(offspring[0], G);
		Output::PrintSolution(offspring[1], G);
		std::cout << "End OrderCrossover\n";//*/
		return offspring;
	}

	template <class T> std::vector<Assignment > ClonedCrossover(Assignment &father1, Assignment &father2, T &G){
		std::cout << "Start ClonedCrossover\n";
		std::vector<Assignment> offspring(2), fathers(2);
		fathers[0] = father1, fathers[1] = father2;
		offspring[0] = father1, offspring[1] = father2;

		std::vector<int> flag(offspring[0].matching.size());
		std::vector<int> different_indexes;
		for(int i = 0; i < offspring[0].matching.size(); i++){
			if(offspring[0].matching[i] != offspring[1].matching[i]){
				different_indexes.push_back(i);
				flag[i] = false;
			}else{
				flag[i] = true;
//				std::cout << "i = "<<i<<" ";
			}
		}
//		std::cout << "\n";
		for(int noff = 0; noff < 2; noff++){//for each son
			for(int dim = 0; dim < G.getNumberOfDimensions(); dim++){
				int next = (noff+1)&1;
				std::next_permutation(different_indexes.begin(), different_indexes.end());		
				for(int i = 0, index_perm = 0; i < offspring[noff].matching.size(); i++){
					if(flag[i] == false){
						offspring[noff].matching[i][dim] = fathers[next].matching[ different_indexes[index_perm] ][dim];
						++index_perm;
					}
				}
			}
		}

		for(int i = 0; i < 2; i++){
			if(!checker(offspring[i].matching)){
				std::cout << "Error at offspring"<<i << "\n";
				exit(0);
			}
			offspring[i].cost = BasicMethods::CalculateCost(offspring[i].matching, G);
		}

/*		Output::PrintSolution(father1, G);
		Output::PrintSolution(father2, G);
		Output::PrintSolution(offspring[0], G);
		Output::PrintSolution(offspring[1], G);
		std::cout << "End ClonedCrossover\n";//*/
		return offspring;
	}

	std::vector<Assignment > FatherSelection(std::vector<Assignment> &individual){
//		cout << "Start FatherSelection for Crossings " << endl;
		std::vector<Assignment > fathers;
		std::set<int> q;
		unsigned int cont;
		for(cont = 0; cont < 2; cont++){
			int index;
			while(q.find(index=rand()%individual.size()) != q.end());
			q.insert(index);
			fathers.push_back(individual[index]);
		}
		q.clear();
//		cout << "End FatherSelection " << endl;
		return fathers;
	}

	template <class T> std::vector<Assignment> SelectedCrossover(const int typeCrossover, 
		double proportionOffspring, const int lsh_type, std::vector<Assignment> &individual, T &graph){
//		std::cout << "Start SelectedCrossover\n";
		std::vector<Assignment > offspring, sons;
//		cout << proportionOffspring << endl;
		int numberOffsprings = (int)floor((individual.size()*proportionOffspring)/100.0);
//		std::cout << "Sons " <<numberOffsprings<< "\n";
		for(int cont = 0; (int)offspring.size() < numberOffsprings && cont < (numberOffsprings*10); cont++){
//			std::cout << "iterando crossover "<<cont << "\n";
			std::vector<Assignment > fathers = FatherSelection(individual);
			switch(typeCrossover){
				case 1:{//
					sons = PartiallyMappedCrossover(fathers[0], fathers[1], graph);
				}break;
				case 2:{//
					sons = CycledCrossover(fathers[0], fathers[1], graph);
				}break;
				case 3:{//
					sons = OrderCrossover(fathers[0], fathers[1], graph);
				}break;
				case 4:{//
					sons = ClonedCrossover(fathers[0], fathers[1], graph);
				}break;
				default:break;
			}
/*			std::cout << "Costo padre1 " << BasicMethods::CalculateCost(fathers[0].matching, graph) << "\n";
			std::cout << "Costo padre2 " << BasicMethods::CalculateCost(fathers[1].matching, graph) << "\n";
			exit(0);//*/
			for(unsigned int j = 0; j < sons.size() && (int)offspring.size() < numberOffsprings; j++){
				sons[j].cost = BasicMethods::CalculateCost(sons[j].matching, graph);
//				cout << "Costo hijo "<< j<<" " << sons[j].cost << " ";				
				LocalMethods::SelectDimensionWiseVariation(lsh_type, sons[j], graph);
				offspring.push_back(sons[j]);//*/
			}
		}
//		std::cout << "End SelectedCrossover "<<offspring.size() << "\n";
		return offspring;
	}
}

#endif