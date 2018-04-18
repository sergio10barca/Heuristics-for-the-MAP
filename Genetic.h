#ifndef GENETIC
#define GENETIC
#include <vector>
#include "BasicMethods.h"
#include "LocalMethods.h"
#include "Mutations.h"
#include "Crossovers.h"
#include "Selections.h"

class Genetic{
	private:
	int typeLSH, typeCrossover, typeMutation, typeSelection;
	int generations, populationSize;
	double probMutation, proportionOffspring;

	public:

	void PrintPopulation(std::vector<Assignment > &individual){
		std::cout << "Start PrintPopulation\n";
		for(unsigned int i = 0; i < individual.size(); i++){
			std::cout << "i = "<<i << " cost = "<<individual[i].cost << "\n";
		}
		std::cout << "End PrintPopulation\n";
		return;
	}

	Genetic(const int typeHeuristic, const int generations, const int populationSize,
			const double probMutation, const double proportionOffspring,
			const int typeCrossover, const int typeMutation, const int typeSelection){

		std::cout << "Start Genetic constructor\n";
		this->typeLSH        = typeHeuristic;
		this->generations    = generations;
		this->populationSize = populationSize;
		this->probMutation   = probMutation;
		this->proportionOffspring = proportionOffspring;
		this->typeCrossover  = typeCrossover;
		this->typeMutation   = typeMutation;
		this->typeSelection  = typeSelection;
		std::cout << "End Genetic constructor\n";
	}

	double fittness(const std::vector<Assignment> &individual){
		double sum = 0;
		for(int i = 0; i < individual.size(); i++){
			sum += individual[i].cost;
		}
		return sum / (double)individual.size();
	}

	template <class T> Assignment executeBasicGenetic(T &graph){
		std::cout << "Start BasicGenetic\n";
		std::vector<Assignment> individual;
		Assignment bestInd;
		for(int i = 0; i < populationSize; i++){
//			std::cout << "individual "<< i <<" "<<populationSize<< "\n";
			Assignment assignment(BasicMethods::RandomizedSolution(graph.getDimensions()), 0);
			assignment.cost = BasicMethods::CalculateCost(assignment.matching, graph);

			LocalMethods::SelectDimensionWiseVariation(this->typeLSH, assignment, graph);

			individual.push_back(assignment);
		}
		std::sort(individual.begin(), individual.end());
		bestInd = individual[0];
//		this->PrintPopulation(individual);//*/
		populationSize = (int)individual.size();
		std::cout << "the populationSize is " << populationSize << "\n";

		clock_t comienzo = clock();
		for(int iteration = 0; (int)BasicMethods::GetRunningTime(comienzo) < generations; iteration++){
//		for(int iteration = 0; iteration < generations; iteration++){
			//crossover
			std::vector<Assignment > offspring =
				Crossovers::SelectedCrossover(this->typeCrossover, this->proportionOffspring,
												this->typeLSH, individual, graph);

			//REALIZO LA SELECCION
//			cout << "Selecting "<<individual.size() << " "<<offspring.size() << endl;
			individual = Selections::Selection(individual, (int)individual.size()-(int)offspring.size(), this->typeSelection);

//			cout << "Selecting offspring from " << individual.size()<<" "<<offspring.size() << endl;
			for(int i = (int)individual.size(), j = 0; i < this->populationSize && j < (int)offspring.size(); i++, j++)
				individual.push_back(offspring[j]);

			//REALIZO LA MUTACION
//			this->PrintPopulation(individual);
			Mutations::SelectedMutation(this->typeMutation, this->probMutation, this->typeLSH, individual, graph);

			std::sort(individual.begin(), individual.end());
//			this->PrintPopulation(individual);//*/
			for(int i = 0; i < (int)individual.size(); i++){
				if(individual[i].cost != BasicMethods::CalculateCost(individual[i].matching, graph)){
					std::cout << "Algun costo fue mal calculado\n";
					exit(0);
				}
			}
			if(individual[0].cost < bestInd.cost){
				bestInd = individual[0];
			}
			if(iteration % 20 == 0){
              std::cout << "Iteration " << iteration <<" ";
			  std::printf("%3.2f %d\n", fittness(individual), bestInd.cost);
			}
		}

/*		for(unsigned int i = 0; i < individual.size(); i++)
			Output::PrintSolution(individual[i].matching);//*/
		std::cout << "End BasicGenetic\n";
		std::cout << "best cost = "<<bestInd.cost << "\n";
		BasicMethods::verify(bestInd, graph, "Genetic");
		return bestInd;
	}

	template <class T> std::vector<Assignment> generateCandidatesPool(const int lsh_type, std::vector<Assignment> individual, T &graph){
		std::vector<Assignment > candidates;
		for(int i = 0; i < 200; i++){
			int a = rand()%individual.size(), b = rand()%individual.size();
			std::vector<Assignment > sons = Crossovers::PartiallyMappedCrossover(individual[a], individual[b], graph);
			for(int j = 0; j < sons.size(); j++){
				LocalMethods::SelectDimensionWiseVariation(lsh_type, sons[j], graph);
				candidates.push_back(sons[j]);
			}
		}
		return candidates;
	}

	template <class T> Assignment hybridGeneticAlgorithm(T &graph){
		std::cout << "Start hybridGeneticAlgorithm\n";
		std::vector<Assignment> individual;
		Assignment bestInd;
		for(int i = 0; i < 100; i++){
//			std::cout << "individual "<< i <<" "<<populationSize<< "\n";
			Assignment assignment(BasicMethods::RandomizedSolution(graph.getDimensions()), 0);
			assignment.cost = BasicMethods::CalculateCost(assignment.matching, graph);

			LocalMethods::SelectDimensionWiseVariation(this->typeLSH, assignment, graph);

			individual.push_back(assignment);
		}
		sort(individual.begin(), individual.end());
		bestInd = individual[0];
		clock_t comienzo = clock();
		int no_success = 0, iteration = 0;

		while(true){
			std::vector<Assignment> candidates = generateCandidatesPool(this->typeLSH, individual, graph);

			sort(candidates.begin(), candidates.end());

			int cont = 0;
			individual[cont++] = candidates[0];
			for(int i = 1; i < candidates.size() && cont < 100; i++){
				if(!(individual[cont-1] == candidates[i])){
					individual[cont++] = candidates[i];
				}
			}

			if(individual[0].cost < bestInd.cost){
				bestInd = individual[0];
				no_success = 0;
			}else{
				no_success++;
			}
			std::printf("%d - %d ", no_success, cont);
			if(no_success > 10 || cont < 100){
				break;
			}
			iteration++;
            std::cout << "Iteration " << iteration <<" ";
			std::printf("%3.2f %d\n", fittness(individual), bestInd.cost);
		}

/*		for(unsigned int i = 0; i < individual.size(); i++)
			Output::PrintSolution(individual[i].matching);//*/
		std::cout << "End hybridGeneticAlgorithm\n";
		std::cout << "best cost = "<<bestInd.cost << "\n";
		BasicMethods::verify(bestInd, graph, "Genetic");
		return bestInd;
	}


};

#endif
