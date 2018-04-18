#ifndef SELECTIONS
#define SELECTIONS
#include <vector>
#include <iostream>
#include "BasicMethods.h"
	
namespace Selections{

	std::vector<Assignment > Elitist(const std::vector<Assignment> &individual, const int resultSize){
//		std::cout << "Start Selection elitist\n";
		std::vector<Assignment> result;
		for(int i = 0; i < resultSize && i < (int)individual.size(); i++)
			result.push_back(individual[i]);
		return result;
	}

	std::vector<Assignment> Tournament(const std::vector<Assignment> &individual, const int resultSize){
//		std::cout << "Start Selection Tournament\n";
		std::vector<Assignment> result;
		std::vector<int> perm(individual.size(), 0);
		for(int i = 0; i < individual.size(); i++)
			perm[i] = i;

		std::random_shuffle(perm.begin(), perm.end());
//		cout << " perm "<<perm.size() << endl;
		for(int i = 0; (int)result.size() < resultSize; i = (i + 2)%individual.size()){
//			cout << "iteracion " << (i/2) << " torneo" << endl;
			int indA = perm[i], indB = perm[(i+1)%individual.size()];
//			cout << indA<<" "<<individual[indA].cost <<" "<<indB<<" "<< individual[indB].cost << endl;			
			if(individual[indA].cost < individual[indB].cost)
				result.push_back(individual[indA]);
			else
				result.push_back(individual[indB]);
		}
	//	std::cout << "End Selection Tournament\n";
		return result;
	}

	std::vector<Assignment> RouletteWheel(const std::vector<Assignment> &individual, const int resultSize){
//		std::cout << "Start Selection RouletteWheel\n";
		std::vector<int> dp(individual.size());
		std::vector<Assignment> result;
		int fitnessSum = 0;
		int maxi = individual[individual.size()-1].cost+1;
		for(int i = 0; i < individual.size(); i++){
//			std::cout << (maxi - individual[i].cost) << "\n";
			fitnessSum += maxi-individual[i].cost;
			dp[i] = fitnessSum;
		}
		for(int i = 0; i < resultSize; i++){
			int alpha = rand()%fitnessSum;
			std::vector<int>::iterator it = lower_bound(dp.begin(), dp.end(), alpha);
			int index = it-dp.begin();
			result.push_back(individual[index]);
		}		
	//	std::cout << "End Selection RouletteWheel\n";
		return result;
	}

	std::vector<Assignment> Selection(const std::vector<Assignment> &individual, int numberOfSurvivors, const int typeSelection){
//		std::cout << "Selection\n";
		switch(typeSelection){
			case 1:{ 
				return Selections::Elitist(individual, numberOfSurvivors);
			}break;
			case 2:{ 
				return Selections::Tournament(individual, numberOfSurvivors);
			}break;
			case 3:{ 
				return Selections::RouletteWheel(individual, numberOfSurvivors);
			}break;
			default:break;
		}
		return std::vector<Assignment> (0);
	}
};
#endif
