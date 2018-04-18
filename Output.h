#ifndef OUTPUT
#define OUTPUT
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include "Structures.h"

namespace Output{

	/*
		This function prints the current matching in time O(r)
	*/
	void PrintMatching(std::vector<std::pair<std::pair<int, int>, int> > matching){
		printf("Current Matching\nwoman->man\n");
		for(unsigned int i = 0; i < matching.size(); i++){
			std::cout << matching[i].first.first << " "<<matching[i].first.second<<" "<<matching[i].second << "\n";
		}
		return;
	}

	/*
		Created By: Sergio Perez
		Modified by:
		Revised By:
		Created on: March, 2013
		Last Update on:
		This function writes the matching of minimum cost in a file in adjacency list format
		Input:
			**output: FILE; the reference of the file where will be printed the matching
	*/
	void Adjacency_list_output(FILE **output, std::vector<std::pair<std::pair<int, int>, int> > matching){
		long long int matching_cost = 0;
		for(unsigned int index = 0; index < matching.size(); index++)
			matching_cost += matching[index].second;
		std::cout << "Minimum cost = " << matching_cost << " Matching size = "<<matching.size()<<"\n";
		fprintf(*output, "%lld %u\n", matching_cost, matching.size());
		for(unsigned int index = 0; index < matching.size(); index++)
			fprintf(*output, "%d %d %d\n", matching[index].first.first, matching[index].first.second, matching[index].second);
		return;
	}

	/*
		Created By: Sergio Perez
		Modified by:
		Revised By:
		Created on: March, 2013
		Last Update on:
		This function prints the minimum cost matching in the formats indicated by the user from the input parameters
		Input:
			**output: FILE; the reference of the file where will be printed the matching
	*/
	void Print_output(const int in_type, int out_type, char *output_name, std::vector<std::pair<std::pair<int, int>, int> > matching){
		std::string s[3] = {"_ady.out", "_tex.out", "_bin.out"};
		std::string mode[3] = {"w", "w", "wb"};
		char file_name[1000];
		int i, j, cont = 1;
		for(cont = 0; out_type > 0 && cont < 3; cont++, out_type>>=1){
			if(out_type&1){
				for(i = 0; (file_name[i] = output_name[i]) != '\0';i++);
				for(j = 0; (file_name[i] = s[cont].c_str()[j]) != '\0'; i++, j++);
	//			strcat(file_name, s[cont].c_str());
				FILE *output = fopen(file_name, mode[cont].c_str());
				if(output == NULL){
					std::cout << "The specified output file does not exist.";
					exit(30);
				}
				switch(cont){
					case 0:{
						Adjacency_list_output(&output, matching);
					}break;
					default:break;
				}
				fclose(output);
			}
		}
		return;
	}

	void PrintDTuple(dTuple u){
		for(int i = 0; i < u.size(); i++)
			std::cout << u[i] << char((i+1) == u.size() ? '\n' : ' ');
		return;
	}

	void PrintSolution(std::vector<dTuple> solution){
		std::cout << "Assignment = \n";
		for(unsigned int i = 0; i < solution.size(); i++)
			for(unsigned int j = 0; j < solution[i].size(); j++)
				std::cout << solution[i][j] << char((j+1) == solution[i].size() ? '\n' : ' ');
		return;
	}

	template <class T> void PrintSolution(Assignment &solution, T &G){
		std::vector<dTuple> &assignment = solution.matching;
		edgeCost totalCost = 0;
		std::cout << "Current solution\n";
		std::cout << assignment.size() << " " << solution.cost << "\n";
		for(unsigned int i = 0; i < assignment.size(); i++){
			for(unsigned int dimA = 0; dimA < assignment[i].size(); dimA++)
				std::cout << assignment[i][dimA] << " ";
			edgeCost cost = G.getEdgeCost(assignment[i]);
			std::cout << cost << "\n";
			totalCost += cost;
		}
/*		if(totalCost != solution.cost){
			std::cout << "Error by calculating assignment cost "<<totalCost << " "<<solution.cost << "\n";
			exit(0);
		}//*/
		return;
	}//*/

	template <class T> void PrintSolution(std::ofstream &output, Assignment &solution, T &G){
		std::vector<dTuple> &assignment = solution.matching;
		edgeCost totalCost = 0;
		output << assignment.size() << " " << solution.cost << "\n";
		for(unsigned int i = 0; i < assignment.size(); i++){
			for(unsigned int dimA = 0; dimA < assignment[i].size(); dimA++)
				output << assignment[i][dimA] << " ";
			edgeCost cost = G.getEdgeCost(assignment[i]);
			output << cost << "\n";
			totalCost += cost;
		}
		if(totalCost != solution.cost){
			std::cout << "Error by calculating assignment cost "<<totalCost << " "<<solution.cost << "\n";
			output.close();
			exit(0);
		}//*/
		return;
	}

}

#endif
