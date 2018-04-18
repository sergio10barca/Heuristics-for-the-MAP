#include "Structures.h"
#include "BasicMethods.h"
#include "Input.h"
#include "Output.h"
#include "ExactMethods.h"
#include "LocalMethods.h"
#include "Genetic.h"

void solve(const int inputType, const int typeHeuristic, const int iterations,
	const bool seconds, const int K_value, const char *inputNameFile, const char *outputNameFile,
	const int selection = 0, const int cross = 0, const int mutation = 0, const int population = 0){
	std::cout << "Start solve\n";

	std::ofstream output(outputNameFile, std::ofstream::out);
	if(!output.is_open()){
		std::cout << "Error: The output file cannot be open.\n";
		exit(0);
	}

	Assignment assignment;
	switch(inputType){
		case 1:{
			HyperGraph hyperGraph;
			Input::ReadAsIndependent(inputNameFile, hyperGraph);

			if(1 <= typeHeuristic && typeHeuristic < 100){
				assignment.matching = BasicMethods::RandomizedSolution(hyperGraph.getDimensions());
				assignment.cost = BasicMethods::CalculateCost(assignment.matching, hyperGraph);
				LocalMethods::LocalSearch(inputNameFile, typeHeuristic, iterations, seconds,
					K_value, assignment, hyperGraph);
				Output::PrintSolution(output, assignment, hyperGraph);
			}
/*			if(70 < typeHeuristic && typeHeuristic < 100){
				assignment = ExactMethods::ExactMethods(typeHeuristic, hyperGraph);
			}//*/
			if(100 <= typeHeuristic){
				Genetic gen(typeHeuristic, iterations, population, 20.0, 50.0, cross, mutation, selection);
				if(100 <= typeHeuristic && typeHeuristic < 110){
					assignment = gen.executeBasicGenetic(hyperGraph);
				}else{
					assignment = gen.hybridGeneticAlgorithm(hyperGraph);
				}
				Output::PrintSolution(output, assignment, hyperGraph);
			}
		}break;

		case 2: case 3:{
			std::cout << "Este es el caso\n";
			CompleteHyperGraph completeHyperGraph;
			if(inputType == 2)
				Input::ReadBalasSaltzman(inputNameFile, completeHyperGraph);
			else if(inputType == 3)
				Input::ReadBinaryAsKarapetyanGutin(inputNameFile, completeHyperGraph);

			if(1 <= typeHeuristic && typeHeuristic < 100){
				assignment.matching = BasicMethods::RandomizedSolution(completeHyperGraph.getDimensions());
				assignment.cost = BasicMethods::CalculateCost(assignment.matching, completeHyperGraph);
				LocalMethods::LocalSearch(inputNameFile, typeHeuristic, iterations, seconds,
					K_value, assignment, completeHyperGraph);
			}

			if(100 <= typeHeuristic){
				Genetic gen(typeHeuristic, iterations, population, 10.0, 50.0, cross, mutation, selection);
				if(100 <= typeHeuristic && typeHeuristic < 110){
					assignment = gen.executeBasicGenetic(completeHyperGraph);
				}else{
					assignment = gen.hybridGeneticAlgorithm(completeHyperGraph);
				}
			}

			Output::PrintSolution(output, assignment, completeHyperGraph);
		}break;
		default:
		break;
	}
	return;
}

void use(){
	std::cout << "\nUse:Executable <inputType> <inputFile> <outputFile> <localSearch> <iterations> <seed> <seconds>\n\n";
	std::cout << "<inputType>: One of the next cases:\n";
	std::cout << "\tA. Input as independent costs (Hyper Graph).\n";
	std::cout << "\t\t1. General Multipartite Hyper Graphs.\n";
	std::cout << "\t\t2. Complete Multipartite Hyper Graphs (Integer format by Balas and Saltzman).\n";
	std::cout << "\t\t3. Complete Multipartite Hyper Graphs (Binary format by Karapetyan and Gutin).\n";
	std::cout << "<inputFile>: Name of data input file.\n\n";
	std::cout << "<outputFile>: Name of data output file.\n\n";
	std::cout << "<localSearch>: One of the next cases:\n";
	std::cout << "\t1. Simple 2-Opt.\n";
	std::cout << "\t2. Inversion.\n";
	std::cout << "\t3. CircularRotation.\n";
	std::cout << "\t4. K-vertex contiguous permutated.\n";
	std::cout << "\t5. Select All Basic (1-5).\n";
	std::cout << "\t6. Multi Start Simple 2-Opt.\n";
	std::cout << "\t7. Multi Start Inversion.\n";
	std::cout << "\t8. Multi Start Circular Rotation.\n";
	std::cout << "\t9. Multi Start K-vertex.\n";
	std::cout << "\t10. Multi Start All basic.\n";
	std::cout << "\t11. Multi Start DVH3+2-opt.\n";
	std::cout << "\t12. Multi Start DVH3+3-opt.\n";

	std::cout << "\t16. Iterated Local Search Simple 2-Opt.\n";
	std::cout << "\t17. Iterated Local Search Inversion.\n";
	std::cout << "\t18. Iterated Local Search Circular Rotation.\n";
	std::cout << "\t19. Iterated Local Search K-vertex.\n";
	std::cout << "\t20. Iterated Local Search All basic.\n";

	std::cout << "\t21. Multi k swap.\n";
	std::cout << "\t22. Multi k swap extended.\n";
	std::cout << "\t31. Round Local DimensionWiseVariation3D.\n";
	std::cout << "\t32. Round Local MultiDimensionwiseVariation3D.\n";
	std::cout << "\t37. K-Opt (Gurobi).\n";
	std::cout << "\t41. Complete DimensionWiseVariation3D.\n";
	std::cout << "\t42. Complete MultiDimensionwiseVariation3D.\n";
	std::cout << "\t46. Multi Start DimensionWiseVariation3D.\n";
	std::cout << "\t47. Multi Start MultiDimensionwiseVariation3D.\n";
	std::cout << "\t48. Multi Start 2-Opt (Gurobi) heuristic.\n";
	std::cout << "\t49. Multi Start 3-Opt (Gurobi) heuristic.\n";
	std::cout << "\t50. Multi Start K-Opt (Gurobi) heuristic.\n";

	std::cout << "\t51. DimensionwiseVariation DV (Gutin and Karapetyan).\n";
	std::cout << "\t52. MultiDimensionwiseVariation (Gutin and Karapetyan).\n";
	std::cout << "\t53. TwoOpt (Gutin and Karapetyan).\n";
	std::cout << "\t54. ThreeOpt (Gutin and Karapetyan).\n";
	std::cout << "\t55. K-Opt (Based on Gutin and Karapetyan ideas).\n";
	std::cout << "\t56. VOpt (Gutin and Karapetyan).\n";
	std::cout << "\t57. Multi start SDV2 (SDV2+2-Opt).\n";
	std::cout << "\t58. Multi start DVH2 (DVH2+2-Opt).\n";
	std::cout << "\t59. Multi start DVHVopt (MDV+Vopt).\n";
	std::cout << "\t60. Multi start DVH2 (DVH2+3-Opt).\n";
	std::cout << "\t61. Complete DimensionwiseVariation DV (Gutin and Karapetyan).\n";
	std::cout << "\t62. Complete MultiDimensionwiseVariation (Gutin and Karapetyan).\n";
	std::cout << "\t63. Complete Two Opt (Gutin and Karapetyan).\n";
	std::cout << "\t64. Complete Three Opt (Gutin and Karapetyan).\n";
	std::cout << "\t66. Multi Start DimensionwiseVariation DV (Gutin and Karapetyan).\n";
	std::cout << "\t67. Multi Start MultiDimensionwiseVariation (Gutin and Karapetyan).\n";
	std::cout << "\t68. Multi Start Two Opt (Gutin and Karapetyan).\n";
	std::cout << "\t69. Multi Start Three Opt (Gutin and Karapetyan).\n";
	std::cout << "\t70. Multi Start K Opt (Gutin and Karapetyan).\n";
	std::cout << "\t71. Exact Method By Permutations (until 8 vertices for 3AP).\n";
	std::cout << "\t72. Exact Method By Dynamic Programming (until 14 vertices for 3AP).\n";
	std::cout << "\t73. Exact Method By Permutations and Hungarian (until 10 vertices for 3AP).\n";
	std::cout << "\t74. Exact Method By Gurobi Solver (until 70 vertices for 3AP in 30 seconds).\n";
	std::cout << "\t75. Exact Method By Gurobi Solver (until 20 vertices for 4AP ).\n";
	std::cout << "\t76. Exact Method By Gurobi Solver (for MAP).\n";
	std::cout << "\t101. Basic Genetic + Single DVH2. Extra parameters are required. Selection, Crossover  and Mutation types.\n";
	std::cout << "\t102. Basic Genetic + Single MDVH2. Extra parameters are required. Selection, Crossover  and Mutation types.\n";
	std::cout << "\t103. Basic Genetic + Complete DVH2. Extra parameters are required. Selection, Crossover  and Mutation types.\n";
	std::cout << "\t104. Basic Genetic + Complete MDVH2. Extra parameters are required. Selection, Crossover  and Mutation types.\n";
	std::cout << "\t105. Basic Genetic + Complete MDVH2+3OPT. Extra parameters are required. Selection, Crossover  and Mutation types.\n";
	std::cout << "\t106. Basic Genetic + Single DVH3. Extra parameters are required. Selection, Crossover  and Mutation types.\n";
	std::cout << "\t107. Basic Genetic + Single MDVH3. Extra parameters are required. Selection, Crossover  and Mutation types.\n";
	std::cout << "\t108. Basic Genetic + Complete MDVH2+3OPT. Extra parameters are required. Selection, Crossover  and Mutation types.\n";
	std::cout << "\t109. Basic Genetic + Complete MDVH2+3OPT. Extra parameters are required. Selection, Crossover  and Mutation types.\n";

	std::cout << "\t111. Hybrid Genetic(Huang and Lim). LSH based on Hungarian.\n";
	std::cout << "\t112. Hybrid Genetic(Huang and Lim). LSH based on FlowAssign.\n";
	std::cout << "\t113. Hybrid Genetic(Huang and Lim). LSH based on Auction.\n";

	std::cout << "\n<iterations>: The number of iterations to execute.\n";
	std::cout << "\n<seed>: The seed to use, a negative value denotes a random seed.\n";
	std::cout << "\n<seconds>: True if seconds should be considered.\n";
	std::cout << "\n<K or population size for memetic>(optional): The value of K for those methods that required it.\n";

	std::cout << "\n\n For the case of the of the memetic algorithm 11 parameters are required, the previous and the next ones:\n";
	std::cout << "\n<Selection>(optional): The selection method (1. Elitist, 2. Tournament, 3. Roulette).\n";
	std::cout << "\n<Crossover>(optional): The crossover method (1. PMX, 2. Cycled, 3. Ordered).\n";
	std::cout << "\n<Mutation>(optional): The mutation method (1. Swapping, 2. Inversion. 3 Permutation).\n";
	exit(0);
}

int main(int argc, char **argv){
	int K_value = 2;
	clock_t comienzo = clock();
	if(argc < 8 || atoi(argv[4]) >= 100 && argc < 12){
		use();
	}
	int seed = atoi(argv[6]);
	srand(seed < 0 ? (unsigned int)time(0) : seed);
	int inputType = atoi(argv[1]);
	int localSearch = atoi(argv[4]);
	int iterations = atoi(argv[5]);
	bool seconds = atoi(argv[7]);
	int mutation, cross, selection, population;
	if(argc >= 9)
		K_value = atoi(argv[8]);
	if((localSearch == 50 || localSearch == 70) && argc < 9){
		printf("The k value is missing\n");
		exit(0);
	}
	if(localSearch >= 100){
		if(argc < 12){
			use();
		}
		selection = atoi(argv[9]);
		cross = atoi(argv[10]);
		mutation = atoi(argv[11]);
		population = atoi(argv[8]);
	}
	printf("k_value %d\n", K_value);
	solve(inputType, localSearch, iterations, seconds, K_value, argv[2], argv[3], selection, cross, mutation, population);
	std::cout << "The running total time was "<<((clock()-comienzo)/(double)CLOCKS_PER_SEC)<<" seconds.\n";
	return 0;
}
