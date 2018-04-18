#ifndef INPUT
#define INPUT
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#define REDEFINED_MATRIX "REDEFINED_MATRIX"

namespace Input{

	void verify(std::ifstream &input, const char *fileName){
		if(!input.is_open()){
			std::cout << "File " << fileName << " can not be opened or it does not exist.\n";
			exit(0);
		}
		return;
	}

	void verify(FILE *input, const char *name){
		if(input == NULL){
			printf("File '%s' can not be opened or it does not exist.\n", name);
			exit(0);
		}
		return;
	}

	int getInteger(unsigned char s[]){
		int value = 0;
		for(int i = 3; i >= 0; i--){
//			printf("%d\n", s[i]);
			value<<=8;
			value += s[i];
		}
		return value;
	}

	void ReadBalasSaltzman(const char *nameFile, CompleteHyperGraph &G){
		std::cout << "Start ReadBalasandSaltzman\n";
		FILE *input = fopen(nameFile, "r");
		
		verify(input, nameFile);

		int numberOfDimensions = 0, vertices;
		if(fscanf(input, "%d%d", &numberOfDimensions, &vertices) != 2 || !(2 < numberOfDimensions && numberOfDimensions < 15)){
			std::cout << "Error: the input cannot be obtained or dimensions are wrong dim = " << numberOfDimensions << "\n";
//			std::exit(0);
		}//*/
		std::vector<int> dimensions(numberOfDimensions, vertices);
		G.setDimensions(dimensions);

		printf("dim %d vertices %d\n", numberOfDimensions, vertices);

		edgeCost cost;
		for(int i = 0; i < G.getNumberOfEdges(); i++){
			fscanf(input, "%d", &cost);
			G.addEdge(cost);
		}
		std::cout << "End ReadBalasandSaltzman\n";
		return;
	}

	void ReadBinaryAsKarapetyanGutin(const char *nameFile, CompleteHyperGraph &G){
		std::cout << "Start ReadBinaryAsKarapetyanGutin\n";
		FILE *input = fopen(nameFile, "rb");
		
		verify(input, nameFile);

		unsigned char buffer[10];

		int numberOfDimensions = 0, vertices;
		size_t result;
		result = fread (buffer,1,4,input);
		numberOfDimensions = getInteger(buffer);
		result = fread (buffer,1,4,input);
		vertices = getInteger(buffer);
		printf("%d %d\n", numberOfDimensions, vertices);
		if(!((3 <= numberOfDimensions && numberOfDimensions <= 12) && (1 <= vertices && vertices <= 1000))){
			std::cout << "Instance is too big s = "<<numberOfDimensions<<" edges = "<<vertices<<"\n";
			exit(0);
		}
		std::vector<int> dimensions (numberOfDimensions, vertices);

		G.setDimensions(dimensions);

		printf("dim %d vertices %d\n", numberOfDimensions, vertices);

		for(int i = 0; i < G.getNumberOfEdges(); i++){
			result = fread (buffer,1,4,input);
			G.addEdge(getInteger(buffer));
		}
		std::cout << "End ReadBinaryAsKarapetyanGutin\n";
		return;
	}


	void ReadAsIndependent(const char *nameFile, HyperGraph &G){
		std::cout << "Start ReadAsIndependent\n";
		int numberOfDimensions, nTuples;
		edgeCost cost;	

		FILE *input = fopen(nameFile, "r");

		verify(input, nameFile);

		fscanf(input, "%d", &numberOfDimensions);
		std::vector<int> dimensions(numberOfDimensions, 0);
		
		for(unsigned int dim = 0; dim < dimensions.size(); dim++){
			fscanf(input, "%d", &dimensions[dim]);
		}

		G.setDimensions(dimensions);

		fscanf(input, "%d", &nTuples);

		printf("dimensiones %d\n", G.getNumberOfDimensions());
		std::vector<unsigned long long int> vertices(G.getNumberOfDimensions(), 0);
		for(int i = 0; i < nTuples; i++){			
			for(int j = 0; j < G.getNumberOfDimensions(); j++)//reading the involved vertices at the hyperedge
				fscanf(input, "%llu", &vertices[j]);
			fscanf(input, "%d", &cost);
			G.addEdge(vertices, cost);
		}
		std::cout << "End ReadAsIndependent\n";
		return;
	}

}
#endif