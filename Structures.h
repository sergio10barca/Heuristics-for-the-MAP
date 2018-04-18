#ifndef STRUCTURES
#define STRUCTURES
#include <map>
#include <vector>
#include <utility>
#include <iostream>
#include <cstdlib>
#include <climits>
#include "LinearAssignment/StructureAP.h"
typedef int edgeCost;

/*
dTuple:: each four bits are equal to the index of a vertex,
we can represent as much as 8 vertices (then dimensions) by word.
*/
typedef std::vector<std::vector<edgeCost> > matriz;
typedef std::vector<int> dTuple;
//typedef unsigned long long int dataTypeHyperEdge;//for a maximum of 8 dimensions of size 8 by each
class HyperGraph{
	typedef unsigned long long int dataTypeHyperEdge;//for a maximum of 8 dimensions of size 8 by each
private:
	std::map<dataTypeHyperEdge, edgeCost> hyperGraph;
	std::vector<int> dimensions;
	std::vector<unsigned long long int> offset;
public:
	HyperGraph(){

	}
	~HyperGraph(){
		std::cout << "entra aqui y destruye\n";
		hyperGraph.clear();
		dimensions.clear();
		offset.clear();
	}
	unsigned long long int getOffset(int dim){
		unsigned long long int powerOfTwo = 1, offset = 0;
		for(; powerOfTwo < dim; powerOfTwo<<=1, offset++);
		return offset;
	}

	void getOffsetVector(std::vector<int> dimensions){
		for(unsigned int dim = 0; dim < dimensions.size(); dim++)
			offset.push_back(getOffset(dimensions[dim]));
		return;
	}

	void setDimensions(const std::vector<int> &dimensions){
		this->dimensions = dimensions;
		this->getOffsetVector(dimensions);
		return;
	}

	std::vector<int> getDimensions(){
		return this->dimensions;
	}

	int getDimension(int i){
		return (int)this->dimensions.size() > i ? this->dimensions[i] : -1;
	}

	int getNumberOfDimensions(){
		return (int)this->dimensions.size();
	}

	int getNumberOfEdges(){
		return (int)this->hyperGraph.size();
	}

	void addEdge(std::vector<unsigned long long int> &v, int cost){
		if(v.size() != this->getNumberOfDimensions()){
			std::cout << "Error: the vertex array doesnot have the required size\n";
			exit(0);
		}
		dataTypeHyperEdge edge = 0;
		for(unsigned int i = 0; i < v.size(); i++){
			edge <<= offset[i];
			edge |= v[i];
		}
		this->hyperGraph[edge] = cost;
		return;
	}

	dataTypeHyperEdge getHyperEdge(const dTuple &tuple){
		dataTypeHyperEdge hyperEdge = 0;
		for(unsigned int i = 0; i < tuple.size(); i++){
			hyperEdge <<= offset[i];
			hyperEdge |= tuple[i];
		}
		return hyperEdge;
	}

	dTuple getEdge(int index){
		dTuple tuple(this->dimensions.size(), 0);
		return tuple;
	}	

	edgeCost getEdgeCost(const int hyperEdge){
		std::map<dataTypeHyperEdge, edgeCost>::iterator it = hyperGraph.find(hyperEdge);
		if(it == hyperGraph.end()){
			std::cout << "Error: the edge doesnot exit\n";
			exit(0);
		}
		return it->second;
	}

	edgeCost getEdgeCost(const dTuple &edge){
		std::map<dataTypeHyperEdge, edgeCost>::iterator it = hyperGraph.find(this->getHyperEdge(edge));
		if(it == hyperGraph.end()){
			std::cout << "Error: the edge doesnot exit\n";
			exit(0);
		}
		return it->second;
	}

	edgeCost getEdgeCost(const int dim1, const int dim2, const int dim3){
		dTuple edge(3, 0);
		edge[0] = dim1;
		edge[1] = dim2;
		edge[2] = dim3;
		return this->getEdgeCost(edge);
	}

	edgeCost getEdgeCost(const int dim1, const int dim2, const int dim3, const int dim4){
		dTuple edge(4, 0);
		edge[0] = dim1;
		edge[1] = dim2;
		edge[2] = dim3;
		edge[3] = dim4;
		return this->getEdgeCost(edge);
	}

	bool isEdge(dTuple &edge){
		return hyperGraph.find(this->getHyperEdge(edge)) != hyperGraph.end() ? true: false;
	}
};

class CompleteHyperGraph{
private:
	int numberOfEdges, next;
	edgeCost *completeGraph;
	std::vector<int> dimensions;
public:

	CompleteHyperGraph(){
		this->completeGraph = NULL;
	}

	CompleteHyperGraph(const int dimensions, const int n){
		this->dimensions = std::vector<int>(dimensions, n);
		this->setDimensions(this->dimensions);
	}

	~CompleteHyperGraph(){
		if(this->completeGraph != NULL)
			delete[] this->completeGraph;
		this->dimensions.clear();
	}

	void reStart(){
		this->next = 0;
		return;
	}

	void setDimensions(const int s, const int n){
		this->dimensions.clear();
		this->numberOfEdges = 1;
		for(int i = 0; i < s; i++){
			this->dimensions.push_back(n);
			this->numberOfEdges *= n;
		}
		this->completeGraph = new edgeCost[this->numberOfEdges];
		this->next = 0;
		return;
	}

	void setDimensions(const std::vector<int> &dimensions){
		this->dimensions = dimensions;
		this->numberOfEdges = 1;
		for(unsigned int i = 0; i < dimensions.size(); i++)
			this->numberOfEdges *= dimensions[i];
		this->completeGraph = new edgeCost[this->numberOfEdges];
		this->next = 0;
		return;
	}

	std::vector<int> getDimensions(){
		return this->dimensions;
	}

	int getDimension(int i){
		return (int)this->dimensions.size() > i ? this->dimensions[i] : -1;
	}

	int getNumberOfDimensions(){
		return (int)this->dimensions.size();
	}

	int getNumberOfEdges(){
		return this->numberOfEdges;
	}

	void addEdge(edgeCost cost){
//		cout << this->next << " " << cost << endl;
		this->completeGraph[this->next++] = cost;
	}

	dTuple getEdge(int index){
		dTuple tuple(this->dimensions.size(), 0);
		for(int i = (int)this->dimensions.size()-1; i >= 0; i--){
			tuple[i] = index%(this->dimensions[i]);
			index /= this->dimensions[i];
		}
		return tuple;
	}

	edgeCost getEdgeCost(const int index){
		return this->completeGraph[index];
	}

	edgeCost getEdgeCost(const dTuple &edge){
		int index = 0, prod = 1;
		for(int i = (int)edge.size()-1; i >= 0; i--){
			index += edge[i]*prod;
			prod *= this->getDimension(i);
		}
		return this->completeGraph[index];
	}

	edgeCost getEdgeCost(const int dim1, const int dim2, const int dim3){
		return this->completeGraph[dim1*this->getDimension(1)*this->getDimension(2)+dim2*this->getDimension(2)+dim3];
	}

	edgeCost getEdgeCost(const int dim1, const int dim2, const int dim3, const int dim4){
		return this->completeGraph[dim1*this->getDimension(1)*this->getDimension(2)*this->getDimension(3)
			+dim2*this->getDimension(2)*this->getDimension(3)+dim3*this->getDimension(3)+dim4];
	}

	bool isEdge(dTuple &edge){
		int index = 0, prod = 1;
		for(int i = (int)edge.size()-1; i >= 0; i--){
			index += edge[i]*prod;
			prod *= this->getDimension(i);
		}
		return index < this->numberOfEdges;
	}

};

class Assignment{
	public:
	std::vector<dTuple> matching;
	edgeCost cost;

	Assignment(){

	}
	Assignment(const std::vector<dTuple> &matching, const int cost){
		this->matching = matching;
		this->cost     = cost;
	}
	Assignment(int cost){
		this->cost = cost;
	}
	~Assignment(){
		this->matching.clear();
	}

	bool operator<(const Assignment &assignment) const {
		return this->cost<assignment.cost;
	}
	bool operator>(const Assignment &assignment) const {
		return this->cost>assignment.cost;
	}
	bool operator==(const Assignment &assignment) const {
		if(this->cost!=assignment.cost){
			return false;
		}

		if(this->matching.size() != assignment.matching.size()){
			return false;
		}

		for(unsigned int i = 0; i < this->matching.size(); i++){
			if(this->matching[i].size() != assignment.matching[i].size()){
				return false;
			}
			for(unsigned int j = 0; j < this->matching[i].size(); j++){
				if(this->matching[i][j] != assignment.matching[i][j]){
					return false;
				}
			}
		}

		return true;
	}
};

//*/
#endif
