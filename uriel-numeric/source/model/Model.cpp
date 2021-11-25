#include "pch.h";
#include "model/Model.h";

Model::Model(int nStates, int nAlgs, int nParams) {
	this->nStates = nStates;
	this->nAlgs = nAlgs;
	this->nParams = nParams;
	this->descr = new ModelVariableDescription[nStates];
}

void Model::setStateDescr(int index, std::string name, std::string unit) {
	descr[index].name = name;
	descr[index].unit = unit;
}

std::string Model::getStateName(int index) {
	return descr[index].name;
}

std::string Model::getStateUnit(int index) {
	return descr[index].unit;
}