// MyEvalData.cpp: implementation of the MyEvalData class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MyEvalData.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MyEvalData::MyEvalData() : GAEvalData(){
	seedStatus = DEUnknown;
}

MyEvalData::MyEvalData(const GAEvalData& ed): GAEvalData(ed){
	const MyEvalData* pmed = dynamic_cast<const MyEvalData*>(&ed);
	seedStatus = pmed->seedStatus;
}

MyEvalData::~MyEvalData() {

}

void MyEvalData::copy(const GAEvalData& ed) {
	const MyEvalData* pmed = dynamic_cast<const MyEvalData*>(&ed);
	seedStatus = pmed->seedStatus;
}


GAEvalData* MyEvalData::clone() const {
	MyEvalData* ned = new MyEvalData();
	ned->seedStatus = seedStatus;
	return ned;
}

void MyEvalData::setSeedStatus(SeedStatus status) {
	seedStatus = status;
}

SeedStatus MyEvalData::getSeedStatus( ){
	return seedStatus;
}

