// MyEvalData.h: interface for the MyEvalData class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MYEVALDATA_H__0AB2E77B_2A27_4B1D_A9EA_261DAE225932__INCLUDED_)
#define AFX_MYEVALDATA_H__0AB2E77B_2A27_4B1D_A9EA_261DAE225932__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <ga/GAEvalData.h>
#include <dataStruct.h>


class MyEvalData : public GAEvalData {
private:
	SeedStatus seedStatus;

public:
	MyEvalData();
	MyEvalData(const GAEvalData& ed);
	virtual ~MyEvalData();
	GAEvalData* clone() const;
	void copy(const GAEvalData& ed);

public:
	void setSeedStatus(SeedStatus status);
	SeedStatus getSeedStatus( );
};


#endif // !defined(AFX_MYEVALDATA_H__0AB2E77B_2A27_4B1D_A9EA_261DAE225932__INCLUDED_)
