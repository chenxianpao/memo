
#if !defined(ZGJ_INFO_ZJUT_PDB_ATOM_ENTRY_INCLUDED)
#define ZGJ_INFO_ZJUT_PDB_ATOM_ENTRY_INCLUDED

#if _MSC_VER > 1000
#pragma once
#endif 

#include "string.h"
#include "stdlib.h"

typedef enum {
	H2N_H2 = 0,
	
	TYR_N,	
	TYR_NH2,  
	TYR_CA, 
	TYR_HA, 
	TYR_CB, 
	TYR_C, 
	TYR_O, 
	TYR_HB1, 
	TYR_HB2, 
	TYR_CG, 
	TYR_CD1, 
	TYR_CD2,
	TYR_HD1,
	TYR_CE1,
	TYR_CE2,
	TYR_HD2,
	TYR_HE1, 
	TYR_CZ, 
	TYR_HE2, 
	TYR_OH, 
	TYR_HH,
	
	GLY1_N, 
	GLY1_HN, 
	GLY1_CA, 
	GLY1_HA1, 
	GLY1_HA2, 
	GLY1_C, 
	GLY1_O,
	
	GLY2_N, 
	GLY2_HN, 
	GLY2_CA, 
	GLY2_HA1, 
	GLY2_HA2, 
	GLY2_C, 
	GLY2_O,
	
	PHE_N,	 
	PHE_HN,  
	PHE_CA,  
	PHE_HA,  
	PHE_CB,  
	PHE_C,   
	PHE_O,  
	PHE_HB1, 
	PHE_HB2, 
	PHE_CG,
	PHE_CD1, 
	PHE_CD2, 
	PHE_HD1, 
	PHE_CE1, 
	PHE_CE2, 
	PHE_HD2, 
	PHE_HE1, 
	PHE_CZ, 
	PHE_HE2, 
	PHE_HZ,
	
	MET_N,   
	MET_HN,  
	MET_CA,  
	MET_HA,  
	MET_CB,  
	MET_C,   
	MET_O,  
	MET_HB1, 
	MET_HB2, 
	MET_CG,
	MET_HG1, 
	MET_HG2, 
	MET_SD,  
	MET_CE,  
	MET_HE1, 
	MET_HE2, 
	MET_HE3,
	
	CXH_O,   
	CXH_H	
} ATOMTYPE;



class CAtomEntry  
{
public:
	CAtomEntry(){ 
		strncpy(recName,    "\0", 8);
		strncpy(atomName,   "\0", 8);
		strncpy(residueName,"\0", 8);
		strncpy(altLoc,		"\0", 2);
		strncpy(chainID,	"\0", 2);
		strncpy(insertCode, "\0", 2);
		x = y = z = occupancy = tempFactor = 0.0;
	}
	virtual ~CAtomEntry() {};

private:
	char	recName[8];		// record name. Format: A6; From-End: 0-5; Description: a literal "ATOM  " (note two trailing spaces). 
	int		serial;			// atom serial number. Format: I5; From-End: 6-10;
	char	atomName[8];	// atom role name. Format: A4; From-End: 12-15;
	char	altLoc[2];		// alternate location indicator. Format:A1; From-End: 16-16;
	char	residueName[8];	// residue name. Format:A3; From-End:17-19;
	char	chainID[2];		// chain ID. Format:A1; From-End:21-21;
	int		seqNo;			// residue sequence number. Format:I4; From-End:22-25;
	char	insertCode[2];	// insertion code. Format:A1; From-End: 26-26;
	double	x;				// atom X coordinate. Format:F8.3; Start-End: 30-37;
	double	y;				// atom Y coordinate. Format:F8.3; Start-End: 38-45;
	double	z;				// atom Z coordinate. Format:F8.3; Start-End: 46-53;
	double	occupancy;		// atom occupancy. Format:F6.2; Start-End: 54-59;
	double	tempFactor;		// B value or temperature factor. Format:F6.2; Start-End:60-65;

public:
	void	setRecName(char* ch) { strncpy(recName, ch, 4);}
	char*	getRecName( ) {return recName; }

	void	setSerial(char* ch) {serial = atoi(ch);}
	int		getSerial( ) { return serial;}

	void	setAtomName(char *ch) {strncpy(atomName, ch, 4);}
	char*	getAtomName( ) { return atomName;}

	void	setAltLoc(char *ch) {strncpy(altLoc, ch, 1);}
	char*	getAltLoc( ) {return altLoc;}

	void	setResName(char* ch) {strncpy(residueName, ch, 3);}
	char*	getResName( ) {return residueName;}

	void	setChainID(char* ch) {strncpy(chainID, ch, 1);}
	char*	getChainID( ) {return chainID;};

	void	setSeqNo(char* ch) {seqNo = atoi(ch); }
	int		getSeqNo( ) {return seqNo;}

	void	setInsertCode(char* ch) {strncpy(insertCode,ch,1);}
	char*	getInsertCode( ) {return insertCode;}

	void	setX(char* ch) {x = atof(ch);}
	double	getX( ) {return x;}

	void	setY(char* ch) {y = atof(ch);}
	double	getY( ) {return y;}

	void	setZ(char* ch) {z = atof(ch);}
	double	getZ( ) {return z;}

	void	setOccupancy(char* ch) {occupancy = atof(ch);}
	double	getOccupancy( ) {return occupancy;}

	void	setTempFactor(char* ch) {tempFactor = atof(ch);}
	double	getTempFactor( ) {return tempFactor;}
};
#endif 
