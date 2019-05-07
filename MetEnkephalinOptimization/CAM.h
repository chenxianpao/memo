// CAM.h: interface for the CCAM class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CAM_H__758F551F_B8CF_4C75_8299_26D6EADB033C__INCLUDED_)
#define AFX_CAM_H__758F551F_B8CF_4C75_8299_26D6EADB033C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "config.h"


/**********************************************************************************************/
// ֧��������
// �ڴ�ռ���������������� �麯��ָ��4B+��Ա����12B = 16B
//			   ��չ������ sizeof(double)��m_dim			 
/**********************************************************************************************/
class Cl {
public:
	Cl(long id, double* x, int n) {			
		m_id = id;
		m_dim = n;
		m_element = new double[n];
		for(int i=0; i<n;i++) {m_element[i]=x[i];}
	}
	
	Cl(const Cl& arg) {
		m_id = arg.m_id;
		m_dim = arg.m_dim;
		m_element = new double[m_dim];
		for(int i=0; i<m_dim; i++) {m_element[i] = arg.m_element[i];}
	}

	virtual void copy(const Cl& arg) {
		m_id = arg.m_id;
		m_dim = arg.m_dim;
		if (m_element !=0) {delete [] m_element; m_element = 0;	}
		m_element = new double [m_dim];
		for(int i=0; i<m_dim; i++) {m_element[i] = arg.m_element[i];}
	}

	virtual ~Cl( ){
		if (m_element!=0) {
			delete [] m_element;
			m_element=0;
		}
	}

	Cl & operator=(const Cl & arg) {copy(arg); return (*this); }

protected:
	long    m_id;		// ������ʶ������1��ʼ
	int	    m_dim;		// ����Ԫ�ظ�����
	double  *m_element;	// ����Ԫ�����飻

public:
	long getId( )  { return m_id; }
	int  getDim( ) { return m_dim;}
	double* getElement( ) {
		return m_element;
	}
	double  getElementData(int index) {
		return m_element[index];
	}
	double  setElementData(int index, double e) {
		return m_element[index] = e;
	}
};

/**********************************************************************************************/
// ClArray: ֧�����������ࡣΪ���������ڴ����ģ�ֻ��ClArray�б���֧�ű���������Ϊָ������÷�ʽ
// �ڴ�ռ���������������� �麯��ָ��4B+��Ա����20B = 24B
/**********************************************************************************************/
class ClArray {
public:
	ClArray(int n) {
		m_dim = n;
		m_count = 0;
		m_chunksize = 200;
		m_allocatedsize = m_chunksize;
		m_l = new Cl* [m_allocatedsize];
	}
	virtual ~ClArray( ) {
		if (m_l !=0) {
			for(int i=0; i<m_count; i++) {
				delete [] m_l[i];
			}
			delete [] m_l;
			m_l = 0;
		}
	}

public:
	int		m_dim;					// ����ά��
	long	m_count;				// ֧��������
	long	m_chunksize;			// ÿ�οռ������
	long	m_allocatedsize;		// �Ѿ�����ռ�
	Cl		**m_l;					// ֧������ָ������

public:
	long	getDim( )			{ return m_dim;}
	long	getCount( )			{ return m_count;}
	long	getChunkSize( )		{ return m_chunksize;}
	long	getAllocatedSize( ) { return m_allocatedsize;}
	Cl**	getClArray( )		{ return m_l;}

	void add(Cl& arg) {
		grow(m_count);
		m_l[m_count] = new Cl(arg);
		m_count ++;
	}
protected:
	long grow(long s) {
		if (s<m_allocatedsize) return m_allocatedsize;
		
		long oldsize = m_allocatedsize;
		while(m_allocatedsize<=s) m_allocatedsize += m_chunksize;
		
		Cl **temp = m_l;
		m_l = new Cl * [m_allocatedsize];
		for(int i=0; i<m_count; i++) {
			m_l[i]  = temp[i];
		}
		delete [] temp;
		return m_allocatedsize;
	}
};

/**********************************************************************************************/
// ��Ҷ��: �����������Ӧ֧������ָ�룬��������֧������ʵ�壬�û�Ӧ�ñ�����Ӧ֧����������������
/**********************************************************************************************/
class CL {
public:
	CL(Cl** sl, int dim) {
		b_calculated = false; 
		b_condition1 = false;
		m_condition1 = false; 
		m_min		 = -1; 
		m_dim		 = dim;
		m_l			 = new Cl* [m_dim];
		for(int i=0; i<m_dim; i++) {m_l[i] = sl[i];}
	}
	CL(const CL& arg) {
		m_dim		 = arg.m_dim;
		m_min		 = arg.m_min;
		m_condition1 = arg.m_condition1;
		b_condition1 = arg.b_condition1;
		b_calculated = arg.b_calculated;
		m_l			 = new Cl* [m_dim];
		for(int i=0; i<m_dim; i++) {m_l[i] = arg.m_l[i];}
	}

	~CL( ){					// ��������ָ��Ӧ��֧����������������������ָ�롣
		delete [] m_l;
	}

public:
	Cl		**m_l;			// ����ָ��֧��������ָ�룻
	double	m_min;			// ��Ӧ�ɳ�������Сֵ��
	bool	b_calculated;	// �Ƿ��Ѽ�������ֵ��
	bool	m_condition1;	// �Ƿ���������I�����Խ�֧����(diagonal dominance)��
	bool	b_condition1;	// ����I�Ƿ���ԣ�
	int		m_dim;			// ����֧������ά���� 

public:
	double getMinimum( ) {
		if (b_calculated) {
			return m_min;
		}
		m_min = 0;
		for(int i=0; i<m_dim; i++) {
			m_min += m_l[i]->getElementData(i);	// ����Խ�Ԫ��֮��
		}
		m_min = 1/m_min;
		
		b_calculated = true;
		return m_min;
	}

	bool getCondition1( ) {
		if (b_condition1) { return m_condition1;}

		m_condition1 = true;
		for(int i=0; i<m_dim; i++) {
			double a = m_l[i]->getElementData(i);
			for (int j=0; j<m_dim; j++) {
				if (a < m_l[j]->getElementData(i) && i!=j) {
					m_condition1 = false;
					break;
				}
			}
			if (!m_condition1) break;
		}

		b_condition1 = true;

		return m_condition1;
	}

	double getDialogElement(int i) {
		return m_l[i]->getElementData(i);
	}
	
	long getClLable(int i) {
		return m_l[i]->getId();
	}
};


class CLLeafArray {
public:
	int			 m_dim;				// ����ά����
	long		 m_count;			// ��Ҷ����Ŀ;
	CL			 **m_arrayCL;		// ָ�����飬ָ����Ҷ;
	CL			 *m_minimumCL;		// ָ����СҶ
	unsigned int m_chunksize;		// ÿ�οռ������Ŀ;
	unsigned int m_allocatedsize;	// �Ѿ�����ռ��С;

public:
	CLLeafArray(int n) {
		m_dim =n;
		m_count = 0;
		m_minimumCL = 0;
		m_chunksize = 200;
		m_allocatedsize = m_chunksize;
		m_arrayCL = new CL* [m_allocatedsize];
	}
	virtual ~CLLeafArray( ) {
		delete [] m_arrayCL;
	}

public:
	long getCount( ) { return m_count;}
	int  getDim( ) { return m_dim;}

	unsigned int grow(unsigned int s) {
		if (s<m_allocatedsize) return m_allocatedsize;

		unsigned int oldsize = m_allocatedsize;
		while(m_allocatedsize <= s) m_allocatedsize += m_chunksize;

		CL **temp = m_arrayCL;
		m_arrayCL = new CL * [m_allocatedsize];
		for(int i=0; i<m_count; i++) {
			m_arrayCL[i] = temp[i];
		}
		delete [] temp;

		return m_allocatedsize;
	}

	void add(CL &L) {
		grow(m_count);
		m_arrayCL[m_count] = new CL(L);
		m_count++;
	}

	CL* findMinimumCL( ) {
		m_minimumCL = m_arrayCL[0];
		for(int i=1; i<m_count; i++) {
			if (m_minimumCL->getMinimum( ) > m_arrayCL[i]->getMinimum( )) {
				m_minimumCL=m_arrayCL[i];
			}
		}
		return m_minimumCL;
	}

	void remove(long i) {
		if(m_count==1) {m_count--; return;}
		for(long j=i; j<m_count-1; j++) {
			m_arrayCL[j] = m_arrayCL[j+1];
		}
		m_count--;
	}

	void QuickSortAscending(CL **c, int l, int r) {
		int i,j; double v; CL* t;
		if(r > l){
			v = c[r]->getMinimum(); i=l-1; j=r;
			for(;;){
				while(c[++i]->getMinimum() < v && i <= r);
				while(c[--j]->getMinimum() > v && j > 0);
				if(i >= j) break;
				t = c[i]; c[i] = c[j]; c[j] = t;
			}
			t = c[i]; c[i] = c[r]; c[r] = t;
			QuickSortAscending(c,l,i-1);
			QuickSortAscending(c,i+1,r);
		}
	}

	CL* findMinimumLeaf( ) {
		QuickSortAscending(m_arrayCL, 0, m_count-1);
		return m_arrayCL[0];
	}
};


class CCAM  
{
public:
	typedef double (*Evaluator)(double *x, int dim);
	static double function(double*x, int dim);
public:
	double	energyFunction(double *x, int dim);

	double	singleEnergyMinimize(double *pGEOM, int length);
	void	writeSingleEnergyMinimizeSourceFile(const char* sourceFile, double* pGEOM, int length);
	double	readSingleEnergyMinimizeResultETOT(const char* resultFile, double*pGEOM, int length);
	void	clearSingleEnergyMinimizeResultFiles( );
public:
	CCAM();
	virtual ~CCAM();
public:
	void updateLeafArray(CLLeafArray& leafArray, Cl& l);
	bool illposed(CL& arg);
public:
	void optimize( );
};

#endif // !defined(AFX_CAM_H__758F551F_B8CF_4C75_8299_26D6EADB033C__INCLUDED_)
