
#ifndef GRID_H
#define GRID_H
 
#include<cmath>
#include<cstdlib>
#include<vector>
#include<time.h>
#include<iostream>
#include<cassert>
#include<set>

#include "Common.h"
#include "pcgsolver/pcg_solver.h"

#define NX							100
#define NY							100
#define RHO							1.0
#define DX							0.1
#define DT							0.005
#define ETA							0.0003
#define NUMBER_OF_PARTICLES			0
#define FRAMES_FOR_PHI_RECALC		10

using namespace std;

typedef vector<double> doubleVector;				//1D array of doubles
typedef vector<doubleVector> doubleMatrix;			//2D array of doubles

typedef vector<bool> boolVector;					//1D array of bools
typedef vector<boolVector> boolMatrix;				//2D array of bools

typedef vec<int> intVec;							//2-int object
typedef vector<intVec> intVecVector;				//1D array of 2-int objects

typedef vec<double> doubleVec;						//2-double object
typedef vector<doubleVec> doubleVecVector;			//1D array of 2-double objects
typedef vector<doubleVecVector> doubleVecMatrix;	//2D array of 2-double objects

const double g = -9.81; 							//gravity

struct indexPair {
	int i, j;

	inline bool operator==(indexPair ip) const {
		return ((i == ip.i) and (j == ip.j));
	}
};

template <typename T> class bag {
	private:
		vector<T> elements;

	public:
		bool has(T elem) {
			for(unsigned int n=0; n<elements.size(); n++)
				if (elements[n] == elem)
					return true;
			return false;
		}
		bool insert(T elem) {
			for(unsigned int n=0; n<elements.size(); n++)
				if (elements[n] == elem)
					return false;
			elements.push_back(elem);
			return true;
		}
		unsigned int size() {
			return elements.size();
		}
		const T& element(int n) const {
			return elements[n];
		}
		inline const T& operator[](int n) const {
			return elements[n];
		}
};

class Grid {
	static const vec<double> uOffset;
	static const vec<double> vOffset;

	private:
		void loadParameters();
		void initialConditions();
		void printGlobal();
		template<typename T> void printMatrix(vector<vector<T> >& q);

		int nx, ny, timeStep;
		double rho, dx, dt, eta;
		doubleMatrix p, u, v;
		doubleMatrix unew, vnew;

		//integration routines
		indexPair backToBox(int i, int j, int field);
		bool checkInbox(int i, int j, int field);
		double linear1D(double q1, double q2, double alpha);
		double catmull1D(double q0, double q1, double q2, double q3, double a);
		double catmullInterpolate(int iP, int jP, double alpha, double beta, const doubleMatrix& q, int field);
		double linearInterpolate(int iP, int jP, double alpha, double beta, const doubleMatrix& q, int field);
		double interpolate(const vec<double>& rP, const vec<double>& offset, const doubleMatrix& q, int field);
		double advectField(int i, int j, double _dt, const doubleMatrix& q, int field, const vec<double>& offset);
		void advect(double _dt);
		void applyGravity();

		//checkers
		bool CFLcondition();
		bool insidePhiBox(int i, int j);
		bool insideUBox(int i, int j);
		bool insideVBox(int i, int j);
		const vec<double> cellCenter(int i, int j);

		//level set
		doubleVecMatrix closestPoint;
		boolMatrix known;
		doubleMatrix phi, phinew;
		void initDistanceCalculation();
		void calcNeighbourMinDistance(int i, int j);
		void recalcSignedDistance(int N);

		//projection
		SparseMatrix<double> matA;
		PCGSolver<double> pcgSolver;
		doubleVector rhs, vecP;
		int index(int i, int j);
		void calcAMatrix();
		void calcRHS();
		void pressureUpdate();
		void project();

		//viscosity
		void viscosity();
		
		//cell types and boundaries
		boolMatrix activeU, activeV;
		boolMatrix activeUnew, activeVnew;
		vector<indexPair> solidCells;
		vector<indexPair> fluidCells;
		boolMatrix isSolid;
		void clearNonFluidFaces();
		void solidBoundaries();
		void boxBoundaries();
		bool isFluid(int i, int j);
		bool hasFluidNeighbour(int i, int j);

		//extrapolation
		void defineFluidCellsForPaint();
		void defineActiveCells();
		bool extrapolateU(int i, int j);
		bool extrapolateV(int i, int j);
		void extrapolateCell(int i, int j);
		void extrapolateEasy(int unsigned N);

		//particles
		doubleVecVector particlePosition;
		int numberOfParticles;
		bool particlesAreVisible;
		void particleInit(int n);
		void particleAdvection();

	public:
		Grid();
		~Grid();
		void resetGrid();
		void update();

		//level set
		double getPhi(int i, int j);
		bool cellIsKnown(int i, int j);
		vec<double> getClosestPoint(int i, int j);

		//cell types
		const indexPair& getSolidCell(int n);
		const indexPair& getFluidCell(int n);
		int getNumberOfSolidCells();
		int getNumberOfFluidCells();
		bool isSolidCell(int i, int j);
		void setSolidCell(int i, int j);

		//particles
		const vec<double>& getParticlePosition(int n);
		int getNumberOfParticles();
		void setParticlesVisible(bool b);

		//accessors
		double maxU();
		double energy();
		double fluidArea();
		vec<double> massCenter();
		double divergence(int i, int j);
		int getTimeStep() { return timeStep; }
		int getNX() { return nx; }
		int getNY() { return ny; }
		double getDX() { return dx; }
		double getDT() { return dt; }
		void setU(int i, int j, double _u, double _v);
		vec<double> getU(int i, int j);
};

#endif
