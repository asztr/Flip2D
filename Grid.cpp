#include "Grid.h"

const vec<double> Grid::uOffset(0.0, 0.5);
const vec<double> Grid::vOffset(0.5, 0.0);

const indexPair nbours4[5] = {{1,0}, {-1,0}, {0,1}, {0,-1}, {0,0}};
const indexPair nbours8[9] = {{-1,1}, {0,1}, {1,1}, {-1,0}, {1,0}, {-1,-1}, {0,-1}, {1,-1}, {0,0}};

//cosas importantes para tener en cuenta sobre los bordes:
//recordar que phi[i][j] solo define adentro y afuera del fluido.
//pero un solido puede tener phi negativo (o sea ser considerado
//como parte interna del fluido. no sé a qué se debe eso.
//o sea que phi<0 no asegura que sea fluido la celda.
//por ahi se podria pensar en unas funciones isFluid, isSolid, isAir.
//para la extrapolación tuve que desencapsular isFluid y poner
//phi < 0, porque se trata de extrapolar adentro de un solido.

//tambien hay que tener en cuenta que muchas veces para condicionar
//las velocidades, hay que fijarse en las condiciones de ambas celdas
//que limitan con la cara de la velocidad.

//todavia queda algun que otro problema con los bordes.
//por ejemplo, el liquido parece quedarse pegado a las paredes un poco,
//y sobre todo al techo.

//otro problema: si pongo todo lleno de agua, o todo lleno de aire,
//se rompe todo. parece un problema del solver de presión vinculado
//con los bordes. pero se soluciona desactivando la extrapolación o
//el recálculo de distancia.

//escribí una rutina de viscosidad, pero todavía necesita unos cuantos
//retoques. habría que lidiar con los bordes, y en realidad armar mejor
//la matriz. por ejemplo, hay que tener en cuenta que la matriz A para
//u es distinta a la matriz para v. y ambas son distintas a la de la
//presion. por ejemplo, index se calcula distinto en cada situación.

//initialization
Grid::~Grid() {
}
Grid::Grid() {
	loadParameters();
	u = doubleMatrix(nx+1, doubleVector(ny+1));
	v = doubleMatrix(nx+1, doubleVector(ny+1));
	unew = doubleMatrix(nx+1, doubleVector(ny+1));
	vnew = doubleMatrix(nx+1, doubleVector(ny+1));
	p = doubleMatrix(nx+1, doubleVector(ny+1));
	isSolid = boolMatrix(nx+1, boolVector(ny+1));

	//projection
	matA = SparseMatrix<double>((nx+1)*(ny+1));
	rhs = doubleVector((nx+1)*(ny+1));
	vecP = doubleVector((nx+1)*(ny+1));

	//level set
	closestPoint = doubleVecMatrix(nx+1, doubleVecVector(ny+1));
	known = boolMatrix(nx+1, boolVector(ny+1));
	phi = doubleMatrix(nx+1, doubleVector(ny+1));
	phinew = doubleMatrix(nx+1, doubleVector(ny+1));

	//extrapolation
	activeU = boolMatrix(nx+1, boolVector(ny+1));
	activeV = boolMatrix(nx+1, boolVector(ny+1));
	activeUnew = boolMatrix(nx+1, boolVector(ny+1));
	activeVnew = boolMatrix(nx+1, boolVector(ny+1));

	//particles
	if (numberOfParticles > 0)
		particlePosition = doubleVecVector(numberOfParticles);
	particlesAreVisible = false;

	//initial conditions
	srand(time(NULL));
	resetGrid();
	initialConditions();
}
void Grid::loadParameters() {
	nx = NX;
	ny = NY;
	rho = RHO;
	dx = DX;
	dt = DT;
	eta = ETA;
	numberOfParticles = NUMBER_OF_PARTICLES;
}
void Grid::initialConditions() {
	//initial settings for level set
	const double R = 1.0;
	const vec<double> rc(nx*dx/2.0, ny*dx/2.0);
	for(int i=1; i<nx-1; i++)
		for(int j=1; j<ny-1; j++) {
			vec<double> r = cellCenter(i,j);
			double d = (r-rc).norm();
			phi[i][j] = d-R;
// 			if (phi[i][j] <= 0.0)
// 				u[i][j] = 5.0;
		}
 	for(int i=0; i<nx; i++)
 		for(int j=0; j<ny/20; j++)
 			phi[i][j] = -1.0;
	//printMatrix(phi);
}
void Grid::resetGrid() {
	timeStep = 0;

	//clear cell border arrays
	for(int i=0; i<nx+1; i++)
		for(int j=0; j<ny+1; j++) {
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			unew[i][j] = 0.0;
			vnew[i][j] = 0.0;
			activeU[i][j] = false;
			activeV[i][j] = false;
		}

	//clear cell center arrays
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++) {
			p[i][j] = 0.0;
			phi[i][j] = 1.0;
			known[i][j] = false;
			closestPoint[i][j].define(0.0, 0.0);
		}

	//set rigid boundaries
	for(int j=0; j<ny; j++) {
		setSolidCell(0, j);
		setSolidCell(nx-1, j);
	}
	for(int i=0; i<nx; i++) {
		setSolidCell(i, 0);
		setSolidCell(i, ny-1);
	}

	//particles
	for (int n=0; n<numberOfParticles; n++)
		particleInit(n);
}

void Grid::update() {
	if (timeStep % FRAMES_FOR_PHI_RECALC == 0)
		recalcSignedDistance(2);

	defineActiveCells();
	extrapolateEasy(5);
	advect(dt);
	applyGravity();
	project();
	//viscosity();
	//project();
	solidBoundaries();
	boxBoundaries();
	defineFluidCellsForPaint();

	timeStep++;
}

//advection
indexPair Grid::backToBox(int i, int j, int field) {
	indexPair ret = {0, 0};
	ret.i = max(0, i);
	ret.j = max(0, j);

	if ((i >= nx) and (field != 0))
		ret.i = nx-1;
	if ((i >= nx) and (field == 0))
		ret.i = nx;

	if ((j >= ny) and (field != 1))
		ret.j = ny-1;
	if ((j >= ny) and (field == 1))
		ret.j = ny;
	return ret;
}
bool Grid::checkInbox(int i, int j, int field) {
	if (field == 0)
		return insideUBox(i, j);
	else if (field == 1)
		return insideVBox(i, j);
	else return insidePhiBox(i, j);
}
double Grid::linear1D(double q1, double q2, double alpha) {
	return alpha*q2 + (1.0-alpha)*q1;
}
double Grid::catmull1D(double q0, double q1, double q2, double q3, double a) {
	double a2 = a*a;
	double a3 = a2*a;
	double ret = q0*(-0.5*a+a2-0.5*a3) + q1*(1.0-2.5*a2+1.5*a3) + q2*(0.5*a+2.0*a2-1.5*a3) + q3*(-0.5*a2+0.5*a3);
	/*if (ret > max(q1,q2))
		return max(q1,q2);
	if (ret < min(q1,q2))
		return min(q1,q2);*/
	return ret;
}
double Grid::linearInterpolate(int iP, int jP, double alpha, double beta, const doubleMatrix& q, int field) {
		double qaux[2][2];
		double raux[2];
		for(int i=0; i<2; i++)
			for(int j=0; j<2; j++)
				if (checkInbox(iP+i, jP+j, field))
					qaux[i][j] = q[iP+i][jP+j];
				else {
					indexPair back = backToBox(iP+i, jP+j, field);
					qaux[i][j] = q[back.i][back.j];
				}

		for(int j=0; j<2; j++)
			raux[j] = linear1D(qaux[0][j], qaux[1][j], alpha);
		return linear1D(raux[0], raux[1], beta);
}
double Grid::catmullInterpolate(int iP, int jP, double alpha, double beta, const doubleMatrix& q, int field) {
		double qaux[4][4];
		double raux[4];
		for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
				if (checkInbox(iP-1+i, jP-1+j, field))
					qaux[i][j] = q[iP-1+i][jP-1+j];
				else
					qaux[i][j] = 0.0;

		for(int j=0; j<4; j++)
			raux[j] = catmull1D(qaux[0][j], qaux[1][j], qaux[2][j], qaux[3][j], alpha);
		return catmull1D(raux[0], raux[1], raux[2], raux[3], beta);
}
double Grid::interpolate(const vec<double>& rP, const vec<double>& offset, const doubleMatrix& q, int field) {
	vec<double> aux = rP/dx-offset;
	int iP = floor(aux.x);
	int jP = floor(aux.y);

	double alpha = aux.x - iP;
	double beta = aux.y - jP;
	return linearInterpolate(iP, jP, alpha, beta, q, field);
}
double Grid::advectField(int i, int j, double _dt, const doubleMatrix& q, int field, const vec<double>& offset = vec<double>(0.5, 0.5)) { //generic (but obfuscated) advection routine
	vec<double> rG((i + offset.x)*dx, (j + offset.y)*dx);
 	vec<double> uG(interpolate(rG, uOffset, u, 0), interpolate(rG, vOffset, v, 1));
	vec<double> rMid = rG - uG*0.5*_dt;
	vec<double> uMid(interpolate(rMid, uOffset, u, 0), interpolate(rMid, vOffset, v, 1));
	vec<double> rP = rG - uMid*_dt;
	//vec<double> rP = rG - uG*_dt;
	return interpolate(rP, offset, q, field);
}
void Grid::advect(double _dt) {
	for(int i=0; i<nx; i++) //check limits
		for(int j=0; j<ny; j++)
			phinew[i][j]  = advectField(i, j, _dt, phi, 2);

	for(int i=0; i<nx; i++) //check limits
		for(int j=0; j<ny; j++)
			phi[i][j] = phinew[i][j];

	for(int i=0; i<nx+1; i++) //check limits
		for(int j=0; j<ny+1; j++) {
			if (insideUBox(i,j))
				unew[i][j] = advectField(i, j, _dt, u, 0, uOffset);
			if (insideVBox(i,j))
				vnew[i][j] = advectField(i, j, _dt, v, 1, vOffset);
		}

	for(int i=0; i<nx+1; i++) //check limits
		for(int j=0; j<ny+1; j++) {
			if (insideUBox(i,j))
				u[i][j] = unew[i][j];
			if (insideVBox(i,j))
				v[i][j] = vnew[i][j];
		}
}
void Grid::applyGravity() {
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			if (isFluid(i,j)) //or (hasFluidNeighbour(i,j)))
				v[i][j] += g*dt;
}

//projection
int Grid::index(int i, int j) {
	return nx*j + i; //this is true for a grid with nx X ny elements (not for u and v for example).
}
void Grid::calcRHS() {
	double scaleRHS = 1.0/dx;

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++) {
			int idx = index(i, j);
			if (isFluid(i,j)) {
				rhs[idx] = -scaleRHS*(-u[i][j]-v[i][j]);
				if (insideUBox(i+1,j))
					rhs[idx] -= scaleRHS*u[i+1][j];
				if (insideVBox(i,j+1))
					rhs[idx] -= scaleRHS*v[i][j+1];

				if (insidePhiBox(i-1,j) and (isSolid[i-1][j]))
					rhs[idx] -= scaleRHS*(u[i][j]-0.0); //uSolid[i][j]
				if (insidePhiBox(i+1,j) and (isSolid[i+1][j]))
					rhs[idx] += scaleRHS*(u[i+1][j]-0.0); //uSolid[i+1][j]
				if (insidePhiBox(i,j-1) and (isSolid[i][j-1]))
					rhs[idx] -= scaleRHS*(v[i][j]-0.0); //vSolid[i][j]
				if (insidePhiBox(i,j+1) and (isSolid[i][j+1]))
					rhs[idx] += scaleRHS*(v[i][j+1]-0.0); //vSolid[i][j+1]
			}
			else
				rhs[idx] = 0.0;
		}
}
void Grid::calcAMatrix() {
	matA.zero();
	double scaleA = dt/(rho*dx*dx);

	for(int i=1; i<nx-1; i++)
		for(int j=1; j<ny-1; j++) {
			//if (isFluid(i,j)) {
				int idx = index(i, j);
				int diag = 4;
				int Aplusi = -1; //these are not Bridson's pseudocode variables
				int Aplusj = -1;
				int Aminusi = -1;
				int Aminusj = -1;

				if ((isSolid[i+1][j])) {
					Aplusi = 0;
					diag--;
				}
				if (not isFluid(i+1,j))
					Aplusi = 0;

				if (isSolid[i-1][j]) { //no entiendo por qué esto no tira segfault.
					Aminusi = 0;
					diag--;
				}
				if (not isFluid(i-1,j))
					Aminusi = 0;

				if (isSolid[i][j+1]) {
					Aplusj = 0;
					diag--;
				}
				if (not isFluid(i,j+1))
					Aplusj = 0;

				if (isSolid[i][j-1]) {
					Aminusj = 0;
					diag--;
				}
				if (not isFluid(i,j-1))
					Aminusj = 0;

				if ((idx-1 >= 0) and (Aminusi != 0))
					matA.set_element(idx-1, idx, scaleA*Aminusi); //sparseMatrix seems to lack element existence check
				if ((idx+1 <= nx*ny-1) and (Aplusi != 0))
					matA.set_element(idx+1, idx, scaleA*Aplusi);
				if ((idx-nx >= 0) and (Aminusj != 0))
					matA.set_element(idx-nx, idx, scaleA*Aminusj);
				if ((idx+nx <= nx*ny-1) and (Aplusj != 0))
					matA.set_element(idx+nx, idx, scaleA*Aplusj);
				if (diag != 0)
					matA.set_element(idx, idx, diag*scaleA);
			}
}
void Grid::pressureUpdate() {
	double scale = dt/(rho*dx);
	for(int i=0; i<nx+1; i++)
		for(int j=0; j<ny+1; j++)
			if (isFluid(i,j)) {
				if (insideUBox(i+1,j)) {
					if (isFluid(i+1,j))
						u[i+1][j] -= scale*(p[i+1][j] - p[i][j]);
					else
						u[i+1][j] -= scale*(p[i+1][j]*phi[i+1][j]/phi[i][j] - p[i][j]);
				}
				if (insideVBox(i,j+1)) {
					if (isFluid(i,j+1))
						v[i][j+1] -= scale*(p[i][j+1] - p[i][j]);
					else
						v[i][j+1] -= scale*(p[i][j+1]*phi[i][j+1]/phi[i][j] - p[i][j]);
				}
			}
}
void Grid::project() {
	double _residual;
	int _iterations;

	calcRHS();
	calcAMatrix();
	pcgSolver.solve(matA, rhs, vecP, _residual, _iterations); //I removed assert. solve does not converge inmediately after inserting a solid cell.
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			p[i][j] = vecP[index(i,j)];
	pressureUpdate();
}

//viscosity
void Grid::viscosity() {
	matA.zero();
	double scaleA = eta*dt/(rho*dx*dx);

	for(int i=0; i<nx+1; i++)
	for(int j=0; j<ny+1; j++) {
		//if (insidePhiBox(i,j) and (not isSolid[i][j])) {
		//if (insidePhiBox(i,j) and (isFluid(i,j))) {
		//if (insidePhiBox(i-1,j) and isFluid(i-1,j)) {
			int idx = index(i, j);
			int diag = 4;
			int Aplusi = -1; //these are not Bridson's pseudocode variables
			int Aplusj = -1;
			int Aminusi = -1;
			int Aminusj = -1;

			/*if (not isFluid(i+1,j))
				Aplusi = 0;
			if (not isFluid(i-1,j))
				Aminusi = 0;
			if (not isFluid(i,j+1))
				Aplusj = 0;
			if (not isFluid(i,j-1))
				Aminusj = 0;*/

			if ((idx-1 >= 0) and (Aminusi != 0))
				matA.set_element(idx-1, idx, scaleA*Aminusi); //sparseMatrix seems to lack element existence check
			if ((idx+1 <= nx*ny-1) and (Aplusi != 0))
				matA.set_element(idx+1, idx, scaleA*Aplusi);
			if ((idx-nx >= 0) and (Aminusj != 0))
				matA.set_element(idx-nx, idx, scaleA*Aminusj);
			if ((idx+nx <= nx*ny-1) and (Aplusj != 0))
				matA.set_element(idx+nx, idx, scaleA*Aplusj);
			if (diag != 0)
				matA.set_element(idx, idx, 1.0 + diag*scaleA);
		}

	double _residual;
	int _iterations;

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			if (not isSolid[i][j])
				rhs[index(i,j)] = u[i][j];
			else
				rhs[index(i,j)] = 0.0;
	pcgSolver.solve(matA, rhs, vecP, _residual, _iterations);
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			u[i][j] = vecP[index(i,j)];

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			if (not isSolid[i][j])
				rhs[index(i,j)] = v[i][j];
			else
				rhs[index(i,j)] = 0.0;
	pcgSolver.solve(matA, rhs, vecP, _residual, _iterations);
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			v[i][j] = vecP[index(i,j)];
}

//level set
void Grid::calcNeighbourMinDistance(int i, int j) {
	if (known[i][j] == false)
		for(int n = 0; n<4; n++) {
			int i2 = i + nbours4[n].i;
			int j2 = j + nbours4[n].j;
			if (insidePhiBox(i2,j2)) {//and (known[i2][j2])) { //i2,j2 is known neighbour
// 				known[i][j] = true;
				vec<double> r = cellCenter(i,j);
				double sqrDistance = (r-closestPoint[i2][j2]).norm2();
// 				if (sqrDistance < pow(phi[i2][j2],2)) //check if unknown distance is less than known distance
// 					known[i2][j2] = false;
				double distance = sqrt(sqrDistance);
 				if (distance < abs(phi[i][j])) { //check if unknown distance is the minimum of the neighbourhood
					phi[i][j] = distance*sign(phi[i2][j2]); //not sure about copying the sign from the neighbour
					closestPoint[i][j] = closestPoint[i2][j2]; //not sure about this
 				}
			}
		}
}
void Grid::initDistanceCalculation() {
	//first step: calc closestPoints and signed distances of cells near the surface. set as unknown the other cells.
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++) {
			known[i][j] = false;
			closestPoint[i][j].define(INFINITY, INFINITY);
		}

	for(int i=0; i<nx; i++) //deal with this later
		for(int j=0; j<ny; j++) {
			phinew[i][j] = phi[i][j];
		}

	for(int j=0; j<ny; j++)
		for(int i=0; i<nx; i++) {
			for(int n = 0; n<4; n++) {
				int i2 = i + nbours4[n].i;
				int j2 = j + nbours4[n].j;
				if (insidePhiBox(i2, j2) and (sign(phi[i][j]) != sign(phi[i2][j2]))) {
					double alpha = phi[i][j]/(phi[i][j] - phi[i2][j2]);
					vec<double> r = cellCenter(i,j);
					vec<double> r2 = cellCenter(i2,j2);
					vec<double> r0 = r*(1.0-alpha) + r2*alpha;
					double sqrDistance = (r0-r).norm2();
					if ((known[i][j] == false) or ((known[i][j] == true) and (sqrDistance <= pow(phinew[i][j],2)))) {
						phinew[i][j] = sign(phi[i][j])*sqrt(sqrDistance);
						closestPoint[i][j] = r0;
					}
					known[i][j] = true;
				}
			}
			if (not known[i][j])
				phinew[i][j] = INFINITY*sign(phi[i][j]); //ver signo
		}

	for(int i=0; i<nx; i++) //deal with this later
		for(int j=0; j<ny; j++) {
			phi[i][j] = phinew[i][j];
		}
}
void Grid::recalcSignedDistance(int N) {
	initDistanceCalculation();

	//fast sweeping method
	for(int n=0; n<N; n++) {
		for(int i=nx-1; i>=0; i--) //r-l d-u
			for(int j=0; j<ny; j++)
				//upwind(i,j);
				calcNeighbourMinDistance(i, j);

		for(int i=nx-1; i>=0; i--) //r-l u-d
			for(int j=ny-1; j>=0; j--)
				//upwind(i,j);
				calcNeighbourMinDistance(i, j);

		for(int i=0; i<nx; i++) //l-r d-u
			for(int j=0; j<ny; j++)
				//upwind(i,j);
				calcNeighbourMinDistance(i, j);

		for(int i=0; i<nx; i++) //l-r u-d
			for(int j=ny-1; j>=0; j--)
				//upwind(i,j);
				calcNeighbourMinDistance(i, j);
	}
}

//solid-air
bool Grid::isFluid(int i, int j) {
	return ((phi[i][j] <= 0.0) and (not isSolid[i][j]));
}
bool Grid::hasFluidNeighbour(int i, int j) {
	for(int n=0; n<8; n++) {
		int i2 = i + nbours8[n].i;
		int j2 = j + nbours8[n].j;
		if (insidePhiBox(i2,j2) and isFluid(i2,j2))
			return true;
	}
	return false;
}
void Grid::boxBoundaries() {
	for(int j=0; j<ny; j++) {
		u[0][j] = 0.0;
		u[nx][j] = 0.0;
	}
	for(int i=0; i<nx; i++) {
		v[i][0] = 0.0;
		v[i][ny] = 0.0;
	}
}
void Grid::solidBoundaries() {
	for(unsigned int n=0; n<solidCells.size(); n++) {
		int i = solidCells[n].i;
		int j = solidCells[n].j;

		if (insidePhiBox(i-1,j) and (not isSolid[i-1][j]))
			u[i][j] = min(u[i][j], 0.0);
		if (insidePhiBox(i+1,j) and (not isSolid[i+1][j]))
			u[i+1][j] = max(u[i+1][j], 0.0);
		if (insidePhiBox(i,j-1) and (not isSolid[i][j-1]))
			v[i][j] = min(v[i][j], 0.0);
		if (insidePhiBox(i,j+1) and (not isSolid[i][j+1]))
			v[i][j+1] = max(v[i][j+1], 0.0);

		p[i][j] = 0.0;
	}
}
void Grid::defineFluidCellsForPaint() {
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			if ((isFluid(i,j)) or (hasFluidNeighbour(i,j)))
				fluidCells.push_back({i, j});
}
void Grid::defineActiveCells() {
	fluidCells = vector<indexPair>(0);
	for(int i=0; i<nx+1; i++)
		for(int j=0; j<ny+1; j++) {
			activeU[i][j] = false;
			activeV[i][j] = false;
		}

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			//if (isFluid(i,j)) {
			if (phi[i][j] <= 0.0) { //this was isFluid, but i changed it
				if ((insideUBox(i+1,j)) and (isFluid(i+1,j)))
					activeU[i+1][j] = true;
				if ((insideVBox(i,j+1)) and (isFluid(i,j+1)))
					activeV[i][j+1] = true;
			}
}

//extrapolation
void Grid::clearNonFluidFaces() {
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			if (not isFluid(i,j)) {
				if ((insideUBox(i+1,j)) and (not isFluid(i+1,j)))
					u[i+1][j] = 0.0;
				if ((insideVBox(i,j+1)) and (not isFluid(i,j+1)))
					v[i][j+1] = 0.0;
			}
}
bool Grid::extrapolateU(int i, int j) {
	double ret = 0.0;
	int numberOfActiveNeighbours = 0;
	for(unsigned int n=0; n<8; n++) {
		int i2 = i + nbours8[n].i;
		int j2 = j + nbours8[n].j;
		if (insideUBox(i2,j2) and (activeU[i2][j2])) {
			ret += u[i2][j2];
			numberOfActiveNeighbours++;
		}
	}
	if (numberOfActiveNeighbours > 0) {
		unew[i][j] = ret/numberOfActiveNeighbours;
		return true;
	}
	return false;
}
bool Grid::extrapolateV(int i, int j) {
	double ret = 0.0;
	int numberOfActiveNeighbours = 0;
	for(unsigned int n=0; n<8; n++) {
		int i2 = i + nbours8[n].i;
		int j2 = j + nbours8[n].j;
		if (insideVBox(i2,j2) and (activeV[i2][j2])) {
			ret += v[i2][j2];
			numberOfActiveNeighbours++;
		}
	}
	if (numberOfActiveNeighbours > 0) {
		vnew[i][j] = ret/numberOfActiveNeighbours;
		return true;
	}
	return false;
}
void Grid::extrapolateCell(int i, int j) {
	if ((insideUBox(i,j)) and (not activeU[i][j]))
		activeUnew[i][j] = extrapolateU(i,j);
	if ((insideVBox(i,j)) and (not activeV[i][j]))
		activeVnew[i][j] = extrapolateV(i,j);
}
void Grid::extrapolateEasy(unsigned int N) {
	clearNonFluidFaces();

	for(unsigned int n=0; n<N; n++) {
		for(int i=0; i<nx+1; i++)
			for(int j=0; j<ny+1; j++) {
				if (insideUBox(i,j))
					unew[i][j] = u[i][j];
				if (insideVBox(i,j))
					vnew[i][j] = v[i][j];
			}

		for(int i=0; i<nx+1; i++)
			for(int j=0; j<ny+1; j++)
				extrapolateCell(i,j);

		for(int i=0; i<nx+1; i++)
			for(int j=0; j<ny+1; j++) {
				if (insideUBox(i,j))
					u[i][j] = unew[i][j];
				if (insideVBox(i,j))
					v[i][j] = vnew[i][j];
				activeU[i][j] = activeUnew[i][j]; //((activeU[i][j]) or (activeUnew[i][j]));
				activeV[i][j] = activeVnew[i][j]; //((activeV[i][j]) or (activeVnew[i][j]));
			}
	}
}

//checkers
bool Grid::insidePhiBox(int i, int j) {
	return ((i >= 0) and (j >= 0) and (i < nx) and (j < ny)); //check limits
}
bool Grid::insideVBox(int i, int j) {
	return ((i >= 0) and (j >= 0) and (i < nx) and (j <= ny)); //check limits
}
bool Grid::insideUBox(int i, int j) {
	return ((i >= 0) and (j >= 0) and (i <= nx) and (j < ny)); //check limits
}
bool Grid::CFLcondition() {
	return (dt < 5*dx/maxU());
}

//particles
const vec<double>& Grid::getParticlePosition(int n) {
	return particlePosition[n];
}
int Grid::getNumberOfParticles() {
	return numberOfParticles;
}
void Grid::particleAdvection() {
	for (int n = 0; n<numberOfParticles; n++) {
		int i = floor(particlePosition[n].x/dx); //check
		int j = floor(particlePosition[n].y/dx);
		if ((not insidePhiBox(i, j)) or (not isFluid(i,j)))
			particleInit(n);
		else
			particlePosition[n] += getU(i, j)*dt;
	}
}
void Grid::particleInit(int n) {
	do {
		double x = ((double)rand()/RAND_MAX)*(dx*nx);
		double y = ((double)rand()/RAND_MAX)*(dx*ny);
		int i = floor(x/dx);
		int j = floor(y/dx);
		if (isFluid(i,j)) {
			particlePosition[n] = {x, y};
			break;
		}
	} while (true);
}
void Grid::setParticlesVisible(bool b) {
	particlesAreVisible = b;
}

//accessors
double Grid::maxU() {
	double ret = 0.0;
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			ret = max(ret, sqrt(pow(u[i][j],2.0) + pow(v[i][j],2.0)));
	return ret;
}
double Grid::fluidArea() {
	double ret = 0.0;
	for(unsigned int n=0; n<fluidCells.size(); n++) {
		int i = fluidCells[n].i;
		int j = fluidCells[n].j;
		if (isFluid(i,j))
			ret += dx*dx;
	}
	return ret;
}
double Grid::divergence(int i, int j) {
	return (u[i+1][j]-u[i][j]+v[i][j+1]-v[i][j])/dx;
}
double Grid::energy() {
	double ret = 0.0;
	for(unsigned int n=0; n<fluidCells.size(); n++) {
		int i = fluidCells[n].i;
		int j = fluidCells[n].j;
		ret += pow(u[i][j], 2) + pow(v[i][j], 2);
	}
	return 0.5*ret;
}
vec<double> Grid::massCenter() {
	vec<double> ret(0.0, 0.0);
	int numberOfFluidCells = 0.0;
	for(unsigned int n=0; n<fluidCells.size(); n++) {
		int i = fluidCells[n].i;
		int j = fluidCells[n].j;
		if (isFluid(i,j)) {
			ret += cellCenter(i,j);
			numberOfFluidCells++;
		}
	}
	return ret/numberOfFluidCells;
}
vec<double> Grid::getU(int i, int j) {
	
	return (vec<double>(u[i][j], v[i][j]));
}
void Grid::setU(int i, int j, double _u, double _v) {
	u[i][j] = _u;
	v[i][j] = _v;
}
const vec<double> Grid::cellCenter(int i, int j) {
	return vec<double>(i+0.5, j+0.5)*dx;
}
bool Grid::cellIsKnown(int i, int j) {
	return known[i][j];
}
double Grid::getPhi(int i, int j) {
	return phi[i][j];
}
vec<double> Grid::getClosestPoint(int i, int j) {
	return closestPoint[i][j];
}
int Grid::getNumberOfFluidCells() {
	return fluidCells.size();
}
int Grid::getNumberOfSolidCells() {
	return solidCells.size();
}
const indexPair& Grid::getFluidCell(int n) {
	return fluidCells[n];
}
const indexPair& Grid::getSolidCell(int n) {
	return solidCells[n];
}
bool Grid::isSolidCell(int i, int j) {
	return isSolid[i][j];
}
void Grid::setSolidCell(int i, int j) {
	solidCells.push_back({i, j});
	isSolid[i][j] = true;
}

//prints
void Grid::printGlobal() {
	//cout << "Energy = " << energy() << endl;
	//cout << "Area = " << fluidArea() << endl;
	cout << timeStep*dt << " " << fluidArea() << " " << energy() << endl;
}
template<typename T> void Grid::printMatrix(vector<vector<T> >& q) {
	cout.precision(1);
	for(int j=ny-1; j>=0; j--) {
		for(int i=0; i<nx; i++) {
			cout << q[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

