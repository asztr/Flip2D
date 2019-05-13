#include "glwidget.h"

const double toDegrees = 180.0 / (4 * atan(1.0));

const QColor GLWidget::fluidColor = QColor(0, 85, 255); //blue fluid
const QColor GLWidget::fluidborderColor = QColor(0, 0, 0); //white level set surface
const QColor GLWidget::solidColor = QColor(0, 0, 0); //black solids
const QColor GLWidget::airColor = QColor(255, 255, 255); //white air
const QColor GLWidget::particleColor = QColor(200, 200, 200); //gray particles
const QColor GLWidget::arrowColor = QColor(206, 51, 53); //dark red arrows
const QColor GLWidget::gridColor = QColor(225, 225, 225); //light gray grid
const QBrush GLWidget::backgroundBrush = QBrush(airColor); //white air background
const QBrush GLWidget::fluidBrush = QBrush(fluidColor, Qt::SolidPattern);

GLWidget::GLWidget(QWidget *parent) : QGLWidget(QGLFormat(QGL::SampleBuffers), parent) {
	//grid initialization
	timePerFrame = TIME_PER_FRAME;
	drawMode = testParticles;
	maxU = grid.maxU();

	//booleans
	videoCapture = VIDEO_CAPTURE;
	pause = false;
	gridVisible = true;
	massCenterVisible = false;
	closestPointsVisible = false;
	closestLinesVisible = false;

	//visualization settings (hardcoded for maximized window on ubuntu 1680x1050)
	//setMinimumSize(1657, 950);
	setMinimumSize(1250, 1250);
	setAutoFillBackground(false);
	grabKeyboard();

	arrowPen = QPen(arrowColor);
	particlePen = QPen(particleColor);
	fluidborderPen = QPen(fluidborderColor);
	gridPen = QPen(gridColor);
	arrowPen.setWidth(1);
	particlePen.setWidth(3);
	fluidborderPen.setWidth(6);
	gridPen.setWidth(1);
}

void GLWidget::captureScreen() {
	grabFrameBuffer(false).save(QString("capture/screenshot%1.bmp").arg(grid.getTimeStep()+1000), "bmp", 100);
}

void GLWidget::defineCellPixelSize() {
	cellPixelSize.setX((double)size().width()/grid.getNX());
	cellPixelSize.setY((double)size().height()/grid.getNY());
}
QPoint GLWidget::phys2Window(const QPointF& physCoord) {
	double dx = grid.getDX();
	int x = floor(size().width()*physCoord.x()/(grid.getNX()*dx)); //round(physCoord.x()*size().width()/(grid.getNX()*dx));
	int y = floor(size().height()*(1.0 - physCoord.y()/(grid.getNY()*dx)) );
	return QPoint(x, y);
}
QPointF GLWidget::cellPhysCenter(int i, int j) {
	return QPointF(i+0.5, j+0.5)*grid.getDX();
}
QPointF GLWidget::cellPhysCorner(int i, int j) {
	return QPointF(i,j)*grid.getDX();
}
QPoint GLWidget::cell2Window(int i, int j) {
	return phys2Window(cellPhysCorner(i,j));
}

void GLWidget::paint(QPainter *painter, QPaintEvent *event) {
	painter->fillRect(event->rect(), backgroundBrush);
	defineCellPixelSize();
	if ((gridVisible) and (drawMode != phiField)) paintGrid(painter);
	paintFluid(painter);
	paintSolidCells(painter);
	//paintFluidCells(painter);
	switch(drawMode) {
		case vectorField:
			//paintVectorField(painter);
			paintVectorModulus(painter);
			break;
		case vectorFieldRenormalized:
			if (drawMode == vectorFieldRenormalized)
				maxU = grid.maxU();
			paintVectorField(painter);
			break;
		case testParticles:
			paintParticles(painter);
			break;
		case phiField:
			paintPhiField(painter);
			if (gridVisible) paintGrid(painter);
			//paintDivergence(painter);
			break;
		default:
			break;
	}
	if (closestLinesVisible) paintClosestLines(painter);
	if (closestPointsVisible) paintClosestPoints(painter);
	if (massCenterVisible) paintMassCenter(painter);
}
void GLWidget::paintCellBorder(QPainter *painter, int i, int j, QColor color) {
	painter->setPen(color);
	QPoint corner = cell2Window(i,j);
	painter->drawRect(corner.x(), (int)(corner.y()-cellPixelSize.y()), (int)cellPixelSize.x(), (int)cellPixelSize.y());
}
void GLWidget::paintClosestLines(QPainter *painter) {
	painter->setPen(arrowPen);
	for (int i = 0; i<grid.getNX(); i++)
		for (int j = 0; j<grid.getNY(); j++) {
			//if (grid.cellIsKnown(i,j)) {
				QPointF physPos(grid.getClosestPoint(i,j).x, grid.getClosestPoint(i,j).y);
				painter->drawLine(phys2Window(cellPhysCenter(i,j)), phys2Window(physPos));
			}
}
void GLWidget::paintClosestPoints(QPainter *painter) {
	painter->setPen(arrowPen);
	for (int i = 0; i<grid.getNX(); i++)
		for (int j = 0; j<grid.getNY(); j++) {
			//if (grid.cellIsKnown(i,j)) {
				QPointF physPos(grid.getClosestPoint(i,j).x, grid.getClosestPoint(i,j).y);
				painter->drawEllipse(phys2Window(physPos), 1, 1);
			}
}
void GLWidget::paintParticles(QPainter *painter) {
	painter->setPen(particlePen);
	for (int n=0; n<grid.getNumberOfParticles(); n++) {
		QPointF _physPos(grid.getParticlePosition(n).x, grid.getParticlePosition(n).y); //in physical units
		painter->drawPoint(phys2Window(_physPos));
	}
}
void GLWidget::paintVectorField(QPainter *painter) {
	//toPixels = cellPixelSize.length()/(sqrt(2.0)*maxU);
	toPixels = 0.5*cellPixelSize.manhattanLength()/(sqrt(2.0)*maxU); //check length
	painter->setPen(arrowPen);
	for (int i = 0; i<grid.getNX(); i++)
		for (int j = 0; j<grid.getNY(); j++) {
			QPointF scale = grid.getU(i,j).norm()*toPixels*QPointF(1./10, 1./12);
			//if (scale.manhattanLength() > .5) {
			double angle = grid.getU(i,j).angle()*toDegrees;
			QPoint _windowPos(phys2Window(cellPhysCenter(i, j)));
			drawArrow(painter, _windowPos, scale, angle);
			//}
		}
}
void GLWidget::paintVectorModulus(QPainter *painter) {
	//toPixels = 0.5*(pow(cellPixelSize.x(),2) + pow(cellPixelSize.y(),2))/(sqrt(2.0)*maxU); //check length
	for (int i = 0; i<grid.getNX(); i++)
		for (int j = 0; j<grid.getNY(); j++) {
			double scale = grid.getU(i,j).norm()*200;
			if (scale > 255)
				scale = 255;
			QPoint _windowPos(phys2Window(cellPhysCenter(i, j)));
			paintCell(painter, i, j, QColor(scale, scale, scale));
		}

}
void GLWidget::paintFluid(QPainter *painter) {
	for(int n=0; n<grid.getNumberOfFluidCells(); n++) {
		indexPair idx = grid.getFluidCell(n);
		paintFluidSquare(idx.i, idx.j, painter);
	}
}
void GLWidget::paintPhiField(QPainter *painter) {
	for (int i = 0; i<grid.getNX(); i++)
		for (int j = 0; j<grid.getNY(); j++)
			if (not grid.isSolidCell(i,j)) {
				double _phi = 30*grid.getPhi(i,j);
				if (abs(_phi) < 255) {
					//if (grid.cellIsKnown(i,j) == true) {
						if (_phi < 0.0)
							paintCell(painter, i, j, QColor(100-_phi, 0.0, 0.0)); //celdas conocidas de fluido (rojo)
						else if (_phi > 0.0) 
							paintCell(painter, i, j, QColor(_phi, _phi, _phi)); //celdas conocidas de aire (gris)
						else
							paintCell(painter, i, j, QColor(0.0, 0.0, 200)); //celdas con phi = 0
					//}
					//else
						//paintCell(painter, i, j, QColor(0.0, 200, 0.0)); //celdas desconocidas (verde)
				}
			}
}
void GLWidget::paintSolidCells(QPainter *painter) {
	for(int n=0; n<grid.getNumberOfSolidCells(); n++) {
		indexPair idx = grid.getSolidCell(n);
		paintCell(painter, idx.i, idx.j, solidColor);
	}
}
void GLWidget::paintCell(QPainter *painter, int i, int j, QColor color) {
	painter->fillRect(QRect(cell2Window(i, j), cell2Window(i+1, j+1) -QPoint(0,1)), color);
}
void GLWidget::paintMassCenter(QPainter *painter) {
	QPen massCenterPen(Qt::black, 2);
	QBrush massCenterBrush(Qt::green);
	painter->setPen(massCenterPen);
	painter->setBrush(massCenterBrush);
	vec<double> cm = grid.massCenter();
	QPointF pnt(cm.x, cm.y);
	painter->drawEllipse(phys2Window(pnt), (int)cellPixelSize.x(), (int)cellPixelSize.y()); 
}
void GLWidget::paintDivergence(QPainter *painter) {
	for (int i = 0; i<grid.getNX(); i++)
		for (int j = 0; j<grid.getNY(); j++)
			if (not grid.isSolidCell(i,j)) {
				double div = abs(grid.divergence(i,j));
				paintCell(painter, i, j, QColor(div, div, div));
			}
}
void GLWidget::drawArrow(QPainter *painter, QPoint p0, QPointF scale, double angle) {
	painter->save();
	painter->translate(p0);
	painter->rotate(-angle);
	painter->scale(scale.x(), scale.y());
	painter->drawLine(0, 0, 10, 0);
	painter->drawLine(10, 0, 7, -3);
	painter->drawLine(10, 0, 7, 3);
	painter->restore();
}
void GLWidget::paintGrid(QPainter *painter) {
	painter->setPen(gridPen);
	for (int i = 0; i<grid.getNX(); i++)
		painter->drawLine(i*cellPixelSize.x(), 0, i*cellPixelSize.x(), size().height());
	for (int j = 0; j<grid.getNY(); j++)
		painter->drawLine(0, j*cellPixelSize.y(), size().width(), j*cellPixelSize.y());
}

void GLWidget::animate() {
	/*if (drawMode == testParticles)
		grid.setParticlesVisible(true);
	else
		grid.setParticlesVisible(false);*/
	
		double dt = grid.getDT();
		int it = grid.getTimeStep();
		double t = it*dt;

		if ((pause) or (it % (int)(timePerFrame/dt) == 0)) {
			repaint();
			if (videoCapture == true)
				captureScreen();
		}
		//if ((t >= TOTAL_TIME) and (TOTAL_TIME > 0.0))
			//exit(0);
	if (not pause)
		grid.update();
}
void GLWidget::mouseMoveEvent(QMouseEvent *event) {
	if (event->buttons() != Qt::NoButton) {
		mouse.setPosition(event->pos());
		QPoint pressedCell = mouse.pressedCell(size(), cellPixelSize);
		if (event->buttons() == Qt::LeftButton) {
			QPointF delta = mouse.deltaPosition(); //in pixels
			if (not delta.isNull())
				grid.setU(pressedCell.x(), pressedCell.y(), delta.x(), delta.y());
		}
		if (event->buttons() == Qt::RightButton)
			grid.setSolidCell(pressedCell.x(), pressedCell.y());
	}
}
void GLWidget::keyPressEvent(QKeyEvent *event) {
	if (event->key() == Qt::Key_L)
		closestLinesVisible = !closestLinesVisible;
	if (event->key() == Qt::Key_C)
		closestPointsVisible = !closestPointsVisible;
	if (event->key() == Qt::Key_M)
		massCenterVisible = !massCenterVisible;
	if (event->key() == Qt::Key_G)
		gridVisible = !gridVisible;
	if (event->key() == Qt::Key_P)
		pause = !pause;
	if (event->key() == Qt::Key_R)
		grid.resetGrid();
	if (event->key() == Qt::Key_D)
		drawMode = (drawModes)switchInteger(drawMode, 1, numberOfDrawModes); //drawMode cycles through 1, 2, 3
}
void GLWidget::paintEvent(QPaintEvent *event) {
	QPainter painter;
	painter.begin(this);
	painter.setRenderHint(QPainter::Antialiasing);
	paint(&painter, event);
	painter.end();
}
