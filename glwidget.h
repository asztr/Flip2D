
#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QtGui>
#include <QGLWidget>
#include <QBrush>
#include <QFont>
#include <QPen>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QImage>
#include <QRect>

#include <cmath>
#include "Grid.h"
#include "Mouse.h"
#include "Common.h"

#define TOTAL_TIME	(0*DT)
#define TIME_PER_FRAME	(1*DT)
#define VIDEO_CAPTURE	false
#define FILL_FLUID	true
#define FLUID_BORDER	false

class QPainter;
class QPaintEvent;
class QWidget;

enum drawModes { testParticles = 1, vectorField, vectorFieldRenormalized, phiField, lastDrawMode };
const int numberOfDrawModes = lastDrawMode-testParticles;

class GLWidget : public QGLWidget {
	Q_OBJECT

	public:
		GLWidget(QWidget *parent);
		void mouseMoveEvent(QMouseEvent *event);
		void keyPressEvent(QKeyEvent *event);

	public slots:
		void animate();

	protected:
		void paintEvent(QPaintEvent *event);

	private:
		void paint(QPainter *painter, QPaintEvent *event);
		void paintCellBorder(QPainter *painter, int i, int j, QColor color);
		void paintClosestLines(QPainter *painter);
		void paintClosestPoints(QPainter *painter);
		void paintParticles(QPainter *painter);
		void paintFluid(QPainter *painter);
		void paintSolidCells(QPainter *painter);
		void paintCell(QPainter *painter, int i, int j, QColor color);
		void paintVectorField(QPainter *painter);
		void paintVectorModulus(QPainter *painter);
		void paintPhiField(QPainter *painter);
		void paintDivergence(QPainter *painter);
		void paintMassCenter(QPainter *painter);
		void paintGrid(QPainter *painter);
		void drawArrow(QPainter *painter, QPoint p0, QPointF scale, double angle);

		void captureScreen();
		void defineCellPixelSize();
		QPoint phys2Window(const QPointF& physCoord);
		QPoint cell2Window(int i, int j);
		QPointF cellPhysCenter(int i, int j);
		QPointF cellPhysCorner(int i, int j);

		QPointF cellPixelSize;
		drawModes drawMode;
		double toPixels;
		double maxU;
		double t, timePerFrame;

		bool videoCapture;
		bool pause;
		bool gridVisible;
		bool massCenterVisible;
		bool closestPointsVisible;
		bool closestLinesVisible;

		Mouse mouse;
		Grid grid;

		//colors, brushes and pens
		static const QColor fluidborderColor;
		static const QColor fluidColor;
		static const QColor solidColor;
		static const QColor airColor;
		static const QColor particleColor;
		static const QColor arrowColor;
		static const QColor gridColor;
		static const QBrush backgroundBrush;
		static const QBrush fluidBrush;
		QPen arrowPen;
		QPen particlePen;
		QPen fluidborderPen;
		QPen gridPen;

		//level set
		QPointF interPoint(int i, int j, int i2, int j2);
		void drawPhysLine(QPointF p1, QPointF p2, QPainter *painter);
		void paintFluidSquare(int i, int j, QPainter *painter);
		void fill1(int i, int j, QPainter *painter);
		void fill2(int i, int j, QPainter *painter);
		void fill3(int i, int j, QPainter *painter);
		void fill4(int i, int j, QPainter *painter);
		void fill5(int i, int j, QPainter *painter);
		void fill6(int i, int j, QPainter *painter);
		void fill7(int i, int j, QPainter *painter);
		void fill8(int i, int j, QPainter *painter);
		void fill9(int i, int j, QPainter *painter);
		void fill10(int i, int j, QPainter *painter);
		void fill11(int i, int j, QPainter *painter);
		void fill12(int i, int j, QPainter *painter);
		void fill13(int i, int j, QPainter *painter);
		void fill14(int i, int j, QPainter *painter);
		void fill15(int i, int j, QPainter *painter);

		void line7(int i, int j, QPainter *painter);
		void line4(int i, int j, QPainter *painter);
		void line2(int i, int j, QPainter *painter);
		void line1(int i, int j, QPainter *painter);
		void lineV(int i, int j, QPainter *painter);
		void lineH(int i, int j, QPainter *painter);
};

#endif
