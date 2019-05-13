
#include "glwidget.h"

QPointF GLWidget::interPoint(int i, int j, int i2, int j2) {
	double phi1 = grid.getPhi(i,j);
	double phi2 = grid.getPhi(i2,j2);
	double alpha = phi1/(phi1-phi2);
	return cellPhysCenter(i,j)*(1.0-alpha) + cellPhysCenter(i2,j2)*alpha;
}
void GLWidget::drawPhysLine(QPointF p1, QPointF p2, QPainter *painter) {
	painter->drawLine(phys2Window(p1), phys2Window(p2));
}
void GLWidget::fill1(int i, int j, QPainter *painter) {
	const QPointF points[3] = { phys2Window(cellPhysCenter(i,j)),
								phys2Window(interPoint(i,j,i,j+1)),
								phys2Window(interPoint(i,j,i+1,j)) };
	painter->drawConvexPolygon(points, 3);
}
void GLWidget::fill2(int i, int j, QPainter *painter) {
	const QPointF points[3] = { phys2Window(cellPhysCenter(i+1,j)),
								phys2Window(interPoint(i+1,j,i+1,j+1)),
								phys2Window(interPoint(i,j,i+1,j)) };
	painter->drawConvexPolygon(points, 3);
}
void GLWidget::fill3(int i, int j, QPainter *painter) {
	const QPointF points[4] = { phys2Window(cellPhysCenter(i,j)),
								phys2Window(interPoint(i,j,i,j+1)),
								phys2Window(interPoint(i+1,j,i+1,j+1)),
								phys2Window(cellPhysCenter(i+1,j))};
	painter->drawConvexPolygon(points, 4);
}
void GLWidget::fill4(int i, int j, QPainter *painter) {
	const QPointF points[3] = { phys2Window(cellPhysCenter(i+1,j+1)),
								phys2Window(interPoint(i+1,j+1,i,j+1)),
								phys2Window(interPoint(i+1,j+1,i+1,j)) };
	painter->drawConvexPolygon(points, 3);
}
void GLWidget::fill5(int i, int j, QPainter *painter) {
	const QPointF points[6] = { phys2Window(cellPhysCenter(i,j)),
								phys2Window(interPoint(i,j,i,j+1)),
								phys2Window(interPoint(i,j+1,i+1,j+1)),
								phys2Window(cellPhysCenter(i+1,j+1)),
								phys2Window(interPoint(i+1,j+1,i+1,j)),
								phys2Window(interPoint(i,j,i+1,j)) };
	painter->drawConvexPolygon(points, 6);
}
void GLWidget::fill6(int i, int j, QPainter *painter) {
	const QPointF points[4] = { phys2Window(cellPhysCenter(i+1,j)),
								phys2Window(interPoint(i+1,j,i,j)),
								phys2Window(interPoint(i+1,j+1,i,j+1)),
								phys2Window(cellPhysCenter(i+1,j+1))};
	painter->drawConvexPolygon(points, 4);
}
void GLWidget::fill7(int i, int j, QPainter *painter) {
	const QPointF points[5] = { phys2Window(cellPhysCenter(i,j)),
								phys2Window(cellPhysCenter(i+1,j)),
								phys2Window(cellPhysCenter(i+1,j+1)),
								phys2Window(interPoint(i+1,j+1,i,j+1)),
								phys2Window(interPoint(i,j,i,j+1)) };
	painter->drawConvexPolygon(points, 5);
}
void GLWidget::fill8(int i, int j, QPainter *painter) {
	const QPointF points[3] = { phys2Window(cellPhysCenter(i,j+1)),
								phys2Window(interPoint(i,j+1,i+1,j+1)),
								phys2Window(interPoint(i,j+1,i,j)) };
	painter->drawConvexPolygon(points, 3);
}
void GLWidget::fill9(int i, int j, QPainter *painter) {
	const QPointF points[4] = { phys2Window(cellPhysCenter(i,j)),
								phys2Window(cellPhysCenter(i,j+1)),
								phys2Window(interPoint(i,j+1,i+1,j+1)),
								phys2Window(interPoint(i,j,i+1,j)) };
	painter->drawConvexPolygon(points, 4);
}
void GLWidget::fill10(int i, int j, QPainter *painter) {
	const QPointF points[6] = { phys2Window(cellPhysCenter(i,j+1)),
								phys2Window(interPoint(i,j+1,i+1,j+1)),
								phys2Window(interPoint(i+1,j+1,i+1,j)),
								phys2Window(cellPhysCenter(i+1,j)),
								phys2Window(interPoint(i,j,i+1,j)),
								phys2Window(interPoint(i,j,i,j+1)) };
	painter->drawConvexPolygon(points, 6);
}
void GLWidget::fill11(int i, int j, QPainter *painter) {
	const QPointF points[5] = { phys2Window(cellPhysCenter(i,j)),
								phys2Window(cellPhysCenter(i,j+1)),
								phys2Window(interPoint(i,j+1,i+1,j+1)),
								phys2Window(interPoint(i+1,j+1,i+1,j)),
								phys2Window(cellPhysCenter(i+1,j)) };
	painter->drawConvexPolygon(points, 5);
}
void GLWidget::fill12(int i, int j, QPainter *painter) {
	const QPointF points[4] = { phys2Window(cellPhysCenter(i,j+1)),
								phys2Window(cellPhysCenter(i+1,j+1)),
								phys2Window(interPoint(i+1,j+1,i+1,j)),
								phys2Window(interPoint(i,j,i,j+1)) };
	painter->drawConvexPolygon(points, 4);
}
void GLWidget::fill13(int i, int j, QPainter *painter) {
	const QPointF points[5] = { phys2Window(cellPhysCenter(i,j)),
								phys2Window(cellPhysCenter(i,j+1)),
								phys2Window(cellPhysCenter(i+1,j+1)),
								phys2Window(interPoint(i+1,j+1,i+1,j)),
								phys2Window(interPoint(i+1,j,i,j)) };
	painter->drawConvexPolygon(points, 5);
}
void GLWidget::fill14(int i, int j, QPainter *painter) {
	const QPointF points[5] = { phys2Window(cellPhysCenter(i,j+1)),
								phys2Window(cellPhysCenter(i+1,j+1)),
								phys2Window(cellPhysCenter(i+1,j)),
								phys2Window(interPoint(i,j,i+1,j)),
								phys2Window(interPoint(i,j,i,j+1)) };
	painter->drawConvexPolygon(points, 5);
}
void GLWidget::fill15(int i, int j, QPainter *painter) {
	painter->fillRect(QRect(phys2Window(cellPhysCenter(i,j)), phys2Window(cellPhysCenter(i+1,j+1))-QPoint(0,1)), fluidColor);
}

void GLWidget::line7(int i, int j, QPainter *painter) { drawPhysLine(interPoint(i,j,i,j+1), interPoint(i,j+1,i+1,j+1), painter); }
void GLWidget::line4(int i, int j, QPainter *painter) { drawPhysLine(interPoint(i,j+1,i+1,j+1), interPoint(i+1,j+1,i+1,j), painter); }
void GLWidget::line2(int i, int j, QPainter *painter) { drawPhysLine(interPoint(i+1,j+1,i+1,j), interPoint(i+1,j,i,j), painter); }
void GLWidget::line1(int i, int j, QPainter *painter) { drawPhysLine(interPoint(i+1,j,i,j), interPoint(i,j,i,j+1), painter); }
void GLWidget::lineV(int i, int j, QPainter *painter) { drawPhysLine(interPoint(i,j+1,i+1,j+1), interPoint(i,j,i+1,j), painter); }
void GLWidget::lineH(int i, int j, QPainter *painter) { drawPhysLine(interPoint(i,j,i,j+1), interPoint(i+1,j,i+1,j+1), painter); }

void GLWidget::paintFluidSquare(int i, int j, QPainter *painter) {
	double d8 = grid.getPhi(i, j+1);
	double d4 = grid.getPhi(i+1, j+1);
	double d2 = grid.getPhi(i+1, j);
	double d1 = grid.getPhi(i, j);
	int caseIndex = 8*(d8<0) + 4*(d4<0) + 2*(d2<0) + 1*(d1<0);

	if (FILL_FLUID == true) {
		painter->setBrush(fluidBrush);
		painter->setPen(QPen(fluidColor));
		switch(caseIndex) {
			case 0: break;
			case 1: fill1(i,j,painter); break;
			case 2: fill2(i,j,painter); break;
			case 3: fill3(i,j,painter); break;
			case 4: fill4(i,j,painter); break; 
			case 5: fill5(i,j,painter); break;
			case 6: fill6(i,j,painter); break;
			case 7: fill7(i,j,painter); break;
			case 8: fill8(i,j,painter); break;
			case 9: fill9(i,j,painter); break;
			case 10: fill10(i,j,painter); break; 
			case 11: fill11(i,j,painter); break;
			case 12: fill12(i,j,painter); break;
			case 13: fill13(i,j,painter); break;
			case 14: fill14(i,j,painter); break;
			case 15: fill15(i,j,painter); break;
		}
	}

	if (FLUID_BORDER == true) {
		painter->setPen(fluidborderPen);
		switch(caseIndex) {
			case 0: break;
			case 1: line1(i,j,painter); break;
			case 2: line2(i,j,painter); break;
			case 3: lineH(i,j,painter); break;
			case 4: line4(i,j,painter); break;
			case 5: line7(i,j,painter); line2(i,j,painter); break;
			case 6: lineV(i,j,painter); break;
			case 7: line7(i,j,painter); break;
			case 8: line7(i,j,painter); break;
			case 9: lineV(i,j,painter); break;
			case 10: line4(i,j,painter); line1(i,j,painter); break;
			case 11: line4(i,j,painter); break;
			case 12: lineH(i,j,painter); break;
			case 13: line2(i,j,painter); break;
			case 14: line1(i,j,painter); break;
			case 15: break;
		}
	}
}
