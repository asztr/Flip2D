
#ifndef MOUSE_H
#define MOUSE_H

#include<QVector2D>
#include<QPoint>
#include<QSize>

class Mouse {
	private:
		QPoint oldPosition;
		QPoint newPosition;
		QPoint _pressedCell;
		QPointF _deltaPosition;

	public:
		void setPosition(const QPoint& _pos);
		QPoint& pressedCell(const QSize& _size, const QPointF& _cellPixelSize);
		QPointF& deltaPosition();
};

#endif