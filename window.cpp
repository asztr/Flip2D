#include <QtGui>
#include "glwidget.h"
#include "window.h"

Window::Window() : QWidget() {
	GLWidget *openGL = new GLWidget(this);

	QGridLayout *layout = new QGridLayout;
	layout->addWidget(openGL, 0, 1);
	setLayout(layout);

	QTimer *timer = new QTimer(this);
	connect(timer, SIGNAL(timeout()), openGL, SLOT(animate()));
	timer->start(50);

	setWindowTitle(tr("Split Navier-Stokes"));
}