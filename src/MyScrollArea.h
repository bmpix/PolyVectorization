#pragma once

#include <QtWidgets>

class MyScrollArea : public QScrollArea
{
protected:
	void wheelEvent(QWheelEvent * event);
};
