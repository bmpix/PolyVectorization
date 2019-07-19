#include "stdafx.h"
#include "MyScrollArea.h"

void MyScrollArea::wheelEvent(QWheelEvent * event)
{
	auto km = QApplication::keyboardModifiers();
	if (!(km & Qt::KeyboardModifier::AltModifier))
		QScrollArea::wheelEvent(event);
}
