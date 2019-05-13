#include "Common.h"

int sign(double r) {
	if (r > 0.0)
		return 1;
	else
		return -1;
}

int switchInteger(int value, int firstValue, int lastValue) {
	return (firstValue + (value+lastValue) % lastValue); //(1 + ((int)drawMode+3) % 3);
}
