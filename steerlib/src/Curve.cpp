//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
// Copyright (c) 2015 Mahyar Khayatkhoei
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
	{
		std::cerr << "Error: DrawCurve does not have enough points" << std::endl;
		return;
	}

	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	Point EP;
	for (unsigned int x = 0; x < controlPoints.size()-1; x++)
	{
		for (float CurrTime = controlPoints[x].time; CurrTime <= controlPoints[x + 1].time; CurrTime += (float)window)
		{
			if (type == hermiteCurve)
			{
				EP = useHermiteCurve(x + 1, CurrTime);
			}
			else if (type == catmullCurve)
			{
				EP = useCatmullCurve(x + 1, CurrTime);
			}

			DrawLib::drawLine(controlPoints[x].position, EP, curveColor, curveThickness);
			controlPoints[x].position = EP;
		}
	}
	return;
#endif
}

bool compare(CurvePoint x, CurvePoint y) {
	return x.time < y.time;
}
// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	std::sort(Curve::controlPoints.begin(), Curve::controlPoints.end(), compare);
	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
// Note that this function should return false if the end of the curve is reached, or no next point can be found
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	// Note that nextPoint is an integer containing the index of the next control point
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve given the next control point (nextPoint)
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
	if (controlPoints.size() < 2)
		return false;
	return true;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	if (controlPoints[0].time == time) {
		nextPoint = 0;
		return true;
	}
	for (int x = 0;x < controlPoints.size();x++) {
		if ((controlPoints[x].time < time) && (controlPoints[x + 1].time >= time)) {
			nextPoint = x + 1;
			return true;
		}
	}
	return false;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	float intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	float normalTime = (time - controlPoints[nextPoint - 1].time) / (intervalTime);

	// Calculate position at t = time on Hermite curve
	Point a = (2 * pow(normalTime, 3) - 3 * pow(normalTime, 2) + 1) * controlPoints[nextPoint - 1].position;
	Vector b = (pow(normalTime, 3) - 2 * pow(normalTime, 2) + normalTime) * controlPoints[nextPoint - 1].tangent * intervalTime;
	Point c = (-2 * pow(normalTime, 3) + 3 * pow(normalTime, 2)) * controlPoints[nextPoint].position;
	Vector d = (pow(normalTime, 3) - pow(normalTime, 2)) * controlPoints[nextPoint].tangent * intervalTime;
	// Return result
	return a + b + c + d;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float intervalTime, normalTime;
	Vector s;
	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	normalTime = (time - controlPoints[nextPoint - 1].time) / (intervalTime);

	// For start point
	if (nextPoint - 1 == 0)
	{
		s = ((controlPoints[nextPoint + 1].time - controlPoints[nextPoint - 1].time) / (controlPoints[nextPoint + 1].time - controlPoints[nextPoint].time) *
			(controlPoints[nextPoint].position - controlPoints[nextPoint - 1].position) / (intervalTime)) -
			((intervalTime) / (controlPoints[nextPoint + 1].time - controlPoints[nextPoint].time) *
			(controlPoints[nextPoint + 1].position - controlPoints[nextPoint - 1].position) / (controlPoints[nextPoint + 1].time - controlPoints[nextPoint - 1].time));
	}
	// Middle Points
	else if (nextPoint != controlPoints.size() - 1)
	{
		s = ((intervalTime) / (controlPoints[nextPoint + 1].time - controlPoints[nextPoint - 1].time) *
			(controlPoints[nextPoint + 1].position - controlPoints[nextPoint].position) / (controlPoints[nextPoint + 1].time - controlPoints[nextPoint].time)) +
			((controlPoints[nextPoint + 1].time - controlPoints[nextPoint].time) / (controlPoints[nextPoint + 1].time - controlPoints[nextPoint - 1].time) *
			(controlPoints[nextPoint].position - controlPoints[nextPoint - 1].position) / (intervalTime));;
	}
	// End point
	else
	{
		s = ((controlPoints[nextPoint].time - controlPoints[nextPoint - 2].time) / (controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time) *
			(controlPoints[nextPoint].position - controlPoints[nextPoint - 1].position) / (controlPoints[nextPoint - 1].time - controlPoints[nextPoint - 2].time)) -
			((controlPoints[nextPoint - 1].time - controlPoints[nextPoint - 2].time) / (controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time) *
			(controlPoints[nextPoint].position - controlPoints[nextPoint - 2].position) / (controlPoints[nextPoint].time - controlPoints[nextPoint - 2].time));
	}

	newPosition = controlPoints[nextPoint - 1].position + (s*normalTime);
	return newPosition;
}