#include "CircularSegment.h"
#include <assert.h>
#include <math.h>
CircularSegment::CircularSegment (double begin, double end, double modulus)
	:begin_(begin),end_(end),modulus_(modulus)
{

}

std::vector<CircularSegment> CircularSegment::unionWith (const CircularSegment& other, double tolerance) const
{
	CircularSegment intersection;
	if (!intersectionWith(other,intersection, tolerance))
	{
		std::vector<CircularSegment> result;
		result.push_back(*this);
		result.push_back(other);
		return result;
	}

	std::vector<CircularSegment> result;
	result.push_back(*this);
	result[0].extendToContainPoint(other.begin());
	result[0].extendToContainPoint(other.end());
	return result;
}

double CircularSegment::pointAtT(double t) const
{
	if ((modulus_ < 0) || (end_+1e-10 > begin_))
		return begin_*(1-t)+end_*t;
	else
	{
		double result = begin_*(1-t) + t*(end_ + modulus_);
		if (result > modulus_)
			result -= modulus_;
		return result;
	}
}

std::vector<CircularSegment> CircularSegment::complement() const
{
	std::vector<CircularSegment> result;
	CircularSegment s1 (0,begin_,modulus_), s2(end_,modulus_,modulus_);
	result.push_back(s1);
	result.push_back(s2);
	return result;
}

bool CircularSegment::intersectionWith (const CircularSegment& other, CircularSegment& intersection, double tolerance) const
{
	assert(fabs(other.modulus()- modulus())<1e-10);
	
	if (isInside(other.begin(),tolerance))
	{
		if (isInside(other.end(),tolerance))
		{
			intersection = other;
			return true;
		}

		intersection = CircularSegment(other.begin(),end(),modulus());
		return true;
	}
	
	if (isInside(other.end(),tolerance))
	{
		intersection = CircularSegment(begin(),other.end(),modulus());
		return true;
	}

	if (other.isInside(end(),tolerance))
	{
		if (other.isInside(begin(),tolerance))
		{
			intersection = *this;
			return true;
		}
		intersection = CircularSegment(other.begin(),end(),modulus());
		return true;
	}

	if (other.isInside(begin(),tolerance))
	{
		intersection = CircularSegment(begin(),other.end(),modulus());
		return true;
	}
	return false;
}

void CircularSegment::extendToContainPoint (double p)
{
	if (isInside (p))
		return;

	if (modulus_ > 0)
	{
		//update left or right
		CircularSegment leftExtension (p,end(),modulus());
		CircularSegment rightExtension (begin(),p,modulus());
		if (leftExtension.length()<rightExtension.length())
		{
			*this = leftExtension;
		}
		else
		{
			*this = rightExtension;
		}
	}
	else
	{
		if (p < begin_)
			*this = CircularSegment(p, end_);
		else if (p > end_)
			*this = CircularSegment(begin_, p);
	}
}

bool CircularSegment::isInside (double p, double tolerance) const
{
	if (modulus_>0)
	{
		double diff = fabs(CircularSegment(begin_,p,modulus_).length() + CircularSegment(p,end_,modulus_).length() - length());
		return (diff < tolerance);
	}
	else
	{
		return ((begin_ - tolerance < p) && (end_ + tolerance > p));
	}
}

double CircularSegment::length() const
{
	double s = end_ - begin_;
	if (s < 0)
	{
		if (modulus_ > 0)
			s+=modulus_;
		else
			return 1e100;
	}
	return s;
}