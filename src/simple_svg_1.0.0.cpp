#include "stdafx.h"
#include "simple_svg_1.0.0.hpp"
namespace svg
{
	std::string elemStart(std::string const & element_name)
	{
		return "\t<" + element_name + " ";
	}
	std::string elemEnd(std::string const & element_name)
	{
		return "</" + element_name + ">\n";
	}
	std::string emptyElemEnd()
	{
		return "/>\n";
	}


	optional<Point> getMinPoint(std::vector<Point> const & points)
	{
		if (points.empty())
			return optional<Point>();

		Point min = points[0];
		for (unsigned i = 0; i < points.size(); ++i) {
			if (points[i].x < min.x)
				min.x = points[i].x;
			if (points[i].y < min.y)
				min.y = points[i].y;
		}
		return optional<Point>(min);
	}

	optional<Point> getMaxPoint(std::vector<Point> const & points)
	{
		if (points.empty())
			return optional<Point>();

		Point max = points[0];
		for (unsigned i = 0; i < points.size(); ++i) {
			if (points[i].x > max.x)
				max.x = points[i].x;
			if (points[i].y > max.y)
				max.y = points[i].y;
		}
		return optional<Point>(max);
	}

	// Convert coordinates in user space to SVG native space.
	double translateX(double x, Layout const & layout)
	{
		if (layout.origin == Layout::BottomRight || layout.origin == Layout::TopRight)
			return layout.dimensions.width - ((x + layout.origin_offset.x) * layout.scale);
		else
			return (layout.origin_offset.x + x) * layout.scale;
	}

	double translateY(double y, Layout const & layout)
	{
		if (layout.origin == Layout::BottomLeft || layout.origin == Layout::BottomRight)
			return layout.dimensions.height - ((y + layout.origin_offset.y) * layout.scale);
		else
			return (layout.origin_offset.y + y) * layout.scale;
	}
	double translateScale(double dimension, Layout const & layout)
	{
		return dimension * layout.scale;
	}
}