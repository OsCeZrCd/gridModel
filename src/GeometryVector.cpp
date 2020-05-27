#include "GeometryVector.h"
#include <iostream>


void GeometryVector::WriteBinary(std::ostream & ofile, DimensionType dimension) const
{
	ofile.write( (char *)(this->x), sizeof(double)*dimension);
}
void GeometryVector::ReadBinary(std::istream & ifile, DimensionType dimension)
{
	ifile.read( (char *)(this->x), sizeof(double)*dimension);
}
