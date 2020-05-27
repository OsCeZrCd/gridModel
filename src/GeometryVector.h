#ifndef GEOMETRYVECTOR_INCLUDED
#define GEOMETRYVECTOR_INCLUDED

typedef int DimensionType;


//since we are using GeometryVector to represent different components of the deviatoric strain
//this should be set to the number of degrees of freedom in the deviatoric strain.
const DimensionType MaxDimension=2;

#include <cstring>
#include <cassert>
#include <fstream>
#include <cmath>



//In debug mode, this class also records the dimension of the GeometryVector.
class GeometryVector
{
public:
	double x[ ::MaxDimension];
	GeometryVector()
	{
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]=0;
	}
	//do not change this "int" to "DimensionType"!
	//otherwize GeometryVector(2) would be ambiguous.
	GeometryVector(int dimension)
	{
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]=0;
	}
	GeometryVector(double xx)
	{
		this->x[0]=xx;
		for(int i=1; i< ::MaxDimension; i++)
			this->x[i]=0;
	}
	GeometryVector(double xx, double yy)
	{
		this->x[0]=xx;
		this->x[1]=yy;
		for(int i=2; i< ::MaxDimension; i++)
			this->x[i]=0;
	}
	GeometryVector(const GeometryVector & src)
	{
		#pragma ivdep
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]=src.x[i];
	}
	void AddFrom(const GeometryVector & right)
	{
		#pragma ivdep
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]+=right.x[i];
	}
	void MinusFrom(const GeometryVector & right)
	{
		#pragma ivdep
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]-=right.x[i];
	}
	void MultiplyFrom(const double & right)
	{
		#pragma ivdep
		for(int i=0; i< ::MaxDimension; i++)
			this->x[i]*=right;
	}
	double Dot(const GeometryVector & right) const
	{
		double result=0;
		for(int i=0; i< ::MaxDimension; i++)
			result+=this->x[i]*right.x[i];
		return result;
	}
	double Modulus2(void) const
	{
		return this->Dot(*this);
	}
	bool IsEqual(const GeometryVector & right) const
	{
		for(int i=0; i< ::MaxDimension; i++)
			if(this->x[i]!=right.x[i])
				return false;
		return true;
	}
	void OutputCoordinate(std::ostream & os, DimensionType dim) const
	{
		for(int i=0; i< dim; i++)
			os<<this->x[i]<<" \t";
	}
	void InputCoordinate(std::istream & is, DimensionType dim)
	{
		for(int i=0; i< dim; i++)
			is>>this->x[i];
	}
	void WriteBinary(std::ostream & ofile, DimensionType dimension) const;
	void ReadBinary(std::istream & ifile, DimensionType dimension);

	friend inline GeometryVector SameDimensionZeroVector(const GeometryVector & src);
};

inline std::ostream & operator << (std::ostream & os, const GeometryVector & a)
{
	//os<<a.Dimension<<" \t";
	a.OutputCoordinate(os, ::MaxDimension);
	os<<'\n';
	return os;
}
inline std::istream & operator >> (std::istream & is, GeometryVector & a)
{
	a.InputCoordinate(is, ::MaxDimension);
	return is;
}
inline GeometryVector operator + (const GeometryVector & left, const GeometryVector & right)
{
	GeometryVector result(left);
	result.AddFrom(right);
	return result;
}
inline GeometryVector operator - (const GeometryVector & left, const GeometryVector & right)
{
	GeometryVector result(left);
	result.MinusFrom(right);
	return result;
}
inline GeometryVector operator * (const GeometryVector & left, const double & right)
{
	GeometryVector result(left);
	result.MultiplyFrom(right);
	return result;
}
inline GeometryVector operator * (const double & left, const GeometryVector & right)
{
	GeometryVector result(right);
	result.MultiplyFrom(left);
	return result;
}
inline bool operator == (const GeometryVector & left, const GeometryVector & right)
{
	return left.IsEqual(right);
}
inline bool operator != (const GeometryVector & left, const GeometryVector & right)
{
	return !left.IsEqual(right);
}

#endif
