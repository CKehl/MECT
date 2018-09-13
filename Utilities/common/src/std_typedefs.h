/**
 * created on Jan 9 2018
 * author: Christian Kehl
 */
#ifndef STD_TYPEDEFS_H_
#define STD_TYPEDEFS_H_

#include <stdint.h>
#include <stddef.h>
#include <string>

#ifndef COMPILE_WITH_CUDA

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long ulong;

typedef struct _uchar2
{
	uchar x,y;
} uchar2;

typedef struct _uchar3
{
	uchar x,y,z;
} uchar3;

typedef struct _uchar4
{
	uchar w,x,y,z;
} uchar4;

typedef struct _int2
{
	int x,y;
} int2;

typedef struct _int3
{
	int x,y,z;
} int3;

typedef struct _int4
{
	int x,y,z,w;
} int4;

typedef struct _uint2
{
	uint x,y;
} uint2;

inline uint distanceSquared(_uint2 in1, _uint2 in2) {
	_uint2 dV;
	dV.x = in2.x-in1.x;
	dV.y = in2.y-in1.y;
	return (dV.x*dV.x+dV.y*dV.y);
}

inline bool operator==(const _uint2& lhs, const _uint2& rhs)
{
    return ((lhs.x==rhs.x)&&(lhs.y==rhs.y));
}

typedef struct _uint3
{
	uint x,y,z;
} uint3;

inline uint distanceSquared(_uint3 in1, _uint3 in2) {
	_uint3 dV;
	dV.x = in2.x-in1.x;
	dV.y = in2.y-in1.y;
	dV.z = in2.z-in1.z;
	return (dV.x*dV.x+dV.y*dV.y+dV.z*dV.z);
}

inline bool operator==(const _uint3& lhs, const _uint3& rhs)
{
    return ((lhs.x==rhs.x)&&(lhs.y==rhs.y)&&(lhs.z==rhs.z));
}

typedef struct _uint4
{
	uint w,x,y,z;
} uint4;

typedef struct _float2
{
	float x,y;
} float2;

inline float distanceSquared(_float2 in1, _float2 in2) {
	float2 dV;
	dV.x = in2.x-in1.x;
	dV.y = in2.y-in1.y;
	return (dV.x*dV.x+dV.y*dV.y);
}

inline bool operator==(const _float2& lhs, const _float2& rhs)
{
    return ((lhs.x==rhs.x)&&(lhs.y==rhs.y));
}

typedef struct _float3
{
	float x,y,z;
} float3;

inline uint distanceSquared(_float3 in1, _float3 in2) {
	_float3 dV;
	dV.x = in2.x-in1.x;
	dV.y = in2.y-in1.y;
	dV.z = in2.z-in1.z;
	return (dV.x*dV.x+dV.y*dV.y+dV.z*dV.z);
}

inline bool operator==(const _float3& lhs, const _float3& rhs)
{
    return ((lhs.x==rhs.x)&&(lhs.y==rhs.y)&&(lhs.z==rhs.z));
}

typedef struct _float4
{
	float x,y,z,w;
} float4;

typedef struct _double2
{
	double x,y;
} double2;

typedef struct _double3
{
	double x,y,z;
} double3;

typedef struct _double4
{
	double x,y,z,w;
} double4;

#else
typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long ulong;
#include <vector_types.h>
#endif

/**
 * short x=value, short a=limLow, short b=limHigh
 */
template<typename Type>
inline Type clamp(Type x, Type a, Type b) {
	return x < a ? a : (x > b ? b : x);
}

/**
 * short x=value, short a=limLow, short b=limHigh
 * wrapper-friendly
 */
template<typename Type>
inline Type clampByValue(Type x, Type a, Type b) {
	return x < a ? a : (x > b ? b : x);
}

/**
 * ushort x=value, ushort a=limLow, ushort b=limHigh
 */
template<typename Type>
inline Type clamp(Type& x, Type& a, Type& b) {
	return x < a ? a : (x > b ? b : x);
}

/**
 * ushort x=value, ushort a=limLow, ushort b=limHigh
 * wrapper-friendly
 */
template<typename Type>
inline Type clampByRef(Type& x, Type& a, Type& b) {
	return x < a ? a : (x > b ? b : x);
}

typedef struct _MinMaxCH {
	char min_x, min_y, min_z;
	char max_x, max_y, max_z;

	bool operator==(const _MinMaxCH& rhs) const {
		return (min_x==rhs.min_x && min_y==rhs.min_y && min_z==rhs.min_z && max_x==rhs.max_x && max_y==rhs.max_y && max_z==rhs.max_z);
	}

	bool operator!=(const _MinMaxCH& rhs) const {
		return (min_x!=rhs.min_x || min_y!=rhs.min_y || min_z!=rhs.min_z || max_x!=rhs.max_x || max_y!=rhs.max_y || max_z!=rhs.max_z);
	}
} MinMaxCH;

typedef struct _MinMaxUC {
	uchar min_x, min_y, min_z;
	uchar max_x, max_y, max_z;

	bool operator==(const _MinMaxUC& rhs) const {
		return (min_x==rhs.min_x && min_y==rhs.min_y && min_z==rhs.min_z && max_x==rhs.max_x && max_y==rhs.max_y && max_z==rhs.max_z);
	}

	bool operator!=(const _MinMaxUC& rhs) const {
		return (min_x!=rhs.min_x || min_y!=rhs.min_y || min_z!=rhs.min_z || max_x!=rhs.max_x || max_y!=rhs.max_y || max_z!=rhs.max_z);
	}
} MinMaxUC;

typedef struct _MinMaxS {
	short min_x, min_y, min_z;
	short max_x, max_y, max_z;

	bool operator==(const _MinMaxS& rhs) const {
		return (min_x==rhs.min_x && min_y==rhs.min_y && min_z==rhs.min_z && max_x==rhs.max_x && max_y==rhs.max_y && max_z==rhs.max_z);
	}

	bool operator!=(const _MinMaxS& rhs) const {
		return (min_x!=rhs.min_x || min_y!=rhs.min_y || min_z!=rhs.min_z || max_x!=rhs.max_x || max_y!=rhs.max_y || max_z!=rhs.max_z);
	}
} MinMaxS;

typedef struct _MinMaxUS {
	ushort min_x, min_y, min_z;
	ushort max_x, max_y, max_z;

	bool operator==(const _MinMaxUS& rhs) const {
		return (min_x==rhs.min_x && min_y==rhs.min_y && min_z==rhs.min_z && max_x==rhs.max_x && max_y==rhs.max_y && max_z==rhs.max_z);
	}

	bool operator!=(const _MinMaxUS& rhs) const {
		return (min_x!=rhs.min_x || min_y!=rhs.min_y || min_z!=rhs.min_z || max_x!=rhs.max_x || max_y!=rhs.max_y || max_z!=rhs.max_z);
	}
} MinMaxUS;

typedef struct _MinMaxI {
	int min_x, min_y, min_z;
	int max_x, max_y, max_z;

	bool operator==(const _MinMaxI& rhs) const {
		return (min_x==rhs.min_x && min_y==rhs.min_y && min_z==rhs.min_z && max_x==rhs.max_x && max_y==rhs.max_y && max_z==rhs.max_z);
	}

	bool operator!=(const _MinMaxI& rhs) const {
		return (min_x!=rhs.min_x || min_y!=rhs.min_y || min_z!=rhs.min_z || max_x!=rhs.max_x || max_y!=rhs.max_y || max_z!=rhs.max_z);
	}
} MinMaxI;

typedef struct _MinMaxUI
{
	uint min_x, min_y, min_z;
	uint max_x, max_y, max_z;

	bool operator==(const _MinMaxUI& rhs) const {
		return (min_x==rhs.min_x && min_y==rhs.min_y && min_z==rhs.min_z && max_x==rhs.max_x && max_y==rhs.max_y && max_z==rhs.max_z);
	}

	bool operator!=(const _MinMaxUI& rhs) const {
		return (min_x!=rhs.min_x || min_y!=rhs.min_y || min_z!=rhs.min_z || max_x!=rhs.max_x || max_y!=rhs.max_y || max_z!=rhs.max_z);
	}
} MinMaxUI;

typedef struct _MinMaxL
{
	long min_x, min_y, min_z;
	long max_x, max_y, max_z;

	bool operator==(const _MinMaxL& rhs) const {
		return (min_x==rhs.min_x && min_y==rhs.min_y && min_z==rhs.min_z && max_x==rhs.max_x && max_y==rhs.max_y && max_z==rhs.max_z);
	}

	bool operator!=(const _MinMaxL& rhs) const {
		return (min_x!=rhs.min_x || min_y!=rhs.min_y || min_z!=rhs.min_z || max_x!=rhs.max_x || max_y!=rhs.max_y || max_z!=rhs.max_z);
	}
} MinMaxL;

typedef struct _MinMaxUL
{
	ulong min_x, min_y, min_z;
	ulong max_x, max_y, max_z;

	bool operator==(const _MinMaxUL& rhs) const {
		return (min_x==rhs.min_x && min_y==rhs.min_y && min_z==rhs.min_z && max_x==rhs.max_x && max_y==rhs.max_y && max_z==rhs.max_z);
	}

	bool operator!=(const _MinMaxUL& rhs) const {
		return (min_x!=rhs.min_x || min_y!=rhs.min_y || min_z!=rhs.min_z || max_x!=rhs.max_x || max_y!=rhs.max_y || max_z!=rhs.max_z);
	}
} MinMaxUL;

typedef struct _MinMaxD
{
	double min_x, min_y, min_z;
	double max_x, max_y, max_z;

	bool operator==(const _MinMaxD& rhs) const {
		return (min_x==rhs.min_x && min_y==rhs.min_y && min_z==rhs.min_z && max_x==rhs.max_x && max_y==rhs.max_y && max_z==rhs.max_z);
	}

	bool operator!=(const _MinMaxD& rhs) const {
		return (min_x!=rhs.min_x || min_y!=rhs.min_y || min_z!=rhs.min_z || max_x!=rhs.max_x || max_y!=rhs.max_y || max_z!=rhs.max_z);
	}
} MinMaxD;

typedef struct _MinMaxF
{
	float min_x, min_y, min_z;
	float max_x, max_y, max_z;

	bool operator==(const _MinMaxF& rhs) const {
		return (min_x==rhs.min_x && min_y==rhs.min_y && min_z==rhs.min_z && max_x==rhs.max_x && max_y==rhs.max_y && max_z==rhs.max_z);
	}

	bool operator!=(const _MinMaxF& rhs) const {
		return (min_x!=rhs.min_x || min_y!=rhs.min_y || min_z!=rhs.min_z || max_x!=rhs.max_x || max_y!=rhs.max_y || max_z!=rhs.max_z);
	}
} MinMaxF;

union int_float_bits_32 {
    int32_t int_bits;
    float float_bits;
};

union int_float_bits_64 {
    int64_t int_bits;
    float float_bits;
};

typedef enum _dtype {
	NONE = 0,
	CHAR = 1,
	UCHAR = 2,
	SHORT = 3,
	USHORT = 4,
	INT = 5,
	UINT = 6,
	LONG = 7,
	ULONG = 8,
	FLOAT = 9,
	DOUBLE = 10
} dtype;

inline std::string getDTypeString(dtype typeval) {
	switch(typeval) {
	case NONE: return "NONE";
	case CHAR: return "CHAR";
	case UCHAR: return "UCHAR";
	case SHORT: return "SHORT";
	case USHORT: return "USHORT";
	case INT: return "INT";
	case UINT: return "UINT";
	case LONG: return "LONG";
	case ULONG: return "ULONG";
	case FLOAT: return "FLOAT";
	case DOUBLE: return "DOUBLE";
	default: return "";
	}
	return "";
}

typedef enum _dorder {
	UNDEF_ORDER = 0,
	ROW_MAJOR = 1,
	COLUMN_MAJOR =2
} dorder;

inline std::string getDOrderString(dorder orderval) {
	switch(orderval) {
	case UNDEF_ORDER: return "UNDEFINED";
	case ROW_MAJOR: return "ROW_MAJOR";
	case COLUMN_MAJOR: return "COLUMN_MAJOR";
	default: return "";
	}
	return "";
}

#endif /* STD_TYPEDEFS_H_ */
