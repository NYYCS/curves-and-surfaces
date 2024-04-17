#include "curve.h"
#include "vertexrecorder.h"
using namespace std;

const float c_pi = 3.14159265358979323846f;

namespace
{
// Approximately equal to.  We don't want to use == because of
// precision issues with floating point.
inline bool approx(const Vector3f& lhs, const Vector3f& rhs)
{
	const float eps = 1e-8f;
	return (lhs - rhs).absSquared() < eps;
}


}

const Matrix4f Mbez = {
	1.0f, -3.0f,  3.0f, -1.0f,
	0.0f,  3.0f, -6.0f,  3.0f,
	0.0f,  0.0f,  3.0f, -3.0f,
	0.0f,  0.0f,  0.0f,  1.0f,
};

const Matrix4f MbezInv = Mbez.inverse();

const Matrix4f Mb =  (1.0f / 6.0f) * Matrix4f({
	1.0f, -3.0f,  3.0f, -1.0f,
	4.0f,  0.0f, -6.0f,  3.0f,
	1.0f,  3.0f,  3.0f, -3.0f,
	0.0f,  0.0f,  0.0f,  1.0f,
});

Vector3f bezier(Matrix4f Gbez, float t)
{	
	Vector4f V = Gbez * Mbez * Vector4f(1, t, t*t, t*t*t);
	return { V[0], V[1], V[2] };
}

Vector3f bezierD1(Matrix4f Gbez, float t)
{
	Vector4f V = Gbez * Mbez * Vector4f(0, 1, 2*t, 3*t*t);
	return { V[0], V[1], V[2] };
};

Curve evalBezier(const vector< Vector3f >& P, unsigned steps)
{
	// Check
	if (P.size() < 4 || P.size() % 3 != 1)
	{
		cerr << "evalBezier must be called with 3n+1 control points." << endl;
		exit(0);
	}

	unsigned n = P.size() - 3;

	// TODO:
	// You should implement this function so that it returns a Curve
	// (e.g., a vector< CurvePoint >).  The variable "steps" tells you
	// the number of points to generate on each piece of the spline.
	// At least, that's how the sample solution is implemented and how
	// the SWP files are written.  But you are free to interpret this
	// variable however you want, so long as you can control the
	// "resolution" of the discretized spline curve with it.

	// Make sure that this function computes all the appropriate
	// Vector3fs for each CurvePoint: V,T,N,B.
	// [NBT] should be unit and orthogonal.

	// Also note that you may assume that all Bezier curves that you
	// receive have G1 continuity.  Otherwise, the TNB will not be
	// be defined at points where this does not hold.

	cerr << "\t>>> evalBezier has been called with the following input:" << endl;

	cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;

	for (int i = 0; i < (int)P.size(); ++i)
	{
		cerr << "\t>>> " << P[i] << endl;
	}

	cerr << "\t>>> Steps (type steps): " << steps << endl;

	Curve C(n * steps + 1);
	Vector3f B0 = { 0.0f, 0.0f, 1.0f };

	for (unsigned i = 0; i < n; ++i)
	{

		Matrix4f Gbez = {
			Vector4f(P[i], 0.0f),
			Vector4f(P[i+1], 0.0f),
			Vector4f(P[i+2], 0.0f),
			Vector4f(P[i+3], 0.0f),
		};

		for (unsigned j = 0; j <= steps; j++)
		{
			unsigned k = i * steps + j;
			float t = float(j) / steps;

			C[k].V = bezier(Gbez, t);

			C[k].T = bezierD1(Gbez, t).normalized();

			C[k].N = Vector3f::cross(B0, C[k].T).normalized();
			
			C[k].B = Vector3f::cross(C[k].T, C[k].N).normalized();

			B0 = C[k].B;
		}
	}

	return C;
}

Curve evalBspline(const vector< Vector3f >& P, unsigned steps)
{
	// Check
	if (P.size() < 4)
	{
		cerr << "evalBspline must be called with 4 or more control points." << endl;
		exit(0);
	}

	unsigned n = P.size() - 3;

	// TODO:
	// It is suggested that you implement this function by changing
	// basis from B-spline to Bezier.  That way, you can just call
	// your evalBezier function.

	cerr << "\t>>> evalBSpline has been called with the following input:" << endl;
	cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;

	for (int i = 0; i < (int)P.size(); ++i)
	{
		cerr << "\t>>> " << P[i] << endl;
	}

	Curve C(n * steps + 1);
	Vector3f B0 = { 0.0f, 0.0f, 1.0f };

	for (int i = 0; i < n; ++i)
	{

		Matrix4f Gb = {
			Vector4f(P[i], 0.0f),
			Vector4f(P[i+1], 0.0f),
			Vector4f(P[i+2], 0.0f),
			Vector4f(P[i+3], 0.0f),
		};

		Matrix4f Gbez = Gb * Mb * MbezInv;

		for (unsigned j = 0; j <= steps; j++)
		{
			unsigned k = i * steps + j;
			float t = float(j) / steps;

			C[k].V = bezier(Gbez, t);
			
			C[k].T = bezierD1(Gbez, t).normalized();

			C[k].N = Vector3f::cross(B0, C[k].T).normalized();
			
			C[k].B = Vector3f::cross(C[k].T, C[k].N).normalized();

			B0 = C[k].B;
		}
	}

	return C;
}

Curve evalCircle(float radius, unsigned steps)
{
	// This is a sample function on how to properly initialize a Curve
	// (which is a vector< CurvePoint >).

	// Preallocate a curve with steps+1 CurvePoints
	Curve R(steps + 1);

	// Fill it in counterclockwise
	for (unsigned i = 0; i <= steps; ++i)
	{
		// step from 0 to 2pi
		float t = 2.0f * c_pi * float(i) / steps;

		// Initialize position
		// We're pivoting counterclockwise around the y-axis
		R[i].V = radius * Vector3f(cos(t), sin(t), 0);

		// Tangent vector is first derivative
		R[i].T = Vector3f(-sin(t), cos(t), 0);

		// Normal vector is second derivative
		R[i].N = Vector3f(-cos(t), -sin(t), 0);

		// Finally, binormal is facing up.
		R[i].B = Vector3f(0, 0, 1);
	}

	return R;
}

void recordCurve(const Curve& curve, VertexRecorder* recorder)
{
	const Vector3f WHITE(1, 1, 1);

	for (int i = 0; i < (int)curve.size() - 1; ++i)
	{
		recorder->record_poscolor(curve[i].V, WHITE);
		recorder->record_poscolor(curve[i + 1].V, WHITE);
	}
}
void recordCurveFrames(const Curve& curve, VertexRecorder* recorder, float framesize)
{
	Matrix4f T;
	const Vector3f RED(1, 0, 0);
	const Vector3f GREEN(0, 1, 0);
	const Vector3f BLUE(0, 0, 1);
	
	const Vector4f ORGN(0, 0, 0, 1);
	const Vector4f AXISX(framesize, 0, 0, 1);
	const Vector4f AXISY(0, framesize, 0, 1);
	const Vector4f AXISZ(0, 0, framesize, 1);

	for (int i = 0; i < (int)curve.size(); ++i)
	{
		T.setCol(0, Vector4f(curve[i].N, 0));
		T.setCol(1, Vector4f(curve[i].B, 0));
		T.setCol(2, Vector4f(curve[i].T, 0));
		T.setCol(3, Vector4f(curve[i].V, 1));
 
		// Transform orthogonal frames into model space
		Vector4f MORGN  = T * ORGN;
		Vector4f MAXISX = T * AXISX;
		Vector4f MAXISY = T * AXISY;
		Vector4f MAXISZ = T * AXISZ;

		// Record in model space
		recorder->record_poscolor(MORGN.xyz(), RED);
		recorder->record_poscolor(MAXISX.xyz(), RED);

		recorder->record_poscolor(MORGN.xyz(), GREEN);
		recorder->record_poscolor(MAXISY.xyz(), GREEN);

		recorder->record_poscolor(MORGN.xyz(), BLUE);
		recorder->record_poscolor(MAXISZ.xyz(), BLUE);
	}
}

