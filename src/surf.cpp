#include "surf.h"
#include "vertexrecorder.h"
using namespace std;

const float c_pi = 3.14159265358979323846f;


namespace
{

    inline bool approx(const Vector3f& lhs, const Vector3f& rhs)
    {
	    const float eps = 1e-8f;
	    return (lhs - rhs).absSquared() < eps;
    }
    // We're only implenting swept surfaces where the profile curve is
    // flat on the xy-plane.  This is a check function.
    static bool checkFlat(const Curve &profile)
    {
        for (unsigned i=0; i<profile.size(); i++)
            if (profile[i].V[2] != 0.0 ||
                profile[i].T[2] != 0.0 ||
                profile[i].N[2] != 0.0)
                return false;
    
        return true;
    }
}



// DEBUG HELPER
Surface quad() { 
	Surface ret;
	ret.VV.push_back(Vector3f(-1, -1, 0));
	ret.VV.push_back(Vector3f(+1, -1, 0));
	ret.VV.push_back(Vector3f(+1, +1, 0));
	ret.VV.push_back(Vector3f(-1, +1, 0));

	ret.VN.push_back(Vector3f(0, 0, 1));
	ret.VN.push_back(Vector3f(0, 0, 1));
	ret.VN.push_back(Vector3f(0, 0, 1));
	ret.VN.push_back(Vector3f(0, 0, 1));

	ret.VF.push_back(Tup3u(0, 1, 2));
	ret.VF.push_back(Tup3u(0, 2, 3));
	return ret;
}

void makeTriangleMeshes(Surface &surface, unsigned m, unsigned n)
{   
    for (unsigned i = 0; i < m; ++i) 
    {
        for (unsigned j = 0; j <= n; ++j)
        {
            unsigned k = i * n + j;
            unsigned t = k + n;

            surface.VF.push_back({ k, t + 1, k + 1 });
            surface.VF.push_back({ k + 1, t + 1, t + 2 });
        }
    }
}

Surface makeSurfRev(const Curve &profile, unsigned steps)
{
    Surface surface;
    
    if (!checkFlat(profile))
    {
        cerr << "surfRev profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // TODO: Here you should build the surface.  See surf.h for details.

    for (unsigned i = 0; i < profile.size(); ++i)
    {
        for (unsigned j = 0; j <= steps; ++j)
        {
            float t = 2.0f * c_pi * float(j) / steps;

            Matrix4f M = Matrix4f::rotateY(t);
            Vector4f P = M * Vector4f(profile[i].V, 1.0f);
            Vector4f N = -(M.inverse().transposed() * Vector4f(profile[i].N, 0.0f)).normalized();

            surface.VV.push_back(P.xyz());
            surface.VN.push_back(N.xyz());
        }
    }

    makeTriangleMeshes(surface, profile.size(), steps);
 
    return surface;
}

Surface makeGenCyl(const Curve &profile, const Curve &sweep )
{
    Surface surface;

    if (!checkFlat(profile))
    {
        cerr << "genCyl profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    float a = -acos(Vector3f::dot(sweep[0].N, sweep[sweep.size() - 1].N));
    float t = a / sweep.size();

    for (unsigned i = 0; i < profile.size(); ++i)
    {   
        for (unsigned j = 0; j < sweep.size(); ++j)
        {   

            Matrix4f M = {
                Vector4f(cos(j * t) * sweep[j].N + sin(j * t) * sweep[j].B, 0.0f),
                Vector4f(cos(j * t) * sweep[j].B - sin(j * t) * sweep[j].N, 0.0f),
                Vector4f(sweep[j].T, 0.0f),
                Vector4f(sweep[j].V, 1.0f),
            };

            Vector4f P = M * Vector4f(profile[i].V, 1.0f);
            Vector4f N = -(M.inverse().transposed() * Vector4f(profile[i].N, 0.0f)).normalized();

            surface.VV.push_back(P.xyz());
            surface.VN.push_back(N.xyz());

        }
    }

    // TODO: Here you should build the surface.  See surf.h for details.

    makeTriangleMeshes(surface, profile.size(), sweep.size());

    return surface;
}

void recordSurface(const Surface &surface, VertexRecorder* recorder) {
	const Vector3f WIRECOLOR(0.4f, 0.4f, 0.4f);
    for (int i=0; i<(int)surface.VF.size(); i++)
    {
		recorder->record(surface.VV[surface.VF[i][0]], surface.VN[surface.VF[i][0]], WIRECOLOR);
		recorder->record(surface.VV[surface.VF[i][1]], surface.VN[surface.VF[i][1]], WIRECOLOR);
		recorder->record(surface.VV[surface.VF[i][2]], surface.VN[surface.VF[i][2]], WIRECOLOR);
    }
}

void recordNormals(const Surface &surface, VertexRecorder* recorder, float len)
{
	const Vector3f NORMALCOLOR(0, 1, 1);
    for (int i=0; i<(int)surface.VV.size(); i++)
    {
		recorder->record_poscolor(surface.VV[i], NORMALCOLOR);
		recorder->record_poscolor(surface.VV[i] + surface.VN[i] * len, NORMALCOLOR);
    }
}

void outputObjFile(ostream &out, const Surface &surface)
{
    
    for (int i=0; i<(int)surface.VV.size(); i++)
        out << "v  "
            << surface.VV[i][0] << " "
            << surface.VV[i][1] << " "
            << surface.VV[i][2] << endl;

    for (int i=0; i<(int)surface.VN.size(); i++)
        out << "vn "
            << surface.VN[i][0] << " "
            << surface.VN[i][1] << " "
            << surface.VN[i][2] << endl;

    out << "vt  0 0 0" << endl;
    
    for (int i=0; i<(int)surface.VF.size(); i++)
    {
        out << "f  ";
        for (unsigned j=0; j<3; j++)
        {
            unsigned a = surface.VF[i][j]+1;
            out << a << "/" << "1" << "/" << a << " ";
        }
        out << endl;
    }
}
