using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using Unity.VisualScripting;
using UnityEditor;
using UnityEngine;
using static UnityEditor.Searcher.SearcherWindow.Alignment;

public class CoordConversionDouble
{
    // Start is called before the first frame update
    void Start()
    {

    }

    // Update is called once per frame
    void Update()
    {

    }

    //Conformal Dymaxion properties
    protected static readonly double VECTOR_SCALE_FACTOR = 1.0f / 1.1473979730192934f;
    protected static readonly int SIDE_LENGTH = 256;

    protected static readonly double[][] vx = Enumerable.Range(0, SIDE_LENGTH + 1)
    .Select(i => new double[SIDE_LENGTH + 1 - i])
    .ToArray();

    protected static readonly double[][] vy = Enumerable.Range(0, SIDE_LENGTH + 1)
    .Select(i => new double[SIDE_LENGTH + 1 - i])
    .ToArray();

    //--------------------------------------------------------------------------------

    protected static readonly double THETA = -150 * Math.PI / 180;
    protected static readonly double SIN_THETA = Math.Sin(THETA);
    protected static readonly double COS_THETA = Math.Cos(THETA);
    protected static readonly double BERING_X = -0.3420420960118339f;//-0.3282152608138795;
    protected static readonly double BERING_Y = -0.322211064085279f;//-0.3281491467713469;
    protected static readonly double ARCTIC_Y = -0.2f;//-0.3281491467713469;
    protected static readonly double ARCTIC_M = (ARCTIC_Y - MathUtilsDouble.ROOT3 * ARC / 4f) / (BERING_X - -0.5f * ARC);
    protected static readonly double ARCTIC_B = ARCTIC_Y - ARCTIC_M * BERING_X;
    protected static readonly double ALEUTIAN_Y = -0.5000446805492526f;//-0.5127463765943157;
    protected static readonly double ALEUTIAN_XL = -0.5149231279757507f;//-0.4957832938238718;
    protected static readonly double ALEUTIAN_XR = -0.45f;
    protected static readonly double ALEUTIAN_M = (BERING_Y - ALEUTIAN_Y) / (BERING_X - ALEUTIAN_XR);
    protected static readonly double ALEUTIAN_B = BERING_Y - ALEUTIAN_M * BERING_X;

    protected static double ARC = 2 * Math.Asin(Math.Sqrt(5 - Math.Sqrt(5)) / Math.Sqrt(10));
    protected static readonly double Z = Math.Sqrt(5 + 2 * Math.Sqrt(5)) / Math.Sqrt(15);
    protected static readonly double EL = Math.Sqrt(8) / Math.Sqrt(5 + Math.Sqrt(5));
    protected static readonly double EL6 = EL / 6;
    protected static readonly double DVE = Math.Sqrt(3 + Math.Sqrt(5)) / Math.Sqrt(5 + Math.Sqrt(5));
    protected static readonly double R = -3 * EL6 / DVE;

    /*
     * Number of iterations for Newton's method
     */

    private static readonly int NEWTON = 5;

    /**
     * This contains the vertices of the icosahedron,
     * identified by their geographic longitude and latitude in degrees.
     * When the class is loaded, a static block below converts all these coordinates
     * to the equivalent spherical coordinates (longitude and colatitude), in radians.
     *
     * @see <a href="https://en.wikipedia.org/wiki/Regular_icosahedron#Spherical_coordinates">Wikipedia</a>
     */
    protected static readonly double[][] VERTICES = {
            new double[] { 10.536199f, 64.700000f },
            new double[] { -5.245390f, 2.300882f },
            new double[] { 58.157706f, 10.447378f },
            new double[] { 122.300000f, 39.100000f },
            new double[] { -143.478490f, 50.103201f },
            new double[] { -67.132330f, 23.717925f },
            new double[] { 36.521510f, -50.103200f },
            new double[] { 112.867673f, -23.717930f },
            new double[] { 174.754610f, -2.300882f },
            new double[] { -121.842290f, -10.447350f },
            new double[] { -57.700000f, -39.100000f },
            new double[] { -169.463800f, -64.700000f },
    };

    /**
     * Indicates the vertices forming each face of the icosahedron.
     * Each entry refers to the index of a vertex in {@link #VERTICES}
     */
    protected static readonly int[,] ISO = {
            { 2, 1, 6 },
            { 1, 0, 2 },
            { 0, 1, 5 },
            { 1, 5, 10 },
            { 1, 6, 10 },
            { 7, 2, 6 },
            { 2, 3, 7 },
            { 3, 0, 2 },
            { 0, 3, 4 },
            { 4, 0, 5 }, //9, qubec
            { 5, 4, 9 },
            { 9, 5, 10 },
            { 10, 9, 11 },
            { 11, 6, 10 },
            { 6, 7, 11 },
            { 8, 3, 7 },
            { 8, 3, 4 },
            { 8, 4, 9 },
            { 9, 8, 11 },
            { 7, 8, 11 },
            { 11, 6, 7 }, //child of 14
            { 3, 7, 8 } //child of 15
    };

    protected static readonly double[,] CENTER_MAP = {
        { -3, 7 },
        { -2, 5 },
        { -1, 7 },
        { 2, 5 },
        { 4, 5 },
        { -4, 1 },
        { -3, -1 },
        { -2, 1 },
        { -1, -1 },
        { 0, 1 },
        { 1, -1 },
        { 2, 1 },
        { 3, -1 },
        { 4, 1 },
        { 5, -1 }, //14, left side, right to be cut
        { -3, -5 },
        { -1, -5 },
        { 1, -5 },
        { 2, -7 },
        { -4, -7 },
        { -5, -5 }, //20, pseudo triangle, child of 14
        { -2, -7 } //21 , pseudo triangle, child of 15
    };

    /**
    * Indicates for each face if it needs to be flipped after projecting
    */
    protected static readonly bool[] FLIP_TRIANGLE = {
            true, false, true, false, false,
            true, false, true, false, true, false, true, false, true, false,
            true, true, true, false, false,
            true, false
    };

    /**
    * This contains the Cartesian coordinates the centroid
    * of each face of the icosahedron.
    */
    protected static readonly double[][] CENTROIDS = Array.ConvertAll(new double[22], x => new double[3]);

    /**
    * Rotation matrices to move the triangles to the reference coordinates from the original positions.
    * Indexed by the face's indices.
    */
    protected static readonly double[][][] ROTATION_MATRICES = Array.ConvertAll(new double[22], x =>
                                                              Array.ConvertAll(new double[3], y => new double[3]));

    /**
    * Rotation matrices to move the triangles from the reference coordinates to their original positions.
    * Indexed by the face's indices.
    */
    protected static readonly double[][][] INVERSE_ROTATION_MATRICES = Array.ConvertAll(new double[22], x =>
                                                                      Array.ConvertAll(new double[3], y => new double[3]));
    //protected static readonly double[][][] INVERSE_ROTATION_MATRICES = new double[22][][];

    protected static readonly int[] FACE_ON_GRID = {
            -1, -1, 0, 1, 2, -1, -1, 3, -1, 4, -1,
            -1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
            20, 19, 15, 21, 16, -1, 17, 18, -1, -1, -1,
    };

    static CoordConversionDouble()
    {

        for (int i = 0; i < 22; i++)
        {
            CENTER_MAP[i, 0] *= 0.5f * ARC;
            CENTER_MAP[i, 1] *= ARC * MathUtilsDouble.ROOT3 / 12f;
        }

        // Will contain the list of vertices in Cartesian coordinates
        double[][] verticesCartesian = new double[VERTICES.Length][];
        for (int i = 0; i < verticesCartesian.Length; i++)
        {
            verticesCartesian[i] = new double[3];
        }

        // Convert the geographic vertices to spherical in radians
        for (int i = 0; i < VERTICES.Length; i++)
        {
            double[] vertexSpherical = MathUtilsDouble.Geo2Spherical(VERTICES[i]);
            double[] vertex = MathUtilsDouble.Spherical2Cartesian(vertexSpherical);
            verticesCartesian[i] = vertex;
            VERTICES[i] = vertexSpherical;
        }

        for (int i = 0; i < 22; i++)
        {

            // Vertices of the current face
            double[] vec1 = verticesCartesian[ISO[i, 0]];
            double[] vec2 = verticesCartesian[ISO[i, 1]];
            double[] vec3 = verticesCartesian[ISO[i, 2]];

            // Find the centroid's projection onto the sphere
            double xsum = vec1[0] + vec2[0] + vec3[0];
            double ysum = vec1[1] + vec2[1] + vec3[1];
            double zsum = vec1[2] + vec2[2] + vec3[2];
            double mag = Math.Sqrt(xsum * xsum + ysum * ysum + zsum * zsum);
            CENTROIDS[i] = new double[] { xsum / mag, ysum / mag, zsum / mag };

            double[] centroidSpherical = MathUtilsDouble.Cartesian2Spherical(CENTROIDS[i]);
            double centroidLambda = centroidSpherical[0];
            double centroidPhi = centroidSpherical[1];

            double[] vertex = VERTICES[ISO[i, 0]];
            double[] v = { vertex[0] - centroidLambda, vertex[1] };
            v = YRot(v, -centroidPhi);

            ROTATION_MATRICES[i] = MathUtilsDouble.ProduceZYZRotationMatrix(-centroidLambda, -centroidPhi, (Math.PI / 2) - v[0]);
            INVERSE_ROTATION_MATRICES[i] = MathUtilsDouble.ProduceZYZRotationMatrix(v[0] - (Math.PI / 2), centroidPhi, centroidLambda);

            //Conformal Dymaxion InvertableVectorField START----------------------------------------------

            // Specify the path to the extracted LZMA file
            string assetsPath = Application.dataPath;
            string lzmaFilePath = assetsPath + "/conformal";

            using FileStream fileStream = new FileStream(lzmaFilePath, FileMode.Open, FileAccess.Read);
            using BinaryReader reader = new BinaryReader(fileStream);
            for (int v1 = 0; v1 < SIDE_LENGTH + 1; v1++)
            {
                for (int u = 0; u < SIDE_LENGTH + 1 - v1; u++)
                {
                    byte[] rawVxBytes = reader.ReadBytes(sizeof(double));
                    if (BitConverter.IsLittleEndian)
                        Array.Reverse(rawVxBytes);

                    double rawVx = BitConverter.ToDouble(rawVxBytes, 0);

                    byte[] rawVyBytes = reader.ReadBytes(sizeof(double));
                    if (BitConverter.IsLittleEndian)
                        Array.Reverse(rawVyBytes);

                    double rawVy = BitConverter.ToDouble(rawVyBytes, 0);

                    vx[u][v1] = (double)rawVx * VECTOR_SCALE_FACTOR;
                    vy[u][v1] = (double)rawVy * VECTOR_SCALE_FACTOR;
                }
            }

        }
    }

    protected static bool IsEurasianPart(double x, double y)
    {

        //catch vast majority of cases in not near boundary
        if (x > 0)
        {
            return false;
        }
        if (x < -0.5 * ARC)
        {
            return true;
        }

        if (y > MathUtilsDouble.ROOT3 * ARC / 4) //above arctic ocean
        {
            return x < 0;
        }

        if (y < ALEUTIAN_Y) //below bering sea
        {
            return y < (ALEUTIAN_Y + ALEUTIAN_XL) - x;
        }

        if (y > BERING_Y)
        { //boundary across arctic ocean

            if (y < ARCTIC_Y)
            {
                return x < BERING_X; //in strait
            }

            return y < ARCTIC_M * x + ARCTIC_B; //above strait
        }

        return y > ALEUTIAN_M * x + ALEUTIAN_B;
    }

    private static int FindTriangleGrid(double x, double y)
    {
        //cast equilateral triangles to 45 degrees right triangles (side length of root2)
        double xp = x / ARC;
        double yp = y / (ARC * MathUtilsDouble.ROOT3);

        int row;
        if (yp > -0.25f)
        {
            if (yp < 0.25f)
            { //middle
                row = 1;
            }
            else if (yp <= 0.75)
            { //top
                row = 0;
                yp = 0.5f - yp; //translate to middle and flip
            }
            else
            {
                return -1;
            }
        }
        else if (yp >= -0.75f)
        { //bottom
            row = 2;
            yp = -yp - 0.5f; //translate to middle and flip
        }
        else
        {
            return -1;
        }

        yp += 0.25f; //change origin to vertex 4, to allow grids to align

        //rotate coords 45 degrees so left and right sides of the triangle become the x/y axies (also side lengths are now 1)
        double xr = xp - yp;
        double yr = xp + yp;

        //assign a order to what grid along the y=x line it is
        int gx = (int)Math.Floor(xr);
        int gy = (int)Math.Floor(yr);

        int col = 2 * gx + (gy != gx ? 1 : 0) + 6;

        //out of bounds
        if (col < 0 || col >= 11)
        {
            return -1;
        }

        return FACE_ON_GRID[row * 11 + col]; //get face at this position
    }

    protected static double[] YRot(double[] spherical, double rot)
    {
        double[] c = MathUtilsDouble.Spherical2Cartesian(spherical);

        double x = c[0];
        c[0] = c[2] * Math.Sin(rot) + x * Math.Cos(rot);
        c[2] = c[2] * Math.Cos(rot) - x * Math.Sin(rot);

        double mag = Math.Sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
        c[0] /= mag;
        c[1] /= mag;
        c[2] /= mag;

        return new double[]{
                Math.Atan2(c[1], c[0]),
                Math.Atan2(Math.Sqrt(c[0] * c[0] + c[1] * c[1]), c[2])
        };
    }

    /**
     * Finds the face of the icosahedron on which to project a point.
     * In practice, it works by finding the face with the closest centroid to the point.
     *
     * @param vector - position vector as double array of length 3, using Cartesian coordinates
     * @return an integer identifying the face on which to project the point
     */
    protected static int FindTriangle(double[] vector)
    {

        double min = double.MaxValue;
        int face = 0;

        for (int i = 0; i < 20; i++)
        {
            double xd = CENTROIDS[i][0] - vector[0];
            double yd = CENTROIDS[i][1] - vector[1];
            double zd = CENTROIDS[i][2] - vector[2];

            double dissq = xd * xd + yd * yd + zd * zd;
            if (dissq < min)
            {

                if (dissq < 0.1) //TODO: enlarge radius
                {
                    return i;
                }

                face = i;
                min = dissq;
            }
        }

        return face;
    }

    private static double[] InverseTriangleTransformNewton(double xpp, double ypp)
    {
        //a & b are linearly related to c, so using the tan of sum formula we know: tan(c+off) = (tanc + tanoff)/(1-tanc*tanoff)
        double tanaoff = Math.Tan(MathUtilsDouble.ROOT3 * ypp + xpp); // a = c + root3*y'' + x''
        double tanboff = Math.Tan(2 * xpp); // b = c + 2x''

        double anumer = tanaoff * tanaoff + 1;
        double bnumer = tanboff * tanboff + 1;

        //we will be solving for tanc, starting at t=0, tan(0) = 0
        double tana = tanaoff;
        double tanb = tanboff;
        double tanc = 0;

        double adenom = 1;
        double bdenom = 1;

        //double fp = anumer + bnumer + 1; //derivative relative to tanc

        //int i = newton;
        for (int i = 0; i < NEWTON; i++)
        {
            double f = tana + tanb + tanc - R; //R = tana + tanb + tanc
            double fp = anumer * adenom * adenom + bnumer * bdenom * bdenom + 1; //derivative relative to tanc

            //TODO: fp could be simplified on first loop: 1 + anumer + bnumer

            tanc -= f / fp;

            adenom = 1 / (1 - tanc * tanaoff);
            bdenom = 1 / (1 - tanc * tanboff);

            tana = (tanc + tanaoff) * adenom;
            tanb = (tanc + tanboff) * bdenom;
        }

        //simple reversal algebra based on tan values
        double yp = MathUtilsDouble.ROOT3 * (DVE * tana + EL6) / 2;
        double xp = DVE * tanb + yp / MathUtilsDouble.ROOT3 + EL6;

        //x = z*xp/Z, y = z*yp/Z, x^2 + y^2 + z^2 = 1
        double xpoZ = xp / Z;
        double ypoZ = yp / Z;

        double z = 1 / Math.Sqrt(1 + xpoZ * xpoZ + ypoZ * ypoZ);

        return new double[] { z * xpoZ, z * ypoZ, z };
    }

    protected static double[] TriangleTransform(double[] vec)
    {

        double S = Z / vec[2];

        double xp = S * vec[0];
        double yp = S * vec[1];

        double a = Math.Atan((2 * yp / MathUtilsDouble.ROOT3 - EL6) / DVE); //ARC/2 terms cancel
        double b = Math.Atan((xp - yp / MathUtilsDouble.ROOT3 - EL6) / DVE);
        double c = Math.Atan((-xp - yp / MathUtilsDouble.ROOT3 - EL6) / DVE);

        return new double[] { 0.5f * (b - c), (2 * a - b - c) / (2 * MathUtilsDouble.ROOT3) };
    }

    protected static double[] InverseTriangleTransform(double x, double y)
    {
        x /= ARC;
        y /= ARC;

        x += 0.5f;
        y += MathUtilsDouble.ROOT3 / 6;

        double[] c = GetInterpolatedVector(x, y);
        return InverseTriangleTransformNewton(c[0], c[1]);
    }

    public static double[] FromGeo(double longitude, double latitude)
    {
        //DymaxionProjection START----------------------------------------------
        double[] vector = MathUtilsDouble.Spherical2Cartesian(MathUtilsDouble.Geo2Spherical(new double[] { longitude, latitude }));

        int face = FindTriangle(vector);

        //apply rotation matrix (move triangle onto template triangle)
        double[] pvec = MathUtilsDouble.MatVecProdD(ROTATION_MATRICES[face], vector);
        double[] projectedVec = TriangleTransform(pvec);

        //flip triangle to correct orientation
        if (FLIP_TRIANGLE[face])
        {
            projectedVec[0] = -projectedVec[0];
            projectedVec[1] = -projectedVec[1];
        }

        vector[0] = projectedVec[0];
        //deal with special snowflakes (child faces 20, 21)
        if (((face == 15 && vector[0] > projectedVec[1] * MathUtilsDouble.ROOT3) || face == 14) && vector[0] > 0)
        {
            projectedVec[0] = 0.5f * vector[0] - 0.5f * MathUtilsDouble.ROOT3 * projectedVec[1];
            projectedVec[1] = 0.5f * MathUtilsDouble.ROOT3 * vector[0] + 0.5f * projectedVec[1];
            face += 6; //shift 14->20 & 15->21
        }

        projectedVec[0] += CENTER_MAP[face, 0];
        projectedVec[1] += CENTER_MAP[face, 1];

        //BTEDymaxionProjection START----------------------------------------------
        double[] c = projectedVec;
        double x = c[0];
        double y = c[1];

        bool easia = IsEurasianPart(x, y);

        y -= 0.75f * ARC * MathUtilsDouble.ROOT3;

        if (easia)
        {
            x += ARC;

            double t = x;
            x = COS_THETA * x - SIN_THETA * y;
            y = SIN_THETA * t + COS_THETA * y;

        }
        else
        {
            x -= ARC;
        }

        c[0] = y;
        c[1] = -x;
        //FlipVerticalProjectionTransform START----------------------------------------------
        c[1] = -c[1];
        //BTEDymaxionProjection END----------------------------------------------
        c[0] *= 7318261.522857145f;
        c[1] *= 7318261.522857145f;
        return c;
    }


    public static double[] ToGeo(double x, double y)
    {
        //ScaleProjectionTransform START----------------------------------------------
        x /= 7318261.522857145f;
        y /= 7318261.522857145f;
        //FlipVerticalProjectionTransform START----------------------------------------------
        y = -y;
        //BTEDymaxionProjection START----------------------------------------------
        bool easia;
        if (y < 0)
        {
            easia = x > 0;
        }
        else if (y > ARC / 2)
        {
            easia = x > -MathUtilsDouble.ROOT3 * ARC / 2;
        }
        else
        {
            easia = y * -MathUtilsDouble.ROOT3 < x;
        }

        double t = x;
        x = -y;
        y = t;

        if (easia)
        {
            t = x;
            x = COS_THETA * x + SIN_THETA * y;
            y = COS_THETA * y - SIN_THETA * t;
            x -= ARC;

        }
        else
        {
            x += ARC;
        }

        y += 0.75f * ARC * MathUtilsDouble.ROOT3;

        //DymaxionProjection START----------------------------------------------
        int face = FindTriangleGrid(x, y);

        x -= CENTER_MAP[face, 0];
        y -= CENTER_MAP[face, 1]; //da cu round mai tare

        //flip triangle to upright orientation (if not already)
        if (FLIP_TRIANGLE[face])
        {
            x = -x;
            y = -y;
        }

        //invert triangle transform
        double[] c = InverseTriangleTransform(x, y);
        x = c[0];
        y = c[1];
        double z = c[2];

        double[] vec = { x, y, z };
        //apply inverse rotation matrix (move triangle from template triangle to correct position on globe)
        double[] vecp = MathUtilsDouble.MatVecProdD(INVERSE_ROTATION_MATRICES[face], vec);

        //convert back to geo coordinates
        return MathUtilsDouble.Spherical2Geo(MathUtilsDouble.Cartesian2Spherical(vecp));
    }

    public static double[] GetInterpolatedVector(double x, double y)
    {
        //scale up triangle to be triangleSize across
        x *= SIDE_LENGTH;
        y *= SIDE_LENGTH;

        //convert to triangle units
        double v = 2 * y / MathUtilsDouble.ROOT3;
        double u = x - v * 0.5f;

        int u1 = (int)u;
        int v1 = (int)v;

        if (u1 < 0)
        {
            u1 = 0;
        }
        else if (u1 >= SIDE_LENGTH)
        {
            u1 = SIDE_LENGTH - 1;
        }

        if (v1 < 0)
        {
            v1 = 0;
        }
        else if (v1 >= SIDE_LENGTH - u1)
        {
            v1 = SIDE_LENGTH - u1 - 1;
        }

        double valx1;
        double valy1;
        double valx2;
        double valy2;
        double valx3;
        double valy3;
        double y3;
        double x3;

        double flip = 1;

        if (y < -MathUtilsDouble.ROOT3 * (x - u1 - v1 - 1) || v1 == SIDE_LENGTH - u1 - 1)
        {
            valx1 = vx[u1][v1];
            valy1 = vy[u1][v1];
            valx2 = vx[u1][v1 + 1];
            valy2 = vy[u1][v1 + 1];
            valx3 = vx[u1 + 1][v1];
            valy3 = vy[u1 + 1][v1];

            y3 = 0.5f * MathUtilsDouble.ROOT3 * v1;
            x3 = (u1 + 1) + 0.5f * v1;
        }
        else
        {
            valx1 = vx[u1][v1 + 1];
            valy1 = vy[u1][v1 + 1];
            valx2 = vx[u1 + 1][v1];
            valy2 = vy[u1 + 1][v1];
            valx3 = vx[u1 + 1][v1 + 1];
            valy3 = vy[u1 + 1][v1 + 1];

            flip = -1;
            y = -y;

            y3 = -(0.5f * MathUtilsDouble.ROOT3 * (v1 + 1));
            x3 = (u1 + 1) + 0.5f * (v1 + 1);
        }

        //TODO: not sure if weights are right (but weirdly mirrors stuff so there may be simplifcation yet)
        double w1 = -(y - y3) / MathUtilsDouble.ROOT3 - (x - x3);
        double w2 = 2 * (y - y3) / MathUtilsDouble.ROOT3;
        double w3 = 1 - w1 - w2;

        return new double[]{ valx1 * w1 + valx2 * w2 + valx3 * w3, valy1 * w1 + valy2 * w2 + valy3 * w3,
                    (valx3 - valx1) * SIDE_LENGTH, SIDE_LENGTH * flip * (2 * valx2 - valx1 - valx3) / MathUtilsDouble.ROOT3,
                    (valy3 - valy1) * SIDE_LENGTH, SIDE_LENGTH * flip * (2 * valy2 - valy1 - valy3) / MathUtilsDouble.ROOT3 };
    }
}

