using System;
using System.Collections;
using UnityEngine;

public class MathUtilsDouble : MonoBehaviour
{
    public static double ROOT3 = Math.Sqrt(3);

    /**
     * Converts geographic latitude and longitude coordinates to spherical coordinates on a sphere of radius 1.
     *
     * @param geo - geographic coordinates as a double array of length 2, {longitude, latitude}, in degrees
     * @return the corresponding spherical coordinates in radians: {longitude, colatitude}
     */
    public static double[] Geo2Spherical(double[] geo)
    {
        double lambda = geo[0] * Math.PI / 180;
        double phi = (90 - geo[1]) * Math.PI / 180;
        return new double[] { lambda, phi };
    }

    /**
     * Converts spherical coordinates to geographic coordinates on a sphere of radius 1.
     *
     * @param spherical - spherical coordinates in radians as a double array of length 2: {longitude, colatitude}
     * @return the corresponding geographic coordinates in degrees: {longitude, latitude}
     */
    public static double[] Spherical2Geo(double[] spherical)
    {
        double lon = spherical[0] * 180 / Math.PI;
        double lat = 90 - spherical[1] * 180 / Math.PI;
        return new double[] { lon, lat };
    }

    /**
     * Converts spherical coordinates to Cartesian coordinates on a sphere of radius 1.
     *
     * @param spherical - spherical coordinates in radians as a double array of length 2: {longitude, colatitude}
     * @return the corresponding Cartesian coordinates: {x, y, z}
     */
    public static double[] Spherical2Cartesian(double[] spherical)
    {
        double sinphi = Math.Sin(spherical[1]);
        double x = sinphi * Math.Cos(spherical[0]);
        double y = sinphi * Math.Sin(spherical[0]);
        double z = Math.Cos(spherical[1]);
        return new double[] { x, y, z };
    }

    /**
     * Converts Cartesian coordinates to spherical coordinates on a sphere of radius 1.
     *
     * @param cartesian - Cartesian coordinates as double array of length 3: {x, y, z}
     * @return the spherical coordinates of the corresponding normalized vector
     */
    public static double[] Cartesian2Spherical(double[] cartesian)
    {
        double lambda = Math.Atan2(cartesian[1], cartesian[0]);
        double phi = Math.Atan2(Math.Sqrt(cartesian[0] * cartesian[0] + cartesian[1] * cartesian[1]), cartesian[2]);
        return new double[] { lambda, phi };
    }

    /**
     * Multiples the given matrix with the given vector.
     * The matrix is assumed to be square and the vector is assumed to be of the same dimension as the matrix.
     *
     * @param matrix - the matrix as a n*n double array
     * @param vector - the vector as double array of length n
     * @return the result of the multiplication as an array of double on length n
     */

    /**
 * TODO produceZYZRotationMatrix javadoc
 *
 * @param a
 * @param b
 * @param c
 * @return
 */
    public static double[][] ProduceZYZRotationMatrix(double a, double b, double c)
    {

        double sina = Math.Sin(a);
        double cosa = Math.Cos(a);
        double sinb = Math.Sin(b);
        double cosb = Math.Cos(b);
        double sinc = Math.Sin(c);
        double cosc = Math.Cos(c);

        double[][] mat = new double[3][];
        for (int i = 0; i < mat.Length; i++)
        {
            mat[i] = new double[3];
        }
        mat[0][0] = cosa * cosb * cosc - sinc * sina;
        mat[0][1] = -sina * cosb * cosc - sinc * cosa;
        mat[0][2] = cosc * sinb;

        mat[1][0] = sinc * cosb * cosa + cosc * sina;
        mat[1][1] = cosc * cosa - sinc * cosb * sina;
        mat[1][2] = sinc * sinb;

        mat[2][0] = -sinb * cosa;
        mat[2][1] = sinb * sina;
        mat[2][2] = cosb;

        return mat;
    }

    public static double[] MatVecProdD(double[][] matrix, double[] vector)
    {
        double[] result = new double[vector.Length];
        for (int i = 0; i < result.Length; i++)
        {
            for (int j = 0; j < matrix[i].Length; j++)
            {
                result[i] += matrix[i][j] * vector[j];
            }
        }
        return result;
    }
}
