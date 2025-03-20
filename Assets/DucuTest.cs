using CesiumForUnity;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Runtime.CompilerServices;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;
using static Unity.VisualScripting.Member;
using static UnityEditor.ShaderGraph.Internal.KeywordDependentCollection;

public class DucuTest : MonoBehaviour
{
    #region Fields

    public GameObject voxel;
    public ComputeShader computeShader;
    [SerializeField]
    private float2 _pos1 = new(0, 0);
    [SerializeField]
    private float2 _pos2 = new(0, 0);
    [SerializeField]
    private int _areaMaxAltitude = 160; // <-- safe altitude to start scanning the area from
    [SerializeField]
    private int _areaGeoidHeight = 35; // <-- geoid height in the area

    #endregion

    public float2 pos1
    {
        get => _pos1;
        set => _pos1 = value;
    }
    public float2 pos2
    {
        get => _pos2;
        set => _pos2 = value;
    }
    public int areaMaxAltitude
    {
        get => _areaMaxAltitude;
        set => _areaMaxAltitude = value;
    }
    public int areaGeoidHeight
    {
        get => _areaGeoidHeight;
        set => _areaGeoidHeight = value;
    }

    private void Update()
    {
        if (Input.GetKeyDown(KeyCode.Space))
        {
            AreaTest3();
        }
    }

    private void CoordTest()
    {
        float x = 4294245.516f;
        float z = -4458512.401f; // aici pui varfu primariei sau idk xd
        float lon = 21.930719999771775f;
        float lat = 47.057053001891255f;//astea is expected values
        float[] geo = CoordConversion.ToGeo(x, z);
        GameObject instance = Instantiate(new GameObject("Test CLONE"), GameObject.Find("CesiumGeoreference").GetComponent<Transform>());
        instance.AddComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geo[0], geo[1], 207);
    }

    private void AreaTest1()
    {
        float[] pos1 = new float[] { 4580409.50f, -4107347.50f };
        float[] pos2 = new float[] { 4580422.50f, -4107332.50f };
        int index1 = 0;
        int index2 = 0;

        int areaMaxAltitude = 160; // <-- safe altitude to start scanning the area from

        int areaGeoidHeight = 35; // <-- geoid height in the area 

        int startHeight = areaGeoidHeight + areaMaxAltitude;

        double originHeight = 167.8; // <-- this is the height of the first point in the area

        float originYRotation = -9.5f; // <-- this is the rotation of the first point in the area

        //instantiate a Collision test cube
        GameObject collisionTestCube = Instantiate(voxel, GameObject.Find("CesiumGeoreference").GetComponent<Transform>());
        collisionTestCube.AddComponent<CesiumGlobeAnchor>();
        collisionTestCube.name = "Collision Test Cube";

        //calculate the angle of the generated grid
        double[] geoAngle1 = CoordConversionDouble.ToGeo(pos1[0], pos1[1]);
        collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geoAngle1[0], geoAngle1[1], startHeight);
        Vector3 pos1Vector3 = collisionTestCube.transform.position;
        double[] geoAngle2 = CoordConversionDouble.ToGeo(pos2[0], pos1[1]);
        collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geoAngle2[0], geoAngle2[1], startHeight);
        Vector3 pos2Vector3 = collisionTestCube.transform.position;
        float angle = Vector3.Angle(pos1Vector3, pos2Vector3);

        originYRotation = angle;
        Destroy(collisionTestCube);

        ////set cesium georeference origin to the first point in the area and hardcoded height
        //float[] geoOrigin = CoordConversion.ToGeo(pos1[0], pos1[1]);
        //GameObject.Find("CesiumGeoreference").GetComponent<CesiumGeoreference>().longitude = geoOrigin[0];
        //GameObject.Find("CesiumGeoreference").GetComponent<CesiumGeoreference>().latitude = geoOrigin[1];
        //GameObject.Find("CesiumGeoreference").GetComponent<CesiumGeoreference>().height = originHeight;

        ////move dynamic camera to cesium georeference
        //GameObject.Find("DynamicCamera").GetComponent<Transform>().position = GameObject.Find("CesiumGeoreference").GetComponent<Transform>().position;

        //for each point in the area spawn a cube at that position and set its height with linecast
        for (float i = pos1[0]; i < pos2[0]; i += 1)
        {   
            index2 = 0;
            for (float j = pos1[1]; j < pos2[1]; j += 1)
            {
                //double method
                double[] geo2 = CoordConversionDouble.ToGeo(i, j);
                GameObject instance2 = Instantiate(voxel, GameObject.Find("CesiumGeoreference").GetComponent<Transform>());
                instance2.tag = "Voxel";
                instance2.AddComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geo2[0], geo2[1], 210);
                instance2.AddComponent<SchematicCoords>();
                instance2.GetComponent<SchematicCoords>().SetX(index1);
                instance2.GetComponent<SchematicCoords>().SetZ(index2);
                //set layer to voxel
                instance2.layer = LayerMask.NameToLayer("Voxel");

                index2++;
            }
            index1++;
        }

        int minimumHeight = 10000;

        //for each cube in the area set its height with linecast
        foreach (GameObject voxel in GameObject.FindGameObjectsWithTag("Voxel"))
        {
            if (Physics.Raycast(voxel.transform.position, Vector3.down, out RaycastHit hit, 1000))
            {
                //round the height to the nearest whole
                hit.point = new Vector3(hit.point.x, Mathf.Round(hit.point.y), hit.point.z);
                voxel.transform.position = hit.point;
                //set gameobject SchematicCoords to the rounded height
                voxel.GetComponent<SchematicCoords>().SetY((int)Mathf.Round(hit.point.y));
                voxel.name = "Voxel x:" + voxel.GetComponent<SchematicCoords>().GetX() + ", y: " + voxel.GetComponent<SchematicCoords>().GetY() + ", z:" + voxel.GetComponent<SchematicCoords>().GetZ();

                //if the height is lower than the minimum height set the minimum height to the height of the current voxel
                if (voxel.GetComponent<SchematicCoords>().GetY() < minimumHeight)
                {
                    minimumHeight = voxel.GetComponent<SchematicCoords>().GetY();
                }
            }
        }

        foreach (GameObject topVoxel in GameObject.FindGameObjectsWithTag("Voxel"))
        {
            int topVoxelX = topVoxel.GetComponent<SchematicCoords>().GetX();
            int topVoxelY = topVoxel.GetComponent<SchematicCoords>().GetY();
            int topVoxelZ = topVoxel.GetComponent<SchematicCoords>().GetZ();

            Vector3 posTopVoxel = topVoxel.transform.position;
            Quaternion rotTopVoxel = topVoxel.transform.rotation;

            int currentY = topVoxelY;
            while(currentY >= minimumHeight)
            {
                currentY--;
                Vector3 posNewVoxel = new(posTopVoxel.x, currentY, posTopVoxel.z);
                //if there is a voxel at the current position determined with BoxOverlap instantiate a voxel at the current position
                Collider[] collisions = Physics.OverlapBox(posNewVoxel, new Vector3(0.5f, 0.5f, 0.5f));
                if(collisions.Length>0 )
                {
                    GameObject newVoxel = Instantiate(voxel, GameObject.Find("CesiumGeoreference").GetComponent<Transform>());
                    newVoxel.tag = "Voxel";
                    newVoxel.transform.position = posNewVoxel;
                    newVoxel.AddComponent<SchematicCoords>();
                    newVoxel.GetComponent<SchematicCoords>().SetX(topVoxelX);
                    newVoxel.GetComponent<SchematicCoords>().SetY(currentY);
                    newVoxel.GetComponent<SchematicCoords>().SetZ(topVoxelZ);
                    //set layer to voxel
                    newVoxel.layer = LayerMask.NameToLayer("Voxel2");
                    //color it blue
                    newVoxel.GetComponent<Renderer>().material.color = UnityEngine.Color.blue;
                }
            }
        }
        string path = @"output.txt";
        if (System.IO.File.Exists(path))
        {
            System.IO.File.Delete(path);
        }

        //for each cube in the area output its x, y, z into a txt file
        foreach (GameObject voxel in GameObject.FindGameObjectsWithTag("Voxel"))
        {
            System.IO.File.AppendAllText(@"output.txt", voxel.GetComponent<SchematicCoords>().GetX() + ", " + voxel.GetComponent<SchematicCoords>().GetY() + ", " + voxel.GetComponent<SchematicCoords>().GetZ() + "\n");
        }
    }

    //New method for area test avoiding instancing
    private void AreaTest2()
    {
        float[] pos1 = new float[] { 4580409.50f, -4107347.50f };
        float[] pos2 = new float[] { 4580422.50f, -4107332.50f };

        int areaMaxAltitude = 160; // <-- safe altitude to start scanning the area from

        int areaGeoidHeight = 35; // <-- geoid height in the area 

        int startHeight = areaGeoidHeight + areaMaxAltitude;

        double originHeight = 167.8; // <-- this is the height of the first point in the area

        //instantiate a Collision test cube
        GameObject collisionTestCube = Instantiate(voxel, GameObject.Find("CesiumGeoreference").GetComponent<Transform>());
        collisionTestCube.AddComponent<CesiumGlobeAnchor>();
        collisionTestCube.name = "Collision Test Cube";

        //calculate the angle of the generated grid
        double[] geoAngle1 = CoordConversionDouble.ToGeo(pos1[0], pos1[1]);
        collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geoAngle1[0], geoAngle1[1], startHeight);
        Vector3 pos1Vector3 = collisionTestCube.transform.position;
        double[] geoAngle2 = CoordConversionDouble.ToGeo(pos2[0], pos1[1]);
        collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geoAngle2[0], geoAngle2[1], startHeight);
        Vector3 pos2Vector3 = collisionTestCube.transform.position;
        float angle = Vector3.Angle(pos1Vector3, pos2Vector3);

        float originYRotation = angle;
        Quaternion originRotation = new(0, originYRotation, 0, 0);
        collisionTestCube.transform.rotation = new Quaternion(0, originYRotation, 0, 0);

        //set cesium georeference origin to the first point in the area and hardcoded height
        float[] geoOrigin = CoordConversion.ToGeo(pos1[0], pos1[1]);
        GameObject.Find("CesiumGeoreference").GetComponent<CesiumGeoreference>().longitude = geoOrigin[0];
        GameObject.Find("CesiumGeoreference").GetComponent<CesiumGeoreference>().latitude = geoOrigin[1];
        GameObject.Find("CesiumGeoreference").GetComponent<CesiumGeoreference>().height = originHeight;

        //move dynamic camera to cesium georeference 10 meters above
        GameObject.Find("DynamicCamera").GetComponent<Transform>().position = GameObject.Find("CesiumGeoreference").GetComponent<Transform>().position + new Vector3(0, 10, 0);

        //map each point in the area to a geo coordinate
        Dictionary<Tuple<float, float>, Tuple<double, double>> points = new();
        for (float i = pos1[0]; i < pos2[0]; i += 1)
        {
            for (float j = pos1[1]; j < pos2[1]; j += 1)
            {
                Tuple<float, float> pos = new((float)i, (float)j);
                double[] geo = CoordConversionDouble.ToGeo(i, j);
                points.Add(pos, new Tuple<double, double>(geo[0], geo[1]));
            }
        }

        int startingHeight = areaGeoidHeight + areaMaxAltitude;

        Dictionary<Tuple<float, float>, int> maxHeights = new();

        int minFoundHeight = 10000;
        int maxFoundHeight = -1000;
        int schemHeight = 0;
        
        //for each point in the area find the height with linecast
        for(float i = pos1[0]; i < pos2[0]; i+=1)
        {
            for(float j = pos1[1]; j < pos2[1]; j+=1)
            {
                //get the current key and value
                Tuple<float, float> pos = new(i, j);
                double[] geo = new double[] { points[pos].Item1, points[pos].Item2 };

                //set collision test cube position to the current geo coordinate
                collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geo[0], geo[1], startingHeight);

                //linecast down from the collision test cube
                if (Physics.Raycast(collisionTestCube.transform.position, Vector3.down, out RaycastHit hit, 1000))
                {
                    //round the height to the nearest whole
                    hit.point = new Vector3(hit.point.x, Mathf.Round(hit.point.y), hit.point.z);
                    //set the height of the current point to the rounded height
                    int height = (int)Mathf.Round(hit.point.y);
                    //add the current point and its height to the dictionary
                    maxHeights.Add(pos, height);

                    //update schemHeight and minFoundHeight
                    if (i == 0 && j == 0)
                    {
                        schemHeight = height;
                    }
                    if (height > maxFoundHeight)
                    {
                        maxFoundHeight = height;
                    }

                    if (height < minFoundHeight)
                    {
                        minFoundHeight = height;
                    }
                }
            }
        }

        //create a dynamic array to store found voxel positions 
        List<int[]> schemList = new();

        //add top voxels to list
        foreach (KeyValuePair <Tuple<float, float>, int> entry in maxHeights)
        {
            float[] pos = new float[] { entry.Key.Item1, entry.Key.Item2 };
            int height = entry.Value;

            int[] schemPos = new int[] { (int)(pos[0] - pos1[0]), height - schemHeight, (int)(pos[1] - pos1[1]) };
            schemList.Add(schemPos);
            Debug.Log("Top Voxel Added at x: " + pos[0] + ", y: " + height + ", z: " + pos[1]);
        }

        //for each point in the area traverse down from its height
        for (float i = pos1[0]; i < pos2[0]; i += 1)
        {
            for (float j = pos1[1]; j < pos2[1]; j += 1)
            {
                //get the current key and value
                Tuple<float, float> pos = new(i, j);
                double[] geo = new double[] { points[pos].Item1, points[pos].Item2 };
                int currentHeight = maxFoundHeight;

                collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geo[0], geo[1], startingHeight);
                Vector3 currentPosition = new(collisionTestCube.transform.position.x, currentHeight, collisionTestCube.transform.position.z);
                Vector3 collisionTestCubePos = collisionTestCube.transform.position;

                while (currentHeight >= minFoundHeight)
                {
                    currentHeight--;

                    // set collision test cube position to the current geo coordinate
                    collisionTestCube.transform.position.Set(collisionTestCubePos.x, currentHeight, collisionTestCubePos.z);

                    //boxcast at the collision test cube position
                    Collider[] collisions = Physics.OverlapBox(collisionTestCube.transform.position, new Vector3(0.5f, 0.5f, 0.5f), originRotation);
                    if(collisions.Length>0)
                    { 
                        //add the current position to the dynamic array
                        int[] schemPos = new int[] { (int)(i - pos1[0]), currentHeight - schemHeight, (int)(j - pos1[1]) };
                        Debug.Log("Added voxel at x: " + i + ", y: " + currentHeight + ", z: " + j);
                        schemList.Add(schemPos);
                    }
                }
            }
        }

        string path = @"output.txt";
        if (System.IO.File.Exists(path))
        {
            System.IO.File.Delete(path);
        }

        //for each point in the dynamic array output its x, y, z into a txt file
        foreach (int[] pos in schemList)
        {
            System.IO.File.AppendAllText(path, Mathf.Floor(pos[0]) + ", " + Mathf.Floor(pos[1]) + ", " + Mathf.Floor(pos[2]) + "\n");
        }
    }

    void AreaTest3()
    {
        float[] pos1 = new float[] { _pos1.x, _pos1.y };
        float[] pos2 = new float[] { _pos2.x, _pos2.y };

        int areaMaxAltitude = _areaMaxAltitude; // <-- safe altitude to start scanning the area from

        int areaGeoidHeight = _areaGeoidHeight; // <-- geoid height in the area

        int startHeight = areaGeoidHeight + areaMaxAltitude;

        //instantiate a Collision test cube
        GameObject collisionTestCube = Instantiate(voxel, GameObject.Find("CesiumGeoreference").GetComponent<Transform>());
        collisionTestCube.AddComponent<CesiumGlobeAnchor>();
        collisionTestCube.name = "Collision Test Cube";

        //calculate the angle of the generated grid
        double[] geoAngle1 = CoordConversionDouble.ToGeo(pos1[0] + 0.5f, pos1[1] + 0.5f);
        collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geoAngle1[0], geoAngle1[1], startHeight);
        Vector3 pos1Vector3 = collisionTestCube.transform.position;
        double[] geoAngle2 = CoordConversionDouble.ToGeo(pos2[0] + 0.5f, pos1[1] + 0.5f);
        collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geoAngle2[0], geoAngle2[1], startHeight);
        Vector3 pos2Vector3 = collisionTestCube.transform.position;
        float angle = Vector3.Angle(pos1Vector3, pos2Vector3);

        float originYRotation = angle;
        collisionTestCube.transform.rotation = new Quaternion(0, originYRotation, 0, 0);

        List<VoxelData> TopDataList = new();

        int minimumHeight = 10000;
        int maximumHeight = -10000;

        int index1 = 0;
        //for each point in the area spawn a cube at that position and set its height with linecast
        for (float i = pos1[0]; i <= pos2[0]; i += 1)
        {
            int index2 = 0;
            for (float j = pos1[1]; j <= pos2[1]; j += 1)
            {
                double[] geo = CoordConversionDouble.ToGeo(i + 0.5f, j + 0.5f);
                collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geo[0], geo[1], startHeight);
                if (Physics.Raycast(collisionTestCube.transform.position, Vector3.down, out RaycastHit hit, 1000))
                {
                    //round the height to the nearest whole
                    int height = (int)Mathf.Round(hit.point.y);
                    //set gameobject SchematicCoords to the rounded height
                    VoxelData voxelData = new VoxelData(i, height, j);
                    TopDataList.Add(voxelData);
                    Debug.Log("Top Voxel Added at x: " + i + ", y: " + height + ", z: " + j + "\n" + "Coords: " + geo[0] + ", " + geo[1]);

                    //if the height is lower than the minimum height set the minimum height to the height of the current voxel
                    if (height < minimumHeight)
                    {
                        minimumHeight = height;
                    }
                    if (height > maximumHeight)
                    {
                        maximumHeight = height;
                    }
                }
                index2++;
            }
            index1++;
        }

        List<VoxelData> BottomDataList = new();

        foreach (VoxelData voxelData in TopDataList)
        {
            double[] geo = CoordConversionDouble.ToGeo(voxelData.GetX() + 0.5f, voxelData.GetZ() + 0.5f);

            collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geo[0], geo[1], maximumHeight);

            Vector3 posTopVoxel = collisionTestCube.transform.position;
            Quaternion rotTopVoxel = collisionTestCube.transform.rotation;

            int currentY = voxelData.GetY();
            while (currentY >= minimumHeight)
            {
                currentY--;
                Vector3 posNewVoxel = new(posTopVoxel.x, currentY, posTopVoxel.z);
                //if there is a voxel at the current position determined with BoxOverlap 
                Collider[] collisions = Physics.OverlapBox(posNewVoxel, new Vector3(0.75f, 0.5f, 0.75f), rotTopVoxel);
                if (collisions.Length > 0)
                {
                    VoxelData newVoxelData = new VoxelData(voxelData.GetX(), currentY, voxelData.GetZ());
                    BottomDataList.Add(newVoxelData);
                    Debug.Log("Bottom Voxel Added at x: " + voxelData.GetX() + ", y: " + currentY + ", z: " + voxelData.GetZ());
                }
            }
        }

        string path = @"output.txt";
        if(System.IO.File.Exists(path))
        {
            System.IO.File.Delete(path);
        }
        else
        {
            System.IO.File.Create(path);
        }

        //for each cube in the area output its x, y, z into a txt file
        foreach (VoxelData voxelData in TopDataList)
        {
            int x = (int)(voxelData.GetX() - pos1[0]);
            int y = voxelData.GetY()-2; //hardcoded  -2
            int z = (int)(voxelData.GetZ() - pos1[1]);
            System.IO.File.AppendAllText(path, x + ", " + y + ", " + z + "\n");
        }
        foreach (VoxelData voxelData in BottomDataList)
        {
            int x = (int)(voxelData.GetX() - pos1[0]);
            int y = voxelData.GetY()-2; //hardcoded -2
            int z = (int)(voxelData.GetZ() - pos1[1]);
            System.IO.File.AppendAllText(path, x + ", " + y + ", " + z + "\n");
        }
    }

    void AreaTest3Color()
    {

        float[] pos1 = new float[] { _pos1.x, _pos1.y};
        float[] pos2 = new float[] { _pos2.x, _pos2.y};

        Debug.Log("pos1: " + pos1[0] + ", " + pos1[1]);
        Debug.Log("pos2: " + pos2[0] + ", " + pos2[1]);

        //int areaMaxAltitude = 200; // <-- safe altitude to start scanning the area from

        //int areaGeoidHeight = 43; // <-- geoid height in the area 

        int startHeight = _areaGeoidHeight + _areaMaxAltitude;

        //instantiate a Collision test cube
        GameObject collisionTestCube = Instantiate(voxel, GameObject.Find("CesiumGeoreference").GetComponent<Transform>());
        collisionTestCube.AddComponent<CesiumGlobeAnchor>();
        collisionTestCube.name = "Collision Test Cube";

        //calculate the angle of the generated grid
        double[] geoAngle1 = CoordConversionDouble.ToGeo(pos1[0] + 0.5f, pos1[1] + 0.5f);
        collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geoAngle1[0], geoAngle1[1], startHeight);
        Vector3 pos1Vector3 = collisionTestCube.transform.position;
        double[] geoAngle2 = CoordConversionDouble.ToGeo(pos2[0] + 0.5f, pos1[1] + 0.5f);
        collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geoAngle2[0], geoAngle2[1], startHeight);
        Vector3 pos2Vector3 = collisionTestCube.transform.position;
        float angle = Vector3.Angle(pos1Vector3, pos2Vector3);

        float originYRotation = angle;
        collisionTestCube.transform.rotation = new Quaternion(0, originYRotation, 0, 0);

        List<VoxelData> TopDataList = new();

        //compute shader initialization
        ComputeBuffer buffer;

        int minimumHeight = 10000;

        int index1 = 0;
        //for each point in the area spawn a cube at that position and set its height with linecast
        for (float i = pos1[0]; i <= pos2[0]; i += 1)
        {
            int index2 = 0;
            for (float j = pos1[1]; j <= pos2[1]; j += 1)
            {
                double[] geo = CoordConversionDouble.ToGeo(i + 0.5f, j + 0.5f);
                collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geo[0], geo[1], startHeight);
                if (Physics.Raycast(collisionTestCube.transform.position, Vector3.down, out RaycastHit hit, 1000))
                {
                    //if collision is with DynamicCamera ignore it
                    if (hit.collider.gameObject.name == "DynamicCamera")
                    {
                        continue;
                    }

                    //round the height to the nearest whole
                    int height = (int)Mathf.Round(hit.point.y);
                    
                    //get uv position from hit point
                    Vector2 uv = hit.textureCoord;
                    //get pixel color from uv position
                    Texture2D source = hit.collider.gameObject.GetComponent<MeshRenderer>().material.GetTexture("_baseColorTexture") as Texture2D;
                    //scale uv to pixel position
                    uv.x *= source.width;
                    uv.y *= source.height;

                    int kernelHandle = computeShader.FindKernel("CSMain");
                    buffer = new ComputeBuffer(1, 3*sizeof(float));
                    float[] readout = new float[3];
                    computeShader.SetBuffer(kernelHandle, "outputBuffer", buffer);
                    computeShader.SetTexture(kernelHandle, "inputTexture", source);
                    computeShader.SetFloat("uv_x", uv.x);
                    computeShader.SetFloat("uv_y", uv.y);

                    computeShader.Dispatch(kernelHandle, 1, 1, 1);

                    buffer.GetData(readout);
                    buffer.Release();

                    Debug.Log("Color at x: " + i + ", y: " + height + ", z: " + j + " is " + "r: " + readout[0] + ", g: " + readout[1] + ", b:" + readout[2]);

                    //set gameobject SchematicCoords to the rounded height
                    VoxelData voxelData = new VoxelData(i, height, j);
                    TopDataList.Add(voxelData);
                    //Debug.Log("Top Voxel Added at x: " + i + ", y: " + height + ", z: " + j);

                    //if the height is lower than the minimum height set the minimum height to the height of the current voxel
                    if (height < minimumHeight)
                    {
                        minimumHeight = height;
                    }
                }

                index2++;
            }
            index1++;
        }

        List<VoxelData> BottomDataList = new();

        foreach (VoxelData voxelData in TopDataList)
        {
            double[] geo = CoordConversionDouble.ToGeo(voxelData.GetX() + 0.5f, voxelData.GetZ() + 0.5f);

            collisionTestCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight = new Unity.Mathematics.double3(geo[0], geo[1], startHeight);

            Vector3 posTopVoxel = collisionTestCube.transform.position;
            Quaternion rotTopVoxel = collisionTestCube.transform.rotation;

            int currentY = voxelData.GetY();
            while (currentY >= minimumHeight)
            {
                currentY--;
                Vector3 posNewVoxel = new(posTopVoxel.x, currentY, posTopVoxel.z);
                //if there is a voxel at the current position determined with BoxOverlap 
                Collider[] collisions = Physics.OverlapBox(posNewVoxel, new Vector3(0.7f, 0.5f, 0.7f), rotTopVoxel);
                if (collisions.Length > 0)
                {
                    //get average color at current position form mesh material shader
                    
                    //Debug.Log("Color at x: " + voxelData.GetX() + ", y: " + currentY + ", z: " + voxelData.GetZ() + " is " + color);

                    VoxelData newVoxelData = new VoxelData(voxelData.GetX(), currentY, voxelData.GetZ());
                    BottomDataList.Add(newVoxelData);
                    //Debug.Log("Bottom Voxel Added at x: " + voxelData.GetX() + ", y: " + currentY + ", z: " + voxelData.GetZ());
                }
            }
        }

        string path = @"output.txt";
        if (System.IO.File.Exists(path))
        {
            System.IO.File.Delete(path);
        }

        //for each cube in the area output its x, y, z into a txt file
        foreach (VoxelData voxelData in TopDataList)
        {
            int x = (int)(voxelData.GetX() - pos1[0]);
            int y = voxelData.GetY();
            int z = (int)(voxelData.GetZ() - pos1[1]);
            System.IO.File.AppendAllText(path, x + ", " + y + ", " + z + "\n");
        }
        foreach (VoxelData voxelData in BottomDataList)
        {
            int x = (int)(voxelData.GetX() - pos1[0]);
            int y = voxelData.GetY();
            int z = (int)(voxelData.GetZ() - pos1[1]);
            System.IO.File.AppendAllText(path, x + ", " + y + ", " + z + "\n");
        }
    }

    public class VoxelData
    {
        public float x;
        public int y;
        public float z;
        public VoxelData(float x, int y, float z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
        public float GetX()
        {
            return x;
        }
        public int GetY()
        {
            return y;
        }
        public float GetZ()
        {
            return z;
        }
        public void SetX(float x)
        {
               this.x = x;
        }
        public void SetY(int y)
        {
            this.y = y;
        }
        public void SetZ(float z)
        {
            this.z = z;
        }
    }
}
