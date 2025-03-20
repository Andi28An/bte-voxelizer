using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using Unity.Burst;
using UnityEngine;
using UnityEngine.Jobs;
using System.Collections.Generic;
using CesiumForUnity;
using System;

public class VoxelizeColorJobs : MonoBehaviour
{
    #region Fields
    public GameObject voxel;
    public ComputeShader computeShader;
    [SerializeField] private float2 _pos1 = new(0, 0);
    [SerializeField] private float2 _pos2 = new(0, 0);
    [SerializeField] private int _areaMaxAltitude = 160;
    [SerializeField] private int _areaGeoidHeight = 35;
    [SerializeField] private int _jobBatchSize = 32;
    [SerializeField] private float _boxCastSize = 0.75f;

    private NativeList<VoxelData> _topVoxels;
    private NativeList<VoxelData> _bottomVoxels;
    private int _minimumHeight;
    private int _maximumHeight;
    private bool _isProcessing = false;
    #endregion

    // Rest of the structs remain the same...
    public struct VoxelData
    {
        public float3 Position;
        public float Height;
        public quaternion Rotation;

        public VoxelData(float3 position, float height, quaternion rotation)
        {
            Position = position;
            Height = height;
            Rotation = rotation;
        }
    }

    public struct RaycastResult
    {
        public bool hit;
        public float3 point;
        public float height;
    }

    [BurstCompile(Debug = true)]
    private struct TerrainSamplingJob : IJobParallelFor
    {
        [ReadOnly] public float2 StartPos;
        [ReadOnly] public float2 EndPos;
        [ReadOnly] public float StartHeight;
        [ReadOnly] public int BatchIndex;
        [ReadOnly] public int BatchSize;
        public NativeArray<RaycastCommand> RaycastCommands;

        public void Execute(int index)
        {
            int globalIndex = BatchIndex * BatchSize + index;
            int x = globalIndex / (int)(math.abs(EndPos.x - StartPos.x));
            int z = globalIndex % (int)(math.abs(EndPos.y - StartPos.y));

            float3 position = new float3(
                StartPos.x + x,
                StartHeight,
                StartPos.y + z
            );

            var parameters = new QueryParameters
            {
                hitMultipleFaces = false,
                layerMask = -1,
                hitTriggers = QueryTriggerInteraction.UseGlobal
            };

            RaycastCommands[index] = new RaycastCommand(position, Vector3.down, parameters, 1000f);
        }
    }

    [BurstCompile(Debug = true)]
    private struct ProcessTerrainDataJob : IJobParallelFor
    {
        [ReadOnly] public NativeArray<RaycastHit> RaycastHits;
        [ReadOnly] public float2 StartPos;
        [ReadOnly] public float2 EndPos;
        [ReadOnly] public int BatchIndex;
        [ReadOnly] public int BatchSize;
        public NativeList<RaycastResult>.ParallelWriter Results;

        public void Execute(int index)
        {
            var hit = RaycastHits[index];

            if (hit.distance > 0)
            {
                int globalIndex = BatchIndex * BatchSize + index;
                int x = globalIndex / (int)(math.abs(EndPos.x - StartPos.x));
                int z = globalIndex % (int)(math.abs(EndPos.y - StartPos.y));

                var result = new RaycastResult
                {
                    hit = true,
                    point = hit.point,
                    height = hit.point.y
                };

                Results.AddNoResize(result);
            }
        }
    }

    [BurstCompile(Debug = true)]
    private struct ProcessHeightBoundsJob : IJob
    {
        [ReadOnly] public NativeList<RaycastResult> Results;
        public NativeReference<int> MinHeight;
        public NativeReference<int> MaxHeight;

        public void Execute()
        {
            int minH = int.MaxValue;
            int maxH = int.MinValue;

            for (int i = 0; i < Results.Length; i++)
            {
                if (Results[i].hit)
                {
                    int height = (int)math.round(Results[i].height);
                    minH = math.min(minH, height);
                    maxH = math.max(maxH, height);
                }
            }

            MinHeight.Value = minH;
            MaxHeight.Value = maxH;
        }
    }

    private void Update()
    {
        if (Input.GetKeyDown(KeyCode.Space) && !_isProcessing)
        {
            StartVoxelization();
        }
    }

    private void StartVoxelization()
    {
        if (_isProcessing) return;
        _isProcessing = true;

        int width = (int)math.abs(_pos2.x - _pos1.x);
        int height = (int)math.abs(_pos2.y - _pos1.y);
        int totalVoxels = width * height;

        _topVoxels = new NativeList<VoxelData>(totalVoxels, Allocator.Persistent);
        _bottomVoxels = new NativeList<VoxelData>(totalVoxels * 2, Allocator.Persistent);

        StartCoroutine(ProcessVoxelsCoroutine(totalVoxels));
    }

    private System.Collections.IEnumerator ProcessVoxelsCoroutine(int totalVoxels)
    {
        var minHeight = new NativeReference<int>(Allocator.TempJob);
        var maxHeight = new NativeReference<int>(Allocator.TempJob);
        var allResults = new NativeList<RaycastResult>(totalVoxels, Allocator.TempJob);

        try
        {
            int numBatches = (totalVoxels + _jobBatchSize - 1) / _jobBatchSize;

            for (int batchIndex = 0; batchIndex < numBatches; batchIndex++)
            {
                int batchSize = math.min(_jobBatchSize, totalVoxels - batchIndex * _jobBatchSize);

                var raycastCommands = new NativeArray<RaycastCommand>(batchSize, Allocator.TempJob);
                var raycastHits = new NativeArray<RaycastHit>(batchSize, Allocator.TempJob);

                var samplingJob = new TerrainSamplingJob
                {
                    StartPos = _pos1,
                    EndPos = _pos2,
                    StartHeight = _areaMaxAltitude + _areaGeoidHeight,
                    BatchIndex = batchIndex,
                    BatchSize = _jobBatchSize,
                    RaycastCommands = raycastCommands
                };

                var samplingHandle = samplingJob.Schedule(batchSize, 32);
                var raycastHandle = RaycastCommand.ScheduleBatch(raycastCommands, raycastHits, 32, samplingHandle);

                var processJob = new ProcessTerrainDataJob
                {
                    RaycastHits = raycastHits,
                    StartPos = _pos1,
                    EndPos = _pos2,
                    BatchIndex = batchIndex,
                    BatchSize = _jobBatchSize,
                    Results = allResults.AsParallelWriter()
                };

                var processHandle = processJob.Schedule(batchSize, 32, raycastHandle);
                processHandle.Complete();

                raycastCommands.Dispose();
                raycastHits.Dispose();

                // Yield after each batch to prevent frame drops
                yield return null;
            }

            // Process height bounds
            var heightBoundsJob = new ProcessHeightBoundsJob
            {
                Results = allResults,
                MinHeight = minHeight,
                MaxHeight = maxHeight
            };

            var heightBoundsHandle = heightBoundsJob.Schedule();
            heightBoundsHandle.Complete();

            _minimumHeight = minHeight.Value;
            _maximumHeight = maxHeight.Value;

            // Process the results
            ProcessResults(allResults);
            WriteResultsToFile();
        }
        finally
        {
            // Cleanup
            if (minHeight.IsCreated) minHeight.Dispose();
            if (maxHeight.IsCreated) maxHeight.Dispose();
            if (allResults.IsCreated) allResults.Dispose();
            _isProcessing = false;
        }
    }

    private void ProcessResults(NativeList<RaycastResult> results)
    {
        var testCube = Instantiate(voxel, GameObject.Find("CesiumGeoreference").GetComponent<Transform>());
        testCube.AddComponent<CesiumGlobeAnchor>();

        foreach (var result in results)
        {
            if (!result.hit) continue;

            var voxelData = new VoxelData(result.point, result.height, Quaternion.identity);
            _topVoxels.Add(voxelData);

            double[] geo = CoordConversionDouble.ToGeo(
                voxelData.Position.x + 0.5f,
                voxelData.Position.z + 0.5f
            );

            testCube.GetComponent<CesiumGlobeAnchor>().longitudeLatitudeHeight =
                new double3(geo[0], geo[1], voxelData.Position.y);

            // Check for bottom voxels
            var currentPos = testCube.transform.position;
            var currentRot = testCube.transform.rotation;

            for (int y = (int)voxelData.Position.y - 1; y >= _minimumHeight; y--)
            {
                var checkPos = new Vector3(currentPos.x, y, currentPos.z);
                if (Physics.CheckBox(checkPos, new Vector3(_boxCastSize, 0.5f, _boxCastSize), currentRot))
                {
                    _bottomVoxels.Add(new VoxelData(
                        new float3(voxelData.Position.x, y, voxelData.Position.z),
                        y,
                        currentRot
                    ));
                }
            }
        }

        Destroy(testCube);
    }

    private void WriteResultsToFile()
    {
        string path = "output.txt";
        if (System.IO.File.Exists(path))
        {
            System.IO.File.Delete(path);
        }

        using (var writer = new System.IO.StreamWriter(path, true))
        {
            foreach (var voxel in _topVoxels)
            {
                int x = (int)(voxel.Position.x - _pos1.x);
                int y = (int)(voxel.Position.y - 2);
                int z = (int)(voxel.Position.z - _pos1.y);
                writer.WriteLine($"{x}, {y}, {z}");
            }

            foreach (var voxel in _bottomVoxels)
            {
                int x = (int)(voxel.Position.x - _pos1.x);
                int y = (int)(voxel.Position.y - 2);
                int z = (int)(voxel.Position.z - _pos1.y);
                writer.WriteLine($"{x}, {y}, {z}");
            }
        }
    }

    private void OnDisable()
    {
        if (_topVoxels.IsCreated) _topVoxels.Dispose();
        if (_bottomVoxels.IsCreated) _bottomVoxels.Dispose();
    }
}