using UnityEngine.Events;
using UnityEngine;

public class SchematicCoords : MonoBehaviour
{

    int _x = 0;
    int _y = 0;
    int _z = 0;
    
    void Start()
    {
        _x = 0;
        _y = 0;
        _z = 0;
    }

    public void SetX(int x)
    {
        _x = x;
    }
    public void SetY(int y)
    {
        _y = y;
    }
    public void SetZ(int z)
    {
        _z = z;
    }
    public int GetX()
    {
        return _x;
    }
    public int GetY()
    {
        return _y;
    }
    public int GetZ()
    {
        return _z;
    }
}