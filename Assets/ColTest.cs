using UnityEngine;

public class Example : MonoBehaviour
{
    private float speed = 2f;

    //Moves this GameObject 2 units a second in the forward direction
    void Update()
    {
        transform.Translate(speed * Time.deltaTime * Vector3.forward);
    }

    //Upon collision with another GameObject, this GameObject will reverse direction
    private void OnTriggerEnter(Collider other)
    {
        speed *= -1;
    }
}