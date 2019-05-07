#include <Python.h>
int main(int argc, char** argv)
{
    Py_Initialize();
    if(!Py_Initialize())
    {
        printf"false in init";
        return -1;
    }


}
