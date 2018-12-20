__kernel void adder(__global const unsigned int* img, __global unsigned int* result)
{
        int idx = get_global_id(0);
        int index = 256 * (idx % 3);
        atomic_add(&result[index + img[idx]], 1);
        //result[index + img[idx]]+=1;
}