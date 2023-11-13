#include <vector>
#include <iostream>

using std::vector;

template <typename T>
std::ostream& operator << (std::ostream& os, const vector<T>& vector)
{
    // Printing all the elements using <<
    int len = vector.size();
    for (int index = 0; index < len - 1; index++) 
    {
        os << vector[index] << " ";
    }
    if (len != 0)
    {
        os << vector[len-1];
    } 
    return os;
}

template <typename T> class table
{
    private:
        vector<T> __data;
        int __n_dimensions;
        int __n_values;
        vector<int> __size;

    public:
        table(int x_dim)
        {
            __n_dimensions = 1;
            __n_values = x_dim;
            __size.resize(__n_dimensions);
            __data.resize(__n_values);
            __size[0] = x_dim;
        }

        table(int x_dim, int y_dim)
        {
            __n_dimensions = 2;
            __size.resize(__n_dimensions);
            __n_values = x_dim * y_dim;
            __data.resize(__n_values);
            __size[0] = x_dim;
            __size[1] = y_dim;
        }

        table(int x_dim, int y_dim, int z_dim)
        {
            __n_dimensions = 3;
            __size.resize(__n_dimensions);
            __n_values = x_dim * y_dim * z_dim;
            __data.resize(__n_values);
            __size[0] = x_dim;
            __size[1] = y_dim;
            __size[1] = z_dim;
        }

        T& operator()(int x)
        {
            if (__n_dimensions != 1)
            {
                std::cout << "incorrect size\n";
                exit(1);
            }
            if (x >= __size[0])
            {
                std::cout << "Segementation_fault\n";
                exit(2);
            }
            return (__data[x]);
        }

        T& operator()(int x, int y)
        {
            if (__n_dimensions != 2)
            {
                std::cout << "incorrect size\n";
                exit(1);
            }
            if ((x >= __size[0])||(y >= __size[1]))
            {
                std::cout << "Segementation_fault\n";
                exit(2);
            }
            return (__data[x*__size[1] + y]);
        }

        T& operator()(int x, int y, int z)
        {
            if (__n_dimensions != 3)
            {
                std::cout << "incorrect size\n";
                exit(1);
            }
            if ((x >= __size[0])||(y >= __size[1])||(z >= __size[2]))
            {
                std::cout << "Segementation_fault\n";
                exit(2);
            }
            return (__data[x*__size[2]*__size[1] + y*__size[2] + z]);
        }

        T& operator () (vector<int> coordinates)
        {
            int dim = coordinates.size();
            if (__n_dimensions != dim)
            {
                std::cout << "incorrect size\n";
                exit(1);
            }
            for (int i = 0; i < dim; i++)
            {
                if (coordinates[i] >= __size[i])
                {
                    std::cout << "Segementation_fault\n";
                    exit(2);
                }
            }
            int index = 0;
            int pre_factor = __n_values;
            for (int i = 0; i < dim; i++)
            {
                if (pre_factor % __size[i] != 0)
                {
                    std::cout << "Something went wrong with indexing\n";
                    exit(3);
                }
                pre_factor = pre_factor / __size[i];
                index += coordinates[i]*pre_factor;
            }
            return (__data[index]);
        }

        void operator = (table<T> the_other)
        {
            if (__n_dimensions != the_other.__n_dimensions)
            {
                std::cout << "Incompatible assignment\n";
                exit(4);
            }
            if (__size != the_other.__size)
            {
                std::cout << "Incompatible assignment\n";
                exit(4);
            }
            __data = the_other.__data;
        }

        bool operator == (table<T> the_other)
        {
            if (__n_dimensions != the_other.__n_dimensions)
            {
                return false;
            }
            if (__size != the_other.__size)
            {
                return false;
            }
            return (__data == the_other.__data);
        }

        void print()
        {
            int index = 0;
            while (index < __n_values)
            {
                for (int i = 0; i < __size[__n_dimensions-1]; i++)
                {
                    std::cout << __data[index] << " ";
                    index++;
                }
                std::cout << std::endl;
            }
        }

};

/*
int main()
{
    vector<int> v = {1, 2, 3};
    std::cout << v << std::endl;
    return 0;
}
*/