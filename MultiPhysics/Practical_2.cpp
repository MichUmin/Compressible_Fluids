#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include  <string>

using std::min;
using std::cout;
using std::cin;
using std::vector;
using std::string;

const double pi = 3.14159265358979311600;

int nPoints = 400;
int padding = 3;
double x0 = 0;
double x1 = 1;
double tStart = 0.0;
double tStop;
double dx = (x1 - x0)/nPoints;
double C = 0.8;
vector<double> gas_coef = {1.4, 1.4};
int num_dimmensions = 1;
int num_variables = 2 + num_dimmensions;
int step = 0;
double v_l, v_r, rho_l, rho_r, p_l, p_r;
bool use_halfstep = true;

vector<vector<double>> u1; // u[0] and u[nPoints+1] are padding
vector<vector<double>> u2; // u holds rho, rho_v, E

vector<vector<double>> f; //f[i] holds the flux from u[i] to u[i+1]
vector<vector<double>> limiters;
vector<vector<double>> u_minus;
vector<vector<double>> u_plus;
vector<double> level_set;

vector<double> conservative_to_primitive(int, const vector<double> &);

vector<double> primitive_to_conservative(int, const vector<double> &);

void boundary_condition(vector<vector<double>> &values)
{
    for (int i = 0; i < padding; i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            values[i][j] = values[padding][j];
            values[padding + nPoints + i][j] = values[padding + nPoints - 1][j];
        }
    }
}

vector<double> (flux)(int, const vector<double> &);

vector<double> (*nummerical_flux)(int, const vector<double>&, const vector<double>&, double, double);

vector<double> Lax_Friedrichs_flux(int material, const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> result;
    result.resize(num_variables);
    vector<double> f_left = flux(material, velocity);
    vector<double> f_right = flux(material, velocity_next);
    for (int i = 0; i < num_variables; i++)
    {
        result[i] = (d_x / d_t)*(velocity[i] - velocity_next[i]);
        result[i] += f_left[i];
        result[i] += f_right[i];
        result[i] = result[i] / 2.0;
    }
    return result;
}

vector<double> Richtmyer_flux(int material, const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> velocity_half_step;
    velocity_half_step.resize(num_variables);
    vector<double> f_left = flux(material, velocity);
    vector<double> f_right = flux(material, velocity_next);
    for (int i = 0; i < num_variables; i++)
    {
        velocity_half_step[i]= 0.5*(velocity[i] + velocity_next[i]);
        velocity_half_step[i] += (-0.5)*(d_t/d_x)*(f_right[i]-f_left[i]);
    }
    return flux(material, velocity_half_step);
}

vector<double> FORCE_flux(int material, const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> result = Lax_Friedrichs_flux(material, velocity, velocity_next, d_x, d_t);
    vector<double> R = Richtmyer_flux(material, velocity, velocity_next, d_x, d_t);
    for (int i = 0; i < num_variables; i++)
    {
        result[i] = 0.5*(result[i] + R[i]);
    }
    return result;
}
 

void print_result(const vector<vector<double>> & data_inside, const vector<vector<double>> & data_outside)
{
    for (int i = 0; i < nPoints; i++)
    {
        std::cout << level_set[i+padding] << " ";
        if (level_set[i+padding] < 0.0)
        {
            for (int j = 0; j < num_variables; j++)
            {
                std::cout << data_inside[i+padding][j] << " ";
            }
        }
        else
        {
            for (int j = 0; j < num_variables; j++)
            {
                std::cout << data_outside[i+padding][j] << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout <<std::endl << std::endl;
}

void debug_result(int material, const vector<vector<double>> & data)
{
    vector<double> help_primitive;
    for (auto i = data.begin(); i != data.end(); i++)
    {
        help_primitive = conservative_to_primitive(material, (*i));
        for (int j = 0; j < num_variables; j++)
        {
            std::cout << help_primitive[j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout <<std::endl << std::endl;
}

double vec_len(const vector<double> & v)
{
    int num = v.size();
    double result = 0;
    for (int i = 0; i < num; i++)
    {
        result += (v[i]*v[i]);
    }
    return (sqrt(result));
}

double power(double base, double exponent)
{
    if (base < 0.0)
    {
        std::cout << "negative power" << base << "^" << exponent << "\n";
        std::cout << step << std::endl;
        exit(2);
    }
    else
    return (pow(base, exponent));
}

double pressure(int material, const vector<double> & state)
{
    double rho = state[0];
    double E = state[num_variables-1];
    double v;
    double p;
    if (rho == 0.0)
    {
        /*
        std::cout << "Zero density when calculating the pressure\n";
        exit(2);
        */
       p = ((gas_coef[material] - 1)*E);
    }
    else
    {
        v = state[1]/rho;
        p = ((gas_coef[material] - 1)*(E - (v*v*rho*0.5)));
    }
    /*
    if (p < 0.0)
    {
        std::cout << "Negative pressure\n";
        exit(3);
    }
    */
    return p;
}


double speed_of_sound(int material, const vector<double> & state)
{
    double rho = state[0];
    if (rho == 0)
    {
        std::cout << "Zero density when calculating the speed of sound\n";
        std::cout << rho << "\n";
        exit(2);
    }
    double p = pressure(material, state);
    if (p < 0.0)
    {
        /*
        debug_result(u);
        std::cout << "Negative pressure when calculating the speed of sound\n";
        std::cout << p << "\n";
        exit(3);
        */
       return 0.0;
    }
    else
    {
        return (sqrt(gas_coef[material]*p/rho));
    }
}

double speed(const vector<double> & state)
{
    if (state[0] == 0.0)
    {
        // zero density
        return 0;
    }
    else
    {
        // momentum divided by density
        return (state[1] / state[0]);
    }
}

vector<double> conservative_to_primitive(int material, const vector<double> & state)
{
    vector<double> result;
    result.resize(num_variables);
    result[0] = state[0]; //rho
    result[1] = speed(state);
    result[2] = pressure(material, state);
    return result;
}

vector<double> primitive_to_conservative(int material, const vector<double> & state)
{
    vector<double> result;
    result.resize(num_variables);
    result[0] = state[0]; //rho
    result[1] = state[0]*state[1]; // rho*v
    result[2] = (state[2]/(gas_coef[material] - 1)) + (state[0]*state[1]*state[1]*0.5); // rho*epsilon + 0.5*rho*v^2
    return result;
}

void EvaluateGhostCells(int index, vector<vector<double>>& values_left, int material_left, vector<vector<double>>& values_right, int material_right)
{
    // index holds the position of the last cell in the left material
    if ((index < padding)||(index >= nPoints+padding))
    {
        std::cout << "You cannot insert ghost cells here: " << index << std::endl;
    }
    vector<double> last_conservative;
    last_conservative = values_left[index];
    //std::cout << "copied state\n";
    vector<double> last_primitive;
    last_primitive = conservative_to_primitive(material_left, last_conservative);
    //std::cout << "converted\n";
    vector<double> real, fake;
    real.resize(3);
    fake.resize(3);
    double entropy = last_primitive[2] / power(last_primitive[0], gas_coef[material_left]);
    for (int i = 1; i <= padding; i++)
    {
        if (level_set[index]*level_set[index+i] <= 0.0)
        {
            //std::cout << "Adding a ghost cell\n";
            real = conservative_to_primitive(material_right, values_right[index+i]);
            fake[1] = real[1];
            fake[2] = real[2];
            fake[0] = fake[2] / entropy;
            fake[0] = pow(fake[0], 1.0 / gas_coef[material_left]);
            values_left[index+i] = primitive_to_conservative(material_left, fake);
        }
        else
        {
            std::cout << "You hit the same material again\n";
        }
    }
    last_conservative = values_right[index+1];
    last_primitive = conservative_to_primitive(material_right, last_conservative);
    entropy = last_primitive[2] / pow(last_primitive[0], gas_coef[material_right]);
    for (int i = 0; i < padding; i++)
    {
        if (level_set[index+1]*level_set[index-i] <= 0.0)
        {
            //std::cout << "Adding a ghost cell\n";
            real = conservative_to_primitive(material_left, values_left[index-i]);
            fake[1] = real[1];
            fake[2] = real[2];
            fake[0] = fake[2] / entropy;
            fake[0] = power(fake[0], 1.0 / gas_coef[material_right]);
            values_right[index-i] = primitive_to_conservative(material_right, fake);
        }
        else
        {
            std::cout << "You hit the same material again\n";
        }
    }
}

void FindGhostCells(vector<vector<double>>& values_inside, int material_inside, vector<vector<double>>& values_outside, int material_outside)
{
    for (int index = padding; index < nPoints+padding-1; index++)
    {
        if (level_set[index]*level_set[index+1] <= 0.0)
        {
            //std::cout << "Located boundary at " << index << "\n";
            if (level_set[index] <= 0.0)
            {
                EvaluateGhostCells(index, values_inside, material_inside, values_outside, material_outside);
            }
            else
            {
                EvaluateGhostCells(index, values_outside, material_outside, values_inside, material_inside);
            }
        }
    }

}

double maximal_v(const vector<vector<double>>& values_inside, int material_inside, const vector<vector<double>>& values_outside, int material_outside)
{
    int number = values_inside.size();
    double result = 0;
    double v_c;
    for (int i = padding; i < number-padding; i++)
    {
        if (level_set[i] < 0.0)
        {
            v_c = speed(values_inside[i]);
            v_c = abs(v_c) + speed_of_sound(material_inside, values_inside[i]);
        }
        else
        {
            v_c = speed(values_outside[i]);
            v_c = abs(v_c) + speed_of_sound(material_outside, values_outside[i]);
        }
        if (v_c > result)
        {
            result = v_c;
        }
        if (level_set[i]*level_set[i+1] < 0.0)
        {
            for (int j = i - padding + 1; j <= i+padding; j++)
            {
                v_c = speed(values_outside[j]);
                v_c = abs(v_c) + speed_of_sound(material_outside, values_outside[j]);
                if (v_c > result)
                {
                    result = v_c;
                }
                v_c = speed(values_inside[j]);
                v_c = abs(v_c) + speed_of_sound(material_inside, values_inside[j]);
                if (v_c > result)
                {
                    result = v_c;
                }
            }
        }
    }
    return result;
}

vector<double> flux(int material, const vector<double> & state)
{
    vector<double> result;
    result.resize(3);
    double rho = state[0];
    double rho_v = state[1];
    double E = state[2];
    double p = pressure(material, state);
    double v;

    result[0] = rho_v;
    if (rho == 0.0)
    {
        /*
        std::cout << "Zero density\n";
        std::cout << rho << " " << rho_v << " " << E << "\n";
        exit(2);
        */
       result[1] = p;
       result[2] = 0;
    }
    else
    {
        v = rho_v / rho;
        result[1] = p + rho_v*v;
        result[2]= (E + p)*v;
    }
    
    return result;
}

double (*limiter)(double r);

double minbee(double r)
{
    if (r <= 0.0)
    {
        return 0.0;
    }
    else
    {
        if ((r <= 1.0)&&(r > 0.0))
        {
            return r;
        }
        else
        {
            return (2.0 / (1.0 + r));
        }
    }
}

double zero_limiter(double r)
{
    return 0.0;
}

double superbee(double r)
{
    if (r <= 0.0)
    {
        return 0.0;
    }
    else
    {
        if (r <= 0.5)
        {
            return 2.0*r;
        }
        else
        {
            if (r <= 1)
            {
                return 1.0;
            }
            else
            {
                return (2.0 / (1.0 + r));
            }
        }
    }
}

double Van_Leer(double r)
{
    if (r <= 0.0)
    {
        return 0.0;
    }
    else
    {
        if (r <= 1.0)
        {
            return ((2.0 * r) / (1.0 + r));
        }
        else
        {
            return (2.0 / (1.0 + r));
        }
    }
}
double (*initial_level_set)(double);

double inside_on_the_left(double x)
{
    return (x - 0.5);
}

double helium_inside(double x)
{
    if (x <= 0.5)
    {
        return (0.4 - x);
    }
    else
    {
        return (x - 0.6);
    }
}

void set_initial_helium()
{
    gas_coef = {1.67, 1.4};
    double x;
    vector<double> left, right, helium, help;
    help.resize(3);
    help[0] = 1.3765;
    help[1] = 0.3948;
    help[2] = 1.57;
    left = primitive_to_conservative(1, help);
    help[0] = 1.0;
    help[1] = 0.0;
    help[2] = 1.0;
    right = primitive_to_conservative(1, help);
    help[0] = 0.138;
    helium = primitive_to_conservative(0, help);

    for (int index = 0; index < nPoints +2*padding; index++)
    {
        x = (index - padding)*dx + (0.5*dx);
        level_set[index] = helium_inside(x);
        if (x < 0.25)
        {
            u2[index] = left;
        }
        else
        {
            u2[index] = right;
        }
        u1[index] = helium;
    }
}

void update_level_set_function(const vector<vector<double>> & inside_state, const vector<vector<double>> & outside_state, double space_step, double time_step)
{
    double v, difference;
    vector<double> result;
    result.resize(nPoints + 2*padding);
    for (int index = 1; index < (nPoints + 2*padding -1); index++)
    {
        if (level_set[index] < 0.0)
        {
            v = speed(inside_state[index]);
        }
        else
        {
            v = speed(outside_state[index]);
        }
        if (v < 0.0)
        {
            difference = level_set[index+1] - level_set[index];
        }
        else
        {
            difference = level_set[index] - level_set[index-1];
        }
        result[index] = level_set[index] - (time_step / space_step) * v * difference; 
    }
    if (result[1] > 0.0)
    {
        result[0] = result[1] + space_step;
    }
    else
    {
        result[0] = result[1] - space_step;
    }
    if (result[nPoints + 2*padding - 2] > 0.0)
    {
        result[nPoints + 2*padding-1] = result[nPoints + 2*padding - 2] - space_step;
    }
    else
    {
        result[nPoints + 2*padding - 1] = result[nPoints + 2*padding - 2] + space_step;
    }
    level_set = result;
}

void SLICK_Step(int material, vector<vector<double>> &current, double space_step, double time_step)
{
    
    int next_x = 1;
    // if (direction != 1)
    // {
    //     std::cout << "Invalid direction\n";
    //     exit(3);
    // }
    u_minus = current;
    u_plus = current;
    double r;
    double lim;
    double grad;
    for (int index_x = 1; index_x < (nPoints + 2*padding - 1); index_x++)
    {
        if (current[index_x][num_variables-1] != current[index_x + next_x][num_variables-1])
        {
            r = ((current[index_x][num_variables-1] - current[index_x - next_x][num_variables-1]) / (current[index_x + next_x][num_variables-1] - current[index_x][num_variables-1]));
            lim = limiter(r);
        }
        else
        {
            lim = 0.0;
        }
        for (int j = 0; j < num_variables; j++)
        {
            
            grad = (current[index_x+1][j] - current[index_x-1][j])*0.5;
            u_minus[index_x][j] -= grad*lim*0.5;
            u_plus[index_x][j] += grad*lim*0.5;
        }
    }
    boundary_condition(u_minus);
    boundary_condition(u_plus);
    if (use_halfstep)
    {
        vector<double> flux_plus, flux_minus;
        double change;
        for (int index_x = padding; index_x < nPoints+ padding; index_x++)
        {
            flux_minus = flux(material, u_minus[index_x]);
            flux_plus = flux(material, u_plus[index_x]);
            for (int j = 0; j < num_variables; j++)
            {
                change = (flux_plus[j] - flux_minus[j]) * time_step / space_step / 2.0;
                u_minus[index_x][j] -= change;
                u_plus[index_x][j] -= change;
            }
        }
        boundary_condition(u_minus);
        boundary_condition(u_plus);
    }

    for (int i = 0; i < (nPoints + (2*padding) - 1); i++)
    {
        f[i] = FORCE_flux(material, u_plus[i], u_minus[i + 1], space_step, time_step);
    }
    
    for (int index_x = padding; index_x < (nPoints + padding); index_x++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            current[index_x][j] = current[index_x][j] - ((time_step / space_step)*(f[index_x][j] - f[index_x - 1][j]));
        }
    }
    boundary_condition(current);
}


int main(int argc, char *argv[])
{
    
    //limiter = zero_limiter;
    limiter = Van_Leer;

    dx = (x1 - x0) / nPoints;
    u1.resize(nPoints + 2*padding);
    u2.resize(nPoints + 2*padding);
    u_plus.resize(nPoints + 2*padding);
    u_minus.resize(nPoints + 2*padding);
    level_set.resize(nPoints + 2*padding);
    f.resize(nPoints + 2*padding - 1);

    for (int i = 0; i < (nPoints+(2*padding)-1); i++)
    {
        u1[i].resize(num_variables);
        u2[i].resize(num_variables);
        f[i].resize(num_variables);
        u_plus[i].resize(num_variables);
        u_minus[i].resize(num_variables);
    }
    int last_index = nPoints+(2*padding)-1;
    u1[last_index].resize(num_variables);
    u2[last_index].resize(num_variables);
    u_plus[last_index].resize(num_variables);
    u_minus[last_index].resize(num_variables);

    if (argc < 2)
    {
        std::cout << "Remember to specify the test case\n";
        exit(1);
    }
    if (argv[1][0] == 'H') // Helium
    {
        set_initial_helium();
        tStop = 0.3;
    }
    else
    {
        if (argv[1][0] = 'S') //Simple
        {

        }
        else
        {
            std::cout << "Invalid test case\n";
            exit(1);
        }
    }

    print_result(u1, u2);
    std::cout << "\t\n";

    double t = tStart;
    double dt;
    double a_max;
    step = 1;
    vector<double> flux_left, flux_right;
    flux_left.resize(num_variables);
    flux_right.resize(num_variables);
    while (t < tStop)
    {
        a_max = maximal_v(u1, 0, u2, 1);
        dt = C * dx / a_max;
        // avoid overshooting the final time
        if ((t + dt) > tStop)
        {
            dt = tStop - t;
        }
        //std::cout << "Found time step\n";     

        FindGhostCells(u1, 0, u2, 1);
        //std::cout << "Ghost cells done\n";
        //debug_result(0, u1);
        //debug_result(1, u2);
        update_level_set_function(u1, u2, dx, dt);
        SLICK_Step(0, u1, dx, dt);
        SLICK_Step(1, u2, dx, dt);

        //print_result(u1, u2);
        t = t + dt;
        step += 1;
        //if (step > 2){exit(0);}
    }
    
    print_result(u1, u2);
    //debug_result(u);
    return 0;
}
