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

int nPoints = 100;
int padding = 1;
int nGhost = 3*padding;
double x0 = 0;
double x1 = 1;
double tStart = 0.0;
double tStop;
double a = 1;
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
        if (level_set[i+padding] < 0.0)
        {
            for (int j = 0; j < num_variables; j++)
            {
                //std::cout << j << " ";
                std::cout << data_inside[i+padding][j] << " ";
            }
        }
        else
        {
            for (int j = 0; j < num_variables; j++)
            {
                //std::cout << j << " ";
                std::cout << data_outside[i+padding][j] << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout <<std::endl << std::endl;
}

void debug_result(const vector<vector<double>> & data)
{
    for (auto i = data.begin(); i != data.end(); i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            std::cout << (*i)[j] << " ";
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

double pressure(int material, vector<double> & state)
{
    double rho = state[0];
    double E = state[num_variables-1];
    double v;
    v.resize(3);
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
        p = ((gas_coef[material] - 1)*(E - (square(v))*rho*0.5));
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
    double p = pressure(state);
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

double maximal_v(const vector<vector<double>>& values_inside, const vector<vector<double>>& values_outside)
{
    int number = values.size();
    double result = 0;
    double v_c;
    for (int i = padding; i < number-padding; i++)
    {
        if (level_set[i] < 0.0)
        {
            v_c = speed(values_inside[i]);
            v_c = abs(v_c) + speed_of_sound(values_inside[i]);
        }
        else
        {
            v_c = speed(values_outside[i]);
            v_c = abs(v_c) + speed_of_sound(values_outside[i]);
        }
        if (v_c > result)
        {
            result = v_c;
        }
        if (level_set[i]*level_set[i+1] < 0.0)
        {
            for (int j = i - nGhost + 1; j <= i+nGhost; j++)
            {
                v_c = speed(values_outside[j]);
                v_c = abs(v_c) + speed_of_sound(values_outside[j]);
                if (v_c > result)
                {
                    result = v_c;
                }
                v_c = speed(values_inside[j]);
                v_c = abs(v_c) + speed_of_sound(values_inside[j]);
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
void update_level_set_function(const vector<vector<double>> & inside_state, const vector<vector<double>> & outside_state, double time_step, double space_step)
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
        result[index] = level_set[index] - (time_step / space_step) * difference; 
    }
    if (result[1] > 0.0)
    {
        result[0] = result[1] + space_step;
    }
    else
    {
        result[0] = result[1] - space_step;
    }
    if (result[nPoints + 2*padding - 1] > 0.0)
    {
        result[nPoints + 2*padding] = result[nPoints + 2*padding - 1] - space_step;
    }
    else
    {
        result[nPoints + 2*padding] = result[nPoints + 2*padding - 1] + space_step;
    }
    level_set = result;
}

void SLICK_Step(int direction, int material, vector<vector<double>> &current, double space_step, double time_step)
{
    
    int next_x = 1;
    if (direction != 1)
    {
        std::cout << "Invalid direction\n";
        exit(3);
    }
    u_minus = current;
    u_plus = current;
    double r;
    double lim;
    double grad;
    for (int index_x = 1; index_x < (nPoints + 2*padding - 1); index_x++)
    {
        if (current[index_x][num_variables-1] != current[index_x + next_x][num`-1])
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
    boundary_condition(u_minus, padding);
    boundary_condition(u_plus, padding);
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
        boundary_condition(u_minus, padding);
        boundary_condition(u_plus, padding);
    }

    for (int i = 0; i < (nPoints + (2*padding) - 1); i++)
    {
        f[i] = nummerical_flux(material, u_plus[i], u_minus[i + 1], space_step, time_step);
    }
    
    for (int index_x = padding; index_x < (nPoints + padding); index_x++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            current[index_x][j] = current[index_x][j] - ((time_step / space_step)*(fs[index_x][j] - f[index_x - 1][j]));
        }
    }
    boundary_condition(current, padding);
}


int main(int argc, char *argv[])
{
    nummerical_flux = FORCE_flux;
    //limiter = zero_limiter;
    limiter = Van_Leer;
    find_limiters = find_limiters_energy;

    if (argc < 2)
    {
        std::cout << "Remember to specify the test case you are running\n";
        exit(1);
    }
    // the first argument is the test case
    if (argv[1][4] == '1') // Test1
    {
        rho_l = 1.0;
        v_l = 0.0;
        p_l = 1.0;
        rho_r = 0.125;
        v_r = 0.0;
        p_r = 0.1;
        tStop = 0.25;
    }
    else
    {
        if (argv[1][4] == '2') // Test2
        {
            rho_l = 1.0;
            v_l = -2.0;
            p_l = 0.4;
            rho_r = 1.0;
            v_r = 2.0;
            p_r = 0.4;
            tStop = 0.15;
        }
        else
        {
            if (argv[1][4] == '3') // Test3
            {
                rho_l = 1.0;
                v_l = 0.0;
                p_l = 1000.0;
                rho_r = 1.0;
                v_r = 0.0;
                p_r = 0.01;
                tStop = 0.012;
            }
            else
            {
                if (argv[1][4] == '4') // Test4
                {
                    rho_l = 1.0;
                    v_l = 0.0;
                    p_l = 0.01;
                    rho_r = 1.0;
                    v_r = 0.0;
                    p_r = 100.0;
                    tStop = 0.035;
                }
                else
                {
                    if (argv[1][4] == '5') // Test5
                    {
                        rho_l =  5.99924;
                        v_l =  19.5975;
                        p_l =  460.894;
                        rho_r = 5.99242;
                        v_r = -6.19633;
                        p_r =  46.095;
                        tStop = 0.035;
                    }
                    else
                    {
                        std::cout << "invalid test case\n";
                        exit(1);
                    }
                }
            }
        }
    }
    

    dx = (x1 - x0) / nPoints;
    u.resize(nPoints + 2*padding);
    u_left.resize(nPoints + 2*padding);
    u_right.resize(nPoints + 2*padding);
    uPlusOne_left.resize(nPoints + 2*padding);
    uPlusOne_right.resize(nPoints + 2*padding);
    f.resize(nPoints + 2*padding - 1);
    center_deltas.resize(nPoints + 2*padding);
    limiters.resize(nPoints + 2*padding);
    between_deltas.resize(nPoints + 2*padding - 1);

    for (int i = 0; i < (nPoints+(2*padding)-1); i++)
    {
        u[i].resize(num_variables);
        f[i].resize(num_variables);
        between_deltas[i].resize(num_variables);
        limiters[i].resize(num_variables);
        center_deltas[i].resize(num_variables);
        u_left[i].resize(num_variables);
        u_right[i].resize(num_variables);
        uPlusOne_left[i].resize(num_variables);
        uPlusOne_right[i].resize(num_variables);
    }
    int last_index = nPoints+(2*padding)-1;
    u[last_index].resize(num_variables);
    limiters[last_index].resize(num_variables);
    center_deltas[last_index].resize(num_variables);
    u_left[last_index].resize(num_variables);
    u_right[last_index].resize(num_variables);
    uPlusOne_left[last_index].resize(num_variables);
    uPlusOne_right[last_index].resize(num_variables);

    set_initial_value(u, padding);
    print_result(u);
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
        a_max = maximal_v(u);
        dt = C * dx / a_max;
        // avoid overshooting the final time
        if ((t + dt) > tStop)
        {
            dt = tStop - t;
        }     
        //std::cout << dt << std::endl; 

        find_deltas(u, between_deltas, center_deltas);
        //std::cout << "Deltas ";
        //find_limiters(between_deltas, limiters);
        find_limiters(u, limiters);
        //std::cout << "limiters\n";
        //debug_result(limiters);

        double  change;
        for (int i = 0; i < (nPoints + (2*padding)); i++)
        {
            for (int j = 0; j < num_variables; j++)
            {
                change = 0.5*(limiters[i][j] * center_deltas[i][j]);
                //std::cout << change << " ";
                u_left[i][j] = u[i][j] - change; 
                u_right[i][j] = u[i][j] + change;
            }
            //std::cout << std::endl;
        }
        if (use_halfstep)
        {
            for (int i = 0; i < nPoints; i++)
            {
                flux_left = flux(u_left[i + padding]);
                flux_right = flux(u_right[i + padding]);
                for (int j = 0; j < num_variables; j++)
                {
                    change = (flux_right[j] - flux_left[j]) * dt / dx / 2.0;
                    uPlusOne_left[i+padding][j] = u_left[i+padding][j] - change;
                    uPlusOne_right[i+padding][j] = u_right[i+padding][j] - change;
                }
            }
            for (int i = 0; i < padding; i++)
            {
                for (int j = 0; j < num_variables; j++)
                {
                    uPlusOne_left[i][j] = u_left[i][j];
                    uPlusOne_right[i][j] = u_right[i][j];
                    uPlusOne_right[nPoints+padding + i][j] = u_right[nPoints+padding+i][j];
                    uPlusOne_left[nPoints+padding + i][j] = u_left[nPoints+padding+i][j];
                }
            }
            u_left = uPlusOne_left;
            u_right = uPlusOne_right;  
        }
        //boundary_condition(u_left, padding);
        //boundary_condition(u_right, padding);
        //debug_result(u_left);
        //debug_result(u_right);

        evaluate_flux(u_left, u_right, f, dx, dt);

        //debug_result(f);

        //std::cout << "Flux evaluated\n";

        u = PDE_Step(u, f, dt, padding);

        //print_result(u);
        //std::cout << "PDE Step done\n";

        t = t + dt;
        step += 1;
    }

    print_result(u);
    //debug_result(u);
    return 0;
}
