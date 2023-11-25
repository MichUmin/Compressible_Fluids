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
double x0 = 0;
double x1 = 1;
double tStart = 0.0;
double tStop;
double a = 1;
double dx = (x1 - x0)/nPoints;
double C = 0.8;
double gas_coef = 1.4;
int num_dimmensions = 1;
int num_variables = 8;
int step = 0;
double vx_l, vx_r, vy_l, vy_r, vz_l, vz_r, rho_l, rho_r, p_l, p_r, Bx_l, Bx_r, By_l, By_r, Bz_l, Bz_r;
bool use_halfstep = true;

vector<vector<double>> u; // u[0] and u[nPoints+1] are padding
                          // u holds rho, rho_vx, rho_vy, rho_vz, Bx, By, Bz, U
vector<vector<double>> f; //f[i] holds the flux from u[i] to u[i+1]
vector<vector<double>> limiters;
vector<vector<double>> u_minus;
vector<vector<double>> u_plus;


void boundary_condition(vector<vector<double>> &values, int pad)
{
    for (int i = 0; i < pad; i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            values[i][j] = values[pad][j];
            values[pad + nPoints + i][j] = values[pad + nPoints - 1][j];
        }
    }
}

vector<double> (flux)(int, const vector<double> &);

vector<double> (*nummerical_flux)(const vector<double>&, const vector<double>&, double, double);

vector<double> Lax_Friedrichs_flux(const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> result;
    result.resize(num_variables);
    vector<double> f_left = flux(velocity);
    vector<double> f_right = flux(velocity_next);
    for (int i = 0; i < num_variables; i++)
    {
        result[i] = (d_x / d_t)*(velocity[i] - velocity_next[i]);
        result[i] += f_left[i];
        result[i] += f_right[i];
        result[i] = result[i] / 2.0;
    }
    return result;
}

vector<double> Richtmyer_flux(const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> velocity_half_step;
    velocity_half_step.resize(num_variables);
    vector<double> f_left = flux(velocity);
    vector<double> f_right = flux(velocity_next);
    for (int i = 0; i < num_variables; i++)
    {
        velocity_half_step[i]= 0.5*(velocity[i] + velocity_next[i]);
        velocity_half_step[i] += (-0.5)*(d_t/d_x)*(f_right[i]-f_left[i]);
    }
    return flux(velocity_half_step);
}

vector<double> FORCE_flux(const vector<double> & velocity, const vector<double> & velocity_next, double d_x, double d_t)
{
    vector<double> result = Lax_Friedrichs_flux(velocity, velocity_next, d_x, d_t);
    vector<double> R = Richtmyer_flux(velocity, velocity_next, d_x, d_t);
    for (int i = 0; i < num_variables; i++)
    {
        result[i] = 0.5*(result[i] + R[i]);
    }
    return result;
}

void set_initial_value(vector<vector<double>> &vel, int pad)
{
    double rho_vx_l = rho_l*vx_l;
    double rho_vy_l = rho_l*vy_l;
    double rho_vz_l = rho_l*vz_l;
    double E_l = p_l/(gas_coef - 1.0);
    E_l += (rho_l*vx_l*vx_l*0.5);
    E_l += (rho_l*vy_l*vy_l*0.5);
    E_l += (rho_l*vz_l*vz_l*0.5);
    E_l += (Bx_l*Bx_l*0.5);
    E_l += (By_l*By_l*0.5);
    E_l += (Bz_l*Bz_l*0.5);
    double rho_vx_r = rho_l*vx_r;
    double rho_vy_r = rho_l*vy_r;
    double rho_vz_r = rho_l*vz_r;
    double E_r = p_r/(gas_coef - 1.0);
    E_r += (rho_r*vx_r*vx_r*0.5);
    E_r += (rho_r*vy_r*vy_r*0.5);
    E_r += (rho_r*vz_r*vz_r*0.5);
    E_r += (Bx_r*Bx_r*0.5);
    E_r += (By_r*By_r*0.5);
    E_r += (Bz_r*Bz_r*0.5);
    for (int i = 0; i <= nPoints/2 ; i ++)
    {
        vel[i+pad][0] = rho_l;
        vel[i+pad][1] = rho_vx_l;
        vel[i+pad][2] = rho_vy_l;
        vel[i+pad][3] = rho_vz_l;
        vel[i+pad][4] = Bx_l;
        vel[i+pad][5] = By_l;
        vel[i+pad][6] = Bz_l;
        vel[i+pad][7] = E_l;
        vel[nPoints-+pad][0] = rho_r;
        vel[nPoints-i+pad][1] = rho_vx_r;
        vel[nPoints-i+pad][2] = rho_vy_r;
        vel[nPoints-i+pad][3] = rho_vz_r;
        vel[nPoints-i+pad][4] = Bx_r;
        vel[nPoints-i+pad][5] = By_r;
        vel[nPoints-i+pad][6] = Bz_r;
        vel[nPoints-i+pad][7] = E_r;
    }
    if (nPoints % 2 == 1)
    {
        vel[(nPoints/2)+1+pad][0] = 0.5*(rho_r + rho_l);
        vel[(nPoints/2)+1+pad][1] = 0.5*(rho_vx_r + rho_vx_l);
        vel[(nPoints/2)+1+pad][2] = 0.5*(rho_vy_r + rho_vy_l);
        vel[(nPoints/2)+1+pad][3] = 0.5*(rho_vz_r + rho_vz_l);
        vel[(nPoints/2)+1+pad][4] = 0.5*(Bx_r + Bx_l);
        vel[(nPoints/2)+1+pad][5] = 0.5*(By_r + By_l);
        vel[(nPoints/2)+1+pad][6] = 0.5*(Bz_r + Bz_l);
        vel[(nPoints/2)+1+pad][7] = 0.5*(E_r + E_l);
    }
    boundary_condition(vel, pad);
}

void print_result(const vector<vector<double>> & data)
{
    for (int i = 0; i < nPoints; i++)
    {
        //std::cout << i << " ";
        for (int j = 0; j < num_variables; j++)
        {
            //std::cout << j << " ";
            std::cout << data[i+padding][j] << " ";
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

double pressure(const vector<double> & state)
{
    double rho = state[0];
    double E = state[num_variables-1];
    vector<double> v;
    v.resize(3);
    double p;
    double B_2;
    if (rho == 0.0)
    {
        /*
        std::cout << "Zero density when calculating the pressure\n";
        exit(2);
        */
       p = ((gas_coef - 1)*E);
    }
    else
    {
        B_2 = 0.5*(state[4]*state[4] + state[5]*state[5] + state[6]*state[6]);
        for (int j = 0; j < 3; j++)
        {
            v[j] = state[1+j]/rho;
        }
        p = ((gas_coef - 1)*(E - (square(vec_len(v))*rho*0.5) - B_2));
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


double speed_of_sound(const vector<double> & state)
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
        return (sqrt(gas_coef*p/rho));
    }
}

double maximum(double a, double b, double c)
{
    double result = abs(a);
    if (abs(b) > result)
    {
        result = abs(b);
    }
    if (abs(c) > result)
    {
        result = abs(c);
    }
    return result;
}

double maximal_v(const vector<vector<double>>& values)
{
    int number = values.size();
    double rho;
    double result = 0;
    double v_c, c_s, c_a, c_f;
    vector<double> v;
    v.resize(3);
    for (int i = 0; i < number; i++)
    {
        rho = values[i][0];
        c_s = speed_of_sound(values[i]);
        c_a = values[i][4] / sqrt(values[i][0]);
        c_f = sqrt(0.5*(c_s*c_s + c_a*c_a + abs(c_s*c_s - c_a*c_a)));

        for (int j = 0; j < 3; j++)
        {
            v[j] = values[i][1+j]/rho;
        }
        v_c = vec_len(v); // rho_v / rho
        v_c += c_f;
        if (v_c > result)
        {
            result = v_c;
        }
    }
    return result;
}

vector<double> flux(int direction, const vector<double> & state)
{
    if ((direction < 1)||(direction > num_dimmensions))
    {
        std::cout << "Invalid direction when evaluating flux\n"
        exit(4);
    }
    vector<double> result;
    result.resize(num_variables);
    double rho = state[0];
    double p = pressure(state);
    double B_2 = 0.5*(state[4]*state[4] + state[5]*state[5] + state[6]*state[6]);
    vector<double> v;
    v.resize(3);

    result[0] = state[direction];
    if (rho == 0.0)
    {
        std::cout << "Zero density\n";
        std::cout << state << std::endl;
        exit(2);
    }
    else
    {
        v[0] = state[1] / rho;
        v[1] = state[2] / rho;
        v[2] = state[3] / rho;
        for (int i = 0; i < 3; i++)
        {
            result[i+1] = state[direction]*v[i] - state[direction+3]*state[i+4];
            result[i+4] = state[i+4]*state[direction] - state[direction+3]*state[i+1];
        }
        result[direction] += p + B_2;
        result[num_variables-1] = (state[num_variables-1] + p + B_2) * v[direction-1];
        for (int i = 0; i < 3; i++)
        {
            result[num_variables-1] -= state[direction + 3]*v[i]*state[i+4];
    }
    
    return result;
}

void evaluate_flux(const vector<vector<double>> & values_left, const vector<vector<double>> & values_right, vector<vector<double>> &flu, double space_step, double time_step)
{
    for (int i = 0; i < (nPoints + (2*padding) - 1); i++)
    {
        //std::cout << i << "\n";
        //std::cout << values_left[i+1][0] << "\n";
        flu[i] = nummerical_flux(values_right[i], values_left[i+1], space_step, time_step);
    }
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

void (*find_limiters)(const vector<vector<double>> & between, vector<vector<double>> & limiters);

void find_limiters_energy(const vector<vector<double>> & values, vector<vector<double>> & limiters)
{
    double r;
    double val;
    for (int i = padding; i < nPoints+padding; i++)
    {
        if (values[i+1][num_variables-1] != values[i][num_variables-1])
            {
                //std::cout << between[i + padding][j] << std::endl;
                r = ((values[i][num_variables] - values[i-1][num_variables-1]) / (values[i+1][num_variables-1] - values[i][num_variables-1]));
                val = limiter(r);
            }
            else
            {
                val = 0.0;
            }
        for (int j = 0; j < num_variables; j++)
        {
            limiters[i][j] = val;
        }
    }
    for (int i = 0; i < padding; i++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            limiters[i][j] = 0.0;
            limiters[i + nPoints + padding][j] = 0.0;
        }
    }
}


void SLICK_Step(int direction, vector<vector<double>> &current, vector<vector<double>> &fluxes, double space_step, double time_step)
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
                flux_minus = flux(direction, u_minus[index_x]);
                flux_plus = flux(direction, u_plus[index_x]);
                for (int j = 0; j < num_variables; j++)
                {
                    change = (flux_plus[j] - flux_minus[j]) * time_step / space_step / 2.0;
                    u_minus[index_x][j] -= change;
                    u_plus[index_x][j] -= change;
                }
            }
        }
        boundary_condition(u_minus, padding);
        boundary_condition(u_plus, padding);
    }

    for (int i = 0; i < (nPoints + (2*padding) - 1); i++)
    {
        fluxes[i] = nummerical_flux(direction, u_plus[i], u_minus[i + 1], space_step, time_step);
    }
    
    for (int index_x = padding; index_x < (nPoints + padding); index_x++)
    {
        for (int j = 0; j < num_variables; j++)
        {
            current[index_x][j] = current[index_x][j] - ((time_step / space_step)*(fluxes[index_x][j] - fluxes[index_x - 1][j]));
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
