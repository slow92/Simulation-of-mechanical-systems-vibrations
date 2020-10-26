#include "Simulations_v1.h"
#include "Templates.h"

#include <cmath>
#include <fstream>



Simulations_v1::Simulations_v1() {}

Simulations_v1::~Simulations_v1() {}


bool Simulations_v1::one_freedom_degree()
{
    reset_vector();

    double m2,Au1,c21,b21,p;
    long N;

    Templates::insert_user_data<double>("Insert mass [kg]: ",m2,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert amplitude of unit jump [m]: ",Au1,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert elasticity coefficient: ",c21,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert dumping factor: ",b21,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert dumping exponent: ",p,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert integration step [s]: ",dt,"Error: data have be floating point");
    Templates::insert_user_data<long>("Insert number of integration steps: ",N,"Error: data have be floating point");

    init_vector(1,N);

    double Fs21,Ft21;
    double u2,v2=0,v1=0,a2;
    double Q2=-m2*Gravity;

    double u1=Au1;
    u2=Q2/c21;

    for(int i=0;i<N;i++)
    {
        Fs21 = -c21*(u2-u1);
        Ft21 = -b21*pow(abs(v2-v1),p)*(std::signbit(v2-v1)? -1:1);

        a2=(Ft21+Fs21+Q2)/m2;
        v2=v2+a2*dt;
        u2=u2+v2*dt;

        v_u2.push_back(u2);
        v_v2.push_back(u2);
        v_a2.push_back(u2);
    }

    return true;
}


bool Simulations_v1::two_freedom_degree()
{
    reset_vector();

    double Au1,p;
    double m2,c21,b21;
    double m3,c32,b32;
    long N;

    Templates::insert_user_data<double>("Insert mass m2 [kg]: ",m2,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert mass m3 [kg]: ",m3,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert amplitude of unit jump [m]: ",Au1,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert elasticity coefficient c21: ",c21,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert elasticity coefficient c32: ",c32,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert dumping factor b21: ",b21,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert dumping exponent: ",p,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert dumping factor b32: ",b32,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert integration step [s]: ",dt,"Error: data have be floating point");
    Templates::insert_user_data<long>("Insert number of integration steps: ",N,"Error: data have be floating point");

    init_vector(2,N);

    double Fs21, Ft21;
    double Fs32, Ft32;
    double v1=0;
    double v2=0,a2;
    double v3=0,a3;

    double Q2=-m2*Gravity;
    double Q3=-m3*Gravity;

    double u1=Au1;
    double u2=(Q2+Q3)/c21;
    double u3=u2+Q3/c32;

    for(int i=0;i<N;i++)
    {
        Fs21 = -c21*(u2-u1);
        Fs32 = -c32*(u3-u2);
        Ft21 = -b21*pow(abs(v2-v1),p)*(std::signbit(v2-v1)? -1:1);
        Ft32 = -b32*(v3-v2);

        a2=(Ft21+Fs21+Q2-Ft32-Fs32)/m2;
        a3=(Ft32+Fs32+Q3)/m3;
        v2=v2+a2*dt;
        v3=v3+a3*dt;
        u2=u2+v2*dt;
        u3=u3+v3*dt;

        v_u2.push_back(u2);
        v_v2.push_back(u2);
        v_a2.push_back(u2);
        v_u3.push_back(u2);
        v_v3.push_back(u2);
        v_a3.push_back(u2);
    }

    return true;
}


bool Simulations_v1::tree_freedom_degree()
{
    reset_vector();

    double Au1,p;
    double m2,c21,b21;
    double m3,c32,b32;
    double m4,c43,b43;
    long N;

    Templates::insert_user_data<double>("Insert mass m2 [kg]: ",m2,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert mass m3 [kg]: ",m3,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert mass m4 [kg]: ",m4,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert amplitude of unit jump [m]: ",Au1,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert elasticity coefficient c21: ",c21,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert elasticity coefficient c32: ",c32,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert elasticity coefficient c43: ",c43,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert dumping factor b21: ",b21,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert dumping exponent: ",p,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert dumping factor b32: ",b32,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert dumping factor b43: ",b43,"Error: data have be floating point");
    Templates::insert_user_data<double>("Insert integration step [s]: ",dt,"Error: data have be floating point");
    Templates::insert_user_data<long>("Insert number of integration steps: ",N,"Error: data have be floating point");

    init_vector(3,N);

    double Fs21, Ft21;
    double Fs32, Ft32;
    double Fs43, Ft43;
    double v1=0;
    double v2=0,a2;
    double v3=0,a3;
    double v4=0,a4;

    double Q2=-m2*Gravity;
    double Q3=-m3*Gravity;
    double Q4=-m4*Gravity;

    double u1=Au1;
    double u2=(Q2+Q3+Q4)/c21;
    double u3=u2+(Q3+Q4)/c32;
    double u4=u3+Q4/c43;

    for(int i=0;i<N;i++)
    {
        Fs21 = -c21*(u2-u1);
        Fs32 = -c32*(u3-u2);
        Fs43 = -c43*(u4-u3)*(std::signbit(u3-u4)? 0:1);

        Ft21 = -b21*pow(abs(v2-v1),p)*(std::signbit(v2-v1)? -1:1);
        Ft32 = -b32*(v3-v2);
        Ft43 = -b43*(v4-v3)*(std::signbit(v3-v4)? 0:1);

        a2=(Ft21+Fs21+Q2-Ft32-Fs32)/m2;
        a3=(Ft32+Fs32+Q3-Ft43-Fs43)/m3;
        a4=(Ft43+Fs43+Q4)/m4;
        v2=v2+a2*dt;
        v3=v3+a3*dt;
        v4=v4+a4*dt;
        u2=u2+v2*dt;
        u3=u3+v3*dt;
        u4=u4+v4*dt;

        v_u2.push_back(u2);
        v_v2.push_back(u2);
        v_a2.push_back(u2);
        v_u3.push_back(u2);
        v_v3.push_back(u2);
        v_a3.push_back(u2);
        v_u4.push_back(u2);
        v_v4.push_back(u2);
        v_a4.push_back(u2);
    }

    return true;
}


bool Simulations_v1::save_data()
{
    if(v_u2.empty())
    {
        std::cout<<"No data to save!"<<std::endl;
        return false;
    }
    std::fstream file;
    std::string file_name;
    Templates::insert_user_data<std::string>("Save as: ",file_name,"");
    file.open(file_name.c_str(),std::ios::out);

    if((!v_u2.empty())&&v_u3.empty()&&v_u4.empty())
    {
        file<<"time[s] u2[m] v2[m/s] a2[m/s2]"<<std::endl;
        for(int i=0;i<static_cast<int>(v_u2.size());i++)
        {
            file<<i*dt<<" "<<v_u2[i]<<" "<<v_v2[i]<<" "<<v_a2[i]<<std::endl;
        }
    }
    if((!v_u2.empty())&&(!v_u3.empty())&&v_u4.empty())
    {
        file<<"time[s] u2[m] u3[m] v2[m/s] v3[m/s] a2[m/s2] a3[m/s2]"<<std::endl;
        for(int i=0;i<static_cast<int>(v_u2.size());i++)
        {
            file<<i*dt
            <<" "<<v_u2[i]<<" "<<v_u3[i]
            <<" "<<v_v2[i]<<" "<<v_v3[i]
            <<" "<<v_a2[i]<<" "<<v_a3[i]<<std::endl;
        }
    }
    if((!v_u2.empty())&&(!v_u3.empty())&&(!v_u4.empty()))
    {
        file<<"time[s] u2[m] u3[m] u4[m] v2[m/s] v3[m/s] v4[m/s] a2[m/s2] a3[m/s2] a4[m/s2]"<<std::endl;
        for(int i=0;i<static_cast<int>(v_u2.size());i++)
        {
            file<<i*dt
            <<" "<<v_u2[i]<<" "<<v_u3[i]<<" "<<v_u4[i]
            <<" "<<v_v2[i]<<" "<<v_v3[i]<<" "<<v_v4[i]
            <<" "<<v_a2[i]<<" "<<v_a3[i]<<" "<<v_a4[i]<<std::endl;
        }
    }

    file.close();
return true;
}


bool Simulations_v1::reset_vector()
{
    if(!v_u2.empty()) v_u2.clear();
    if(!v_u3.empty()) v_u3.clear();
    if(!v_u4.empty()) v_u4.clear();
    if(!v_v2.empty()) v_v2.clear();
    if(!v_v3.empty()) v_v3.clear();
    if(!v_v4.empty()) v_v4.clear();
    if(!v_a2.empty()) v_a2.clear();
    if(!v_a3.empty()) v_a3.clear();
    if(!v_a4.empty()) v_a4.clear();
    return true;
}


bool Simulations_v1::init_vector(const int& freedom_degree_number, const long& iteration_number)
{
    if(freedom_degree_number>=1)
    {
        v_u2.reserve(iteration_number);
        v_v2.reserve(iteration_number);
        v_a2.reserve(iteration_number);
    }
    if(freedom_degree_number>2)
    {
        v_u3.reserve(iteration_number);
        v_v3.reserve(iteration_number);
        v_a3.reserve(iteration_number);
    }
    if(freedom_degree_number>=3)
    {
        v_u4.reserve(iteration_number);
        v_v4.reserve(iteration_number);
        v_a4.reserve(iteration_number);
    }
    return true;
}
